/**
 * @file ParallelSketcher.cc (v8.6.2)
 * @author CM
 * @brief Parallel implementation of a mergeable Sketch (Reduce local sketches into a final one)
 * 
 * 
 * 1) Initially, all processes get their data sample from a common file (VALIDATION MODE)
 * or self-generated distribution with a different seed for each process (RUN MODE).
 * 
 * @note: Since we don't plan to use a parallel file system, use canonical I/O
 * @note: To perform scalability tests, don't use a memory local_buffer to store generated data
 * 
 * 2) Then, each process populates a local Sketch.
 * @note: In order to take advantage of keys ordering, we still make use of a std::map<int, int>
 * 
 * In VALIDATION MODE, to take into account positive and negative items, 2 sketches are used:
 * - one for positive values and the special bucket B* for 0 and near-0 values
 * - one for negative ones
 * 
 * @note In the negative sketch, keys are calculated wrt the absolute value of each item:
 * keys monotonicity is inverse of that of the values 
 * (i.e. item = -100 will have a key greater than that of item = -10, 
 * lowest (negative) items have greatest (positive) keys) 
 *  
 * @note In RUN MODE, only positive values are generated and processed (1 positive sketch only)
 * 
 * @note Since in RUN MODE we only want to deal with positive items, while in 
 * VALIDATION MODE we will deal also with negative items, we can 
 * - handle a pair of sketches 
 * - or have only a positive one, managing the values generated in RUN Mode to be sure items are ALL positives
 * 
 * \begin{verbatim}

    std::map<int,int> posiSketch;               
    int posipop = 0;                            
    int posibins = 0;                           

    #ifdef VALIDATE
        std::map<int,int> negaSketch;              
        int negapop = 0;                         
        int negabins = 0;                          
    #endif
  
 * \end{verbatim}
 * 
 * @note When the sum of positive and negative bins exceeds the fixed sketch bound, a uniform collapse procedure is invoked.
 * The collapse works separately on the two sketches and the sum of bins is taken into account at each iteration.
 *
 * @note For the DDOG version, collapse takes care of both sketches, but it merges only pair of buckets.
 * If we collapse left, we start from the negative one, merging the 2 leftmost bins, until the number of bins is the negative sketch is equal to 1.
 * Then we start collapsing the leftmost 2 bins in the positive sketch (but not the B* bucket)
 *  Reversly when we collapse right.
 * 
 * 
 * 3) When sketches are filled, a merge procedure is started.
 * 
 * Merge requires 2 steps:
 * - in the first steps all sketches are collapsed to the common MAX ALPHA
 * - in the second step, the merge of buckets is performed
 * 
 * To perform Reduction, an MPI_PACKED data is serialized from the original local sketches.
 * Since in each Reduction step, the merge of posi and nega sketches can exceed the sketch bound
 * a collapse procedure must be implemented to work on array of Buckets instead of std::map<,>
 * 
 * @note In the final phase, keys in the negative sketch must be sorted in the reverse order before querying quantiles.
 * 
 * 4) P0 performs queries to estimate quantiles on the unique global reduced sketch
 * 
 * 5) P0 logs to stderr in csv format the running results
 * 
 * 
 * @note - version 2
 *          In version 2, specialized the VALIDATION and RUN mode (respectively -DCHECK and -DTEST)
 *          - VALIDATION to compare sequential results with parallel results
 *          - RUN mode to check scalability and throughput
 *          In RUN mode, added the possibility to generate input distribution on the fly and to consider only positive values
 * 
 * @note - version 6
 *          reverted specialized computation to use both sketch pairs (posi and nega)
 * 
 * @note - version 7
 *          Added the original DDSketch of Datadog
 * 
 * @note - version 8
 *          Try to improve code performances
 * 
 * @note - version 8.1
 *          In RUN Mode items are generated but not allocated in a memory local_buffer
 * 
 * @note - version 8.3
 *          Testing VALIDATION Mode and added exact quantiles computation
 * 
 * @note - version 8.4
 *          improved merging of DDOG sketches in MPI_Reduce()
 * 
 * @note - version 8.6.1
 *          .) added computation of the final alpha value for the thrash bin in DDOG (ONLY FOR LOW BINS COLLAPSE PROCEDURE)
 *          .) added command line seed configuration: (if) the user defined seed is used to initialize the random generator 
 *              (otherwise, the current timestamp is used)
 *          .) the Beta distribution is added
 *  
 * @note - version 8.6.2: added the LogNormal distribution 
 * 
 * 
 * @note 
 * NULLBOUND (the bound used to fill the special bucket B* for 0 and near-0 items) is evaluated at the beginning and never changed during collapses 
 * 
 * @todo: extend the DDOG alphaErr computation for the HIGH BINS COLLAPSE PROCEDURE
 * 
 * 
 */

#include "Merger.h"
#include "Quantiles.h"                          // for quickselect


#include <limits>                               // to compute lowest item value (DDOG)
#include <cstddef>


char VERSION[] = "ParallelSketcher (v8.6.1)";   // code version for logging

//*** QQ contains the quantiles to estimate
double QQ[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0};
int nQ = sizeof(QQ)/sizeof(double);             // size of the QQ array


double NULLBOUND;                               // NULLBOUND is the bound (on the right) of the special bucket B* for 0 and near-0 values 
int Gbuffer_size = 0;                           // the size of the MPI_PACKED exchanged during MPI_Reduce

typedef struct CC{
    int collapse;
    long items;
}CollapseCounter;


CollapseCounter monitorCollapses[1];                // for each collapse, the number of received items until now is recorded
int m_idx = 0;                                      // current position in the monitorCollapses array




int main(int argc, char *argv[]) {

    //*************** 1. Check the minimum parameters are specified
    if (argc < MANDATORY_PARAMS) {   
        //at least (-a -b -f) OR (-a -b -d)
        std::cout << "Bad command line configuration\n" << std::endl;
        usage(argv[0]);
        exit(EPARAMS);
    }//fi argc
    
    //*************** 2. Read command line 

    Params configParam; // configParam collects command line arguments                     

    int res = check_params(argc, argv, &configParam);
    if (res) {
        std::cout << "Unrecognized command line configuration\n" << std::endl;
        usage(argv[0]);
        exit(EPARAMS);
    }//fi res

    #ifdef VERBOSE
        debugParams(&configParam);
    #endif

    //*************** 3. Command line configuration is ok, check that mode is coherent with params

    #ifdef RUN
        if (configParam.dtype == UNKNOWN) {
            std::cout << "\tATTENTION:\n\tYou are in RUN mode, input is expected to be self-generated from a given distribution\n" << std::endl;
            exit(EDATA);
        }
        if (configParam.processInputLen <= 0) {
            std::cout << "\tATTENTION:\n\tYou are in RUN mode, you must specify the number of items to generate from the given distribution\n" << std::endl;
            exit(EDATA);
        }
    #else //VALIDATE
        if (!configParam.fileflag) {
            std::cout << "\tATTENTION:\n\tYou are in VALIDATION mode, input is expected from given a binary file\n" << std::endl;
            exit(EDATA);
        }


        //*** 3.1 in VALIDATION mode exact quantiles computation is performed for relative error tests (if -e is specified on the command line)
        int exactLen = 0;           // to compute the exact quantiles
        double *exactStats = NULL;  // to local_buffer the entire file items

    #endif //RUN vs VALIDATE

    
    //*************** 4. Set up logfile for each of the running process (based on its pid)

    #ifdef LOG
        FILE *f_log = logProcess(&configParam);
    #else
        FILE *f_log = NULL;
    #endif
    
    

    //*************** 5. MPI data and MPI initialization
    int ecode = MPI_SUCCESS;            // MPI error code returned by the APIs

    int comm_size;                      // communicator size
    int process_rank;                   // rank of the process 
    int root_process = 0;               // process root, managing the reduce 
    int pidRoot = 0;                    // pid of the root process   

    ecode = warming_up_MPI(&argc, &argv, &comm_size, &process_rank);    
    if (ecode != MPI_SUCCESS) {
        if (f_log) {
            fprintf(f_log, "Error initializing MPI... shutting down\n");
            fclose(f_log);
        }
        exit(ecode);
    }//fi ecode
    
    
    //*************** 6. DDSketch data
    
    std::map<int, long> posiSketch;             // sketch for positive and null items
    int posibins = 0;                           // number of bins in the positive sketch respectively
    long posipop = 0;                           // population in the positive sketch
    
    std::map<int, long> negaSketch;             // sketch for negative items
    int negabins = 0;                           // number of bins in the negative sketch 
    long negapop = 0;                           // population in the negative sketch
    
    double p_alpha = configParam.initial_alpha; // process's actual α  
    double gamma = (1+p_alpha)/(1-p_alpha);     // γ = (1+α)/(1-α),  γ> 1 ∀α
    double base = std::log10(gamma);            // logγ
    NULLBOUND = pow(gamma, -MIN_KEY);           // define here NULLBOUND since we don't want it to change in collapses
    
    
    int process_collapses = 0;                  // total number of collapses executed on the overall sketch by current process
    
    double xmin = 0.0;                          // minimum value among the items in the stream
    xmin = std::numeric_limits<double>::max();  // initialize to the max double value
    

    //*************** 7. Recap process configuration and starting run
    
    if (!process_rank) {
        pidRoot = getpid();
        logStartup(&configParam, comm_size, pidRoot);
    }//fi 

    
    

    //*************** 8. Data Decomposition or Data generation
    
    double *local_buffer = NULL;    // input buffer to store the items to process (VALIDATION mode only)
    long local_bsize = 0;           // length of the data sample in charge of the process
    long totalInputSize = 0;        // total length of the input (=comm_size * local_bsize)

    double inputTime = 0;           // monitoring elapsed seconds to get input data (from file or self-generated)
    
    #ifdef VALIDATE
        
        //*** Step 8.1: (VALIDATION MODE)
        // each process reads the binary file and loads its correspondent data sample 
        // On return, totalInputSize contains the size of the _entire_ binary file in input
        
        if (configParam.processInputLen > 0) {
            //*** Limit the number of rows read from the file
            //fprintf(stderr, "Get only the first %ld items (if they are more than this)\n", configParam.processInputLen);
            inputTime -= MPI_Wtime();   
                local_buffer = (double *)readDataSlot(configParam.filename, comm_size, process_rank, &local_bsize, f_log, &totalInputSize, &configParam.processInputLen);
            inputTime += MPI_Wtime();
        }
        else 
        {   
            //*** Read all the file entries
            //fprintf(stderr, "Get all the binary file items\n");
            inputTime -= MPI_Wtime();   
                local_buffer = (double *)readDataSlot(configParam.filename, comm_size, process_rank, &local_bsize, f_log, &totalInputSize, NULL);
            inputTime += MPI_Wtime();
        }

        if (!local_buffer) { 
            std::cerr << "Error on getting input data\n" << std::endl; 
            if (f_log) {
                fprintf(f_log, "Error on generating the local dataset .. shutting down\n");
                fclose(f_log);
            }
            MPI_Finalize();
            exit(EDATA);
        } else {
            
            #ifdef LOG
                if (f_log) {
                    fprintf(f_log, "%.6f [seconds] elapsed to load data sample of %ld items\n\n", inputTime, local_bsize);
                }//fi 
            #endif

            //*** Step 8.2: (VALIDATION MODE) only process 0 performs the exact computations
            if (configParam.exactFlag && !process_rank) {

                long GLen = totalInputSize;
                double *Glocal_buffer = (double *) malloc(GLen * sizeof(double));
                
                if (!Glocal_buffer) {   
                    fprintf(stderr, "Error allocating the memory to perform exact computation\n");
                } else {
                    
                    FILE *fp_in = fopen(configParam.filename, "rb");
                    if (fp_in == NULL) {
                        fprintf(stderr, "Error opening %s\n", configParam.filename);
                    } else {
                        
                        fread(Glocal_buffer, sizeof(double)*GLen, 1, fp_in);
                        fclose(fp_in);

                        //*** Step 8.3: compute exact statistics 
                        
                        //fprintf(stderr, "Starting Exact Quantiles computation\n");
                        exactLen = 0;
                        exactStats = getExactQuantiles(Glocal_buffer, GLen, &exactLen, f_log);

                        if (Glocal_buffer) {
                            free(Glocal_buffer);
                        }
                    }//fi (fp_in)
                }//fi Glocal_buffer
            }//fi configParam.exactFlag       

        }//fi local_buffer
    

    #else //RUN
        //*** Step 8.4: (RUN MODE)
        // each process self-generates its data sample starting with a custom seed
        // each local shard has size given by the processInputLen
        // @note: generateData will fill the local_buffer with (items > NULLBOUND)
        // @note: instead of local_buffering items, they are generated on the fly

        local_bsize = configParam.processInputLen;
        totalInputSize = configParam.processInputLen * comm_size;
        
        //*** Step 8.5: local_buffering self-generated data items
        // inputTime -= MPI_Wtime();                    
        // local_buffer = (double *)generateData(&configParam, comm_size, process_rank);
        // inputTime += MPI_Wtime();

        //*** Step 8.6: Define the stream distribution params 
        std::uniform_real_distribution<double> udistribution(configParam.x, configParam.y);
        std::exponential_distribution<double> edistribution(configParam.x);
        std::normal_distribution<double> ndistribution(configParam.x, configParam.y);
        std::lognormal_distribution<double> ldistribution(configParam.x, configParam.y);
        
        std::default_random_engine generator;
        std::function<double()> randomizer;


        /*
        Beta distribution can be obtained as
        Z = X/(X+Y), where X is Gamma(x,1.0) and Y is Gamma(y, 1.0)
        so Z is Beta(x,y)
        */

        std::gamma_distribution<double> gammaX(configParam.x, 1.0);
        std::gamma_distribution<double> gammaY(configParam.y, 1.0);

        std::default_random_engine generateX;
        std::default_random_engine generateY;

        std::function<double()> randomizeX;
        std::function<double()> randomizeY;

        //generator.seed(std::chrono::system_clock::now().time_since_epoch().count() + process_rank);
        if (configParam.seed != -1){
            //fprintf(stdout, "%d) seed (%d)\n",process_rank, configParam.seed);

            generator.seed(configParam.seed + process_rank);

            generateX.seed(configParam.seed + process_rank);
            generateY.seed(configParam.seed + comm_size + process_rank);
        } else {
            generator.seed(std::chrono::system_clock::now().time_since_epoch().count() + process_rank);
            
            generateX.seed(std::chrono::system_clock::now().time_since_epoch().count() + process_rank);
            generateY.seed(std::chrono::system_clock::now().time_since_epoch().count() + comm_size + process_rank);
        }//fi


        switch(configParam.dtype) 
        {
            
            case UNIFORM:
                randomizer = std::bind(udistribution, generator);
                break;

            case EXPONENTIAL:
                randomizer = std::bind(edistribution, generator);
                break;

            case NORMAL:
                randomizer = std::bind(ndistribution, generator);
                break;

            case BETA:
                randomizeX = std::bind(gammaX, generateX);
                randomizeY = std::bind(gammaY, generateY);
                break;
            
            case LOGN:
                randomizer = std::bind(ldistribution, generator);
                break;

            default:
                randomizer = std::bind(ndistribution, generator);    
                break;
        }//switch
        
        //*** Step 8.6: instead of local_buffer[i] we have now a single var
        double current_item = 0.0; 
    #endif //RUN vs VALIDATE


    //*************** 9. Local Sketch construction

    /* 
     For a couple of sketches, the key estimation for each items is:
      key = ⎡log(x)⎤    iif  x   > 0
      key = ⎡log(-x)⎤   iif  x   < 0
      key = B*          iif |x|  < ε
    */


    double NCKEY = 0.75;        // DDOG: a value to distinguish not collapsed leftmost and rightmost keys
    double kNega = NCKEY;       // DDOG: bucket key of the collapsed bin in the negative sketch
    double kPosi = NCKEY;       // DDOG: bucket key of the collapsed bin in the positive sketch
    int trashable = 0;          // DDOG: items that falls directly into the trash bin



    int key = 0;                        // bucket key for items

    //***--> START MONITORING PARALLEL PROCESSING TIME

    double totTime = 0.0;               // time (seconds) spent between the BARRIERs
    
    double fillTime = 0.0;              // time (seconds) to fill the local sketch with the data sample
    double additionalTime = 0.0;        // time (seconds) spent for the additional collapses procedure before the Reduce 
    double mergeTime = 0.0;             // time (seconds) spent for the Reduce procedure

    double maxSketchTime = 0.0;             // max of times (seconds) spent for filling the local sketches
    double maxAdditionalTime = 0.0;         // max of times (seconds) spent for performing additional Collapses
    double maxMergeTime = 0.0;              // max of times (seconds) spent for the Reduce procedure

    
    
    MPI_Barrier(MPI_COMM_WORLD);
    totTime -= MPI_Wtime();
    


    //*** Step 9.1: Eeach process fills its local sketch
    fillTime = totTime;
    for(long i = 0; i < local_bsize; ++i) {
        
        #ifdef RUN   
            //*** Instead of filling the local_buffer, items are generated and processed without local_buffering in the Heap
            
            // current_item = randomizer();
            // while (current_item <= NULLBOUND) {
            //     current_item = randomizer();
            // }
            
            if (configParam.dtype != BETA) { 
                current_item = randomizer();
            
                while (current_item <= NULLBOUND) {
                    current_item = randomizer();
                }//wend

            } else {
                double X = randomizeX();
                double Y = randomizeY();
                current_item = X/(X+Y);
            }//fi dtype


            // keep minimum value for Low Bins Collapse in DDOG
            if (xmin > current_item) {
                xmin = current_item;
            }//fi 


            /* Add to positive sketch only
             * because we are sure that each item is >> 0
             * @note see generateData() */
            key = std::ceil(std::log10(current_item)/base); //std::ceil(std::log10(local_buffer[i])/base);
            
            #ifdef LowBins
                // if kPosi is no more equal to NCKEY, we have that the leftmost bin is collapsed
                // so, instead of adding and collapsing a new bin, add directly to the collapsed one
                if (process_collapses > 0 && key < (int)kPosi ) {
                    posibins += addKeyToSketch(posiSketch, (int)kPosi);
                    ++posipop;
                    ++trashable;
                    
                    #ifdef LOGCollapses
                        fprintf(f_log, "key %d in %d (%d)\n", key, (int)kPosi, trashable);
                    #endif

                } else {
                    posibins += addKeyToSketch(posiSketch, key);
                    ++posipop;
                    #ifdef LOGCollapses
                        fprintf(f_log, "New bin with key %d\n", key);
                    #endif
                }//fi
            
            #else //HighBins

                // if kPosi is no more equal to NCKEY, we have that the rightmost bin is collapsed
                // so, instead of adding and collapsing a new bin, add directly to the collapsed one
                if (process_collapses > 0 && key > (int)kPosi) {
                    posibins += addKeyToSketch(posiSketch, (int)kPosi);
                    ++posipop;

                    ++trashable;
                    #ifdef LOGCollapses
                        fprintf(f_log, "key %d in %d (%d)\n", key, (int)kPosi, trashable);
                    #endif
                } else {
                    posibins += addKeyToSketch(posiSketch, key);
                    ++posipop;
                    #ifdef LOGCollapses
                        fprintf(f_log, "New bin with key %d\n", key);
                    #endif
                }//fi
            #endif //BINS

        #else // VALIDATION
            
            /*
             * In VALIDATION MODE consider both positive and negative items: 
             * items are parsed from a file that was allocated in Heap's local_buffer
             * and we don't know neither their sign nor their magnitude */

            if (local_buffer[i] > NULLBOUND) 
            {   
                key = std::ceil(std::log10(local_buffer[i])/base);   
    
                #ifdef LowBins
                    // if we have collapsed yet, we have to check in which sketch is the collapsed bin:
                    // since we are adding to the positive sketch, we must verify that it already contains the leftmost collapsed bin 
                    if (process_collapses > 0 && ( (kPosi != NCKEY) && (key < (int)kPosi) ) )  {
                        posibins += addKeyToSketch(posiSketch, (int)kPosi);
                        ++posipop;

                        ++trashable;
                        #ifdef LOGCollapses
                            fprintf(f_log, "key %d in %d (%d)\n", key, (int)kPosi, trashable);
                        #endif

                    } else {
                        posibins += addKeyToSketch(posiSketch, key);
                        ++posipop;
                        #ifdef LOGCollapses
                            fprintf(f_log, "New bin with key %d\n", key);
                        #endif
                    }//fi

                #else //High bins
                    
                    // if we have collapsed yet, we have to check in which sketch is the collapsed bin:
                    // since we are adding to the positive sketch, we must verify that it already contains the rightmost collapsed bin
                    if (process_collapses > 0 && ( (kPosi != NCKEY) && (key > (int)kPosi) ) )  {
                        posibins += addKeyToSketch(posiSketch, (int)kPosi);
                        ++posipop;

                        ++trashable;
                        #ifdef LOGCollapses
                            fprintf(f_log, "key %d in %d (%d)\n", key, (int)kPosi, trashable);
                        #endif
                    
                    } else {
                        posibins += addKeyToSketch(posiSketch, key);
                        ++posipop;
                        #ifdef LOGCollapses
                            fprintf(f_log, "New bin with key %d\n", key);
                        #endif
                    }//fi
                #endif //Bins
            } 
            else if ( -NULLBOUND <= local_buffer[i] && local_buffer[i] <= NULLBOUND) 
            {  //  |local_buffer[i]| <= NULLBOUND         
            
                key = -MIN_KEY; //special bucket B* (for 0 and near-0 values): -(2^30)
                posibins += addKeyToSketch(posiSketch, key);
                ++posipop;
            } else {   
                // local_buffer[i] < -NULLBOUND
                //add to negative sketch
                key = std::ceil(std::log10(-local_buffer[i])/base);
            
                #ifdef LowBins
                    // since we are adding to the negative sketch, if it has been collapsed
                    // and the new item falls in a destined to collapse bin, add it directly to the leftmost bin
                    // that corresponds to the biggest key in the negative sketch
                    if (process_collapses > 0 && ( (kNega != NCKEY) && (key > (int)kNega) ) )  {
                        negabins += addKeyToSketch(negaSketch, (int)kNega);
                        ++negapop;

                        ++trashable;
                        #ifdef LOGCollapses
                            fprintf(f_log, "key %d in %d (%d)\n", key, (int)kNega, trashable);
                        #endif
                    } else {
                        negabins += addKeyToSketch(negaSketch, key);
                        ++negapop;
                        #ifdef LOGCollapses
                        fprintf(f_log, "New bin with key %d\n", key);
                        #endif
                    }//fi
                
                #else //High bins
                    // if the negative sketch has been collapsed
                    // and the new item falls in a destined to collapse bin, add it directly to the rightmost bin
                    // that corresponds to the smallest key in the negative sketch (near -0 values)
                    if (process_collapses > 0 && ( (kNega != NCKEY) && (key < (int)kNega) ) )  {
                        negabins += addKeyToSketch(negaSketch, (int)kNega);
                        ++negapop;

                        ++trashable;
                        #ifdef LOGCollapses
                            fprintf(f_log, "key %d in %d (%d)\n", key, (int)kNega, trashable);
                        #endif
                    } else {
                        negabins += addKeyToSketch(negaSketch, key);
                        ++negapop;
                        #ifdef LOGCollapses
                            fprintf(f_log, "New bin with key %d\n", key);
                        #endif
                    }//fi
                #endif //Bins            
            }//fi buffer[i]

        #endif //RUN vs VALIDATION

        //*** Step 9.2: Each process collapses its sketches for BOUND - if necessary
        process_collapses += OriginalPairCollapse(posiSketch, &posibins, negaSketch, &negabins, configParam.sketch_bound, &kNega, &kPosi, f_log);
        #ifdef LOGCollapses
            fprintf(stderr, "%ld) DDOG collapses: %d, trashed items: %d \n", i, process_collapses, trashable);
        #endif

        //snapshot of the collapse step
        monitorCollapses[m_idx].collapse = process_collapses;
        monitorCollapses[m_idx].items = posipop + negapop;
        //++m_idx;

    }//for fill local sketch
    fillTime += MPI_Wtime();
    
    // max time to fill the sketch among all the processes
    MPI_Allreduce(&fillTime, &maxSketchTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    
    #ifdef LOG
        if (f_log) {
            fprintf(f_log, "P%d) α: %.6f, collapses: %d, fillTime: %.6f\n", process_rank, p_alpha, process_collapses, fillTime);
            fprintf(f_log, "DDOG trashed items: %d \n", trashable);
        }
    #endif

    //*** Step 9.5: free the local_buffer of items (loaded from input file or self-generated)
    if (local_buffer) {
        free(local_buffer);
        local_buffer = NULL;
    }

    #ifndef NDEBUG
        //check population in the sketches
        {
            long psum = 0;
            for (auto it=posiSketch.begin(); it != posiSketch.end(); ++it) {
                psum += it->second;
            }//for
            
            long nsum = 0;
            for (auto it=negaSketch.begin(); it != negaSketch.end(); ++it) {
                nsum += it->second;
            }//for
            
            assert((psum+nsum) == (posipop+negapop));
        }
    #endif


    //*************** 10. CUSTOM MERGE of local sketches

    int collapsesA = 0;       //loop counter
    
    additionalTime -= MPI_Wtime();
        //*** Step 10.3: DDOG collapses before REDUCE
        // int DCollapses_bR = 0;          // DDOG collapses before reduce executed by P0
        // if (!process_rank) {
        //     DCollapses_bR = process_collapses;
        // }

        int DDog_MaxCollapses = 0;      // the max among all the collapses executed over all the MPI processes
        MPI_Allreduce(&process_collapses, &DDog_MaxCollapses, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        #ifdef LOG
            if (f_log) {
                fprintf(f_log, "\tP%d) local α=%.6f, local collapses=%d, ALL_MAX %d\n", process_rank, p_alpha, process_collapses, DDog_MaxCollapses);
            }
        #endif   
    additionalTime += MPI_Wtime(); 

    //get the max of the times for executing an ALLReduce() 
    // MPI_Allreduce(&additionalTime, &maxAdditionalTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 

    //*************** 11. Custom Packing of local sketches

    //*** Step 11.1: pack the pair of sketches into a local local_buffer
    Gbuffer_size = 0;  // the global variable is initialized here and used throughout the code during the Reduce phase
    
    mergeTime -= MPI_Wtime();

        char *localSketch = pack_SketchesPair(p_alpha, 0, configParam.sketch_bound, posiSketch, posibins, posipop, negaSketch, negabins, negapop, &Gbuffer_size); 


        if (!localSketch) {
            std::cout << "Could not allocate space for packing the local sketches\n";
            
            if (local_buffer) {
                free(local_buffer);
                local_buffer = NULL;
            }

            // close log file
            if (f_log) {
                fprintf(f_log,"Process computation ended!!!\n");
                fclose(f_log);
            }

            MPI_Finalize();
            exit(EXIT_FAILURE);
        }//fi localSketch



        //*** Step 11.2: allocate space for the global sketch
        char *finalSketch = (char *)malloc(Gbuffer_size * sizeof(char)); 
        if (!finalSketch) {
            std::cout << "Could not allocate space for packing the merged sketch\n";
            
            if (localSketch) {
                free(localSketch);
            }
            
            if (local_buffer) {
                free(local_buffer);
                local_buffer = NULL;
            }

            // close log file
            if (f_log) {
                fprintf(f_log,"Process computation ended!!!\n");
                fclose(f_log);
            }
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }//fi finalSketch
        memset(finalSketch, 0, Gbuffer_size * sizeof(char)); 
        


        //*** Step 11.3: Define custom data type

        MPI_Datatype Sketchtype;
        MPI_Type_contiguous(Gbuffer_size, MPI_BYTE, &Sketchtype);
        MPI_Type_commit(&Sketchtype);

        //*** Step 11.4: Define custom MPI_Reduce()
        int commute = 1;         // reduction is commutative
        MPI_Op custom_merge;

        MPI_Op_create(reduceOriginal, commute, &custom_merge);      //*** DDOG

            
            
        //*** Step 11.5: execute the MPI_Reduce()
        //    double reduceTime = 0.0;
        //reduceTime -= MPI_Wtime();
        MPI_Reduce(localSketch, finalSketch, 1, Sketchtype, custom_merge, root_process, MPI_COMM_WORLD);
        //reduceTime += MPI_Wtime();
        

    mergeTime += MPI_Wtime();    
    // save the time spent in performing the Reduce step
    MPI_Allreduce(&mergeTime, &maxMergeTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


    //***--> STOP MONITORING PARALLEL PROCESSING TIME
    MPI_Barrier(MPI_COMM_WORLD);
    totTime +=MPI_Wtime();
    

    double maxParallelTime = 0.0;
    MPI_Reduce(&totTime, &maxParallelTime, 1, MPI_DOUBLE, MPI_MAX, root_process, MPI_COMM_WORLD);

    

    //*** Step 11.8: Destroy custom merge and associated resources
    MPI_Op_free(&custom_merge);
    MPI_Type_free(&Sketchtype);

    // each process destroy its localSketch
    if (localSketch) {
        free(localSketch);
    }
    
    
    //all P's log to stdout the snapshot of the collapses executed
    fprintf(stdout, "P%d) collapses: %d, items: %ld (additionals: %d)\n", process_rank, monitorCollapses[m_idx].collapse, monitorCollapses[m_idx].items, collapsesA);
    


    //*************** 12. Queries on the Global Merged Sketch

    /*
        for DDOG only in the Low Bins collapse schema:
        keep the minimum item value in the processed stream
        to compute the alpha-accuracy for the collapsed bin at the end of the merge procedure
    */
    int local_i_min = std::ceil(std::log10(xmin)/base);
    int i_min = 0;
    MPI_Allreduce(&local_i_min, &i_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    //fprintf(stderr, "%d) %.6f: %d\n", process_rank, xmin, i_min);


    // only the process with rank 0 executes these analytics
    if (!process_rank) 
    {
        struct PackedPair *globalSketch = deserializeGlobalSketchPair(finalSketch, Gbuffer_size);

        #ifndef NDEBUG
            //check population in the final sketches
            {    
                fprintf(stderr, "Negabins: %d, npop: %ld\n", globalSketch->negabins, globalSketch->negapop);
                long dsumn = 0;
                if (globalSketch->negabins>0 && globalSketch->nega) {
                    for (int k=0; k<globalSketch->negabins;++k) {
                        dsumn += globalSketch->nega[k].count;
                    }
                }//fi nega

                fprintf(stderr, "Posibins: %d, ppop: %ld\n", globalSketch->posibins, globalSketch->posipop);
                long dsump = 0;
                if (globalSketch->posibins>0 && globalSketch->posi) {
                    for (int k=0; k<globalSketch->posibins; ++k) {
                        dsump += globalSketch->posi[k].count;
                    }
                }//fi posi

                std::cout<<"\nTotal population in final sketch " << dsump+dsumn << std::endl;
                std::cout<<"Total Input Size "<<totalInputSize << std::endl;
                assert((dsump+dsumn) == totalInputSize);
            }
        #endif
        
        //*** Throughput in million of update/s (upd(10^6)/s)
        
        // totTime is elapsed from the fill of the sketches until the end of the Reduces
        double throughput = totalInputSize/(totTime*1000000.0); //udp(10^6)/s
        
        // spentTime is the sum of the major phases of the algo: 
        // fill of the local sketches, additional collapses before Reduce, packing the local sketches, Reduce time
        //double spentTime = maxSketchTime + maxAdditionalTime + maxPackTime + maxMergeTime;
        double spentTime = maxSketchTime + maxAdditionalTime + maxMergeTime;
        double realThroughput = totalInputSize/(spentTime*1000000.0); //udp(10^6)/s

        //*** Compute quantiles (wrt the global array QQ)
        double estimates[nQ];             // where approximated quantiles are saved
        struct Bucket BKeys[nQ];          // where the bucket indexes corresponding to quantiles are saved 
        
        double quantilesTime = 0;
        quantilesTime -= MPI_Wtime();    
        for (int h = 0; h < nQ; ++h) {
            estimates[h] = PairQuantile(QQ[h], globalSketch, &BKeys[h].key, &BKeys[h].count);
        }//for h
        quantilesTime += MPI_Wtime();


        //*************** 13. Log results to csv file 
        //std::cout<< "Logging to stderr to redirect to csv file" << std::endl;        
        fprintf(stdout, "Distribution,param1:param2,StreamLen,#P,SketchBound,iAlpha,");
        fprintf(stdout, "fAlpha,Collapses,q0,q1,Ppop,Npop,Tpop,Pbins,Nbins,Tbins,");
        fprintf(stdout, "Throughput,RunningTime,QueryTime,realThroughput,SpentTime,MaxSketchTime,MaxAdditionalTime,MaxMergeTime");
        for(int c=0; c<nQ; ++c){
            fprintf(stdout, ",%.2f,Estimate,binKey,binCount",QQ[c]);
        }
        fprintf(stdout, "\n");



        #ifdef RUN
            char params[BSIZE];
            snprintf(params, BSIZE-1, "%.3f:%.3f", configParam.x, configParam.y);
            
            //Distribution, dist's Params, (ex: UNIFORM,1.000:1000000.000,)
            fprintf(stderr,"%s,%s,",configParam.dname, params);
        #else
            //Distribution, dist's Params, (ex: binFiles/uniformN.dat,-,) 
            fprintf(stderr,"%s,%s,",configParam.filename, "-");
        #endif
        
        //datalen, #processes, Sketch Bound, initial Alpha,
        fprintf(stderr,"%ld,%d,%d,%.6f,", totalInputSize, comm_size, configParam.sketch_bound, configParam.initial_alpha);

            
        // DDOG is q0-q1 accurate 
        double q0B = 0.0, q1B = 1.0;
        
        double gammaErr = 0.0, alphaErr = 0.0;

        // collapsed during the MPI_Reduce merge of pair of sketches
        int collapsesMerge = globalSketch->collapses;

        //force this check always
        if (1 || collapsesMerge || DDog_MaxCollapses) {

            int bound = configParam.sketch_bound;
            int pb = globalSketch->posibins;
            int nb = globalSketch->negabins;
            long ppop = globalSketch->posipop;
            long npop = globalSketch->negapop;
        
            #ifdef LowBins

                long tCount = 0;
                int tFirst = 0;

                if (!npop) {
                    // we must only deal with the positive sketch  

                    // the first bucket is "trash"
                    tCount = globalSketch->posi[0].count;
                    tFirst = globalSketch->posi[0].key; //this is the key of the first bin 

                    if (globalSketch->posi[0].key == -MIN_KEY) {
                        tCount += globalSketch->posi[1].count;
                        tFirst = globalSketch->posi[1].key;
                    }

                    // so the min accurate quantile q0 is:
                    q0B = tCount/(1.0*ppop);

                    // to compute the error for the first bin (that is collapsed in this scenario)
                    gammaErr = pow(gamma, tFirst - i_min + 1);      // this is the gamma representing the range of values in the first collapsed bin
                    alphaErr = (gammaErr - 1.0)/(gammaErr + 1.0);   // this is the corresponding alpha value
                    //fprintf(stdout, "%d, %d, %.6f\n", tFirst, i_min, gammaErr);
                    
                } else {

                    //we have also the negative sketch
                    if ((bound - pb) > 1) {
                        //only the first bucket in the negative sketch is the trash
                        tCount = globalSketch->nega[0].count;

                        // so the min accurate quantile q0 is:
                        q0B = tCount/(1.0*(ppop+npop));
                    } else {

                        /** @note
                         * If we have only a single 1 bin in the negative sketch, 
                         * for sure it's a trash bucket
                         * but ** at this time ** we can't say 
                         * if also the first bin in the positive final sketch 
                         * is trash bin at the end of the MPI_Reduce 
                         * (we know this in the sequential phase thanks to the fencing keys instead)
                        */

                        // so, let say that in this scenario, the first positive bin is trash
                        tCount = globalSketch->nega[0].count;
                        tCount += globalSketch->posi[0].count;
                        if (globalSketch->posi[0].key == -MIN_KEY) {
                            tCount += globalSketch->posi[1].count;
                        }

                        // so the min accurate quantile q0 is:
                        q0B = tCount/(1.0*(ppop+npop));
                    }//fi sizen
                }//fi npop

            #else //HighBins

                long tCount = 0;

                if (ppop) {

                    if ((pb > 2) || (pb==2 && globalSketch->posi[0].key != -MIN_KEY) ) {
                        tCount = globalSketch->posi[pb-1].count;
                        q1B = 1.0 - (tCount/(1.0*(ppop+npop)));
                    } else {
                        //pb == 2 and == -MIN_KEY
                        //pb == 1
                        
                        tCount = globalSketch->posi[0].count; //B* count
                        if (pb == 2) {
                            tCount += globalSketch->posi[1].count;
                        }

                        //and add the last bucket of the negative sketch
                        //@note same as in the LowBins indeterminable scenario 
                        tCount += globalSketch->nega[nb-1].count;
                        q1B = 1.0 - (tCount/(1.0*(ppop+npop)));     
                    }//fi pb>2
                
                } else {
                    //only negative sketch
                    tCount = globalSketch->nega[nb-1].count;
                    q1B = 1.0 - (tCount/(1.0*npop));     
                }//fi ppop
            #endif
        }//fi q0-q1 check
        
        // keep the max collapses among all the processes
        // finalAlpha on Collapsed bin, Tot Collapses, q0-bound, q1-bound
        fprintf(stderr,"%.6f,%d,%.6f,%.6f,", alphaErr, DDog_MaxCollapses, q0B, q1B);


        //Ppop, Npop, Tpop,
        fprintf(stderr, "%ld,%ld,%ld,", globalSketch->posipop, globalSketch->negapop, globalSketch->posipop + globalSketch->negapop);
            
        //Pbins,Nbins,Tbins,
        fprintf(stderr, "%d,%d,%d,", globalSketch->posibins, globalSketch->negabins, globalSketch->posibins + globalSketch->negabins);
            
        //Throughput, total running time [s], P0 QueryTime, 
        fprintf(stderr, "%.16f,%.16f,%.16f,", throughput, totTime, quantilesTime);
        

        //effective time spent in computations: fillTime, collapseTime, packTime, reduceTime [s], 
        //fprintf(stderr, "%.16f,%.16f,%.16f,%.16f,%.16f,%.16f", realThroughput, spentTime, maxSketchTime, maxAdditionalTime, maxPackTime, maxMergeTime);
        fprintf(stderr, "%.16f,%.16f,%.16f,%.16f,%.16f", realThroughput, spentTime, maxSketchTime, maxAdditionalTime, maxMergeTime);
        

        
        //q0 approximate quantiles
        #ifdef VALIDATE   
            if (exactStats) {
                for (int h = 0; h < nQ; ++h ) {
                    //double rel_err = std::abs((estimates[h] - exactStats[h])/exactStats[h]);
                    double rel_err = std::abs(estimates[h] - exactStats[h]);
                    double norm = std::abs(exactStats[h]);
                    if (norm > 0) {
                        rel_err /= norm;
                    }
                    fprintf(stderr, ",%.2f,%.6f,%.6f,%.6f", QQ[h], estimates[h], exactStats[h], rel_err);
                    //fprintf(stderr, ",%.2f,%.6f,%.6f,%d,%ld", QQ[h], estimates[h], rel_err, BKeys[h].key, BKeys[h].count);
                }//for h
                fprintf(stderr, "\n");
                free(exactStats);
                exactStats = NULL;
            }//fi exactStats
        #else //RUN
            for (int h = 0; h < nQ; ++h ) {
                fprintf(stderr, ",%.2f,%.6f,%d,%ld", QQ[h], estimates[h], BKeys[h].key, BKeys[h].count);
                //fprintf(stderr, ",%.2f,%.6f", QQ[h], estimates[h]);
            }//for h
            fprintf(stderr, "\n");
        #endif

        if (globalSketch) {

            if (globalSketch->posi){
                free(globalSketch->posi);
            }

            if (globalSketch->nega){
                free(globalSketch->nega);
            }

            free(globalSketch);
        }//fi globalSketch

    }//fi P0

    //************************************** 14. On exit
    
    if (finalSketch) {
        free(finalSketch);
    }

    // close log file
    if (f_log) {
        fprintf(f_log,"Process computation ended!!!\n");
        fclose(f_log);
    }

    // shutting down MPI    
    MPI_Finalize();

    // if (!process_rank) {
    //     #ifdef VALIDATE
    //         std::cout << "\nVALIDATION test ended!\n"<< std::endl;
    //     #else
    //         std::cout << "\nRUN test ended!\n"<< std::endl;
    //     #endif
    // }

    return 0;
}
