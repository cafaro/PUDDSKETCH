/**
 * @file Utility.cc
 * @author CM
*/

#include "Utility.h"

#include <unistd.h>
#include <string.h>
#include <stdlib.h>

extern char VERSION[];
extern double NULLBOUND;


//************************************************************ Command Line Management

void usage(char *msg) {
    
    // {[-f path-to-file]} | {[-d input-distribution] [ -n input_len ] [-x param1] [-y param2]}
    // [ -a initial-alpha ] [-b max_sketch_bound] 

   	std::cerr << "\n\tUsage: " << msg;
    std::cerr << " { ([-f path-to-file]) || ([-d input-distribution] [-x distribution param1] [-y distribution param2] [ -n input_len ]) } ";
    std::cerr << "[ -s initial-seed ] ";
    std::cerr << "[ -a initial-alpha ] [-b max_sketch_bound] -e [if specified, exact quantiles comparison is done]" << std::endl;
    
    std::cerr << "\n -d can be: \n";
    std::cerr << " : uniform distribution [1], so -x is for param 'a' and -y for param 'b' of Uniform [a,b)\n";
    std::cerr << " : exponential distribution [2], so -x is for 'λ' of Exponential\n";
    std::cerr << " : normal distribution [3], so -x is for param 'mean' and -y for param 'stddev' of Normal\n";
    std::cerr << " : beta distribution [4], so -x and -y are the shape params of Beta\n";
    std::cerr << " : lognormal distribution [5], so -x and -y are the shape params of LogNormal\n";

    std::cerr << "\n -n input_len must be expressed in term of Millions: \n";
    std::cerr << " : 0.5\t= (0.5*1M)\t= 500.000, \n";
    std::cerr << " : 1\t= (1*1M)\t= 1.000.000, \n";
    std::cerr << " : 2\t= (2*1M)\t= 2.000.000, ...\n"; 
}



/**
 * @brief Maps enumeration to corresponding string
 */
void getDistName(int dtype, char *dname, int len) {

    switch (dtype){
        case UNIFORM:
            memcpy(dname, "UNIFORM", sizeof(char)*len);
            break;

        case EXPONENTIAL:
            memcpy(dname, "EXPONENTIAL", sizeof(char)*len);
            break;

        case NORMAL:
            memcpy(dname, "NORMAL", sizeof(char)*len);
            break;
        
        case BETA:
            memcpy(dname, "BETA", sizeof(char)*len);
            break;

        case LOGN:
            memcpy(dname, "LOGNORMAL", sizeof(char)*len);
            break;
        
        default:
            memcpy(dname, "UNKNOWN", sizeof(char)*len);
    }//switch
}


void debugParams(Params *p) {

    if (p) {
        fprintf(stderr, "sketch bound : %d\n", p->sketch_bound);            
        fprintf(stderr, "alpha: %.6f\n", p->initial_alpha);           
        fprintf(stderr, "len of local data for current process: %ld\n", p->processInputLen);
        if (p->fileflag){
            fprintf(stderr, "Filename: %s\n", p->filename);
        }
        else{
            fprintf(stderr, "distribution type choosen by user: %d (%s)\n", p->dtype, p->dname);
            fprintf(stderr, "distribution param 1: %6f\n", p->x ); 
            fprintf(stderr, "distribution param 2: %6f\n", p->y ); 
            //fprintf(stderr, "number of items generated: %ld\n", p->processInputLen);
        }
        if (p->exactFlag){
            fprintf(stderr, "exact quantiles processing is enabled\n");
        }
    }//fi p
}



int check_params(int argc, char *argv[], Params *p) {
    
    int res = 0;
    
    if (!p) {
        std::cout << "Failed to pass Params"<<std::endl;
        exit(EPARAMS);
    }

    bool fbound = false;            // sketch bound is mandatory
    bool falpha = false;            // alpha is mandatory

    p->sketch_bound = 0;            // bound on the number of bins in the sketch
    p->initial_alpha = 0;           // initial α for each DDSketch

    // version 8.3: put a limit on the number of rows gathered from a long input file
    p->processInputLen = 0;         // len of local data for current process
    //p->processInputLen = INPUTLEN;

    p->seed = -1;                    // invalid initial seed

    p->dtype = UNKNOWN;             // distribution type choosen by user
    p->x = 0;                       // distribution param 1
    p->y = 0;                       // distribution param 2
    p->fileflag = false;            // true if input is from file, false if self-generated by the process 

    p->exactFlag = false;           // by default, exact quantiles computation is disabled

    bool fdist = false;             // fdist is true if input distribution is given
    int etype = -1;                 // to map distribution type

    int c=0;    
    while ( (c = getopt(argc, argv, "f:b:a:d:n:x:y:s:e")) != -1) {
    
        // {[-f path-to-file]} | {[-d input-distribution] [ -n input_len ] [-x param1] [-y param2]}
        // [ -a initial-alpha ] [-b max_sketch_bound] [-s seed]
        
        switch (c) {  
        
            case 'f':
                // input filepath
                if (strlen(optarg) <= BSIZE) {
                    snprintf(p->filename, BSIZE-1,"%s", optarg);
                    p->fileflag = true;
                }
                break;
            
            case 'b':
                // sketch memory bound
                p->sketch_bound = atoi(optarg);
                fbound = true;
                break;

            case 'a':
                // initial alpha for DDSketch
                p->initial_alpha = strtod(optarg, NULL); 
                falpha = true;
                break;

            case 'd':
                // distribution generator
                etype = atoi(optarg); 
                fdist = true;
                break;        

            case 'x':
                //lambda param for EXPONENTIAL
                //a param for UNIFORM
                //mu param for NORMAL
                //alpha for BETA
                p->x = strtod(optarg, NULL);
                break;
            
            case 'y':
                //b param for UNIFORM
                //stddev param for NORMAL
                //beta for BETA
                p->y = strtod(optarg, NULL);
                break;

            case 'n':
                // # of items to generate for each process, in millions (10^6)
                // version 8.3: put a limit on the number of rows from a file
                p->processInputLen = INPUTLEN * strtod(optarg, NULL); 
                //p->processInputLen *= strtod(optarg, NULL); //strtol(optarg, NULL, 10);
                break;

            case 'e':
                //enable exact quantiles computation
                p->exactFlag = true;
                break;

            case 's':
                // initial seed
                p->seed = atoi(optarg);
                break;

            default:
                fprintf(stderr, "Incorrect parameter specification for: '%c'\n", optopt);
                break;
        }// switch
    }//wend

    if (fdist) {
        fdist = false;

        switch (etype) {

            case UNIFORM:
                fdist = true;
                p->dtype = UNIFORM;

                if ( (p->x==0 && p->y==0) || (p->y < p->x) ) {
                    std::cout<<"Error on setting the range [a,b) for Uniform distribution\n" << std::endl;
                    return EPARAMS;
                }
                break;

            case EXPONENTIAL:
                fdist = true;
                p->dtype = EXPONENTIAL;
                if (p->x==0){
                    std::cout<<"Error on setting the lambda (λ) for Exponential distribution\n" << std::endl;
                    return EPARAMS;
                }
                break;

            case NORMAL:
                fdist = true;
                p->dtype = NORMAL;
                if (p->x==0 && p->y==0){
                    std::cout<<"Error on setting mean (μ) and stddev (σ) for Normal distribution\n" <<std::endl;
                    return EPARAMS;
                }
                break;

            case BETA:
                fdist = true;
                p->dtype = BETA;
                if (p->x==0.0 && p->y==0.0){
                    std::cout<<"Error on setting alpha (α) and beta (β) for Beta distribution\n" <<std::endl;
                    return EPARAMS;
                }
                break;
            

            case LOGN:
                fdist = true;
                p->dtype = LOGN;
                if (p->x==0 && p->y==0){
                    std::cout<<"Error on setting mean (μ) and stddev (σ) for LogNormal distribution\n" <<std::endl;
                    return EPARAMS;
                }
                break;

            default:
                p->dtype = UNKNOWN;
                break;
        }//switch

        getDistName(p->dtype, p->dname, 32);     // get distribution name
    }//fi fdist

    if (!fdist && !p->fileflag) {
        std::cout << "You must specify either the filename or the kind of distribution to process!\n"<<std::endl;
        return EPARAMS;
    }

    if (!fbound) {
        std::cout << "You must specify a bound on the sketch size!\n"<<std::endl;
        return EPARAMS;
    }

    if (!falpha) {
        //@note: set to default alpha ???
        std::cout << "You don't specify a value for Alpha, setting to default: "<< ALPHA << std::endl;
        p->initial_alpha = ALPHA; 
    }
    
return res;
}



void logStartup(Params *p, int nprocs, int pid) {
    //std::cout << "\n\tStarting "<< VERSION << std::endl;

    #ifdef UDD
        std::cout << "\n\tStarting "<< VERSION << " for Uniform-DDSketch"<< std::endl;
    #else
        std::cout << "\n\tStarting "<< VERSION << " for DataDog's original DDSketch"<< std::endl;
    #endif

    std::cout << "\twith parameters:\n";
    std::cout << "\t : MPI initialized on " << nprocs << " process(es)" << std::endl;
    std::cout << "\t : Root process's pid " << pid <<std::endl;    
    std::cout << "\t : Initial α: "<< p->initial_alpha << std::endl;      // alpha << std::endl;
    std::cout << "\t : Sketch bound: "<< p->sketch_bound << std::endl;    // bound << std::endl;
    std::cout << "\t : Random Seed: "<< p->seed << std::endl;    
    
    if (p->fileflag)
    {
        std::cout << "\t : From input file: "<< p->filename; // << std::endl;     
        if (p->processInputLen>0) {
            std::cout << " read only "<< p->processInputLen << " rows "<< std::endl;      
        } else {
            std::cout << " read all lines "<< std::endl;      
        }

        if (p->exactFlag){
            std::cout << "\t : Enabled exact quantiles computation "<< std::endl;     
        }
    } 
    else
    {
        std::cout << "\t : Distribution type: "<< p->dname <<" with param: "<< p->x;
        if (p->dtype != EXPONENTIAL){
            std::cout << ","<<p->y;
        } 
        std::cout << ", # items: "<< p->processInputLen << " per process"<<std::endl;  
    }//fi distr
    
    #ifdef VALIDATE
        std::cout << "\tVALIDATION mode"<< std::endl;
    #else
        std::cout << "\tRUN mode"<< std::endl;
    #endif
    
    #ifdef LOG
        std::cout << "...Enabled log to file for each process\n"<< std::endl;
    #endif
}



//************************************************************ Process Logging

FILE *logProcess(Params *p) {

    int pid = getpid();

    char logfile[BSIZE];

    if (p->fileflag){
   
        std::string name = p->filename;
        std::size_t from_p = name.find_last_of("/"); 
        std::string distr = name.substr(from_p+1, 4); 
        snprintf(logfile, BSIZE-1, "./PLogs/%d-%s_%.6f_%d_%ld.log", pid, distr.c_str(), p->initial_alpha, p->sketch_bound, p->processInputLen/1000);
    } else {
   
        snprintf(logfile, BSIZE-1, "./PLogs/%d-%s_%.6f_%d_%ld.log", pid, p->dname, p->initial_alpha, p->sketch_bound, p->processInputLen/1000);
    }//fi filename
    

    FILE *log = fopen(logfile, "w");
    if (log == NULL) {
        std::cerr<< "Process with "<<pid<<": cannot open " << logfile << std::endl;
        exit(ELOGFILE);
    } 
    else {
        fprintf(log, "Log for process %d\n", pid);
        #ifdef UDD
        fprintf(log, "UDD Sketch\t");
        #else
        fprintf(log, "DDog Sketch\t");
        #endif
        
        #ifdef RUN
        fprintf(log, "RUN Mode\n");
        #else
        fprintf(log, "VALIDATION Mode\n");
        #endif
    }

    return log;
}




/**
 * @brief to save buffer of doubles to binary file
 */
FILE *logBinary(Params *p) {

    char logfile[BSIZE];
    int pid = getpid();
    if (p->fileflag){
        std::string name = p->filename;
        std::size_t from_p = name.find_last_of("/"); 
        std::string distr = name.substr(from_p+1, 4); 
        snprintf(logfile, BSIZE-1, "./BinaryLogs/%d-%s-%.6f-%d-%ld.log", pid, distr.c_str(),p->initial_alpha, p->sketch_bound, p->processInputLen/1000);
    } else {
        snprintf(logfile, BSIZE-1, "./BinaryLogs/%d-%s-%.6f-%d-%ld.log", pid, p->dname,p->initial_alpha, p->sketch_bound, p->processInputLen/1000);
    }//fi filename
    
    FILE *log = fopen(logfile, "wb");
    if (log == NULL){
        std::cerr<< "Process with "<<pid<<": cannot open " << logfile << std::endl;
        exit(ELOGFILE);
    } 
    return log;
}

//************************************************************ MPI Initialization 

int warming_up_MPI(int *argc, char **argv[], int *comm_size, int *process_rank) {
    
    int ecode = MPI_SUCCESS;    // MPI error code returned by the APIs

    ecode = MPI_Init(argc, argv);
    if (ecode != MPI_SUCCESS){
        std::cerr << "MPI Error code " << ecode << std::endl;
        MPI_Finalize();
        return ecode;
    }

    ecode = MPI_Comm_size(MPI_COMM_WORLD, comm_size);
    if (ecode != MPI_SUCCESS){
        std::cerr << "MPI Error code " << ecode << std::endl;
        MPI_Finalize();
        return ecode;
    }

    ecode = MPI_Comm_rank(MPI_COMM_WORLD, process_rank);
    if (ecode != MPI_SUCCESS){
        std::cerr << "MPI Error code " << ecode << std::endl;
        MPI_Finalize();
        return ecode;
    }

return ecode;
}



//************************************************************ Process Termination

void endExecution(int rank, FILE *logfd, int ecode) {
    
    // shutting down MPI environment
    std::cerr << "shutting down MPI environment for P" << rank << std::endl;
    MPI_Finalize();
    
    //close log file
    if(logfd) {
        fclose(logfd);
    }
    
    exit(ecode);
}



//************************************************************ Data Decomposition

/**
 * @brief This function reads a binary file into a data slot
 * based on the file size in byte and on the process rank
 */
double *readDataSlot(char *binfilename, int comm_size, int process_rank, long *data_size, FILE *log, long *nitems, long *maxrows) {    
    
    FILE *fp_in = fopen(binfilename, "rb");         // binary file
    if (fp_in == NULL) {
        std::cerr << "Error opening " << binfilename << std::endl;
        return NULL;
    }

    fseek(fp_in, 0, SEEK_END);
    long fileSize = ftell(fp_in);                   // file size in bytes
    long ndouble = fileSize/sizeof(double);         // # of items (double type)

    // version 8.3: specify a limit on the number of rows
    if (maxrows && *maxrows < ndouble) {
        ndouble = *maxrows;
    }

    /* 
    Partitioning of the data set:
    computing ⎡i*N/p⎤, where:
    - N is the # of items 
    - p the # of processes
    - i the process's rank
    */
    long start_pos = std::ceil((double)(process_rank * ndouble)/comm_size);          // ⎡i*N/p⎤
    long end_pos = std::ceil((double)((process_rank + 1) * ndouble)/comm_size) - 1;  // ⎡(i+1)*N/p⎤ - 1
    
    *data_size = (end_pos - start_pos + 1);                                          // local shard size
    
    double *buffer = (double *)malloc(sizeof(double)*(*data_size));

    long offset = start_pos * sizeof(double);                                        //offset in bytes
    fseek(fp_in, offset, SEEK_SET);
    fread(buffer, sizeof(double)*(*data_size), 1, fp_in);                            // read *data_size doubles to local buffer
    fclose(fp_in);
    
    #ifdef VERBOSE
        if (log) {
            fprintf(log, "Binary file size %ld (in bytes)\t # items: %ld\n", fileSize, ndouble);
            fprintf(log, "Partitioning the shard among %d processes\n", comm_size);
            fprintf(log, "P%d reads from item %ld-th to item %ld-th (%ld items)\n", process_rank, start_pos, end_pos, *data_size);
        }//fi log
    #endif

    *nitems = ndouble; //return the number of items in the processed input file
    return buffer;
}



void debugLocalShard(int rank, double *buffer, long size, FILE *fd){
    fprintf(fd, "\nProcess %d local shard\n", rank);
    for(long i = 0; i < size; ++i){
        fprintf(fd, "%.32f\n", buffer[i]);
    }//for
    fprintf(fd, "\n");
}



/**
 * @brief This function generate a data slot from a given distribution
 * @note ATTENTION: we make sure that items are all >> 0
 * 
 *  
 * Random number distribution that produces floating-point values according to:
 * 
 * 1) Normal distribution
    Distribution parameters: mean (μ) and stddev (σ), (set on construction).
    std::normal_distribution<double> distribution((μ),(σ));

    This distribution produces random numbers around the distribution mean (μ) with a specific standard deviation (σ).
    The normal distribution is a common distribution used for many kind of processes, 
    since it is the distribution that the aggregation of a large number of independent random variables approximates to, 
    when all follow the same distribution (no matter which distribution).
        
 * 2) Uniform real distribution
    Distribution parameters: a and b - range [a,b) - (set on construction).
    std::uniform_real_distribution<double> distribution((a),(b));
        
    This distribution (also know as rectangular distribution) produces random numbers in a range [a,b) 
    where all intervals of the same length within it are equally probable.


 * 3) Exponential distribution 
    Distribution parameter: lambda (λ) - (set on construction).
    std::exponential_distribution<double> distribution((λ));
        
    This distribution produces random numbers where each value represents the interval between two random events 
    that are independent but statistically defined by a constant average rate of occurrence (its lambda, λ).
    Its analogous discrete distribution is the geometric_distribution.

 * 4) chi_squared distribution 
    Distribution parameter: degrees of freedom (n) - (set on construction).
    td::chi_squared_distribution<double> distribution((n));
    
    This distribution produces random numbers as if the square of n independent standard normal 
    random variables (Normal with μ=0.0 and σ=1.0) were aggregated, 
    where n is the distribution parameter, known as degrees of freedom.
 *
 */
double *generateData(Params *p, int comm_size, int process_rank) {

    if (!p){
        std::cout<<"Error passing configuration parameters"<<std::endl;
        return NULL;
    }    

    //*** Define the stream distribution params 

    //std::uniform_real_distribution<double> distribution((a),(b));
    std::uniform_real_distribution<double> udistribution(p->x, p->y);

    //std::exponential_distribution<double> distribution((λ));
    std::exponential_distribution<double> edistribution(p->x);
    
    //std::normal_distribution<double> distribution((μ),(σ));
    std::normal_distribution<double> ndistribution(p->x, p->y);
    
    std::default_random_engine generator;
    std::function<double()> randomizer;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count() + process_rank );

    switch(p->dtype){

        case UNIFORM:
            randomizer = std::bind(udistribution, generator);
            break;

        case EXPONENTIAL:
            randomizer = std::bind(edistribution, generator);
            break;

        case NORMAL:
            randomizer = std::bind(ndistribution, generator);
            break;
        
        default:
            randomizer = std::bind(ndistribution, generator);    
            break;
    }//switch

    //*** fill the buffer
    double *buffer = (double *)malloc(sizeof(double)*(p->processInputLen));
    int i = 0;
    while (i < p->processInputLen) {
        buffer[i] = randomizer();
        // @note make sure buffer[i] >>0 
        if (buffer[i] > NULLBOUND) {
           ++i;
        } //fi 
    }//wend

    return buffer;
}



//************************************************************ DEBUGS

void debugPack(struct PackedPair *p, FILE *log) {
    
    if (log == NULL) {
        log = stderr;
    }

    fprintf(log, "DebugPack--->\n");
    fprintf(log, "Alpha: %.6f, Bound: %d, Collapses: %d\n", p->alpha, p->bound, p->collapses);
    
    fprintf(log, "PosiBins: %d, PosiPop: %ld\n", p->posibins, p->posipop);
    
    if (p->posibins>0 && p->posi != NULL){
        for(int i = 0; i < p->posibins; ++i) {
            fprintf(log, "\tkey: %d, count: %ld\n", p->posi[i].key, p->posi[i].count);
        }
    }

    fprintf(log, "NegaBins: %d, NegaPop: %ld\n", p->negabins, p->negapop);
    
    if (p->negabins>0 && p->nega != NULL){
        for(int i = 0; i < p->negabins; ++i) {
            fprintf(log, "\tkey: %d, count: %ld\n",p->nega[i].key, p->nega[i].count);
        }
    }
    fprintf(log, "***\n");
}


//************************************************************ SKETCHT Utilities

void initSketchT(SketchT *S, double alpha, int bound) {

    /*
     * In RUN MODE we only want to deal with positive items
     * In VALIDATION MODE we will deal also with negative items
     * 
     * We can 
     * - handle a pair of sketches 
     * - or have only a positive one, managing the values generated 
     *      in RUN Mode to be sure items are ALL positives
     * 
     * \begin{verbatim}
     * 
            std::map<int,int> posiSketch;               
            int posipop = 0;                            
            int posibins = 0;                           
            #ifdef VALIDATE
                std::map<int,int> negaSketch;              
                int negapop = 0;                         
                int negabins = 0;                          
            #endif
     * 
     * \end{verbatim}
     * 
     */

    S->posipop = 0;
    S->posibins = 0;

    S->negapop = 0;
    S->negabins = 0;

    S->process_collapses = 0;
    S->mbound = bound;

    S->p_alpha = alpha;
    S->gamma = (1 + S->p_alpha)/(1 - S->p_alpha);
    S->base = std::log10(S->gamma);
    S->NULLBOUND = pow(S->gamma, -MIN_KEY);
}
