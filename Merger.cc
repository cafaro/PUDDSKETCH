/**
 * @file Merger.cc 
 * @author CM
*/

#include "Merger.h"

extern int Gbuffer_size;

//****************************************************** 

/**
 * This procedure reverses keys ordering in the negative sketch to simplify quantile queries 
 */
struct PackedPair *deserializeGlobalSketchPair(void *buffer, int size) {
    
    int pos = 0;
    
    struct PackedPair *res = (struct PackedPair *) malloc( sizeof(struct PackedPair) );
    memset(res, 0, sizeof(struct PackedPair));

    MPI_Unpack(buffer, size, &pos, &(res->alpha), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer, size, &pos, &(res->collapses), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, size, &pos, &(res->bound), 1, MPI_INT, MPI_COMM_WORLD);
    
    MPI_Unpack(buffer, size, &pos, &(res->posibins), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, size, &pos, &(res->posipop), 1, MPI_LONG, MPI_COMM_WORLD);
    //MPI_Unpack(buffer, size, &pos, &(res->posipop), 1, MPI_INT, MPI_COMM_WORLD);


    #ifndef NDEBUG
        fprintf(stderr, "deserializeGlobalSketchPair => res->posibins: %d, posipop: %ld\n", res->posibins, res->posipop);
        long psum = 0;
    #endif

    
    if (res->posibins > 0){
        
        res->posi = (struct Bucket *)malloc(sizeof(struct Bucket)*(res->posibins));
        
        //MPI_Unpack(buffer, size, &pos, res->posi, 2*(res->posibins), MPI_INT, MPI_COMM_WORLD);
        for (int b = 0; b < res->posibins; ++b) {
            MPI_Unpack(buffer, size, &pos, &(res->posi[b].key), 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, size, &pos, &(res->posi[b].count), 1, MPI_LONG, MPI_COMM_WORLD);

            #ifndef NDEBUG
                psum += res->posi[b].count;
            #endif
        }
    }//fi posibins
    
    #ifndef NDEBUG
        fprintf(stdout, "posipop: %ld, sum: %ld\n", res->posipop, psum);
        assert(psum == res->posipop);
    #endif

    MPI_Unpack(buffer, size, &pos, &(res->negabins), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, size, &pos, &(res->negapop), 1, MPI_LONG, MPI_COMM_WORLD);
    //MPI_Unpack(buffer, size, &pos, &(res->negapop), 1, MPI_INT, MPI_COMM_WORLD);

    #ifndef NDEBUG
        fprintf(stdout, "deserializeGlobalSketchPair => res->negabins: %d, negapop: %ld\n", res->negabins, res->negapop);
        long nsum = 0, nn = 0;
    #endif

    if (res->negabins > 0) {
        
        // as discussed, for quantile queries we have to reverse keys order in the negative sketch         
        
        res->nega = (struct Bucket *) malloc( sizeof(struct Bucket)*res->negabins);
        //MPI_Unpack(buffer, size, &pos, res->nega, 2*(res->negabins), MPI_INT, MPI_COMM_WORLD);
        
        //struct Bucket negaValues[res->negabins];
        struct Bucket *negaValues = (struct Bucket *) malloc( sizeof(struct Bucket)*res->negabins);
        if (!negaValues){
            fprintf(stderr, "unable to allocate buffer in heap\n");
        }

        //MPI_Unpack(buffer, size, &pos, negaValues, 2*(res->negabins), MPI_INT, MPI_COMM_WORLD);
        for (int b = 0; b < res->negabins; ++b) {
            MPI_Unpack(buffer, size, &pos, &(negaValues[b].key), 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, size, &pos, &(negaValues[b].count), 1, MPI_LONG, MPI_COMM_WORLD);

            #ifndef NDEBUG
            nn += negaValues[b].count;
            #endif
        }

        for (int i = 0; i < res->negabins; ++i) { //for (int s = res->negabins-1; s>=0; --s){
            int s = (res->negabins-1)-i;
            res->nega[s].key = negaValues[i].key;       //res->nega[res->negabins-1-i].key = negaValues[i].key;
            res->nega[s].count = negaValues[i].count;   //res->nega[res->negabins-1-i].count = negaValues[i].count;
            
            #ifndef NDEBUG
            nsum += res->nega[s].count;
            #endif
        }

        if (negaValues){
            free(negaValues);
        }
        
    }//fi negabins

    #ifndef NDEBUG
        fprintf(stderr, "negapop: %ld, sum: %ld, nn: %ld\n", res->negapop, nsum, nn);
        assert(nsum == res->negapop);
    #endif
    
    return res;
}



//****************************************************** DEBUG

void debugGlobalSketchPair(struct PackedPair *p, FILE *log) {
    
    if (log == NULL){
        log = stderr;
    }

    double gamma = (1+p->alpha)/(1-p->alpha);
    fprintf(log, " ---> Debugging Global Sketch --->\n");
    fprintf(log, "Alpha: %.6f, Collapses: %d, Bound: %d\n", p->alpha, p->collapses, p->bound);
    fprintf(log, "PosiBins: %d, PosiPop: %ld\n", p->posibins, p->posipop);
    fprintf(log, "NegaBins: %d, NegaPop: %ld\n", p->negabins, p->negapop);
    
    if (p->negabins > 0 && p->nega != NULL) {
        
        for(int i = 0; i < p->negabins; ++i) {
            fprintf(log, "\tkey: %d,\tcount: %ld,\t",p->nega[i].key, p->nega[i].count);
            double estimate = -(2 * pow(gamma, p->nega[i].key) / (gamma + 1));
            fprintf(log, "estimate: %.6f\n", estimate);
        }
    }
    if (p->posibins>0 && p->posi != NULL){
        for(int i = 0; i <p->posibins; ++i){
            fprintf(log, "\tkey: %d,\tcount: %ld,\t", p->posi[i].key, p->posi[i].count);
            double estimate = 2 * pow(gamma, p->posi[i].key) / (gamma + 1);
            fprintf(log, "estimate: %.6f\n", estimate);
        }
    }
    fprintf(log, "\n");
}




//****************************************************** 

struct Bucket *mergeBins(struct Bucket *in, int len1, long pop1, struct Bucket *out, int len2, long pop2, int *count, long *pop) {
 
    int m = len1+len2;                  // at most all keys (len1+len2) are different
    struct Bucket *res = (struct Bucket *) malloc (sizeof(struct Bucket) * m);
    memset(res, 0, sizeof(struct Bucket) * m);
    
    if (len1 == 0) {
        
        if (len2 != 0) {

            #ifdef VERBOSE 
                std::cerr << "mergeBins(): only outbuffer has a positive/negative sketch: get it\n";
            #endif

            memcpy(res, out, len2 * sizeof(struct Bucket));
            (*count) = len2;
            (*pop) = pop2;
            return res;

        } else {
            // else neither ... res remains 0
            
            (*count) = 0;
            (*pop) = 0;
            return res;
        }

    } else if (len2 == 0) {
        
        #ifdef VERBOSE 
            std::cerr << "mergeBins(): only inbuffer has a positive/negative sketch: get it\n"; 
        #endif
        
        memcpy(res, in, sizeof(struct Bucket) * len1);
        (*count) = len1;
        (*pop) = pop1;
        return res;

    } else {

        //effective merge
        (*pop) = pop1 + pop2;
        
        int l = 0, r = 0;
        bool flag1 = true, flag2 = true;
        int p = 0;                          // index of res (distinct keys in res)

        while(m>0){

            if ( flag1 && ( (!flag2) || (in[l].key < out[r].key) ) ) {
                // get from inbuffer only
                // because its key is < than out.key 
                // or because we reached the end of outbuffer

                res[p].key = in[l].key;
                res[p].count = in[l].count;

                #ifdef VERBOSE 
                    std::cerr << p << ")\tMerged.key=" << res[p].key << "\tMerged.count=" << res[p].count << " (inbuffer)"<< std::endl;
                #endif
                
                ++p;            
                m -= 1;
                ++l;
            
            } else if ( (flag1 && flag2) && (in[l].key == out[r].key) ) {
                // equal keys

                res[p].key = in[l].key;
                res[p].count = in[l].count + out[r].count;

                #ifdef VERBOSE 
                    std::cerr << p << ")\tMerged.key=" << res[p].key << "\tMerged.count=" << res[p].count << " (merged)"<< std::endl;
                #endif
                
                ++p;            
                m -= 2;
                ++l;++r;

            } else if (flag2 && ( !flag1 || (in[l].key > out[r].key) ) ) {
                // get from outbuffer
                // because we reached the end of inbuffer
                // or because out.key < in.key

                res[p].key = out[r].key;
                res[p].count = out[r].count;

                #ifdef VERBOSE 
                    std::cerr << p << ")\tMerged.key=" << res[p].key << "\tMerged.count=" << res[p].count << " (outbuffer)"<< std::endl;
                #endif

                ++p;            
                m -= 1;
                ++r;
            
            }//fi keys

            if (l >= len1) {
                flag1 = false;
            }

            if (r >= len2) {
                flag2 = false;
            }

        }//wend m

        (*count) = p; // number of distinct bins in the merged array

        #ifdef VERBOSE 
            long sum = 0;
            for (int t = 0; t < p; ++t) {
                std::cerr << "(" << res[t].key <<", " << res[t].count << ") " << std::endl;
                sum += res[t].count;
            }
            if (sum != (*pop)){
                std::cerr << "Error in mergeBins(): population differs from counted items" << std::endl;
            }
        #endif

        return res;
    }//fi merge 
   
   //to warning off...
   return res;
}


//****************************************************** The MPI_Reduce() custom operator

//*** Reduce positive and negative sketches
void reduceSketchPair(void *inbuf, void *outbuf, int *buffer_size, MPI_Datatype *dptr) {

    //Gbuffer_size is a global variable with value equal to the size of the Reduce msg that never change
    struct PackedPair *inpack = deserializeMsgPair(inbuf, Gbuffer_size);
    struct PackedPair *outpack = deserializeMsgPair(outbuf, Gbuffer_size);
    
    if (!inpack || !outpack){
        std::cout << "Error in Reducing sketches\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    #ifdef VERBOSE
        int r = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &r);
        fprintf(stdout, "reduceSketchPair => %d) PBins: %d, NBins: %d, outPB %d, outNB %d\n", r, inpack->posibins, inpack->negabins, outpack->posibins, outpack->negabins);
    #endif


    struct PackedPair result;
    memset(&result, 0, sizeof(struct PackedPair));

    #ifdef VERBOSE
        std::cerr << "Reducing sketches: a1: "<< inpack->alpha << ", a2: " << outpack->alpha << std::endl;
        debugPack(inpack, NULL);
        debugPack(outpack, NULL);
    #endif

    //*** 1. Get the sketch bound
    result.bound = inpack->bound; //this never changes

    int collapsesA = 0, collapsesB = 0;     // the additional collapses executed before merge (A) and after merge (B)
    
    //*** 2. Get initial alpha
    //if (inpack->alpha != outpack->alpha) {
    if (inpack->collapses != outpack->collapses) {
        
        #ifdef VERBOSE
            std::cerr << "Merge of sketches with different alphas...\n"; 
        #endif

        //if (inpack->alpha < outpack->alpha){
            //collapsesA =collapsePairToAlpha(inpack->posi, &inpack->posibins, inpack->nega, &inpack->negabins, inpack->bound, &inpack->alpha, outpack->alpha);

        if (inpack->collapses < outpack->collapses) {

            collapsesA = collapsePairToMaxCollapses(inpack->posi, &inpack->posibins, inpack->nega, &inpack->negabins, inpack->bound, &inpack->alpha, &inpack->collapses, outpack->collapses);

            #ifdef VERBOSE
                std::cerr << "collapsing inpack: additional collapse "<< collapsesA << std::endl;
            #endif
            
            #ifndef NDEBUG
                assert(inpack->collapses == outpack->collapses);
            #endif
            
            result.alpha = outpack->alpha;
            result.collapses = outpack->collapses;

        } else {

            //collapsesA =collapsePairToAlpha(outpack->posi, &outpack->posibins, outpack->nega, &outpack->negabins, outpack->bound, &outpack->alpha, inpack->alpha);
            collapsesA = collapsePairToMaxCollapses(outpack->posi, &outpack->posibins, outpack->nega, &outpack->negabins, outpack->bound, &outpack->alpha, &outpack->collapses, inpack->collapses);

            #ifdef VERBOSE
                std::cerr << "collapsing outpack: additional collapse "<< collapsesA << std::endl;
            #endif

            #ifndef NDEBUG
                assert(inpack->collapses == outpack->collapses);
            #endif
            
            result.alpha = inpack->alpha;
            result.collapses = inpack->collapses;
        }//fi collapsing

    } else {
        //std::cerr << "No need to collapse incoming sketches\n";
        result.alpha = inpack->alpha;
        result.collapses = inpack->collapses;
    }//fi alpha

    //*** Merge sketches

    //*** 3.1 Merge of positive sketches
    result.posi = mergeBins(inpack->posi, inpack->posibins, inpack->posipop, outpack->posi, outpack->posibins, outpack->posipop, &result.posibins, &result.posipop);

    //*** 3.2 Merge of negative sketches
    result.nega = mergeBins(inpack->nega, inpack->negabins, inpack->negapop, outpack->nega, outpack->negabins, outpack->negapop, &result.negabins, &result.negapop);

    //*** 4. Check sketch bound and update final alpha
    collapsesB = collapseArrayPair(result.posi, &result.posibins, result.nega, &result.negabins, result.bound, &(result.alpha) );
    
    #ifdef VERBOSE 
        std::cerr << "Additional collapses: " << collapsesB << std::endl;
        std::cerr << "New alpha: " << result.alpha << std::endl;
        debugPack(&result, NULL);
    #endif

    //*** 5.Update the final collapses
    result.collapses += collapsesB;

    #ifndef NDEBUG 
    {
        long respop = result.posipop+result.negapop;
        long sump = 0;
        if (result.negapop > 0 && result.negabins > 0 && result.nega) {
            for (int k=0; k<result.negabins;++k){
                sump += result.nega[k].count;
            }
        }//fi nega
        if (result.posibins>0 && result.posi) {
            for (int k=0; k<result.posibins; ++k){
                sump += result.posi[k].count;
            }
        }//fi posi
        assert(sump==respop);
    }
    #endif
    

    //*** 6. Pack the result: serialize and copy in the output buffer
    memset(outbuf, 0, Gbuffer_size * sizeof(char));
    int buffer_position = 0;

    MPI_Pack(&result.alpha, 1, MPI_DOUBLE, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&result.collapses, 1, MPI_INT, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&result.bound, 1, MPI_INT, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);

    MPI_Pack(&result.posibins, 1, MPI_INT, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&result.posipop, 1, MPI_LONG, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    //MPI_Pack(&result.posipop, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);

    for (int i = 0; i < result.posibins; ++i) {
        //MPI_Pack(&result.posi[i].count, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&result.posi[i].key, 1, MPI_INT, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&result.posi[i].count, 1, MPI_LONG, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    }

    MPI_Pack(&result.negabins, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&result.negapop, 1, MPI_LONG, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    //MPI_Pack(&result.negapop, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    
    for (int i = 0; i < result.negabins; ++i) {
        //MPI_Pack(&result.nega[i].count, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&result.nega[i].key, 1, MPI_INT, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&result.nega[i].count, 1, MPI_LONG, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    }

    //*** 7. Free other resources 

    if (result.posi) {
        free(result.posi);
    }

    if (result.nega) {
        free(result.nega);
    }

    if (inpack){
        if (inpack->posi){
            free(inpack->posi);
        }
        if (inpack->nega){
            free(inpack->nega);
        }
        free(inpack);
    }

    if (outpack){
        if (outpack->posi){
            free(outpack->posi);
        }
        if (outpack->nega){
            free(outpack->nega);
        }
        free(outpack);
    }

    //Reduce ended
}






//*************************************** The MPI_Reduce() custom operator for DataDog Sketches

void reduceOriginal(void *inbuf, void *outbuf, int *buffer_size, MPI_Datatype *dptr) {
    
    static int reduce_round = 0;
    int p_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    
    struct PackedPair *inpack = deserializeMsgPair(inbuf, Gbuffer_size);
    struct PackedPair *outpack = deserializeMsgPair(outbuf, Gbuffer_size);
    struct PackedPair result;
    memset(&result, 0, sizeof(struct PackedPair));

    //*** 1. Get the sketch bound, the alpha, the collapses 
    result.bound = inpack->bound;                               // this never changes
    result.alpha = inpack->alpha;                               // this never changes in DDog
    result.collapses = inpack->collapses+outpack->collapses;    // initially these are 0

    #ifdef LOGCollapses
        ++reduce_round;
        fprintf(stderr, "\tOriginal MPI_Reduce => round %d) => P %d)\n", reduce_round, p_rank);
        fprintf(stderr, "P%d) inpack->collapses: %d, outpack->collapses: %d, result: %d\n", p_rank, inpack->collapses, outpack->collapses, result.collapses);
    #endif
    

    //*** 1.1 Check that messages are ok
    #ifndef NDEBUG 
    {
        long tpop = inpack->posipop+inpack->negapop;
        long sum = 0;
        if (inpack->negapop>0 && inpack->negabins>0 && inpack->nega) {
            for (int k=0; k<inpack->negabins;++k){
                sum += inpack->nega[k].count;
            }
        }//fi nega

        if (inpack->posibins>0 && inpack->posi) {
            for (int k=0; k<inpack->posibins; ++k){
                sum += inpack->posi[k].count;
            }
        }//fi posi
       
        assert(sum == tpop);

        tpop = outpack->posipop+outpack->negapop;
        sum = 0;
        if (outpack->negapop>0 && outpack->negabins>0 && outpack->nega) {
            for (int k=0; k<outpack->negabins;++k){
                sum += outpack->nega[k].count;
            }
        }//fi nega
        if (outpack->posibins>0 && outpack->posi) {
            for (int k=0; k<outpack->posibins; ++k){
                sum += outpack->posi[k].count;
            }
        }//fi posi
        
        assert(sum == tpop);
    }
    #endif
    
    
    int pmerge = 0, nmerge = 0; //to take count of collapses on the positive and negative arrays

    //*** 2. Merge sketches
    #ifdef LowBins
        
        /*  collapse on the left, 
            so that if there are negative bins, they are collapsed first
            while B* remains not collapsed
            and the leftmost/lowest keys in the positive sketch are merged into one single bin
        */
        int bound = result.bound;
        bool negaF = false;
        if (inpack->negapop > 0 || outpack->negapop > 0) {
            negaF = true;
            bound -= 1;
        }

        //*** 2.1 Merge of positive sketches
        
        result.posi = optMergeBinsPosi(inpack->posi, inpack->posibins, inpack->posipop, outpack->posi, outpack->posibins, outpack->posipop, bound, &result.posibins, &result.posipop, &pmerge);
        
        #ifdef LOGCollapses
            fprintf(stderr, "DDOG LowBins posi: bound=%d => R: posibins: %d, posipop: %ld, collapses :%d\n", bound, result.posibins, result.posipop, pmerge);
        #endif

        if (negaF) {

            //*** 2.2 Merge negative sketches into the remaining space
            bound = result.bound - result.posibins; // 1 or more

            result.nega = optMergeBinsNega(inpack->nega, inpack->negabins, inpack->negapop, outpack->nega, outpack->negabins, outpack->negapop, bound, &result.negabins, &result.negapop, &nmerge);

            #ifdef LOGCollapses
                fprintf(stderr, "DDOG LowBins nega: bound=%d => R: negabins: %d, negapop: %ld, collapses :%d\n", bound, result.negabins, result.negapop, nmerge);
            #endif
        }//fi negaF

    #else //High Bins

        int bound = result.bound;
        bool posiF = false;

        if (inpack->posibins > 0 || outpack->posibins > 0) {
            posiF = true;
            bound -= 1;
            
            #ifdef LOGCollapses
                fprintf(stderr, "\nreduceOriginal High Bins => inpack->posibins %d, outpack->posibins %d, bound: %d\n", inpack->posibins, outpack->posibins, bound);
            #endif

            if ((inpack->posibins > 1 && inpack->posi[0].key == -MIN_KEY) || (outpack->posibins > 1 &&  outpack->posi[0].key == -MIN_KEY) ) {
                bound -= 1;
            }
        }//fi posibins

        #ifdef LOGCollapses
            fprintf(stderr, "\nreduceOriginal High Bins => at max bound: %d for negabins\n", bound);
        #endif

        //*** 2.1 Merge negative sketches (from highest positive keys to big negatives keys)
        result.nega = optMergeBinsNega(inpack->nega, inpack->negabins, inpack->negapop, outpack->nega, outpack->negabins, outpack->negapop, bound, &result.negabins, &result.negapop, &nmerge);
        
        #ifdef LOGCollapses
            fprintf(stderr, "\nreduceOriginal High Bins Nega => round %d) => P%d) optMergeBinsNega bins = %d, collapses: %d\n", reduce_round, p_rank, result.negabins, nmerge);
        #endif

        bound = result.bound - result.negabins; // 1, 2 or more

        //*** 2.2 Merge of positive sketches
        if (posiF) {
            result.posi = optMergeBinsPosi(inpack->posi, inpack->posibins, inpack->posipop, outpack->posi, outpack->posibins, outpack->posipop, bound, &result.posibins, &result.posipop, &pmerge);
            
            #ifdef LOGCollapses
                fprintf(stderr, "\nreduceOriginal High Bins => P%d) (bound %d) optMergeBinsPosi bins = %d, collapses: %d\n", p_rank, bound, result.posibins, pmerge);
            #endif
        }//fi posiF
    #endif
    

    #ifndef NDEBUG
    {
        long sump = 0;
        if (result.posibins>0 && result.posi) {

            for (int k=0; k<result.posibins; ++k) {
                sump += result.posi[k].count;
            }
        }//fi posi
        //std::cerr << "resultP: "<<result.posipop << " vs " << sump <<std::endl;
        assert(result.posipop == sump);
    
        long sumn = 0;
        if (result.negapop>0 && result.nega) {
            for (int k=0; k<result.negabins; ++k){
                sumn += result.nega[k].count;
            }
        }//fi posi
        //std::cerr << "resultN: "<<result.negapop<<" vs "<<sumn <<std::endl;
        assert(result.negapop == sumn);
        assert((result.posipop + result.negapop)==(sumn+sump));
    }    
    #endif

    //*** 3. Check sketch bound and update final alpha
    int collapsesB = originalCollapseArrayPair(result.posi, &result.posibins, result.nega, &result.negabins, result.bound);
    #ifdef LOGCollapses 
        std::cerr << "Additional collapses: " << collapsesB << std::endl;
        //debugPack(&result, NULL);
    #endif

    //*** 5.Update the final collapses (only on P0)
    //result.collapses += collapsesB;
    //if (!p_rank) {
        result.collapses += collapsesB + pmerge + nmerge;
    //}

    #ifdef LOGCollapses 
        fprintf(stderr, "final value for result.collapses: %d\n", result.collapses);
    #endif

    #ifndef NDEBUG 
    {
        long respop = result.posipop+result.negapop;
        long sumn = 0;
        long sump = 0;
        
        if (result.negapop>0 && result.negabins>0 && result.nega) {
        
            for (int k=0; k < result.negabins;++k){
                sumn += result.nega[k].count;
            }
            //std::cout << "Merger(): Pnega: "<< sumn << std::endl;
        }//fi nega

        if (result.posibins>0 && result.posi) {
            for (int k=0; k < result.posibins; ++k) {
                sump += result.posi[k].count;
                //std::cout<<"\t "<<result.posi[k].key<< " : "<<result.posi[k].count << " sum: "<< sump <<std::endl;
            }
        }//fi posi
        assert((sumn + sump) == respop);
    }
    #endif
    
    //*** 6. Pack the result: serialize and copy in the output buffer
    memset(outbuf, 0, Gbuffer_size * sizeof(char));
    int buffer_position = 0;

    MPI_Pack(&result.alpha, 1, MPI_DOUBLE, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&result.collapses, 1, MPI_INT, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&result.bound, 1, MPI_INT, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);

    MPI_Pack(&result.posibins, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&result.posipop, 1, MPI_LONG, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    //MPI_Pack(&result.posipop, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    
    for (int i = 0; i < result.posibins; ++i) {
        //MPI_Pack(&result.posi[i].count, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&result.posi[i].key, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&result.posi[i].count, 1, MPI_LONG, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    }

    MPI_Pack(&result.negabins, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&result.negapop, 1, MPI_LONG, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    //MPI_Pack(&result.negapop, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    
    for (int i = 0; i < result.negabins; ++i) {
        //MPI_Pack(&result.nega[i].count, 1, MPI_INT, outbuf,Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&result.nega[i].key, 1, MPI_INT, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&result.nega[i].count, 1, MPI_LONG, outbuf, Gbuffer_size, &buffer_position, MPI_COMM_WORLD);
    }


    //*** 7. Free other resources
    if (result.posi) {
        free(result.posi);
    }

    if (result.nega){
        free(result.nega);
    }

    if (inpack) {
        if (inpack->posi){
            free(inpack->posi);
        }

        if (inpack->nega){
            free(inpack->nega);
        }

        free(inpack);
    }

    if (outpack){
        if (outpack->posi){
            free(outpack->posi);
        }

        if (outpack->nega){
            free(outpack->nega);
        }

        free(outpack);
    }

    //Reduce ended
}


//*************************************** Merge bins optimized for DataDog Sketches


struct Bucket *optMergeBinsPosi(struct Bucket *in, int len1, long pop1, struct Bucket *out, int len2, long pop2, int bound, int *count, long *pop, int* collapses) {

    //in = inpack->posi, len1 = inpack->posibins, pop1 = inpack->posipop, 
    //out = outpack->posi, len2 = outpack->posibins, pop2 = outpack->posipop

    #ifdef LOGCollapses
        fprintf(stderr, "optMergeBinsPosi) posibins1=%d, posibins2=%d, bound: %d\n", len1, len2, bound);
    #endif

    // we can fill at most _bound_ bins (but they can be less)
    struct Bucket *sketch = (struct Bucket *) malloc(sizeof(struct Bucket) * bound);
    memset(sketch, 0, sizeof(struct Bucket) * bound);

    //effective merge
    (*pop) = pop1 + pop2;               // total population from the pair of positive sketches
    long r_pop1 = pop1, r_pop2 = pop2;  // residual of population from each sketch 
    long tpop = 0;                      // # of items in the merged sketch
    int collapsed = 0;

    #ifdef LowBins
        
        //collapse lowest/leftmost keys on the positive sketch (sx)
        //so start merging from the right
        
        int idx = bound - 1;                    // index of the current position in the merged sketch (of distinct keys)
        int idx_low = 0;                        // min index available in the merged sketch (to take into account the B* bucket)
        
        bool flag1 = true, flag2 = true;        // flags are true while we are in the range of the arrays _struct Bucket_

        // if _in_ or _out_ arrays contain the B*, this will stay un-collapsed
        long bIn = 0, bOut = 0;

        int lmin = 0;                           // min position of _in_ array        
        int l = (len1 - 1);                     // max position of _in_ array
        if (l < 0) {
            flag1 = false; //for empty array
        }

        if (flag1 && (in[0].key == -MIN_KEY)) {
            bIn = in[0].count;      // the count of items of B* in the _in_ array is saved
            lmin = 1;               // the min index of the _in_ array is 1 and not 0

            if (l < lmin){
                flag1 = false;
            }
        }
        
        int rmin = 0;                           // min position of _out_ array
        int r = (len2 - 1);                     // max position of _out_ array
        if (r < 0) {
            flag2 = false;
        }

        if (flag2 && (out[0].key == -MIN_KEY)) {
            bOut = out[0].count;    // the count of items of B* in the _out_ array is saved
            rmin = 1;               // the min index of the _out_ array is 1 and not 0

            if (r < rmin){
                flag2 = false;
            }
        }

        //don't collapse bucket B*
        if (bIn || bOut) {
            
            // in the merged sketch, the B* is saved first as the lowest possible key
            sketch[0].key = -MIN_KEY;   
            
            // count is the sum of the original bins
            sketch[0].count = bIn + bOut;
            tpop += sketch[0].count;

            // from the remaining population delete the items in B*
            r_pop1 -= bIn; 
            r_pop2 -= bOut;

            // in the merged sketch, the min key is now at position 1
            idx_low = 1;

            #ifdef LOGCollapses
                fprintf(stderr, "Bucket B* has %ld items (%ld + %ld) - bound: %d\n", sketch[0].count, bIn, bOut, bound);
            #endif
        }

        while ((idx >= idx_low) && (flag1 && flag2)) {  
            
            #ifdef LOGCollapses
                fprintf(stderr, "idx %d, in[%d]: %d, out[%d]: %d\n", idx, l, in[l].key, r, out[r].key);
            #endif

            if (out[r].key > in[l].key ) { 
                sketch[idx].key = out[r].key;
                sketch[idx].count = out[r].count;
                r_pop2 -= out[r].count;

                --r;
            } else if (in[l].key == out[r].key)  {
                    
                sketch[idx].key = out[r].key;
                sketch[idx].count = out[r].count + in[l].count;
                r_pop1 -= in[l].count;
                r_pop2 -= out[r].count;

                --l;
                --r;
            } else { // (in[l].key > out[r].key) 
                sketch[idx].key = in[l].key;
                sketch[idx].count = in[l].count;
                r_pop1 -= in[l].count;

                --l;
            }//fi

            tpop += sketch[idx].count;
            --idx;

            if (l < lmin) {
                flag1 = false;
            }
            if (r < rmin) {
                flag2 = false;
            }
        }//wend

        
        while ((idx >= idx_low) && flag1) {

            #ifdef LOGCollapses
                fprintf(stderr, "filling with sketch left; idx %d, in[%d]: %d\n", idx, l, in[l].key);
            #endif
            // there is still room in the merged sketch
            // and
            // only keys in _in_ are left
            sketch[idx].key = in[l].key;
            sketch[idx].count = in[l].count;
            
            r_pop1 -= in[l].count;
            tpop += sketch[idx].count;
            
            --idx;
            --l;
            if (l<lmin){
                flag1 = false;
            }
        }//wend flag1


        while ((idx >= idx_low) && flag2) {

            #ifdef LOGCollapses
                fprintf(stderr, "filling with sketch right; idx %d, out[%d]: %d\n",idx, r, out[r].key);
            #endif

            // there is still room in the merged sketch
            // and
            // only keys in _out_ are left
            sketch[idx].key = out[r].key;
            sketch[idx].count = out[r].count;
            tpop += sketch[idx].count;
            r_pop2 -= out[r].count;
            
            --idx;
            --r;
            if (r < rmin) {
                flag2 = false;
            }
        }//wend 

        if (flag1 || flag2) {

            #ifdef LOGCollapses
                fprintf(stderr, "space exhausted, at the end idx is %d\n", idx);
            #endif
            
            /*  while() ended because space exhausted
                now idx = -1 OR idx = 0 
                and, all the remaning keys go into this final 'trash' bin
                So, we get the remaining maximum key and sum up all the remaining items
            */     
            
            ++idx; //into the last filled bucket went all the remaining items
            sketch[idx].count += r_pop1 + r_pop2;
            tpop += r_pop1 + r_pop2;
            ++collapsed;

            #ifdef LOGCollapses
                fprintf(stderr, "after while, into first (trash) bucket (idx is %d) are %ld items\n", idx, sketch[idx].count);
            #endif

            #ifndef NDEBUG
            {
                // fprintf(stderr, "Total ppop = %ld, pop: %ld\n", tpop, *pop);
                assert(tpop == *pop);

                long check = 0;
                for (int i = 0; i < bound; ++i) {
                    check += sketch[i].count;
                }
                // fprintf(stderr, "Total check = %ld\n", check);
                assert(tpop == check);
            }
            #endif

            (*count) = bound; // number of distinct bins in the merged array
            *collapses = collapsed;
            return sketch;

        } else {

            #ifdef LOGCollapses
                fprintf(stderr, "(idx is %d), no more keys to merge\n", idx);
            #endif

            /*  while() ended because keys exhausted
                so we have to resize the returned pointer...
            */
            if (idx_low == 1) { //we have B* in the merged sketch and we have to slide it right
                //move B* forward
                sketch[idx].key = -MIN_KEY;
                sketch[idx].count = sketch[0].count;
            } else { 
                ++idx;
            }

            (*count) = bound - idx;
            
            struct Bucket *Fsketch = (struct Bucket *) malloc((*count) * sizeof(struct Bucket));
            memset(Fsketch, 0, (*count) * sizeof(struct Bucket));
            memcpy(Fsketch, &sketch[idx], (*count) * sizeof(struct Bucket));

            #ifndef NDEBUG
            {
                // fprintf(stderr, "Total ppop = %ld, pop: %ld\n", tpop, *pop);
                assert(tpop == *pop);

                long check = 0;
                for (int i = 0; i < (*count); ++i) {
                    check += Fsketch[i].count;
                }
                // for (int i = idx; i < bound; ++i) {
                //     check += sketch[i].count;
                // }
                
                // fprintf(stderr, "Total check = %ld\n", check);
                assert(tpop == check);
            }
            #endif

            *collapses = collapsed;
            
            //return &sketch[idx];
            free(sketch);
            return Fsketch;
        }
    
    #else // High Bins

        //collapse highest/rightmost keys (dx)
        //so start merging from the left
        
        int idx = 0;                            // index of sketch (distinct keys)
        int l = 0;
        int r = 0;
        
        long bIn = 0, bOut = 0;

        if ((len2 >= 1 ) && (out[0].key == -MIN_KEY)) {
            bOut = out[0].count;
            r_pop2 -= out[0].count;

            ++r;
        }

        if ((len1 >=  1) && (in[0].key == -MIN_KEY)) {
            bIn = in[0].count;
            r_pop1 -= in[0].count;

            ++l;
        }

        if (bIn || bOut) {
            //don't collapse bucket B*
            sketch[0].key = -MIN_KEY;
            sketch[0].count = bIn + bOut;
            idx = 1;

            tpop += sketch[0].count;
        }

        bool flag1 = true, flag2 = true;
        if (l >= len1) {
            flag1 = false;
        }

        if (r >= len2) {
            flag2 = false;
        }

        while ((flag1 && flag2) && (idx < bound)) {
            
            if (in[l].key < out[r].key) {   
                sketch[idx].key = in[l].key;
                sketch[idx].count = in[l].count;
                r_pop1 -= in[l].count;
                
                ++l;
            } else if (in[l].key == out[r].key) {
                sketch[idx].key = out[r].key;
                sketch[idx].count = out[r].count + in[l].count;
                r_pop1 -= in[l].count;
                r_pop2 -= out[r].count;

                ++r;++l;
            } else {// (in[l].key > out[r].key) 
                sketch[idx].key = out[r].key;
                sketch[idx].count = out[r].count;
                r_pop2 -= out[r].count;
                
                ++r;
            }//fi

            tpop += sketch[idx].count;
            ++idx;
            
            if (l >= len1) {
                flag1 = false;
            }

            if (r >= len2) {
                flag2 = false;
            }
        }//wend

        while (flag1 && idx < bound){
            sketch[idx].key = in[l].key;
            sketch[idx].count = in[l].count;
            r_pop1 -= in[l].count;
            
            ++l;
            if (l >= len1) {
                flag1 = false;
            }

            tpop += sketch[idx].count;
            ++idx;
        }//wend

        while (flag2 && idx < bound){
            sketch[idx].key = out[r].key;
            sketch[idx].count = out[r].count;
            r_pop2 -= out[r].count;
            
            ++r;
            if (r >= len2) {
                flag2 = false;
            }

            tpop += sketch[idx].count;
            ++idx;
        }//wend


        if (flag1 || flag2) {

            /* while() ended because space exhausted 
            now idx = bound
            and, all the remaning keys go into this final 'trash' bin
            So, we get the remaining minimum key and sum up all the remaining items
            */     
           --idx; 
            sketch[idx].count += r_pop1 + r_pop2;
            tpop += r_pop1 + r_pop2; 
            ++collapsed;

            *count = idx + 1; //bound;
            *collapses = collapsed;

            #ifndef NDEBUG
            {
                // fprintf(stderr, "Space exhausted, Total ppop = %ld, pop: %ld\n", tpop, *pop);
                assert(tpop == *pop);

                long check = 0;
                for (int i = 0; i <= idx; ++i) {
                    check += sketch[i].count;
                }
                // fprintf(stderr, "Total check = %ld\n", check);
                assert(tpop == check);
            }
            #endif

            return sketch;
        } else {
            //keys exhausted, idx >= bound
            *count = idx; // number of distinct bins in the merged array

             #ifndef NDEBUG
            {
                // fprintf(stderr, "Total ppop = %ld, pop: %ld\n", tpop, *pop);
                assert(tpop == *pop);

                long check = 0;
                for (int i = 0; i < *count; ++i) {
                    check += sketch[i].count;
                }
                // fprintf(stderr, "Total check = %ld\n", check);
                assert(tpop == check);
            }
            #endif

            return sketch;
        }
    #endif //LowBins vs HighBins
   
   //to warning off...
   return sketch;
}



struct Bucket *optMergeBinsNega(struct Bucket *in, int len1, long pop1, struct Bucket *out, int len2, long pop2, int bound, int *count, long *pop, int* collapses) {

    //in = inpack->nega, len1 = inpack->negabins, pop1 = inpack->negapop, 
    //out = outpack->nega, len2 = outpack->negabins, pop2 = outpack->negapop
    
    // we can fill at most bound bins (but they can be less)
    struct Bucket *sketch = (struct Bucket *) malloc (sizeof(struct Bucket) * bound);
    memset(sketch, 0, sizeof(struct Bucket) * bound);

    //effective merge
    (*pop) = pop1 + pop2;
    long r_pop1 = pop1, r_pop2 = pop2;
    int Ncollapses = 0;
    long tpop = 0;

    #ifdef LowBins
        
        //highest keys on the negative sketch should collapse (be merged)
        //so start merging from 0 
        int l = 0;
        int r = 0;
        int idx = 0;    // index of sketch (distinct keys)
        
        bool flag1 = true, flag2 = true;
        if (!len1) {
            flag1 = false;
        }
        if (!len2) {
            flag2 = false;
        }

        while ((flag1 && flag2) && (idx < bound) ) {
            
            if (in[l].key < out[r].key ) {   
                
                sketch[idx].key = in[l].key;
                sketch[idx].count = in[l].count;
                r_pop1 -= in[l].count;
                
                ++l;
            } else if (in[l].key == out[r].key) {
                
                sketch[idx].key = out[r].key;
                sketch[idx].count = out[r].count + in[l].count;

                r_pop1 -= in[l].count;
                r_pop2 -= out[r].count;
                
                ++r;++l;
            } else { //(in[l].key > out[r].key)
                
                sketch[idx].key = out[r].key;
                sketch[idx].count = out[r].count;
                r_pop2 -= out[r].count;
                
                ++r;
            }//fi

            tpop += sketch[idx].count;
            ++idx;

            if (l >= len1 ) {
                flag1 = false;
            }

            if (r >= len2 ) {
                flag2 = false;
            }
        }//wend

        while (flag1 && idx < bound) {
            // fill the merged sketch with _in_ items
            sketch[idx].key = in[l].key;
            sketch[idx].count = in[l].count;
            r_pop1 -= in[l].count;
                
            ++l;

            tpop += sketch[idx].count;
            ++idx;
        }//wend

        while (flag2 && idx < bound) {
            // fill the merged sketch with _out_ items
            sketch[idx].key = out[r].key;
            sketch[idx].count = out[r].count;
            r_pop2 -= out[r].count;
                
            ++r;

            tpop += sketch[idx].count;
            ++idx;
        }//wend

        if (flag1 || flag2) {

            /*  while() ended because space exhausted
                idx==bound

                all the remaning keys go into this final 'trash' bin
                So, we get the remaining maximum key and sum up all the remaining items
            */  

            --idx; //in the last bucket add the remaining items 
            sketch[idx].count += r_pop1 + r_pop2;
            tpop += r_pop1 + r_pop2;
            ++Ncollapses;

            #ifndef NDEBUG
            {
                // fprintf(stderr, "Total tpop = %ld, pop: %ld\n", tpop, *pop);
                assert(tpop == *pop);

                long check = 0;
                for (int i = 0; i < bound; ++i) {
                    check += sketch[i].count;
                }
                // fprintf(stderr, "Total check = %ld\n", check);
                assert(tpop == check);
            }
            #endif


            *collapses = Ncollapses;
            (*count) = bound; // number of distinct bins in the merged array
            return sketch;

        } else {
            /*  while() ended because keys exhausted
                so we have to resize the returned pointer...
            */
           (*count) = idx + 1;

            #ifndef NDEBUG
            {
                // fprintf(stderr, "Total tpop = %ld, pop: %ld\n", tpop, *pop);
                assert(tpop == *pop);

                long check = 0;
                for (int i = 0; i < (*count); ++i) {
                    check += sketch[i].count;
                }
                // fprintf(stderr, "Total check = %ld\n", check);
                assert(tpop == check);
            }
            #endif
            
            return sketch;
        }

    #else //High bins
        
        // start reading from the biggest keys and merging the lowest ones
        
        int idx = bound - 1;    // index of sketch (distinct keys)
        int l = len1 - 1;
        int r = len2 - 1;
        bool flag1 = true, flag2 = true;
        if (!len1){
            flag1 = false;
        }
        if (!len2){
            flag2 = false;
        }

        while ((flag1 && flag2) && (idx >= 0)) {
            
            if (out[r].key > in[l].key) {   
                sketch[idx].key = out[r].key;
                sketch[idx].count = out[r].count;
                r_pop2 -= out[r].count;

                --r;
            } else if (in[l].key == out[r].key) {
                
                sketch[idx].key = out[r].key;
                sketch[idx].count = out[r].count + in[l].count;
                r_pop1 -= in[l].count;
                r_pop2 -= out[r].count;

                --l;
                --r;
            } else { // (in[l].key > out[r].key)
                
                sketch[idx].key = in[l].key;
                sketch[idx].count = in[l].count;                
                r_pop1 -= in[l].count;
                
                --l;
            }//fi

            tpop += sketch[idx].count;
            --idx;
            
            if (l < 0) {
                flag1 = false;
            }

            if (r < 0) {
                flag2 = false;
            }
        }//wend

        while (flag1 && idx >= 0) {
            sketch[idx].key = in[l].key;
            sketch[idx].count = in[l].count;                
            r_pop1 -= in[l].count;
            
            --l;
            if (l < 0) {
                flag1 = false;
            }

            tpop += sketch[idx].count;
            --idx;
        }//wend


        while (flag2 && idx>= 0) {
            sketch[idx].key = out[r].key;
            sketch[idx].count = out[r].count;                
            r_pop2 -= out[r].count;
            
            --r;
            if (r < 0) {
                flag2 = false;
            }
            
            tpop += sketch[idx].count;
            --idx;
        }//wend

        if (flag1 || flag2) { 

            /* while() ended because space exhausted 
            now idx = -1
            and, all the remaning keys go into this final 'trash' bin
            So, we get the remaining minimum key and sum up all the remaining items
            */ 

            ++idx; //to go back to the last bin
            sketch[idx].count += r_pop1 + r_pop2;
            tpop += r_pop1 + r_pop2;
            Ncollapses += 1;

            *count = bound - idx;

            #ifdef LOGCollapses
                fprintf(stderr, "space exhausted, negabins: %d\n", *count);
            #endif

            #ifndef NDEBUG
            {
                // fprintf(stderr, "Total tpop = %ld, pop: %ld\n", tpop, *pop);
                assert(tpop == *pop);

                long check = 0;
                for (int i = 0; i < *count; ++i) {
                    check += sketch[i].count;
                }
                // fprintf(stderr, "Total check = %ld\n", check);
                assert(tpop == check);
            }
            #endif

            *collapses = Ncollapses;
            return sketch;
        } else {
            /* while() ended because keys exhausted 
                so idx >= 0
            */
            ++idx; // step back to stay on the last filled bin

            *count = (bound - idx); // number of distinct bins in the merged array

            #ifdef LOGCollapses
                fprintf(stderr, "keys exhausted, negabins: %d, resizing\n", *count);
            #endif

            struct Bucket *Fsketch = (struct Bucket *)malloc((*count) * sizeof(struct Bucket));
            memset(Fsketch, 0, (*count) * sizeof(struct Bucket));
            memcpy(Fsketch, &sketch[idx], (*count) * sizeof(struct Bucket));
            free(sketch);

            #ifndef NDEBUG
            {
                // fprintf(stderr, "Total tpop = %ld, pop: %ld\n", tpop, *pop);
                assert(tpop == *pop);

                long check = 0;
                for (int i = 0; i < (*count); ++i) {
                    check += Fsketch[i].count;
                }
                // fprintf(stderr, "Total check = %ld\n", check);
                assert(tpop == check);
            }
            #endif

            //return &sketch[idx];
            *collapses = Ncollapses;
            return Fsketch;
        }
    #endif //LowBins vs HighBins


   //to warning off...
   return sketch;
}

