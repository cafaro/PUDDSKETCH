/**
 * @file ArraySketch.cc
 * @author CM
*/

#include "ArraySketch.h"

//************************************************************ QUANTILES ESTIMATIONS

double PairQuantile(double q, struct PackedPair *p, int *index, long *bcount) {
    // if ( (q < 1) && (q >= 0) ) 
    
    double gamma = (1+p->alpha)/(1-p->alpha);
    long n = p->posipop + p->negapop;
    long threshold = q*(n-1);

    long count = 0;             // sum of items in the sketch
    int bkey = 0;               // bkey, bucket key
    long current_count = 0;     // count of the bucket where quantile falls
    int sign = 0;               // -1 or +1, for nega and posi sketch
    int i_p = -1, i_n=-1;

    if (p->negabins>0 && p->nega != NULL) {
        i_n = 0;

        bkey = p->nega[i_n].key;
        count = p->nega[i_n].count;
        sign = -1;
        ++i_n;

        current_count = count;
    }//fi nega
        

    if (p->posibins>0 && p->posi != NULL) {
            i_p = 0;
        
            // If there are no bins in the negative sketch, 
            // initialize bkey and count with the first positive bucket
            if (sign==0) {
                bkey = p->posi[i_p].key;
                count = p->posi[i_p].count;
                sign = 1;
                ++i_p;

                current_count = count;
            }//fi sign
    }//fi posi
    
    
    //@note: to be sure this will not become an infinite loop, add a second check condition on i_p
    while (count <= threshold && (i_p < p->posibins)) {

        if (i_n >= 0 && i_n < p->negabins){
            
            sign = -1;
            bkey = p->nega[i_n].key;
            count += p->nega[i_n].count;
            
            current_count = p->nega[i_n].count; //current bucket's count
        
        ++i_n;
        } else if (i_p >= 0 && i_p < p->posibins){
        
            sign = +1;
            bkey = p->posi[i_p].key;
            count += p->posi[i_p].count;
            
            current_count = p->posi[i_p].count; //current bucket's count

        ++i_p;
        } 
    }//wend count

    double quantileEstimate = (sign * 2 * pow(gamma, bkey))/(gamma + 1);
    (*index) = bkey;
    (*bcount) = current_count;
    return quantileEstimate; 
}







//****************************************************** ARRAY-BASED sketches and operations


//*** used in Collapsing an already merged sketch 
struct Bucket *collapseBins(struct Bucket *sketch, int count, int *size) {

    struct Bucket *result = (struct Bucket *) malloc( sizeof(struct Bucket)*count );
    int j = 0;
    int i = 0;

    while (i < count) {
    
        int key = sketch[i].key; 
        int count = sketch[i].count;
        
        if (key == -MIN_KEY) {
            //B* is never collapsed
            result[j].key = key;
            result[j].count = count;
            ++j;
            ++i;
        } else {

            if (key%2 != 0) {
                
                // key is odd
                // if key is positive, we look for the next even key (i.e. 3 -> 4)
                // if key is negative, we look for the next even key (i.e. -5 -> -4)
                int nextk = key+1;


                // the new 'collapsed' key will be 1/2 of the even key (i.e. 4 -> 2, -4 -> -2)
                // while with positive key, the new key 2 has already been 'collapsed' itself, -2 has not yet been processed...
                // so we use an auxiliary sketch
                int newK = (key+1)/2;

                if (sketch[i+1].key != nextk) {
                    // the successive even key doesn't exist 
                    result[j].key = newK;
                    result[j].count = count;
                    i += 1;
                    ++j;
                } else {
                    // the successive even key exists, these two buckets must 'sum up' in the new 'collapsed' key 
                    result[j].key = newK;
                    result[j].count = count + sketch[i+1].count;
                    i += 2;
                    ++j;
                }
                    
            }//fi key odd 
            else 
            {   /* 
                KEY IS EVEN : (-6, -4, 0, 2, .. 8 ...)
                Prev_key = key-1 :  (-7, -5, -1, 1, .. 7...)

                Se il bucket dispari precedente - corrispondente alla key (key-1) - esistesse, 
                avremmo già fatto il collapse di quello col bucket pari corrente. 
                Evidentemente non c'è, quindi modifichiamo solo la key di questo bucket...

                ATTENZIONE: mentre con key positive la riduzione verso 0 sovrascrive key già processate,
                con key negative procedendo allo stesso modo si va a sovrascrivere bucket non ancora processati
                (a meno di usare floor(key-1/2) e tendere a sinistra anche per key negative)
                */
                
                int newK = (key)/2; 
                
                result[j].key = newK;
                result[j].count = count;
                ++i;        
                ++j;
            }//fi key even
        }//fi -MIN_KEY

    }//wend

    *size = j;
    return result;
}




//*** Used in reduceSketchPair(): collapses a pair of positive and negative sketches until the sketch bound is met
int collapseArrayPair(struct Bucket *posi, int *count_p, struct Bucket *nega, int *count_n, int bound, double *alpha) {

    struct Bucket *p = NULL;
    struct Bucket *n = NULL;
    int sizep = (*count_p);
    int sizen = (*count_n);
    double alpha_out = *alpha;
    int sizep_out = 0, sizen_out = 0;

    int collapses = 0;
    while ( (sizep + sizen) > bound) {
        
        p = collapseBins(posi, sizep, &sizep_out);
        n = collapseBins(nega, sizen, &sizen_out);
        
        memset(posi, 0, sizeof(struct Bucket)*sizep);
        memcpy(posi, p, sizeof(struct Bucket)*sizep_out);
        sizep = sizep_out;
        free(p);

        memset(nega, 0, sizeof(struct Bucket)*sizen);
        memcpy(nega, n, sizeof(struct Bucket)*sizen_out);
        sizen = sizen_out;
        free(n);
        
        alpha_out = (2*alpha_out)/(1 + pow(alpha_out,2));
        
        ++collapses;
    }//wend

    (*alpha) = alpha_out;
    (*count_p) = sizep;
    (*count_n) = sizen;

    return collapses;
}



//******************************************* DDog version

int originalCollapseArrayPair(struct Bucket *posi, int *count_p, struct Bucket *nega, int *count_n, int bound) {
    
    int sizeP = (*count_p);
    int sizeN = (*count_n);
    
    struct Bucket pp[sizeP];
    memcpy(pp, posi, sizeof(struct Bucket)*sizeP);
    
    struct Bucket nn[sizeN];
    memcpy(nn, nega, sizeof(struct Bucket)*sizeN);
    
    int lowP = 0;
    int BB = 0;

    if (pp[lowP].key == -MIN_KEY) {
        BB = pp[lowP].count;
        ++lowP;

        #ifdef LOGCollapses
            std::cout<<"originalCollapseArrayPair(): B* contains "<< BB << " items and is not collapsed" <<std::endl;
        #endif 
    }//fi don't collapse B*

    int maxP = sizeP-1;
    
    int lowN = sizeN-1;
    int maxN = 0;

    int Ppop = 0;
    int Npop = 0;

    int collapses = 0;

    while ( (sizeP + sizeN) > bound) {

        // #ifdef LOGCollapses
        //     std::cout<<"originalCollapseArrayPair "<< collapses+1 <<") "<<sizeP << " + " <<sizeN <<std::endl;
        // #endif 
 
        #ifdef LowBins

            if (sizeN > 1)
            {
                //collapse Negative starting from highest keys
                int i = lowN; //max key
                int j = i-1;

                nn[j].count += nn[i].count; //ultimo nel penultimo 

                #ifdef LOGCollapses
                    std::cout << collapses+1 <<") Negative merge: key " << nn[i].key << " -> " << nn[j].key << std::endl;
                #endif

                nn[i].key = 0;
                nn[i].count = 0;
                lowN = j;

                sizeN -= 1;
                ++collapses;

                //std::cout<<"originalCollapseArrayPair(): Negative collapse: bucket "<<nn[j].key<< " count is: "<<nn[j].count <<std::endl;
                // #ifndef NDEBUG
                // {
                //     int sum = 0;
                //     for(int c = 0; c <= lowN; ++c){
                //         sum += nn[c].count;
                //         std::cout<<"\t "<<nn[c].key<< " : "<<nn[c].count << " sum: "<< sum <<std::endl;
                //     }
                //     std::cout<<"Nega pop after collapse: "<<sum<<std::endl;
                // }
                // #endif

            } 
            else 
            {   //collapse Positive starting from lowest keys
                
                int i = lowP;
                int j = i+1;

                pp[j].count += pp[i].count;

                #ifdef LOGCollapses
                    std::cout <<"\t"<< collapses+1 <<") Positive merge: key " << pp[i].key << " -> " << pp[j].key << std::endl;
                #endif

                pp[i].key = 0;
                pp[i].count = 0;
                lowP = j;

                sizeP -= 1;
                ++collapses;

                //std::cout<<"originalCollapseArrayPair(): Positive collapse: bucket "<<pp[j].key<< " count is: "<<pp[j].count <<std::endl;
                // #ifndef NDEBUG
                // {
                // int sum = BB;
                // for(int c = lowP; c < (*count_p); ++c){
                //     sum += pp[c].count;
                //     std::cout<<"\t "<<pp[c].key<< " : "<<pp[c].count << " sum: "<< sum <<std::endl;
                // }
                // std::cout<<"Posi pop after collapse: "<<sum<<std::endl;
                // }
                // #endif
            }
        #else
            if ( (pp[0].key == -MIN_KEY && sizeP <=2 ) || (pp[0].key != -MIN_KEY && sizeP < 2) ) 
            {
                //collapse negative sketch
                int i = maxN;
                int j = i+1;

                nn[j].count += nn[i].count;

                #ifdef LOGCollapses
                    std::cout << collapses+1 <<") Negative merge: key " << nn[i].key << " -> " << nn[j].key << std::endl;
                #endif

                nn[i].key = 0;
                nn[i].count = 0;
                maxN = j;

                sizeN -= 1;
                ++collapses;

                //std::cout<<"originalCollapseArrayPair(): Negative collapse: bucket "<<nn[j].key<< " count is: "<<nn[j].count <<std::endl;
                // #ifndef NDEBUG
                // {
                // int sum = 0;
                // for(int c = maxN; c < (*count_n); ++c){
                //     sum += nn[c].count;
                //     std::cout<<"\t "<<nn[c].key<< " : "<<nn[c].count << " sum: "<< sum <<std::endl;
                // }
                // std::cout<<"Nega pop after collapse: "<<sum<<std::endl;
                // }
                // #endif
            } 
            else 
            {   //collapse Positive sketch

                int i = maxP;
                int j = i-1;

                pp[j].count += pp[i].count;

                #ifdef LOGCollapses
                    std::cout <<"\t"<< collapses+1 <<") Positive merge: key " << pp[i].key << " -> " << pp[j].key << std::endl;
                #endif
               
                pp[i].key = 0;
                pp[i].count = 0;
                maxP = j;

                sizeP -= 1;
                ++collapses;

                //std::cout<<"originalCollapseArrayPair(): Positive collapse: bucket "<<pp[j].key<< " count is: "<<pp[j].count <<std::endl;
                // #ifndef NDEBUG
                // {
                // int sum = 0;
                // for(int c = 0; c < sizeP; ++c){
                //     sum += pp[c].count;
                // }
                // std::cout<<"Posi pop after collapse: "<<sum<<std::endl;
                // }
                // #endif

            }//fi
        #endif
    }//wend
    
    memset(posi, 0, sizeof(struct Bucket)*(*count_p));
    memset(nega, 0, sizeof(struct Bucket)*(*count_n));

    #ifdef LowBins

        //copy back negative sketch
        for(int i = 0; i <= lowN; ++i) {
            nega[i].key = nn[i].key;
            nega[i].count = nn[i].count;
        }
        //memcpy(nega, nn, sizeof(struct Bucket)*sizeN);


        //copy back positive sketch
        int i = 0;
        //int sum = 0;
        //B* bucket remains unchanged
        if (pp[0].key == -MIN_KEY) {
            posi[i].key = pp[0].key;
            posi[i].count = pp[0].count;
            //sum += posi[i].count;
        ++i;
        }//fi
        //std::cout<< "Copy back from "<<lowP <<" with key: "<<pp[lowP].key << std::endl;   
        
        for(int k = lowP; k < (*count_p); ++k) {
            posi[i].key = pp[k].key;
            posi[i].count = pp[k].count;
            ++i;
            //sum +=posi[i].count; 
        }
        //memcpy(&posi[1], &pp[lowP], sizeof(struct Bucket)*sizeP);
    
    #else
        //copy back negative sketch
        int som = 0;
        int j = 0;
        //std::cout<< "Copy back Nega from "<<maxN <<" with key: "<<nn[maxN].key;   
        for(int h=maxN; h < (*count_n); ++h){
            nega[j].key = nn[h].key;
            nega[j].count = nn[h].count;
            ++j;
            //som += nega[j].count;
        }
        //std::cout<<  " and sum "<<som <<std::endl;

        //memcpy(nega, &nn[maxN], sizeof(struct Bucket)*sizeN);
        
        //copy back positive sketch
        for(int h = 0; h <= maxP; ++h) {
            posi[h].key = pp[h].key;
            posi[h].count = pp[h].count;
        }
        //memcpy(posi, pp, sizeof(struct Bucket)*sizeP);

    #endif

    (*count_p) = sizeP;
    (*count_n) = sizeN;

    return collapses;
}




//****************************************************** COLLAPSE ARRAY-BASED 




/**
 * In the Reduce, before the merge of pair of sketches, 
 * we have to check that both sketches have the same alpha,
 * so we need to collapse the sketch that has alpha less than the other's alpha
 */
int collapsePairToAlpha(struct Bucket *posi, int *count_p, struct Bucket *nega, int *count_n, int bound, double *alpha, double maxAlpha) {

    struct Bucket *p = NULL;
    struct Bucket *n = NULL;
    int sizep = (*count_p);
    int sizen = (*count_n);
    double alpha_out = *alpha;
    int collapse = 0;
    int sizep_out = 0, sizen_out = 0;

    while (alpha_out < maxAlpha) {
        
        p = collapseBins(posi, sizep, &sizep_out);
        n = collapseBins(nega, sizen, &sizen_out);

        // p = collapsePositiveBuckets(posi, sizep, &sizep_out);
        // n = collapseNegativeBuckets(nega, sizen, &sizen_out);

        memset(posi, 0, sizeof(struct Bucket)*sizep);
        memcpy(posi, p, sizeof(struct Bucket)*sizep_out);
        sizep = sizep_out;
        free(p);

        memset(nega, 0, sizeof(struct Bucket)*sizen);
        memcpy(nega, n, sizeof(struct Bucket)*sizen_out);
        sizen = sizen_out;
        free(n);

        alpha_out = 2*(alpha_out/(1 + pow(alpha_out,2)));

        ++collapse;
    }//wend

    *alpha = alpha_out;
    (*count_p) = sizep;
    (*count_n) = sizen;

    return collapse;
}



/**
 * In the Reduce, before the merge of pair of sketches, 
 * we have to check that both sketches have been collapsed the same number of times
 * so we need additional collapse until the other's max
 */
int collapsePairToMaxCollapses(struct Bucket *posi, int *count_p, struct Bucket *nega, int *count_n, int bound, double *alpha, int *mincollapses, int maxCollapses) {

    struct Bucket *p = NULL;
    struct Bucket *n = NULL;
    int sizep = (*count_p);
    int sizen = (*count_n);
    int sizep_out = 0, sizen_out = 0;
    
    double out_alpha = (*alpha);
    int collapses = (*mincollapses);
    int loop = 0; //additional collapses

    while (collapses < maxCollapses) {
        
        p = collapseBins(posi, sizep, &sizep_out);
        n = collapseBins(nega, sizen, &sizen_out);

        memset(posi, 0, sizeof(struct Bucket)*sizep);
        memcpy(posi, p, sizeof(struct Bucket)*sizep_out);
        sizep = sizep_out;
        free(p);

        memset(nega, 0, sizeof(struct Bucket)*sizen);
        memcpy(nega, n, sizeof(struct Bucket)*sizen_out);
        sizen = sizen_out;
        free(n);

        out_alpha = (2*out_alpha)/(1 + pow(out_alpha,2));

        ++collapses;
        ++loop;
    }//wend

    *alpha = out_alpha;
    (*mincollapses) = collapses;
    (*count_p) = sizep;
    (*count_n) = sizen;
    
    return loop;
}


