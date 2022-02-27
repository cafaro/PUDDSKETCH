/**
 * @file MapSketch.cc
 * @author CM
 */

#include "MapSketch.h"
extern int pid0;


//****************************************************** UTILITY FUNCTIONS

int addKeyToSketch(std::map<int,long>&sketch, int key) {
    int res = 0;

    std::map<int,long>::iterator it = sketch.find(key);
    if ( it == sketch.end() ) {   
        //key not in sketch, add new key
        sketch[key] = 1;
        res = 1;
    } else {
        //key already inserted in the sketch, increase count
        sketch[key] += 1;
        res = 0;
    }//fi find key

    return res;
}


void debugSketch(std::map<int,long>& mySketch, int rank, FILE *log){    
    for(auto it = mySketch.begin(); it != mySketch.end(); ++it){
        fprintf(log, "%d, %ld\n", it->first , it->second);
    }
    fprintf(log, "\n");
}




//****************************************************** COLLAPSE for MAP-BASED SKETCHES

//*** Uniform Collapse for both positive and negative sketches
void collapseUniformly(std::map<int,long>& mySketch) {
    
    // key = ⎡log(x)/log(γ) ⎤, ⎡log(γ)/log(γ) ⎤ = 1, ⎡log(1)/log(γ) ⎤ = 0, ⎡log(1/γ)/log(γ) ⎤ = -1
    // when collapsing:
    // key = 1 and key = 2 go into k*=1
    // key = 3 and key = 4 go into k*=2
    // ...
    // key = -1 and key = 0 go into k*=0
    // key = -3 and key = -2 go into k*=-1
    // ...
    // @note: - std::ceil(x) != std::ceil(-x)
    
    std::map<int,long> newSketch; // where new collapsed buckets will be stored    
    std::map<int,long>::iterator it;

    // for(it = mySketch.begin(); it != mySketch.end(); ++it){
    //     if (it->first == -MIN_KEY) { // -MIN_KEY: bucket B* (for 0 and near-0 values): -(2^30)      
    //         newSketch[-MIN_KEY] += it->second;
    //     } else {
    //         double k = it->first;
    //         int k_new = std::ceil(k/2);
    //         newSketch[k_new] += it->second;
    //     }//fi
    // }//for

    //this IF is met only for positive sketches and only for the first bucket
    it = mySketch.begin();
    if (it->first == -MIN_KEY) { // -MIN_KEY: bucket B* (for 0 and near-0 values): -(2^30)
        newSketch[-MIN_KEY] += it->second;
        ++it;
    }
    while(it != mySketch.end() ){
        double k = it->first;
        int k_new = std::ceil(k/2);
        newSketch[k_new] += it->second;
        ++it;
    }//wend

    //copy back
    mySketch.swap(newSketch);
}

//**************************************** DataDog Collapse

int OriginalPairCollapse(std::map<int,long>& posiSketch, int *posibins, std::map<int,long>& negaSketch, int *negabins, int bound, double *BNega, double *BPosi, FILE *log) {
    
    int collapses = 0;
    
    int sizep = *posibins;
    int sizen = *negabins;
    
    #ifdef LOGCollapses
        fprintf(stderr, "OriginalPairCollapse: sizep: %d, sizen: %d\n", sizep, sizen);
    #endif

    //*** collapse sketches for BOUND - if necessary
    while ( (sizep + sizen) > bound) 
    {
        
        #ifdef LowBins
            if (sizen > 1)
            {   //collapse the first two bins (with the highest keys) in the negative Sketch
                
                std::map<int,long>::reverse_iterator rp1, rp2;
                rp1 = negaSketch.rbegin();//last bucket
                rp2 = rp1;//before last
                
                ++rp2;
                rp2->second += rp1->second;//ultimo nel penultimo
                collapses += 1;

                #ifdef LOGCollapses
                    if (log) {
                        fprintf(log, "DDOG Low Bins, negative sketch, biggest key: ");
                        fprintf(log, "%d in %d\n", rp1->first, rp2->first);
                    } 
                //std::cout << "Key "<<rp1->first << " goes into "<<rp2->first<<std::endl;
                #endif

                negaSketch.erase(rp1->first); //delete last bucket
            
                #ifndef NDEBUG
                assert(negaSketch.size()==sizen-1);
                #endif
                
                sizen = negaSketch.size(); //get new size
                
                *BNega = rp2->first;
            } 
            else 
            {   
                // 0 or 1 bin in the negative Sketch => collapse Positive Sketch
                //collapse the first two bins on the lowest keys (except for B*)

                std::map<int,long>::iterator p1, p2;
                p1 = posiSketch.begin();
                
                if (p1->first == -MIN_KEY)
                {
                    ++p1;
                }//fi don't collapse B*

                p2 = p1;
                ++p2;

                p2->second += p1->second;
                collapses += 1;

                #ifdef LOGCollapses
                    if (log) {
                        fprintf(log, "DDOG Low Bins, positive sketch, lowest key: ");
                        fprintf(log, "%d in %d\n", p1->first, p2->first);
                    }
                #endif

                posiSketch.erase(p1->first);
                
                #ifndef NDEBUG
                assert(posiSketch.size()==sizep-1);
                #endif

                sizep = posiSketch.size();
                *BPosi = p2->first;
            }//fi sizen

        #else //High bins collapsing

            //@note 8.4 correction, don't merge B*
            int hasZero = posiSketch.count(-MIN_KEY); //1 if key exists, 0 otherwise

            if ((hasZero && sizep > 2) || (!hasZero && sizep > 1)) {  

                #ifdef LOGCollapses
                    fprintf(stderr, "high collapse: hasZero: %d, sizep: %d\n", hasZero, sizep);
                #endif
                //collapse Positive Sketch
                //collapse the last two bins on the highest keys

                std::map<int,long>::reverse_iterator rp1, rp2;
                rp1 = posiSketch.rbegin();

                rp2 = rp1;
                ++rp2;

                //@note 8.4 correction, don't merge B*
                //if (rp2->first != -MIN_KEY) {
                    
                    rp2->second += rp1->second;
                    collapses += 1;
                    
                    #ifdef LOGCollapses
                        if (log) {
                            fprintf(log, "DDOG High Bins, positive sketch, biggest key: ");
                            fprintf(log, "%d in %d\n", rp1->first, rp2->first);
                        }
                    //std::cout<<"count is "<<rp2->second<<std::endl;
                    #endif

                    posiSketch.erase(rp1->first);

                    //std::map<int,int>::reverse_iterator c = posiSketch.rbegin();
                    //std::cout<< "new last bins: " << c->first<<std::endl;

                    sizep = posiSketch.size();
                    //std::cout<<"Collapsing Positive Sketch: final Size: "<<sizep<<std::endl;
                    
                    *BPosi = rp2->first;
                //}//fi don't collapse B*
            }//fi hasZero
            else
            {   //collapse negative Sketch (because we have 0 or 1 bin in the positive Sketch)
                //collapse the last two bins on the lowest keys
                #ifdef LOGCollapses
                    fprintf(stderr, "high collapse: sizen: %d\n", sizen);
                #endif
                std::map<int,long>::iterator p1, p2;
                p1 = negaSketch.begin();

                p2 = p1;
                ++p2;
                
                p2->second += p1->second;
                collapses += 1;
                
                #ifdef LOGCollapses
                    if (log) {
                        fprintf(log, "DDOG High Bins, negative sketch, lowest key: ");
                        fprintf(log, "%d in %d\n", p1->first, p2->first);
                    }
                #endif
                
                negaSketch.erase(p1->first);

                sizen = negaSketch.size();
                *BNega = p2->first;
            }
        #endif

    }//wend collapses
    
    *posibins = sizep;
    *negabins = sizen;

return collapses;
}


//**************************************** DataDog

void OriginalCollapse(SketchT *S) {

    //*** collapse sketches for BOUND - if necessary
    while ((S->posibins + S->negabins) > S->mbound) 
    {
        #ifdef LowBins
            if (S->negabins > 1)
            {   //collapse negative Sketch
                //collapse the first two bins on the highest keys 
                
                std::map<int,long>::reverse_iterator rp1, rp2;
                rp1 = S->negaSketch.rbegin();
                rp2 = rp1;
                ++rp2;
                rp2->second += rp1->second;

                S->negaSketch.erase(rp1->first);
                S->negabins = S->negaSketch.size();
                S->process_collapses += 1;
            } 
            else
            {   //collapse positive Sketch
                //collapse the first two bins on the lowest keys (except for B*)
                
                std::map<int,long>::iterator p1, p2;
                p1 = S->posiSketch.begin();
                
                if (p1->first == -MIN_KEY){
                    ++p1;
                }//fi don't collapse B*
                p2 = p1;
                ++p2;
                
                p2->second += p1->second;
                S->posiSketch.erase(p1->first);

                S->posibins = S->posiSketch.size();
                S->process_collapses += 1;
            }
        #else
            if (S->posibins > 1)
            {   //collapse Positive Sketch
                //collapse the last two bins on the highest keys

                std::map<int,long>::reverse_iterator rp1, rp2;
                rp1 = S->posiSketch.rbegin();
                rp2 = rp1;
                ++rp2;
                
                rp2->second += rp1->second;
                S->posiSketch.erase(rp1->first);

                S->posibins = S->posiSketch.size();
                S->process_collapses += 1;
            }
            else
            {   //collapse negative Sketch
                //collapse the last two bins on the lowest keys
                
                std::map<int,long>::iterator p1, p2;
                p1 = S->posiSketch.begin();
                p2 = p1;
                ++p2;
                
                p2->second += p1->second;
                S->posiSketch.erase(p1->first);

                S->negabins = S->negaSketch.size();
                S->process_collapses += 1;

            }//fi
        #endif
    }//wend collapses
}




void fillSketchForDDog(double *buffer, int bsize, SketchT *S) {

    int key = 0;                            // bucket key for items
    /* estimate the key:
     * key = ⎡log(x)⎤   if x > 0
     * key = ⎡log(-x)⎤  if x < 0
     * key = B*         if |x| < ε
     */

    //*** Step 1: each process fills its local sketch
    for(int i = 0; i < bsize; ++i) 
    {
        #ifdef RUN

            /* Add to positive sketch only
             * because we are sure that each item is >> 0
             * @note see generateData()
             */

            //if (buffer[i] > NULLBOUND) {
            
                key = std::ceil(std::log10(buffer[i])/S->base);
                S->posibins += addKeyToSketch(S->posiSketch, key);
                S->posipop += 1;

            //}//fi (1)

            // /* In RUN MODE we consider only positive items,
            //  * so every negative item or near 0 falls in the NULL Bucket B* (-∞)
            //  * 
            //  * buffer[i] <= NULLBOUND
            //  */
            // else { 
            //     key = -MIN_KEY; //special bucket B* (for 0 and near-0 values): -(2^30)
            //     posibins += addKeyToSketch(posiSketch, key);
            //     ++posipop;
            // }//fi (2) 

        #else
            
            /* In VALIDATION MODE consider both positive and negative items: 
             * items are parsed from a file 
             * and we don't know their sign and magnitude
             */

            if (buffer[i] > S->NULLBOUND) 
            {   //add to positive sketch
               
                key = std::ceil(std::log10(buffer[i])/S->base);
                S->posibins += addKeyToSketch(S->posiSketch, key);
                S->posipop += 1;
            
            }//fi (1)
            else if ( -S->NULLBOUND <= buffer[i] && buffer[i] <= S->NULLBOUND) 
            {  //  |buffer[i]| <= NULLBOUND         
                
                key = -MIN_KEY; //special bucket B* (for 0 and near-0 values): -(2^30)
                S->posibins += addKeyToSketch(S->posiSketch, key);
                S->posipop +=1;
            
            } 
            else 
            {   // buffer[i] < -NULLBOUND
                //add to negative sketch
                key = std::ceil(std::log10(-buffer[i])/S->base);
                S->negabins += addKeyToSketch(S->negaSketch, key);
                S->negapop += 1;

            }//posi vs nega

        #endif 
        
        //check for the memory bound
        OriginalCollapse(S);

    }//for fill local sketch
}



