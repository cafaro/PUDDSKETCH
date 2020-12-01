/**
 * @file ArraySketch.h 
 * @author CM
 * @brief Functions that operate on sketch based on array 
 * The array is derived from a std::map structure, so that keys in the array are inserted into ascending order
*/

#ifndef __ARRAYSKETCH_H__
#define __ARRAYSKETCH_H__


#include "Header.h"


/**
 * @brief get quantile approximate value from negative and positive sketches
 */
double PairQuantile(double q, struct PackedPair *p, int *index, long *bcount);


/**
 * @brief It collapse an array representing a sketch
 */
struct Bucket *collapseBins(struct Bucket *sketch, int count, int *size);


/**
 * @brief The function collapses positive and negative array representing the sketches
 * 
 * Used in reduceSketchPair(): collapses a pair of positive and negative sketches until the sketch bound is met
 */
int collapseArrayPair(struct Bucket *posi, int *count_p, struct Bucket *nega, int *count_n, int bound, double *alpha);



/**
 * In the input to the Reduce, we can merge sketches if both have the same alpha,
 * so we need to collapse the sketch that has alpha less than the other's alpha
 * 
 */
int collapsePairToAlpha(struct Bucket *posi, int *count_p, struct Bucket *nega, int *count_n, int bound, double *alpha, double maxAlpha);


/** @brief In the Reduce, before the merge of pair of sketches, 
 * we have to check that both sketches have been collapsed the same number of times
 * so we need additional collapse until the other's max
 */
int collapsePairToMaxCollapses(struct Bucket *posi, int *count_p, struct Bucket *nega, int *count_n, int bound, double *alpha, int *mincollapses, int maxCollapses);



//**************************************** DataDog


int originalCollapseArrayPair(struct Bucket *posi, int *count_p, struct Bucket *nega, int *count_n, int bound);


#endif //__ARRAYSKETCH_H__

