/**
 * @file Quantiles.h
 * @author CM
 * @brief QuickSelect for order statistics
 * 
 * (Some sample C code for the quickselect algorithm, taken from Numerical Recipes in C):
 * http://www.stat.cmu.edu/~ryantibs/median/quickselect.c
 * 
 */

#ifndef __QUANTILES_H__
#define __QUANTILES_H__


#include "Header.h"


#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

/**
 * @brief In an array of data, of type double, with length 'len', the statistic of order 'pos' is given by:
 * 
 * @param data, the dataset/window
 * @param len, the len of the dataset/window
 * @param pos, the order statistic to calculate (passed as order-1)
 * 
 * @return the pos-th element of the array
 * 
 * @note: 
 * the min of the data array is for pos = 0 (order = first, 1st)
 * the max is order len-1 (order = len-th)
 * the median is at (len-1)/2 (if len is odd)
 * and so on...
 */
double quickselect(double *data, int len, int pos);




/**
 * @brief Computes exact order statististics over the aggregated dataset
 * 
 * @param gdataset, the dataset 
 * @param dlen, the len of the dataset
 * @param eLen, the out params with the number of order statistics calculated 
 * @param fp, the log file
 * 
 * @return the array of exact stats
 */
double *getExactQuantiles(double *gdataset, long dlen, int *eLen, FILE *fp);

#endif // __QUANTILES_H__
