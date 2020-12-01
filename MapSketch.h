/**
 * @file MapSketch.h
 * @author CM
 * 
 * @brief Implementation of a Sketch with uniform collapse using a std::map container
 * 
 */


#ifndef __MAPSKETCH_H__
#define __MAPSKETCH_H__


#include "Header.h"

/**
 * @brief Given a key, the function adds the key to the sketch 
 * it returns 1 if a new bucket was created
 * it returns 0 if key was added to an existing bucket 
 */
int addKeyToSketch(std::map<int,long>&sketch, int key);


/**
 * @brief Simplified uniform collapse of the sketch
 */
void collapseUniformly(std::map<int,long>& mySketch);


/**
 * @brief Logging of the sketch 
 */
void debugSketch(std::map<int,long>& mySketch, int rank, FILE *log);

//****************************************************** DataDog vs UDD


//**************************************** UDD

void UDDCollapse(SketchT *S);

int UDDPairCollapse(std::map<int,long>& posiSketch, int *posibins, std::map<int,long>& negaSketch, int *negabins, int bound, double *alpha, double *gamma, double *base);

void fillSketchForUDD(double *buffer, int bsize, SketchT *S);



//**************************************** DataDog

void OriginalCollapse(SketchT *S);

int OriginalPairCollapse(std::map<int,long>& posiSketch, int *posibins, std::map<int,long>& negaSketch, int *negabins, int bound, double *BNega, double *BPosi, FILE *log);

void fillSketchForDDog(double *buffer, int bsize, SketchT *S);

#endif //__MAPSKETCH_H__

