/**
 * @file Summary.h
 * @author CM
 * 
 * @brief Implementation of a Sketch with uniform collapse  
 * 
 */


#ifndef __SUMMARY_H__
#define __SUMMARY_H__

#include "Header.h"

/**
 * @brief Sketches serialization for Reduce
 */
char *pack_SketchesPair(double alpha, int collapses, int bound, std::map<int, long>& posiSketch, int posibins, long posipop, std::map<int, long>& negaSketch, int negabins, long negapop, int *final_size);


/**
 * @brief Sketches deserialization for Reduce
 */
struct PackedPair *deserializeMsgPair(void *buffer, int size);


#endif //__SUMMARY_H__
