/**
 * @file Merger.h 
 * @author CM
 * @brief The Parallel Reduction Merge function definition and all utility routines
*/

#ifndef __MERGE_H__
#define __MERGE_H__

#include "Summary.h"
#include "ArraySketch.h"
#include "MapSketch.h"
#include "Utility.h"
#include "Header.h"


/**
 * @brief The MPI_Reduce() custom-defined operator
 */
void reduceSketchPair(void *inbuf, void *outbuf, int *buffer_size, MPI_Datatype *dptr);


/**
 * @brief The function merges 2 arrays of type struct Bucket into a final one
 */
struct Bucket *mergeBins(struct Bucket *in, int len1, long pop1, struct Bucket *out, int len2, long pop2, int *count, long *pop);
 

/**
 * @brief Deserialize the MPI_PACKED resulting at the end of the MPI_Reduce()
 */
struct PackedPair *deserializeGlobalSketchPair(void *buffer, int size);

/**
 * @brief Work on the final Global sketch (into which the negative sketch keys permutation was in the reversed order)
 */
void debugGlobalSketchPair(struct PackedPair *p, FILE *log);


//**************************************** DataDog

/**
 * @brief The MPI_Reduce() custom-defined operator for the original DataDog sketch collapse procedure
 */
void reduceOriginal(void *inbuf, void *outbuf, int *buffer_size, MPI_Datatype *dptr);

//struct Bucket *optMergeBinsPosi(struct Bucket *in, int len1, long pop1, struct Bucket *out, int len2, long pop2, int bound, int *count, long *pop);
struct Bucket *optMergeBinsPosi(struct Bucket *in, int len1, long pop1, struct Bucket *out, int len2, long pop2, int bound, int *count, long *pop, int* collapses);

struct Bucket *optMergeBinsNega(struct Bucket *in, int len1, long pop1, struct Bucket *out, int len2, long pop2, int bound, int *count, long *pop, int* collapses);


#endif //__MERGE_H__

