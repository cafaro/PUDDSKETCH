/**
 * @file Utility.h 
 * @author CM
 * @brief Utility functions: command line arguments, input file partitioning, debug, ...
*/

#ifndef __UTILITY_H__
#define __UTILITY_H__

#include "Header.h"

/**
 * @brief Recap command line configuration
 */
void usage(char *msg);


/**
 * @brief It parses command line to configure running instance
 * Returns 0 if configuration is ok, EPARAMS if command line arguments are not correctly set
 */
int check_params(int argc, char *argv[], Params *p);


/**
 * @brief Recap the running conditions
 */
void logStartup(Params *p, int nprocs, int pid);

void debugParams(Params *p);



/**
 * @brief Open a log file for each running process
 */
FILE *logProcess(Params *p);


/**
 * @brief Open a binary file in write mode to log the dataset of the running process
 */
FILE *logBinary(Params *p);


/**
 * @brief Initialize MPI library...
 */
int warming_up_MPI(int *argc, char **argv[], int *comm_size, int *process_rank);


/**
 * @brief Exit function to deallocate MPI and to close log file
 */
// void endExecution(int rank, FILE *logfd, int ecode);


/**
 * @brief Read the input file and load the local data sample (shard)
 * 
 * This function reads a binary file into a data slot
 * based on the file size in byte and on the process rank
 */
double *readDataSlot(char *binfilename, int comm_size, int process_rank, long *data_size, FILE *log, long *nitems, long *maxrows);

 

/**
 * @brief This function generate a data slot from a given distribution
 */
double *generateData(Params *p, int comm_size, int process_rank);


/**
 * @brief Log the local buffer to the log file
 */
void debugLocalShard(int rank, double *buffer, long size, FILE *fd);


/**
 * @brief Log an MPI_PACKED message
 */
void debugPack(struct PackedPair *p, FILE *log);


/**
 * @brief Used to set-up a SketchT structure to hold a pair of sketches
 */
void initSketchT(SketchT *S, double alpha, int bound);


#endif //__UTILITY_H__

