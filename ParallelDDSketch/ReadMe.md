
# Parallel DDSketch
Based on the DDSketch implementation, we derived a Parallel DDSketch with Uniform Collapse procedure

The DDOG version has a LowBins collapse (when collapse is performed on the leftmost buckets,
starting from the Negative Sketch and jumping to the Positive one when in the former there are less than 2 bins),
and a HighBins collapse (when collapse is performed on the rightmost buckets, starting from the Positive Sketch)
In each case the special bucket handling 0 and near-0 values is never collapsed.


## Folder arrangement
Each folder has a special purpose:
- PLogs contains the logs of the MPI processes
- BinaryLogs contains, for each process, a binary file with the processed items
- binFiles contains the input binary files

## Modes
The code can be executed in 2 modes:
- VALIDATE, that is to compare results against the sequential version, processing and input binary file
- RUN, in which each process self-generates its local data from a given distribution

# Compile and Run
To speed up the testing phase, I've specialized the Makefile to let each script compile its own version
of binary file and to execute its run.


## How to test
Code can be tested in:
- sequential mode
- parallel mode

If you want to test the original DataDog sketch, compile its binary:

$ make -f MakefileO.mk

In Sequential mode simply launch (this is in VALIDATION mode since it's getting its input from file)

$ ./SketcherDog -f binFiles/uniformN.dat -a 0.001 -b 512

For the Parallel test (this is in RUN mode since it's generating its own input)

$ mpirun -np 8 ./SketcherDog -d 1 -x 1 -y 10000 -n 1.25 -a 0.001 -b 500

### ATTENTION
When in RUN mode, the value of -n is intended to be multiplied for 10^6.
In the example, the buffer len is 10^6*1.25=1250000 items for each process.
