#
# ************************************************************
# @file Makefile
# @author CM
# @brief compiling Parallel DDSketch,
#			DDOG version, with collapse on the left
# ************************************************************
#

NAME=$(shell uname -s)

ifeq ($(NAME),Linux)
	CC=mpicxx
	CFLAGS=-std=c++11 -O3
	EXE=mpirun
else
	CC=mpicxx
	CFLAGS=-std=c++11 -Os
	EXE=mpirun
endif

LFLAGS=#-L.
LIBS=#-l

SRCDIR=./src
IFLAGS=-I$(SRCDIR)

DEPS=$(SRCDIR)/ParallelSketcher.cc $(SRCDIR)/Quantiles.cc $(SRCDIR)/Merger.cc $(SRCDIR)/ArraySketch.cc $(SRCDIR)/MapSketch.cc $(SRCDIR)/Summary.cc $(SRCDIR)/Utility.cc

# MODE: 
# VALIDATE, uses input from file: used to compare results with the sequential version 
# RUN, uses self-generated distribution: used to test parallel running times
MODE=-DRUN 

# -DLOG: used to enable log to file for each process
# -DOFF: to disable log to file
LOG=-DLOG# -DOFF#   


# if LOGB is enabled, a binary file containing the processed items for each process is saved in BinaryLogs
COPY=#-DLOGB

TYPE=-DDDOG

#what bins to collapse? lowest or highest?
BINS=-DLowBins# -DHighBins # 

V=-DLOGCollapses# -DVERBOSE 

TARGET=SketcherDogL# SketcherDogH# 

all: $(TARGET)

# define NDEBUG to disable assert checks -DNDEBUG  
# $(CC) $(CFLAGS) -o $@ $(DEPS) $(IFLAGS) -DNDEBUG $(MODE) $(TYPE) $(BINS) $(LOG) $(COPY) $(LFLAGS) $(LIBS) $(V)

$(TARGET):$(DEPS)
	@echo "Compiling for" $(NAME)
	$(CC) $(CFLAGS) -o $@ $(DEPS) $(IFLAGS) -DNDEBUG $(MODE) $(TYPE) $(BINS) $(LOG) $(COPY) $(LFLAGS) $(LIBS)


clean:
	rm -f $(TARGET)
	rm -f PLogs/*.* BinaryLogs/*.*

