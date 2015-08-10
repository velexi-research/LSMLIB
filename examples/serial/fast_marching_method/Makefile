##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for fast marching method example code
##

BUILD_DIR = ../../..
include $(BUILD_DIR)/config/Makefile.config

LIB_DIRS     = -L$(LSMLIB_LIB_DIR)

PROGS        = computeDistanceFunction2d        \
               computeExtensionFields2d         \
               computeDistanceFunction3d        \
               solveEikonalEquation2d           \

all:  $(PROGS)

computeDistanceFunction2d: computeDistanceFunction2d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) $< $(LSMLIB_LIBS) -o $@

computeExtensionFields2d: computeExtensionFields2d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) $< $(LSMLIB_LIBS) -o $@

computeDistanceFunction3d: computeDistanceFunction3d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) $< $(LSMLIB_LIBS) -o $@

solveEikonalEquation2d: solveEikonalEquation2d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) $< $(LSMLIB_LIBS) -o $@

clean:
	rm -f *.o 
	rm -f *.dat
	rm -f $(PROGS)
		
spotless:  clean
	rm -f *.log
