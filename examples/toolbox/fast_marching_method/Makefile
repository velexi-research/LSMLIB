##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for fast marching method example/test code
##

BUILD_DIR = ../../..
include $(BUILD_DIR)/config/Makefile.config

LIB_DIRS     = -L$(LSMLIB_LIB_DIR)

PROGS        = fmm_heap_test                         \

all:  $(PROGS)

fmm_heap_test:  fmm_heap_test.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

clean:
	rm -f *.o 
	rm -f *.dat
	rm -f $(PROGS)
		
spotless:  clean
	rm -f *.log
