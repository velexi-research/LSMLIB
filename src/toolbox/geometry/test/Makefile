##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for geometry test programs
##

SRC_DIR = .

BUILD_DIR = ../../../..
include $(BUILD_DIR)/config/Makefile.config

CFLAGS_EXTRA = -I$(LSMLIB_INCLUDE)

TEST_PROGRAMS = test_find_line_in_tetrahedron_1       \
                test_find_line_in_tetrahedron_2

all:   $(TEST_PROGRAMS)

test_find_line_in_tetrahedron_1:  test_find_line_in_tetrahedron_1.o
	gcc -o $@ $^ ../lsm_geometry3d_c.o -lm

test_find_line_in_tetrahedron_2:  test_find_line_in_tetrahedron_2.o
	gcc -o $@ $^ ../lsm_geometry3d_c.o -lm

clean:
	rm -f $(TEST_PROGRAMS)
	rm -f *.o

spotless:  clean

