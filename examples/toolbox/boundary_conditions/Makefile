##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for boundary condition example code
##

BUILD_DIR = ../../..
include $(BUILD_DIR)/config/Makefile.config

LIB_DIRS     = -L$(LSMLIB_LIB_DIR)

PROGS        =                                      \
               neumann_bc_ENO1_1d                   \
               neumann_bc_ENO2_1d                   \
               neumann_bc_ENO3_1d                   \
               neumann_bc_WENO5_1d                  \
               neumann_bc_ENO1_2d                   \
               neumann_bc_ENO2_2d                   \
               neumann_bc_ENO3_2d                   \
               neumann_bc_WENO5_2d                  \
               neumann_bc_ENO1_3d                   \
               neumann_bc_ENO2_3d                   \
               neumann_bc_ENO3_3d                   \
               neumann_bc_WENO5_3d                  \

all:  $(PROGS)

neumann_bc_ENO1_1d:  neumann_bc_ENO1_1d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_ENO2_1d:  neumann_bc_ENO2_1d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_ENO3_1d:  neumann_bc_ENO3_1d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_WENO5_1d:  neumann_bc_WENO5_1d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_ENO2_2d:  neumann_bc_ENO2_2d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_ENO3_2d:  neumann_bc_ENO3_2d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_WENO5_2d:  neumann_bc_WENO5_2d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_ENO1_2d:  neumann_bc_ENO1_2d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_ENO2_3d:  neumann_bc_ENO2_3d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_ENO1_3d:  neumann_bc_ENO1_3d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_ENO3_3d:  neumann_bc_ENO3_3d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

neumann_bc_WENO5_3d:  neumann_bc_WENO5_3d.o
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) -o $@ $^ $(LSMLIB_LIBS)

clean:
	rm -f *.o 
	rm -f *.dat
	rm -f $(PROGS)
		
spotless:  clean
	rm -f *.log
