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

LIB_DIRS     = -L$(LSMLIB_LIB_DIR) -L.

PROGS   = curvature_model
	
OBJLIST = lsm_options.o curvature_model_top.o curvature_model3d.o \
          curvature_model3d_local.o curvature_model.o

all:  $(PROGS)

curvature_model: $(OBJLIST)
	gcc -O3 -funroll-loops -fPIC -g -O2 $(LIB_DIRS) $(INCLUDE) $(OBJLIST) \
	$(LSMLIB_LIBS) -o $@

clean:
	$(RM) *.o 
	$(RM) $(PROGS)

spotless:  clean
	rm -f data_final*
	rm -f data_init*
	rm -f grid*
	rm -f mask*
	rm -f out_file*

