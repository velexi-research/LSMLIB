##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for utility subroutines
##

SRC_DIR = .

BUILD_DIR = ../../..
include $(BUILD_DIR)/config/Makefile.config

all:  includes library

includes:  
	cp -f $(SRC_DIR)/lsm_calculus_toolbox.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_calculus_toolbox2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_calculus_toolbox2d_local.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_calculus_toolbox3d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_utilities1d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_utilities2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_utilities2d_local.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_utilities3d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_utilities3d_local.h $(BUILD_DIR)/include/

library:  lsm_calculus_toolbox.o         \
	  lsm_calculus_toolbox.o         \
	  lsm_calculus_toolbox2d.o       \
	  lsm_calculus_toolbox2d_local.o \
	  lsm_calculus_toolbox3d.o       \
          lsm_utilities1d.o              \
          lsm_utilities2d.o              \
	  lsm_utilities2d_local.o        \
          lsm_utilities3d.o              \
	  lsm_utilities3d_local.o        \

clean:
	rm -f *.o 

