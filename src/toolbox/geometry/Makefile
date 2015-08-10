##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for geometry subroutines
##

SRC_DIR = .

BUILD_DIR = ../../..
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE


all:  includes library

includes:  
	cp -f $(SRC_DIR)/lsm_geometry1d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_geometry2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_geometry2d_local.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_geometry3d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_curvature2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_curvature2d_local.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_curvature3d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/lsm_curvature3d_local.h $(BUILD_DIR)/include/

library:  lsm_curvature2d.o                 \
          lsm_curvature2d_local.o           \
          lsm_curvature3d.o                 \
          lsm_curvature3d_local.o           \
          lsm_geometry1d.o                  \
          lsm_geometry2d.o                  \
          lsm_geometry2d_local.o            \
          lsm_geometry3d_fort.o             \
          lsm_geometry3d_c.o                \

clean:
	rm -f *.o 
	cd test; make clean || exit 1

