##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for building MATLAB MEX-files for spatial derivatives
##

SRC_DIR = .

BUILD_DIR = ../../..
include $(BUILD_DIR)/config/Makefile.config

# include MATLAB configuration
include $(BUILD_DIR)/config/Makefile.config.MATLAB


library: mex
	if [ ! -d "$(BUILD_DIR)/lib/matlab" ]; then mkdir $(BUILD_DIR)/lib/matlab; fi
	cp -f *.m $(BUILD_DIR)/lib/matlab
	cp -f *. $(BUILD_DIR)/lib/matlab
	if [ -d "$(MATLAB_INSTALL_DIR)/LSMLIB" ]; then cp -f *.m $(MATLAB_INSTALL_DIR)/LSMLIB/; fi
	if [ -d "$(MATLAB_INSTALL_DIR)/LSMLIB" ]; then cp -f *. $(MATLAB_INSTALL_DIR)/LSMLIB/; fi

mex:  matlab

matlab:  HJ_ENO1_1D.            \
         HJ_ENO2_1D.            \
         HJ_ENO3_1D.            \
         HJ_WENO5_1D.           \
         HJ_ENO1_2D.            \
         HJ_ENO2_2D.            \
         HJ_ENO3_2D.            \
         HJ_WENO5_2D.           \
         HJ_ENO1_3D.            \
         HJ_ENO2_3D.            \
         HJ_ENO3_3D.            \
         HJ_WENO5_3D.           \
         UPWIND_HJ_ENO1_1D.     \
         UPWIND_HJ_ENO2_1D.     \
         UPWIND_HJ_ENO3_1D.     \
         UPWIND_HJ_WENO5_1D.    \
         UPWIND_HJ_ENO1_2D.     \
         UPWIND_HJ_ENO2_2D.     \
         UPWIND_HJ_ENO3_2D.     \
         UPWIND_HJ_WENO5_2D.    \
         UPWIND_HJ_ENO1_3D.     \
         UPWIND_HJ_ENO2_3D.     \
         UPWIND_HJ_ENO3_3D.     \
         UPWIND_HJ_WENO5_3D.

clean:
		rm -f *.
		rm -f *.o 

