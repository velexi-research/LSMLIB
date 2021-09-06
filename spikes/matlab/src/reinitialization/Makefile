##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for building MATLAB MEX-files for reinitialization
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

matlab:  COMPUTE_REINIT_EQN_RHS_2D.   \
         COMPUTE_REINIT_EQN_RHS_3D.   \

clean:
		rm -f *.
		rm -f *.o 

