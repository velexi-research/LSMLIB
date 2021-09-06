##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for installing MATLAB scripts for TVD-RK time 
##              integration

SRC_DIR = .

LSMLIB_BUILD_DIR = ../../..
include $(LSMLIB_BUILD_DIR)/config/Makefile.config

# include MATLAB configuration
include $(LSMLIB_BUILD_DIR)/config/Makefile.config.MATLAB

library:  
	if [ ! -d "$(LSMLIB_BUILD_DIR)/lib/matlab" ]; then mkdir $(LSMLIB_BUILD_DIR)/lib/matlab; fi
	cp -f *.m $(LSMLIB_BUILD_DIR)/lib/matlab
	if [ -d "$(MATLAB_INSTALL_DIR)/LSMLIB" ]; then cp -f *.m $(MATLAB_INSTALL_DIR)/LSMLIB/; fi

matlab:

clean:
		rm -f *.
		rm -f *.o 

