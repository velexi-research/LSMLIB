##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: top-level makefile for source tree (does NOT build library 
##              archive)
##

SRC_DIR = .

BUILD_DIR = ..
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE



SUBDIRS = toolbox      \
          serial       \
          parallel


all:  includes library

includes: 
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

library: 
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

clean:
	rm -f *.o
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done
	cd matlab; make clean   # clean out matlab directory

spotless: clean
	cd matlab; make clean   # clean out matlab directory
