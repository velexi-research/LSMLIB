##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description:	makefile for parallel level set method example programs
##

BUILD_DIR = ../../..
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE


# configuration for SAMRAI library
include /work/02697/jrigelo/opt/SAMRAI/include/../config/Makefile.config

SUBDIRS = advection field_extension orthogonalization toolbox

all:
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

clean:
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

cleanrun:
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

spotless:	clean cleanrun
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

