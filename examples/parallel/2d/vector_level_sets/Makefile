##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for parallel level set method example programs
##

BUILD_DIR = ../../../..

# configuration for SAMRAI library
include /work/02697/jrigelo/opt/SAMRAI/include/../config/Makefile.config

# NOTE:  LSMLIB's Makefile.config MUST be included last because it
#        redefines the pattern rules for compiling source code
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE


CPPFLAGS_EXTRA = -I$(FORTRAN)

all: main

main:                             \
    main.o                        \
    VelocityFieldModule.o         \
    PatchModule.o                 \
    patchmodule_fort.o            \
    velocityfield_fort.o
	$(CXX) -O3 -funroll-loops -fPIC -g -O2 $(CXXFLAGS) $(LDFLAGS) -L$(LSMLIB_LIB_DIR) $^          \
        $(LSMLIB_LIBS) $(LDLIBS_EXTRA)                                     \
        $(LIBSAMRAI3D) $(LIBSAMRAI2D) $(LIBSAMRAI1D) $(LIBSAMRAI)          \
        $(LDLIBS) -o main
     
clean:
		rm -f main
		rm -f *.o core 
		rm -f PI*

cleanrun:
		rm -f PI*
		rm -f *.log
		rm -f *.log.*
		rm -f -r *.visit
		rm -f -r *.restart

spotless: clean cleanrun

include Makefile.depend

FORTRAN = fortran

patchmodule_fort.o: $(FORTRAN)/patchmodule_fort.f
	@F77@ -O3 -fPIC  $(FFLAGS) $(FFLAGS_EXTRA) -c $^ 

velocityfield_fort.o: $(FORTRAN)/velocityfield_fort.f
	@F77@ -O3 -fPIC  $(FFLAGS) $(FFLAGS_EXTRA) -c $^ 

