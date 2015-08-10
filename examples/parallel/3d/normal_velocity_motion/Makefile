##
## File:        Makefile.in
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: makefile for 3d parallel level set method example programs
##

BUILD_DIR = ../../../..

# configuration for SAMRAI library
include /work/02697/jrigelo/opt/SAMRAI/include/../config/Makefile.config

# NOTE:  LSMLIB's Makefile.config MUST be included last because it
#        redefines the pattern rules for compiling source code
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE


CPPFLAGS_EXTRA = -I$(FORTRAN)

all: TestLSM_3d_main

TestLSM_3d_main:                             \
    TestLSM_3d_main.o                        \
    TestLSM_3d_VelocityFieldModule.o         \
    TestLSM_3d_PatchModule.o                 \
    testlsm_3d_patchmodule_fort.o            \
    testlsm_3d_velocityfield_fort.o
	$(CXX) -O3 -funroll-loops -fPIC -g -O2 $(CXXFLAGS) $(LDFLAGS) -L$(LSMLIB_LIB_DIR) $^          \
        $(LSMLIB_LIBS) $(LDLIBS_EXTRA)                                     \
        $(LIBSAMRAI3D) $(LIBSAMRAI2D) $(LIBSAMRAI1D) $(LIBSAMRAI)          \
        $(LDLIBS) -o TestLSM_3d_main
     
clean:
		rm -f TestLSM_3d_main
		rm -f *.o core 

cleanrun:
		rm -f PI*
		rm -f *.log
		rm -f *.log.*
		rm -f -r *.visit
		rm -f -r *.restart

spotless: clean cleanrun

include Makefile.depend

FORTRAN = fortran

testlsm_3d_patchmodule_fort.o: $(FORTRAN)/testlsm_3d_patchmodule_fort.f
	@F77@ -O3 -fPIC  $(FFLAGS) $(FFLAGS_EXTRA) -c $^ 

testlsm_3d_velocityfield_fort.o: $(FORTRAN)/testlsm_3d_velocityfield_fort.f
	@F77@ -O3 -fPIC  $(FFLAGS) $(FFLAGS_EXTRA) -c $^ 

