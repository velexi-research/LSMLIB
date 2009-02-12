#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# 
# Autoconf input scripts for lsmlib.

#===========================================================================
# INITIALIZATION
#===========================================================================
AC_PREREQ(2.57)
AC_INIT([LSMLIB], 1.0.0, [ktchu@serendipityresearch.org])
AC_REVISION($Revision$)
AC_CONFIG_SRCDIR([src/parallel/LevelSetMethodAlgorithm.cc])
AC_CONFIG_HEADER([config/LSMLIB_config.h])
AC_CONFIG_AUX_DIR([config/autoconf])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
case $target in
  *powerpc*darwin* | *ppc*darwin*)
    AC_SUBST([arch],darwin)
  ;;
  *86*darwin*)
    AC_SUBST([arch],darwin86)
  ;;
  *cygwin*)
    AC_SUBST([arch],linux32)
  ;;
  *32*linux*)
    AC_SUBST([arch],linux32)
  ;;
  *64*linux*)
    AC_SUBST([arch],linux64)
  ;;
  *bgl*)
    AC_SUBST([arch],bgl)
  ;;
  *)
    AC_SUBST([arch],linux32)
    AC_MSG_WARN([m4_text_wrap([Unknown target architecture...assuming linux...])])
  ;;
esac

# Clear the default compilation flags
CFLAGS=
CXXFLAGS=
FFLAGS=


#===========================================================================
# Check for programs
#===========================================================================

# compilers and linkers
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_F77

# set Fortran compiler flags
case $F77 in
  xlf) 
    FFLAGS="$FFLAGS -qextname"
    COMPILER_FULL_PATH=`which xlf`
    COMPILER_BIN_PATH=`dirname ${COMPILER_FULL_PATH}`
    LIB_DIR_EXTRA="-L`dirname ${COMPILER_BIN_PATH}`/lib"
    LIBS_EXTRA="-lxlfmath -lxlf90 -lxl"
    LDFLAGS="${LIB_DIR_EXTRA} ${LIBS_EXTRA}"
  ;;
  *)
esac
AC_SUBST([FFLAGS])
AC_SUBST([LIB_DIR_EXTRA])
AC_SUBST([LIBS_EXTRA])

# utilities
AC_CHECK_PROGS(CP, cp) CP="$CP -f" AC_SUBST([CP])
AC_CHECK_PROGS(RM, rm) RM="$RM -f" AC_SUBST([RM])
AC_CHECK_PROGS(MV, mv) MV="$MV -f" AC_SUBST([MV])
AC_CHECK_PROGS(MKDIR, mkdir)
AC_CHECK_PROGS(AR, ar)
AC_CHECK_PROGS(XARGS, xargs)
AC_CHECK_PROGS(SED, sed)
AC_CHECK_PROGS(AWK, awk)
AC_CHECK_PROGS(FIND, find)
AC_CHECK_PROGS(MAKE, make)
AC_CHECK_PROGS(INSTALL, install)
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AS_IF([test $arch = "darwin"],
      [RANLIB="$RANLIB -s"
      AC_SUBST([RANLIB])]
     )


#==========================================================================
# Checks for libraries.
#==========================================================================

# check for SAMRAI library
samrai_dir=""  # initialize samrai_dir
AC_MSG_CHECKING([SAMRAI library])
AC_ARG_WITH([samrai],
            [AC_HELP_STRING([--with-samrai=DIR],
            [Specify the location of SAMRAI (default [NO])])],
            [], [withval=no])
AS_IF([test $withval = no],
      [AC_MSG_RESULT([m4_text_wrap([configured without SAMRAI...parallel functionality unavailable])])],
      [samrai_dir=$with_samrai 
       AC_MSG_RESULT([$samrai_dir])])
AS_IF([! test -z "$samrai_dir"],
      [AC_CHECK_FILE([$samrai_dir/include/SAMRAI_config.h],
                     [samrai_inc="$samrai_dir/include"
                      samrai_lib="$samrai_dir/lib"
                      AC_SUBST(samrai_dir)
                      AC_SUBST(samrai_inc)
                      AC_SUBST(samrai_lib)],
                     [AC_MSG_ERROR([SAMRAI library appears to be incomplete or invalid])
                     ])
      ])

# Check for mpi if SAMRAI has been configured
AS_IF([! test -z "$samrai_dir"], 
      [CHECK_MPI("LSMLIB")]) 

# set PARALLEL_PACKAGE if configured with SAMRAI
AS_IF([test -z "$samrai_dir"], 
      [PARALLEL_PACKAGE="" LSMLIB_PARALLEL=""],
      [PARALLEL_PACKAGE="parallel" LSMLIB_PARALLEL="-llsm_parallel"])
AC_SUBST([PARALLEL_PACKAGE])
AC_SUBST([LSMLIB_PARALLEL])

# check for MATLAB libraries
matlab_dir=""  # initialize matlab_dir
AC_MSG_CHECKING([MATLAB configuration])
AC_ARG_WITH([matlab],
            [AC_HELP_STRING([--with-matlab=DIR],
            [Specify the MATLAB patch (default [NO])])],
            [], [withval=no])
AS_IF([test $withval = no], 
      [AC_MSG_RESULT([m4_text_wrap([configured without MATLAB...MEX-file generation unavailable])])],
      [matlab_dir=$with_matlab
       AC_MSG_RESULT([$matlab_dir])]
)
AC_SUBST(matlab_dir)

# set mex flags
AS_IF([! test -z "$matlab_dir"],
      [AC_CHECK_FILE([$matlab_dir/bin/mex],
                     [mex_flags='-I$(BUILD_DIR)/include -L$(BUILD_DIR)/lib -llsm_serial -llsm_toolbox'],
                     [AC_MSG_ERROR([MATLAB installation appears to be incomplete or invalid])])
])
AC_SUBST(mex_flags)

# set mex-file suffix if library configured with MATLAB
AS_IF([! test -z "$matlab_dir"],
      [case $arch in
         darwin)
           AC_SUBST([mex_extension],mexmac)
         ;;
         darwin86)
           AC_SUBST([mex_extension],mexmaci)
         ;;
         linux32)
           AC_SUBST([mex_extension],mexglx)
         ;;
         linux64)
           AC_SUBST([mex_extension],mexa64)
         ;;
         windows)
           AC_SUBST([mex_extension],dll)
         ;;
         solaris)
           AC_SUBST([mex_extension],mexsol)
         ;;
         *)
           AC_MSG_ERROR([m4_text_wrap([Unknown target architecture...])])
         ;;
       esac])
AC_SUBST(mex_extension)

# check for MATLAB install directory
matlab_install_dir=""  # initialize matlab_install_dir
AC_ARG_WITH([matlab_install_dir],
            [AC_HELP_STRING([--with-matlab-install-dir=DIR],
            [Specify the location where MATLAB files should be installed (default [NO])])], 
            [], [withval=no])
AS_IF([! test -z "$matlab_dir"],
      [AC_MSG_CHECKING([MATLAB install directory])]
       AS_IF([test $withval = no], [],
             [matlab_install_dir=$with_matlab_install_dir])
       AC_MSG_RESULT([$matlab_install_dir])
)
AC_SUBST(matlab_install_dir)


#==========================================================================
# Checks for header files.
#==========================================================================
AC_CHECK_HEADERS([float.h limits.h])


#==========================================================================
# Checks for typedefs, structures, and compiler characteristics.
#==========================================================================
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE


#--------------------------------------------------------------------------
# Determine if allow implicit-template-instantiation
#--------------------------------------------------------------------------
AC_MSG_CHECKING([whether to use implicit-template-instantiation])
AC_ARG_ENABLE([implicit-template-instantiation],
              [AC_HELP_STRING([--enable-implicit-template-instantiation],
              [Enable implicit template instantiation (default [NO])])],
              [], [enable_implicit_template_instantiation=no])
AC_MSG_RESULT([$enable_implicit_template_instantiation])
AS_IF([test $enable_implicit_template_instantiation = no],
      [case "$CXX" in
         g++)
           CXXFLAGS="-fno-implicit-templates $CXXFLAGS"
         ;;
         xlC)
           CXXFLAGS="-qnotempinc $CXXFLAGS"
         ;;
         icpc)
           :
         ;;
         cc)
           CXXFLAGS="-instances=explicit $CXXFLAGS"
        ;;
      esac])


#==========================================================================
# Process remaining configure options.
#==========================================================================
AC_MSG_CHECKING([whether to build single-precision library])
AC_ARG_ENABLE([float],
              [AC_HELP_STRING([--enable-float],
              [Compile LSMLIB for single-precision calculations (default NO)])],
              [], [enable_float=no])
AC_MSG_RESULT([$enable_float])
AS_IF([test $enable_float = no],
      [case "$F77" in
         g77)
           FFLAGS="-fdefault-real-8 -fdefault-double-8 $FFLAGS"
         ;;
         gfortran)
           FFLAGS="-fdefault-real-8 -fdefault-double-8 $FFLAGS"
         ;;
         ifort)
           FFLAGS="-r8 $FFLAGS"
         ;;
         xlf)
           FFLAGS="-qrealsize=8 $FFLAGS"
         ;;
         blrts_xlf)
           FFLAGS="-qrealsize=8 $FFLAGS"
         ;;
      esac])
AS_IF([test $enable_float = no],
      [AC_DEFINE(LSMLIB_DOUBLE_PRECISION)],
      [])
AS_IF([test $enable_float = no],
      [AC_DEFINE(LSMLIB_REAL,double)],
      [AC_DEFINE(LSMLIB_REAL,float)])
AS_IF([test $enable_float = no],
      [AC_DEFINE(LSMLIB_ZERO_TOL,1.e-11)],
      [AC_DEFINE(LSMLIB_ZERO_TOL,1.e-5)])
AS_IF([test $enable_float = no],
      [AC_SUBST(lsmlib_zero_tol,1.d-11)],
      [AC_SUBST(lsmlib_zero_tol,1.e-5)])
AS_IF([test $enable_float = no],
      [AC_SUBST(tiny_nonzero_number,1.d-99)],
      [AC_SUBST(tiny_nonzero_number,1.d-35)])
AS_IF([test $enable_float = no],
      [AC_DEFINE(LSMLIB_REAL_MAX,DBL_MAX)],
      [AC_DEFINE(LSMLIB_REAL_MAX,FLT_MAX)])
AS_IF([test $enable_float = no],
      [AC_DEFINE(LSMLIB_REAL_MIN,DBL_MIN)],
      [AC_DEFINE(LSMLIB_REAL_MIN,FLT_MIN)])
AS_IF([test $enable_float = no],
      [AC_DEFINE(LSMLIB_REAL_EPSILON,DBL_EPSILON)],
      [AC_DEFINE(LSMLIB_REAL_EPSILON,FLT_EPSILON)])


#--------------------------------------------------------------------------
# Determine build mode; either "optimize", "debug", or "profile".
#--------------------------------------------------------------------------
SET_BUILD_MODE

# set extra options for specified build mode
case $build_mode in
  optimize) : ;;

  debug) :

    # set up LSMLIB_config.h
    AC_DEFINE(LSMLIB_DEBUG_NO_INLINE)

  ;;

  profile)
    :
  ;;

esac

# substitute compiler flags
AC_SUBST([CFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([FFLAGS])


#==========================================================================
# Checks for library functions.
#==========================================================================
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt floor pow])


#=======================================================================
# Generate the output files.
#=======================================================================

AC_OUTPUT([config/Makefile.config
           config/Makefile.config.MATLAB
           Makefile
           src/Makefile
           src/matlab/Makefile
           src/matlab/fast_marching_method/Makefile
           src/matlab/level_set_evolution/Makefile
           src/matlab/reinitialization/Makefile
           src/matlab/spatial_derivatives/Makefile
           src/matlab/time_integration/Makefile
           src/parallel/Makefile
           src/parallel/Makefile.depend
           src/parallel/fortran/Makefile
           src/parallel/templates/Makefile
           src/parallel/templates/Makefile.depend
           src/serial/Makefile
           src/serial/Makefile.depend
           src/toolbox/Makefile
           src/toolbox/boundary_conditions/Makefile
           src/toolbox/fast_marching_method/Makefile
           src/toolbox/fast_marching_method/Makefile.depend
           src/toolbox/field_extension/Makefile
           src/toolbox/geometry/Makefile
           src/toolbox/geometry/lsm_curvature2d.f
           src/toolbox/geometry/lsm_curvature2d_local.f
           src/toolbox/geometry/lsm_curvature3d.f
           src/toolbox/geometry/lsm_curvature3d_local.f
           src/toolbox/geometry/lsm_geometry1d.f
           src/toolbox/geometry/lsm_geometry2d.f
           src/toolbox/geometry/lsm_geometry2d_local.f
           src/toolbox/geometry/lsm_geometry3d_fort.f
           src/toolbox/geometry/test/Makefile
           src/toolbox/level_set_evolution/Makefile
           src/toolbox/level_set_evolution/lsm_level_set_evolution1d.f
           src/toolbox/level_set_evolution/lsm_level_set_evolution2d.f
           src/toolbox/level_set_evolution/lsm_level_set_evolution2d_local.f
           src/toolbox/level_set_evolution/lsm_level_set_evolution3d.f
           src/toolbox/level_set_evolution/lsm_level_set_evolution3d_local.f
           src/toolbox/localization/Makefile
           src/toolbox/reinitialization/Makefile
           src/toolbox/reinitialization/lsm_reinitialization1d.f
           src/toolbox/reinitialization/lsm_reinitialization2d.f
           src/toolbox/reinitialization/lsm_reinitialization2d_local.f
           src/toolbox/reinitialization/lsm_reinitialization3d.f
           src/toolbox/reinitialization/lsm_reinitialization3d_local.f
           src/toolbox/spatial_derivatives/Makefile
           src/toolbox/spatial_derivatives/lsm_spatial_derivatives1d.f
           src/toolbox/spatial_derivatives/lsm_spatial_derivatives2d.f
           src/toolbox/spatial_derivatives/lsm_spatial_derivatives2d_local.f
           src/toolbox/spatial_derivatives/lsm_spatial_derivatives3d.f
           src/toolbox/spatial_derivatives/lsm_spatial_derivatives3d_local.f
           src/toolbox/time_integration/Makefile
           src/toolbox/utilities/Makefile
           src/toolbox/utilities/lsm_calculus_toolbox2d.f
           src/toolbox/utilities/lsm_calculus_toolbox2d_local.f
           src/toolbox/utilities/lsm_calculus_toolbox3d.f
           src/toolbox/utilities/lsm_utilities1d.f
           src/toolbox/utilities/lsm_utilities2d.f
           src/toolbox/utilities/lsm_utilities3d.f
           src/toolbox/utilities/lsm_utilities2d_local.f
           src/toolbox/utilities/lsm_utilities3d_local.f
           examples/Makefile
           examples/parallel/Makefile
           examples/parallel/2d/Makefile
           examples/parallel/2d/advection/Makefile
           examples/parallel/2d/field_extension/Makefile
           examples/parallel/2d/normal_velocity_motion/Makefile
           examples/parallel/2d/toolbox/Makefile
           examples/parallel/2d/vector_level_sets/Makefile
           examples/parallel/3d/Makefile
           examples/parallel/3d/advection/Makefile
           examples/parallel/3d/field_extension/Makefile
           examples/parallel/3d/normal_velocity_motion/Makefile
           examples/parallel/3d/orthogonalization/Makefile
           examples/parallel/3d/toolbox/Makefile
           examples/serial/Makefile
           examples/serial/curvature_example/Makefile
           examples/serial/fast_marching_method/Makefile
           examples/serial/reinitialization_example/Makefile
           examples/toolbox/Makefile
           examples/toolbox/boundary_conditions/Makefile
           examples/toolbox/fast_marching_method/Makefile])


#==========================================================================
# Make includes directory
#==========================================================================
make includes

#==========================================================================
# Copy MATLAB scripts to build directory
#==========================================================================
if test "$srcdir" != "."; then
  echo "Copying MATLAB scripts to build directory..."

  MATLAB_SRC_DIRS="fast_marching_method level_set_evolution 
                   reinitialization spatial_derivatives time_integration"
  for dir in ${MATLAB_SRC_DIRS}; do
    cp -f $srcdir/src/matlab/$dir/*.m src/matlab/$dir
  done
fi

#==========================================================================
# Copy example code to build directory
#==========================================================================
if test "$srcdir" != "."; then
  echo "Copying examples to build directory..."
  cp -r $srcdir/examples/* examples/
  configure_input_files=`find examples -name "*.in"`
  for file in $configure_input_files; do rm -f $file; done
  svn_directories=`find examples -name ".svn"`
  for dir in $svn_directories; do rm -rf $dir; done
fi

#==========================================================================
# Print out configuration result
#==========================================================================
echo""
echo "  Configuration Summmary"
echo "  Compiling Options:"
echo "              Compilation Mode: $build_mode"
if test "$enable_float" = "no"; then
echo "      Floating-Point Precision: double"
else 
echo "      Floating-Point Precision: single"
fi
echo "                C Compiler(CC): $CC"
echo "                        CFLAGS: $CFLAGS"
echo "             C++ Compiler(CXX): $CXX"
echo "                      CXXFLAGS: $CXXFLAGS"
echo "                       LDFLAGS: $LDFLAGS"
echo "      Fortran 77 Compiler(F77): $F77"
echo "                        FFLAGS: $FFLAGS"
if test -z "$samrai_dir"; then
echo "                        SAMRAI: not configured"
echo "                           MPI: not needed"
else
echo "           SAMRAI include path: $samrai_inc"
echo "           SAMRAI library path: $samrai_lib"
echo "                       MPI LIB: $mpi_type"
echo "              MPI include path: $mpi_include_dirs"
echo "              MPI library path: $mpi_lib_dir"
echo "                 MPI libraries: $mpi_libs"
fi
if test -z "$matlab_dir"; then
echo "                        MATLAB: not configured"
else
echo "                   MATLAB path: $matlab_dir"
echo "      MATLAB install directory: $matlab_install_dir"
echo "             Install directory: $prefix"
fi
echo ""
echo ""
echo "If you are happy with the above configuration, type 'make' "
echo "to compile the LSMLIB library followed by 'make install' to "
echo "install the LSMLIB library. Please send comments or bugs "
echo "to <ktchu AT serendipityresearch DOT org>"
