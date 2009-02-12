# generated automatically by aclocal 1.10 -*- Autoconf -*-

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004,
# 2005, 2006  Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.


#
# check for mpi library
#
AC_DEFUN([CHECK_MPI],[

PACKAGE_NAME_LOCAL_VAR="$1"

AC_MSG_CHECKING([mpi library]) 

# determine which MPI library to use
mpi_type="" # initialize mpi_type
AC_ARG_WITH([lammpi], [AC_HELP_STRING([--with-lammpi=DIR],
            [Specify the lammpi library path (default [NO])])],
            [], [withval=no])
AS_IF([test $withval != no], [mpi_type=lammpi])

AC_ARG_WITH([mpich], [AC_HELP_STRING([--with-mpich=DIR],
            [Specify the mpich library path (default [NO])])],
            [], [withval=no])
AS_IF([test $withval != no], [mpi_type=mpich])

AC_ARG_WITH([openmpi], [AC_HELP_STRING([--with-openmpi=DIR],
            [Specify the openmpi library path (default [NO])])],
            [], [withval=no])
AS_IF([test $withval != no], [mpi_type=openmpi])

AC_ARG_WITH([ibmmpi], [AC_HELP_STRING([--with-ibmmpi=DIR],
            [Specify the ibmmpi library path (default [NO])])],
            [], [withval=no])
AS_IF([test $withval != no], [mpi_type=ibmmpi])

AC_ARG_WITH([intelmpi], [AC_HELP_STRING([--with-intelmpi=DIR],
                [Specify the intelmpi library path (default [NO])])],
        [], [withval=no])
AS_IF([test $withval != no],
        [mpi_type=intelmpi])

AC_ARG_WITH([sgimpi], [AC_HELP_STRING([--with-sgimpi=DIR],
            [Specify the sgimpi library path (default [NO])])],
            [], [withval=no])
AS_IF([test $withval != no], [mpi_type=sgimpi])

# output result to stdout
AS_IF([test -z "$mpi_type"], 
      AC_MSG_ERROR([m4_text_wrap([$PACKAGE_NAME_LOCAL_VAR requires MPI])]),
      AC_MSG_RESULT($mpi_type)
     )

# set MPI variables
case $mpi_type in
  lammpi)
    AC_PROGRAM_PATH(MPICC,mpicc)
    AC_PROGRAM_PATH(MPICXX,mpiCC)
    AC_PROGRAM_PATH(MPIF77,mpif77)
    AC_PROGRAM_PATH(MPIF90,mpif90)
    MPILIBNAME="lammpi"
    MPIBOOT="lamboot"
    MPIUNBOOT="wipe"
    MPIRUN="mpirun"
    mpi_dir=$with_lammpi
    mpi_include_dirs=$with_lammpi/include
    mpi_lib_dir=$with_lammpi/lib
    mpi_libs="-lmpi -llam -llammpi++ -llamf77mpi"
    ;;
  mpich)
    AC_PROGRAM_PATH(MPICC,mpicc)
    AC_PROGRAM_PATH(MPICXX,mpiCC)
    AC_PROGRAM_PATH(MPIF77,mpif77)
    AC_PROGRAM_PATH(MPIF90,mpif90)
    MPILIBNAME="mpich"
    AC_PATH_PROG(MPIRUN,mpirun)
    AC_PATH_PROG(MPIBOOT,mpichboot)
    AC_PATH_PROG(MPIUNBOOT,mpichstop)
    mpi_dir=$with_mpich
    mpi_include_dirs=$with_mpich/include
    mpi_lib_dir=$with_mpich/lib
    mpi_libs="-lmpich -lpmpich"
    ;;
  openmpi)
    AC_PROGRAM_PATH(MPICC,mpicc)
    AC_PROGRAM_PATH(MPICXX,mpiCC)
    AC_PROGRAM_PATH(MPIF77,mpif77)
    AC_PROGRAM_PATH(MPIF90,mpif90)
    MPILIBNAME="OpenMPI"
    AC_PATH_PROG(MPIEXEC,mpiexec)
    mpi_dir=$with_openmpi
    mpi_include_dirs="$with_openmpi/include $with_openmpi/include/openmpi/ompi"
    mpi_lib_dir=$with_openmpi/lib
    mpi_libs="-lmpi_f90 -lmpi_cxx -lmpi"
    ;;
  ibmmpi)
    AC_PROGRAM_PATH(MPICC,mpcc)
    AC_PROGRAM_PATH(MPIF77,mpxlf)
    mpi_dir=$with_ibmmpi
    mpi_include_dirs=$with_ibmmpi/include
    mpi_lib_dir=$with_ibmmpi/lib
    mpi_libs="??"
    ;;
  intelmpi)
    mpi_dir="$with_intelmpi"
    mpi_include_dirs="-I$with_intelmpi/include"
    mpi_lib_dir="-L$with_intelmpi/lib" 
    mpi_libs="-lmpi -lmpi++"
    ;;
  sgimpi)
    mpi_dir=$with_sgimpi
    mpi_include_dirs=$with_sgimpi/include
    mpi_lib_dir=$with_sgimpi/lib
    mpi_libs="??"
    ;;
  *) AC_MSG_ERROR([m4_text_wrap(
      [unrecognized mpi library --with-$mpi_type;
      use `lammpi', `mpich', `ibmmpi' `intelmpi' or `sgimpi'], [    ], [[]], [60])]) ;;
esac

AC_SUBST([mpi_dir])
AC_SUBST([mpi_include_dirs])
AC_SUBST([mpi_lib_dir])
AC_SUBST([mpi_libs])

]) # end of CHECK_MPI function


#
# set the build mode
#
AC_DEFUN([SET_BUILD_MODE],[

AC_MSG_CHECKING([build mode])

AC_ARG_ENABLE([opt], [AC_HELP_STRING([--enable-opt],
              [build with optimizations enabled (default [YES]);
               same as --enable-mode=opt])],
              [enable_opt=yes], [enable_opt=no])

AC_ARG_ENABLE([debug], [AC_HELP_STRING([--enable-debug],
              [build with debugging information (default [NO]);
               same as --enable-mode=debug])],
              [enable_debug=yes], [enable_debug=no])

AC_ARG_ENABLE([profile], [AC_HELP_STRING([--enable-profile],
              [build with profiling information (default [NO]);
               same as --enable-mode=profile])],
              [enable_profile=yes], [enable_profile=no])

AC_ARG_ENABLE([mode], [AC_HELP_STRING([--enable-mode=mode],
              [set build mode; recognized modes are `opt', `debug', `profile' (default [OPTIMIZE])])],
              [case $enableval in
                 opt|debug|profile) build_mode=$enableval ;;
                 *) AC_MSG_ERROR([m4_text_wrap([unrecognized mode --enable-mode=$enableval; use `opt', `debug', or `profile'], [    ], [[]], [60])]) ;;
               esac])

# build mode priorities:  profile, debug, optimize
# default build mode:     optimize
AS_IF([test $enable_profile = yes],
      [build_mode=profile],
      [AS_IF([test $enable_debug = yes],
             [build_mode=debug],
             [AS_IF([test $enable_opt = yes],
                    [build_mode=optimize],
                    [build_mode=optimize])
            ])
      ])

# Set default compilers to avoid errors
AS_IF([test -z "$CC"], [CC=gcc])
AS_IF([test -z "$CXX"], [CXX=g++])
AS_IF([test -z "$F77"], [F77=gfortran])

# set compiler flags to build position-independent code
case "$CC" in
  gcc)
    CFLAGS="-fPIC $CFLAGS"
  ;;
  xlc)
    CFLAGS="-fPIC $CFLAGS"
  ;;
  icc)
    CFLAGS="-fPIC $CFLAGS"
  ;;
  blrts_xlc)
    # no position-independent code for BlueGene/L
  ;;
  *)
    # do not use compiler flags for position-independent code 
esac 

case "$CXX" in
  g++)
    CXXFLAGS="-fPIC $CXXFLAGS"
  ;;
  xlC)
    CXXFLAGS="-fPIC $CXXFLAGS"
  ;;
  icpc)
    CXXFLAGS="-fPIC $CXXFLAGS"
  ;;
  blrts_xlC|blrts_xlc++)
    # no position-independent code for BlueGene/L
  ;;
  *)
    # do not use compiler flags for position-independent code 
esac 

case "$F77" in
  g77)
    FFLAGS="-fPIC $FFLAGS"
  ;;
  gfortran)
    FFLAGS="-fPIC $FFLAGS"
  ;;
  xlf)
    FFLAGS="-fPIC $FFLAGS"
  ;;
  ifort)
    FFLAGS="-fPIC $FFLAGS"
  ;;
  blrts_xlf)
    # no position-independent code for BlueGene/L
  ;;
  *)
    # do not use compiler flags for position-independent code 
esac 


case $build_mode in
  optimize)
    case "$CC" in
      gcc)
        CFLAGS="-O3 -funroll-loops $CFLAGS"
      ;;
      xlc)
        CFLAGS="-O $CFLAGS"
      ;;
      icc)
        CFLAGS="-O3 $CFLAGS"
      ;;
      blrts_xlc)
        CFLAGS="-O3 -qmaxmem=64000 -qarch=440 $CFLAGS"
      ;;
      *)
        CFLAGS="-O $CFLAGS"
        AC_MSG_WARN([m4_text_wrap([unrecognized C compiler...trying '-O' compiler flag...])])
    esac 

    case "$CXX" in
      g++)
        CXXFLAGS="-O3 -funroll-loops $CXXFLAGS"
      ;;
      xlC)
        CXXFLAGS="-O $CXXFLAGS"
      ;;
      icpc)
        CXXFLAGS="-O3 $CXXFLAGS"
      ;;
      blrts_xlC|blrts_xlc++)
        CXXFLAGS="-O3 -qmaxmem=64000 -qarch=440 $CXXFLAGS"
      ;;
      *)
        CXXFLAGS="-O $CXXFLAGS"
        AC_MSG_WARN([m4_text_wrap([unrecognized C++ compiler...trying '-O' compiler flag...])])
    esac 

    case "$F77" in
      g77)
        FFLAGS="-O3 $FFLAGS"
      ;;
      gfortran)
        FFLAGS="-O3 -funroll-loops $FFLAGS"
      ;;
      xlf)
        FFLAGS="-O $FFLAGS"
      ;;
      ifort)
        FFLAGS="-O3 $FFLAGS"
      ;;
      blrts_xlf)
        FFLAGS="-O3 -qmaxmem=64000 -qarch=440 $FFLAGS"
      ;;
      *)
        FFLAGS="-O $FFLAGS"
        AC_MSG_WARN([m4_text_wrap([unrecognized Fortran compiler...trying '-O' compiler flag...])])
    esac 
  ;;

  debug)
    case "$CC" in
      gcc)
        CFLAGS="-g -Wall $CFLAGS"
      ;;
      xlc)
        CFLAGS="-g $CFLAGS"
      ;;
      icc)
        CFLAGS="-g $CFLAGS"
      ;;
      blrts_xlc)
        CFLAGS="-g $CFLAGS"
      ;;
      *)
        CFLAGS="-g $CFLAGS"
        AC_MSG_WARN([m4_text_wrap([unrecognized C compiler...trying '-g' compiler flag...])])
    esac 

    case "$CXX" in
      g++)
        CXXFLAGS="-g -Wall $CXXFLAGS"
      ;;
      xlC)
        CXXFLAGS="-g $CXXFLAGS"
      ;;
      icpc)
        CXXFLAGS="-g $CXXFLAGS"
      ;;
      blrts_xlC|blrts_xlc++)
        CXXFLAGS="-g $CXXFLAGS"
      ;;
      *)
        CXXFLAGS="-g $CXXFLAGS"
        AC_MSG_WARN([m4_text_wrap([unrecognized C++ compiler...trying '-g' compiler flag...])])
    esac 

    case "$F77" in
      g77)
        FFLAGS="-g -Wall $FFLAGS"
      ;;
      gfortran)
        FFLAGS="-g -Wall $FFLAGS"
      ;;
      xlf)
        FFLAGS="-g $FFLAGS"
      ;;
      ifort)
        FFLAGS="-g $FFLAGS"
      ;;
      blrts_xlf)
        FFLAGS="-g $FFLAGS"
      ;;
      *)
        FFLAGS="-g $FFLAGS"
        AC_MSG_WARN([m4_text_wrap([unrecognized Fortran compiler...trying '-g' compiler flag...])])
    esac 
  ;;

  profile)
    case "$CC" in
      gcc)
        CFLAGS="-pg $CFLAGS"
      ;;
      xlc)
        CFLAGS="-pg $CFLAGS"
      ;;
      icc)
        CFLAGS="-p $CFLAGS"
      ;;
      blrts_xlc)
        # no profiling for BlueGene/L
      ;;
      *)
        CFLAGS="-pg $CFLAGS"
        AC_MSG_WARN([m4_text_wrap([unrecognized C compiler...trying '-pg' compiler flag...])])
    esac 

    case "$CXX" in
      g++)
        CXXFLAGS="-pg $CXXFLAGS"
      ;;
      xlC)
        CXXFLAGS="-pg $CXXFLAGS"
      ;;
      icpc)
        CXXFLAGS="-p $CXXFLAGS"
      ;;
      blrts_xlC|blrts_xlc++)
        # no profiling for BlueGene/L
      ;;
      *)
        CXXFLAGS="-pg $CXXFLAGS"
        AC_MSG_WARN([m4_text_wrap([unrecognized C++ compiler...trying '-pg' compiler flag...])])
    esac 

    case "$F77" in
      g77)
        FFLAGS="-pg $FFLAGS"
      ;;
      gfortran)
        FFLAGS="-pg $FFLAGS"
      ;;
      xlf)
        FFLAGS="-pg $FFLAGS"
      ;;
      ifort)
        FFLAGS="-p $FFLAGS"
      ;;
      blrts_xlf)
        # no profiling for BlueGene/L
      ;;
      *)
        FFLAGS="-pg $FFLAGS"
        AC_MSG_WARN([m4_text_wrap([unrecognized Fortran compiler...trying '-pg' compiler flag...])])
    esac 
  ;;

esac

AC_MSG_RESULT([$build_mode])

]) # end of SET_BUILD_MODE function

