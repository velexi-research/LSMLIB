#
## File:        check_mpi.m4
## Copyright:   (c) 2005-2006 Princeton University
## Authors:     Zhaoxuan Wu and Kevin T. Chu
## Revision:    $Revision: 1.3 $
## Modified:    $Date: 2007/05/08 19:16:50 $
## Description: autoconf macro for checking mpi 
##

#
# check for mpi library
#
AC_DEFUN(CHECK_MPI,[

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
    mpi_libs="-lmpi_f90 -lmpi_cxx -lmpi -lorte -lopal"
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
