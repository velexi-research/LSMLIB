##
## File:        set_build_mode.m4
## Copyright:   (c) 2005-2006 Princeton University
## Authors:     Zhaoxuan Wu and Kevin T. Chu
## Revision:    $Revision: 1.7 $
## Modified:    $Date: 2007/05/08 19:16:50 $
## Description: autoconf macro for setting the build mode 
##

#
# set the build mode
#
AC_DEFUN(SET_BUILD_MODE,[

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
AS_IF([test -z "$F77"], [F77=g77])

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
        CFLAGS="-O3 $CFLAGS"
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
        CXXFLAGS="-O3 $CXXFLAGS"
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
