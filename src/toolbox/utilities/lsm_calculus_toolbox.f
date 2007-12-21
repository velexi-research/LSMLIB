c***********************************************************************
c
c  File:        lsm_calculus_toolbox.f
c  Copyright:   (c) 2005-2006 Kevin T. Chu
c  Revision:    $Revision: 1.8 $
c  Modified:    $Date: 2006/04/20 12:22:18 $
c  Description: F77 routines for several common level set method
c               calculus calculations
c
c***********************************************************************

c***********************************************************************
c
c  lsmHeaviside() computes the value of the standard smoothed Heaviside 
c  function for level set method calculations
c
c  Arguments:
c    H (out):       value of Heaviside function at specified position
c    x (in):        spatial position Heaviside function evaluated at
c    epsilon (in):  width of numerical smoothing
c
c***********************************************************************
      subroutine lsmHeaviside(
     &  H,
     &  x,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real H
      real x
      real epsilon
      real one_over_epsilon
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)

      if (x .lt. -epsilon) then
        H = 0.0d0
      elseif (x .gt. epsilon) then
        H = 1.0d0
      else
        one_over_epsilon = 1.d0/epsilon
        H = 0.5d0*( 1.d0+x*one_over_epsilon
     &                  +one_over_pi*sin(pi*x*one_over_epsilon) )
      endif

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsmDeltaFunction() computes the value of the standard smoothed delta 
c  function for level set method calculations
c
c  Arguments:
c    delta (out):   value of delta-function at specified position
c    x (in):        spatial position delta-function evaluated at
c    epsilon (in):  width of numerical smoothing
c
c***********************************************************************
      subroutine lsmDeltaFunction(
     &  delta,
     &  x,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real delta
      real x
      real epsilon
      real one_over_epsilon
      real pi
      parameter (pi=3.14159265358979323846d0)

      if (abs(x) .gt. epsilon) then
        delta = 0.0d0
      else
        one_over_epsilon = 1.d0/epsilon
        delta = 0.5d0*one_over_epsilon
     &               *( 1.d0+cos(pi*x*one_over_epsilon) )
      endif

      return
      end
c } end subroutine
c***********************************************************************
