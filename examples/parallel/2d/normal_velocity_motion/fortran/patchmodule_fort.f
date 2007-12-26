c***********************************************************************
c
c  File:        patchmodule_fort.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision: 1.2 $
c  Modified:    $Date: 2006/01/24 21:45:59 $
c  Description: F77 patch module routines for 2d LSM example problem
c
c***********************************************************************

c***********************************************************************
      subroutine initcircle(
     &  level_set,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  x_lower, 
     &  dx,
     &  center,
     &  radius)
c***********************************************************************
      implicit none

c     _gb refers to ghost box
c     _fb refers to fill box
      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real level_set(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real x_lower(0:1)
      real dx(0:1)
      real x,y
      integer i,j
      real radius
      real center(0:1)

c     loop over grid {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          x = x_lower(0) + dx(0)*(i-ilo_fb+0.5)
          y = x_lower(1) + dx(1)*(j-jlo_fb+0.5)

          level_set(i,j) = sqrt((x-center(0))**2 + (y-center(1))**2) 
     &                   - radius

        enddo
      enddo
c     } end loop over grid 

      return
      end
c***********************************************************************

c***********************************************************************
      subroutine initlobes(
     &  level_set,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  x_lower, 
     &  dx,
     &  center,
     &  radius,
     &  num_lobes)
c***********************************************************************
      implicit none

c     _gb refers to ghost box
c     _fb refers to fill box
      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real level_set(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real x_lower(0:1)
      real dx(0:1)
      real x,y
      real r, theta
      integer i,j
      real radius
      real center(0:1)
      integer num_lobes
      real PI
      parameter (PI=3.14159265358979d0)

c     loop over grid {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          x = x_lower(0) + dx(0)*(i-ilo_fb+0.5)
          y = x_lower(1) + dx(1)*(j-jlo_fb+0.5)

          r = sqrt((x-center(0))**2 + (y-center(1))**2) 
          theta = atan( (y-center(1))/(x-center(0)) )
          if ( (x-center(0)) .lt. 0.d0 ) then
            theta = theta + PI
          endif

          level_set(i,j) = r 
     &                   - radius * ( cos(num_lobes*theta) + 2 ) /3.d0

        enddo
      enddo
c     } end loop over grid 

      return
      end
c***********************************************************************
