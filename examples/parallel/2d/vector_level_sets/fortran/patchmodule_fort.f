c***********************************************************************
c
c  File:        patchmodule_fort.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision: 1.2 $
c  Modified:    $Date: 2006/01/24 21:46:05 $
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
