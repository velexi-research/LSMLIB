c***********************************************************************
c
c  File:        testlsm_3d_patchmodule.f
c  Copyright:   (c) 2005-2006 Kevin T. Chu
c  Revision:    $Revision: 1.2 $
c  Modified:    $Date: 2006/01/24 21:46:09 $
c  Description: F77 patch module routines for 3d LSM test problem
c
c***********************************************************************
c***********************************************************************
      subroutine initsphere(
     &  level_set,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb, 
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb, 
     &  x_lower, 
     &  dx,
     &  center,
     &  radius)
c***********************************************************************
      implicit none

c     _gb refers to ghost box
c     _fb refers to fill box
      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      double precision level_set(ilo_gb:ihi_gb,
     &                           jlo_gb:jhi_gb,
     &                           klo_gb:khi_gb)
      double precision x_lower(0:2)
      double precision dx(0:2)
      double precision x,y,z
      integer i,j,k
      double precision radius
      double precision center(0:2)

c     loop over grid {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            x = x_lower(0) + dx(0)*(i-ilo_fb+0.5)
            y = x_lower(1) + dx(1)*(j-jlo_fb+0.5)
            z = x_lower(2) + dx(2)*(k-klo_fb+0.5)
  
            level_set(i,j,k) = sqrt( (x-center(0))**2 
     &                             + (y-center(1))**2
     &                             + (z-center(2))**2) 
     &                       - radius
  
          enddo
        enddo
      enddo
c     } end loop over grid 

      return
      end
c***********************************************************************
