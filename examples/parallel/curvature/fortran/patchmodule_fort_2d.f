c***********************************************************************
c
c  File:        patchmodule_fort.f
c  Description: F77 patch module routines for 2d LSM example problem
c
c***********************************************************************
c***********************************************************************
      subroutine initcircle(
     &  phi,
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
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real x_lower(0:1)
      real dx(0:1)
      real x,y
      integer i,j
      real center(0:1)
      real radius

c     loop over grid {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          x = x_lower(0) + dx(0)*(i-ilo_fb+0.5)
          y = x_lower(1) + dx(1)*(j-jlo_fb+0.5)

          phi(i,j) = sqrt((x-center(0))**2 + (y-center(1))**2) - radius

        enddo
      enddo
c     } end loop over grid

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine initstar(
     &  phi,
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
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real x_lower(0:1)
      real dx(0:1)
      real x,y
      integer i,j
      real center(0:1)
      real radius
      real theta

c     loop over grid {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          x = x_lower(0) + dx(0)*(i-ilo_fb+0.5) - center(0)
          y = x_lower(1) + dx(1)*(j-jlo_fb+0.5) - center(1)

          theta = atan2(y,x)
          phi(i,j) = sqrt(x**2+y**2) - radius*(1+0.5*sin(5*theta));

        enddo
      enddo
c     } end loop over grid

      return
      end
c***********************************************************************
