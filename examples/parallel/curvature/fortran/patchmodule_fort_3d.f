c***********************************************************************
c
c  File:        patchmodule_fort_3d.f
c  Description: F77 patch module routines for 3D LSM example problem
c
c***********************************************************************
c***********************************************************************
      subroutine initsphere(
     &  phi,
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
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      real x_lower(0:2)
      real dx(0:2)
      real x,y,z
      integer i,j,k
      real center(0:2)
      real radius

c     loop over grid {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            x = x_lower(0) + dx(0)*(i-ilo_fb+0.5) - center(0)
            y = x_lower(1) + dx(1)*(j-jlo_fb+0.5) - center(1)
            z = x_lower(2) + dx(2)*(k-klo_fb+0.5) - center(2)

            phi(i,j,k) = sqrt(x**2 + y**2 + z**2) - radius

          enddo
        enddo
      enddo
c     } end loop over grid

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine initbumpysphere(
     &  phi,
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
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      real x_lower(0:2)
      real dx(0:2)
      real x,y,z,rho
      integer i,j,k
      real center(0:2)
      real radius
      real theta, azimuth

c     loop over grid {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            x = x_lower(0) + dx(0)*(i-ilo_fb+0.5) - center(0)
            y = x_lower(1) + dx(1)*(j-jlo_fb+0.5) - center(1)
            z = x_lower(2) + dx(2)*(k-klo_fb+0.5) - center(2)

            theta = atan2(y,x)
            rho = sqrt(x**2 + y**2 + z**2)
            if (rho .gt. 1e-7) then
                azimuth = acos(z / rho)
            else
                azimuth = 0.0
            endif

            phi(i,j,k) = rho
     &          - radius*(1 + 0.25 * sin(4*azimuth) * cos(5*theta))

          enddo
        enddo
      enddo
c     } end loop over grid

      return
      end
c***********************************************************************
