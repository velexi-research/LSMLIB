c***********************************************************************
c
c  File:        velocityfield_fort.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 velocity field routines for 3d LSM example problem
c
c***********************************************************************
c***********************************************************************
c  Uniform velocity in x-direction with magnitude 1:  U = (1,0,0)
c***********************************************************************
      subroutine uniformvelx(
     &  u,v,w,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real u(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real v(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real w(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      integer i,j,k
      real zero,one
      parameter (zero=0.0)
      parameter (one=1.0)

c     loop over box {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb
  
            u(i,j,k) = one
            v(i,j,k) = zero
            w(i,j,k) = zero

          enddo
        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
c***********************************************************************
c  Uniform velocity in y-direction with magnitude 1:  U = (0,1,0)
c***********************************************************************
      subroutine uniformvely(
     &  u,v,w,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real u(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real v(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real w(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      integer i,j,k
      real zero,one
      parameter (zero=0.0)
      parameter (one=1.0)

c     loop over box {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb
  
            u(i,j,k) = zero
            v(i,j,k) = one
            w(i,j,k) = zero

          enddo
        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
c***********************************************************************
c  Uniform velocity in (1,1,0)-direction with magnitude sqrt(2):  
c    U = (1,1,0)
c***********************************************************************
      subroutine uniformvelxy(
     &  u,v,w,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real u(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real v(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real w(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      integer i,j,k
      real one
      parameter (one=1.0)

c     loop over box {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb
  
            u(i,j,k) = one
            v(i,j,k) = one
            w(i,j,k) = 0.0

          enddo
        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
c***********************************************************************
c Pure rotation velocity field with angular velocity 1: 
c   U = (-y,x,0)
c***********************************************************************
      subroutine rotatingvel(
     &  u,v,w,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dx,
     &  x_lower)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real u(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real v(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real w(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      integer i,j,k
      real dx(0:2)
      real x_lower(0:2)
      real x,y,z

c     loop over box {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb
 
            x = x_lower(0) + dx(0)*(0.5+i-ilo_gb)
            y = x_lower(1) + dx(1)*(0.5+j-jlo_gb)
            z = x_lower(2) + dx(2)*(0.5+k-klo_gb)
            u(i,j,k) = -y
            v(i,j,k) = x
            w(i,j,k) = 0.0

          enddo
        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
c***********************************************************************
c Pure expansion/compression velocity field oscillating in time:
c   U = speed*cos(omega*t)
c***********************************************************************
      subroutine expandingvel(
     &  u,v,w,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dx,
     &  x_lower,
     &  speed,
     &  omega,
     &  time)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real u(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real v(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      real w(ilo_gb:ihi_gb,
     &       jlo_gb:jhi_gb,
     &       klo_gb:khi_gb)
      integer i,j,k
      real dx(0:2)
      real x_lower(0:2)
      real time
      real speed
      real omega
      real x,y,z
      real r

c     loop over box {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb
 
            x = x_lower(0) + dx(0)*(0.5+i-ilo_fb)
            y = x_lower(1) + dx(1)*(0.5+j-jlo_fb)
            z = x_lower(2) + dx(2)*(0.5+k-klo_fb)
            r = sqrt(x**2 + y**2 + z**2)
            if (r .ne. 0) then
              u(i,j,k) = speed*cos(omega*time)*x/r
              v(i,j,k) = speed*cos(omega*time)*y/r
              w(i,j,k) = speed*cos(omega*time)*z/r
            else
              u(i,j,k) = 0.0
              v(i,j,k) = 0.0
              w(i,j,k) = 0.0
            endif

          enddo
        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
