c***********************************************************************
c
c  File:        testlsm_2d_velocityfield.f
c  Copyright:   (c) 2005-2006 Kevin T. Chu
c  Revision:    $Revision: 1.2 $
c  Modified:    $Date: 2006/01/24 21:46:02 $
c  Description: F77 velocity field routines for 2d LSM test problem
c
c***********************************************************************
c***********************************************************************
c  Uniform velocity in x-direction with magnitude 1:  U = (1,0)
c***********************************************************************
      subroutine uniformvelx(
     &  u,v,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      double precision u(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      double precision v(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      integer i,j
      double precision zero,one
      parameter (zero=0.0)
      parameter (one=1.0)

c     loop over box {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb
  
          u(i,j) = one
          v(i,j) = zero

        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
c***********************************************************************
c  Uniform velocity in y-direction with magnitude 1:  U = (0,1)
c***********************************************************************
      subroutine uniformvely(
     &  u,v,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      double precision u(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      double precision v(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      integer i,j
      double precision zero,one
      parameter (zero=0.0)
      parameter (one=1.0)

c     loop over box {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb
  
          u(i,j) = zero
          v(i,j) = one

        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
c***********************************************************************
c  Uniform velocity in (1,1)-direction with magnitude sqrt(2):  
c    U = (1,1)
c***********************************************************************
      subroutine uniformvelxy(
     &  u,v,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      double precision u(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      double precision v(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      integer i,j
      double precision one
      parameter (one=1.0)

c     loop over box {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb
  
          u(i,j) = one
          v(i,j) = one

        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
c***********************************************************************
c Pure rotation velocity field with angular velocity 1: 
c   U = (-y,x)
c***********************************************************************
      subroutine rotatingvel(
     &  u,v,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx,
     &  x_lower)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      double precision u(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      double precision v(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      integer i,j
      double precision dx(0:1)
      double precision x_lower(0:1)
      double precision x,y

c     loop over box {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb
 
          x = x_lower(0) + dx(0)*(0.5+i-ilo_gb)
          y = x_lower(1) + dx(1)*(0.5+j-jlo_gb)
          u(i,j) = -y
          v(i,j) = x

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
     &  u,v,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx,
     &  x_lower,
     &  speed,
     &  omega,
     &  time)
c***********************************************************************
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      double precision u(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      double precision v(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      integer i,j
      double precision dx(0:1)
      double precision x_lower(0:1)
      double precision time
      double precision speed
      double precision omega
      double precision x,y
      double precision r

c     loop over box {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb
 
          x = x_lower(0) + dx(0)*(0.5+i-ilo_fb)
          y = x_lower(1) + dx(1)*(0.5+j-jlo_fb)
          r = sqrt(x**2 + y**2)
          if (r .ne. 0) then
            u(i,j) = speed*cos(omega*time)*x/r
            v(i,j) = speed*cos(omega*time)*y/r
          else
            u(i,j) = 0.0
            v(i,j) = 0.0
          endif

        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
