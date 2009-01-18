c***********************************************************************
c
c  File:        velocityfield.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 normal velocity field routines for 2d LSM example 
c               problem
c
c***********************************************************************

c***********************************************************************
c Pure expansion/compression velocity field oscillating in time:
c   V_n = speed*cos(omega*t)
c***********************************************************************
      subroutine oscillatingExpansionVelocity (
     &  vel_n,
     &  ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  x_lower,
     &  speed,
     &  omega,
     &  time)
c***********************************************************************
      implicit none

      integer ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real vel_n(ilo_vel_gb:ihi_vel_gb, jlo_vel_gb:jhi_vel_gb)
      integer i,j
      real x_lower(0:1)
      real time
      real speed
      real omega
      real speed_n

      speed_n = speed*cos(omega*time)

c     loop over box {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          vel_n(i,j) = speed_n

        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************

c***********************************************************************
c Motion due to curvature; velocity field proportional to minus the
c mean curvature.
c   V_n = -speed*kappa
c***********************************************************************
      subroutine meanCurvatureVelocity (
     &  vel_n,
     &  ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy,
     &  speed)
c***********************************************************************
      implicit none

      integer ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real vel_n(ilo_vel_gb:ihi_vel_gb, jlo_vel_gb:jhi_vel_gb)
      real phi(ilo_phi_gb:ihi_phi_gb, jlo_phi_gb:jhi_phi_gb)
      integer i,j
      real dx, dy
      real speed
      real kappa
      real phi_x
      real phi_y
      real phi_xx
      real phi_xy
      real phi_yy
      real norm_grad_phi

c     loop over box {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb
 
c          approximate mean curvature 
c          (only valid for signed distance functions)
c          kappa = ( phi(i+1,j)+phi(i-1,j)-2*phi(i,j) )/dx/dx
c     &          + ( phi(i,j+1)+phi(i,j-1)-2*phi(i,j) )/dy/dy

          phi_x = (phi(i+1,j)-phi(i-1,j))/dx
          phi_y = (phi(i,j+1)-phi(i,j-1))/dy
          phi_xx = (phi(i+1,j)+phi(i-1,j)-2*phi(i,j))/dx/dx
          phi_yy = (phi(i,j+1)+phi(i,j-1)-2*phi(i,j))/dy/dy
          phi_xy = ( (phi(i+1,j+1)-phi(i-1,j+1))
     &             - (phi(i+1,j-1)-phi(i-1,j-1)) )/dx/dy
          norm_grad_phi = sqrt(phi_x*phi_x + phi_y*phi_y)

          kappa = ( phi_x*phi_x*phi_yy + phi_y*phi_y*phi_xx
     &            - 2*phi_x*phi_y*phi_xy )
     &          / norm_grad_phi/norm_grad_phi/norm_grad_phi
          vel_n(i,j) = -speed*kappa

        enddo
      enddo
c     } end loop over box

      return
      end
c***********************************************************************
