c***********************************************************************
c
c  File:        testlsm_2d_normalvelocityfield.f
c  Copyright:   (c) 2005-2006 Kevin T. Chu
c  Revision:    $Revision: 1.2 $
c  Modified:    $Date: 2006/01/24 21:45:59 $
c  Description: F77 normal velocity field routines for 2d LSM test problem
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
     &  dx,
     &  x_lower,
     &  speed,
     &  omega,
     &  time)
c***********************************************************************
      implicit none

      integer ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      double precision vel_n(ilo_vel_gb:ihi_vel_gb,
     &                       jlo_vel_gb:jhi_vel_gb)
      integer bdry_loc
      integer i,j
      double precision dx(0:1)
      double precision x_lower(0:1)
      double precision time
      double precision speed
      double precision omega
      double precision x,y
      double precision r
      double precision speed_n

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
     &  x_lower, y_lower,
     &  speed)
c***********************************************************************
      implicit none

      integer ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      double precision vel_n(ilo_vel_gb:ihi_vel_gb,
     &                       jlo_vel_gb:jhi_vel_gb)
      double precision phi(ilo_phi_gb:ihi_phi_gb,
     &                     jlo_phi_gb:jhi_phi_gb)
      integer bdry_loc
      integer i,j
      double precision dx, dy
      double precision x_lower, y_lower
      double precision speed
      double precision x,y
      double precision r
      double precision kappa
      double precision phi_x
      double precision phi_y
      double precision phi_xx
      double precision phi_xy
      double precision phi_yy
      double precision norm_grad_phi

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
