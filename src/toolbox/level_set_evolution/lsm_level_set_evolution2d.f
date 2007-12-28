c***********************************************************************
c
c  File:        lsm_level_set_evolution2d.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision: 1.5 $
c  Modified:    $Date: 2007/05/06 23:47:28 $
c  Description: F77 subroutines for 2D level set evolution equation
c
c***********************************************************************

c***********************************************************************
c
c  lsm2dZeroOutLevelSetEqnRHS() zeros out the right-hand side of the
c  level set equation when it is written in the form:
c
c    phi_t = ...
c
c  Arguments:
c    lse_rhs (in/out):  right-hand of level set equation
c    *_gb (in):         index range for ghostbox
c
c***********************************************************************
      subroutine lsm2dZeroOutLevelSetEqnRHS(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      integer i,j

c     { begin loop over grid
      do j=jlo_lse_rhs_gb,jhi_lse_rhs_gb
        do i=ilo_lse_rhs_gb,ihi_lse_rhs_gb
        
          lse_rhs(i,j) = 0.d0
      
        enddo 
      enddo 
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dAddAdvectionTermToLSERHS adds the contribution of an advection
c  term (external vector velocity field) to the right-hand side of the 
c  level set equation when it is written in the form:
c
c    phi_t = -vel dot grad(phi) + ...
c
c  Arguments:
c    lse_rhs (in/out):  right-hand
c    phi_* (in):        components of grad(phi) at t = t_cur
c    vel_* (in):        components of velocity at t = t_cur
c    *_gb (in):         index range for ghostbox
c    *_fb (in):         index range for fillbox
c
c***********************************************************************
      subroutine lsm2dAddAdvectionTermToLSERHS(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  vel_x, vel_y,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      integer i,j

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb
        
          lse_rhs(i,j) = lse_rhs(i,j) - ( vel_x(i,j)*phi_x(i,j)
     &                                  + vel_y(i,j)*phi_y(i,j) )
      
        enddo 
      enddo 
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dAddNormalVelTermToLSERHS() adds the contribution of a normal 
c  (scalar) velocity term to the right-hand side of the level set 
c  equation when it is written in the form:
c
c    phi_t = -V_n |grad(phi)| + ...
c
c  Arguments:
c    lse_rhs (in/out):  right-hand of level set equation
c    phi_*_plus (in):   components of forward approx to grad(phi) at 
c                       t = t_cur
c    phi_*_minus (in):  components of backward approx to grad(phi) at 
c                       t = t_cur
c    vel_n (in):        normal velocity at t = t_cur
c    *_gb (in):         index range for ghostbox
c    *_fb (in):         index range for fillbox
c
c***********************************************************************
      subroutine lsm2dAddNormalVelTermToLSERHS(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  vel_n,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real vel_n(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      integer i,j
      real vel_n_cur
      real norm_grad_phi_sq
      real zero_tol
      parameter (zero_tol=@lsmlib_zero_tol@)

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          vel_n_cur = vel_n(i,j)
          if (abs(vel_n_cur) .ge. zero_tol) then

c           { begin Godunov selection of grad_phi

            if (vel_n_cur .gt. 0.d0) then
              norm_grad_phi_sq = max(max(phi_x_minus(i,j),0.d0)**2,
     &                               min(phi_x_plus(i,j),0.d0)**2 )
     &                         + max(max(phi_y_minus(i,j),0.d0)**2,
     &                               min(phi_y_plus(i,j),0.d0)**2 )
            else
              norm_grad_phi_sq = max(min(phi_x_minus(i,j),0.d0)**2,
     &                               max(phi_x_plus(i,j),0.d0)**2 )
     &                         + max(min(phi_y_minus(i,j),0.d0)**2,
     &                               max(phi_y_plus(i,j),0.d0)**2 )
            endif

c           } end Godunov selection of grad_phi


c           compute contribution to lse_rhs(i,j) 
            lse_rhs(i,j) = lse_rhs(i,j) 
     &                   - vel_n_cur*sqrt(norm_grad_phi_sq)

          endif
     
        enddo 
      enddo 
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm2dAddConstNormalVelTermToLSERHS() adds the contribution of a normal 
c  (scalar, constant for all grid points) velocity term to the right-hand 
c  side of the level set equation when it is written in the form:
c
c    phi_t = -V_n |grad(phi)| + ...
c
c  Arguments:
c    lse_rhs (in/out):  right-hand of level set equation
c    phi_*_plus (in):   components of forward approx to grad(phi) at 
c                       t = t_cur
c    phi_*_minus (in):  components of backward approx to grad(phi) at 
c                       t = t_cur
c    vel_n (in):        normal velocity at t = t_cur, a scalar value
c    *_gb (in):         index range for ghostbox
c    *_fb (in):         index range for fillbox
c
c***********************************************************************
      subroutine lsm2dAddConstNormalVelTermToLSERHS(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  vel_n,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real vel_n
      
      integer i,j
      real norm_grad_phi_sq
      real zero_tol
      parameter (zero_tol=@lsmlib_zero_tol@)

      if (abs(vel_n) .ge. zero_tol) then

c       { begin loop over grid
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb


c           { begin Godunov selection of grad_phi

            if (vel_n .gt. 0.d0) then
              norm_grad_phi_sq = max(max(phi_x_minus(i,j),0.d0)**2,
     &                               min(phi_x_plus(i,j),0.d0)**2 )
     &                         + max(max(phi_y_minus(i,j),0.d0)**2,
     &                               min(phi_y_plus(i,j),0.d0)**2 )
            else
              norm_grad_phi_sq = max(min(phi_x_minus(i,j),0.d0)**2,
     &                               max(phi_x_plus(i,j),0.d0)**2 )
     &                         + max(min(phi_y_minus(i,j),0.d0)**2,
     &                               max(phi_y_plus(i,j),0.d0)**2 )
            endif

c           } end Godunov selection of grad_phi


c           compute contribution to lse_rhs(i,j) 
              lse_rhs(i,j) = lse_rhs(i,j) 
     &                     - vel_n*sqrt(norm_grad_phi_sq)

          enddo 
        enddo 
c       } end loop over grid

      endif
     

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm2dAddConstCurvTermToLSERHS() adds the contribution of a curvature 
c  term to the right-hand side of the level set equation when it is 
c  written in the form:
c
c    phi_t = b*kappa*|grad(phi)| + ...
c  
c  kappa (mean curvature, div ( grad*(phi) / |grad(phi)| ) will be computed 
c  from second order derivatives.
c
c  Arguments:
c    lse_rhs (in/out):  right-hand of level set equation
c    phi_*      (in):   derivatives (the 1st and 2nd order)   
c    b     (in):        scalar curvature term component 
c    *_gb (in):         index range for ghostbox
c    *_fb (in):         index range for fillbox
c
c***********************************************************************
      subroutine lsm2dAddConstCurvTermToLSERHS(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  phi_xx, phi_xy, phi_yy,
     &  ilo_grad2_phi_gb, ihi_grad2_phi_gb,
     &  jlo_grad2_phi_gb, jhi_grad2_phi_gb,
     &  b,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_grad2_phi_gb, ihi_grad2_phi_gb
      integer jlo_grad2_phi_gb, jhi_grad2_phi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_xx(ilo_grad2_phi_gb:ihi_grad2_phi_gb,
     &            jlo_grad2_phi_gb:jhi_grad2_phi_gb)
      real phi_yy(ilo_grad2_phi_gb:ihi_grad2_phi_gb,
     &            jlo_grad2_phi_gb:jhi_grad2_phi_gb)
      real phi_xy(ilo_grad2_phi_gb:ihi_grad2_phi_gb,
     &            jlo_grad2_phi_gb:jhi_grad2_phi_gb)
      real b
      
      integer i,j
      real grad_mag2, curv
      real zero_tol
      parameter (zero_tol=@lsmlib_zero_tol@)

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb
	
c         compute squared magnitude of gradient
	  grad_mag2 = phi_x(i,j) * phi_x(i,j) + phi_y(i,j) * phi_y(i,j) 

	  if(grad_mag2 .lt. zero_tol) then
	      curv = 0.d0
	  else
	    curv = phi_xx(i,j)*phi_y(i,j)*phi_y(i,j)  
     &	       +   phi_yy(i,j)*phi_x(i,j)*phi_x(i,j)  
     &	       - 2*phi_xy(i,j)*phi_x(i,j)*phi_y(i,j)  
	    curv = curv / grad_mag2 
	  endif

	  lse_rhs(i,j) = lse_rhs(i,j) + b*curv
      
        enddo 
      enddo 
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
