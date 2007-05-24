c***********************************************************************
c
c  File:        lsm_reinitialization1d.f
c  Copyright:   (c) 2005-2006 Kevin T. Chu
c  Revision:    $Revision: 1.8 $
c  Modified:    $Date: 2006/01/24 21:46:42 $
c  Description: F77 routines for reinitialization of 1d level set functions
c
c***********************************************************************

c***********************************************************************
c The algorithms and notation in these subroutines closely follows
c the discussion in Osher & Fedkiw (2003).
c***********************************************************************

c***********************************************************************
c
c  lsm1dComputeReinitializationEqnRHS() computes the right-hand side of 
c  the reinitialization equation using a Godunov scheme to select the 
c  numerical discretization of the sgn(phi) |grad(phi)| term.  
c  Forward (plus) and backward (minus) spatial derivatives used in 
c  the Godunov calculation must be supplied by the user.
c
c  Arguments:
c    reinit_rhs (out):       right-hand side of reinitialization 
c                            equation
c    phi (in):               level set function at current iteration
c                            of reinitialization process
c    phi0 (in):              level set function at initial iteration
c                            iteration of reinitialization process
c    phi_*_plus (in):        forward spatial derivatives for grad(phi)
c    phi_*_minus (in):       backward spatial derivatives for grad(phi)
c    use_phi0_for_sgn (in):  flag to specify whether phi0 should be
c                            used in the computation of sgn(phi).
c                              0 = use phi (do NOT use phi0)
c                              1 = use phi0
c    *_gb (in):              index range for ghostbox
c    *_fb (in):              index range for fillbox
c
c  NOTES:
c   (1) if use_phi0_for_sgn is not equal to 0 or 1, the default
c       behavior of lsm1dComputeReinitializationEqnRHS() is to use 
c       phi (i.e. equivalent to setting use_phi0_for_sgn to 0)
c
c***********************************************************************
      subroutine lsm1dComputeReinitializationEqnRHS(
     &  reinit_rhs,
     &  ilo_rhs_gb, ihi_rhs_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  phi0,
     &  ilo_phi0_gb, ihi_phi0_gb, 
     &  phi_x_plus, 
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  phi_x_minus, 
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  ilo_fb, ihi_fb, 
     &  dx,
     &  use_phi0_for_sgn)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_rhs_gb, ihi_rhs_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_phi0_gb, ihi_phi0_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer ilo_fb, ihi_fb
      double precision reinit_rhs(ilo_rhs_gb:ihi_rhs_gb)
      double precision phi(ilo_phi_gb:ihi_phi_gb)
      double precision phi0(ilo_phi0_gb:ihi_phi0_gb)
      double precision phi_x_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb)
      double precision phi_x_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb)
      double precision dx
      integer use_phi0_for_sgn
      double precision phi_cur
      double precision grad_phi_plus_cur
      double precision grad_phi_minus_cur
      double precision grad_phi_star
      integer i
      double precision sgn_phi
      double precision norm_grad_phi
      double precision dx_sq
      double precision tol
      parameter (tol=1.d-13)
      double precision one
      parameter (one=1.d0)

c     set value of dx_sq to be square of dx
      dx_sq = dx*dx

c----------------------------------------------------
c      compute RHS of reinitialization equation
c      using Godunov's method
c----------------------------------------------------
c     { begin condition on use_phi0_for_sgn
      if (use_phi0_for_sgn .ne. 1) then

c       -----------------------------------------------
c       use phi in computation of smoothed sgn function
c       -----------------------------------------------
c       { begin loop over grid
        do i=ilo_fb,ihi_fb

c         cache phi and spatial derivative approximations
          phi_cur = phi(i)
          grad_phi_plus_cur = phi_x_plus(i)
          grad_phi_minus_cur = phi_x_minus(i)

c         { begin Godunov selection of grad_phi

          if (phi_cur .gt. 0.d0) then
            grad_phi_plus_cur = max(-grad_phi_plus_cur,0.d0)
            grad_phi_minus_cur = max(grad_phi_minus_cur,0.d0)
          else
            grad_phi_plus_cur = max(grad_phi_plus_cur,0.d0)
            grad_phi_minus_cur = max(-grad_phi_minus_cur,0.d0)
          endif

          grad_phi_star = max(grad_phi_plus_cur,grad_phi_minus_cur)

c         } end Godunov selection of grad_phi


c         compute reinit_rhs(i) using smoothed sgn(phi)
          if (abs(phi_cur) .ge. tol) then
            norm_grad_phi = abs(grad_phi_star)
            sgn_phi = phi_cur/sqrt( phi_cur*phi_cur 
     &                            + norm_grad_phi*norm_grad_phi*dx_sq)
            reinit_rhs(i) = sgn_phi*(one - norm_grad_phi)
          else
            reinit_rhs(i) = 0.d0
          endif

        enddo
c       } end loop over grid

      else

c       ------------------------------------------------
c       use phi0 in computation of smoothed sgn function 
c       ------------------------------------------------
c       { begin loop over grid
        do i=ilo_fb,ihi_fb

c         cache phi and spatial derivative approximations
          phi_cur = phi0(i)
          grad_phi_plus_cur = phi_x_plus(i)
          grad_phi_minus_cur = phi_x_minus(i)

c         { begin Godunov selection of grad_phi

          if (phi_cur .gt. 0.d0) then
            grad_phi_plus_cur = max(-grad_phi_plus_cur,0.d0)
            grad_phi_minus_cur = max(grad_phi_minus_cur,0.d0)
          else
            grad_phi_plus_cur = max(grad_phi_plus_cur,0.d0)
            grad_phi_minus_cur = max(-grad_phi_minus_cur,0.d0)
          endif

          grad_phi_star = max(grad_phi_plus_cur,grad_phi_minus_cur)

c         } end Godunov selection of grad_phi

c         compute reinit_rhs(i) using smoothed sgn(phi)
          if (abs(phi_cur) .ge. tol) then
            norm_grad_phi = abs(grad_phi_star)
            sgn_phi = phi_cur / sqrt(phi_cur*phi_cur + dx_sq)
            reinit_rhs(i) = sgn_phi*(one - norm_grad_phi)
          else
            reinit_rhs(i) = 0.d0
          endif

        enddo
c       } end loop over grid

      endif
c     } end condition on use_phi0_for_sgn

      return
      end
c } end subroutine
c***********************************************************************

