c***********************************************************************
c
c  File:        lsm_reinitialization2d.f
c  Copyright:   (c) 2005-2006 Kevin T. Chu
c  Revision:    $Revision: 1.8 $
c  Modified:    $Date: 2006/01/24 21:46:42 $
c  Description: F77 routines for reinitialization of 2d level set functions
c
c***********************************************************************

c***********************************************************************
c The algorithms and notation in these subroutines closely follows
c the discussion in Osher & Fedkiw (2003).
c***********************************************************************

c***********************************************************************
c
c  lsm2dComputeReinitializationEqnRHS() computes the right-hand side of 
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
c       behavior of lsm2dComputeReinitializationEqnRHS() is to use 
c       phi (i.e. equivalent to setting use_phi0_for_sgn to 0)
c
c***********************************************************************
      subroutine lsm2dComputeReinitializationEqnRHS(
     &  reinit_rhs,
     &  ilo_rhs_gb, ihi_rhs_gb, jlo_rhs_gb, jhi_rhs_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  phi0,
     &  ilo_phi0_gb, ihi_phi0_gb, jlo_phi0_gb, jhi_phi0_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy,
     &  use_phi0_for_sgn)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_rhs_gb, ihi_rhs_gb, jlo_rhs_gb, jhi_rhs_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_phi0_gb, ihi_phi0_gb, jlo_phi0_gb, jhi_phi0_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real reinit_rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                            jlo_rhs_gb:jhi_rhs_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &                     jlo_phi_gb:jhi_phi_gb)
      real phi0(ilo_phi0_gb:ihi_phi0_gb,
     &                      jlo_phi0_gb:jhi_phi0_gb)
      real phi_x_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real dx, dy
      integer use_phi0_for_sgn
      real phi_cur
      integer DIM
      parameter (DIM=2)
      real grad_phi_plus_cur(1:DIM)
      real grad_phi_minus_cur(1:DIM)
      real grad_phi_star(1:DIM)
      integer i,j
      integer dir
      real sgn_phi
      real norm_grad_phi_sq
      real dx_sq
      real tol
      parameter (tol=1.d-13)
      real one
      parameter (one=1.d0)

c     set value of dx_sq to be square of max{dx,dy}
      dx_sq = max(dx,dy)
      dx_sq = dx_sq*dx_sq

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
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

c           cache phi and spatial derivative approximations
            phi_cur = phi(i,j)
            grad_phi_plus_cur(1) = phi_x_plus(i,j)
            grad_phi_plus_cur(2) = phi_y_plus(i,j)
            grad_phi_minus_cur(1) = phi_x_minus(i,j)
            grad_phi_minus_cur(2) = phi_y_minus(i,j)

c           { begin Godunov selection of grad_phi
            do dir=1,DIM

              if (phi_cur .gt. 0.d0) then
                grad_phi_plus_cur(dir) = max(-grad_phi_plus_cur(dir),
     &                                       0.d0)
                grad_phi_minus_cur(dir) = max(grad_phi_minus_cur(dir),
     &                                        0.d0)
              else
                grad_phi_plus_cur(dir) = max(grad_phi_plus_cur(dir),
     &                                       0.d0)
                grad_phi_minus_cur(dir) = max(-grad_phi_minus_cur(dir),
     &                                        0.d0)
              endif

              grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                 grad_phi_minus_cur(dir)) 

            enddo
c           } end Godunov selection of grad_phi

c           compute reinit_rhs(i,j) using smoothed sgn(phi)
            if (abs(phi_cur) .gt. tol) then 
              norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                         + grad_phi_star(2)*grad_phi_star(2)
              sgn_phi = phi_cur
     &                / sqrt(phi_cur*phi_cur + norm_grad_phi_sq*dx_sq)
              reinit_rhs(i,j) = sgn_phi*(one - sqrt(norm_grad_phi_sq))
            else
              reinit_rhs(i,j) = 0.d0
            endif

          enddo
        enddo
c       } end loop over grid

      else

c       ------------------------------------------------
c       use phi0 in computation of smoothed sgn function 
c       ------------------------------------------------
c       { begin loop over grid
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

c           cache phi and spatial derivative approximations
            phi_cur = phi0(i,j)
            grad_phi_plus_cur(1) = phi_x_plus(i,j)
            grad_phi_plus_cur(2) = phi_y_plus(i,j)
            grad_phi_minus_cur(1) = phi_x_minus(i,j)
            grad_phi_minus_cur(2) = phi_y_minus(i,j)

c           { begin Godunov selection of grad_phi
            do dir=1,DIM

              if (phi_cur .gt. 0.d0) then
                grad_phi_plus_cur(dir) = max(-grad_phi_plus_cur(dir),
     &                                       0.d0)
                grad_phi_minus_cur(dir) = max(grad_phi_minus_cur(dir),
     &                                        0.d0)
              else
                grad_phi_plus_cur(dir) = max(grad_phi_plus_cur(dir),
     &                                       0.d0)
                grad_phi_minus_cur(dir) = max(-grad_phi_minus_cur(dir),
     &                                        0.d0)
              endif

              grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                 grad_phi_minus_cur(dir)) 

            enddo
c           } end Godunov selection of phi_y in y-direction

c           compute reinit_rhs(i,j) using smoothed sgn(phi)
            if (abs(phi_cur) .gt. tol) then
              norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                         + grad_phi_star(2)*grad_phi_star(2)
              sgn_phi = phi_cur / sqrt(phi_cur*phi_cur + dx_sq)
              reinit_rhs(i,j) = sgn_phi*(one - sqrt(norm_grad_phi_sq))
            else
              reinit_rhs(i,j) = 0.d0
            endif

          enddo
        enddo
c       } end loop over grid

      endif
c     } end condition on use_phi0_for_sgn

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dComputeOrthogonalizationEqnRHS() computes the right-hand side 
c  of the orthogonalization equation:
c
c    phi_t + grad(phi) dot { sgn(psi)/|grad(psi)| grad(psi) } = 0
c
c  Upwinding is used to select whether the forward (plus) or backward
c  (minus) spatial derivative should be used for grad(phi).  grad(psi)
c  is computed by averaging the forward and backward spatial derivatives
c  for grad(psi).  Forward and backward spatial derivatives used in the 
c  calculation must be supplied by the user.
c
c  Arguments:
c    othro_rhs (out):        right-hand side of orthogonalization
c                            equation
c    psi (in):               data array for psi
c    phi_*_plus (in):        forward spatial derivatives for grad(phi)
c    phi_*_minus (in):       backward spatial derivatives for grad(phi)
c    psi_*_plus (in):        forward spatial derivatives for grad(psi)
c    psi_*_minus (in):       backward spatial derivatives for grad(psi)
c    *_gb (in):              index range for ghostbox
c    *_fb (in):              index range for fillbox
c
c***********************************************************************
      subroutine lsm2dComputeOrthogonalizationEqnRHS(
     &  ortho_rhs,
     &  ilo_rhs_gb, ihi_rhs_gb, 
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  psi,
     &  ilo_psi_gb, ihi_psi_gb, 
     &  jlo_psi_gb, jhi_psi_gb,
     &  psi_x_plus, psi_y_plus,
     &  ilo_grad_psi_plus_gb, ihi_grad_psi_plus_gb, 
     &  jlo_grad_psi_plus_gb, jhi_grad_psi_plus_gb,
     &  psi_x_minus, psi_y_minus,
     &  ilo_grad_psi_minus_gb, ihi_grad_psi_minus_gb, 
     &  jlo_grad_psi_minus_gb, jhi_grad_psi_minus_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, 
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_psi_gb, ihi_psi_gb
      integer jlo_psi_gb, jhi_psi_gb
      integer ilo_grad_psi_plus_gb, ihi_grad_psi_plus_gb
      integer jlo_grad_psi_plus_gb, jhi_grad_psi_plus_gb
      integer ilo_grad_psi_minus_gb, ihi_grad_psi_minus_gb
      integer jlo_grad_psi_minus_gb, jhi_grad_psi_minus_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real ortho_rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                           jlo_rhs_gb:jhi_rhs_gb)
      real psi(ilo_psi_gb:ihi_psi_gb,
     &                     jlo_psi_gb:jhi_psi_gb)
      real phi_x_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real psi_x_plus(
     &                    ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                    jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb)
      real psi_y_plus(
     &                    ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                    jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb)
      real psi_x_minus(
     &                    ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                    jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb)
      real psi_y_minus(
     &                    ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                    jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb)
      real dx, dy
      real dx_sq
      integer DIM
      parameter (DIM=2)
      real grad_psi_star(1:DIM)
      real norm_grad_psi
      real sgn_psi
      integer i,j
      real psi_tol, grad_psi_tol
      parameter (psi_tol=1.d-13,grad_psi_tol=1.d-8)
      real one, half
      parameter (one=1.d0,half=0.5d0)

c     set value of dx_sq to be square of max{dx,dy}
      dx_sq = max(dx,dy)
      dx_sq = dx_sq*dx_sq

c----------------------------------------------------
c      compute RHS of orthogonalization equation
c      using upwinding
c----------------------------------------------------

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         average forward and backward derivatives of 
c         grad(psi) to get grad_psi_star
          grad_psi_star(1) = half 
     &                     * (psi_x_plus(i,j)+ psi_x_minus(i,j))
          grad_psi_star(2) = half 
     &                     * (psi_y_plus(i,j)+ psi_y_minus(i,j))
 
c         compute norm of grad(psi)
          norm_grad_psi = grad_psi_star(1)*grad_psi_star(1)
     &                  + grad_psi_star(2)*grad_psi_star(2)
          norm_grad_psi = sqrt(norm_grad_psi)

c         compute ortho_rhs(i,j) using upwinding on sgn_psi*grad(psi) 
          if ( (abs(psi(i,j)) .gt. psi_tol) .and.
     &         (norm_grad_psi .gt. grad_psi_tol) ) then

c           CASE: nontrivial psi and grad(psi) 

            sgn_psi = psi(i,j)/sqrt(psi(i,j)*psi(i,j) + dx_sq)

            ortho_rhs(i,j) = -1.d0/norm_grad_psi 
     &        * ( max(sgn_psi*grad_psi_star(1),0.d0)
     &                       *phi_x_minus(i,j)
     &          + min(sgn_psi*grad_psi_star(1),0.d0)
     &                        *phi_x_plus(i,j)
     &          + max(sgn_psi*grad_psi_star(2),0.d0)
     &                        *phi_y_minus(i,j)
     &          + min(sgn_psi*grad_psi_star(2),0.d0)
     &                        *phi_y_plus(i,j) )
          else

c           CASE: grad(psi) = 0 CASE: psi = 0 

            ortho_rhs(i,j) = 0.d0

          endif

        enddo
      enddo
c     } end loop over grid


      return
      end
c } end subroutine
c***********************************************************************
