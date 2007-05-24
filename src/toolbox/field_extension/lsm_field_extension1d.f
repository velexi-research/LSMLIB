c***********************************************************************
c
c  File:        lsm_field_extension1d.f
c  Copyright:   (c) 2005-2006 Kevin T. Chu
c  Revision:    $Revision: 1.3 $
c  Modified:    $Date: 2006/05/18 01:50:44 $
c  Description: 2D F77 routines for extending fields off of the 
c               zero level set 
c
c***********************************************************************

c***********************************************************************
c
c  lsm1dComputeFieldExtensionEqnRHS() computes right-hand side of the 
c  field extension equation when it is written in the form:
c
c  S_t = -sgn(phi) N dot grad(S)
c
c  Arguments:
c    rhs (out):             right-hand side of field extension equation
c    S (in):                field to be extended off of the zero level set
c    phi (in):              level set function used to compute normal vector
c    S_*_upwind (in):       upwind spatial derivatives for grad(S)
c    signed_normal_* (in):  signed normal 
c    *_gb (in):             index range for ghostbox
c    *_fb (in):             index range for fillbox
c
c***********************************************************************
      subroutine lsm1dComputeFieldExtensionEqnRHS(
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  S,
     &  ilo_S_gb, ihi_S_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  S_x_upwind, 
     &  ilo_grad_S_upwind_gb, ihi_grad_S_upwind_gb,
     &  signed_normal_x,
     &  ilo_signed_normal_gb, ihi_signed_normal_gb,
     &  ilo_fb, ihi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fill-box

      integer ilo_rhs_gb, ihi_rhs_gb
      integer ilo_S_gb, ihi_S_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_grad_S_upwind_gb, ihi_grad_S_upwind_gb
      integer ilo_signed_normal_gb, ihi_signed_normal_gb
      integer ilo_fb, ihi_fb
      double precision rhs(ilo_rhs_gb:ihi_rhs_gb)
      double precision S(ilo_S_gb:ihi_S_gb)
      double precision phi(ilo_phi_gb:ihi_phi_gb)
      double precision S_x_upwind(
     &                   ilo_grad_S_upwind_gb:ihi_grad_S_upwind_gb)
      double precision signed_normal_x(
     &                   ilo_signed_normal_gb:ihi_signed_normal_gb)
      double precision dx
      integer i
      double precision zero
      parameter (zero=0.0d0)
      double precision zero_level_set_cutoff

c     set zero_level_set_cutoff to 3*dx
      zero_level_set_cutoff = 3.0d0*dx

c     compute RHS
c     { begin loop over grid
      do i=ilo_fb,ihi_fb

        if ( abs(phi(i)) .gt. zero_level_set_cutoff ) then

          rhs(i) = -signed_normal_x(i)*S_x_upwind(i)

        else

          if ( (phi(i+1)*phi(i) .le. zero) .and. 
     &         (phi(i-1)*phi(i) .le. zero) ) then

            rhs(i) = -signed_normal_x(i)*S_x_upwind(i)

          else

            rhs(i) = zero

          endif
        endif

      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
