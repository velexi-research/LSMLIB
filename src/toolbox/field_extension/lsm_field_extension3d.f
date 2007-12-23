c***********************************************************************
c
c  File:        lsm_field_extension3d.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision: 1.3 $
c  Modified:    $Date: 2006/05/18 01:50:44 $
c  Description: 3D F77 routines for extending fields off of the 
c               zero level set 
c
c***********************************************************************

c***********************************************************************
c
c  lsm3dComputeFieldExtensionEqnRHS() computes right-hand side of the 
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
      subroutine lsm3dComputeFieldExtensionEqnRHS(
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  S,
     &  ilo_S_gb, ihi_S_gb,
     &  jlo_S_gb, jhi_S_gb,
     &  klo_S_gb, khi_S_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  S_x_upwind, S_y_upwind, S_z_upwind,
     &  ilo_grad_S_upwind_gb, ihi_grad_S_upwind_gb,
     &  jlo_grad_S_upwind_gb, jhi_grad_S_upwind_gb,
     &  klo_grad_S_upwind_gb, khi_grad_S_upwind_gb,
     &  signed_normal_x, signed_normal_y, signed_normal_z,
     &  ilo_signed_normal_gb, ihi_signed_normal_gb,
     &  jlo_signed_normal_gb, jhi_signed_normal_gb,
     &  klo_signed_normal_gb, khi_signed_normal_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  klo_fb, khi_fb,
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fill-box

      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      integer ilo_S_gb, ihi_S_gb
      integer jlo_S_gb, jhi_S_gb
      integer klo_S_gb, khi_S_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_grad_S_upwind_gb, ihi_grad_S_upwind_gb
      integer jlo_grad_S_upwind_gb, jhi_grad_S_upwind_gb
      integer klo_grad_S_upwind_gb, khi_grad_S_upwind_gb
      integer ilo_signed_normal_gb, ihi_signed_normal_gb
      integer jlo_signed_normal_gb, jhi_signed_normal_gb
      integer klo_signed_normal_gb, khi_signed_normal_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      integer klo_fb, khi_fb
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &         jlo_rhs_gb:jhi_rhs_gb,
     &         klo_rhs_gb:khi_rhs_gb)
      real S(ilo_S_gb:ihi_S_gb,
     &       jlo_S_gb:jhi_S_gb,
     &       klo_S_gb:khi_S_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real S_x_upwind(ilo_grad_S_upwind_gb:ihi_grad_S_upwind_gb,
     &                jlo_grad_S_upwind_gb:jhi_grad_S_upwind_gb,
     &                klo_grad_S_upwind_gb:khi_grad_S_upwind_gb)
      real S_y_upwind(ilo_grad_S_upwind_gb:ihi_grad_S_upwind_gb,
     &                jlo_grad_S_upwind_gb:jhi_grad_S_upwind_gb,
     &                klo_grad_S_upwind_gb:khi_grad_S_upwind_gb)
      real S_z_upwind(ilo_grad_S_upwind_gb:ihi_grad_S_upwind_gb,
     &                jlo_grad_S_upwind_gb:jhi_grad_S_upwind_gb,
     &                klo_grad_S_upwind_gb:khi_grad_S_upwind_gb)
      real signed_normal_x(ilo_signed_normal_gb:ihi_signed_normal_gb,
     &                     jlo_signed_normal_gb:jhi_signed_normal_gb,
     &                     klo_signed_normal_gb:khi_signed_normal_gb)
      real signed_normal_y(ilo_signed_normal_gb:ihi_signed_normal_gb,
     &                     jlo_signed_normal_gb:jhi_signed_normal_gb,
     &                     klo_signed_normal_gb:khi_signed_normal_gb)
      real signed_normal_z(ilo_signed_normal_gb:ihi_signed_normal_gb,
     &                     jlo_signed_normal_gb:jhi_signed_normal_gb,
     &                     klo_signed_normal_gb:khi_signed_normal_gb)
      real dx, dy, dz
      integer i,j,k
      real zero
      parameter (zero=0.0d0)
      real zero_level_set_cutoff

c     set zero_level_set_cutoff to 3*max(dx,dy,dz)
      zero_level_set_cutoff = 3.0d0*max(dx,dy,dz)

c     compute RHS
c     { begin loop over grid
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            if ( abs(phi(i,j,k)) .gt. zero_level_set_cutoff ) then

              rhs(i,j,k) = -( signed_normal_x(i,j,k)*S_x_upwind(i,j,k)
     &                      + signed_normal_y(i,j,k)*S_y_upwind(i,j,k)
     &                      + signed_normal_z(i,j,k)*S_z_upwind(i,j,k) )

            else

              if ( (phi(i,j,k)*phi(i-1,j,k) .gt. zero) .and. 
     &             (phi(i,j,k)*phi(i+1,j,k) .gt. zero) .and.
     &             (phi(i,j,k)*phi(i,j-1,k) .gt. zero) .and.
     &             (phi(i,j,k)*phi(i,j+1,k) .gt. zero) .and.
     &             (phi(i,j,k)*phi(i,j,k-1) .gt. zero) .and.
     &             (phi(i,j,k)*phi(i,j,k+1) .gt. zero) ) then

                rhs(i,j,k) = 
     &             -( signed_normal_x(i,j,k)*S_x_upwind(i,j,k)
     &              + signed_normal_y(i,j,k)*S_y_upwind(i,j,k)
     &              + signed_normal_z(i,j,k)*S_z_upwind(i,j,k) )

              else

                rhs(i,j,k) = zero

              endif
            endif

          enddo
        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
