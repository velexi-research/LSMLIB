c***********************************************************************
c
c  File:        patchmodule_fort.f
c  Copyright:   (c) 2005-2006 Kevin T. Chu
c  Revision:    $Revision: 1.1 $
c  Modified:    $Date: 2006/04/10 15:25:52 $
c  Description: F77 patch module routines for LSMLIB example program
c
c***********************************************************************

c***********************************************************************
c  initializePeriodicArrayOfLines() initializes phi and psi for a
c  straight line with direction (sqrt(2/3),sqrt(1/3)) using level
c  set functions that are NOT orthogonal to begin with.
c
c  Arguments:
c    phi (out):               phi data
c    psi (out):               psi data
c    *_gb (in):               ghostbox for phi and psi data
c    *_fb (in):               fillbox for phi and psi data
c    x_lower (in):            physical coordinates of lower corner of box
c    dx (in):                 grid spacing
c
c***********************************************************************
      subroutine initializePeriodicArrayOfLines(
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb, 
     &  klo_phi_gb, khi_phi_gb,
     &  psi,
     &  ilo_psi_gb, ihi_psi_gb, 
     &  jlo_psi_gb, jhi_psi_gb, 
     &  klo_psi_gb, khi_psi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  x_lower, 
     &  dx)
c***********************************************************************
      implicit none

c     _gb refers to ghost box
c     _fb refers to fill box
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_psi_gb, ihi_psi_gb
      integer jlo_psi_gb, jhi_psi_gb
      integer klo_psi_gb, khi_psi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real psi(ilo_psi_gb:ihi_psi_gb,
     &         jlo_psi_gb:jhi_psi_gb,
     &         klo_psi_gb:khi_psi_gb)
      real x_lower(1:3)
      real dx(1:3)
      integer i,j,k
      real x,y,z
      real n_x, n_y, n_z
      real r, theta
      real temp1, temp2
      real tol
      parameter (tol = 1.0d-8)

c     normal vector for dislocation line
      n_x = -1/sqrt(3.0)
      n_y = sqrt(2.0/3.0)
      n_z = 0.0

c     loop over grid {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            x = x_lower(1) + dx(1)*(i-ilo_fb+0.5)
            y = x_lower(2) + dx(2)*(j-jlo_fb+0.5)
            z = x_lower(3) + dx(3)*(k-klo_fb+0.5)

            temp1 = n_x*(x-0.75) + n_y*y
            temp2 = n_x*(x+1.25) + n_y*y
            theta = atan2(-z,temp1) + atan2(-z,temp2)
            r = min(sqrt(z**2+temp1**2),sqrt(z**2+temp2**2))
            psi(i,j,k) = r * cos(theta)

            theta = atan2(z,n_x*x + n_y*y)
            phi(i,j,k) = sin(2*theta)

          enddo
        enddo
      enddo
c     } end loop over grid 

      return
      end
c***********************************************************************

