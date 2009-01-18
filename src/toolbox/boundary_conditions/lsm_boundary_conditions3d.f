c***********************************************************************
c
c  File:        lsm_boundary_conditions3d.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 routines for applying boundary conditions in 3D
c
c***********************************************************************

c***********************************************************************
c  
c The boundary location index is used to identify the location of the
c boundary relative to the computational domain.  In 3D, the boundary
c location index conventions are:
c     
c   x_lo: 0
c   x_hi: 1
c   y_lo: 2
c   y_hi: 3
c   z_lo: 4
c   z_hi: 5
c
c***********************************************************************

c***********************************************************************
c
c lsm3dLinearExtrapolation() extrapolates 3D data in from the index
c range of the fillbox into cells in ghostbox at the specified boundary
c location using linear extrapolation.
c
c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c 
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 3D, then no
c    ghostcell values are set
c 
c***********************************************************************
      subroutine lsm3dLinearExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
c     local variables       
      integer i,j,k
      real dist, slope

      if (bdry_location_idx .eq. 0) then
c     { extrapolate data in x-direction at lower end

c       { begin k,j loop
        do k = klo_gb, khi_gb
          do j = jlo_gb, jhi_gb
            slope = phi(ilo_fb,j,k) - phi(ilo_fb+1,j,k)
            do i = ilo_gb, ilo_fb-1
              dist = ilo_fb - i
              phi(i,j,k) = phi(ilo_fb,j,k) + slope*dist
            enddo
          enddo
        enddo
c       } end k,j loop

c     } end extrapolate data in x-direction at lower end

      elseif (bdry_location_idx .eq. 1) then
c     { extrapolate data in x-direction at upper end
  
c       { begin k,j loop
        do k = klo_gb, khi_gb
          do j = jlo_gb, jhi_gb
            slope = phi(ihi_fb,j,k) - phi(ihi_fb-1,j,k)
            do i = ihi_fb+1, ihi_gb
              dist = i - ihi_fb
              phi(i,j,k) = phi(ihi_fb,j,k) + slope*dist
            enddo 
          enddo
        enddo
c       } end k,j loop

c     } extrapolate data in x-direction at upper end

      elseif (bdry_location_idx .eq. 2) then
c     { extrapolate data in y-direction at lower end

c       { begin i,k loop
        do i = ilo_gb, ihi_gb
          do k = klo_gb, khi_gb
            slope = phi(i,jlo_fb,k) - phi(i,jlo_fb+1,k)
            do j = jlo_gb, jlo_fb-1
              dist = jlo_fb - j
              phi(i,j,k) = phi(i,jlo_fb,k) + slope*dist
            enddo
          enddo
        enddo
c       } end i,k loop

c     } end extrapolate data in y-direction at lower end

      elseif (bdry_location_idx .eq. 3) then
c     { extrapolate data in y-direction at upper end

c       { begin i,k loop
        do i = ilo_gb, ihi_gb
          do k = klo_gb, khi_gb
            slope = phi(i,jhi_fb,k) - phi(i,jhi_fb-1,k)
            do j = jhi_fb+1, jhi_gb
              dist = j - jhi_fb
              phi(i,j,k) = phi(i,jhi_fb,k) + slope*dist
            enddo 
          enddo
        enddo
c       } end i,k loop

c     } end extrapolate data in y-direction at upper end

      elseif (bdry_location_idx .eq. 4) then
c     { extrapolate data in z-direction at lower end

c       { begin i,j loop
        do i = ilo_gb, ihi_gb
          do j = jlo_gb, jhi_gb
            slope = phi(i,j,klo_fb) - phi(i,j,klo_fb+1)
            do k = klo_gb, klo_fb-1
              dist = klo_fb - k
              phi(i,j,k) = phi(i,j,klo_fb) + slope*dist
            enddo
          enddo
        enddo
c       } end i,j loop

c     } end extrapolate data in z-direction at lower end

      elseif (bdry_location_idx .eq. 5) then
c     { extrapolate data in z-direction at upper end

c       { begin i,j loop
        do i = ilo_gb, ihi_gb
          do j = jlo_gb, jhi_gb
            slope = phi(i,j,khi_fb) - phi(i,j,khi_fb-1)
            do k = khi_fb+1, khi_gb
              dist = k - khi_fb
              phi(i,j,k) = phi(i,j,khi_fb) + slope*dist
            enddo 
          enddo
        enddo
c       } end i,j loop

c     } end extrapolate data in z-direction at upper end

      endif

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm3dSignedLinearExtrapolation() extrapolates 3D data from the index
c range of the fillbox into cells in ghostbox at the specified boundary
c location using signed linear extrapolation.  Extrapolation is "away" 
c from the zero level set, i.e. the zero level set will not be
c artifically created at the boundary.
c
c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c 
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 3D, then no
c    ghostcell values are set
c 
c***********************************************************************
      subroutine lsm3dSignedLinearExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
c     local variables       
      integer i,j,k
      real s, abs_diff, dist, slope
      real one
      parameter (one = 1.0d0)

      if (bdry_location_idx .eq. 0) then
c     { extrapolate data in x-direction at lower end

c       { begin k,j loop
        do k = klo_gb, khi_gb
          do j = jlo_gb, jhi_gb
            s = sign(one,phi(ilo_fb,j,k))
            abs_diff = abs(phi(ilo_fb,j,k) - phi(ilo_fb+1,j,k))
            slope = s*abs_diff
            do i = ilo_gb, ilo_fb-1
              dist = ilo_fb - i
              phi(i,j,k) = phi(ilo_fb,j,k) + slope*dist
            enddo
          enddo
        enddo
c       } end k,j loop

c     } end extrapolate data in x-direction at lower end

      elseif (bdry_location_idx .eq. 1) then
c     { extrapolate data in x-direction at upper end

c       { begin k,j loop
        do k = klo_gb, khi_gb
          do j = jlo_gb, jhi_gb
            s = sign(one,phi(ihi_fb,j,k))
            abs_diff = abs(phi(ihi_fb,j,k) - phi(ihi_fb-1,j,k))
            slope = s*abs_diff
            do i = ihi_fb+1, ihi_gb
              dist = i - ihi_fb
              phi(i,j,k) = phi(ihi_fb,j,k) + slope*dist
            enddo 
          enddo
        enddo
c       } end k,j loop

c     } end extrapolate data in x-direction at upper end

      elseif (bdry_location_idx .eq. 2) then
c     { extrapolate data in y-direction at lower end

c       { begin i,k loop
        do i = ilo_gb, ihi_gb
          do k = klo_gb, khi_gb
            s = sign(one,phi(i,jlo_fb,k))
            abs_diff = abs(phi(i,jlo_fb,k) - phi(i,jlo_fb+1,k))
            slope = s*abs_diff
            do j = jlo_gb, jlo_fb-1
              dist = jlo_fb - j
              phi(i,j,k) = phi(i,jlo_fb,k) + slope*dist
            enddo
          enddo
        enddo
c       } end i,k loop

c     } end extrapolate data in y-direction at lower end

      elseif (bdry_location_idx .eq. 3) then
c     { extrapolate data in y-direction at upper end

c       { begin i,k loop
        do i = ilo_gb, ihi_gb
          do k = klo_gb, khi_gb
            s = sign(one,phi(i,jhi_fb,k))
            abs_diff = abs(phi(i,jhi_fb,k) - phi(i,jhi_fb-1,k))
            slope = s*abs_diff
            do j = jhi_fb+1, jhi_gb
              dist = j - jhi_fb
              phi(i,j,k) = phi(i,jhi_fb,k) + slope*dist
            enddo 
          enddo
        enddo
c       } end i,k loop

c     } end extrapolate data in y-direction at upper end

      elseif (bdry_location_idx .eq. 4) then
c     { extrapolate data in z-direction at lower end

c       { begin i,j loop
        do i = ilo_gb, ihi_gb
          do j = jlo_gb, jhi_gb
            s = sign(one,phi(i,j,klo_fb))
            abs_diff = abs(phi(i,j,klo_fb) - phi(i,j,klo_fb+1))
            slope = s*abs_diff
            do k = klo_gb, klo_fb-1
              dist = klo_fb - k
              phi(i,j,k) = phi(i,j,klo_fb) + slope*dist
            enddo
          enddo
        enddo
c       } end i,j loop

c     } end extrapolate data in z-direction at lower end

      elseif (bdry_location_idx .eq. 5) then
c     { extrapolate data in z-direction at upper end

c       { begin i,j loop
        do i = ilo_gb, ihi_gb
          do j = jlo_gb, jhi_gb
            s = sign(one,phi(i,j,khi_fb))
            abs_diff = abs(phi(i,j,khi_fb) - phi(i,j,khi_fb-1))
            slope = s*abs_diff
            do k = khi_fb+1, khi_gb
              dist = k - khi_fb
              phi(i,j,k) = phi(i,j,khi_fb) + slope*dist
            enddo 
          enddo
        enddo
c       } end i,j loop

c     } end extrapolate data in z-direction at upper end

      endif
      return
      end
c } end subroutine
c***********************************************************************

***********************************************************************
c
c lsm3dCopyExtrapolation() trivially extrapolates 3D data from the
c index range of the fillbox into cells in ghostbox at the specified
c boundary location by merely copying data from the closest grid cell
c in the fillbox. 
c                
c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 3D, then no
c    ghostcell values are set
c 
c***********************************************************************
      subroutine lsm3dCopyExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
c     local variables       
      integer i,j,k
 
      if (bdry_location_idx .eq. 0) then
c     { copy data in x-direction at lower end

c       { begin j,k loop
        do j = jlo_gb, jhi_gb
          do k = klo_gb, khi_gb
            do i = ilo_gb, ilo_fb-1
              phi(i,j,k) = phi(ilo_fb,j,k)
            enddo
          enddo
        enddo 
c       } end j,k loop

c     } end copy data in x-direction at lower end

      elseif (bdry_location_idx .eq. 1) then
c     { copy data in x-direction at upper end

c       { begin j,k loop
        do j = jlo_gb, jhi_gb
          do k = klo_gb, khi_gb
            do i = ihi_fb+1, ihi_gb
              phi(i,j,k) = phi(ihi_fb,j,k)
            enddo 
          enddo
        enddo 
c       } end j,k loop

c     } end copy data in x-direction at upper end

      elseif (bdry_location_idx .eq. 2) then
c     { copy data in y-direction at lower end

c       { begin i,k loop
        do i = ilo_gb, ihi_gb
          do k = klo_gb, khi_gb
            do j = jlo_gb, jlo_fb-1
              phi(i,j,k) = phi(i,jlo_fb,k)
            enddo
          enddo
        enddo 
c       } end i,k loop

c     } end copy data in y-direction at lower end

      elseif (bdry_location_idx .eq. 3) then
c     { copy data in y-direction at upper end

c       { begin i,k loop
        do i = ilo_gb, ihi_gb
          do k = klo_gb, khi_gb
            do j = jhi_fb+1, jhi_gb
              phi(i,j,k) = phi(i,jhi_fb,k)
            enddo 
          enddo
        enddo 
c       } end i,k loop

c     } end copy data in y-direction at upper end

      elseif (bdry_location_idx .eq. 4) then
c     { copy data in z-direction at lower end

c       { begin i,j loop
        do i = ilo_gb, ihi_gb
          do j = jlo_gb, jhi_gb
            do k = klo_gb, klo_fb-1
              phi(i,j,k) = phi(i,j,klo_fb)
            enddo
          enddo
        enddo 
c       } end i,j loop

c     } end copy data in z-direction at lower end

      elseif (bdry_location_idx .eq. 5) then
c     { copy data in z-direction at upper end

c       { begin i,j loop
        do i = ilo_gb, ihi_gb
          do j = jlo_gb, jhi_gb
            do k = khi_fb+1, khi_gb
              phi(i,j,k) = phi(i,j,khi_fb)
            enddo 
          enddo
        enddo 
c       } end i,j loop

c     } copy data in z-direction at upper end

      endif

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm3dHomogeneousNeumannENO1() sets the values of phi in the
c ghostcells to impose a homogeneous Neumann boundary condition
c at the specified boundary location for an ENO1 discretization
c of the derivative.  In this case, the boundary condition reduces
c to copy extrapolation.

c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c 
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 3D, then no
c    ghostcell values are set
c
c***********************************************************************
      subroutine lsm3dHomogeneousNeumannENO1(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
      call lsm3dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &        ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm3dHomogeneousNeumannENO2() sets the values of phi in the
c ghostcells to impose a homogeneous Neumann boundary condition
c at the specified boundary location for an ENO2 discretization
c of the derivative.  For this case, the boundary condition reduces 
c to copy extrapolation.
c
c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c 
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 3D, then no
c    ghostcell values are set
c
c***********************************************************************
      subroutine lsm3dHomogeneousNeumannENO2(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
      call lsm3dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &        ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm3dHomogeneousNeumannENO3() sets the values of phi in the
c ghostcells to impose a homogeneous Neumann boundary condition
c at the specified boundary location for an ENO3 discretization
c of the derivative.  For this case, the boundary condition reduces 
c to copy extrapolation.

c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c 
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 3D, then no
c    ghostcell values are set
c
c***********************************************************************
      subroutine lsm3dHomogeneousNeumannENO3(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
      call lsm3dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &        ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm3dHomogeneousNeumannWENO5() sets the values of phi in the
c ghostcells to impose a homogeneous Neumann boundary condition
c at the specified boundary location for an WENO5 discretization
c of the derivative.  For this case, the boundary condition reduces 
c to copy extrapolation.
c
c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c 
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 3D, then no
c    ghostcell values are set
c
c***********************************************************************
      subroutine lsm3dHomogeneousNeumannWENO5(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
      call lsm3dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb,
     &        ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************
