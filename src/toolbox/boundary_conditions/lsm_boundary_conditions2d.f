c***********************************************************************
c
c  File:        lsm_boundary_conditions2d.f
c  Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
c  Revision:    $Revision: 1.6 $
c  Modified:    $Date$
c  Description: F77 routines for applying boundary conditions in 2D
c
c***********************************************************************

c***********************************************************************
c
c The boundary location index is used to identify the location of the 
c boundary relative to the computational domain.  In 2D, the boundary 
c location index conventions are:
c
c   x_lo: 0
c   x_hi: 1
c   y_lo: 2
c   y_hi: 3
c
c***********************************************************************

c***********************************************************************
c
c lsm2dLinearExtrapolation() extrapolates 2D data in from the index 
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
c  - if bdry_location_idx is out of the range for 2D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm2dLinearExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
c     local variables       
      integer i,j
      real dist, slope

      if (bdry_location_idx .eq. 0) then
c     { extrapolate data in x-direction at lower end

c       { begin j loop
        do j = jlo_gb, jhi_gb
          slope = phi(ilo_fb,j) - phi(ilo_fb+1,j)
          do i = ilo_gb, ilo_fb-1
            dist = ilo_fb - i
            phi(i,j) = phi(ilo_fb,j) + slope*dist
          enddo
        enddo  
c       } end j loop

c     } end extrapolate data in x-direction at lower end
         
      elseif (bdry_location_idx .eq. 1) then
c     {  extrapolate data in x-direction at upper end

c       { begin j loop
        do j = jlo_gb, jhi_gb
          slope = phi(ihi_fb,j) - phi(ihi_fb-1,j)
          do i = ihi_fb+1, ihi_gb
            dist = i - ihi_fb
            phi(i,j) = phi(ihi_fb,j) + slope*dist
          enddo 
        enddo  
c       } end j loop

c     } end extrapolate data in x-direction at upper end

      elseif (bdry_location_idx .eq. 2) then
c     {  extrapolate data in y-direction at lower end

c       { begin i loop
        do i = ilo_gb, ihi_gb
          slope = phi(i,jlo_fb) - phi(i,jlo_fb+1)
          do j = jlo_gb, jlo_fb-1
            dist = jlo_fb - j 
            phi(i,j) = phi(i,jlo_fb) + slope*dist
          enddo
        enddo
c       } end i loop

c     } end extrapolate data in y-direction at lower end

      elseif (bdry_location_idx .eq. 3) then
c     {  extrapolate data in y-direction at upper end

c       { begin i loop
        do i = ilo_gb, ihi_gb
          slope = phi(i,jhi_fb) - phi(i,jhi_fb-1)
          do j = jhi_fb+1, jhi_gb
            dist = j - jhi_fb
            phi(i,j) = phi(i,jhi_fb) + slope*dist
          enddo 
        enddo 
c       } end i loop

c     } end extrapolate data in y-direction at upper end

      endif

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm2dSignedLinearExtrapolation() extrapolates 2D data from the index 
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
c  - if bdry_location_idx is out of the range for 2D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm2dSignedLinearExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
c     local variables       
      integer i,j, zero
      parameter (zero = 0)
      real s, abs_diff, dist, slope
      real one
      parameter (one = 1.0d0)

      if (bdry_location_idx .eq. 0) then
c     { extrapolate data in x-direction at lower end

c       { begin j loop
        do j = jlo_gb, jhi_gb
          s = sign(one,phi(ilo_fb,j))
          abs_diff = abs(phi(ilo_fb,j) - phi(ilo_fb+1,j))
          slope = s*abs_diff
          do i = ilo_gb, ilo_fb-1
            dist = ilo_fb - i
            phi(i,j) = phi(ilo_fb,j) + slope*dist
          enddo
        enddo  
c       } end j loop

c     } end extrapolate data in x-direction at lower end
       
      elseif (bdry_location_idx .eq. 1) then
c     { extrapolate data in x-direction at upper end

c       { begin j loop
        do j = jlo_gb, jhi_gb
          s = sign(one,phi(ihi_fb,j))
          abs_diff = abs(phi(ihi_fb,j) - phi(ihi_fb-1,j))
          slope = s*abs_diff
          do i = ihi_fb+1, ihi_gb
            dist = i - ihi_fb
            phi(i,j) = phi(ihi_fb,j) + slope*dist
          enddo 
        enddo  
c       } end j loop

c     } end extrapolate data in x-direction at upper end

      elseif (bdry_location_idx .eq. 2) then
c     { extrapolate data in y-direction at lower end

c       { begin i loop
        do i = ilo_gb, ihi_gb
          s = sign(one,phi(i,jlo_fb))
          abs_diff = abs(phi(i,jlo_fb) - phi(i,jlo_fb+1))
          slope = s*abs_diff
          do j = jlo_gb, jlo_fb-1
            dist = jlo_fb - j 
            phi(i,j) = phi(i,jlo_fb) + slope*dist
          enddo
        enddo 
c       } end i loop

c     } end extrapolate data in y-direction at lower end

      elseif (bdry_location_idx .eq. 3) then
c     { extrapolate data in y-direction at upper end

c       { begin i loop
        do i = ilo_gb, ihi_gb
          s = sign(one,phi(i,jhi_fb))
          abs_diff = abs(phi(i,jhi_fb) - phi(i,jhi_fb-1))
          slope = s*abs_diff
          do j = jhi_fb+1, jhi_gb
            dist = j - jhi_fb
            phi(i,j) = phi(i,jhi_fb) + slope*dist
          enddo 
        enddo 
c       } end i loop

c     } end extrapolate data in y-direction at lower end

      endif

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm2dCopyExtrapolation() trivially extrapolates 2D data from the 
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
c  - if bdry_location_idx is out of the range for 2D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm2dCopyExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
c     local variables       
      integer i,j

      if (bdry_location_idx .eq. 0) then
c     { copy data in x-direction at lower end

c       { begin j loop
        do j = jlo_gb, jhi_gb
          do i = ilo_gb, ilo_fb-1
            phi(i,j) = phi(ilo_fb,j)
          enddo
        enddo  
c       } end j loop

c     } end copy data in x-direction at lower end
       
      elseif (bdry_location_idx .eq. 1) then
c     { copy data in x-direction at upper end

c       { begin j loop
        do j = jlo_gb, jhi_gb
          do i = ihi_fb+1, ihi_gb
            phi(i,j) = phi(ihi_fb,j)
          enddo 
        enddo  
c       } end j loop

c     } end copy data in x-direction at upper end

      elseif (bdry_location_idx .eq. 2) then
c     { copy data in y-direction at lower end

c       { begin i loop
        do i = ilo_gb, ihi_gb
          do j = jlo_gb, jlo_fb-1
            phi(i,j) = phi(i,jlo_fb)
          enddo
        enddo 
c       } end i loop

c     } end copy data in y-direction at lower end

      elseif (bdry_location_idx .eq. 3) then
c     { copy data in y-direction at upper end

c       { begin i loop
        do i = ilo_gb, ihi_gb
          do j = jhi_fb+1, jhi_gb
            phi(i,j) = phi(i,jhi_fb)
          enddo 
        enddo 
c       } end i loop

c     } end copy data in y-direction at upper end

      endif

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm2dHomogeneousNeumannENO1() sets the values of phi in the 
c ghostcells to impose a homogeneous Neumann boundary condition 
c at the specified boundary location for an ENO1 discretization 
c of the derivative.  In this case, the boundary condition reduces 
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
c  - if bdry_location_idx is out of the range for 2D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm2dHomogeneousNeumannENO1(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
      call lsm2dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &        ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm2dHomogeneousNeumannENO2() sets the values of phi in the ghostcells
c to impose a homogeneous Neumann boundary condition at the specified
c boundary location for an ENO2 discretization of the derivative.
c For this case, the boundary condition reduces to copy extrapolation.
c                
c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c 
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 2D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm2dHomogeneousNeumannENO2(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
      call lsm2dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &        ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm2dHomogeneousNeumannENO3() sets the values of phi in the ghostcells
c to impose a homogeneous Neumann boundary condition at the specified
c boundary location for an ENO3 discretization of the derivative.
c For this case, the boundary condition reduces to copy extrapolation.
c                
c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c 
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 2D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm2dHomogeneousNeumannENO3(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
      call lsm2dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &        ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm2dHomogeneousNeumannWENO5() sets the values of phi in the ghostcells
c to impose a homogeneous Neumann boundary condition at the specified
c boundary location for an WENO5 discretization of the derivative.
c For this case, the boundary condition reduces to copy extrapolation.
c                
c Arguments:
c   phi (in/out):            phi
c   bdry_location_idx (in):  boundary location index
c   *_gb (in):               index range for ghostbox
c   *_fb (in):               index range for fillbox
c 
c NOTES:
c  - fillbox indices must be a subset of ghostbox indices.
c  - if bdry_location_idx is out of the range for 2D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm2dHomogeneousNeumannWENO5(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
      call lsm2dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &        ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************
