c***********************************************************************
c
c  File:        lsm_boundary_conditions1d.f
c  Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
c  Revision:    $Revision: 1.4 $
c  Modified:    $Date$
c  Description: F77 routines for applying boundary conditions in 1D
c
c***********************************************************************

c***********************************************************************
c
c The boundary location index is used to identify the location of the 
c boundary relative to the computational domain.  In 1D, the boundary 
c location index conventions are:
c
c   x_lo: 0
c   x_hi: 1
c
c***********************************************************************

c***********************************************************************
c
c lsm1dLinearExtrapolation() extrapolates 1D data in from the index 
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
c  - if bdry_location_idx is out of the range for 1D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm1dLinearExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb,
     &  ilo_fb, ihi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer ilo_fb, ihi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb)
      
c     local variables       
      integer i
      real dist, slope

c     lower end extrapolation   
      if (bdry_location_idx .eq. 0) then
        slope = phi(ilo_fb) - phi(ilo_fb+1)
        do i = ilo_gb, ilo_fb-1
          dist = ilo_fb - i
          phi(i) = phi(ilo_fb) + slope*dist
        enddo
      endif
     
c     upper end extrapolation
      if (bdry_location_idx .eq. 1) then
        slope = phi(ihi_fb) - phi(ihi_fb-1)
        do i = ihi_fb+1, ihi_gb
          dist = i - ihi_fb
          phi(i) = phi(ihi_fb) + slope*dist
        enddo 
      endif

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm1dSignedLinearExtrapolation() extrapolates 1D data from the index 
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
c  - if bdry_location_idx is out of the range for 1D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm1dSignedLinearExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb,
     &  ilo_fb, ihi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer ilo_fb, ihi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb)
      
c     local variables       
      integer i
      real s, abs_diff, dist, slope
      real one
      parameter (one = 1.0d0)

c     lower end extrapolation   
      if (bdry_location_idx .eq. 0) then
        s = sign(one,phi(ilo_fb+1))
        abs_diff = abs(phi(ilo_fb) - phi(ilo_fb+1))
        slope = s*abs_diff
        do i = ilo_gb, ilo_fb-1
          dist = ilo_fb - i
          phi(i) = phi(ilo_fb) + slope*dist
        enddo
      endif
   
c     upper end extrapolation
      if (bdry_location_idx .eq. 1) then
        s = sign(one,phi(ihi_fb-1))
        abs_diff = abs(phi(ihi_fb) - phi(ihi_fb-1))
        slope = s*abs_diff
        do i = ihi_fb+1, ihi_gb
          dist = i - ihi_fb
          phi(i) = phi(ihi_fb) + slope*dist
        enddo 
      endif

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm1dCopyExtrapolation() trivially extrapolates 1D data from the 
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
c  - if bdry_location_idx is out of the range for 1D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm1dCopyExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb,
     &  ilo_fb, ihi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer ilo_fb, ihi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb)
      
c     local variables       
      integer i

c     lower end extrapolation
      if (bdry_location_idx .eq. 0) then
        do i = ilo_gb, ilo_fb-1
          phi(i) = phi(ilo_fb)
        enddo
      endif
   
c     upper end extrapolation
      if (bdry_location_idx .eq. 1) then
        do i = ihi_fb+1, ihi_gb
          phi(i) = phi(ihi_fb)
        enddo 
      endif

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm1dHomogeneousNeumannENO1() sets the values of phi in the 
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
c  - if bdry_location_idx is out of the range for 1D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm1dHomogeneousNeumannENO1(
     &  phi,
     &  ilo_gb, ihi_gb,
     &  ilo_fb, ihi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer ilo_fb, ihi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb)
      
      call lsm1dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb,
     &        ilo_fb, ihi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm1dHomogeneousNeumannENO2() sets the values of phi in the ghostcells
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
c  - if bdry_location_idx is out of the range for 1D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm1dHomogeneousNeumannENO2(
     &  phi,
     &  ilo_gb, ihi_gb,
     &  ilo_fb, ihi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer ilo_fb, ihi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb)
      
      call lsm1dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb,
     &        ilo_fb, ihi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm1dHomogeneousNeumannENO3() sets the values of phi in the ghostcells
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
c  - if bdry_location_idx is out of the range for 1D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm1dHomogeneousNeumannENO3(
     &  phi,
     &  ilo_gb, ihi_gb,
     &  ilo_fb, ihi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer ilo_fb, ihi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb)
      
      call lsm1dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb,
     &        ilo_fb, ihi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c lsm1dHomogeneousNeumannWENO5() sets the values of phi in the ghostcells
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
c  - if bdry_location_idx is out of the range for 1D, then no 
c    ghostcell values are set 
c
c***********************************************************************
      subroutine lsm1dHomogeneousNeumannWENO5(
     &  phi,
     &  ilo_gb, ihi_gb,
     &  ilo_fb, ihi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer ilo_fb, ihi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb)
      
      call lsm1dCopyExtrapolation(phi,
     &        ilo_gb, ihi_gb,
     &        ilo_fb, ihi_fb,
     &        bdry_location_idx)

      return
      end
c } end subroutine
c***********************************************************************
