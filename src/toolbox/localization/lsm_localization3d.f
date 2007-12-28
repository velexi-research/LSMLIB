c***********************************************************************
c
c  File:        lsm_localization3d.f
c  Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
c  Revision:    $Revision: 1.8 $
c  Modified:    $Date: 2006/01/24 21:46:42 $
c  Description: F77 routines for 3D narrow-band level set calculations
c
c***********************************************************************

c***********************************************************************
c lsm3dMarkNarrowBandNeighbors() finds neighbors (iteratively +/-1 in each
c coordinate direction, so-called 6-connectivity) of the level 0 narrow band
c up to desired level.
c Narrow band of level 0 - actual narrow band voxels
c Narrow band of level L - voxels +/-L voxels away from the level 0 narrow
c                          band in x-, y- or z- direction
c                          (needed for correct computation of derivatives
c                          at the actual narrow band voxels) 
c
c Arguments:
c    phi(in):          level set functions (assumed signed distance function)
c    narrow_band(out): array with values L+1 for narrow band level L voxels
c                      and 0 otherwise
c    index_[xyz](out): array with [xyz] coordinates of narrow band voxels
c                      indices of level L narrow band stored consecutively
c    n*_index(in):     (allocated) index range of index_[xyz] arrays 
c    n_lo(in/out):     array, n_lo[L] is starting index of the level L narrow
c                      band voxels
c    n_hi(in/out):     array, n_hi[L] is ending index of the level L narrow
c                      band voxels
c    level(in):        number of narrow band levels to mark
c    *_gb (in):        index range for ghostbox
c
c Notes:
c    - narrow_band and index_[xyz] arrays assumed allocated beforehand
c    - narrow band level 0 assumed set already - appropriate voxels in marked
c     as 1 in 'narrow_band', n_lo[0], n_hi[0] known and appropriate indices
c     in that range saved in arrays index_x,index_y,index_z
c    - voxels that are outside fillbox ARE still INCLUDED in the narrow 
c     band; use lsm3dMarkNarrowBandBoundaryLayer() to distinguish the voxels
c    near volume boundary
c
c***********************************************************************
      subroutine lsm3dMarkNarrowBandNeighbors(
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb,
     &  jlo_nb_gb, jhi_nb_gb,
     &  klo_nb_gb, khi_nb_gb,
     &  index_x,
     &  index_y, 
     &  index_z,
     &  nlo_index, nhi_index,
     &  n_lo, n_hi,
     &  level)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer level
      integer n_lo(0:level), n_hi(0:level)
      
      integer i, j, k, l, m, count
      integer*1 mark
      integer start_level

      count = n_hi(0)+1
      start_level = 1
      
c     { begin loop over all levels      
      do l=start_level,level

        mark = l+1
c       { begin examine points from one level less
        do  m=n_lo(l-1),n_hi(l-1)          
	  i=index_x(m)
	  j=index_y(m)
	  k=index_z(m)

c         check upper x-coordinate neighbor	  
	  if( (i .lt. ihi_gb) .and. (narrow_band(i+1,j,k) .eq. 0) ) then
	    index_x(count) = i+1
	    index_y(count) = j
	    index_z(count) = k
	    narrow_band(i+1,j,k) = mark
	    count = count+1
	  endif

c         check lower x-coordinate neighbor	  
	  if( (i .gt. ilo_gb) .and. (narrow_band(i-1,j,k) .eq. 0) ) then
	    index_x(count) = i-1
	    index_y(count) = j
	    index_z(count) = k
	    narrow_band(i-1,j,k) = mark
	    count = count+1
	  endif

c         check upper y-coordinate neighbor	  
	  if( (j .lt. jhi_gb) .and. (narrow_band(i,j+1,k) .eq. 0) ) then
	    index_x(count) = i
	    index_y(count) = j+1
	    index_z(count) = k
	    narrow_band(i,j+1,k) = mark
	    count = count+1
	  endif

c         check lower y-coordinate neighbor	  
	  if( (j .gt. jlo_gb) .and. (narrow_band(i,j-1,k) .eq. 0) ) then
	    index_x(count) = i
	    index_y(count) = j-1
	    index_z(count) = k
	    narrow_band(i,j-1,k) = mark
	    count = count+1
	  endif

c         check upper z-coordinate neighbor	  
	  if( (k .lt. khi_gb) .and. (narrow_band(i,j,k+1) .eq. 0) ) then
	    index_x(count) = i
	    index_y(count) = j
	    index_z(count) = k+1
	    narrow_band(i,j,k+1) = mark
	    count = count+1
	  endif

c         check lower z-coordinate neighbor	  
	  if( (k .gt. klo_gb) .and. (narrow_band(i,j,k-1) .eq. 0) ) then
	    index_x(count) = i
	    index_y(count) = j
	    index_z(count) = k-1
	    narrow_band(i,j,k-1) = mark
	    count = count+1
	  endif
	  	  	  
        enddo
c       } end examine points one level less
		
	if( count .gt. ( n_hi(l-1) + 1 ) ) then 
	  n_lo(l) = n_hi(l-1) + 1
	  n_hi(l) = count - 1
	else
	  n_lo(l) = -1
	  n_hi(l) = -1
	endif
	
      enddo
c     } end loop over all levels  
  
      return
      end
c } end subroutine
c*********************************************************************** 


c***********************************************************************
c
c  lsm3dDetermineNarrowBand() finds the narrow band voxels around the zero
c  level set of the specified width. Narrow band neighbors (up to the desired
c  level) are marked as well.
c  Narrow band of level 0 - actual narrow band voxels
c  Narrow band of level L - voxels +/-L voxels in x-, y- or z- direction
c                         (needed for correct computation of derivatives
c                         at the actual narrow band voxels) 
c
c  Arguments:
c    phi(in):          level set functions (assumed signed distance function)
c    narrow_band(out): array with values L+1 for narrow band level L voxels
c                      and 0 otherwise
c    index_[xyz](out): array with [xyz] coordinates of narrow band voxels
c                      indices of level L narrow band stored consecutively
c    n*_index(in):     (allocated) index range of index_[xyz] arrays 
c    n_lo(out):        array, n_lo[L] is starting index of the level L narrow
c                      band voxels
c    n_hi(out):        array, n_hi[L] is ending index of the level L narrow
c                      band voxels  
c    level(in):        number of narrow band levels to mark
c    width(in):        narrow band width (distance to the zero level set)
c    width_inner(in):  inner narrow band width
c    index_outer(out): indices of the narrow band voxels such that 
c                      width_inner <= abs(phi) < width
c    n*_plus(out):     index range of 'index_outer'  elements for which
c                      phi values satisfy  0 < width_inner <= phi < width  
c    n*_minus(out):    index range of 'index_outer'  elements for which
c                      phi values satisfy  0> -width_inner >= phi > -width 
c    *_gb (in):        index range for ghostbox
c
c    Notes:
c    - narrow_band, index_[xyz], index_outer arrays assumed allocated beforehand
c    - index_outer stores indices of points in narrow band with positive and
c     negative phi values separately in order to be able to identify change
c     of signs for phi (i.e. when zero level set crosses into the outer layer),
c     see lsm3dCheckOuterNarrowBandLayer() 
c***********************************************************************
      subroutine lsm3dDetermineNarrowBand(
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb,
     &  jlo_nb_gb, jhi_nb_gb,
     &  klo_nb_gb, khi_nb_gb,
     &  index_x,
     &  index_y, 
     &  index_z,
     &  nlo_index, nhi_index,
     &  n_lo, n_hi,
     &  index_outer,
     &  nlo_index_outer, nhi_index_outer,
     &  nlo_outer_plus, nhi_outer_plus,
     &  nlo_outer_minus, nhi_outer_minus,
     &  width,
     &  width_inner,
     &  level)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      real width, width_inner
      integer nlo_index_outer, nhi_index_outer
      integer nlo_outer_plus, nhi_outer_plus
      integer nlo_outer_minus, nhi_outer_minus
      integer index_outer(nlo_index_outer:nhi_index_outer)
      integer level
      integer n_lo(0:level), n_hi(0:level)
      
      integer i,j,k, count, count_outer_minus, count_outer_plus
      real abs_phi_val
      integer*1  one, zero

c     get level 0 narrow band points 
      count = nlo_index
      n_lo(0) = nlo_index
      one = 1
      zero = 0
      
c     index_outer is essentially allocated beforehand to hold all gridpts
c     outer narrow band points with negative phi will be stored at the front
c     and the positive at the end of the array
      
      count_outer_minus = nlo_index_outer
      nlo_outer_minus =   nlo_index_outer
      
      count_outer_plus =  nhi_index_outer
      nhi_outer_plus =    nhi_index_outer
      
c     begin loop over grid      
      do k=klo_gb,khi_gb
	do j=jlo_gb,jhi_gb
          do i=ilo_gb,ihi_gb  
	      
	    abs_phi_val = abs(phi(i,j,k))  
	    if ( abs_phi_val .lt. width ) then
	       index_x(count) = i
	       index_y(count) = j
	       index_z(count) = k
	       narrow_band(i,j,k) = one
	       
	       if( abs_phi_val .ge. width_inner )  then
	       	       
	          if(phi(i,j,k) .le. 0d0 ) then
		      index_outer(count_outer_minus) = count
		      count_outer_minus = count_outer_minus+1
		  else
		      index_outer(count_outer_plus) = count
		      count_outer_plus = count_outer_plus-1
		  endif     
	       
	       endif	   
	       
	       count = count+1       
	    else
	       narrow_band(i,j,k) = zero  
	    endif
	      
          enddo
	enddo
      enddo
c      } end loop over grid 

      n_hi(0) = count-1
      nhi_outer_minus = count_outer_minus - 1
      nlo_outer_plus  = count_outer_plus  + 1
  
      if(level .gt. 0) then    
        call  lsm3dMarkNarrowBandNeighbors(phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
     &  index_x, index_y, index_z,
     &  nlo_index, nhi_index,
     &  n_lo, n_hi,
     &  level)
      endif
          
      return
      end     
c } end subroutine
c***********************************************************************      
  
      
c***********************************************************************
c lsm3dMarkNarrowBandBoundaryLayer() marks planes x = ilo_fb, x = ihi_fb,
c y = jlo_fb, y = jhi_fb, z = klo_fb and z = khi_fb in narrow band array 
c with the specified mark. Potentially used for identification of voxel
c layers near the volume boundary.
c
c  Arguments:
c    narrow_band(out): array with values L+1 for narrow band level L voxels
c                      and 0 otherwise
c    *_gb (in):        index range for ghostbox
c    *_fb (in):        index range for fillbox
c    mark_boundary_layer  distinctive mark (>=120 or so) for boundary layer
c                      i.e. ghost box voxels that do not belong to fill box
c
c***********************************************************************
      subroutine lsm3dMarkNarrowBandBoundaryLayer(
     &  narrow_band,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  klo_fb, khi_fb,
     &  mark_boundary_layer)
c***********************************************************************
c { begin subroutine      
      implicit none
      
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer*1 narrow_band(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      integer klo_fb, khi_fb
      integer*1 mark_boundary_layer
      
c     local variables
      integer i,j,k      
      
      do k=klo_gb,khi_gb
        do j=jlo_gb,jhi_gb
          do i=ilo_fb,ilo_fb
            narrow_band(i,j,k) = mark_boundary_layer
          enddo
        enddo
      enddo

      do k=klo_gb,khi_gb
        do j=jlo_gb,jhi_gb
          do i=ihi_fb,ihi_fb
            narrow_band(i,j,k) = mark_boundary_layer
          enddo
        enddo
      enddo

      do k=klo_gb,khi_gb
        do j=jlo_fb,jlo_fb
          do i=ilo_gb,ihi_gb
            narrow_band(i,j,k) = mark_boundary_layer
          enddo
        enddo
      enddo

      do k=klo_gb,khi_gb
        do j=jhi_fb,jhi_fb
          do i=ilo_gb,ihi_gb
            narrow_band(i,j,k) = mark_boundary_layer
          enddo
        enddo
      enddo

      do k=klo_fb,klo_fb
        do j=jlo_gb,jhi_gb
          do i=ilo_gb,ihi_gb
            narrow_band(i,j,k) = mark_boundary_layer
          enddo
        enddo
      enddo

      do k=khi_fb,khi_fb
        do j=jlo_gb,jhi_gb
          do i=ilo_gb,ihi_gb
            narrow_band(i,j,k) = mark_boundary_layer
          enddo
        enddo
      enddo

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm3dDetermineNarrowBandFromMask() initializes the narrow band voxels as 
c  either positive of negative phase of array mask. Narrow band neighbors 
c  (up to the desired level) are marked as well.
c  Narrow band of level 0 - actual narrow band voxels
c  Narrow band of level L - voxels +/-L voxels in x-, y- or z- direction
c                         (needed for correct computation of derivatives
c                         at the actual narrow band voxels) 
c
c  Arguments:
c    phi(in):          level set functions (assumed signed distance function)
c    mask(in):         level set function one phase of which will determine
c                      narrow band
c    narrow_band(out): array with values L+1 for narrow band level L voxels
c                      and 0 otherwise
c    index_[xyz](out): array with [xyz] coordinates of narrow band voxels
c                      indices of level L narrow band stored consecutively
c    nlo_index, nhi_index: (allocated) index range of index_[xyz] arrays 
c    n_lo(out):        array, n_lo[L] is starting index of the level L narrow
c                      band voxels
c    n_hi(out):        array, n_hi[L] is ending index of the level L narrow
c                      band voxels
c    level(in):        number of narrow band levels to mark
c    width(in):        narrow band width (distance to the zero level set)
c    *_gb (in):        index range for ghostbox
c    use_mask_sign     1 or -1 depending if positive or negative phase of 
c                      mask to be used to initialize the narrow band
c
c    Notes:
c    - narrow_band and index_[xyz] arrays assumed allocated beforehand
c
c***********************************************************************
      subroutine lsm3dDetermineNarrowBandFromMask(
     &  phi,
     &  mask,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb,
     &  jlo_nb_gb, jhi_nb_gb,
     &  klo_nb_gb, khi_nb_gb,
     &  index_x,
     &  index_y, 
     &  index_z,
     &  nlo_index, nhi_index,
     &  n_lo, n_hi,
     &  level,
     &  use_mask_sign)
c***********************************************************************
c { begin subroutine
      implicit none
      
        integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      real mask(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer level
      integer n_lo(0:level), n_hi(0:level)
      integer use_mask_sign     
      
c     local variables      
      integer i,j,k,count

c     get level 0 narrow band points 
      count = nlo_index
      n_lo(0) = nlo_index
      
c     begin loop over grid      
      do k=klo_gb,khi_gb
	do j=jlo_gb,jhi_gb
          do i=ilo_gb,ihi_gb  
	     
	    if( use_mask_sign .gt. 0 ) then
	      if (mask(i,j,k) .gt. 0)  then
		 index_x(count) = i
		 index_y(count) = j
		 index_z(count) = k
		 narrow_band(i,j,k) = 1
		 count = count+1
	      else
		 narrow_band(i,j,k) = 0   
	      endif
	    else
	       if (mask(i,j,k) .lt. 0) then
		 index_x(count) = i
		 index_y(count) = j
		 index_z(count) = k
		 narrow_band(i,j,k) = 1
		 count = count+1
	      else
		 narrow_band(i,j,k) = 0   
	      endif
	    endif
	      
          enddo
	enddo
      enddo
c     } end loop over grid 

      n_hi(0) = count-1

c      write(*,*) level
       
       call  lsm3dMarkNarrowBandNeighbors(phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb, klo_gb, khi_gb, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
     &  index_x, index_y, index_z,
     &  nlo_index, nhi_index,
     &  n_lo, n_hi,
     &  level)  
      
      return
      end
c } end subroutine
c***********************************************************************      


c***********************************************************************
c
c  lsm3dMultiplyCutOffLSERHSLOCAL() multiplies the right hand side of
c  the level set method equation with the cut off function described in
c  Peng at al. '99.
c    phi_t = ...
c
c  Arguments:
c    phi:               level set method function
c    lse_rhs (in/out):  right-hand of level set equation         
c    index_[xyz](out):  array with [xyz] coordinates of narrow band voxels
c                       indices
c    narrow_band(in):    array that marks voxels outside desired fillbox
c    mark_fb(in):        upper limit narrow band value for voxels in 
c                        fillbox
c    nlo_index, nhi_index: index range of index_[xyz] arrays
c    *_gb (in):         index range for ghostbox
c
c***********************************************************************
      subroutine lsm3dMultiplyCutOffLSERHSLOCAL(
     &  phi,
     &  lse_rhs,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  index_x,
     &  index_y, 
     &  index_z, 
     &  nlo_index, nhi_index,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb,
     &  jlo_nb_gb, jhi_nb_gb,
     &  klo_nb_gb, khi_nb_gb,
     &  mark_fb,
     &  beta, gamma)
c**********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      real lse_rhs(ilo_gb:ihi_gb,jlo_gb:jhi_gb,
     &             klo_gb:khi_gb)
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      real beta, gamma
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb
      
c     local variables      
      integer i,j,k,l
      real abs_phi_val, cut_off_coeff

c     { begin loop over indexed points
      do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)
	k=index_z(l)        
	
        if( narrow_band(i,j,k) .le. mark_fb ) then	

	    abs_phi_val = abs(phi(i,j,k))
	      
	    if( abs_phi_val .le. beta ) then
	       cut_off_coeff = 1
	    else if( abs_phi_val .le. gamma ) then
	       cut_off_coeff =  
     &              ((abs_phi_val - gamma)*(abs_phi_val - gamma)
     &             *(2*abs_phi_val + gamma - 3*beta)) / 
     &              ((gamma - beta)*(gamma - beta)*
     &                              (gamma - beta))		 
	    else 
	       cut_off_coeff = 0
	    endif	 
		 		 	   
            lse_rhs(i,j,k) = cut_off_coeff*lse_rhs(i,j,k)
	endif
	    
      enddo 
c     } end loop over indexed points
	
      
      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dImposeMaxPhiMaskLocal() 
c
c  Arguments:
c    phi(in):          level set functions (assumed signed distance function)
c    mask(in):         masking level set function
c    index_[xyz](out): array with [xyz] coordinates of narrow band voxels
c                      indices of level L narrow band stored consecutively
c    nlo_index, nhi_index: (allocated) index range of index_[xyz] arrays 
c    *_gb (in):        index range for ghostbox
c
c***********************************************************************
      subroutine lsm3dImposeMaskLocal(
     &  phi,
     &  mask,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  index_x,
     &  index_y, 
     &  index_z,
     &  nlo_index, nhi_index)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      real mask(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
            
c     local variables      
      integer i,j,k,l

c     { begin loop over indexed points
      do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)
	k=index_z(l)        
	    	      
	if( mask(i,j,k) .gt. phi(i,j,k) ) then
	   phi(i,j,k) = mask(i,j,k)
	endif
	   
      enddo 
c     } end loop over indexed points
    
      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dCheckOuterNarrowBandLayer() checks outer narrow band voxels for changes
c  in the sign of level set function phi values (indicates that the interface
c  represented by the zero level set has crossed over into the outer layer.
c
c  Arguments:
c    phi(in):          level set functions (assumed signed distance function)
c    index_[xyz](out): array with [xyz] coordinates of narrow band voxels
c                      indices of level L narrow band stored consecutively
c    n*_index(in):     index range of points to loop over of index_[xy] arrays 
c    *_gb (in):        index range for ghostbox
c
c***********************************************************************
      subroutine lsm3dCheckOuterNarrowBandLayer(
     &  change_sign, 
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  index_x,
     &  index_y, 
     &  index_z,
     &  nlo_index, nhi_index,
     &  index_outer,
     &  nlo_index_outer, nhi_index_outer,
     &  nlo_outer_plus, nhi_outer_plus,
     &  nlo_outer_minus, nhi_outer_minus)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer change_sign
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer nlo_index_outer, nhi_index_outer
      integer nlo_outer_plus, nhi_outer_plus
      integer nlo_outer_minus, nhi_outer_minus
      integer index_outer(nlo_index_outer:nhi_index_outer) 
         
c     local variables      
      integer i,j,k,l,m
      integer sign_plus, sign_minus
      
      sign_plus = 0
      sign_minus = 0      
      change_sign = 0
      
c     { begin loop over indexed points
      do m=nlo_outer_plus, nhi_outer_plus
        l=index_outer(m)
	i=index_x(l)
	j=index_y(l)
	k=index_z(l)
	      	    	      
	if( phi(i,j,k) .gt. 0d0 ) then
	   sign_plus = 1
	else 
	   sign_minus = 1
	endif

        if( (sign_plus .eq. 1) .and. (sign_minus .eq. 1) ) then 
	   change_sign = 1
	   exit
        endif
			     
      enddo
      
      if(change_sign .eq. 0 ) then
        sign_plus = 0
        sign_minus = 0      
        do m=nlo_outer_minus, nhi_outer_minus
          l=index_outer(m)
	  i=index_x(l)
	  j=index_y(l)
	  k=index_z(l)
	    	    	      
	  if( phi(i,j,k) .gt. 0d0 ) then
	     sign_plus = 1
	  else 
	     sign_minus = 1
	  endif

          if( (sign_plus .eq. 1) .and. (sign_minus .eq. 1) ) then 
	     change_sign = 1
	     exit
          endif 
	  	     
        enddo
      endif                
c     } end loop over indexed points     
	        
      return
      end     
c } end subroutine
c***********************************************************************
