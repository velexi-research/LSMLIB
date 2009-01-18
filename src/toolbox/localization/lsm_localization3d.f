c***********************************************************************
c
c  File:        lsm_localization3d.f
c  Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
c  Revision:    $Revision: 1.8 $
c  Modified:    $Date$
c  Description: F77 routines for 3D narrow-band level set calculations
c
c***********************************************************************

c***********************************************************************
c FUNCTION INTERNAL TO THIS FILE- DO NOT DELETE THIS DOCUMENTATION!
c
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
      
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
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
      
c     { begin loop over all narrow band levels      
      do l=start_level,level

        mark = l+1
c       { begin examine points from one level less
        do  m=n_lo(l-1),n_hi(l-1)          
	  i=index_x(m)
	  j=index_y(m)
	  k=index_z(m)

c         check upper x-coordinate neighbor 
	  if( (i .lt. ihi_nb_gb) .and. 
     &        (narrow_band(i+1,j,k) .eq. 0) ) then
	    index_x(count) = i+1
	    index_y(count) = j
	    index_z(count) = k
	    narrow_band(i+1,j,k) = mark
	    count = count+1
	  endif

c         check lower x-coordinate neighbor  
	  if( (i .gt. ilo_nb_gb) .and. 
     &        (narrow_band(i-1,j,k) .eq. 0) ) then
	    index_x(count) = i-1
	    index_y(count) = j
	    index_z(count) = k
	    narrow_band(i-1,j,k) = mark
	    count = count+1
	  endif

c         check upper y-coordinate neighbor
	  if( (j .lt. jhi_nb_gb) .and. 
     &       (narrow_band(i,j+1,k) .eq. 0) ) then
	    index_x(count) = i
	    index_y(count) = j+1
	    index_z(count) = k
	    narrow_band(i,j+1,k) = mark
	    count = count+1
	  endif

c         check lower y-coordinate neighbor  
	  if( (j .gt. jlo_nb_gb) .and.
     &        (narrow_band(i,j-1,k) .eq. 0) ) then
	    index_x(count) = i
	    index_y(count) = j-1
	    index_z(count) = k
	    narrow_band(i,j-1,k) = mark
	    count = count+1
	  endif

c         check upper z-coordinate neighbor  
	  if( (k .lt. khi_nb_gb) .and. 
     &        (narrow_band(i,j,k+1) .eq. 0) ) then
	    index_x(count) = i
	    index_y(count) = j
	    index_z(count) = k+1
	    narrow_band(i,j,k+1) = mark
	    count = count+1
	  endif

c         check lower z-coordinate neighbor  
	  if( (k .gt. klo_nb_gb) .and.
     &        (narrow_band(i,j,k-1) .eq. 0) ) then
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

      if( count .gt. nlo_index ) then
      
         n_hi(0) = count-1
         nhi_outer_minus = count_outer_minus - 1
         nlo_outer_plus  = count_outer_plus  + 1
  
         call  lsm3dMarkNarrowBandNeighbors(
     &   narrow_band,
     &   ilo_nb_gb, ihi_nb_gb, jlo_nb_gb, jhi_nb_gb, 
     &   klo_nb_gb, khi_nb_gb,
     &   index_x, index_y, index_z,
     &   nlo_index, nhi_index,
     &   n_lo, n_hi,
     &   level)
         
      else
        n_hi(0) = count;
        nhi_outer_minus = count_outer_minus
        nlo_outer_plus  = count_outer_plus 
      endif
       
          
      return
      end     
c } end subroutine
c***********************************************************************      
 

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
      subroutine lsm3dDetermineNarrowBandFromMask(
     &  mask,
     &  ilo_mask_gb, ihi_mask_gb,
     &  jlo_mask_gb, jhi_mask_gb,
     &  klo_mask_gb, khi_mask_gb,
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
      
      integer ilo_mask_gb, ihi_mask_gb
      integer jlo_mask_gb, jhi_mask_gb
      integer klo_mask_gb, khi_mask_gb
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      real mask(ilo_mask_gb:ihi_mask_gb,
     &          jlo_mask_gb:jhi_mask_gb,
     &          klo_mask_gb:khi_mask_gb)
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
      do k=klo_mask_gb,khi_mask_gb
	do j=jlo_mask_gb,jhi_mask_gb
          do i=ilo_mask_gb,ihi_mask_gb  
	     
	    if( use_mask_sign .gt. 0 ) then
	      if (mask(i,j,k) .ge. 0)  then
		 index_x(count) = i
		 index_y(count) = j
		 index_z(count) = k
		 narrow_band(i,j,k) = 1
		 count = count+1
	      else
		 narrow_band(i,j,k) = 0   
	      endif
	    else
	       if (mask(i,j,k) .le. 0) then
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

       call  lsm3dMarkNarrowBandNeighbors(
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
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb
      
c     local variables      
      integer i,j,k,l
      real abs_phi_val, cut_off_coeff
      real gb_const1, gb_const2, temp
      real tol
      parameter (tol=1.d-13)
      
      gb_const1 = gamma - 3*beta;
      gb_const2 = (gamma - beta);
      gb_const2 = gb_const2*gb_const2*gb_const2;
      
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
	       temp = (abs_phi_val - gamma);
	       cut_off_coeff = ( temp * temp
     &               *(2*abs_phi_val + gb_const1) ) / gb_const2
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



c***********************************************************************
      subroutine lsm3dImposeMaskLocal(
     &  dest,
     &  ilo_dest_gb, ihi_dest_gb,
     &  jlo_dest_gb, jhi_dest_gb,
     &  klo_dest_gb, khi_dest_gb,
     &  src,
     &  ilo_src_gb, ihi_src_gb,
     &  jlo_src_gb, jhi_src_gb,
     &  klo_src_gb, khi_src_gb,
     &  mask,
     &  ilo_mask_gb, ihi_mask_gb,
     &  jlo_mask_gb, jhi_mask_gb,
     &  klo_mask_gb, khi_mask_gb,
     &  index_x,
     &  index_y, 
     &  index_z,
     &  nlo_index, nhi_index)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_dest_gb, ihi_dest_gb
      integer jlo_dest_gb, jhi_dest_gb
      integer klo_dest_gb, khi_dest_gb
      integer ilo_src_gb,  ihi_src_gb
      integer jlo_src_gb,  jhi_src_gb
      integer klo_src_gb,  khi_src_gb
      integer ilo_mask_gb, ihi_mask_gb
      integer jlo_mask_gb, jhi_mask_gb
      integer klo_mask_gb, khi_mask_gb
      real dest(ilo_dest_gb:ihi_dest_gb,
     &                      jlo_dest_gb:jhi_dest_gb,
     &                      klo_dest_gb:khi_dest_gb)
      real src(ilo_src_gb:ihi_src_gb,
     &                     jlo_src_gb:jhi_src_gb,
     &                     klo_src_gb:khi_src_gb)
      real mask(ilo_mask_gb:ihi_mask_gb,
     &                      jlo_mask_gb:jhi_mask_gb,     
     &                      klo_mask_gb:khi_mask_gb)
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
	   
	dest(i,j,k) = max(mask(i,j,k),src(i,j,k))
	   
      enddo 
c     } end loop over indexed points
    
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dCopyDataLocal(
     &  dest,
     &  ilo_dest_gb, ihi_dest_gb,
     &  jlo_dest_gb, jhi_dest_gb,
     &  klo_dest_gb, khi_dest_gb,
     &  src,
     &  ilo_src_gb, ihi_src_gb,
     &  jlo_src_gb, jhi_src_gb,
     &  klo_src_gb, khi_src_gb,
     &  index_x,
     &  index_y,
     &  index_z,
     &  nlo_index, nhi_index)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox
      integer ilo_dest_gb, ihi_dest_gb
      integer jlo_dest_gb, jhi_dest_gb
      integer klo_dest_gb, khi_dest_gb
      integer ilo_src_gb, ihi_src_gb
      integer jlo_src_gb, jhi_src_gb
      integer klo_src_gb, khi_src_gb
      real dest(ilo_dest_gb:ihi_dest_gb,
     &                      jlo_dest_gb:jhi_dest_gb,
     &                      klo_dest_gb:khi_dest_gb)
      real src(ilo_src_gb:ihi_src_gb,
     &                     jlo_src_gb:jhi_src_gb,
     &                     klo_src_gb:khi_src_gb)
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
         dest(i,j,k) = src(i,j,k)

        enddo
c       } end loop over indexed points
      
      return
      end
c } end subroutine
c***********************************************************************
