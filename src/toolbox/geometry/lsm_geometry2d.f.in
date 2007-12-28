c***********************************************************************
c
c  File:        lsm_geometry2d.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision: 1.12 $
c  Modified:    $Date: 2006/11/01 00:25:17 $
c  Description: F77 routines for several common level set method
c               geometry calculations
c
c***********************************************************************

c***********************************************************************
c
c  lsm2dComputeUnitNormal() computes the unit normal vector to the 
c  interface from grad(phi). 
c
c  Arguments:
c    normal_* (out):  components of unit normal vector
c    phi_* (in):      components of grad(phi) 
c    dx, dy (in):     grid spacing
c    *_gb (in):       index range for ghostbox
c    *_fb (in):       index range for fillbox
c
c***********************************************************************
      subroutine lsm2dComputeUnitNormal(
     &  normal_x, normal_y,
     &  ilo_normal_gb, ihi_normal_gb,
     &  jlo_normal_gb, jhi_normal_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostboxes 
c     _fb refers to fill-box for normal data

      integer ilo_normal_gb, ihi_normal_gb
      integer jlo_normal_gb, jhi_normal_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real normal_x(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb)
      real normal_y(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real norm_grad_phi, inv_norm_grad_phi
      integer i,j
      real half
      parameter (half=0.5d0)
      real zero_tol
      parameter (zero_tol=@lsmlib_zero_tol@)

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         compute unit normal 

          norm_grad_phi = sqrt( phi_x(i,j)*phi_x(i,j)
     &                        + phi_y(i,j)*phi_y(i,j) )

          if (norm_grad_phi .ge. zero_tol) then
            inv_norm_grad_phi = 1.0d0/norm_grad_phi
            normal_x(i,j) = phi_x(i,j)*inv_norm_grad_phi
            normal_y(i,j) = phi_y(i,j)*inv_norm_grad_phi
          else
            normal_x(i,j) = 1.0d0
            normal_y(i,j) = 0.0d0
          endif

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dComputeSignedUnitNormal() computes the signed unit normal 
c  vector (sgn(phi)*normal) to the interface from grad(phi) using 
c  the following smoothed sgn function 
c
c    sgn(phi) = phi/sqrt( phi^2 + |grad(phi)|^2 * dx^2 )
c
c  Arguments:
c    normal_* (out):     components of unit normal vector
c    phi_* (in):         components of grad(phi) 
c    phi (in):           level set function
c    dx, dy (in):        grid spacing
c    *_gb (in):          index range for ghostbox
c    *_fb (in):          index range for fillbox
c
c***********************************************************************
      subroutine lsm2dComputeSignedUnitNormal(
     &  normal_x, normal_y, 
     &  ilo_normal_gb, ihi_normal_gb,
     &  jlo_normal_gb, jhi_normal_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostboxes 
c     _fb refers to fill-box for normal data

      integer ilo_normal_gb, ihi_normal_gb
      integer jlo_normal_gb, jhi_normal_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real normal_x(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb)
      real normal_y(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      real phi_cur
      real sgn_phi
      real norm_grad_phi_sq, inv_norm_grad_phi
      real dx_sq
      integer i,j
      real half
      parameter (half=0.5d0)
      real zero_tol
      parameter (zero_tol=@lsmlib_zero_tol@)

c     set value of dx_sq to be square of max{dx,dy}
      dx_sq = max(dx,dy)
      dx_sq = dx_sq*dx_sq

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         cache phi_cur
          phi_cur = phi(i,j)

c         compute sgn(phi)*normal
          if (abs(phi_cur) .gt. zero_tol) then

            norm_grad_phi_sq = phi_x(i,j)*phi_x(i,j)
     &                       + phi_y(i,j)*phi_y(i,j)

            if (norm_grad_phi_sq .ge. zero_tol) then
              sgn_phi = phi_cur
     &                / sqrt(phi_cur*phi_cur + norm_grad_phi_sq*dx_sq)

              inv_norm_grad_phi = 1.d0/sqrt(norm_grad_phi_sq)

              normal_x(i,j) = sgn_phi*phi_x(i,j)*inv_norm_grad_phi
              normal_y(i,j) = sgn_phi*phi_y(i,j)*inv_norm_grad_phi
            else
              normal_x(i,j) = 1.0d0
              normal_y(i,j) = 0.0d0
            endif

          else

            normal_x(i,j) = 0.0d0
            normal_y(i,j) = 0.0d0

          endif

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dAreaRegionPhiLessThanZero() computes the area of the 
c  region where the level set function is less than 0.  
c
c  Arguments:
c    area (out):            area of the region where phi < 0
c    phi (in):              level set function
c    dx, dy (in):           grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm2dAreaRegionPhiLessThanZero(
     &  area,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy, 
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real area

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real dA
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dA = dx * dy 
      dA = dx * dy 

c     initialize area to zero
      area = 0.0d0

c      loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon

              if (phi_cur .lt. -epsilon) then
                area = area + dA
              elseif (phi_cur .lt. epsilon) then
                one_minus_H = 0.5d0*( 1 - phi_cur_over_epsilon
     &                                  - one_over_pi
     &                                  * sin(pi*phi_cur_over_epsilon) )
                area = area + one_minus_H*dA
              endif
    
          enddo
        enddo
c       } end loop over grid
      
      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dAreaRegionPhiGreaterThanZero() computes the area of the 
c  region where the level set function is greater than 0.  
c
c  Arguments:
c    area (out):            area of the region where phi > 0
c    phi (in):              level set function
c    dx, dy (in):           grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm2dAreaRegionPhiGreaterThanZero(
     &  area,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real area

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real dA
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dA = dx * dy 
      dA = dx * dy 

c     initialize area to zero
      area = 0.0d0

c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
  
              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon
   
              if (phi_cur .gt. epsilon) then
                area = area + dA
              elseif (phi_cur .gt. -epsilon) then
                H = 0.5d0*( 1 + phi_cur_over_epsilon 
     &                      + one_over_pi*sin(pi*phi_cur_over_epsilon) )
                area = area + H*dA
              endif

          enddo
        enddo
c       } end loop over grid
      
      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dPerimeterZeroLevelSet() computes the perimeter of the zero 
c  level set. 
c
c  Arguments:
c    perimeter (out):       perimeter of curve defined by the zero level
c                           set
c    phi (in):              level set function
c    phi_* (in):            components of grad(phi)
c    dx, dy (in):           grid spacing
c    epsilon (in):          width of numerical smoothing to use for
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm2dPerimeterZeroLevelSet(
     &  perimeter,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy, 
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real perimeter

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)  
      real dx,dy
      real epsilon
      real one_over_epsilon
      integer i,j
      real phi_cur
      real delta
      real norm_grad_phi
      real dA
      real pi
      parameter (pi=3.14159265358979323846d0)
      

c     compute dA = dx * dy
      dA = dx * dy

c     compute one_over_epsilon
      one_over_epsilon = 1.d0/epsilon

c     initialize perimeter to zero
      perimeter = 0.0d0

c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
             
             phi_cur = phi(i,j)
  
             if (abs(phi_cur) .lt. epsilon) then
                delta = 0.5d0*one_over_epsilon
     &                * ( 1+cos(pi*phi_cur*one_over_epsilon) ) 

                norm_grad_phi = sqrt(
     &              phi_x(i,j)*phi_x(i,j)
     &            + phi_y(i,j)*phi_y(i,j) )

                perimeter = perimeter + delta*norm_grad_phi*dA
             endif

         enddo
        enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dAreaRegionPhiLessThanZeroControlVolume() computes the area
c  of the region of the computational domain where the level set
c  function is less than 0.  The computational domain contains only
c  those cells that are included by the control volume data.
c
c  Arguments:
c    area (out):            area of the region where phi < 0
c    phi (in):              level set function
c    control_vol (in):      control volume data (used to exclude cells
c                           from the integral calculation)
c    control_vol_sgn (in):  1 (-1) if positive (negative) control volume
c                           points should be used
c    dx, dy (in):           grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm2dAreaRegionPhiLessThanZeroControlVolume(
     &  area,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy, 
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real area

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real dA
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dA = dx * dy 
      dA = dx * dy 

c     initialize area to zero
      area = 0.0d0

      if (control_vol_sgn .gt. 0) then
        
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
  
c           only include cell in max norm calculation if it has a 
c           positive control volume
            if (control_vol(i,j) .gt. 0.d0) then

              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon

              if (phi_cur .lt. -epsilon) then
                area = area + dA
              elseif (phi_cur .lt. epsilon) then
                one_minus_H = 0.5d0*( 1 - phi_cur_over_epsilon
     &                                  - one_over_pi
     &                                  * sin(pi*phi_cur_over_epsilon) )
                area = area + one_minus_H*dA
              endif

            endif
      
          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
  
c           only include cell in max norm calculation if it has a 
c           negative control volume
            if (control_vol(i,j) .lt. 0.d0) then

              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon

              if (phi_cur .lt. -epsilon) then
                area = area + dA
              elseif (phi_cur .lt. epsilon) then
                one_minus_H = 0.5d0*( 1 - phi_cur_over_epsilon
     &                                  - one_over_pi
     &                                  * sin(pi*phi_cur_over_epsilon) )
                area = area + one_minus_H*dA
              endif

            endif
      
          enddo
        enddo
c       } end loop over grid

      endif
      
      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dAreaRegionPhiGreaterThanZeroControlVolume() computes the area 
c  of the region of the computational domain where the level set
c  function is greater than 0.  The computational domain contains only
c  those cells that are included by the control volume data.
c
c  Arguments:
c    area (out):            area of the region where phi > 0
c    phi (in):              level set function
c    control_vol (in):      control volume data (used to exclude cells
c                           from the integral calculation)
c    control_vol_sgn (in):  1 (-1) if positive (negative) control volume
c                           points should be used
c    dx, dy (in):           grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm2dAreaRegionPhiGreaterThanZeroControlVolume(
     &  area,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real area

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real dA
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dA = dx * dy 
      dA = dx * dy 

c     initialize area to zero
      area = 0.0d0


      if (control_vol_sgn .gt. 0) then

c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
   
c           only include cell in max norm calculation if it has a 
c           positive control volume
            if (control_vol(i,j) .gt. 0.d0) then
  
              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon
   
              if (phi_cur .gt. epsilon) then
                area = area + dA
              elseif (phi_cur .gt. -epsilon) then
                H = 0.5d0*( 1 + phi_cur_over_epsilon 
     &                      + one_over_pi*sin(pi*phi_cur_over_epsilon) )
                area = area + H*dA
              endif

            endif

          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
   
c           only include cell in max norm calculation if it has a 
c           negative control volume
            if (control_vol(i,j) .lt. 0.d0) then
  
              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon
   
              if (phi_cur .gt. epsilon) then
                area = area + dA
              elseif (phi_cur .gt. -epsilon) then
                H = 0.5d0*( 1 + phi_cur_over_epsilon 
     &                      + one_over_pi*sin(pi*phi_cur_over_epsilon) )
                area = area + H*dA
              endif

            endif

          enddo
        enddo
c       } end loop over grid
      endif
      
      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dPerimeterZeroLevelSetControlVolume() computes the perimeter 
c  of the zero level set within the computational domain.  The 
c  computational domain contains only those cells that are included 
c  by the control volume data.
c
c  Arguments:
c    perimeter (out):       perimeter of curve defined by the zero level
c                           set
c    phi (in):              level set function
c    phi_* (in):            components of grad(phi)
c    control_vol (in):      control volume data (used to exclude cells
c                           from the integral calculation)
c    control_vol_sgn (in):  1 (-1) if positive (negative) control volume
c                           points should be used
c    dx, dy (in):           grid spacing
c    epsilon (in):          width of numerical smoothing to use for
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm2dPerimeterZeroLevelSetControlVolume(
     &  perimeter,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy, 
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real perimeter

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy
      real epsilon
      real one_over_epsilon
      integer i,j
      real phi_cur
      real delta
      real norm_grad_phi
      real dA
      real pi
      parameter (pi=3.14159265358979323846d0)
      

c     compute dA = dx * dy
      dA = dx * dy

c     compute one_over_epsilon
      one_over_epsilon = 1.d0/epsilon

c     initialize perimeter to zero
      perimeter = 0.0d0

      if (control_vol_sgn .gt. 0) then
  
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in max norm calculation if it has a 
c           positive control volume
            if (control_vol(i,j) .gt. 0.d0) then

              phi_cur = phi(i,j)
  
              if (abs(phi_cur) .lt. epsilon) then
                delta = 0.5d0*one_over_epsilon
     &                * ( 1+cos(pi*phi_cur*one_over_epsilon) ) 

                norm_grad_phi = sqrt(
     &              phi_x(i,j)*phi_x(i,j)
     &            + phi_y(i,j)*phi_y(i,j) )

                perimeter = perimeter + delta*norm_grad_phi*dA
              endif

            endif
        
          enddo
        enddo
c       } end loop over grid

      else
      
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in max norm calculation if it has a 
c           negative control volume
            if (control_vol(i,j) .lt. 0.d0) then

              phi_cur = phi(i,j)
  
              if (abs(phi_cur) .lt. epsilon) then
                delta = 0.5d0*one_over_epsilon
     &                * ( 1+cos(pi*phi_cur*one_over_epsilon) ) 

                norm_grad_phi = sqrt(
     &              phi_x(i,j)*phi_x(i,j)
     &            + phi_y(i,j)*phi_y(i,j) )

                perimeter = perimeter + delta*norm_grad_phi*dA
              endif

            endif
        
          enddo
        enddo
c       } end loop over grid

      endif
      
      return
      end
c } end subroutine
c***********************************************************************
