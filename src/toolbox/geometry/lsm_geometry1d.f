c***********************************************************************
c
c  File:        lsm_geometry1d.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision: 1.12 $
c  Modified:    $Date: 2006/11/01 00:25:17 $
c  Description: F77 routines for several common level set method
c               geometry calculations
c
c***********************************************************************

c***********************************************************************
c
c  lsm1dComputeUnitNormal() computes the unit normal vector to the 
c  interface from grad(phi).
c
c  Arguments:
c    normal (out):  unit normal vector
c    phi_* (in):    components of grad(phi) 
c    dx (in):       grid spacing
c    *_gb (in):     index range for ghostbox
c    *_fb (in):     index range for fillbox
c
c***********************************************************************
      subroutine lsm1dComputeUnitNormal(
     &  normal,
     &  ilo_normal_gb, ihi_normal_gb,
     &  phi_x,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  ilo_fb, ihi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fill-box for normal data

      integer ilo_normal_gb, ihi_normal_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer ilo_fb, ihi_fb
      real normal(ilo_normal_gb:ihi_normal_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb)
      integer i
      real abs_phi_x_cur
      real half
      parameter (half=0.5d0)
      real zero_tol
      parameter (zero_tol=@lsmlib_zero_tol@)

c     compute unit normal by taking average of upwind and downwind 
c     derivatives and normalizing
c     { begin loop over grid
      do i=ilo_fb,ihi_fb

        abs_phi_x_cur = abs(phi_x(i))

        if (abs_phi_x_cur .ge. zero_tol) then
          normal(i) = normal(i)/abs_phi_x_cur
        else
          normal(i) = 1.0d0
        endif

      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm1dComputeSignedUnitNormal() computes the signed unit normal 
c  vector (sgn(phi)*normal) to the interface from grad(phi) using 
c  the following smoothed sgn function 
c
c    sgn(phi) = phi/sqrt( phi^2 + |grad(phi)|^2 * dx^2 )
c
c  This expression avoids division by zero in computing the unit
c  normal vector.
c
c  Arguments:
c    normal_* (out):     components of unit normal vector
c    phi_* (in):         components of grad(phi) 
c    phi (in):           level set function
c    dx (in):            grid spacing
c    *_gb (in):          index range for ghostbox
c    *_fb (in):          index range for fillbox
c
c***********************************************************************
      subroutine lsm1dComputeSignedUnitNormal(
     &  normal_x, 
     &  ilo_normal_gb, ihi_normal_gb,
     &  phi_x, 
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  ilo_fb, ihi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostboxes 
c     _fb refers to fill-box for normal data

      integer ilo_normal_gb, ihi_normal_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_fb, ihi_fb
      real normal_x(ilo_normal_gb:ihi_normal_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb)
      real dx
      real phi_cur
      real phi_x_cur
      real sgn_phi
      real norm_grad_phi_sq
      real dx_sq
      integer i
      real half
      parameter (half=0.5d0)
      real zero_tol
      parameter (zero_tol=@lsmlib_zero_tol@)

c     compute dx_sq
      dx_sq = dx**2

c     { begin loop over grid
      do i=ilo_fb,ihi_fb

c       cache phi_cur
        phi_cur = phi(i)
        phi_x_cur = phi_x(i)

c       compute sgn(phi)*normal
        if (abs(phi_cur) .gt. zero_tol) then
          norm_grad_phi_sq = phi_x_cur*phi_x_cur

          if (norm_grad_phi_sq .ge. zero_tol) then
            sgn_phi = phi_cur
     &              / sqrt(phi_cur*phi_cur + norm_grad_phi_sq*dx_sq)

            normal_x(i) = sgn_phi*phi_x_cur/abs(phi_x_cur)

          else

            normal_x(i) = 1.0d0

          endif 

        else

          normal_x(i) = 0.0d0

        endif

      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm1dLengthRegionPhiLessThanZero() computes the length of the 
c  region where the level set function is less than 0.  
c
c  Arguments:
c    length (out):          length of the region where phi < 0
c    phi (in):              level set function
c    dx (in):               grid spacing
c    epsilon (in):          width of numerical smoothing to use for
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm1dLengthRegionPhiLessThanZero(
     &  length,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  dx,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real length

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_ib, ihi_ib
      real phi(ilo_phi_gb:ihi_phi_gb)
      real dx
      real epsilon
      integer i
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     initialize length to zero
      length = 0.0d0

c     loop over included cells {
      do i=ilo_ib,ihi_ib
          phi_cur = phi(i)
          phi_cur_over_epsilon = phi_cur/epsilon

          if (phi_cur .lt. -epsilon) then
            length = length + dx
          elseif (phi_cur .lt. epsilon) then
            one_minus_H = 0.5d0*( 1 - phi_cur_over_epsilon
     &                              - one_over_pi
     &                              * sin(pi*phi_cur_over_epsilon) )
            length = length + one_minus_H*dx
          endif
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm1dLengthRegionPhiGreaterThanZero() computes the length of the 
c  region where the level set function is greater than 0.  
c
c  Arguments:
c    length (out):          length of the region where phi > 0
c    phi (in):              level set function
c    dx (in):               grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm1dLengthRegionPhiGreaterThanZero(
     &  length,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  dx, 
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real length

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_ib, ihi_ib
      real phi(ilo_phi_gb:ihi_phi_gb)
      real dx
      real epsilon
      integer i
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     initialize length to zero
      length = 0.0d0

c     loop over included cells {
      do i=ilo_ib,ihi_ib
          phi_cur = phi(i)
          phi_cur_over_epsilon = phi_cur/epsilon

          if (phi_cur .gt. epsilon) then
            length = length + dx
          elseif (phi_cur .gt. -epsilon) then
            H = 0.5d0*( 1 + phi_cur_over_epsilon 
     &                    + one_over_pi*sin(pi*phi_cur_over_epsilon) )
            length = length + H*dx
          endif
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm1dSizeZeroLevelSet() computes the size of the zero level set. 
c
c  Arguments:
c    size (out):            number of points contained in the zero level
c                           set
c    phi (in):              level set function
c    phi_* (in):            components of grad(phi)
c    dx (in):               grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm1dSizeZeroLevelSet(
     &  size,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  phi_x, 
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  ilo_ib, ihi_ib,
     &  dx, 
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real size

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer ilo_ib, ihi_ib
      real phi(ilo_phi_gb:ihi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb)
      real dx
      real epsilon
      real one_over_epsilon
      integer i
      real phi_cur
      real delta
      real pi
      parameter (pi=3.14159265358979323846d0)
      

c     initialize size to zero
      size = 0.0d0

c     compute one_over_epsilon
      one_over_epsilon = 1.d0/epsilon

c     loop over included cells {
      do i=ilo_ib,ihi_ib
          phi_cur = phi(i)

          if (abs(phi_cur) .lt. epsilon) then
            delta = 0.5d0*one_over_epsilon
     &            * ( 1+cos(pi*phi_cur*one_over_epsilon) ) 

            size = size + delta*abs(phi_x(i))*dx
          endif
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm1dLengthRegionPhiLessThanZeroControlVolume() computes the length 
c  of the region of the computational domain where the level set 
c  function is less than 0.  The computational domain contains only
c  those cells that are included by the control volume data.
c
c  Arguments:
c    length (out):          length of the region where phi < 0
c    phi (in):              level set function
c    control_vol (in):      control volume data (used to exclude cells
c                           from the integral calculation)
c    control_vol_sgn (in):  1 (-1) if positive (negative) control volume
c                           points should be used
c    dx (in):               grid spacing
c    epsilon (in):          width of numerical smoothing to use for
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm1dLengthRegionPhiLessThanZeroControlVolume(
     &  length,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  dx,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real length

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer ilo_ib, ihi_ib
      real phi(ilo_phi_gb:ihi_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb)
      integer control_vol_sgn
      real dx
      real epsilon
      integer i
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     initialize length to zero
      length = 0.0d0

c     loop over included cells {
      do i=ilo_ib,ihi_ib
  
c       only include cell in calculation if it has a 
c       control volume of desired sign
        if (control_vol_sgn*control_vol(i) .gt. 0.d0) then

          phi_cur = phi(i)
          phi_cur_over_epsilon = phi_cur/epsilon

          if (phi_cur .lt. -epsilon) then
            length = length + dx
          elseif (phi_cur .lt. epsilon) then
            one_minus_H = 0.5d0*( 1 - phi_cur_over_epsilon
     &                              - one_over_pi
     &                              * sin(pi*phi_cur_over_epsilon) )
            length = length + one_minus_H*dx
          endif

        endif
      
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm1dLengthRegionPhiGreaterThanZeroControlVolume() computes the 
c  length of the region of the computational domain where the level 
c  set function is greater than 0.  The computational domain contains 
c  only those cells that are included by the control volume data.
c
c  Arguments:
c    length (out):          length of the region where phi > 0
c    phi (in):              level set function
c    control_vol (in):      control volume data (used to exclude cells
c                           from the integral calculation)
c    control_vol_sgn (in):  1 (-1) if positive (negative) control volume
c                           points should be used
c    dx (in):               grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm1dLengthRegionPhiGreaterThanZeroControlVolume(
     &  length,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  dx, 
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real length

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer ilo_ib, ihi_ib
      real phi(ilo_phi_gb:ihi_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb)
      integer control_vol_sgn
      real dx
      real epsilon
      integer i
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     initialize length to zero
      length = 0.0d0

c     loop over included cells {
      do i=ilo_ib,ihi_ib
 
c       only include cell in calculation if it has a 
c       control volume of desired sign
        if (control_vol_sgn*control_vol(i) .gt. 0.d0) then

          phi_cur = phi(i)
          phi_cur_over_epsilon = phi_cur/epsilon

          if (phi_cur .gt. epsilon) then
            length = length + dx
          elseif (phi_cur .gt. -epsilon) then
            H = 0.5d0*( 1 + phi_cur_over_epsilon 
     &                    + one_over_pi*sin(pi*phi_cur_over_epsilon) )
            length = length + H*dx
          endif

        endif

      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm1dSizeZeroLevelSetControlVolume() computes the size of the zero 
c  level set within the computational domain.  The computational domain 
c  contains only those cells that are included by the control volume 
c  data.
c
c  Arguments:
c    size (out):            number of points contained in the zero level
c                           set
c    phi (in):              level set function
c    phi_* (in):            components of grad(phi)
c    control_vol (in):      control volume data (used to exclude cells
c                           from the integral calculation)
c    control_vol_sgn (in):  1 (-1) if positive (negative) control volume
c                           points should be used
c    dx (in):               grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm1dSizeZeroLevelSetControlVolume(
     &  size,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  phi_x, 
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  control_vol,
     &  control_vol_sgn,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  ilo_ib, ihi_ib,
     &  dx, 
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real size

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer ilo_ib, ihi_ib
      real phi(ilo_phi_gb:ihi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb)
      integer control_vol_sgn
      real dx
      real epsilon
      real one_over_epsilon
      integer i
      real phi_cur
      real delta
      real pi
      parameter (pi=3.14159265358979323846d0)
      

c     initialize size to zero
      size = 0.0d0

c     compute one_over_epsilon
      one_over_epsilon = 1.d0/epsilon

c     loop over included cells {
      do i=ilo_ib,ihi_ib

c       only include cell in calculation if it has a 
c       control volume of desired sign
        if (control_vol_sgn*control_vol(i) .gt. 0.d0) then

          phi_cur = phi(i)

          if (abs(phi_cur) .lt. epsilon) then
            delta = 0.5d0*one_over_epsilon
     &            * ( 1+cos(pi*phi_cur*one_over_epsilon) ) 

            size = size + delta*abs(phi_x(i))*dx
          endif

        endif

      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
