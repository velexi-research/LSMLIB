c***********************************************************************
c
c  File:        lsm_level_set_evolution2d_local.f
c  Copyright:   (c) 2005-2008 Masa Prodanovic and Kevin T. Chu
c  Revision:    $Revision: 1.6 $
c  Modified:    $Date: 2007/05/06 23:47:29 $
c  Description: F77 subroutines for 2D, narrow-band level set 
c               evolution equation
c
c***********************************************************************

c***********************************************************************
c
c  lsm2dZeroOutLevelSetEqnRHSLOCAL() zeros out the right-hand side of the
c  level set equation when it is written in the form:
c
c    phi_t = ...
c
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    lse_rhs (in/out):  right-hand of level set equation
c    *_gb (in):         index range for ghostbox
c    index_[xy](in):    [xy] coordinates of local (narrow band) points
c    n*_index(in):      index range of points to loop over in index_*
c
c***********************************************************************
      subroutine lsm2dZeroOutLevelSetEqnRHSLOCAL(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb,
     &  index_x,
     &  index_y, 
     &  nlo_index, nhi_index)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)

c     local variables      
      integer i,j,l

c     { begin loop over indexed points
      do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)
	
	lse_rhs(i,j) = 0.d0
      enddo 
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dAddAdvectionTermToLSERHSLOCAL adds the contribution of an advection
c  term (external vector velocity field) to the right-hand side of the 
c  level set equation when it is written in the form:
c
c    phi_t = -vel dot grad(phi) + ...
c
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    lse_rhs (in/out):  right-hand
c    phi_* (in):        components of grad(phi) at t = t_cur
c    vel_* (in):        components of velocity at t = t_cur
c    *_gb (in):         index range for ghostbox
c    index_[xy](in):    [xy] coordinates of local (narrow band) points
c    n*_index(in):      index range of points to loop over in index_*
c    narrow_band(in):   array that marks voxels outside desired fillbox
c    mark_fb(in):       upper limit narrow band value for voxels in 
c                       fillbox
c
c***********************************************************************
      subroutine lsm2dAddAdvectionTermToLSERHSLOCAL(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  vel_x, vel_y,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  index_x,
     &  index_y, 
     &  nlo_index, nhi_index,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)      
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb
      
      integer i,j,l 

c     { begin loop over indexed points
        do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)

c       include only points with appropriate mark        
	  if( narrow_band(i,j) .le. mark_fb ) then
	      
            lse_rhs(i,j) = lse_rhs(i,j) - ( vel_x(i,j)*phi_x(i,j)
     &                                  + vel_y(i,j)*phi_y(i,j) )
          endif
	  
        enddo
c     } end loop over indexed points

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm2dAddNormalVelTermToLSERHSLOCAL() adds the contribution of a normal 
c  (scalar) velocity term to the right-hand side of the level set 
c  equation when it is written in the form:
c
c    phi_t = -V_n |grad(phi)| + ...
c
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    lse_rhs (in/out):  right-hand of level set equation
c    phi_*_plus (in):   components of forward approx to grad(phi) at 
c                       t = t_cur
c    phi_*_minus (in):  components of backward approx to grad(phi) at 
c                       t = t_cur
c    vel_n (in):        normal velocity at t = t_cur
c    *_gb (in):         index range for ghostbox
c    index_[xy](in):    [xy] coordinates of local (narrow band) points
c    n*_index(in):      index range of points to loop over in index_*
c    narrow_band(in):    array that marks voxels outside desired fillbox
c    mark_fb(in):        upper limit narrow band value for voxels in 
c                        fillbox
c
c***********************************************************************
      subroutine lsm2dAddNormalVelTermToLSERHSLOCAL(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  vel_n,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  index_x,
     &  index_y, 
     &  nlo_index, nhi_index,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real vel_n(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb  
      integer i,j,l
      real vel_n_cur
      real norm_grad_phi_sq
      real tol
      parameter (tol=1.d-13)

c     { begin loop over indexed points
      do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)

        if( narrow_band(i,j) .le. mark_fb ) then
          vel_n_cur = vel_n(i,j)

c         { begin Godunov selection of grad_phi

          if (vel_n_cur .gt. 0.d0) then
            norm_grad_phi_sq = max(max(phi_x_minus(i,j),0.d0)**2,
     &                             min(phi_x_plus(i,j),0.d0)**2 )
     &                       + max(max(phi_y_minus(i,j),0.d0)**2,
     &                             min(phi_y_plus(i,j),0.d0)**2 )
          else
            norm_grad_phi_sq = max(min(phi_x_minus(i,j),0.d0)**2,
     &                             max(phi_x_plus(i,j),0.d0)**2 )
     &                       + max(min(phi_y_minus(i,j),0.d0)**2,
     &                             max(phi_y_plus(i,j),0.d0)**2 )
          endif

c         } end Godunov selection of grad_phi


c         compute contribution to lse_rhs(i,j) 
          if (abs(vel_n_cur) .ge. tol) then
            lse_rhs(i,j) = lse_rhs(i,j) 
     &                   - vel_n_cur*sqrt(norm_grad_phi_sq)
          endif
     
        endif 
      enddo 
c     } end loop over indexed points

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dAddConstNormalVelTermToLSERHSLOCAL() adds the contribution of a normal 
c  (scalar) velocity term to the right-hand side of the level set 
c  equation when it is written in the form:
c
c    phi_t = -V_n |grad(phi)| + ...
c
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    lse_rhs (in/out):  right-hand of level set equation
c    phi_*_plus (in):   components of forward approx to grad(phi) at 
c                       t = t_cur
c    phi_*_minus (in):  components of backward approx to grad(phi) at 
c                       t = t_cur
c    vel_n (in):        scalar normal velocity at t = t_cur
c    *_gb (in):         index range for ghostbox
c    index_[xy](in):    [xy] coordinates of local (narrow band) points
c    n*_index(in):      index range of points to loop over in index_*
c    narrow_band(in):    array that marks voxels outside desired fillbox
c    mark_fb(in):        upper limit narrow band value for voxels in 
c                        fillbox
c
c***********************************************************************
      subroutine lsm2dAddConstNormalVelTermToLSERHSLOCAL(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  vel_n,
     &  index_x,
     &  index_y, 
     &  nlo_index, nhi_index,  
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real vel_n
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb  
      integer i,j,l
      real norm_grad_phi_sq
      real tol
      parameter (tol=1.d-13)

c     { begin loop over indexed points
        do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)
	
c       { begin Godunov selection of grad_phi

       if( narrow_band(i,j) .le. mark_fb ) then	
     
        if (vel_n .gt. 0.d0) then
          norm_grad_phi_sq = max(max(phi_x_minus(i,j),0.d0)**2,
     &                           min(phi_x_plus(i,j),0.d0)**2 )
     &                     + max(max(phi_y_minus(i,j),0.d0)**2,
     &                           min(phi_y_plus(i,j),0.d0)**2 )
        else
          norm_grad_phi_sq = max(min(phi_x_minus(i,j),0.d0)**2,
     &                           max(phi_x_plus(i,j),0.d0)**2 )
     &                     + max(min(phi_y_minus(i,j),0.d0)**2,
     &                           max(phi_y_plus(i,j),0.d0)**2 )
        endif

c       } end Godunov selection of grad_phi


c       compute contribution to lse_rhs(i,j) 
        if (abs(vel_n) .ge. tol) then
          lse_rhs(i,j) = lse_rhs(i,j) 
     &                   - vel_n*sqrt(norm_grad_phi_sq)
        endif
	
	endif
      enddo 
c     } end loop over indexed points
	
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm2dAddConstCurvTermToLSERHSLOCAL() adds the contribution of a curvature 
c  term to the right-hand side of the level set equation when it is 
c  written in the form:
c
c    phi_t = -b*kappa*|grad(phi)| + ...
c  
c  kappa (mean curvature) will be computed from second order derivatives.
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    lse_rhs (in/out):  right-hand of level set equation
c    phi_*      (in):   derivatives (the 1st and 2nd order)   
c    b     (in):        scalar curvature term component 
c    *_gb (in):         index range for ghostbox
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):     index range of points in index_*
c    narrow_band(in):  array that marks voxels outside desired fillbox
c    mark_fb(in):      upper limit narrow band value for voxels in 
c                      fillbox
c
c***********************************************************************
      subroutine lsm2dAddConstCurvTermToLSERHSLOCAL(
     &  lse_rhs,
     &  ilo_lse_rhs_gb, ihi_lse_rhs_gb,
     &  jlo_lse_rhs_gb, jhi_lse_rhs_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  phi_xx, phi_xy, phi_yy, 
     &  ilo_grad2_phi_gb, ihi_grad2_phi_gb,
     &  jlo_grad2_phi_gb, jhi_grad2_phi_gb,
     &  b,
     &  index_x,
     &  index_y, 
     &  nlo_index, nhi_index,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_lse_rhs_gb, ihi_lse_rhs_gb
      integer jlo_lse_rhs_gb, jhi_lse_rhs_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_grad2_phi_gb, ihi_grad2_phi_gb
      integer jlo_grad2_phi_gb, jhi_grad2_phi_gb
      real lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb,
     &             jlo_lse_rhs_gb:jhi_lse_rhs_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_xx(ilo_grad2_phi_gb:ihi_grad2_phi_gb,
     &            jlo_grad2_phi_gb:jhi_grad2_phi_gb)
      real phi_yy(ilo_grad2_phi_gb:ihi_grad2_phi_gb,
     &            jlo_grad2_phi_gb:jhi_grad2_phi_gb)
      real phi_xy(ilo_grad2_phi_gb:ihi_grad2_phi_gb,
     &            jlo_grad2_phi_gb:jhi_grad2_phi_gb)
      real b
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb  
c     local variables      
      integer i,j,l
      real grad_mag2, curv
      real tol
      parameter (tol=1.d-13)

c     { begin loop over indexed points
      do l= nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)
	
       if( narrow_band(i,j) .le. mark_fb ) then	
     
c           compute squared magnitude of gradient
	    grad_mag2 = phi_x(i,j) * phi_x(i,j) + 
     &	                phi_y(i,j) * phi_y(i,j)
	    if(grad_mag2 .lt. tol) then
	      curv = 0.d0
	    else
	      curv = phi_xx(i,j)*phi_y(i,j)*phi_y(i,j)  
     &	         +   phi_yy(i,j)*phi_x(i,j)*phi_x(i,j)  
     &	         - 2*phi_xy(i,j)*phi_x(i,j)*phi_y(i,j)
	      curv = curv / grad_mag2 
	      endif

	      lse_rhs(i,j) = lse_rhs(i,j) + b*curv
	      
	endif      
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

