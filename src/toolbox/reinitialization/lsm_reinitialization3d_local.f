c***********************************************************************
c
c  lsm3dComputeReinitializationEqnRHSLOCAL() computes the right-hand side of 
c  the reinitialization equation using a Godunov scheme to select the 
c  numerical discretization of the sgn(phi) |grad(phi)| term.  
c  Forward (plus) and backward (minus) spatial derivatives used in 
c  the Godunov calculation must be supplied by the user.
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    reinit_rhs (out):       right-hand side of reinitialization 
c                            equation
c    phi (in):               level set function at current iteration
c                            of reinitialization process
c    phi0 (in):              level set function at initial iteration
c                            iteration of reinitialization process
c    phi_*_plus (in):        forward spatial derivatives for grad(phi)
c    phi_*_minus (in):       backward spatial derivatives for grad(phi)
c    
c    use_phi0_for_sgn (in):  flag to specify whether phi0 should be
c                            used in the computation of sgn(phi).
c                              0 = use phi (do NOT use phi0)
c                              1 = use phi0
c    *_gb (in):              index range for ghostbox
c    index_[xyz](in):       [xyz] coordinates of local (narrow band) points
c    n*_index(in):          index range of points to loop over in index_*
c    narrow_band(in):    array that marks voxels outside desired fillbox
c    mark_fb(in):        upper limit narrow band value for voxels in 
c                        fillbox
c
c  NOTES:
c   (1) if use_phi0_for_sgn is not equal to 0 or 1, the default
c       behavior of lsm3dComputeReinitializationEqnRHS() is to use 
c       phi (i.e. equivalent to setting use_phi0_for_sgn to 0)
c
c***********************************************************************
      subroutine lsm3dComputeReinitializationEqnRHSLOCAL(
     &  reinit_rhs,
     &  ilo_rhs_gb, ihi_rhs_gb, 
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  phi0,
     &  ilo_phi0_gb, ihi_phi0_gb, 
     &  jlo_phi0_gb, jhi_phi0_gb,
     &  klo_phi0_gb, khi_phi0_gb,
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb,
     &  dx, dy, dz,
     &  use_phi0_for_sgn,
     &  index_x,
     &  index_y, 
     &  index_z, 
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_phi0_gb, ihi_phi0_gb
      integer jlo_phi0_gb, jhi_phi0_gb
      integer klo_phi0_gb, khi_phi0_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      double precision reinit_rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                            jlo_rhs_gb:jhi_rhs_gb,
     &                            klo_rhs_gb:khi_rhs_gb)
      double precision phi(ilo_phi_gb:ihi_phi_gb,
     &                     jlo_phi_gb:jhi_phi_gb,
     &                     klo_phi_gb:khi_phi_gb)
      double precision phi0(ilo_phi0_gb:ihi_phi0_gb,
     &                      jlo_phi0_gb:jhi_phi0_gb,
     &                      klo_phi0_gb:khi_phi0_gb)
      double precision phi_x_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                    klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      double precision phi_y_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                    klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      double precision phi_z_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                    klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      double precision phi_x_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                    klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      double precision phi_y_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                    klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      double precision phi_z_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                    klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      double precision dx, dy, dz
      integer use_phi0_for_sgn
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb     

c     local variables      
      double precision phi_cur
      integer DIM
      parameter (DIM=3)
      double precision grad_phi_plus_cur(1:DIM)
      double precision grad_phi_minus_cur(1:DIM)
      double precision grad_phi_star(1:DIM)
      integer i,j,k,l
      integer dir
      double precision sgn_phi
      double precision norm_grad_phi_sq
      double precision dx_sq
      double precision tol
      parameter (tol=1.d-13)
      double precision one
      parameter (one=1.d0)

c     set value of dx_sq to be square of max{dx,dy,dz}
      dx_sq = max(dx,dy,dz)
      dx_sq = dx_sq*dx_sq

c----------------------------------------------------
c      compute RHS of reinitialization equation
c      using Godunov's method
c----------------------------------------------------
c     { begin condition on use_phi0_for_sgn
      if (use_phi0_for_sgn .ne. 1) then

c       -----------------------------------------------
c       use phi in computation of smoothed sgn function
c       -----------------------------------------------
c     { begin loop over indexed points
      do l=nlo_index,nhi_index      
        i=index_x(l)
	j=index_y(l)
	k=index_z(l)
	 
        if( narrow_band(i,j,k) .le. mark_fb ) then	
      
c               cache phi and spatial derivative approximations
        	phi_cur = phi(i,j,k)
        	grad_phi_plus_cur(1) = phi_x_plus(i,j,k)
        	grad_phi_plus_cur(2) = phi_y_plus(i,j,k)
        	grad_phi_plus_cur(3) = phi_z_plus(i,j,k)
        	grad_phi_minus_cur(1) = phi_x_minus(i,j,k)
        	grad_phi_minus_cur(2) = phi_y_minus(i,j,k)
        	grad_phi_minus_cur(3) = phi_z_minus(i,j,k)

c               { begin Godunov selection of grad_phi
        	do dir=1,DIM

                  if (phi_cur .gt. 0.d0) then
                    grad_phi_plus_cur(dir) = 
     &                              max(-grad_phi_plus_cur(dir),0.d0)
                    grad_phi_minus_cur(dir) = 
     &                              max(grad_phi_minus_cur(dir),0.d0)
                  else
                    grad_phi_plus_cur(dir) = 
     &                max(grad_phi_plus_cur(dir), 0.d0)
                    grad_phi_minus_cur(dir) = 
     &                max(-grad_phi_minus_cur(dir), 0.d0)
                  endif

                  grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                     grad_phi_minus_cur(dir)) 

        	enddo
c               } end Godunov selection of grad_phi

c               compute reinit_rhs(i,j) using smoothed sgn(phi)
        	if (abs(phi_cur) .ge. tol) then
                  norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                             + grad_phi_star(2)*grad_phi_star(2)
     &                             + grad_phi_star(3)*grad_phi_star(3)
                  sgn_phi = phi_cur
     &                / sqrt(phi_cur*phi_cur + norm_grad_phi_sq*dx_sq)
                  reinit_rhs(i,j,k) = 
     &              sgn_phi*(one - sqrt(norm_grad_phi_sq))
        	else
                  reinit_rhs(i,j,k) = 0.d0
        	endif
          endif		
      enddo
c       } end loop over indexed points

      else

c       ------------------------------------------------
c       use phi0 in computation of smoothed sgn function 
c       ------------------------------------------------
c     { begin loop over indexed points
      do l=nlo_index,nhi_index      
        i=index_x(l)
	j=index_y(l)
	k=index_z(l)

       if( narrow_band(i,j,k) .le. mark_fb ) then	
     
c               cache phi and spatial derivative approximations
        	phi_cur = phi0(i,j,k)
        	grad_phi_plus_cur(1) = phi_x_plus(i,j,k)
        	grad_phi_plus_cur(2) = phi_y_plus(i,j,k)
        	grad_phi_plus_cur(3) = phi_z_plus(i,j,k)
        	grad_phi_minus_cur(1) = phi_x_minus(i,j,k)
        	grad_phi_minus_cur(2) = phi_y_minus(i,j,k)
        	grad_phi_minus_cur(3) = phi_z_minus(i,j,k)

c               { begin Godunov selection of grad_phi
        	do dir=1,DIM

                  if (phi_cur .gt. 0.d0) then
                    grad_phi_plus_cur(dir) = 
     &                             max(-grad_phi_plus_cur(dir),0.d0)
                    grad_phi_minus_cur(dir) = 
     &                             max(grad_phi_minus_cur(dir),0.d0)
                  else
                    grad_phi_plus_cur(dir) = 
     &                max(grad_phi_plus_cur(dir), 0.d0)
                    grad_phi_minus_cur(dir) = 
     &                max(-grad_phi_minus_cur(dir),0.d0)
                  endif

                  grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                     grad_phi_minus_cur(dir)) 

        	enddo
c               } end Godunov selection of grad_phi

c               compute reinit_rhs(i,j) using smoothed sgn(phi)
        	if (abs(phi_cur) .ge. tol) then
                  norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                             + grad_phi_star(2)*grad_phi_star(2)
     &                             + grad_phi_star(3)*grad_phi_star(3)
                  sgn_phi = phi_cur / sqrt(phi_cur*phi_cur + dx_sq)
                  reinit_rhs(i,j,k) = 
     &              sgn_phi*(one - sqrt(norm_grad_phi_sq))
        	else
                  reinit_rhs(i,j,k) = 0.d0
        	endif	    
          endif
       enddo
c      } end loop over indexed points

      endif
c     } end condition on use_phi0_for_sgn

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dComputeOrthogonalizationEqnRHSLOCAL() computes the right-hand side 
c  of the orthogonalization equation:
c
c    phi_t + grad(phi) dot { sgn(psi)/|grad(psi)| grad(psi) } = 0
c
c  Upwinding is used to select whether the forward (plus) or backward
c  (minus) spatial derivative should be used for grad(phi).  grad(psi)
c  is computed by averaging the forward and backward spatial derivatives
c  for grad(psi).  Forward and backward spatial derivatives used in the 
c  calculation must be supplied by the user.
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    othro_rhs (out):        right-hand side of orthogonalization
c                            equation
c    psi (in):               data array for psi
c    phi_*_plus (in):        forward spatial derivatives for grad(phi)
c    phi_*_minus (in):       backward spatial derivatives for grad(phi)
c    psi_*_plus (in):        forward spatial derivatives for grad(psi)
c    psi_*_minus (in):       backward spatial derivatives for grad(psi)
c    *_gb (in):              index range for ghostbox
c    *_fb (in):              index range for fillbox
c    index_[xyz](in):       [xyz] coordinates of local (narrow band) points
c    n*_index(in):          index range of points to loop over in index_*
c    narrow_band(in):       array that marks voxels outside desired fillbox
c    mark_fb(in):           upper limit narrow band value for voxels in 
c                           fillbox
c
c***********************************************************************
      subroutine lsm3dComputeOrthogonalizationEqnRHSLOCAL(
     &  ortho_rhs,
     &  ilo_rhs_gb, ihi_rhs_gb, 
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb,
     &  psi,
     &  ilo_psi_gb, ihi_psi_gb, 
     &  jlo_psi_gb, jhi_psi_gb,
     &  klo_psi_gb, khi_psi_gb,
     &  psi_x_plus, psi_y_plus, psi_z_plus,
     &  ilo_grad_psi_plus_gb, ihi_grad_psi_plus_gb, 
     &  jlo_grad_psi_plus_gb, jhi_grad_psi_plus_gb,
     &  klo_grad_psi_plus_gb, khi_grad_psi_plus_gb,
     &  psi_x_minus, psi_y_minus, psi_z_minus,
     &  ilo_grad_psi_minus_gb, ihi_grad_psi_minus_gb, 
     &  jlo_grad_psi_minus_gb, jhi_grad_psi_minus_gb,
     &  klo_grad_psi_minus_gb, khi_grad_psi_minus_gb,
     &  dx, dy, dz,
     &  index_x,
     &  index_y, 
     &  index_z, 
     &  nlo_index, nhi_index,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      integer ilo_psi_gb, ihi_psi_gb
      integer jlo_psi_gb, jhi_psi_gb
      integer klo_psi_gb, khi_psi_gb
      integer ilo_grad_psi_plus_gb, ihi_grad_psi_plus_gb
      integer jlo_grad_psi_plus_gb, jhi_grad_psi_plus_gb
      integer klo_grad_psi_plus_gb, khi_grad_psi_plus_gb
      integer ilo_grad_psi_minus_gb, ihi_grad_psi_minus_gb
      integer jlo_grad_psi_minus_gb, jhi_grad_psi_minus_gb
      integer klo_grad_psi_minus_gb, khi_grad_psi_minus_gb
      double precision ortho_rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                           jlo_rhs_gb:jhi_rhs_gb,
     &                           klo_rhs_gb:khi_rhs_gb)
      double precision psi(ilo_psi_gb:ihi_psi_gb,
     &                     jlo_psi_gb:jhi_psi_gb,
     &                     klo_psi_gb:khi_psi_gb)
      double precision phi_x_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                    klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      double precision phi_y_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                    klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      double precision phi_z_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                    klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      double precision phi_x_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                    klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      double precision phi_y_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                    klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      double precision phi_z_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                    klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      double precision psi_x_plus(
     &                    ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                    jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb,
     &                    klo_grad_psi_plus_gb:khi_grad_psi_plus_gb)
      double precision psi_y_plus(
     &                    ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                    jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb,
     &                    klo_grad_psi_plus_gb:khi_grad_psi_plus_gb)
      double precision psi_z_plus(
     &                    ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                    jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb,
     &                    klo_grad_psi_plus_gb:khi_grad_psi_plus_gb)
      double precision psi_x_minus(
     &                    ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                    jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb,
     &                    klo_grad_psi_minus_gb:khi_grad_psi_minus_gb)
      double precision psi_y_minus(
     &                    ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                    jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb,
     &                    klo_grad_psi_minus_gb:khi_grad_psi_minus_gb)
      double precision psi_z_minus(
     &                    ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                    jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb,
     &                    klo_grad_psi_minus_gb:khi_grad_psi_minus_gb)
      double precision dx, dy, dz
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb
      
      double precision dx_sq
      integer DIM
      parameter (DIM=3)
      double precision grad_psi_star(1:DIM)
      double precision norm_grad_psi
      double precision sgn_psi
      integer i,j,k,l
      double precision psi_tol, grad_psi_tol
      parameter (psi_tol=1.d-13,grad_psi_tol=1.d-8)
      double precision one, half
      parameter (one=1.d0,half=0.5d0)

c     set value of dx_sq to be square of max{dx,dy,dz}
      dx_sq = max(dx,dy,dz)
      dx_sq = dx_sq*dx_sq

c----------------------------------------------------
c      compute RHS of orthogonalization equation
c      using upwinding
c----------------------------------------------------

c     { begin loop over indexed points
      do l=nlo_index,nhi_index      
        i=index_x(l)
	j=index_y(l)
	k=index_z(l)

       if( narrow_band(i,j,k) .le. mark_fb ) then	

c           average forward and backward derivatives of 
c           grad(psi) to get grad_psi_star
            grad_psi_star(1) = half 
     &                       * (psi_x_plus(i,j,k)+ psi_x_minus(i,j,k))
            grad_psi_star(2) = half 
     &                       * (psi_y_plus(i,j,k)+ psi_y_minus(i,j,k))
            grad_psi_star(3) = half 
     &                       * (psi_z_plus(i,j,k)+ psi_z_minus(i,j,k))
 
c           compute norm of grad(psi)
            norm_grad_psi = grad_psi_star(1)*grad_psi_star(1)
     &                    + grad_psi_star(2)*grad_psi_star(2)
     &                    + grad_psi_star(3)*grad_psi_star(3)
            norm_grad_psi = sqrt(norm_grad_psi)

c           compute ortho_rhs(i,j,k) using upwinding on sgn_psi*grad(psi) 
            if ( (abs(psi(i,j,k)) .gt. psi_tol) .and.
     &           (norm_grad_psi .gt. grad_psi_tol) ) then

c             CASE: nontrivial psi and grad(psi) 

              sgn_psi = psi(i,j,k)/sqrt( psi(i,j,k)*psi(i,j,k) 
     &                                 + dx_sq*norm_grad_psi**2)

              ortho_rhs(i,j,k) = -1.d0/norm_grad_psi 
     &          * ( max(sgn_psi*grad_psi_star(1),0.d0)
     &                         *phi_x_minus(i,j,k)
     &            + min(sgn_psi*grad_psi_star(1),0.d0)
     &                          *phi_x_plus(i,j,k)
     &            + max(sgn_psi*grad_psi_star(2),0.d0)
     &                          *phi_y_minus(i,j,k)
     &            + min(sgn_psi*grad_psi_star(2),0.d0)
     &                          *phi_y_plus(i,j,k)
     &            + max(sgn_psi*grad_psi_star(3),0.d0)
     &                          *phi_z_minus(i,j,k)
     &            + min(sgn_psi*grad_psi_star(3),0.d0)
     &                          *phi_z_plus(i,j,k) )
            else

c             CASE: grad(psi) = 0 or psi = 0 

              ortho_rhs(i,j,k) = 0.d0

            endif

        endif
      enddo
c     } end loop over indexed points


      return
      end
c } end subroutine
c***********************************************************************
