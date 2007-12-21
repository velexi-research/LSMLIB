
c***********************************************************************
c
c  lsm3dRK1StepLOCAL() takes a single first-order Runge-Kutta (i.e. Forward 
c  Euler) step. The routine loops only over local (narrow band) points.
c  The routine loops only over local (narrow band) points.
c  
c  Arguments:
c    u_next (out):     u(t_cur+dt)
c    u_cur (in):       u(t_cur)
c    rhs (in):         right-hand side of time evolution equation
c    dt (in):          step size
c    *_gb (in):        index range for ghostbox
c    *_fb (in):        index range for fillbox
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):     index range of points to loop over in index_*
c    narrow_band(in):  array that marks voxels outside desired fillbox
c    mark_fb(in):      upper limit narrow band value for voxels in 
c                      fillbox 
c
c***********************************************************************
      subroutine lsm3dRK1StepLOCAL(
     &  u_next,
     &  ilo_u_next_gb, ihi_u_next_gb,
     &  jlo_u_next_gb, jhi_u_next_gb,
     &  klo_u_next_gb, khi_u_next_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  klo_u_cur_gb, khi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  dt,
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

      integer ilo_u_next_gb, ihi_u_next_gb
      integer jlo_u_next_gb, jhi_u_next_gb
      integer klo_u_next_gb, khi_u_next_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer klo_u_cur_gb, khi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      real u_next(ilo_u_next_gb:ihi_u_next_gb,
     &                        jlo_u_next_gb:jhi_u_next_gb,
     &                        klo_u_next_gb:khi_u_next_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb,
     &                       klo_u_cur_gb:khi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb,
     &                     klo_rhs_gb:khi_rhs_gb)
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
      integer i,j,k,l     
      real dt
     
c     { begin loop over indexed points
      do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)
	k=index_z(l)

        if( narrow_band(i,j,k) .le. mark_fb ) then	             
	   u_next(i,j,k) = u_cur(i,j,k) + dt*rhs(i,j,k)             
        endif
      enddo
c     } end loop over indexed points

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm3dTVDRK2Stage1LOCAL() advances the solution through first stage of the 
c  second-order TVD Runge-Kutta method. The routine loops only over local 
c  (narrow band) points.
c  The routine loops only over local (narrow band) points.
c  
c  Arguments:
c    u_stage1 (out):   u_approx(t_cur+dt)
c    u_cur (in):       u(t_cur)
c    rhs (in):         right-hand side of time evolution equation
c    dt (in):          step size
c    *_gb (in):        index range for ghostbox
c    *_fb (in):        index range for fillbox
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):     index range of points to loop over in index_*
c    narrow_band(in):  array that marks voxels outside desired fillbox
c    mark_fb(in):      upper limit narrow band value for voxels in 
c                      fillbox 
c
c  NOTES:
c   - the first stage of TVD RK2 is identical to a single RK1 step
c
c***********************************************************************
      subroutine lsm3dTVDRK2Stage1LOCAL(
     &  u_stage1,
     &  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &  klo_u_stage1_gb, khi_u_stage1_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  klo_u_cur_gb, khi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  dt,
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

      integer ilo_u_stage1_gb, ihi_u_stage1_gb
      integer jlo_u_stage1_gb, jhi_u_stage1_gb
      integer klo_u_stage1_gb, khi_u_stage1_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer klo_u_cur_gb, khi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      real u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb,
     &                          jlo_u_stage1_gb:jhi_u_stage1_gb,
     &                          klo_u_stage1_gb:khi_u_stage1_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb,
     &                       klo_u_cur_gb:khi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb,
     &                     klo_rhs_gb:khi_rhs_gb)
      real dt
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

c     use lsm3dRK1StepLOCAL() to compute first stage
      call lsm3dRK1StepLOCAL(u_stage1,
     &                  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &                  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &                  klo_u_stage1_gb, khi_u_stage1_gb,
     &                  u_cur,
     &                  ilo_u_cur_gb, ihi_u_cur_gb,
     &                  jlo_u_cur_gb, jhi_u_cur_gb,
     &                  klo_u_cur_gb, khi_u_cur_gb,
     &                  rhs,
     &                  ilo_rhs_gb, ihi_rhs_gb,
     &                  jlo_rhs_gb, jhi_rhs_gb,
     &                  klo_rhs_gb, khi_rhs_gb,
     &                  dt,
     &                  index_x,
     &                  index_y, 
     &                  index_z, 
     &                  nlo_index, nhi_index,
     &                  narrow_band,
     &                  ilo_nb_gb, ihi_nb_gb, 
     &                  jlo_nb_gb, jhi_nb_gb, 
     &                  klo_nb_gb, khi_nb_gb,
     &                  mark_fb)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dTVDRK2Stage2LOCAL() completes advancing the solution through a 
c  single step of the second-order TVD Runge-Kutta method. The routine 
c  loops only over local (narrow band) points.
c  The routine loops only over local (narrow band) points.
c  
c  Arguments:
c    u_next (out):    u(t_cur+dt)
c    u_stage1 (in):   u_approx(t_cur+dt)
c    u_cur (in):      u(t_cur)
c    rhs (in):        right-hand side of time evolution equation
c    dt (in):         step size
c    *_gb (in):       index range for ghostbox
c    *_fb (in):       index range for fillbox
c    index_[xyz](in): [xyz] coordinates of local (narrow band) points
c    n*_index(in):    index range of points to loop over in index_* 
c
c***********************************************************************
      subroutine lsm3dTVDRK2Stage2LOCAL(
     &  u_next,
     &  ilo_u_next_gb, ihi_u_next_gb,
     &  jlo_u_next_gb, jhi_u_next_gb,
     &  klo_u_next_gb, khi_u_next_gb,
     &  u_stage1,
     &  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &  klo_u_stage1_gb, khi_u_stage1_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  klo_u_cur_gb, khi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  dt,
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

      integer ilo_u_next_gb, ihi_u_next_gb
      integer jlo_u_next_gb, jhi_u_next_gb
      integer klo_u_next_gb, khi_u_next_gb
      integer ilo_u_stage1_gb, ihi_u_stage1_gb
      integer jlo_u_stage1_gb, jhi_u_stage1_gb
      integer klo_u_stage1_gb, khi_u_stage1_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer klo_u_cur_gb, khi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      real u_next(ilo_u_next_gb:ihi_u_next_gb,
     &                        jlo_u_next_gb:jhi_u_next_gb,
     &                        klo_u_next_gb:khi_u_next_gb)
      real u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb,
     &                          jlo_u_stage1_gb:jhi_u_stage1_gb,
     &                          klo_u_stage1_gb:khi_u_stage1_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb,
     &                       klo_u_cur_gb:khi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb,
     &                     klo_rhs_gb:khi_rhs_gb)
      real dt 
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
      integer i, j, k, l

c     { begin loop over indexed points
      do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)
	k=index_z(l)
	
        if( narrow_band(i,j,k) .le. mark_fb ) then	      
          u_next(i,j,k) = 0.5d0*( u_cur(i,j,k) 
     &                        + u_stage1(i,j,k) + dt*rhs(i,j,k) )
        endif
	
      enddo
c     } end loop over indexed points

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm3dTVDRK3Stage1LOCAL() advances the solution through first stage of the 
c  third-order TVD Runge-Kutta step.
c  The routine loops only over local (narrow band) points.
c  
c  Arguments:
c    u_stage1 (out):   u_approx(t_cur+dt)
c    u_cur (in):       u(t_cur)
c    rhs (in):         right-hand side of time evolution equation
c    dt (in):          step size
c    *_gb (in):        index range for ghostbox
c    *_fb (in):        index range for fillbox
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):     index range of points to loop over in index_*
c    narrow_band(in):  array that marks voxels outside desired fillbox
c    mark_fb(in):      upper limit narrow band value for voxels in 
c                      fillbox 
c
c  NOTES:
c   - the first stage of TVD RK3 is identical to a single RK1 step
c
c***********************************************************************
      subroutine lsm3dTVDRK3Stage1LOCAL(
     &  u_stage1,
     &  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &  klo_u_stage1_gb, khi_u_stage1_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  klo_u_cur_gb, khi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  dt,
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

      integer ilo_u_stage1_gb, ihi_u_stage1_gb
      integer jlo_u_stage1_gb, jhi_u_stage1_gb
      integer klo_u_stage1_gb, khi_u_stage1_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer klo_u_cur_gb, khi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      real u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb,
     &                          jlo_u_stage1_gb:jhi_u_stage1_gb,
     &                          klo_u_stage1_gb:khi_u_stage1_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb,
     &                       klo_u_cur_gb:khi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb,
     &                     klo_rhs_gb:khi_rhs_gb)
      real dt
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

c     use lsm3dRK1Step() to compute first stage
      call lsm3dRK1StepLOCAL(u_stage1,
     &                  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &                  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &                  klo_u_stage1_gb, khi_u_stage1_gb,
     &                  u_cur,
     &                  ilo_u_cur_gb, ihi_u_cur_gb,
     &                  jlo_u_cur_gb, jhi_u_cur_gb,
     &                  klo_u_cur_gb, khi_u_cur_gb,
     &                  rhs,
     &                  ilo_rhs_gb, ihi_rhs_gb,
     &                  jlo_rhs_gb, jhi_rhs_gb,
     &                  klo_rhs_gb, khi_rhs_gb,
     &                  dt,
     &                  index_x,
     &                  index_y, 
     &                  index_z, 
     &                  nlo_index, nhi_index,
     &                  narrow_band,
     &                  ilo_nb_gb, ihi_nb_gb, 
     &                  jlo_nb_gb, jhi_nb_gb, 
     &                  klo_nb_gb, khi_nb_gb,
     &                  mark_fb)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dTVDRK3Stage2LOCAL() advances the solution through second stage of 
c  the third-order TVD Runge-Kutta step.
c  The routine loops only over local (narrow band) points.
c  
c  Arguments:
c    u_stage2 (out):  u_approx(t_cur+dt/2)
c    u_stage1 (in):   u_approx(t_cur+dt)
c    u_cur (in):      u(t_cur)
c    rhs (in):        right-hand side of time evolution equation
c    dt (in):         step size
c    *_gb (in):       index range for ghostbox
c    *_fb (in):       index range for fillbox
c    index_[xyz](in): [xyz] coordinates of local (narrow band) points
c    n*_index(in):    index range of points to loop over in index_*
c    narrow_band(in): array that marks voxels outside desired fillbox
c    mark_fb(in):     upper limit narrow band value for voxels in 
c                     fillbox
c
c***********************************************************************
      subroutine lsm3dTVDRK3Stage2LOCAL(
     &  u_stage2,
     &  ilo_u_stage2_gb, ihi_u_stage2_gb,
     &  jlo_u_stage2_gb, jhi_u_stage2_gb,
     &  klo_u_stage2_gb, khi_u_stage2_gb,
     &  u_stage1,
     &  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &  klo_u_stage1_gb, khi_u_stage1_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  klo_u_cur_gb, khi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  dt,
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

      integer ilo_u_stage2_gb, ihi_u_stage2_gb
      integer jlo_u_stage2_gb, jhi_u_stage2_gb
      integer klo_u_stage2_gb, khi_u_stage2_gb
      integer ilo_u_stage1_gb, ihi_u_stage1_gb
      integer jlo_u_stage1_gb, jhi_u_stage1_gb
      integer klo_u_stage1_gb, khi_u_stage1_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer klo_u_cur_gb, khi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      real u_stage2(ilo_u_stage2_gb:ihi_u_stage2_gb,
     &                          jlo_u_stage2_gb:jhi_u_stage2_gb,
     &                          klo_u_stage2_gb:khi_u_stage2_gb)
      real u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb,
     &                          jlo_u_stage1_gb:jhi_u_stage1_gb,
     &                          klo_u_stage1_gb:khi_u_stage1_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb,
     &                       klo_u_cur_gb:khi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb,
     &                     klo_rhs_gb:khi_rhs_gb)
      real dt
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
      integer i,j,k,l
     
c     { begin loop over indexed points
      do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)
	k=index_z(l)

        if( narrow_band(i,j,k) .le. mark_fb ) then	  

            u_stage2(i,j,k) = 0.75d0*u_cur(i,j,k) 
     &                      + 0.25d0*(u_stage1(i,j,k) + dt*rhs(i,j,k))

        endif
      enddo
c     } end loop over indexed points

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dTVDRK3Stage3LOCAL() completes advancing the solution through a 
c  single step of the third-order TVD Runge-Kutta step.
c  The routine loops only over local (narrow band) points.
c  
c  Arguments:
c    u_next (out):     u(t_cur+dt)
c    u_stage2 (in):    u_approx(t_cur+dt/2)
c    u_cur (in):       u(t_cur)
c    rhs (in):         right-hand side of time evolution equation
c    dt (in):          step size
c    *_gb (in):        index range for ghostbox
c    *_fb (in):        index range for fillbox
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):     index range of points to loop over in index_*
c    narrow_band(in):  array that marks voxels outside desired fillbox
c    mark_fb(in):      upper limit narrow band value for voxels in 
c                      fillbox
c
c***********************************************************************
      subroutine lsm3dTVDRK3Stage3LOCAL(
     &  u_next,
     &  ilo_u_next_gb, ihi_u_next_gb,
     &  jlo_u_next_gb, jhi_u_next_gb,
     &  klo_u_next_gb, khi_u_next_gb,
     &  u_stage2,
     &  ilo_u_stage2_gb, ihi_u_stage2_gb,
     &  jlo_u_stage2_gb, jhi_u_stage2_gb,
     &  klo_u_stage2_gb, khi_u_stage2_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  klo_u_cur_gb, khi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  dt,
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

      integer ilo_u_next_gb, ihi_u_next_gb
      integer jlo_u_next_gb, jhi_u_next_gb
      integer klo_u_next_gb, khi_u_next_gb
      integer ilo_u_stage2_gb, ihi_u_stage2_gb
      integer jlo_u_stage2_gb, jhi_u_stage2_gb
      integer klo_u_stage2_gb, khi_u_stage2_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer klo_u_cur_gb, khi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      real u_next(ilo_u_next_gb:ihi_u_next_gb,
     &                        jlo_u_next_gb:jhi_u_next_gb,
     &                        klo_u_next_gb:khi_u_next_gb)
      real u_stage2(ilo_u_stage2_gb:ihi_u_stage2_gb,
     &                          jlo_u_stage2_gb:jhi_u_stage2_gb,
     &                          klo_u_stage2_gb:khi_u_stage2_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb,
     &                       klo_u_cur_gb:khi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb,
     &                     klo_rhs_gb:khi_rhs_gb)
      real dt
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
      integer i,j,k,l
      
      real one_third, two_thirds
      parameter (one_third = 1.d0/3.d0)
      parameter (two_thirds = 2.d0/3.d0)

c     { begin loop over indexed points
      do l=nlo_index, nhi_index      
        i=index_x(l)
	j=index_y(l)
	k=index_z(l)

       if( narrow_band(i,j,k) .le. mark_fb ) then	 
            u_next(i,j,k) = one_third*u_cur(i,j,k)
     &                    + two_thirds*( u_stage2(i,j,k) 
     &                                 + dt*rhs(i,j,k) )

       
        endif
      enddo
c     } end loop over indexed points

      return
      end
c } end subroutine
c***********************************************************************
