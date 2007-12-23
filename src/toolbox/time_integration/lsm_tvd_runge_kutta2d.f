c***********************************************************************
c
c  File:        lsm_tvd_runge_kutta2d.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision: 1.8 $
c  Modified:    $Date: 2007/02/11 17:14:03 $
c  Description: F77 routines for 2D TVD Runge-Kutta time integration 
c
c***********************************************************************

c***********************************************************************
c The TVD Runge-Kutta methods used in these subroutines are discussed 
c in Osher & Fedkiw (2003).
c***********************************************************************

c***********************************************************************
c
c  lsm2dRK1Step() takes a single first-order Runge-Kutta (i.e. Forward 
c  Euler) step.
c  
c  Arguments:
c    u_next (out):  u(t_cur+dt)
c    u_cur (in):    u(t_cur)
c    rhs (in):      right-hand side of time evolution equation
c    dt (in):       step size
c    *_gb (in):     index range for ghostbox
c    *_fb (in):     index range for fillbox
c
c***********************************************************************
      subroutine lsm2dRK1Step(
     &  u_next,
     &  ilo_u_next_gb, ihi_u_next_gb,
     &  jlo_u_next_gb, jhi_u_next_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_next_gb, ihi_u_next_gb
      integer jlo_u_next_gb, jhi_u_next_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real u_next(ilo_u_next_gb:ihi_u_next_gb,
     &                        jlo_u_next_gb:jhi_u_next_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb)
      integer i, j
      real dt

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          u_next(i,j) = u_cur(i,j) + dt*rhs(i,j)

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dTVDRK2Stage1() advances the solution through first stage of the 
c  second-order TVD Runge-Kutta step.
c  
c  Arguments:
c    u_stage1 (out):  u_approx(t_cur+dt)
c    u_cur (in):      u(t_cur)
c    rhs (in):        right-hand side of time evolution equation
c    dt (in):         step size
c    *_gb (in):       index range for ghostbox
c    *_fb (in):       index range for fillbox
c
c  NOTES:
c   - the first stage of TVD RK2 is identical to a single RK1 step
c
c***********************************************************************
      subroutine lsm2dTVDRK2Stage1(
     &  u_stage1,
     &  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_stage1_gb, ihi_u_stage1_gb
      integer jlo_u_stage1_gb, jhi_u_stage1_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb,
     &                          jlo_u_stage1_gb:jhi_u_stage1_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb)
      real dt

c     use lsm2dRK1Step() to compute first stage
      call lsm2dRK1Step(u_stage1, 
     &                  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &                  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &                  u_cur,
     &                  ilo_u_cur_gb, ihi_u_cur_gb,
     &                  jlo_u_cur_gb, jhi_u_cur_gb,
     &                  rhs,
     &                  ilo_rhs_gb, ihi_rhs_gb,
     &                  jlo_rhs_gb, jhi_rhs_gb,
     &                  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &                  dt)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dTVDRK2Stage2() completes advancing the solution through a 
c  single step of the second-order TVD Runge-Kutta step.
c  
c  Arguments:
c    u_next (out):   u(t_cur+dt)
c    u_stage1 (in):  u_approx(t_cur+dt)
c    u_cur (in):     u(t_cur)
c    rhs (in):       right-hand side of time evolution equation
c    dt (in):        step size
c    *_gb (in):      index range for ghostbox
c    *_fb (in):      index range for fillbox
c
c***********************************************************************
      subroutine lsm2dTVDRK2Stage2(
     &  u_next,
     &  ilo_u_next_gb, ihi_u_next_gb,
     &  jlo_u_next_gb, jhi_u_next_gb,
     &  u_stage1,
     &  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_next_gb, ihi_u_next_gb
      integer jlo_u_next_gb, jhi_u_next_gb
      integer ilo_u_stage1_gb, ihi_u_stage1_gb
      integer jlo_u_stage1_gb, jhi_u_stage1_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real u_next(ilo_u_next_gb:ihi_u_next_gb,
     &                        jlo_u_next_gb:jhi_u_next_gb)
      real u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb,
     &                          jlo_u_stage1_gb:jhi_u_stage1_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb)
      integer i, j
      real dt

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          u_next(i,j) = 0.5d0*( u_cur(i,j) 
     &                        + u_stage1(i,j) + dt*rhs(i,j) )

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dTVDRK3Stage1() advances the solution through first stage of the 
c  third-order TVD Runge-Kutta step.
c  
c  Arguments:
c    u_stage1 (out):  u_approx(t_cur+dt)
c    u_cur (in):      u(t_cur)
c    rhs (in):        right-hand side of time evolution equation
c    dt (in):         step size
c    *_gb (in):       index range for ghostbox
c    *_fb (in):       index range for fillbox
c
c  NOTES:
c   - the first stage of TVD RK3 is identical to a single RK1 step
c
c***********************************************************************
      subroutine lsm2dTVDRK3Stage1(
     &  u_stage1,
     &  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_stage1_gb, ihi_u_stage1_gb
      integer jlo_u_stage1_gb, jhi_u_stage1_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb,
     &                          jlo_u_stage1_gb:jhi_u_stage1_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb)
      real dt

c     use lsm2dRK1Step() to compute first stage
      call lsm2dRK1Step(u_stage1, 
     &                  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &                  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &                  u_cur,
     &                  ilo_u_cur_gb, ihi_u_cur_gb,
     &                  jlo_u_cur_gb, jhi_u_cur_gb,
     &                  rhs,
     &                  ilo_rhs_gb, ihi_rhs_gb,
     &                  jlo_rhs_gb, jhi_rhs_gb,
     &                  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &                  dt)

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dTVDRK3Stage2() advances the solution through second stage of 
c  the third-order TVD Runge-Kutta step.
c  
c  Arguments:
c    u_stage2 (out):  u_approx(t_cur+dt/2)
c    u_stage1 (in):   u_approx(t_cur+dt)
c    u_cur (in):      u(t_cur)
c    rhs (in):        right-hand side of time evolution equation
c    dt (in):         step size
c    *_gb (in):       index range for ghostbox
c    *_fb (in):       index range for fillbox
c
c***********************************************************************
      subroutine lsm2dTVDRK3Stage2(
     &  u_stage2,
     &  ilo_u_stage2_gb, ihi_u_stage2_gb,
     &  jlo_u_stage2_gb, jhi_u_stage2_gb,
     &  u_stage1,
     &  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_stage2_gb, ihi_u_stage2_gb
      integer jlo_u_stage2_gb, jhi_u_stage2_gb
      integer ilo_u_stage1_gb, ihi_u_stage1_gb
      integer jlo_u_stage1_gb, jhi_u_stage1_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real u_stage2(ilo_u_stage2_gb:ihi_u_stage2_gb,
     &                          jlo_u_stage2_gb:jhi_u_stage2_gb)
      real u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb,
     &                          jlo_u_stage1_gb:jhi_u_stage1_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb)
      integer i, j
      real dt

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          u_stage2(i,j) = 0.75d0*u_cur(i,j) 
     &                  + 0.25d0*( u_stage1(i,j) + dt*rhs(i,j) )

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dTVDRK3Stage3() completes advancing the solution through a 
c  single step of the third-order TVD Runge-Kutta step.
c  
c  Arguments:
c    u_next (out):   u(t_cur+dt)
c    u_stage2 (in):  u_approx(t_cur+dt/2)
c    u_cur (in):     u(t_cur)
c    rhs (in):       right-hand side of time evolution equation
c    dt (in):        step size
c    *_gb (in):      index range for ghostbox
c    *_fb (in):      index range for fillbox
c
c***********************************************************************
      subroutine lsm2dTVDRK3Stage3(
     &  u_next,
     &  ilo_u_next_gb, ihi_u_next_gb,
     &  jlo_u_next_gb, jhi_u_next_gb,
     &  u_stage2,
     &  ilo_u_stage2_gb, ihi_u_stage2_gb,
     &  jlo_u_stage2_gb, jhi_u_stage2_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_next_gb, ihi_u_next_gb
      integer jlo_u_next_gb, jhi_u_next_gb
      integer ilo_u_stage2_gb, ihi_u_stage2_gb
      integer jlo_u_stage2_gb, jhi_u_stage2_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real u_next(ilo_u_next_gb:ihi_u_next_gb,
     &                        jlo_u_next_gb:jhi_u_next_gb)
      real u_stage2(ilo_u_stage2_gb:ihi_u_stage2_gb,
     &                          jlo_u_stage2_gb:jhi_u_stage2_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb)
      integer i, j
      real dt
      real one_third, two_thirds
      parameter (one_third = 1.d0/3.d0)
      parameter (two_thirds = 2.d0/3.d0)

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          u_next(i,j) = one_third*u_cur(i,j)
     &                + two_thirds*( u_stage2(i,j) + dt*rhs(i,j) )

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
