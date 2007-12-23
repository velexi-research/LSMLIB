c***********************************************************************
c
c  File:        lsm_samrai_f77_utilities.f
c  Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
c  Revision:    $Revision: 1.3 $
c  Modified:    $Date: 2006/02/09 16:43:04 $
c  Description: Utility F77 subroutines for SAMRAI implementation of
c               level set method
c
c***********************************************************************

c***********************************************************************
c
c  lsm1dSAMRAIUtilitiesCopyData() copies data from the source to the 
c  destination for 1D problems.
c
c  Arguments:
c    dst_data (in):     dst data to be copied
c    src_data (out):  src space
c    *_gb (in):           index range for ghostbox
c    *_fb (in):           index range for fillbox
c
c***********************************************************************
      subroutine lsm1dSAMRAIUtilitiesCopyData(
     &  dst_data,
     &  ilo_dst_gb, ihi_dst_gb, 
     &  src_data,
     &  ilo_src_gb, ihi_src_gb, 
     &  ilo_fb, ihi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_dst_gb, ihi_dst_gb
      integer ilo_src_gb, ihi_src_gb
      integer ilo_fb, ihi_fb
      real dst_data(ilo_dst_gb:ihi_dst_gb)
      real src_data(ilo_src_gb:ihi_src_gb)
      integer i

c     loop over cells with sufficient data {
      do i=ilo_fb,ihi_fb

        dst_data(i) = src_data(i)

      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dSAMRAIUtilitiesCopyData() copies data from the source to the 
c  destination for 2D problems.
c
c  Arguments:
c    dst_data (in):     dst data to be copied
c    src_data (out):  src space
c    *_gb (in):           index range for ghostbox
c    *_fb (in):           index range for fillbox
c
c***********************************************************************
      subroutine lsm2dSAMRAIUtilitiesCopyData(
     &  dst_data,
     &  ilo_dst_gb, ihi_dst_gb, 
     &  jlo_dst_gb, jhi_dst_gb, 
     &  src_data,
     &  ilo_src_gb, ihi_src_gb, 
     &  jlo_src_gb, jhi_src_gb, 
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_dst_gb, ihi_dst_gb
      integer jlo_dst_gb, jhi_dst_gb
      integer ilo_src_gb, ihi_src_gb
      integer jlo_src_gb, jhi_src_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real dst_data(ilo_dst_gb:ihi_dst_gb,
     &                          jlo_dst_gb:jhi_dst_gb)
      real src_data(ilo_src_gb:ihi_src_gb,
     &                          jlo_src_gb:jhi_src_gb)
      integer i,j

c     loop over cells with sufficient data {
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          dst_data(i,j) = src_data(i,j)

        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dSAMRAIUtilitiesCopyData() copies data from the source to the 
c  destination for 3D problems.
c
c  Arguments:
c    dst_data (in):     dst data to be copied
c    src_data (out):  src space
c    *_gb (in):           index range for ghostbox
c    *_fb (in):           index range for fillbox
c
c***********************************************************************
      subroutine lsm3dSAMRAIUtilitiesCopyData(
     &  dst_data,
     &  ilo_dst_gb, ihi_dst_gb, 
     &  jlo_dst_gb, jhi_dst_gb, 
     &  klo_dst_gb, khi_dst_gb, 
     &  src_data,
     &  ilo_src_gb, ihi_src_gb, 
     &  jlo_src_gb, jhi_src_gb, 
     &  klo_src_gb, khi_src_gb, 
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  klo_fb, khi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_dst_gb, ihi_dst_gb
      integer jlo_dst_gb, jhi_dst_gb
      integer klo_dst_gb, khi_dst_gb
      integer ilo_src_gb, ihi_src_gb
      integer jlo_src_gb, jhi_src_gb
      integer klo_src_gb, khi_src_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      integer klo_fb, khi_fb
      real dst_data(ilo_dst_gb:ihi_dst_gb,
     &                          jlo_dst_gb:jhi_dst_gb,
     &                          klo_dst_gb:khi_dst_gb)
      real src_data(ilo_src_gb:ihi_src_gb,
     &                          jlo_src_gb:jhi_src_gb,
     &                          klo_src_gb:khi_src_gb)
      integer i,j,k

c     loop over cells with sufficient data {
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            dst_data(i,j,k) = src_data(i,j,k)

          enddo
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************
