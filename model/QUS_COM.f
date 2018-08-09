#include "rundeck_opts.h"
      MODULE SOMTQ_COM
!@sum  SOMTQ_COM contains the arrays containing second order moments
!@auth Gary Russell
      USE QUSDEF
      USE RESOLUTION, only : im,jm,lm
      IMPLICIT NONE
      SAVE
!     REAL*8, DIMENSION(NMOM,IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO ,LM)  
!    &        :: TMOM,QMOM
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TMOM,QMOM

      END MODULE SOMTQ_COM

      SUBROUTINE ALLOC_SMOMTQ(grid)
!@sum  init_smomtq allocates the arrays in this module which
!@+    must now be dynamic for the distributed memory implementation.
!@auth Rosalinda de Fainchtein
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE QUSDEF, ONLY : NMOM
      USE RESOLUTION, ONLY : LM
      USE SOMTQ_COM, ONLY : TMOM,QMOM
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE ( TMOM(NMOM , I_0H:I_1H , J_0H:J_1H , LM),
     &           QMOM(NMOM , I_0H:I_1H , J_0H:J_1H , LM),
     &   STAT=IER )

      TMOM = 0.
      QMOM = 0.

      END SUBROUTINE ALLOC_SMOMTQ

      SUBROUTINE io_somtq(kunit,iaction,ioerr)
!@sum  io_somtq reads and writes second order moments to file
!@auth Gavin Schmidt
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, only : PACK_COLUMN, UNPACK_COLUMN, AM_I_ROOT
      USE SOMTQ_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "QUS01"
      REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: TMOM_GLOB, QMOM_GLOB
      integer :: img, jmg

      IF(AM_I_ROOT()) then
         img = IM
         jmg = JM
      else
         img = 1
         jmg = 1
      end if
      ALLOCATE(TMOM_GLOB(NMOM,img,jmg,LM),QMOM_GLOB(NMOM,img,jmg,LM))
      write (MODULE_HEADER(lhead+1:80),'(a7,i2,a)')
     * 'R8 dim(',nmom,',im,jm,lm):Tmom,Qmom'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to standard restart file
        CALL PACK_COLUMN(grid, TMOM, TMOM_GLOB)
        CALL PACK_COLUMN(grid, QMOM, QMOM_GLOB)
        IF (AM_I_ROOT())
     &      WRITE (KUNIT,ERR=10) MODULE_HEADER, TMOM_GLOB, QMOM_GLOB

      CASE (IOREAD:)            ! input from restart file
        if ( AM_I_ROOT() ) then
          READ (KUNIT,ERR=10) HEADER, TMOM_GLOB, QMOM_GLOB
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if
        CALL UNPACK_COLUMN(grid, TMOM_GLOB, TMOM)
        CALL UNPACK_COLUMN(grid, QMOM_GLOB, QMOM)
      END SELECT

      call freespace
      RETURN
 10   IOERR=1
      call freespace
      RETURN

      contains

      subroutine freespace
      DEALLOCATE(TMOM_GLOB,QMOM_GLOB)
      end subroutine freespace

      END SUBROUTINE io_somtq

#ifdef NEW_IO
      subroutine def_rsf_somtq(fid)
!@sum  def_rsf_somtq defines QUS T/Q array structure in restart files
!@auth M. Kelley
!@ver  beta
      use somtq_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,tmom,'tmom(nmom,dist_im,dist_jm,lm)')
      call defvar(grid,fid,qmom,'qmom(nmom,dist_im,dist_jm,lm)')
      return
      end subroutine def_rsf_somtq

      subroutine new_io_somtq(fid,iaction)
!@sum  new_io_somtq read/write QUS T/Q arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use somtq_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to restart file
        call write_dist_data(grid, fid, 'tmom', tmom, jdim=3)
        call write_dist_data(grid, fid, 'qmom', qmom, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'tmom', tmom, jdim=3)
        call read_dist_data(grid, fid, 'qmom', qmom, jdim=3)
      end select
      return
      end subroutine new_io_somtq
#endif /* NEW_IO */

      subroutine tq_zmom_init(t,q,pmid,pedn)
      USE DOMAIN_DECOMP_ATM, ONLY: grid
      USE SOMTQ_COM
      implicit none
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) :: t,q
      REAL*8 :: pmid(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo)
      REAL*8 :: pedn(lm+1,grid%i_strt_halo:grid%i_stop_halo,
     &                    grid%j_strt_halo:grid%j_stop_halo)
      integer :: i,j,l
      REAL*8 :: rdsig

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S
C****
C**** Extract useful local domain parameters from "grid"
C****
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP

C**** INITIALIZES VERTICAL SLOPES OF T,Q
      DO J=J_0,J_1
        DO I=I_0,I_1
          RDSIG=(PMID(1,I,J)-PEDN(2,I,J))/(PMID(1,I,J)-PMID(2,I,J))
          TMOM(MZ,I,J,1)=(T(I,J,2)-T(I,J,1))*RDSIG
          QMOM(MZ,I,J,1)=(Q(I,J,2)-Q(I,J,1))*RDSIG
          IF(Q(I,J,1)+QMOM(MZ,I,J,1).LT.0.) QMOM(MZ,I,J,1)=-Q(I,J,1)
          DO L=2,LM-1
            RDSIG=(PMID(L,I,J)-PEDN(L+1,I,J))/
     *           (PMID(L-1,I,J)-PMID(L+1,I,J))
            TMOM(MZ,I,J,L)=(T(I,J,L+1)-T(I,J,L-1))*RDSIG
            QMOM(MZ,I,J,L)=(Q(I,J,L+1)-Q(I,J,L-1))*RDSIG
            IF(Q(I,J,L)+QMOM(MZ,I,J,L).LT.0.) QMOM(MZ,I,J,L)=-Q(I,J,L)
          END DO
          RDSIG=(PMID(LM,I,J)-PEDN(LM+1,I,J))/
     *         (PMID(LM-1,I,J)-PMID(LM,I,J))
          TMOM(MZ,I,J,LM)=(T(I,J,LM)-T(I,J,LM-1))*RDSIG
          QMOM(MZ,I,J,LM)=(Q(I,J,LM)-Q(I,J,LM-1))*RDSIG
          IF(Q(I,J,LM)+QMOM(MZ,I,J,LM).LT.0.) QMOM(MZ,I,J,LM)=-Q(I,J,LM)
        END DO
      END DO
      return
      end subroutine tq_zmom_init
