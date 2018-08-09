#include "rundeck_opts.h"

      MODULE PBLCOM
!@sum  PBLCOM contains the arrays used by the Boundary Layer code
!@auth Greg Hartke/Ye Cheng
      USE SOCPBL, only : npbl=>n
      IMPLICIT NONE
      SAVE

!@var ROUGHL log10(zgs/roughness length), prescribed with zgs=30 m.
      REAL*8, allocatable, dimension(:,:) :: roughl

!@var DCLEV     LAYER TO WHICH DRY CONVECTION MIXES (1)
!@var ugeo,vgeo components of geostrophic wind at the top of the BL
!@var bldep     boundary layer depth (m)
      REAL*8, allocatable, dimension(:,:) ::
     &     dclev,pblht,pblptop,ugeo,vgeo,bldep

!@var [tuv]1_after_aturb first-layer temp/winds after ATURB completes
!@+   (used to compute tendencies seen by the PBL code)
      REAL*8, allocatable, dimension(:,:) ::
     &     t1_after_aturb,u1_after_aturb,v1_after_aturb

!@var egcm  3-d turbulent kinetic energy in the whole atmosphere
!@var w2gcm vertical component of egcm
!@var t2gcm 3-d turbulent temperature variance in the whole atmosphere
      real*8, allocatable, dimension(:,:,:) :: egcm,w2gcm,t2gcm

!@ egcm_init_max maximum initial vaule of egcm
      real*8, parameter :: egcm_init_max=0.5d0

      END MODULE PBLCOM

c      SUBROUTINE io_pbl(kunit,iaction,ioerr)
c!@sum  io_pbl reads and writes model variables to file
c!@auth Gavin Schmidt
c      USE MODEL_COM, only : ioread,irsfic,irerun,iowrite,irsficno,lhead
c      USE PBLCOM
c      USE DOMAIN_DECOMP_ATM, only : grid
c      USE DOMAIN_DECOMP_1D, only : GET, AM_I_ROOT
c      USE DOMAIN_DECOMP_1D, only : pack_column, pack_data
c      USE DOMAIN_DECOMP_1D, only : unpack_column, unpack_data
c      USE DOMAIN_DECOMP_1D, only : pack_block , unpack_block
c#ifdef TRACERS_ON
c      use tracer_com, only : trname
c#endif
c      IMPLICIT NONE
c
c      INTEGER kunit   !@var kunit unit number of read/write
c      INTEGER iaction !@var iaction flag for reading or writing to file
c!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
c      INTEGER, INTENT(INOUT) :: IOERR
c      real*8, dimension(:,:,:,:), allocatable :: ! (npbl,4,IM,JM)
c     &     uabl_glob,vabl_glob,tabl_glob,qabl_glob,eabl_glob
c      REAL*8, DIMENSION(:,:,:), allocatable ::   ! (4,IM,JM)
c     &     cmgs_glob, chgs_glob, cqgs_glob
c      INTEGER, DIMENSION(:,:,:), allocatable ::  ! (4,IM,JM)
c     &     ipbl_glob
c      INTEGER :: J_0, J_1, j_0h, j_1h, n, i_0h, i_1h
c!@var HEADER Character string label for individual records
c      CHARACTER*80 :: HEADER, MODULE_HEADER = "PBL01"
c#ifdef TRACERS_ON
c!@var TR_HEADER Character string label for tracer record
c      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TRPBL01"
c      REAL*8, DIMENSION(:,:,:,:), allocatable :: trabl_glob,trabl_loc
c#endif
c      integer :: img, jmg
c
c#ifdef PBL_USES_GCM_TENDENCIES
c      call stop_model('io_pbl: need to save [tuv]1_after_aturb',255)
c#endif
c
c      write (MODULE_HEADER(lhead+1:80),'(a7,i2,a)') 'R8 dim(',npbl,
c     *  ',4,ijm):Ut,Vt,Tt,Qt,Et dim(4,ijm,3):Cmhq, I:Ipb(4,ijm)'
c
c      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
c     &     J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
c      I_0H = grid%I_STRT_HALO
c      I_1H = grid%I_STOP_HALO
c
c      if(am_i_root()) then
c         img = IM
c         jmg = JM
c      else
c         img = 1
c         jmg = 1
c      end if
c      allocate(uabl_glob(npbl,4,img,jmg))
c      allocate(vabl_glob(npbl,4,img,jmg))
c      allocate(tabl_glob(npbl,4,img,jmg))
c      allocate(qabl_glob(npbl,4,img,jmg))
c      allocate(eabl_glob(npbl,4,img,jmg))
c      allocate(cmgs_glob(4,img,jmg))
c      allocate(chgs_glob(4,img,jmg))
c      allocate(cqgs_glob(4,img,jmg))
c      allocate(ipbl_glob(4,img,jmg))
c#ifdef TRACERS_ON
c      allocate(trabl_glob(npbl,4,im,jm))
c#endif
c#ifdef TRACERS_ON
c      allocate(trabl_loc(npbl,4,i_0h:i_1h,j_0h:j_1h))
c#endif
c
c      SELECT CASE (IACTION)
c      CASE (:IOWRITE)            ! output to standard restart file
c        CALL PACK_BLOCK(grid, uabl, uabl_glob)
c        CALL PACK_BLOCK(grid, vabl, vabl_glob)
c        CALL PACK_BLOCK(grid, tabl, tabl_glob)
c        CALL PACK_BLOCK(grid, qabl, qabl_glob)
c        CALL PACK_BLOCK(grid, eabl, eabl_glob)
c
c        CALL PACK_COLUMN(grid, cmgs, cmgs_glob)
c        CALL PACK_COLUMN(grid, chgs, chgs_glob)
c        CALL PACK_COLUMN(grid, cqgs, cqgs_glob)
c        CALL PACK_COLUMN(grid, ipbl, ipbl_glob)
c
c        IF (AM_I_ROOT()) THEN
c          WRITE (KUNIT,ERR=10) MODULE_HEADER,UABL_GLOB,VABL_GLOB
c     *       ,TABL_GLOB,QABL_GLOB,EABL_GLOB,CMGS_GLOB
c     *       ,CHGS_GLOB,CQGS_GLOB,IPBL_GLOB
c        END IF
c#ifdef TRACERS_ON
c        do n=1,ntm
c          trabl_loc(:,:,:,:) = trabl(:,n,:,:,:)
c          CALL PACK_BLOCK(grid, trabl_loc, trabl_glob)
c          write (TR_MODULE_HEADER(lhead+1:80),'(a8,a8,i2,a)')
c     &         trname(n),
c     &         ' R8 dim(',npbl,',4,ijm):TRt'
c          IF (AM_I_ROOT()) write(kunit,err=10)
c     &         TR_MODULE_HEADER,TRABL_GLOB
c        enddo
c#endif
c
c      CASE (IOREAD:)            ! input from restart file or restart
c        if ( AM_I_ROOT() ) then
c          READ (KUNIT,ERR=10) HEADER,UABL_glob,VABL_glob,TABL_glob,
c     &         QABL_glob,EABL_glob,CMGS_glob,CHGS_glob,CQGS_glob,
c     &         IPBL_glob
c          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
c            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
c            GO TO 10
c          END IF
c        end if
c
c        call UNPACK_BLOCK(grid, uabl_glob, uabl)
c        call UNPACK_BLOCK(grid, vabl_glob, vabl)
c        call UNPACK_BLOCK(grid, tabl_glob, tabl)
c        call UNPACK_BLOCK(grid, qabl_glob, qabl)
c        call UNPACK_BLOCK(grid, eabl_glob, eabl)
c
c        call UNPACK_COLUMN(grid, cmgs_glob, cmgs)
c        call UNPACK_COLUMN(grid, chgs_glob, chgs)
c        call UNPACK_COLUMN(grid, cqgs_glob, cqgs)
c        call UNPACK_COLUMN(grid, ipbl_glob, ipbl)
c#ifdef TRACERS_ON
c        SELECT CASE (IACTION)
c        CASE (IOREAD,IRERUN,IRSFIC,IRSFICNO)    ! restarts
c        do n=1,ntm
c          IF (AM_I_ROOT()) then
c            read(kunit,err=10) TR_HEADER,TRABL_GLOB
c            IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
c              PRINT*,"Discrepancy in tracer module version ",TR_HEADER
c     *             ,TR_MODULE_HEADER
c              GO TO 10
c            END IF
c          ENDIF
c          CALL UNPACK_BLOCK(grid, trabl_glob, trabl_loc)
c          trabl(:,n,:,:,:) = trabl_loc(:,:,:,:)
c        enddo
c        END SELECT
c#endif
c      END SELECT
c      call freemem
c      RETURN
c 10   IOERR=1
c      call freemem
c      RETURN
c      contains
c      subroutine freemem
c      deallocate(uabl_glob)
c      deallocate(vabl_glob)
c      deallocate(tabl_glob)
c      deallocate(qabl_glob)
c      deallocate(eabl_glob)
c      deallocate(cmgs_glob)
c      deallocate(chgs_glob)
c      deallocate(cqgs_glob)
c      deallocate(ipbl_glob)
c#ifdef TRACERS_ON
c        deallocate(trabl_glob)
c#endif
c#ifdef TRACERS_ON
c      deallocate(trabl_loc)
c#endif
c      end subroutine freemem
c      END SUBROUTINE io_pbl
c
c      SUBROUTINE io_bldat(kunit,iaction,ioerr)
c!@sum  io_bldat reads and writes boundary layer data to file
c!@auth Gavin Schmidt
c      USE MODEL_COM, only : ioread,iowrite,lhead
c      USE DOMAIN_DECOMP_ATM, only : grid
c      USE DOMAIN_DECOMP_1D, only : GET, AM_I_ROOT
c      USE DOMAIN_DECOMP_1D, only : UNPACK_DATA, UNPACK_COLUMN
c      USE DOMAIN_DECOMP_1D, only : PACK_DATA, PACK_COLUMN
c      USE PBLCOM
c      IMPLICIT NONE
c
c      INTEGER kunit   !@var kunit unit number of read/write
c      INTEGER iaction !@var iaction flag for reading or writing to file
c!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
c      INTEGER, INTENT(INOUT) :: IOERR
c!@var HEADER Character string label for individual records
c      CHARACTER*80 :: HEADER, MODULE_HEADER = "BLD02"
c      INTEGER :: J_0, J_1
c      REAL*8, DIMENSION(IM,JM) :: wsavg_glob,tsavg_glob,
c     *       qsavg_glob,dclev_glob,usavg_glob,
c     *       vsavg_glob,tauavg_glob,qgavg_glob,tgvavg_glob
c      REAL*8, DIMENSION(LM,IM,JM) :: egcm_glob, w2gcm_glob
c      REAL*8, DIMENSION(4,IM,JM) :: ustar_pbl_glob
c
c      MODULE_HEADER(lhead+1:80) = 'R8 dim(ijm):ws,ts,qs,'//
c     *  'LvlDC,us,vs,tau,u*(4,.),ke;w2(lijm),tgv,qg'
c
c      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)
c
c      SELECT CASE (IACTION)
c      CASE (:IOWRITE)            ! output to standard restart file
c        CALL PACK_DATA(grid, wsavg, wsavg_glob)
c        CALL PACK_DATA(grid, tsavg, tsavg_glob)
c        CALL PACK_DATA(grid, qsavg, qsavg_glob)
c        CALL PACK_DATA(grid, dclev, dclev_glob)
c        CALL PACK_DATA(grid, usavg, usavg_glob)
c        CALL PACK_DATA(grid, vsavg, vsavg_glob)
c        CALL PACK_DATA(grid, tauavg, tauavg_glob)
c        CALL PACK_DATA(grid, tgvavg, tgvavg_glob)
c        CALL PACK_DATA(grid, qgavg, qgavg_glob)
c
c        CALL PACK_COLUMN(grid, ustar_pbl, ustar_pbl_glob)
c        CALL PACK_COLUMN(grid, egcm, egcm_glob)
c        CALL PACK_COLUMN(grid, w2gcm, w2gcm_glob)
c
c        IF (AM_I_ROOT()) THEN
c          WRITE (kunit,err=10) MODULE_HEADER,wsavg_glob,tsavg_glob
c     *         ,qsavg_glob,dclev_glob,usavg_glob,vsavg_glob,tauavg_glob
c     *         ,ustar_pbl_glob,egcm_glob,w2gcm_glob,tgvavg_glob
c     *         ,qgavg_glob
c        END IF
c
c      CASE (IOREAD:)            ! input from restart file
c        if ( AM_I_ROOT() ) then
c          READ (kunit,err=10) HEADER,wsavg_glob,tsavg_glob,
c     *       qsavg_glob,dclev_glob,usavg_glob,
c     *       vsavg_glob,tauavg_glob,ustar_pbl_glob,egcm_glob,
c     *       w2gcm_glob,tgvavg_glob,qgavg_glob
c          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
c            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
c            GO TO 10
c          END IF
c        end if
c
c        call UNPACK_DATA(grid, wsavg_glob, wsavg)
c        call UNPACK_DATA(grid, tsavg_glob, tsavg)
c        call UNPACK_DATA(grid, qsavg_glob, qsavg)
c        call UNPACK_DATA(grid, dclev_glob, dclev)
c        call UNPACK_DATA(grid, usavg_glob, usavg)
c        call UNPACK_DATA(grid, vsavg_glob, vsavg)
c        call UNPACK_DATA(grid, tauavg_glob, tauavg)
c        call UNPACK_DATA(grid, tgvavg_glob, tgvavg)
c        call UNPACK_DATA(grid, qgavg_glob, qgavg)
c
c        call UNPACK_COLUMN(grid, ustar_pbl_glob, ustar_pbl)
c        call UNPACK_COLUMN(grid, egcm_glob, egcm)
c        call UNPACK_COLUMN(grid, w2gcm_glob, w2gcm)
c
c      END SELECT
c
c      RETURN
c 10   IOERR=1
c      RETURN
c      END SUBROUTINE io_bldat

#ifdef NEW_IO
      subroutine def_rsf_pbl(fid)
!@sum  def_rsf_pbl defines pbl array structures in restart files
!@auth M. Kelley
!@ver  beta
      use pblcom
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      character(len=29) :: dimstr
      dimstr='(dist_im,dist_jm)'
      call defvar(grid,fid,t1_after_aturb,'t1_after_aturb'//dimstr)
      call defvar(grid,fid,u1_after_aturb,'u1_after_aturb'//dimstr)
      call defvar(grid,fid,v1_after_aturb,'v1_after_aturb'//dimstr)
      return
      end subroutine def_rsf_pbl

      subroutine new_io_pbl(fid,iaction)
!@sum  new_io_pbl read/write pbl arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use pblcom
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 't1_after_aturb',t1_after_aturb)
        call write_dist_data(grid, fid, 'u1_after_aturb',u1_after_aturb)
        call write_dist_data(grid, fid, 'v1_after_aturb',v1_after_aturb)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 't1_after_aturb',t1_after_aturb)
        call read_dist_data(grid, fid, 'u1_after_aturb',u1_after_aturb)
        call read_dist_data(grid, fid, 'v1_after_aturb',v1_after_aturb)
      end select
      return
      end subroutine new_io_pbl

      subroutine def_rsf_bldat(fid)
!@sum  def_rsf_bldat defines bldat array structure in restart files
!@auth M. Kelley
!@ver  beta
      use pblcom
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,dclev,'dclev(dist_im,dist_jm)')
      call defvar(grid,fid,pblht,'pblht(dist_im,dist_jm)')
      call defvar(grid,fid,pblptop,'pblptop(dist_im,dist_jm)')
      call defvar(grid,fid,egcm,'egcm(lm,dist_im,dist_jm)')
      call defvar(grid,fid,w2gcm,'w2gcm(lm,dist_im,dist_jm)')
      return
      end subroutine def_rsf_bldat

      subroutine new_io_bldat(fid,iaction)
!@sum  new_io_bldat read/write bldat arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use pblcom
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'dclev',dclev)
        call write_dist_data(grid,fid,'pblht',pblht)
        call write_dist_data(grid,fid,'pblptop',pblptop)
        call write_dist_data(grid,fid,'egcm',egcm, jdim=3)
        call write_dist_data(grid,fid,'w2gcm',w2gcm, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'dclev',dclev)
        call read_dist_data(grid,fid,'pblht',pblht)
        call read_dist_data(grid,fid,'pblptop',pblptop)
        call read_dist_data(grid,fid,'egcm',egcm, jdim=3)
        call read_dist_data(grid,fid,'w2gcm',w2gcm, jdim=3)
      end select
      return
      end subroutine new_io_bldat
#endif /* NEW_IO */

      SUBROUTINE ALLOC_PBL_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE RESOLUTION, only : lm
      USE CONSTANT, only : by3
      USE PBLCOM
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, getDomainBounds
      USE FLUXES, only : atmocns,atmices,atmglas,atmlnds,asflx
#ifdef GLINT2
      USE FLUXES, only : atmglas_hp
#endif

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER, IP3, K, L

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE(    roughl(I_0H:I_1H,J_0H:J_1H),
     *              dclev(I_0H:I_1H,J_0H:J_1H),
     *              pblht(I_0H:I_1H,J_0H:J_1H),
     *              pblptop(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

C**** SET LAYER THROUGH WHICH DRY CONVECTION MIXES TO 1
      DCLEV(:,:)=1.

      ALLOCATE(    egcm(lm,I_0H:I_1H,J_0H:J_1H),
     *            w2gcm(lm,I_0H:I_1H,J_0H:J_1H),
     *            t2gcm(lm,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(
     &         ugeo(I_0H:I_1H,J_0H:J_1H),
     &         vgeo(I_0H:I_1H,J_0H:J_1H),
     &         bldep(I_0H:I_1H,J_0H:J_1H)
     &     )

      ip3 = 0
      do k=1,size(asflx)
        select case(asflx(k)%itype4)
        case (1)
          call alloc_atmsrf_pbl_vars(grid,
     &         atmocns(1)%atmsrf_xchng_vars,asflx(k))
        case (2)
          call alloc_atmsrf_pbl_vars(grid,
     &         atmices(1)%atmsrf_xchng_vars,asflx(k))
        case (3)
          ip3 = ip3 + 1
#ifdef GLINT2
          call alloc_atmsrf_pbl_vars(grid,
     &         atmglas_hp(ip3)%atmsrf_xchng_vars,asflx(k))
#endif
      ! dest in alloc_atmsrf_pbl_vars() should point to atmglas
      ! in the end, not atmglas_hp
          call alloc_atmsrf_pbl_vars(grid,
     &         atmglas(ip3)%atmsrf_xchng_vars,asflx(k))
        case (4)
          call alloc_atmsrf_pbl_vars(grid,
     &         atmlnds(1)%atmsrf_xchng_vars,asflx(k))
        end select
      enddo

      ALLOCATE(t1_after_aturb(I_0H:I_1H,J_0H:J_1H),
     &         u1_after_aturb(I_0H:I_1H,J_0H:J_1H),
     &         v1_after_aturb(I_0H:I_1H,J_0H:J_1H))
      t1_after_aturb(:,J_0H:J_1H) = 0.
      u1_after_aturb(:,J_0H:J_1H) = 0.
      v1_after_aturb(:,J_0H:J_1H) = 0.

C**** initialize egcm to be used in ATURB.f
      DO L=1,LM
        egcm(l,:,:)=egcm_init_max/(float(l)**2)
        w2gcm(l,:,:)=egcm(l,:,:)*2.*by3
      END DO

      END SUBROUTINE ALLOC_PBL_COM

      subroutine alloc_atmsrf_pbl_vars(grid,this,dest)
      use pblcom
      use exchange_types, only : atmsrf_xchng_vars
      use domain_decomp_1d, only : dist_grid, getDomainBounds
      implicit none
      type (dist_grid), intent(in) :: grid
      type(atmsrf_xchng_vars) :: this,dest
      integer :: i_0h, i_1h, j_1h, j_0h

      call getDomainBounds(grid, j_strt_halo=j_0h, j_stop_halo=j_1h)
      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo

      allocate(
     &     this % uabl(npbl,i_0h:i_1h,j_0h:j_1h)
     &    ,this % vabl(npbl,i_0h:i_1h,j_0h:j_1h)
     &    ,this % tabl(npbl,i_0h:i_1h,j_0h:j_1h)
     &    ,this % qabl(npbl,i_0h:i_1h,j_0h:j_1h)
     &    ,this % eabl(npbl,i_0h:i_1h,j_0h:j_1h)
#ifdef TRACERS_ON
     &    ,this % trabl(npbl,this%ntm,i_0h:i_1h,j_0h:j_1h)
#endif
     &     )
      this % qabl = 0.  ! initialise to make life easier

      dest%uabl => this%uabl
      dest%vabl => this%vabl
      dest%tabl => this%tabl
      dest%qabl => this%qabl
      dest%eabl => this%eabl
#ifdef TRACERS_ON
      dest%trabl => this%trabl
#endif
      return
      end subroutine alloc_atmsrf_pbl_vars
