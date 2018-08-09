#include "rundeck_opts.h"

      subroutine diag_ocean_prep
      use oceanr_dim, only : grid=>ogrid
      implicit none
      if(.not. grid%have_domain) return
      call basin_prep
      call oijl_prep
      return
      end subroutine diag_ocean_prep

      SUBROUTINE diag_OCEAN
#ifndef NEW_IO
!@sum  diag_OCEAN prints out diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
C**** Note this is an incorporation and modification of the stand alone
C**** ocean diagnostic programs from Gary. All diagnostics are on the
C**** ocean grid.
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,amon,modelEclock
      USE OCEAN, only : im,jm,lmo,dts,imaxj,lmm,ze
      USE DIAG_COM, only : qdiag,zoc_pout=>zoc,zoc1_pout=>zoc1
      USE MDIAG_COM, only : acc_period
      USE ODIAG
      USE FILEMANAGER, only : openunit
      IMPLICIT NONE
      INTEGER I,J,L,N,NOL(LMO),JEQ,JDLAT,KXLB
      REAL*8 DLON

      zoc_pout(1:lmo) = zoc(1:lmo)
      zoc1_pout(1:lmo+1) = zoc1(1:lmo+1)

C**** Calculate latitudes
      JEQ = JM/2
      JDLAT = NINT(180./(JM-1))
      DO J=1,JM
        FLAT(J,1)  = JDLAT*(J-JEQ-0.5)  ! primary grid
      END DO
      FLAT(1,2)=-90.
      DO J=1,JM-1
        FLAT(J+1,2)  = JDLAT*(J-JEQ)  ! secondary grid  (shifted by 1)
      END DO
C**** Calculate longitudes
      DLON=360./REAL(IM,KIND=8)
      DO I=1,IM
        FLON(I,1)  = -180.+DLON*(I-0.5)  ! primary grid
      END DO
      DO I=1,IM
        FLON(I,2)  = -180.+DLON*I        ! secondary grid
      END DO

C**** determine label to be added to all titles
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB = TRIM(XLB)//" "//XLABEL(1:KXLB)

C**** Open output files
      IF(QDIAG) then
        call open_ij(trim(acc_period)//'.oij'//XLABEL(1:LRUNID),im,jm)
        call open_jl(trim(acc_period)//'.ojl'//XLABEL(1:LRUNID),jm,lmo,0
     *       ,flat)
        call open_il(trim(acc_period)//'.oil'//XLABEL(1:LRUNID),im,lmo,0
     *       )
        call openunit(trim(acc_period)//'.otj'//XLABEL(1:LRUNID),iu_otj,
     *       .FALSE.,.FALSE.)
      END IF
C**** Call diagnostic calculating routines
      CALL OIJOUT  ! lat-lon diags
      CALL OSFOUT  ! overturning streamfunction diags
      CALL OJLOUT  ! lat-height and lon-height sections and means
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      CALL OTJOUT  ! transports
      CALL STROUT  ! strait diags
C**** Miscellaneous diags for print out

C****
C**** vertical mean diagnostics
C****
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
C**** Calculate number of points in average
      NOL = 0.
      DO J=1,JM
      DO I=1,IMAXJ(J)
      DO L=1,LMM(I,J)
        NOL(L)=NOL(L)+1
      END DO
      END DO
      END DO

      WRITE(6,*) " Ocean Mean quantities:"
      WRITE(6,*) " Level   Rho0       Temp    Salinity"
      DO L=1,LMO
        WRITE(6,'(2X,I4,10F10.3)') L,OL(L,L_RHO)/(NOL(L)*IDACC(1)),
     *           OL(L,L_TEMP)/(NOL(L)*IDACC(1)),
     *       1d3*OL(L,L_SALT)/(NOL(L)*IDACC(1))
      END DO
      call sys_flush(6)
C****
      IF (QDIAG) THEN
        call close_ij
        call close_jl
        call close_il
      END IF
C****
      RETURN
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
#endif
      END SUBROUTINE diag_OCEAN

#ifndef NEW_IO
      SUBROUTINE OIJOUT
!@sum  OIJOUT prints out lat-lon diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : undef,teeny,rhows
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,amon,modelEclock
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM,only:n_water,n_obio,tracerlist,ocn_tracer_entry
#endif
      USE OCEAN, only: im,jm,lmo,focean,dxypo,dts,imaxj,lmm,ze,dxvo,dypo
      USE DIAG_COM, only : qdiag
      USE MDIAG_COM, only :
     &     sname_strlen,units_strlen,lname_strlen
      USE STRAITS, only : nmst,wist,dist,lmst,name_st

      USE OCEAN, only : oDLAT_DG, oLAT_DG, oLON_DG

      USE ODIAG
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM) :: Q,SFIJ,QS,QS700
      REAL*8, DIMENSION(IM,JM,LMO) :: Q3
c for now we are assuming that igrid=jgrid in arguments to pout_ij
      INTEGER I,J,K,L,N,KB,IJGRID,IP1,k1,KK
      INTEGER :: LMINEF=1, LMAXEF=LMO, KVMF(3) = (/ 3, 6, 9/),
     *           LMINMF=1, LMAXMF=1,   KCMF(3) = (/ 3, 6, 9/),
     *           LMINSF=1, LMAXSF=LMO, KVDC(3) = (/ 3, 6, 9/)
     *     ,KCMFfull(13) = (/1,2,3,4,5,6,7,8,9,10,11,12,13/)
      REAL*8 GOS,SOS,FAC,FACST,GSMAX,GSMIN,CKMIN,CKMAX,ACMIN,ACMAX
     *     ,TSUM,TEMGS,QJ(JM),QSUM,MQ,DLON,byiacc,volgs
      CHARACTER NAME(KOLNST)*40,TITLE*80
      CHARACTER(len=lname_strlen) :: lname
      CHARACTER(len=sname_strlen) :: sname
      CHARACTER(len=units_strlen) :: units
      character*50 :: unit_string
      character(len=2), dimension(lmo) :: levstr
      real*8 :: scale_jm2
      type(ocn_tracer_entry), pointer :: entry

      QJ=0.
      QSUM=0.
      IJGRID=1
      Q=0.
      Q3=0.

C**** lat/lon diagnostics
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)

      IF(QDIAG) then

      DO L=1,LMO
        WRITE(LEVSTR(L),'(I2.2)') L
      ENDDO

C****
C**** Ocean Potential Temperature (C)
C****
      K=IJL_PTM
      LNAME=LNAME_OIJL(K)
      UNITS=UNITS_OIJL(K)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)  THEN
          GOS = OIJL(I,J,L,IJL_G0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          SOS = OIJL(I,J,L,IJL_S0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          Q(I,J) = TEMGS(GOS,SOS)
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME = 'oc_temp_L'//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** Ocean Salinity (per mil)
C****
      K=IJL_S0M
      LNAME=LNAME_OIJL(K)
      UNITS=UNITS_OIJL(K)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)
     *       Q(I,J) = OIJL(I,J,L,K)*SCALE_OIJL(K) /
     &       (OIJL(I,J,L,IJL_MO)*DXYPO(J))
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME = 'oc_salt_L'//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** Ocean Potential Density (kg/m^3) (w.r.t. 0m)
C****
      K=IJL_PDM
      LNAME=LNAME_OIJL(K)
      UNITS=UNITS_OIJL(K)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)  THEN
          GOS = OIJL(I,J,L,IJL_G0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          SOS = OIJL(I,J,L,IJL_S0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          Q(I,J) = 1./VOLGS(GOS,SOS) - 1000.
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME = 'oc_pot_den_L'//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
#ifdef TRACERS_OCEAN
C****
C**** Ocean Tracers
C****
      DO N=1,tracerlist%getsize()
      entry=>tracerlist%at(n)
      LNAME="OCEAN "//entry%trname
      if (entry%to_per_mil.gt.0) THEN
        UNITS="per mil"
      ELSE
        UNITS=unit_string(entry%ntrocn,'kg/kg')
      END IF
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.) THEN
          if (entry%to_per_mil.gt.0 .and. TOIJL(I,J,L,TOIJL_CONC
     *         ,n_water).gt.0) THEN
            Q(I,J)=1d3*(TOIJL(I,J,L,TOIJL_CONC,N)/(TOIJL(I,J,L
     *           ,TOIJL_CONC,n_water)*entry%trw0)-1.)
c          Q(I,J)=1d3*(TOIJL(I,J,L,TOIJL_CONC,N)/((OIJL(I,J,L,IJL_MO)
c     *         *DXYPO(J)-OIJL(I,J,L,IJL_S0M))*entry%trw0)-1.)
          else
            Q(I,J)=10.**(-entry%ntrocn)*TOIJL(I,J,L,TOIJL_CONC,N)/
     *           (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          end if
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME='oc_'//trim(entry%trname)//'_L'//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
      END DO
#endif
C****
C**** East-West or North-South Velocities (cm/s)
C****
      K=IJL_MFU
      DO LMINMF=1,LMO
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        LMAXMF=LMINMF
      Q = 0.
      DO J=1,JM
        I=IM
        DO IP1=1,IMAXJ(J)
          MQ = 0.
          DO L=LMINMF,LMAXMF
            MQ = MQ + (OIJL(I,J,L,IJL_MO)+OIJL(IP1,J,L,IJL_MO))
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = SCALE_OIJL(K) * Q(I,J) / (.5*MQ*DYPO(J)+teeny)
          I=IP1
        END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) THEN
        WRITE (TITLE(39:41),'(I3)') LMINMF
        WRITE (LNAME(39:41),'(I3)') LMINMF
      ELSEIF(LMINMF.lt.LMAXMF) THEN
        WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
        WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
      END IF
      SNAME="uvel_L"//LEVSTR(LMINMF)
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO

      K=IJL_MFV
      DO LMINMF=1,LMO
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        LMAXMF=LMINMF
      Q = 0.
      DO J=1,JM-1
        DO I=1,IMAXJ(J)
          MQ = 0.
          DO L=LMINMF,LMAXMF
            MQ = MQ + (OIJL(I,J,L,IJL_MO)+OIJL(I,J+1,L,IJL_MO))
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) =  SCALE_OIJL(K) * Q(I,J) / (.5*MQ*DXVO(J)+teeny)
        END DO
      END DO
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) THEN
        WRITE (TITLE(39:41),'(I3)') LMINMF
        WRITE (LNAME(39:41),'(I3)') LMINMF
      ELSEIF(LMINMF.lt.LMAXMF) THEN
        WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
        WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
      END IF
      SNAME="vvel_L"//LEVSTR(LMINMF)
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO
C****
C**** Vertical Velocity (cm/s)
C****
      K=IJL_MFW
      DO L=1,LMO-1
c     IF(KVMF(K).le.0)  GO TO 370
c        L =KVMF(K)
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME='vert_vel_L'//LEVSTR(L)
        DO J=1,JM
          DO I=1,IMAXJ(J)
            Q(I,J)=OIJL(I,J,L,K)*SCALE_OIJL(K)/(IDACC(1)*DXYPO(J))
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** Vertical Mass Flux^2 (kg2/m4)
C****
      K=IJL_MFW2
      DO L=1,LMO-1
c     IF(KVMF(K).le.0)  GO TO 370
c        L =KVMF(K)
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME='vert_ms_flx_sq_L'//LEVSTR(L)
        DO J=1,JM
          DO I=1,IMAXJ(J)
            Q(I,J)=OIJL(I,J,L,K)*SCALE_OIJL(K)/
     *           (IDACC(1)*DXYPO(J)*DXYPO(J))
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** East-West or North-South Mass fluxes (kg/s)
C****
      K=IJL_MFU
      DO LMINMF=1,LMO
      LNAME="EAST-WEST MASS FLUX"
      UNITS="10^9 kg/s"
        LMAXMF=LMINMF
      Q = 0.
      DO J=1,JM
        I=IM
        DO IP1=1,IMAXJ(J)
          DO L=LMINMF,LMAXMF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = 1d-9* Q(I,J) / (IDACC(1)*DTS+teeny)
          I=IP1
        END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) THEN
        WRITE (TITLE(40:47),'(A5,I3)') "Level",LMINMF ! anl
        WRITE (LNAME(40:47),'(A5,I3)') "Level",LMINMF ! anl
      ELSEIF(LMINMF.lt.LMAXMF) THEN
        WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
        WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
      END IF
      SNAME="umfl_L"//LEVSTR(LMINMF)
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO

      K=IJL_MFV
      DO LMINMF=1,LMO
      LNAME="NORTH-SOUTH MASS FLUX"
      UNITS="10^9 kg/s"
        LMAXMF=LMINMF
      Q = 0.
      DO J=1,JM-1
        DO I=1,IMAXJ(J)
          DO L=LMINMF,LMAXMF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = 1d-9* Q(I,J) / (IDACC(1)*DTS+teeny)
        END DO
      END DO
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) THEN
        WRITE (TITLE(40:47),'(A5,I3)') "Level",LMINMF ! anl
        WRITE (LNAME(40:47),'(A5,I3)') "Level",LMINMF ! anl
      ELSEIF(LMINMF.lt.LMAXMF) THEN
        WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
        WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
      END IF
      SNAME="vmfl_L"//LEVSTR(LMINMF)
c        WRITE (TITLE(40:47),'(A5,I3)') "Level",L ! anl
c        WRITE (LNAME(40:47),'(A5,I3)') "Level",L ! anl
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO
C****
C**** Vertical Mass Flux (kg/s)
C****
      DO L=1,LMO-1
c     IF(KVMF(K).le.0)  GO TO 370
c        L =KVMF(K)
        LNAME="DOWNWARD VERTICAL MASS FLUX"
        UNITS="10^9 kg/s"
        SNAME='vert_mfl_L'//LEVSTR(L)
        DO J=1,JM
          DO I=1,IMAXJ(J)
            Q(I,J) = 1d-9*OIJL(I,J,L,IJL_MFW) / real(IDACC(1)*DTS)
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** East-West or North-South Heat Flux (10^15 W)
C****
      DO K=IJL_GFLX,IJL_GFLX+1
c      IF(.not.QL(K))  GO TO 440
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        IF (K.eq.IJL_GFLX) THEN
          SNAME="ew_hflx"
        ELSE
          SNAME="ns_hflx"
        END IF
        Q = 0.
        DO J=1,JM
        DO I=1,IMAXJ(J)
          DO L=LMINEF,LMAXEF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = SCALE_OIJL(K)*Q(I,J)/IDACC(1)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        IF(LMINEF.eq.LMAXEF) WRITE (TITLE(39:41),'(I3)') LMINEF
        IF(LMINEF.lt.LMAXEF) WRITE (TITLE(39:46),'(I3,A2,I3)') LMINEF,
     *       " -",LMAXEF
        IF(LMINEF.eq.LMAXEF) WRITE (LNAME(39:41),'(I3)') LMINEF
        IF(LMINEF.lt.LMAXEF) WRITE (LNAME(39:46),'(I3,A2,I3)') LMINEF,
     *       " -",LMAXEF
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** East-West or North-South Salt Flux (10^6 kg/s)
C****
      DO K=IJL_SFLX,IJL_SFLX+1
c      IF(.not.QL(K))  GO TO 540
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        IF (K.eq.IJL_SFLX) THEN
          SNAME="ew_sflx"
        ELSE
          SNAME="ns_sflx"
        END IF
        Q = 0.
        DO J=1,JM
        DO I=1,IMAXJ(J)
          DO L=LMINSF,LMAXSF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = SCALE_OIJL(K)*Q(I,J)/IDACC(1)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        IF(LMINSF.eq.LMAXSF) WRITE (TITLE(39:41),'(I3)') LMINSF
        IF(LMINSF.lt.LMAXSF) WRITE (TITLE(39:46),'(I3,A2,I3)') LMINSF,
     *       " -",LMAXSF
        IF(LMINSF.eq.LMAXSF) WRITE (LNAME(39:41),'(I3)') LMINSF
        IF(LMINSF.lt.LMAXSF) WRITE (LNAME(39:46),'(I3,A2,I3)') LMINSF,
     *       " -",LMAXSF
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** Gent-McWilliams fluxes ([J or kg]/s or [J kg]/m2/s for vert)
C****
      DO KK=0,2
        K=KK+IJL_GGMFL
        DO L=1,lmo 
          LNAME=LNAME_OIJL(K)
          UNITS=UNITS_OIJL(K)
          SELECT CASE (KK)
          CASE (0)      ! E-W fluxes
            SNAME='gm_ew_hflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/IDACC(1)
              END DO
            END DO
          CASE (1)  ! N-S fluxes
            SNAME='gm_ns_hflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/IDACC(1)
              END DO
            END DO
          CASE (2)    !  Vertical fluxes
            SNAME='gm_vt_hflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/(IDACC(1)*DXYPO(J))
              END DO
            END DO
          END SELECT
          Q(2:IM,JM)=Q(1,JM)
          Q(2:IM,1)=Q(1,1)
          TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
          WRITE (TITLE(40:47),'(A5,I3)') "Level",L
          WRITE (LNAME(40:47),'(A5,I3)') "Level",L
          TITLE(51:80)=XLB
          CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
        END DO
      END DO
C****
C**** Gent-McWilliams Salt Fluxes
C****
      DO KK=0,2
        K=KK+IJL_SGMFL
        DO L=1,lmo
          LNAME=LNAME_OIJL(K)
          UNITS=UNITS_OIJL(K)
          SELECT CASE (KK)
          CASE (0)      ! E-W fluxes
            SNAME='gm_ew_sflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/IDACC(1)
              END DO
            END DO
          CASE (1)  ! N-S fluxes
            SNAME='gm_ns_sflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/IDACC(1)
              END DO
            END DO
          CASE (2)    !  Vertical fluxes
            SNAME='gm_vt_sflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/(IDACC(1)*DXYPO(J))
              END DO
            END DO
          END SELECT

          Q(2:IM,JM)=Q(1,JM)
          Q(2:IM,1)=Q(1,1)
          TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
          WRITE (TITLE(40:47),'(A5,I3)') "Level",L
          WRITE (LNAME(40:47),'(A5,I3)') "Level",L
          TITLE(51:80)=XLB
          CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
        END DO
      END DO
#ifdef TRACERS_OCEAN
C****
C**** Gent-McWilliams Tracer Fluxes
C****
      do n=1,tracerlist%getsize()
      entry=>tracerlist%at(n)
      DO KK=0,2
        DO L=1,lmo
          SELECT CASE (KK)
          CASE (0)      ! E-W fluxes
            LNAME="GM/EDDY E-W FLUX "//entry%trname
            UNITS=unit_string(entry%ntrocn,'kg/s')
            SNAME="oc_gm_ewflx"//trim(entry%trname)//"_L"//LEVSTR(L)
          CASE (1)  ! N-S fluxes
            LNAME="GM/EDDY N-S FLUX "//entry%trname
            UNITS=unit_string(entry%ntrocn,'kg/s')
            SNAME="oc_gm_nstflx"//trim(entry%trname)//"_L"//LEVSTR(L)
          CASE (2)    !  Vertical fluxes
            LNAME="GM/EDDY VERT. FLUX "//entry%trname
            UNITS=unit_string(entry%ntrocn-6,'kg/m^2 s')
            SNAME="gm_vt_tflx"//trim(entry%trname)//"_L"//LEVSTR(L)
          END SELECT
          if (KK.EQ.3) THEN     ! vert fluxes, scale/divide by area

          DO J=1,JM
            DO I=1,IMAXJ(J)
              IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)
     *             THEN
                if (TOIJL(I,J,L,TOIJL_CONC,N).gt.0) Q(I,J)=10.**
     *                (6-entry%ntrocn)*TOIJL(I,J,L,KK+TOIJL_GMFL,N)/
     *                (IDACC(1)*DTS*DXYPO(J))
              ENDIF
            END DO
          END DO

          ELSE
          
          DO J=1,JM
            DO I=1,IMAXJ(J)
              IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)
     *             THEN
                if (TOIJL(I,J,L,TOIJL_CONC,N).gt.0) Q(I,J)=10.**
     *            (-entry%ntrocn)*TOIJL(I,J,L,KK+TOIJL_GMFL,N)/(IDACC(1)
     *               *DTS)
              ENDIF
            END DO
          END DO

          END IF
          Q(2:IM,JM)=Q(1,JM)
          Q(2:IM,1)=Q(1,1)
          TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
          WRITE (TITLE(40:47),'(A5,I3)') "Level",L
          WRITE (LNAME(40:47),'(A5,I3)') "Level",L
          TITLE(51:80)=XLB
          CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
        END DO
      END DO
      enddo
#endif
C****
C**** Vertical Diffusion Coefficients (cm/s)
C****
      K=IJL_KVM
      DO L=1,lmo-1
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME="kvm_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (OIJL(I,J,L+1,IJL_MO).gt.0) Q(I,J)=OIJL(I,J,L,K)*
     &         SCALE_OIJL(K)/IDACC(1)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO

      K=IJL_KVG
      DO L=1,lmo-1
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME="kvg_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (OIJL(I,J,L+1,IJL_MO).gt.0) Q(I,J)=OIJL(I,J,L,K)*
     &         SCALE_OIJL(K)/IDACC(1)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO

      K=IJL_WGFL
      DO L=1,lmo-1
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME="wgfl_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
          DO I=1,IMAXJ(J)
            IF (OIJL(I,J,L+1,IJL_MO).gt.0) Q(I,J)=OIJL(I,J,L,K)
     &           *SCALE_OIJL(K)/(IDACC(1)*dxypo(j))
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO

      K=IJL_WSFL
      DO L=1,lmo-1
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME="wsfl_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (OIJL(I,J,L+1,IJL_MO).gt.0) Q(I,J)=OIJL(I,J,L,K)
     &         *SCALE_OIJL(K)/(IDACC(1)*dxypo(j))
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO

#ifdef TRACERS_OCEAN
      do n=1,tracerlist%getsize()
      entry=>tracerlist%at(n)
      DO L=1,lmo-1
        LNAME="VERT. DIFF. FLUX "//entry%trname
        UNITS=unit_string(entry%ntrocn-6,'kg/m^2 s')
        SNAME="wtfltr"//trim(entry%trname)//"_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF((OIJL(I,J,L+1,IJL_MO).gt.0).and.(TOIJL(I,J,L,TOIJL_CONC
     *         ,N).gt.0)) Q(I,J)=10.**(6-entry%ntrocn)*TOIJL(I,J,L
     *         ,TOIJL_wtfl,N)/(IDACC(1)*DTS*DXYPO(J))
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO
      enddo
#endif

C****
C**** Simple scaled OIJ diagnostics
C****
      DO K=1,KOIJ
        if(trim(sname_oij(k)).eq.'unused') cycle
        byiacc=1./(IDACC(IA_OIJ(K))+teeny)
        if(igrid_oij(k).eq.1 .and. jgrid_oij(k).eq.1) then
          Q=UNDEF
          DO J=1,JM
            DO I=1,IMAXJ(J)
              IF (FOCEAN(I,J).gt.0.5)
     *             Q(I,J)=SCALE_OIJ(K)*OIJ(I,J,K)*byiacc
            END DO
          END DO
          Q(2:IM,JM)=Q(1,JM)
          Q(2:IM,1)=Q(1,1)
        else ! horizontal fluxes
          Q=SCALE_OIJ(K)*OIJ(:,:,K)*byiacc
        endif
        lname=lname_oij(k)
        TITLE=trim(LNAME)//" ("//trim(UNITS_OIJ(K))//") "
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME_OIJ(K),LNAME_OIJ(K),UNITS_OIJ(K),Q,QJ
     *       ,QSUM,IGRID_OIJ(K),JGRID_OIJ(K))

      END DO

      END IF

C****
C**** Calculate Salt Stream Function and write it
C****
c      FAC   = -1d-6/(IDACC(1)*DTS)
c      FACST = -.5E-6/(IDACC(1)*DTS)
c      CALL STRMIJ (OIJL(1,1,1,IJL_SFLX),FAC,OLNST(1,1,LN_SFLX),FACST
c     *     ,SFIJS)
C**** Subtract .035 times the Mass Stream Function
c      DO J=1,JM
c        DO I=1,IM
c        IF (SFIJS(I,J).ne.UNDEF) SFIJS(I,J) = SFIJS(I,J)-35.*SFIJM(I,J)
c        END DO
c      END DO
c      TITLE = NAME(2)//'  Run '//XLABEL(1:6)
c      WRITE(TITLE(63:80),'(A6,I4)') JMON0,JYEAR0

C**** Output Key diagnostics: Gulf Stream, ACC, Kuroshio
      SFIJ = OIJ(:,:,IJ_SF)/idacc(1)
      WRITE(6,'(A)') " Key horizontal mass stream function diags:"

      GSMAX=0 ; GSMIN=100.
      CKMAX=0 ; CKMIN=100.
      ACMAX=0 ; ACMIN=100.
      DO J=2,JM-1
        DO I=1,IM
          if (oLAT_DG(j,2).ge.24 .and. oLAT_DG(j,2).le.38 .and. oLON_DG
     $         (i,2).ge.-80 .and. oLON_DG(i,2).le.-60) then
            IF (SFIJ(I,J).GT.GSMAX) GSMAX=SFIJ(I,J)
            IF (SFIJ(I,J).LT.GSMIN) GSMIN=SFIJ(I,J)
          end if
          if (oLAT_DG(j,2).ge.24 .and. oLAT_DG(j,2).le.38 .and. oLON_DG
     $         (i,2).ge.135 .and. oLON_DG(i,2).le.155) then
            IF (SFIJ(I,J).GT.CKMAX) CKMAX=SFIJ(I,J)
            IF (SFIJ(I,J).LT.CKMIN) CKMIN=SFIJ(I,J)
          end if
          if (oLAT_DG(j,2).ge.-72 .and. oLAT_DG(j,2).le.-54
     $         .and. oLON_DG(i,2).ge.-65-0.5*oDLAT_DG .and. oLON_DG(i,2)
     $         .le. -65+0.5*oDLAT_DG) then
            IF (SFIJ(I,J).GT.ACMAX) ACMAX=SFIJ(I,J)
            IF (SFIJ(I,J).LT.ACMIN) ACMIN=SFIJ(I,J)
          end if
        END DO
      END DO

      WRITE(6,'(a,F10.3)') " Gulf Stream (MAX in (24-38N,60-80W)):",
     *     GSMAX-GSMIN
      WRITE(6,'(a,F10.3)') " Kuroshio  (MAX in (24-38N,135-150E)):",
     *     CKMAX-CKMIN
      WRITE(6,'(a,F10.3)') " ACC                 (Drakes Passage):",
     *     ACMAX-ACMIN
C****
      IF (QDIAG) THEN
C****
C**** Ocean Heat Content (J/m^2)
C****
      K=IJL_G0M
      LNAME=LNAME_OIJL(K)
      UNITS='10^6 J/m^2' ! UNITS_OIJL(K) is J/kg and SCALE_OIJL(K) is 1
      scale_jm2 = 1d-6
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      QS=0. ; QS700=0.
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)  THEN
          Q(I,J) = SCALE_JM2*OIJL(I,J,L,K) / (DXYPO(J)*IDACC(1))
          QS(I,J) = QS(I,J)+OIJL(I,J,L,K) / DXYPO(J)
          IF (ZE(L-1).LE.700) QS700(I,J)=QS700(I,J)+OIJL(I,J,L,K)
     *         /DXYPO(J)
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME="oc_ht_L"//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C**** sum to 700m
      QS700(2:IM,JM)=QS700(1,JM)
      QS700(2:IM,1)=QS700(1,1)
      QS700 = SCALE_JM2*QS700/REAL(IDACC(1))
      WRITE (LNAME(40:47),'(A8)') 'SUM 700m'
      WRITE (TITLE(40:47),'(A8)') 'SUM 700m'
      SNAME="oc_ht_700m"
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,QS700,QJ,QSUM,IJGRID,IJGRID)
C**** total sum
      QS(2:IM,JM)=QS(1,JM)
      QS(2:IM,1)=QS(1,1)
      QS = SCALE_JM2*QS/REAL(IDACC(1))
      WRITE (LNAME(40:47),'(A8)') '  SUM  '
      WRITE (TITLE(40:47),'(A8)') '  SUM  '
      SNAME="oc_ht_total"
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,QS,QJ,QSUM,IJGRID,IJGRID)
      END IF
C****
      RETURN
      END SUBROUTINE OIJOUT

      SUBROUTINE OSFOUT
!@sum  OSFOUT prints out streamfunction diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : undef
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,amon,modelEclock
      USE OCEAN, only : im,jm,lmo,dxypo,imaxj,ze
      USE DIAG_COM, only : qdiag
      USE MDIAG_COM, only :
     &     sname_strlen,units_strlen,lname_strlen

      USE OCEAN, only : oDLAT_DG, oLAT_DG

      USE ODIAG
      IMPLICIT NONE
      REAL*8, DIMENSION(JM+3,LMO+1) :: XJL
      INTEGER I,J,K,L,NS,N,KB,IP1,LP1,JEQ
      REAL*8 FAC,FACST
      CHARACTER TITLE*80
      CHARACTER(len=lname_strlen) :: lname
      CHARACTER(len=sname_strlen) :: sname
      CHARACTER(len=units_strlen) :: units

      TITLE(51:80)=XLB
      LNAME=""
      XJL=undef
C****
C**** Write the Mass Stream Function
C****
      TITLE(9:50)=" Mass Stream Function (Sv)"
      DO KB=1,4
        TITLE(1:8)=TRIM(BASIN(KB))
        lname(1:31)=TITLE(1:31)
        sname='sf_'//trim(BASIN(KB))
        units='Sv'
        XJL(2:JM,1:LMO+1)=SFM(1:JM-1,0:LMO,KB)
        where(xjl.ne.undef) xjl=xjl/idacc(1)
        LP1=LMO+1
        IF (QDIAG) CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LP1,XJL
     *       ,ZOC1,"Latitude","Depth (m)")
C****
C**** Write data to PRinT file
C****
c        WRITE (6,921) 'Southern Hemisphere',(NINT(FLAT(J,2)),J=1,JEQ)
c        DO L=0,LMO
c          WRITE (6,922) ZOC1(L+1),(NINT(SFM(J,L,KB)),J=1,JEQ)
c        END DO
c        WRITE (6,921) 'Northern Hemisphere',(NINT(FLAT(J,2)),J=JEQ,
c    *     JM-1)
c        DO L=0,LMO
c          WRITE (6,922) ZOC1(L+1),(NINT(SFM(J,L,KB)),J=JEQ,JM-1)
c        END DO
      END DO

C**** Output Key diagnostics
      do L=1,LMO-1
        do J=2,JM-1
C**** North Atl. + North Pac. overturning
          if (oLAT_DG(j+1,2).gt.48-0.5*oDLAT_DG .and. oLAT_DG(j+1,2)
     *         .lt. 48+0.5*oDLAT_DG .and. 0.5*(ZE(L)+ZE(L-1)).le.900 
     *         .and. 0.5*(ZE(L)+ZE(L+1)).gt.900) then 
             WRITE(6,'(A46,F6.2)') 
     *            " North Atlantic overturning: 900m 48N: ",
     &           SFM(j,l,1)/idacc(1)
             WRITE(6,'(A46,F6.2)') 
     *            " North Pacific overturning:  900m 48N: ",
     &            SFM(j,l,2)/idacc(1)
          end if
C**** AABW production
          if (oLAT_DG(j+1,2).ge.-52-0.5*oDLAT_DG .and. oLAT_DG(j+1,2)
     *         .lt. -52+0.5*oDLAT_DG .and. 0.5*(ZE(L)+ZE(L-1)).le.3000
     *         .and. 0.5*(ZE(L)+ZE(L+1)).gt.3000) then
             WRITE(6,'(A46,F6.2)') " Antarctic Bottom Water production:"
     *            //" 3000m 52S: ",SFM(j,l,4)/idacc(1)
          end if
        end do
      end do

C****
C**** Calculate Salt Stream Function and write it
C****
c      FAC   = -1d-6/(IDACC(1)*NDYNO*DTO)
c      FACST = -.5E-6/(IDACC(1)*NDYNO*DTO)
c     CALL STRMJL(OIJL(1,1,1,IJL_SFLX),FAC,OLNST(1,1,LN_SFLX),FACST,SFS)
cC**** Subtract .035 times the Mass Stream Function
c      TITLE(9:50)=" Salt Stream Function (Sv)"
c      DO KB=1,4
c        DO L=0,LMO-1
c        DO J=1,JM-1
c          IF (SFS(J,L,KB).ne.UNDEF) SFS(J,L,KB) = SFS(J,L,KB) -
c     *         35.*SFM(J,L,KB)
c        END DO
c        END DO
c        TITLE(1:8)=TRIM(BASIN(KB))
c        lname(1:31)=TITLE(1:31)
c        sname='salt_sf_'//trim(BASIN(KB))
c        units='Sv'
c        XJL(2:JM,1:LMO+1)=SFS(1:JM-1,0:LMO,KB)
c        IF (QDIAG) CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO+1,XJL
c     *       ,ZOC1,"Latitude","Depth (m)")
c      END DO
C****
      RETURN
  920 FORMAT ('1',A80/)
  921 FORMAT ('0',A20 / '0',6X,23I5 / 8X,23('-----'))
  922 FORMAT (F6.0,2X,23I5)

      END SUBROUTINE OSFOUT

      SUBROUTINE STROUT
!@sum  STROUT prints out strait diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,amon,modelEclock
      USE OCEAN, only : im,jm,lmo,dts,ze
      USE STRAITS, only : nmst,wist,dist,lmst,name_st
      USE ODIAG
      IMPLICIT NONE
      REAL*8, DIMENSION(LMO,NMST) :: AS
      INTEGER I,J,K,L,NS,N,KB,IP1
     *     ,LMSTMIN(LMO),SUMORMN(KOLNST)
      REAL*8 MQ
      CHARACTER :: TITLE*80
C****
C**** Strait diagnostics
C****
      TITLE=""
      TITLE(51:80)=XLB
      SUMORMN(:)=1
      SUMORMN(LN_KVM)=2
      SUMORMN(LN_ICFL)=0
      DO N=1,KOLNST
        if (n.eq.LN_MFLX .or. n.eq.LN_GFLX.or. n.eq.LN_SFLX.or. n.eq
     *       .LN_KVM.or. n.eq.LN_ICFL) THEN
          AS = 0.
          TITLE(1:40) = TRIM(LNAME_OLNST(n))//' ('//
     &         TRIM(UNITS_OLNST(N))//')'
          DO NS=1,NMST
            DO L=1,LMST(NS)
              AS(L,NS) = OLNST(L,NS,N)*SCALE_OLNST(N)/IDACC(IA_OLNST(N))
            END DO
          END DO
          IF (N.eq.LN_ICFL) THEN
c            LMSTMIN(:)=1
c            CALL STABLE (LMSTMIN,name_st,AS,TITLE,SUMORMN(N))
          ELSE
            CALL STABLE (LMST,name_st,AS,TITLE,SUMORMN(N))
          END IF
        END IF
      END DO
C****
      RETURN
      END SUBROUTINE STROUT

      SUBROUTINE STABLE (LMST,STRAIT,AX,TITLE,SUMORMN)
!@sum STABLE produces a layer by strait table on the line printer.
C****
C**** Input:
!@var LMST = number of layers for each strait
!@var STRAIT = names of straits
!@var AX = two dimensional input array
!@var TITLE = title of run
!@var SUMORMN flag for whether sum or mean or neither ar printed
C****
      USE OCEAN, only : lmo
      USE STRAITS, only : nmst
      IMPLICIT NONE
      INTEGER, PARAMETER :: LMAXST=6
      INTEGER, INTENT(IN), DIMENSION(NMST) :: LMST
      INTEGER, INTENT(IN) :: SUMORMN
      REAL*8, DIMENSION(LMO,NMST), INTENT(IN) :: AX
      CHARACTER*20, INTENT(IN), DIMENSION(NMST) :: STRAIT
      CHARACTER*80, INTENT(IN) :: TITLE
      INTEGER, SAVE :: LINECT = 0.
      REAL*8 SUML
      INTEGER NS,L,LMAX
C****
C**** Produce leading title lines
C****
      LINECT = LINECT + 19
      IF(LINECT.LE.63)  THEN
        WRITE (6,900) TITLE
      ELSE
        LINECT = 8
        WRITE (6,901) TITLE
      ENDIF
      WRITE (6,902) ('------',L=1,LMAXST),(L,L=1,LMAXST),
     *              ('------',L=1,LMAXST)
C****
C**** Calculate sums and print the table
C****
      DO NS=1,NMST
        LMAX = LMST(NS)
        SUML = 0.
        DO L=1,LMAX
          SUML = SUML + AX(L,NS)
        END DO
        SELECT CASE (SUMORMN)
        CASE (0)                ! no sum or mean
          WRITE (6,'(1X,A20,9X,10I6)') STRAIT(NS),NINT(AX(1:LMAX,NS))
        CASE (1)
          WRITE (6,'(1X,A20,F8.1,1X,10I6)') STRAIT(NS),SUML
     *         ,NINT(AX(1:LMAX,NS))
        CASE (2)
          SUML=SUML/(LMAX-1)
          WRITE (6,'(1X,A20,F8.1,1X,10I6)') STRAIT(NS),SUML
     *         ,NINT(AX(1:LMAX,NS))
        END SELECT
      END DO
C****
      RETURN
  900 FORMAT ('0'//1X,A80)
  901 FORMAT ('1',A80)
  902 FORMAT ('0',29('-'),6A6 / ' Strait',15X,'Sum/Mean',6I6 /
     *        ' ',29('-'),6A6)
  904 FORMAT (1X,A20,F8.1,1X,10I6)
      END SUBROUTINE STABLE

#endif /* not NEW_IO */

      SUBROUTINE STRMIJ (MFU,FAC,OLNST,FACST,SF)
!@sum STRMIJ calculates the latitude by longitude stream function for a
!@+   given quantity.
C****
C**** Input:
!@var MFU  = west-east and south-north tracer fluxes (kg/s)
!@var   FAC = global scaling factor
!@var OLNST = strait mass flux (kg/s)
!@var FACST = global scaling factor for straits
C**** Output:
!@var    SF = stream function (kg/s)
C****
      Use OCEAN, Only: IM,JM,LMO, oDLAT_DG,oFJEQ=>FJEQ
      USE STRAITS, only : nmst,lmst
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FAC,FACST
      REAL*8, INTENT(IN), DIMENSION(IM,JM) :: MFU
      REAL*8, INTENT(IN), DIMENSION(LMO,NMST) :: OLNST
      REAL*8, INTENT(OUT), DIMENSION(IM,JM) :: SF
      Integer*4 :: I,J,L, NSUM
      REAL*8 TSUM
C****
C**** Integrate up from South Pole
C****
      SF=0
      DO J=2,JM-1
        DO I=1,IM
          SF(I,J) = SF(I,J-1) + MFU(I,J)*FAC
        END DO
C****
C**** Add strait flow from to the Stream Function
C****
        CALL STRMIJ_STRAITS(J,SF,OLNST,FACST)
      EndDo  !  End of Do J=2,JM-1

C****
C**** Recalibrate SF to be 0 over middle of North America
C**** Include UO cells whose centers reside in 110-90 W, 32-40 N
C****
      NSUM = 0
      TSUM = 0
      Do J = Ceiling(oFJEQ+32/oDLAT_DG), Floor(oFJEQ+40/oDLAT_DG)
      Do I = Ceiling(70*IM/360.), Floor(90*IM/360.)
        TSUM = TSUM + SF(I,J)
        NSUM = NSUM + 1  ;  EndDo  ;  EndDo
      SF(:,:) = SF(:,:) - TSUM/NSUM
      RETURN
      END

#ifndef NEW_IO
      SUBROUTINE OJLOUT
!@sum OJLOUT calculates basin means, lat. and long. sections
!@+   and advective tracer fluxes
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : idacc
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : n_waterm tracerlist, ocn_tracer_entry
#endif
      USE OCEAN, only : im,jm,lmo,ze,imaxj,focean,dypo,dts,dxvo
     *     ,dxypo, oDLAT_DG, oDLON_DG
      USE DIAG_COM, only : qdiag
      USE MDIAG_COM, only :
     &     sname_strlen,units_strlen,lname_strlen
      USE ODIAG
      IMPLICIT NONE
      REAL*8, DIMENSION(JM+3,LMO+1) :: XJL
      REAL*8 XB0(JM,LMO,NBAS),X0(IM,LMO+1),XS(IM,LMO+1),
     *     XBG(JM,LMO,NBAS),XBS(JM,LMO,NBAS),XG(IM,LMO+1)
#ifdef TRACERS_OCEAN
      REAL*8 XBT(JM,LMO,NBAS,tracerlist%getsize()),
     &   XT(IM,LMO+1,tracerlist%getsize()),XBTW(JM,LMO,NBAS)
#endif
      CHARACTER TITLE*80,EW*1,NS*1
      CHARACTER(len=lname_strlen) :: lname
      CHARACTER(len=sname_strlen) :: sname
      CHARACTER(len=units_strlen) :: units
      CHARACTER LABI*16,LABJ*16
      character*50 :: unit_string
      INTEGER I,J,L,KB,ISEC,II,ILON,JLAT,N,I1
      REAL*8 GOS,SOS,TEMGS,ASUM(IM),GSUM,ZONAL(LMO)
      type(ocn_tracer_entry), pointer :: entry

      IF (QDIAG) then
      XJL=0.

      XB0 = OJL(:,:,:,JL_M)
      XBG = OJL(:,:,:,JL_PT)
      XBS = OJL(:,:,:,JL_S)

      TITLE(51:80) = XLB
C****
C**** Quantites that are calculated
C****
C**** Basin and Zonal averages of tracers
C****
#ifdef TRACERS_OCEAN
      XBT = 0.
      DO L=1,LMO
        DO J= 1,JM
          DO I=1,IMAXJ(J)
            KB=KBASIN(I,J)
            IF(FOCEAN(I,J).gt..5)  THEN
              XBT(J,L,KB,:) = XBT(J,L,KB,:) + TOIJL(I,J,L,TOIJL_CONC,:)
            END IF
          END DO
        END DO
      END DO
      XBT(:,:,4,:)= SUM(XBT(:,:,1:4,:),DIM=3)
      if (n_water.ne.0) XBTW(:,:,:) =  XBT(:,:,:,n_water)
#endif

      DO KB=1,4
        DO J=1,JM
          DO L=1,LMO
            IF (XB0(J,L,KB).ne.0) THEN
#ifdef TRACERS_OCEAN
              do n=1,tracerlist%getsize()
              entry=>tracerlist%at(n)
              if (entry%to_per_mil.gt.0) then
                if (XBTW(j,l,kb).gt.0) then
                  XBT(j,l,kb,n)= 1d3*(XBT(j,l,kb,n)/
     *                 (XBTW(j,l,kb)*entry%trw0)-1.)
                else
                  XBT(j,l,kb,n)=undef
                end if
c                XBT(j,l,kb,n)= 1d3*(XBT(j,l,kb,n)/
c     *               ((XB0(J,L,KB)-XBS(J,L,KB))*entry%trw0)-1.)
              else
                XBT(j,l,kb,n)= 10.**(-entry%ntrocn)*XBT(j,l,kb,n)/
     *               (XB0(J,L,KB))
              end if
              end do
#endif
              XBG(J,L,KB) = XBG(J,L,KB)/XB0(J,L,KB)
              XBS(J,L,KB) = XBS(J,L,KB)/XB0(J,L,KB)
            ELSE
              XBG(J,L,KB)=undef
              XBS(J,L,KB)=undef
#ifdef TRACERS_OCEAN
              XBT(j,l,kb,:)= undef
#endif
            END IF
          END DO
        END DO
        if (KB.le.3) then
          TITLE(21:50)=" "//TRIM(BASIN(KB))//" Basin Average"
          SNAME(5:30)="_"//BASIN(KB)(1:3)
        else
          TITLE(21:50)=" Zonal Average"
          SNAME(5:30)="_Zonal"
        end if
        LNAME(21:50)=TITLE(21:50)
        TITLE(1:20)="Temperature (C)"
        LNAME(1:20)="Temperature"
        SNAME(1:4)="Temp"
        UNITS="C"
        XJL(1:JM,1:LMO)=XBG(1:JM,1:LMO,KB)
        CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *       ,"Latitude","Depth (m)")
        TITLE(1:20)="Salinity  (ppt)"
        LNAME(1:20)="Salinity"
        SNAME(1:4)="Salt"
        UNITS="ppt"
        XJL(1:JM,1:LMO)=XBS(1:JM,1:LMO,KB)
        CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *       ,"Latitude","Depth (m)")
#ifdef TRACERS_OCEAN
        DO N=1,tracerlist%getsize()
          entry=>tracerlist%at(n)
          if (entry%to_per_mil.gt.0) then
            UNITS="permil"
          else
            UNITS=unit_string(entry%ntrocn,'kg/kg')
          end if
          TITLE(1:20)=trim(entry%trname)//" ("//trim(UNITS)//")"
          LNAME(1:20)=trim(entry%trname)
          SNAME(1:4)=entry%trname(1:4)
          XJL(1:JM,1:LMO)=XBT(1:JM,1:LMO,KB,N)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
        END DO
#endif
      END DO
C****
C**** Longitudinal sections
C****
      DO ISEC=1,NSEC
        XB0 = undef  ; XBS = undef ;  XBG = undef
#ifdef TRACERS_OCEAN
        XBT = undef
#endif
        IF (SEC_LON(ISEC).ne.-1) THEN
          I=nint( .5*im + SEC_LON(ISEC)/odlon_dg)
          IF (I.le.IM/2) THEN
            EW="W"
          ELSE
            EW="E"
          END IF
          DO L=1,LMO
            DO J= 1,JM
              IF (OIJL(I,J,L,IJL_MO).ne.0) THEN
                GOS = OIJL(I,J,L,IJL_G0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                SOS = OIJL(I,J,L,IJL_S0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                XBG(J,L,1)= TEMGS(GOS,SOS)
                XBS(J,L,1)= 1d3*SOS
#ifdef TRACERS_OCEAN
                do n=1,tracerlist%getsize()
                  entry=>tracerlist%at(n)
                  if (entry%to_per_mil.gt.0) then
                    XBT(j,l,1,n)=1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
     *              (TOIJL(I,J,L,TOIJL_CONC,n_water)*entry%trw0)-1.)
c                    XBT(j,l,1,n)=1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
c     *              ((OIJL(I,J,L,IJL_MO)*DXYPO(J)-OIJL(I,J,L,IJL_S0M))
c     *                   *entry%trw0)-1.)
                  else
                    XBT(j,l,1,n)= 10.**(-entry%ntrocn)*TOIJL(I,J,L
     *                   ,TOIJL_CONC,n)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                  end if
                end do
#endif
              ELSE
                XBG(J,L,1)=undef ; XBS(J,L,1)=undef
#ifdef TRACERS_OCEAN
                XBT(J,L,1,:)=undef
#endif
              END IF
            END DO
          END DO
          ILON=ABS(SEC_LON(ISEC)-odlon_dg*0.5)
          WRITE(LABI,'(I3,A1)') ILON,EW
          TITLE(1:50)="Temperature Section         (C)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="Temp_"//adjustl(labi)
          UNITS="C"
          XJL(1:JM,1:LMO)=XBG(1:JM,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          TITLE(1:50)="Salinity Section            (ppt)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="Salt_"//adjustl(labi)
          UNITS="ppt"
          XJL(1:JM,1:LMO)=XBS(1:JM,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
#ifdef TRACERS_OCEAN
          DO N=1,tracerlist%getsize()
          entry=>tracerlist%at(n)
          if (entry%to_per_mil.gt.0) then
            UNITS="permil"
          else
            UNITS=unit_string(entry%ntrocn,'kg/kg')
          end if
          TITLE(1:50)=entry%trname
     &              //" Section          ("//TRIM(UNITS)//")"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME=trim(entry%trname)//"_"//adjustl(labi)
          XJL(1:JM,1:LMO)=XBT(1:JM,1:LMO,1,N)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          END DO
#endif

C**** EW fluxes of mass, heat and salt
C**** Note that ocean fluxes are defined 1:JM-1 and need to be moved up
C**** one place for compatibility with AGCM output format routines
          DO L=1,LMO
            DO J= 1,JM-1
C**** GM fluxes are also saved, so add the GM heat and salt fluxes here
              IF (OIJL(I,J,L,IJL_MFU).ne.0) THEN
                XB0(J,L,1)= 1d0*OIJL(I,J,L,IJL_MFU)/(IDACC(1)*DTS
     *               *DYPO(J) *(ZE(L)-ZE(L-1)))
                XBG(J,L,1)= 1d-6*(OIJL(I,J,L,IJL_GFLX)+OIJL(I,J,L
     *               ,IJL_GGMFL))/(IDACC(1)*DTS*DYPO(J)*(ZE(L)-ZE(L-1)))
                XBS(J,L,1)= 1d-6*(OIJL(I,J,L,IJL_SFLX)+OIJL(I,J,L
     *               ,IJL_SGMFL))/(IDACC(1)*DTS*DYPO(J)*(ZE(L)-ZE(L-1)))
              ELSE
                XB0(J,L,1)=undef
                XBG(J,L,1)=undef
                XBS(J,L,1)=undef
              END IF
            END DO
          END DO
          ILON=ABS(SEC_LON(ISEC))
          WRITE(LABI,'(I3,A1)') ILON,EW
          TITLE(1:50)=" EW Mass Flux            (kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="EWmflx_"//adjustl(labi)
          UNITS="kg/s m^2"
          XJL(2:JM,1:LMO)=XB0(1:JM-1,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          TITLE(1:50)=" EW Heat Flux            (10^6 J/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="EWhflx_"//adjustl(labi)
          UNITS="10^6 J/s m^2"
          XJL(2:JM,1:LMO)=XBG(1:JM-1,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          TITLE(1:50)=" EW Salt Flux            (10^6 kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="EWsflx_"//adjustl(labi)
          UNITS="10^6 kg/s m^2"
          XJL(2:JM,1:LMO)=XBS(1:JM-1,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
        END IF
      END DO
C****
C**** Latitudinal sections
C****
C**** Define starting point for wrap (~20W) (to avoid splitting pacific)
      I1 = 1+200./odlon_dg

      ASUM=0. ; GSUM=0. ; ZONAL=0.
      DO ISEC=1,NSEC
        X0=undef ; XS=undef ; XG=undef
#ifdef TRACERS_OCEAN
        XT=undef
#endif
        IF (SEC_LAT(ISEC).ne.-1) THEN
          IF (SEC_LAT(ISEC).le.0) THEN
            NS="S"
          ELSE
            NS="N"
          END IF
          J=nint(.5*jm+SEC_LAT(ISEC)/odlat_dg )
          DO L=1,LMO
            DO II= 1,IM
              I = II+I1-1
              IF (I.gt.IM) I=I-IM
              IF (OIJL(I,J,L,IJL_MO).ne.0) THEN
                GOS = OIJL(I,J,L,IJL_G0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                SOS = OIJL(I,J,L,IJL_S0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                XG(II,L)= TEMGS(GOS,SOS)
                XS(II,L)= 1d3*SOS
#ifdef TRACERS_OCEAN
                do n=1,tracerlist%getsize()
                  entry=>tracerlist%at(n)
                  if (entry%to_per_mil.gt.0.and.TOIJL(I,J,L,TOIJL_CONC
     *                 ,n_water).gt.0) then
                    XT(II,l,n)= 1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
     *              (TOIJL(I,J,L,TOIJL_CONC,n_water)*entry%trw0)-1.)
c                    XT(II,l,n)= 1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
c     *              ((OIJL(I,J,L,IJL_MO)*DXYPO(J)-OIJL(I,J,L,IJL_S0M))
c     *                   *entry%trw0)-1.)
                  else
                    XT(II,l,n)= 10.**(-entry%ntrocn)*TOIJL(I,J,L,
     &                   TOIJL_CONC,n)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                  end if
                end do
#endif
              ELSE
                XG(II,L)=undef ; XS(II,L)=undef
#ifdef TRACERS_OCEAN
                XT(II,L,:)=undef
#endif
              END IF
            END DO
          END DO
          JLAT=ABS(SEC_LAT(ISEC)-0.5*odlat_dg)
          WRITE(LABJ,'(I3,A1)') JLAT,NS
          TITLE(1:50)="Temperature Section              (C)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="Temp_"//adjustl(labj)
          UNITS="C"
          CALL POUT_IL(TITLE,sname,lname,units,I1,1,LMO,XG
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          TITLE(1:50)="Salinity Section                 (ppt)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="Salt_"//adjustl(labj)
          UNITS="ppt"
          CALL POUT_IL(TITLE,sname,lname,units,I1,1,LMO,XS
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
#ifdef TRACERS_OCEAN
          DO N=1,tracerlist%getsize()
          entry=>tracerlist%at(n)
          if (entry%to_per_mil.gt.0) then
            UNITS="permil"
          else
            UNITS=unit_string(entry%ntrocn,'kg/kg')
          end if
          TITLE(1:50)=entry%trname
     &         //" Section          ("//TRIM(UNITS)//")"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME=trim(entry%trname)//"_"//adjustl(labj)
          CALL POUT_IL(TITLE,sname,lname,units,I1,1,LMO,XT(1,1,N)
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          END DO
#endif
C**** Fluxes
          DO L=1,LMO
            DO II= 1,IM
              I = II+I1-1
              IF (I.gt.IM) I=I-IM
              IF (OIJL(I,J,L,IJL_MFV).ne.0) THEN
                X0(II,L)= OIJL(I,J,L,IJL_MFV)/(IDACC(1)*DTS
     *               *DXVO(J) *(ZE(L)-ZE(L-1)))
                XG(II,L)= 1d-6*(OIJL(I,J,L,IJL_GFLX+1)+OIJL(I,J,L
     *               ,IJL_GGMFL+1))/(IDACC(1)*DTS*DXVO(J)*(ZE(L)-ZE(L-1)
     *               ))
                XS(II,L)= 1d-6*(OIJL(I,J,L,IJL_SFLX+1)+OIJL(I,J,L
     *               ,IJL_SGMFL+1))/(IDACC(1)*DTS*DXVO(J)*(ZE(L)-ZE(L-1)
     *               ))
              ELSE
                X0(II,L)=undef ; XG(II,L)=undef ; XS(II,L)=undef
              END IF
            END DO
          END DO
          JLAT=ABS(SEC_LAT(ISEC))
          WRITE(LABJ,'(I3,A1)') JLAT,NS
          TITLE(1:50)=" NS Mass Flux            (kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="NSmflx_"//adjustl(labj)
          UNITS="kg/s m^2"
          CALL POUT_IL(TITLE,sname,lname,units,I1,2,LMO,X0
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          TITLE(1:50)=" NS Heat Flux            (10^6  J/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="NShflx_"//adjustl(labj)
          UNITS="10^6 W/m^2"
          CALL POUT_IL(TITLE,sname,lname,units,I1,2,LMO,XG
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          TITLE(1:50)=" NS Salt Flux            (10^6 kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="NSsflx_"//adjustl(labj)
          UNITS="10^6 kg/s m^2"
          CALL POUT_IL(TITLE,sname,lname,units,I1,2,LMO,XS
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
        END IF
      END DO
      END IF
C****
      RETURN
      END SUBROUTINE OJLOUT
#endif

      subroutine oijl_prep
c
c Convert oijl accumulations into the desired units
c
      use ocean, only : im,jm,lmo,lmm,imaxj,focean,dxypo,dxvo,dypo,dts
      use ocean, only : nbyzm,i1yzm,i2yzm
      use ocnmeso_com, only : use_tdmix
      use odiag, only : koijl,oijl_out,oijl=>oijl_loc,ijl_area
     &     ,igrid_oijl,jgrid_oijl,sname_oijl
     &     ,ijl_mo,ijl_mou,ijl_mov,ijl_g0m,ijl_s0m,ijl_ptm,ijl_pdm
     &     ,ijl_mfu,ijl_mfv,ijl_mfw,ijl_mfw2,ijl_ggmfl,ijl_sgmfl
#ifdef TDMIX_AUX_DIAGS
     &     ,ijl_gsymmf,ijl_ssymmf
#endif
     &     ,ijl_wgfl,ijl_wsfl,ijl_kvm,ijl_kvg,ijl_kvx,ijl_gflx,ijl_sflx
     &     ,ijl_mfub,ijl_mfvb,ijl_mfwb,ijl_isdm,ijl_pdm2
     &     ,oij=>oij_loc,ij_sf,olnst,ln_mflx
#ifdef OCN_GISS_TURB
     &     ,ijl_ri,ijl_rrho,ijl_bv2,ijl_otke,ijl_kvs,ijl_kvc,ijl_buoy
#endif
#ifdef OCN_GISS_SM
     &     ,ijl_fvb
#endif
#ifdef OCEAN_TENDENCY_DIAGS
     &     ,olnst,ln_gflx,ln_sflx
#endif
      use odiag, only : ia_oijl
#ifdef TRACERS_OCEAN
      use odiag, only :
     &     ktoijlx,toijl_out,divbya_toijl,kn_toijl,toijl_loc,toijl_conc
     &    ,toijl_tflx,toijl_gmfl
      USE OCN_TRACER_COM, only : n_Water, tracerlist, ocn_tracer_entry
#endif
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : am_i_root,halo_update,south
     &     ,hassouthpole,hasnorthpole
     &     ,pack_data,unpack_data ! for horz stream function
      use mdiag_com, only : ia_cpl
      use model_com, only : idacc
      use constant, only : grav
      use kpp_com, only : use_tdiss
#ifdef OCEAN_TENDENCY_DIAGS
      use straits, only: ist,jst,lmst
      use domain_decomp_1d, only : broadcast
#endif
      use straits, only : nmst
      implicit none
      integer i,j,l,k,kk,n
      real*8 mass,gos,sos,temgs,volgs,volgsp,fac,facst,dpr
      real*8 wdenom,xedge,wtdn,wtup,massdn,massup
      integer :: j_0,j_1,j_0s,j_1s
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     mfu,pres
      real*8, dimension(:,:), allocatable :: mfu_glob,sf_glob
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     mfub,mfvb
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,0:lmo) ::
     &     mfwb,dmfwb
      real*8 :: sncor
      integer :: ib
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: entry
#endif
#ifdef OCEAN_TENDENCY_DIAGS
      integer, parameter :: max_num_tends=16 ! maximum number of tendency outputs
      integer :: flxind(max_num_tends)
      character(len=32) :: tendname(max_num_tends)
      real*8 :: sign_end
#endif
      real*8 dumarr(1,1)

      j_0 = grid%j_strt
      j_1 = grid%j_stop
      j_0s = grid%j_strt_skp
      j_1s = grid%j_stop_skp

      oijl_out(:,:,:,:) = 0.

c
c Cell-centered quantities. Some conversions to per square meter
c
      pres = 0. ! need to add atm + sea ice press, but this is a small
                ! effect for typical uses of in-situ density diag
      do l=1,lmo
      do j=j_0,j_1
      do i=1,imaxj(j)
        mass = oijl(i,j,l,ijl_mo)*dxypo(j)
        if(focean(i,j).le..5 .or. mass.le.0.) cycle
        oijl_out(i,j,l,ijl_mo) = mass
        oijl_out(i,j,l,ijl_g0m) = oijl(i,j,l,ijl_g0m)
        oijl_out(i,j,l,ijl_s0m) = oijl(i,j,l,ijl_s0m)
c
c compute potential temperature, potential density, in-situ density
c
        dpr = grav*oijl(i,j,l,ijl_mo)/idacc(ia_oijl(ijl_mo))
        pres(i,j) = pres(i,j) + .5d0*dpr
        gos = oijl(i,j,l,ijl_g0m) / mass
        sos = oijl(i,j,l,ijl_s0m) / mass
        oijl_out(i,j,l,ijl_ptm) = mass*temgs(gos,sos)
        oijl_out(i,j,l,ijl_pdm) = mass*(1d0/volgs(gos,sos)-1000d0) !sigma0
        oijl_out(i,j,l,ijl_pdm2)= mass*(1d0/volgsp(gos,sos,2d7)-1d3) !sigma2
        oijl_out(i,j,l,ijl_isdm) =
     &       mass*(1d0/volgsp(gos,sos,pres(i,j))-1000d0)
        pres(i,j) = pres(i,j) + .5d0*dpr
      enddo
      enddo
      enddo


      dmfwb = 0.
      mfub = oijl(:,:,:,ijl_mfub)
      mfvb = oijl(:,:,:,ijl_mfvb)
      mfwb = 0.
      call halo_update(grid,mfvb,from=south)
      if(use_tdmix==1) then
        ! Calculate dmfwb, the bolus-induced component of the mass
        ! flux used in the remapping (vertical advective) tracer flux.
        ! It is used below for post-hoc repartitioning of accumulated
        ! advective fluxes into resolved- and bolus-velocity components.
        ! See notes in TDMIX.
        mfwb(:,:,1:lmo) = oijl(:,:,1:lmo,ijl_mfwb)
        do l=1,lmo-1
        do j=j_0s,j_1s
          i=1
          if(l.lt.lmm(i,j)) then
            dmfwb(i,j,l) = dmfwb(i,j,l-1) + (
     &         (mfub(im,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &        +(mfwb(i,j,l-1)-mfwb(i,j,l))
     &         )!/dxypo(j)
          endif
          do n=1,nbyzm(j,l+1)
          do i=max(2,i1yzm(n,j,l+1)),i2yzm(n,j,l+1)
            dmfwb(i,j,l) = dmfwb(i,j,l-1) + (
     &         (mfub(i-1,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &        +(mfwb(i,j,l-1)-mfwb(i,j,l))
     &         )!/dxypo(j)
          enddo
          enddo
        enddo
        enddo
      else
c****
c**** derive bolus vertical mass flux from bolus horizontal mass fluxes
c**** for skew-GM
        do l=1,lmo-1
        do j=j_0s,j_1s
          i=1
          if(l.lt.lmm(i,j)) then
            mfwb(i,j,l) = mfwb(i,j,l-1) + (
     &         (mfub(im,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &         )
          endif
          do n=1,nbyzm(j,l+1)
          do i=max(2,i1yzm(n,j,l+1)),i2yzm(n,j,l+1)
            mfwb(i,j,l) = mfwb(i,j,l-1) + (
     &         (mfub(i-1,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &         )
          enddo
          enddo
        enddo
        enddo
        oijl(:,:,:,ijl_mfwb) = mfwb(:,:,1:lmo)
      endif
c
c Vertical fluxes.  Some conversions to per square meter
c
      do l=1,lmo-1
      do j=j_0,j_1
      do n=1,nbyzm(j,l+1)
      do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
        oijl_out(i,j,l,ijl_area) = idacc(ia_cpl)*dxypo(j)
          
        oijl_out(i,j,l,ijl_mfw) = oijl(i,j,l,ijl_mfw)
     &       -dmfwb(i,j,l) ! subtract bolus-induced part of remap flux
        oijl_out(i,j,l,ijl_mfwb) = oijl(i,j,l,ijl_mfwb)
     &       +dmfwb(i,j,l) !      add bolus-induced part of remap flux


        massup = oijl(i,j,l  ,ijl_mo)*dxypo(j)
        massdn = oijl(i,j,l+1,ijl_mo)*dxypo(j)
        wtdn = massup/(massup+massdn)
        wtup = 1d0-wtdn

        ! centered approximation of layer edge g for repartitioning
        xedge =
     &       wtup* (oijl(i,j,l  ,ijl_g0m) / massup)
     &      +wtdn* (oijl(i,j,l+1,ijl_g0m) / massdn)
        oijl_out(i,j,l,ijl_gflx+2) = oijl(i,j,l,ijl_gflx+2)
     &       -dmfwb(i,j,l)*xedge ! subtract bolus-induced part of remap flux
        oijl_out(i,j,l,ijl_ggmfl+2) = oijl(i,j,l,ijl_ggmfl+2)
     &       +dmfwb(i,j,l)*xedge !      add bolus-induced part of remap flux

        ! centered approximation of layer edge s for repartitioning
        xedge =
     &       wtup* (oijl(i,j,l  ,ijl_s0m) / massup)
     &      +wtdn* (oijl(i,j,l+1,ijl_s0m) / massdn)
        oijl_out(i,j,l,ijl_sflx+2) = oijl(i,j,l,ijl_sflx+2)
     &       -dmfwb(i,j,l)*xedge ! subtract bolus-induced part of remap flux
        oijl_out(i,j,l,ijl_sgmfl+2) = oijl(i,j,l,ijl_sgmfl+2)
     &       +dmfwb(i,j,l)*xedge !      add bolus-induced part of remap flux

        oijl_out(i,j,l,ijl_mfw2) = oijl(i,j,l,ijl_mfw2)/dxypo(j)
#ifdef TDMIX_AUX_DIAGS
        oijl_out(i,j,l,ijl_gsymmf+2) = oijl(i,j,l,ijl_gsymmf+2)
        oijl_out(i,j,l,ijl_ssymmf+2) = oijl(i,j,l,ijl_ssymmf+2)
#endif

        oijl_out(i,j,l,ijl_wgfl) = oijl(i,j,l,ijl_wgfl)
        oijl_out(i,j,l,ijl_wsfl) = oijl(i,j,l,ijl_wsfl)
        oijl_out(i,j,l,ijl_kvm) = oijl(i,j,l,ijl_kvm)*dxypo(j)
        oijl_out(i,j,l,ijl_kvg) = oijl(i,j,l,ijl_kvg)*dxypo(j)
        if(use_tdiss==1) then
          oijl_out(i,j,l,ijl_kvx) = oijl(i,j,l,ijl_kvx)*dxypo(j)
        endif
#ifdef OCN_GISS_TURB
        oijl_out(i,j,l,ijl_kvs) = oijl(i,j,l,ijl_kvs)*dxypo(j)
        oijl_out(i,j,l,ijl_kvc) = oijl(i,j,l,ijl_kvc)*dxypo(j)
        oijl_out(i,j,l,ijl_ri) = oijl(i,j,l,ijl_ri)*dxypo(j)
        oijl_out(i,j,l,ijl_rrho) = oijl(i,j,l,ijl_rrho)*dxypo(j)
        oijl_out(i,j,l,ijl_bv2) = oijl(i,j,l,ijl_bv2)*dxypo(j)
        oijl_out(i,j,l,ijl_buoy) = oijl(i,j,l,ijl_buoy)*dxypo(j)
        oijl_out(i,j,l,ijl_otke) = oijl(i,j,l,ijl_otke)*dxypo(j)
#endif
#ifdef OCN_GISS_SM
        oijl_out(i,j,l,ijl_fvb) = oijl(i,j,l,ijl_fvb)*dxypo(j)
#endif
      enddo
      enddo
      enddo
      enddo

c
c Horizontal fluxes.  Some conversions to per meter
c
      do k=1,koijl
        call halo_update(grid,oijl(:,:,:,k),from=south)
      enddo
      do l=1,lmo
        do j=j_0s,j_1s
        do i=1,im
          oijl_out(i,j,l,ijl_mfu) = oijl(i,j,l,ijl_mfu)
          oijl_out(i,j,l,ijl_mfub) = oijl(i,j,l,ijl_mfub)
          oijl_out(i,j,l,ijl_gflx) = oijl(i,j,l,ijl_gflx)
          oijl_out(i,j,l,ijl_sflx) = oijl(i,j,l,ijl_sflx)
          oijl_out(i,j,l,ijl_ggmfl) = oijl(i,j,l,ijl_ggmfl)
          oijl_out(i,j,l,ijl_sgmfl) = oijl(i,j,l,ijl_sgmfl)
#ifdef TDMIX_AUX_DIAGS
          oijl_out(i,j,l,ijl_gsymmf) = oijl(i,j,l,ijl_gsymmf)
          oijl_out(i,j,l,ijl_ssymmf) = oijl(i,j,l,ijl_ssymmf)
#endif
        enddo
        do i=1,im-1
          oijl_out(i,j,l,ijl_mou) =
     &         .5*(oijl(i,j,l,ijl_mo)+oijl(i+1,j,l,ijl_mo))*dypo(j)
        enddo
        i=im
          oijl_out(i,j,l,ijl_mou) =
     &         .5*(oijl(i,j,l,ijl_mo)+oijl(1,j,l,ijl_mo))*dypo(j)
        enddo ! j
        do j=max(2,j_0),j_1
        do i=1,im
          oijl_out(i,j,l,ijl_mfv) = oijl(i,j-1,l,ijl_mfv)
          oijl_out(i,j,l,ijl_mfvb) = oijl(i,j-1,l,ijl_mfvb)
          oijl_out(i,j,l,ijl_mov) =
     &         .5*(oijl(i,j,l,ijl_mo)+oijl(i,j-1,l,ijl_mo))*dxvo(j-1)
          oijl_out(i,j,l,ijl_gflx+1) = oijl(i,j-1,l,ijl_gflx+1)
          oijl_out(i,j,l,ijl_sflx+1) = oijl(i,j-1,l,ijl_sflx+1)
          oijl_out(i,j,l,ijl_ggmfl+1) = oijl(i,j-1,l,ijl_ggmfl+1)
          oijl_out(i,j,l,ijl_sgmfl+1) = oijl(i,j-1,l,ijl_sgmfl+1)
#ifdef TDMIX_AUX_DIAGS
          oijl_out(i,j,l,ijl_gsymmf+1) = oijl(i,j-1,l,ijl_gsymmf+1)
          oijl_out(i,j,l,ijl_ssymmf+1) = oijl(i,j-1,l,ijl_ssymmf+1)
#endif
        enddo
        enddo ! j
      enddo

#ifdef OCEAN_TENDENCY_DIAGS
      ! compute convergences of 3D resolved/parameterized fluxes
      mfwb(:,:,0) = 0.

      k = 0
      k = k+1; flxind(k) = ijl_gflx; tendname(k) = 'g_advtend'
      k = k+1; flxind(k) = ijl_sflx; tendname(k) = 's_advtend'
      k = k+1; flxind(k) = ijl_ggmfl; tendname(k) = 'g_mesotend'
      k = k+1; flxind(k) = ijl_sgmfl; tendname(k) = 's_mesotend'
#ifdef TDMIX_AUX_DIAGS
      k = k+1; flxind(k) = ijl_gsymmf; tendname(k) = 'g_mesotend_sym'
      k = k+1; flxind(k) = ijl_ssymmf; tendname(k) = 's_mesotend_sym'
#endif
      kk = k
      ! currently, convergences are in assumed positions in oijl (xflux+3)
      do k=1,kk ! so sanity-check the names first
        if(trim(sname_oijl(flxind(k)+3)).ne.trim(tendname(k))) then
          call stop_model(
     &         'oijl_prep: OCEAN_TENDENCY_DIAGS name mismatch',255)
        endif
      enddo
      do k=1,kk
        mfwb(:,:,1:lmo) = oijl_out(:,:,1:lmo,flxind(k)+2)
        call do_fluxconv_3d(
     &       oijl(:,:,:,flxind(k)+0),
     &       oijl(:,:,:,flxind(k)+1),
     &       mfwb,
     &       oijl_out(:,:,:,flxind(k)+3) )
      enddo

      ! Deal with the straits
      call broadcast(grid, olnst)
      do n=1,nmst
        sign_end = -.5 ! .5 is the scale factor for olnst
        do k=1,2
          i=ist(n,k)
          j=jst(n,k)
          if(j.ge.j_0 .and. j.le.j_1) then
            do l=1,lmst(n)
              oijl_out(i,j,l,ijl_gflx+3) = oijl_out(i,j,l,ijl_gflx+3)
     &             + olnst(l,n,ln_gflx)*sign_end
              oijl_out(i,j,l,ijl_sflx+3) = oijl_out(i,j,l,ijl_sflx+3)
     &             + olnst(l,n,ln_sflx)*sign_end
            enddo
          endif
          sign_end = -sign_end
        enddo
      enddo
#endif

      ! Fill poles
      do k=1,koijl
        if(jgrid_oijl(k).ne.1 .or. igrid_oijl(k).ne.1) cycle
        do l=1,lmo
          !if(hassouthpole(grid)) then
          !  j = j_0
          !  oijl_out(2:im,j,l,k) = oijl_out(1,j,l,k)
          !endif
          if(hasnorthpole(grid)) then
            j = j_1
            oijl_out(2:im,j,l,k) = oijl_out(1,j,l,k)
          endif
        enddo
      enddo

C****
C**** Calculate Horizontal Mass Stream Function
C****
      do j=j_0s,j_1s
      do i=1,im
        mfu(i,j) = 0.
        do l=1,lmm(i,j)
          mfu(i,j) = mfu(i,j) + oijl(i,j,l,ijl_mfu)
        enddo
      enddo
      enddo
      if(am_i_root()) then
        allocate(mfu_glob(im,jm),sf_glob(im,jm))
      else
        allocate(mfu_glob(1,1),sf_glob(1,1))
      endif
      call pack_data(grid,mfu,mfu_glob)
      if(am_i_root()) then
        FAC   = -1d-9/dts
        FACST = -1d-9/dts
        if(nmst.gt.0) then
          CALL STRMIJ(MFU_GLOB,FAC,OLNST(1,1,LN_MFLX),FACST,SF_GLOB)
        else
          CALL STRMIJ(MFU_GLOB,FAC,DUMARR,FACST,SF_GLOB)
        endif
      endif
      call unpack_data(grid,sf_glob,oij(:,:,ij_sf))
      deallocate(mfu_glob,sf_glob)

#ifdef TRACERS_OCEAN
C****
C**** Tracers
C****
      toijl_out(:,:,:,:) = 0.
      kk = 1
      toijl_out(:,:,:,kk) = oijl_out(:,:,:,ijl_mo)
      do kk=2,ktoijlx
        k = kn_toijl(1,kk)
        n = kn_toijl(2,kk)
        entry=>tracerlist%at(n)
        if(k.le.0 .or. n.le.0) cycle
        toijl_out(:,:,:,kk) = toijl_loc(:,:,:,k,n)
        if(divbya_toijl(kk)) then
          do l=1,lmo; do j=j_0,j_1
            toijl_out(:,j,l,kk) = toijl_out(:,j,l,kk)/dxypo(j)
          enddo; enddo
        endif
        if(entry%to_per_mil>0 .and. n.ne.n_Water) then
          toijl_out(:,:,:,kk) = 1d3*(toijl_out(:,:,:,kk)/entry%trw0
     &         -toijl_loc(:,:,:,TOIJL_conc,n_water))
        endif
        if(use_tdmix==1 .and.
     &       (k.eq.toijl_gmfl+2 .or. k.eq.toijl_tflx+2)) then
          if(k.eq.toijl_gmfl+2) then
            sncor = +1d0
          else
            sncor = -1d0
          endif
          do l=1,lmo-1
          do j=j_0,j_1
          do ib=1,nbyzm(j,l+1)
          do i=i1yzm(ib,j,l+1),i2yzm(ib,j,l+1)
            massup = oijl(i,j,l  ,ijl_mo)*dxypo(j)
            massdn = oijl(i,j,l+1,ijl_mo)*dxypo(j)
            wtdn = massup/(massup+massdn)
            wtup = 1d0-wtdn
            ! centered approximation of layer edge tracer for repartitioning
            xedge =
     &           wtup* (toijl_loc(i,j,l  ,toijl_conc,n) / massup)
     &          +wtdn* (toijl_loc(i,j,l+1,toijl_conc,n) / massdn)
            toijl_out(i,j,l,kk) = toijl_out(i,j,l,kk) 
                 ! add or subtract bolus-induced part of remap flux
     &           +(dmfwb(i,j,l)/dxypo(j))*xedge*sncor
          enddo
          enddo
          enddo
          enddo
        endif
      enddo
#endif

      return

      contains

      subroutine do_fluxconv_3d(mfub,mfvb,mfwb,conv_out)
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     mfub,mfvb,conv_out
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,0:lmo) ::
     &     mfwb

      real*8 :: byim

      !call halo_update(grid,mfvb,from=south)
      do l=1,lmo
        do j=grid%j_strt_skp,grid%j_stop_skp
          i=1
          if(l.le.lmm(i,j)) then
            conv_out(i,j,l) =
     &         (mfub(im,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &        +(mfwb(i,j,l-1)-mfwb(i,j,l))
          endif
          do n=1,nbyzm(j,l)
          do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
            conv_out(i,j,l) =
     &         (mfub(i-1,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &        +(mfwb(i,j,l-1)-mfwb(i,j,l))
          enddo
          enddo
        enddo
      enddo
      if(hasnorthpole(grid)) then
        j = grid%j_stop
        byim = 1d0/real(im,kind=8)
        do l=1,lmm(1,j)
          conv_out(:,j,l) = sum(mfvb(:,j-1,l))*byim
     &         +(mfwb(1,j,l-1)-mfwb(1,j,l))
        enddo
      endif
      end subroutine do_fluxconv_3d

      end subroutine oijl_prep

#ifndef NEW_IO
      SUBROUTINE OTJOUT
!@sum OTJOUT print vertically integrated basin and zonal
!@+   northward transports
!@auth Gavin Schmidt/Gary Russell
      USE MODEL_COM, only : idacc
      USE OCEAN, only : im,jm,lmo
      USE DIAG_COM, only : qdiag
      USE ODIAG
      IMPLICIT NONE
      CHARACTER TITLE*80, YAXIS(3)*8
      REAL*8 VLAT(0:JM)
      INTEGER, PARAMETER :: INC=1+(JM-1)/24
      CHARACTER*50, DIMENSION(11) ::  OTJNAME=(/
     *     'Northward Transport of Mass (10^9 kg/sec)        ',
     *     'North. Trans. of Potential Enthalpy (10^15 W)    ',
     *     'North. Trans. of Salt - .035*Mass (10^6 kg/sec)  ',
     *     'North. Trans. of Heat in Atl. Ocean (10^15 W)    ',
     *     'North. Trans. of Heat in Pac. Ocean (10^15 W)    ',
     *     'North. Trans. of Heat in Indian Ocean (10^15 W)  ',
     *     'North. Trans. of Heat in Global Ocean (10^15 W)  ',
     *     'North. Trans. of Salt in Atl. Ocean (10^6 kg/s)  ',
     *     'North. Trans. of Salt in Pac. Ocean (10^6 kg/s)  ',
     *     'North. Trans. of Salt in Indian Ocean (10^6 kg/s)',
     *     'North. Trans. of Salt in Global Ocean (10^6 kg/s)'/)
      DATA YAXIS /'  Mass  ', 'Enthalpy', '  Salt  '/
      INTEGER J,KQ,KB,IDLAT
      REAL*8 FJEQ
      REAL*8 X(0:JM,4,3),XCOMP(0:JM,4,3,3)

      X = OTJ/idacc(1)
      XCOMP = OTJCOMP/idacc(1)

      IDLAT= NINT(180./(JM-1))
      FJEQ = JM/2.
      DO J=1,JM-1
        VLAT( J) = IDLAT*(J-FJEQ)
      END DO
      VLAT( 0) = -90.
      VLAT(JM) =  90.
C****
C**** Write titles and data to disk.
C****
      TITLE(51:80)=XLB
      DO KQ=1,3
        TITLE(1:50)=OTJNAME(KQ)
        WRITE (6,907) TITLE(1:72)
C**** print out truncated series to PRT file
        WRITE (6,903) (NINT(VLAT(J)),J=JM,0,-INC)
        DO KB=1,4
          WRITE (6,906) BASIN(KB),(X(J,KB,KQ),J=JM,0,-INC)
        END DO
        WRITE(6,904)
        IF (QDIAG) THEN
          WRITE (iu_otj,*) TITLE
          WRITE (iu_otj,*) 'Latitude'
          WRITE (iu_otj,*) YAXIS(KQ)
          WRITE (iu_otj,*)
     *         ' Lat   Atlantic   Pacific    Indian     Global'
          WRITE (iu_otj,971) (VLAT(J),(X(J,KB,KQ),KB=1,4),J=0,JM)
          WRITE (iu_otj,*)
        END IF
      END DO
C****
C**** Write titles and data to disk for northward transports
C**** by components
C****
      DO KQ=2,3
        WRITE(6,908)
        DO KB=1,4
          TITLE(1:50)=OTJNAME(3+KB+4*(KQ-2))
          WRITE (6,907) TITLE(1:72)
C**** print out truncated series to PRT file
          WRITE(6,903) (NINT(VLAT(J)),J=JM,0,-INC)
          WRITE(6,906)BASIN(KB)(1:3)//" Overtrn",
     &         (XCOMP(J,KB,1,KQ),J=JM,0,-INC)
          WRITE(6,906)BASIN(KB)(1:3)//" GM flx ",
     &         (XCOMP(J,KB,2,KQ),J=JM,0,-INC)
          WRITE(6,906)BASIN(KB)(1:3)//" Hor gyr",
     &         (XCOMP(J,KB,3,KQ),J=JM,0,-INC)
          WRITE(6,906)BASIN(KB)(1:3)//" Total  ",
     &         (X(J,KB,KQ),J=JM,0,-INC)
          WRITE(6,904)
          IF (QDIAG) THEN
            WRITE (iu_otj,*) TITLE
            WRITE (iu_otj,*) 'Latitude'
            WRITE (iu_otj,*) YAXIS(2)
            WRITE (iu_otj,*)
     *           ' Lat      Total    Overturn    GM_flux    Hor_Gyre'
            WRITE (iu_otj,976)(VLAT(J),X(J,KB,KQ),
     &           XCOMP(J,KB,1:3,KQ),J=0,JM)
            WRITE (iu_otj,*)
          END IF
        END DO
      END DO
      RETURN
C****
  970 FORMAT (A50,A6,2X,2A4)
  971 FORMAT (F6.0,4F10.3)
  972 FORMAT (A,1P,3E15.4)
  973 FORMAT (A,2I6)
  976 FORMAT (F6.0,4F10.3)
  903 FORMAT (' ',131('-')/,' Latitude  ',24I5)
  904 FORMAT (' ',131('-'))
  906 FORMAT (' ',A11,24F5.1)
  907 FORMAT ('0',A)
  908 FORMAT ('1')

      END SUBROUTINE OTJOUT
#endif

      Subroutine STRMIJ_STRAITS (J,SF,OLNST,FACST)
C****
C**** Add strait flow to IxJ Strean Function
C****
C**** Input:   J = ocean model latitude index
C****      OLNST = ocean strait mass flux (kg/s)
C****      FACST = global scaling factor for ocean straits
C**** Output: SF = IxJ stream function (kg/s)
C****
      Use OCEAN,   Only: IM,JM,LMO
      Use STRAITS, Only: NMST,LMST, IST,JST
      Implicit None

      Integer*4,Intent(In) :: J
      Real*8,Intent(InOut) :: SF(IM,JM)
      Real*8,Intent(In)    :: OLNST(LMO,NMST),FACST

C**** Local variables
      Integer*4 N,LM, I1,J1, I2,J2

      Do 40 N=1,NMST
      I1 = IST(N,1)       ;  I2 = IST(N,2)
      J1 = JST(N,1)       ;  J2 = JST(N,2)
      LM = LMST(N)
      If (J2 - J1) 10,20,30

C**** JST(N,2) < JST(N,1)
   10 If (J==J1 .or. J==J2)  Then
         SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST*.5/(J1-J2)
         GoTo 40  ;  EndIf
      If (J2 < J .and. J < J1)
     *   SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST/(J1-J2)
      GoTo 40

C**** JST(N,2) = JST(N,1)
   20 If (J==J1)
     *   SF(I1:I2-1,J) = SF(I1:I2-1,J) + Sum(OLNST(1:LM,N))*FACST
      GoTo 40

C**** JST(N,2) > JST(N,1)
   30 If (J==J1 .or. J==J2)  Then
         SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST*.5/(J2-J1)
         GoTo 40  ;  EndIf
      If (J1 < J .and. J < J2)
     *   SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST/(J2-J1)
   40 Continue
      Return
      EndSubroutine STRMIJ_STRAITS
