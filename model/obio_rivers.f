(from obio_com.f)
#ifdef OBIO_RUNOFF

#ifdef NITR_RUNOFF
!      real, ALLOCATABLE, DIMENSION(:,:)    :: rnitrmflo_loc      ! riverine nitrate mass
flow rate (kg/s)
      real, ALLOCATABLE, DIMENSION(:,:)    :: rnitrconc_loc      ! riverine nitrate
concentration (kg/kg)
#endif
#ifdef DIC_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rdicconc_loc       ! riverine dic
concentration (kg/kg)
#endif
#ifdef DOC_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rdocconc_loc       ! riverine doc
concentration (kg/kg)
#endif 
#ifdef SILI_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rsiliconc_loc      ! riverine silica
concentration (kg/kg)
#endif
#ifdef IRON_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rironconc_loc      ! riverine iron
concentration (kg/kg)
#endif
#ifdef POC_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rpocconc_loc       ! riverine poc
concentration (kg/kg)
#endif
#ifdef ALK_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: ralkconc_loc       ! riverine alkalinity
concentration (mol/kg)
#endif

#endif
#ifdef OBIO_RUNOFF

#ifdef NITR_RUNOFF
      real rnitrconc_ij
!     .    , rnitrmflo_ij
#endif
#ifdef DIC_RUNOFF
      real rdicconc_ij
#endif
#ifdef DOC_RUNOFF
      real rdocconc_ij
#endif
#ifdef SILI_RUNOFF
      real rsiliconc_ij
#endif
#ifdef IRON_RUNOFF
      real rironconc_ij
#endif
#ifdef POC_RUNOFF
      real rpocconc_ij
#endif
#ifdef ALK_RUNOFF
      real ralkconc_ij
#endif
#endif

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!      ALLOCATE(rnitrmflo_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(rnitrconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef DIC_RUNOFF
      ALLOCATE(rdicconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef DOC_RUNOFF
      ALLOCATE(rdocconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef SILI_RUNOFF
      ALLOCATE(rsiliconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef IRON_RUNOFF
      ALLOCATE(rironconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef POC_RUNOFF
      ALLOCATE(rpocconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef ALK_RUNOFF
      ALLOCATE(ralkconc_loc(i_0:i_1,j_0:j_1))
#endif
#endif

--------------------------------------------------------------------------------
(from obio_incom.f)

#ifdef OBIO_RUNOFF
#ifdef IRON_RUNOFF
      real, parameter ::  estFe = 0.01   ! estuarine retention rate, can vary between 0.2
and 0.01 (daCunha 2007)
#endif
#endif

--------------------------------------------------------------------------------
(from obio_init.f)
#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!     .                    ,rnitrmflo_loc
     .                    ,rnitrconc_loc
#endif
#ifdef DIC_RUNOFF
     .                    ,rdicconc_loc
#endif
#ifdef DOC_RUNOFF
     .                    ,rdocconc_loc
#endif
#ifdef SILI_RUNOFF
     .                    ,rsiliconc_loc
#endif
#ifdef IRON_RUNOFF
     .                    ,rironconc_loc
#endif
#ifdef POC_RUNOFF
     .                    ,rpocconc_loc
#endif
#ifdef ALK_RUNOFF
     .                    ,ralkconc_loc
#endif
#endif

#ifdef OBIO_RUNOFF
! read in nutrient concentrations, already regridded to model grid
        if (AM_I_ROOT()) then
        print*, '    '
        print*, 'reading nutrient runoff data.....'
        print*, '    '
        endif
#ifdef NITR_RUNOFF
!        filename='rnitr_mflo'
        filename='rnitr_conc'
        fid=par_open(ogrid,filename,'read')
!       call read_dist_data(ogrid,fid,'din',rnitrmflo_loc)
        call read_dist_data(ogrid,fid,'din',rnitrconc_loc)
        call par_close(ogrid,fid)
#endif
#ifdef DIC_RUNOFF
        filename='rdic_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'dic',rdicconc_loc)
        call par_close(ogrid,fid)
        write(*,*)'reading dic from',filename
#endif
#ifdef DOC_RUNOFF
        filename='rdoc_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'doc',rdocconc_loc)
        call par_close(ogrid,fid)
#endif
#ifdef SILI_RUNOFF
        filename='rsili_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'sil',rsiliconc_loc)
        call par_close(ogrid,fid)
#endif
#ifdef IRON_RUNOFF
        filename='riron_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'fe',rironconc_loc)
        call par_close(ogrid,fid)
#endif
#ifdef POC_RUNOFF
        filename='rpoc_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'poc',rpocconc_loc)
        call par_close(ogrid,fid)
#endif
#ifdef ALK_RUNOFF
        filename='ralk_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'alk',ralkconc_loc)
        call par_close(ogrid,fid)
#endif
#endif

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!      call add_diag("Nitrate mass flow from rivers", "oij_rnitrmflo",
!     &               "kg/s", IJ_rnitrmflo)
      call add_diag("Nitrate conc in runoff", "oij_rnitrconc",
     &              "kg/kg", .false., IJ_rnitrconc)
#endif
#ifdef DIC_RUNOFF
      call add_diag("DIC conc in runoff", "oij_rdicconc",
     &              "kg/kg", .false., IJ_rdicconc)
#endif
#ifdef DOC_RUNOFF
      call add_diag("DOC conc in runoff", "oij_rdocconc",
     &              "kg/kg", .false., IJ_rdocconc)
#endif
#ifdef SILI_RUNOFF
      call add_diag("silica conc in runoff", "oij_rsiliconc",
     &              "kg/kg", .false., IJ_rsiliconc)
#endif
#ifdef IRON_RUNOFF
      call add_diag("iron conc in runoff", "oij_rironconc",
     &              "kg/kg", .false., IJ_rironconc)
#endif
#ifdef POC_RUNOFF
      call add_diag("poc conc in runoff", "oij_rpocconc",
     &              "kg/kg", .false., IJ_rpocconc)
#endif
#ifdef ALK_RUNOFF
      call add_diag("alkalinity conc in runoff", "oij_ralkconc",
     &              "mol/kg", .false., IJ_ralkconc)
#endif
#endif

--------------------------------------------------------------------------------
(from obio_model)
#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!      use obio_com, only: rnitrmflo_loc
      use obio_com, only: rnitrconc_loc
#endif
#ifdef DIC_RUNOFF
      use obio_com, only: rdicconc_loc
#endif
#ifdef DOC_RUNOFF
      use obio_com, only: rdocconc_loc
#endif
#ifdef SILI_RUNOFF
      use obio_com, only: rsiliconc_loc
#endif
#ifdef IRON_RUNOFF
      use obio_com, only: rironconc_loc
#endif
#ifdef POC_RUNOFF
      use obio_com, only: rpocconc_loc
#endif
#ifdef ALK_RUNOFF
      use obio_com, only: ralkconc_loc
#endif
#endif

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
      USE obio_diag, only: ij_rnitrconc
!     .                 ,ij_rnitrmflo
#endif
#ifdef DIC_RUNOFF
      USE obio_diag, only: ij_rdicconc
#endif
#ifdef DOC_RUNOFF
      USE obio_diag, only:  ij_rdocconc
#endif
#ifdef SILI_RUNOFF
      USE obio_diag, only:  ij_rsiliconc
#endif
#ifdef IRON_RUNOFF
      USE obio_diag, only:  ij_rironconc
#endif
#ifdef POC_RUNOFF
      USE obio_diag, only:  ij_rpocconc
#endif
#ifdef ALK_RUNOFF
      USE obio_diag, only:  ij_ralkconc
#endif
#endif


#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!       OIJ(I,J,IJ_rnitrmflo) = OIJ(I,J,IJ_rnitrmflo)+rnitrmflo_loc(i,j)  ! riverine nitr
mass flow from dC (kg/s)
       OIJ(I,J,IJ_rnitrconc) = OIJ(I,J,IJ_rnitrconc)+rnitrconc_loc(i,j)  ! riverine nitr
conc from dC (kg/kg)
#endif
#ifdef DIC_RUNOFF
       OIJ(I,J,IJ_rdicconc) = OIJ(I,J,IJ_rdicconc)+rdicconc_loc(i,j)     ! riverine dic conc
from dC (kg/kg)
#endif
#ifdef DOC_RUNOFF
       OIJ(I,J,IJ_rdocconc) = OIJ(I,J,IJ_rdocconc)+rdocconc_loc(i,j)     ! riverine doc conc
from dC (kg/kg)
#endif
#ifdef SILI_RUNOFF
       OIJ(I,J,IJ_rsiliconc) = OIJ(I,J,IJ_rsiliconc)+rsiliconc_loc(i,j)  ! riverine silica
conc from dC (kg/kg)
#endif
#ifdef IRON_RUNOFF
       OIJ(I,J,IJ_rironconc) = OIJ(I,J,IJ_rironconc)+rironconc_loc(i,j)  ! riverine iron
conc from dC (kg/kg)
#endif
#ifdef POC_RUNOFF
       OIJ(I,J,IJ_rpocconc) = OIJ(I,J,IJ_rpocconc)+rpocconc_loc(i,j)     ! riverine poc conc
from dC (kg/kg)
#endif
#ifdef ALK_RUNOFF
       OIJ(I,J,IJ_ralkconc) = OIJ(I,J,IJ_ralkconc)+ralkconc_loc(i,j)     ! riverine
alkalinity conc from A-S (mol/kg)
#endif
#endif




--------------------------------------------------------------------------------
(from ptend)

#ifdef OBIO_RUNOFF
#ifdef POC_RUNOFF
     .                      ,rpocconc_loc
#endif
#ifdef NITR_RUNOFF
     .                      ,rnitrconc_loc
!                         ,rnitrfmlo_loc
#endif
#ifdef SILI_RUNOFF
     .                      ,rsiliconc_loc
#endif
#ifdef IRON_RUNOFF
     .                      ,rironconc_loc
      USE obio_incom, only: estFe
#endif
#endif

#ifdef OBIO_RUNOFF
      USE OFLUXES, only:  oFLOWO
!      USE ocean, only:  oxyp
      USE MODEL_COM, only:  dtsrc
#endif


#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!       if (oFLOWO(i,j) .gt. 0.)then
!        rnitr_loc = rnitrmflo_loc
!     .    / (oFLOWO(i,j) * dxypo(j)) ! kg/s => kg,N/kg,w/s
!     .    * 1.d3     ! kg,N to g,N
!     .    * (1./14.) ! g,N to mol,N
!     .    * 1.d3     ! mol,N to mmol,N
!     .    * rho_water ! kg,water to m3 water
!         if (i.eq.169.and. j.eq.59) then
!           write(*,'(/,a,2i5,6e12.4)')'i,j,rnitrmflo, rnitr, oFLOWO,
!     .       dxypo, nitr, rho_water:',i,j,rnitrmflo_loc(i,j),
!     .       rnitr_loc(i,j),
!     .       oFLOWO(i,j),dxypo(j),obio_P(1,1),rho_water
!       else
!         rnitr_loc = 0.
!       endif
         term = rnitrconc_loc(i,j)
     .    * oFLOWO(i,j)/dtsrc         ! kg,N/kg,w => kg,N/m2,w/s
     .    / dp1d(1)                   ! kg,N/m2,w/s => kg,N/m3,w/s
     .    * 1.d6/14.                  ! kg,N/m3,w/s => mmol,N/m3,w/s
     .    * 3600.                     ! mmol,N/m3,w/s => mmol,N/m3,w/hr
         rhs(1,1,17) = term
         P_tend(1,1) = P_tend(1,1) + term

!       if (i.eq.169 .and. j.eq.59) then
!         write(*,'(/,a,2i5,5e12.4)')'i,j,rnitrconc,rnitr,oFLOWO,
!     .       dtsrc,dp1d(1):',i,j,rnitrconc_loc(i,j),rhs(1,1,17),
!     .       oFLOWO(i,j),dtsrc,dp1d(1) 
!                endif

#endif
#ifdef SILI_RUNOFF
        term = rsiliconc_loc(i,j)
     .   * oFLOWO(i,j)/dtsrc         ! kg,S/kg,w => kg,S/m2,w/s
     .   / dp1d(1)                   ! kg,S/m2,w/s => kg,S/m3,w/s
     .   * 1.d6/28.055               ! kg,S/m3,w/s => mmol,S/m3,w/s
     .   * 3600.                     ! mmol,S/m3,w/s => mmol,S/m3,w/hr
        rhs(1,3,17) = term
        P_tend(1,3) = P_tend(1,3) + term
#endif
#ifdef IRON_RUNOFF
        term = rironconc_loc(i,j)
     .   * oFLOWO(i,j)/dtsrc         ! kg,Fe/kg,w => kg,Fe/m2,w/s
     .   / dp1d(1)                   ! kg,Fe/m2,w/s => kg,Fe/m3,w/s
     .   * 1.d9/55.845              ! kg,Fe/m3,w/s => umol,Fe/m3,w/s
     .   * 3600.                     ! umol,Fe/m3,w/s => umol,Fe/m3,w/hr
     .   * estFe                     ! estuarine retention rate
        rhs(1,4,17) = term
        P_tend(1,4) = P_tend(1,4) + term
#endif
#ifdef POC_RUNOFF
        term = rpocconc_loc(i,j)
     .   * oFLOWO(i,j)/dtsrc         ! kg,C/kg,w => kg,C/m2,w/s
     .   / dp1d(1)                   ! kg,C/m2,w/s => kg,C/m3,w/s
     .   * 1.d6                      ! kg,C/m3,w/s => mg,C/m3,w/s
     .   * 3600.                     ! mg,C/m3,w/s => mg,C/m3,w/hr
        rhs(1,10,17) = term
        D_tend(1,1) = D_tend(1,1) + term
#endif
#endif


--------------------------------------------------------------------------------
(from obio_carbon.f)
#ifdef OBIO_RUNOFF
#ifdef DOC_RUNOFF
     .                    ,rdocconc_loc
#endif
#ifdef DIC_RUNOFF
     .                    ,rdicconc_loc
#endif
#endif

#ifdef OBIO_RUNOFF
      USE OFLUXES, only:  oFLOWO
#endif

fdef OBIO_RUNOFF
#ifdef DOC_RUNOFF
        term = rdocconc_loc(i,j)
     .    * oFLOWO(i,j)/dtsrc          ! kg,C/kg,w => kg,C/m2,w/s
     .    / dp1d(1)                    ! kg,C/m2,w/s => kg,C/m3,w/s
     .    * 1.d6/12.                   ! kg,C/m3,w/s => mmol,C/m3,w/s
     .    * 3600.                      ! mmol,C/m3,w/s => mmol,C/m3,w/hr
        rhs(1,13,17) = term
        C_tend(1,1) = C_tend(1,1) + term
#endif
#ifdef DIC_RUNOFF
        term = rdicconc_loc(i,j)
     .    * oFLOWO(i,j)/dtsrc          ! kg,C/kg,w => kg,C/m2,w/s
     .    / dp1d(1)                    ! kg,C/m2,w/s => kg,C/m3,w/s
     .    * 1.d6/12.                   ! kg,C/m3,w/s => mmol,C/m3,w/s
     .    *3600.                       ! mmol,C/m3,w/s => mmol,C/m3,w/hr
        rhs(1,14,17) = term
        C_tend(1,2) = C_tend(1,2) + term


        if (i.eq.169.and.j.eq.59) then
        write(*,'(/,a,2i5,2e12.4)')'i,j,rdicconc,term:',
     .          i,j,rdicconc_loc(i,j),term
        endif
#endif
#endif

--------------------------------------------------------------------------------
(from obio_alkalinity.f)
#ifdef OBIO_RUNOFF
#ifdef ALK_RUNOFF
      use obio_com, only: ralkconc_loc
      USE OFLUXES, only:  oFLOWO
      USE MODEL_COM, only: dtsrc
#endif
#endif

#ifdef OBIO_RUNOFF
#ifdef ALK_RUNOFF
      term = ralkconc_loc(i,j)
     . * oFLOWO(i,j)/dtsrc           ! mol/kg,w => mol/m2,w/s
     . / dp1d(1)                     ! mol/m2,w/s => mol/m3,w/s
     . * 1.d3                        ! mol/m3,w/s => mmol/m3,w/s
!    . * 3600.                       ! mmol/m3,w/s => mmol/m3/hr    commented out July 2016
      rhs(1,15,17) = term
      A_tend(1) = A_tend(1) + term
#endif
#endif 


