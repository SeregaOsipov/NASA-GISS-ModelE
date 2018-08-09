Emaobio_08112017.R GISS Model E  coupled version          larissa   04/15/2010

!! obio_cfc is for NIsurf=1 (U00a=0.72; U00b=1.68)

Emaobio_062717_2: E4F40 coupled to 1x1.25deg 32-layer GISS ocean model
based on E4F40 = modelE as frozen in April 2010:
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean: coupled to 1x1.25deg 32-layer GISS ocean model (Russell - Schmidt)
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define CHECK_OCEAN                 ! needed to compile aux/file CMPE002
#define NEW_IO
#define OCN_LAYERING L32
#define OBIO_ON_GISSocean           ! obio on GISS ocean
#define TRACERS_OCEAN               ! GISS Ocean tracers activated
#define TRACERS_OCEAN_INDEP         ! independently defined ocn tracers
#define TRACERS_OceanBiology
#define pCO2_ONLINE
#define TRACERS_GASEXCH_ocean       ! ANY ocean: special tracers to be passed to ocean
#define TRACERS_GASEXCH_ocean_CO2   ! ANY ocean: special tracers to be passed to ocean
#define OCN_CFC
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
Atm144x90                  ! horizontal resolution is 144x90 -> 2x2.5deg
AtmL40                      ! vertical resolution is 40 layers -> 0.1mb
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform
ORES_1Qx1 OFFT288E                  ! ocean horiz res 1.25x1deg

IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

    ! lat-lon grid specific source codes
AtmRes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
MODEL_COM                           ! calendar, timing variables
MODELE_DRV                          ! ModelE cap
MODELE                              ! initialization and main loop
ATM_COM                             ! main atmospheric variables
ATM_DRV                             ! driver for atmosphere-grid components
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
IRRIGMOD                            ! irrigation module
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV     ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO ocalbedo ! radiation and albedo
DIAG_COM DIAG DEFACC                ! diagnostics
OCN_DRV                             ! driver for ocean-grid components
OLAYERS                             ! ocean layering options
ODIAG_COM  OCEAN_COM  OGEOM         ! dynamic ocean modules
ODIAG_ZONAL
OCNDYN  OCNDYN2  OTIDELL            ! dynamic ocean routines
OCN_Interp                          ! dynamic ocean routines
OSTRAITS_COM  OSTRAITS              ! straits parameterization
OCNKPP                              ! ocean vertical mixing
OCNMESO_DRV OCNTDMIX OCNGM          ! ocean mesoscale mixing
OCEANR_DIM  OFLUXES
ODIAG_PRT                           ! ocean diagnostic print out
OCNFUNTAB                           ! ocean function look up table
SparseCommunicator_mod              ! sparse gather/scatter module
OCNQUS                              ! QUS advection scheme
OCNGISS_TURB                        ! GISS Turbulence vertical mixing scheme
OCNGISS_SM                          ! GISS Sub-mesoscale mixing scheme

! Codes common to atmospheric tracer sets
TRACER_COM                          ! configurable tracer code
TRACERS_DRV |-O0|                   ! O0 speeds compilation and no difference for this file
initTracerGriddedData		    ! grid independent info
initTracerMetadata		    ! misc initialization that used to be mixed in with metadata
TRACERS                             ! generic tracer code
TRDRYDEP                            ! dry deposition of tracers
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
lightning                           ! lightning module
MiscTracersMetadata
sharedTracersMetadata

OCN_Int_LATLON                      ! atm-ocn regrid routines

     ! atmospheric tracers
TRDIAG
TRACER_GASEXCH_CO2                  ! tracer functions needed for gas exch expts
TRACER_GASEXCH_CFC                 ! tracer functions needed for gas exch expts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!   OCEAN TRACERS       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OCN_TRACER_COM
OCN_TRACER


     ! ocean carbon cycle
obio_dim         |$(R8)|
obio_incom       |$(R8)|
obio_com_R       |$(R8)|
obio_forc        |$(R8)|
obio_init        |$(R8)|
obio_bioinit     |$(R8)|
obio_model       |$(R8)|
obio_daysetrad   |$(R8)|
obio_daysetbio   |$(R8)|
obio_sfcirr      |$(R8)|
obio_edeu        |$(R8)|
obio_ptend       |$(R8)|
obio_carbon      |$(R8)|
obio_update      |$(R8)|
obio_sinksettl_R |$(R8)|
!!!obio_alkalinity  |$(R8)|
!!!obio_diffmod|$(R8)|
obio_chkbalances |$(R8)|
obio_conserv_R |$(R8)|

Components:
tracers
shared MPI_Support solvers giss_LSM
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT 

Data input files:
    ! start from the restart file of an earlier run ...                 ISTART=8
! AIC=1....rsfE... ! initial conditions, no GIC needed, use
!! AIC=1JAN1961.rsfE4F40.MXL65m   ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
! AIC=AIC.RES_F40.D771201.nc      ! observed initial conditions for F40 1977/12/01
! AIC=AIC_144x90_DEC01_L96.nc     ! observed initial conditions for F96 1977/12/01
AIC=NCARIC.144x90.D7712010_ext.nc ! AIC for automatic relayering to model vertical grid
GIC=GIC.144X90.DEC01.1.ext_1.nc   ! initial ground conditions
OIC=OIC288X180.D1201.nc             ! Levitus ocean intial conditions
OFTAB=OFTABLE_NEW                   ! ocean function table
KBASIN=KB288X180.modelE.BS1.nc      ! ocean basin designations (1 cell Bering Strait)
TOPO_OC=OZ1QX1N.BS1.nc              ! ocean fraction and topography (1 cell Bering Strait)
TIDES=TIDAL_e_v2_1QX1               ! ocean bottom tidal energy and velocity squared
OSTRAITS=OSTRAITS_288x180.nml       ! parameterized straits info
TOPO=Z2HX2fromZ1QX1N.BS1.nc           ! surface fractions and topography

RVR=RD_modelE_Fa.nc ! river direction file
NAMERVR=RD_modelE_Fa.names.txt ! named river outlets

CDN=CD144X90.ext.nc
VEG=V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.ext.nc
LAIMAX=V144x90_EntMM16_lai_max_trimmed_scaled_ext.nc
HITEent=V144x90_EntMM16_height_trimmed_scaled_ext.nc
LAI=V144x90_EntMM16_lai_trimmed_scaled_ext.nc
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
IRRIG=Irrig144x90_1848to2100_FixedFuture_v3.nc
SOIL=S144X900098M.ext.nc
TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc             ! topographic variation for gwdrag
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
GLMELT=GLMELT_144X90_gas.OCN.nc
RADN1=sgpgxg.table8                           ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800  ! rad.tables and history files
RADN4=LWCorrTables33k                         ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-2012.May13_hdr
RADN8=cloud.epsilon4.72x46
!RADN9=solar.lean2015.ann1610-2014.nc ! need KSOLAR=2
RADN9=solar.CMIP6official.ann1850-2299.nc ! need KSOLAR=2
RADNE=topcld.trscat8

ISCCP=ISCCP.tautables
GHG=GHG.CMIP6.1-2014.txt  !  GreenHouse Gases for CMIP6 runs up to 2014
CO2profile=CO2profile.Jul2017.txt ! scaling of CO2 in stratosphere
! CO2profile=CO2profile.Jan2017.txt ! original 40-layer version
dH2O=dH2O_by_CH4_monthly

DUSTaer=TcadiAR5_aero/144x90/DUST_Tcadi2012_Bauer_kg_m2_144x90x40_unlim.nc
BC_dep=TcadiAR5_aero/144x90/BC_dep_Tcadi2012_Bauer_kg_m2_s_144x90_1850-2100.nc
! updated aerosols need MADAER=3
TAero_SUL=TcadiAR5_aero/144x90/SUL_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100.nc
TAero_SSA=TcadiAR5_aero/144x90/SSA_Tcadi2012_Bauer_kg_m2_144x90x40.nc
TAero_NIT=NIT_Bauer2008_kg_m2_144x90x40_1850-2100h.nc
TAero_OCA=TcadiAR5_aero/144x90/OCA_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100.nc
TAero_BCA=TcadiAR5_aero/144x90/BCA_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100.nc
TAero_BCB=TcadiAR5_aero/144x90/BCB_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100.nc


! O3file soon to be replaced by one from latest chemistry code
O3file=jan2012_o3_shindell_144x90x49x12_1850-2010_ple.nc
Ox_ref=o3_2010_shindell_144x90x49_April1850.nc

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

!!!!!!!!!!!!!!!!!!! obio  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cfle1=abw25b.dat                         ! seawater spectral absorp.
                                         ! and scatt. coefs
cfle2=acbc25b.dat                        ! phytoplankton spectrl absorp.
                                         ! and scatt. coefs
!!!!pco2table=pco2.tbl.asc               ! table to compute pco2 values
                                         ! from sst,sss,dic,alk
                                         ! if not defined pCO2_ONLINE
nitrates_inicond=no3_nodc_annmean_180x288.nc    ! initial cond for nitrates (NODC)
silicate_inicond=sio2_nodc_annmean_180x288.nc   ! initial cond for silicate (NODC)
dic_inicond=dic_glodap_annmean_180x288.nc       ! initial cond for dic (GLODAP)
alk_inicond=alk_glodap_annmean_180x288.nc       ! initial cond/forc for alk(GLODAP)
!!!oasimdirect=oasimdirect_20w_new       ! spectral light components
                                         ! if not def OBIO_RAD_coupling
atmFe_inicond=iron_gocart_1x1mon_180x288.nc     ! GOCART iron flux
atmFedirect1=iron_ron_195x180_20w.asc    ! Ron Miller's dust fluxes
facirr=facirr.asc                        ! factors for mean irrad w/in water
eda_esa_ratios=eda_esa_ratios.asc        ! ratios of rad spectrl components
!!!!!!!!!!!!!!!!!!! obio_rad  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHL_DATA=CHL_WG_2x2.5zavg                !CHL_WG_4x5 or CHL_WG_2x2.5
!! CHL_DATA=CHL_WG_2x2.5                    !CHL_WG_4x5 or CHL_WG_2x2.5
                                         !in RUSSELL'socean grid
                                         !to be used with CHL_from_SeaWIFs
cfcatm_data=CFC_ATM_Hist_2014.txt   !CFC concentrations in atmosphere

Label and Namelist:  (next 2 lines)
Emaobio_08112017 (NIsurf=1; 2x2.5x40, 1850 atm.;  1x1.25x32 ocean)


&&PARAMETERS
ocean_trname = 'OceanAge  Ventilatn WatrMass1  WatrMass2  WatrMass3 aoCFC'
! parameters set for coupled ocean runs:
KOCEAN=1            ! ocn is prognostic
OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)
OTIDE = 0           !  Ocean tides are not used
variable_lk=1
init_flake=1

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000055  ! threshold (1/s) for triggering deformation waves
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.1       ! default is 0.5
CDEF=1.6       ! tuning factor for deformation -> momentum flux
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF


! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.60 ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.00 ! below 850mb and MC regions; tune this last  to get rad.balance

WMUI_multiplier = 1.

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=1850
!crops_yr=1850  ! if -1, crops in VEG-file is used
!s0_yr=1850
!s0_day=182
!ghg_yr=1850
!ghg_day=182
!irrig_yr=1850
volc_yr=-1
!volc_day=182
!aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
!albsn_yr=1850
dalbsnX=.024
!o3_yr=-1850
!aer_int_yr=1850    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=0
!orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols

DTO=112.5
DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)

! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time
SUBDD=' '        ! no sub-daily frequency diags
NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
KCOPY=1          ! 0: no output; 1: save .acc; 2: unused; 3: include ocean data
KRSF=12          ! 0: no output; X: save rsf at the beginning of every X month
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=480
&&END_PARAMETERS

 &INPUTZ
!  YEARI=1900,MONTHI=12,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
!  YEARE=2001,MONTHE=12,DATEE=1,HOURE=0, KDIAG=13*0,
!  ISTART=2,IRANDI=0, YEARE=1900,MONTHE=12,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,

   YEARI=1850,MONTHI=1,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=1850,MONTHE=1,DATEE=2,HOURE=0, KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1850,MONTHE=1,DATEE=2,HOURE=1,IWRITE=1,JWRITE=1,
/
