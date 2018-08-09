obio_cfc.R GISS Model E  coupled version          larissa   04/15/2010

!! obio_cfc is for NIsurf=1 (U00a=0.72; U00b=1.68)

obio_cfc: E4F40 coupled to 1x1.25deg 32-layer GISS ocean model
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

#include "latlon_source_files"
#include "modelE4_source_files"
#include "dynamic_ocn_source_files"
#include "tracer_shared_source_files"
OCN_Int_LATLON                      ! atm-ocn regrid routines

#include "ocarbon_cycle_oR_files" ! both gas exch and ocean tracer oR model

Components:
tracers
#include "E4_components"    /* without "Ent" */
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT /* needed for "Ent" only */

Data input files:
#include "IC_144x90_input_files"
#include "dynamic_ocn_288x180_input_files_AR5"
TOPO=Z2HX2fromZ1QX1N.BS1.nc           ! surface fractions and topography

RVR=RD_modelE_Fa.nc ! river direction file
NAMERVR=RD_modelE_Fa.names.txt ! named river outlets

#include "land144x90_input_files"
#include "rad_input_files"
#include "rad_144x90_input_files"

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
                                         !in GISS ocean grid
                                         !to be used with CHL_from_SeaWIFs
cfcatm_data=CFC_ATM_Hist_2014.txt   !CFC concentrations in atmosphere

Label and Namelist:  (next 2 lines)
obio_cfc (NIsurf=1; 2x2.5x40, 1850 atm.;  1x1.25x32 ocean)


&&PARAMETERS
ocean_trname = 'OceanAge  Ventilatn WatrMass1  WatrMass2  WatrMass3 aoCFC'
#include "dynamic_ocn_params"

#include "sdragF40_params"
#include "gwdragF40_params"

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.60 ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.00 ! below 850mb and MC regions; tune this last  to get rad.balance

WMUI_multiplier = 1.

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

#include "atmCompos_1850_params"
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

#include "diag_params"

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=480
&&END_PARAMETERS

 &INPUTZ
!  YEARI=1900,MONTHI=12,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
!  YEARE=2001,MONTHE=12,DATEE=1,HOURE=0, KDIAG=13*0,
!  ISTART=2,IRANDI=0, YEARE=1900,MONTHE=12,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,

   YEARI=1931,MONTHI=1,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=1931,MONTHE=1,DATEE=2,HOURE=0, KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1931,MONTHE=1,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
/
