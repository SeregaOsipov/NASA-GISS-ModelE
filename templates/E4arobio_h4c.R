E4arobio_h4c.R GISS Model E  coupled version          tnl   05/17/2009

E4arobio_h4c: 2x2.5x40 layers modelE version, 1850 atm.; 
	   1x1x26 layers in the HYCOM ocean

modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)     
atmospheric composition from year 1850                               
ocean: coupled to HYCOM ocean model 
uses turbulence scheme (no dry conv),  grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs  
filters: U,V in E-W direction (after every dynamics time step)              
         sea level pressure (after every physics time step)                 

Preprocessor Options
#define NEW_IO
#define CHECK_OCEAN                 ! needed to compile aux/file CMPE002
#define ATM2x2h             !2x2.5 40 layer atm & 26 layer 1deg hycom (387x360)
#define HYCOM1deg           !2x2.5 40 layer atm & 26 layer 1deg hycom (387x360)
#define TRACERS_OceanBiology
#define pCO2_ONLINE
! #define constCO2
#define TRACERS_ON                  ! include tracers code
#define TRACERS_GASEXCH_ocean       ! ANY ocean: special tracers to be passed to ocean
#define TRACERS_GASEXCH_ocean_CO2   ! ANY ocean: special tracers to be passed to ocean
! #define TRACERS_HYCOM_Ventilation
End Preprocessor Options


Object modules:
     ! resolution-specific source codes
Atm144x90                  ! horizontal resolution is 144x90 -> 2x2.5deg
AtmL40                      ! vertical resolution is 40 layers -> 0.1mb
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform
IO_DRV                              ! new i/o
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV                             ! advection of Q/tracers
TQUS_DRV                            ! advection of Q
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

#include "latlon_source_files"
#include "modelE4_source_files"
#include "hycom_source_files"
#include "ocarbon_cycle_oH_files"   
#include "tracer_shared_source_files"

Components:
#include "E4_components_nc"    /* without "Ent" */
tracers Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB

Data input files:
AIC=AIC.RES_F40.D771201.nc      ! observed init cond (atm. only) ISTART=2
GIC=GIC.144X90.DEC01.1.ext_1.nc ! initial ground conditions      ISTART=2
TOPO=Z144X90N.1deghycom_1.nc
CDN=CD144X90.ext.nc              ! neutral drag coefficient
VEG=V144X90_no_crops.ext.nc      ! vegatation file 
CROPS=CROPS_144X90N_nocasp.ext.nc ! crops
SOIL=S144X900098M.ext.nc         ! soil properties
REG=REG2X2.5                     ! special regions-diag
RVR=RD_modelE_Fa_1deghycom.nc             ! river direction file
NAMERVR=RD_modelE_Fa_1deghycom.names.txt  ! named river outlets
#include "rad_input_files"

#include "rad_144x90_input_files"

TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc             ! topographic variation for gwdrag
MSU_wts=MSU_SSU_RSS_weights.txt
GLMELT=GLMELT_144X90_gas.OCN.nc   ! glacial melt distribution
latlonij=latlon387x360.4bin             ! lat & lon at each i,j
hycomtopo=depth387x360.4bin_1  ! topography used in ocean model, NO Baltic
temp_ini=temp387x360x26jan_hv_z1.txt ! 3-d temperature as initial condition
salt_ini=salt387x360x26jan_hv_z1.txt ! 3-d salinity as initial condition
pout_ini=pout387x360x26jan_hv_z1.txt ! 3-d layer pressure as initial condition
ibasin=ibasin387x360.txt_1  ! basin mask, Baltic
flxa2o=flxa2o387x360.8bin_1 ! coupler weights for flux from atm to ocean, NO Baltic
flxo2a=flxo2a387x360.8bin_1 ! coupler weights for flux from atm to ocean, NO Baltic
taua2o=taua2o387x360.8bin_1 ! coupler weights for vector from atm to ocean
ssta2o=sata2o387x360.8bin_1 ! interp weights for scalar located in the center atm grid to the center hycom grid
ssto2a=ssto2a387x360.8bin_1 ! coupler weights for sst from ocean to atm
e_o2a=e_o2a387x360.8bin_1   ! coupler weights for eastward vel from ocean to atm
n_o2a=n_o2a387x360.8bin_1   ! coupler weights for northward vel from ocean to atm
cososino=cososino387x360.8bin           ! cos/sin of i,j axis angle on ocean grid
kpar=seawifs_kpar_387x360.tbin          ! monthly/annual seawifs_kpar data
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OCMIP_cfc=OCMIP_cfc.dat
cfle1=abw25b.dat                         ! seawater spectral absorp. and scatt.
coefs
cfle2=acbc25b.dat                        ! phytoplankton spectrl absorp.  and
scatt. coefs
!!!pco2table=pco2.tbl.asc                   ! table to compute pco2 vals from
sst,sss,dic,alk
                                            ! if not defined pCO2_ONLINE
nitrates_inicond=no3_nodc_annmean_387x360.nc    ! initial cond for nitrates (NODC)
silicate_inicond=sio2_nodc_annmean_387x360.nc   ! initial cond for silicate (NODC)
dic_inicond=dic_glodap_annmean_387x360.nc       ! initial cond for dic (GLODAP)
alk_inicond=alk_glodap_annmean_387x360.nc       ! initial cond/forcing for alk (GLODAP)
!!!oasimdirect=oasimdirect_20w_new          ! spectral light components
                                            ! if not defined
atmFe_inicond=iron_gocart_1x1mon_387x360.nc     ! GOCART iron flux
atmFedirect1=iron_ron_195x180_20w.asc    ! Ron Miller's dust fluxes
facirr=facirr.asc                        ! factors for mean irradiance w/in
water
eda_esa_ratios=eda_esa_ratios.asc        ! ratios of radiation spectral
components
!!!!!!!!!!!!!!!!!!! obio_rad  input data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHL_DATA=CHL_WG_2x2.5zavg                !CHL_WG_4x5 or CHL_WG_2x2.5

Label and Namelist:
E4arobio_h4c (2x2.5x40, 1850 atm.;  1x1x26 HYCOM ocean)

DTFIX=180.
&&PARAMETERS
#include "dynamic_ocn_params"
#include "sdragF40_params"
#include "gwdragF40_params"

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

U00a=0.73    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.68    ! below 850mb and MC regions; then tune this to get rad.balance
U00ice=.60       ! tune this first to get: glob. ann. mean plan.alb=30%   (U00ice up=>albedo down)
U00wtrX=1.47     ! this to get: glob. ann. mean net heat at surf. = 0   (U00wtrX+.01=>NetHtSrf+.7)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
madaer=3    ! updated aerosols

#include "atmCompos_1850_params"

! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

! parameters that may have to be changed in emergencies:
DTsrc=1800.
DT=225.
NIsurf=1        ! increase as layer 1 gets thinner

#include "diag_params"

isccp_diags=1   ! use =0 to save cpu time if isccp-diags are not essential
nssw=2         ! until diurnal diags are fixed, Nssw has to be even
nssw=48         ! until diurnal diags are fixed, Nssw has to be even
itest=-1        
jtest=-1       
iocnmx=2        ! default is 0
brntop=50.      ! default is 0.
brnbot=200.     ! default is 300.
diapyn=3.e-7    ! default is 3.e-7
diapyc=.5e-4    ! default is 1.e-4
jerlv0=1

!parameters that affect CO2 gas exchange
!atmCO2=368.6     !uatm for year 2000
!atmCO2=289.9      !uatm for preindustrial runs
atmCO2=0.          !prognostic atmCO2
to_volume_MixRat=1    ! for tracer printout
!!!solFe=0.02            ! default iron solubility
solFe=0.05            ! enhanced iron solubility

&&END_PARAMETERS

 &INPUTZ
   YEARI=1800,MONTHI=01,DATEI=01,HOURI=00, !  from default: IYEAR1=YEARI
   YEARE=1800,MONTHE=01,DATEE=02,HOURE=00, KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1800,MONTHE=01,DATEE=01,HOURE=01,IWRITE=1,JWRITE=1,
/

