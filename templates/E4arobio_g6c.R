E4arobio_g6c.R GISS Model E  coupled version          aromanou  08/04/2009

E4arobio_g6c: obio in GISS ocean based on Larissa's E4F40o32.R
E4F40o32: 2x2.5x40 layers modelE version, 1850 atm.;
   ocean: coupled to GISS ocean model (Russell - Schmidt),
   32 vert. layers in the ocean
   NOTE: new ocean initial condition OIC=OIC.WOA98.2HX2.L32.D1201
   uses turbulence scheme (no dry conv), no grav.wave drag
   time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs
   filters: U,V in E-W direction (after every dynamics time step)
   sea level pressure (after every physics time step)

Preprocessor Options
#define NEW_IO
#define TRACERS_ON                  ! include tracers code
#define CHECK_OCEAN                 ! needed to compile aux/file CMPE002
#define OCN_LAYERING L32
#define OBIO_ON_GISSocean           ! obio on GISS ocean
#define TRACERS_OCEAN               ! GISS Ocean tracers activated
#define TRACERS_OCEAN_INDEP         ! independently defined ocn tracers
#define TRACERS_OceanBiology
#define pCO2_ONLINE
#define TRACERS_GASEXCH_ocean       ! ANY ocean: special tracers to be passed to ocean
#define TRACERS_GASEXCH_ocean_CO2   ! ANY ocean: special tracers to be passed to ocean
!!!!#define CHL_from_SeaWIFs - don't use it - use run-time parameter instead: chl_from_seawifs=1
End Preprocessor Options

Object modules: (in order of decreasing priority)
     ! resolution-specific source codes
Atm144x90                  ! horizontal resolution is 144x90 -> 2x2.5deg
AtmL40                      ! vertical resolution is 40 layers -> 0.1mb
DIAG_RES_F                          ! diagnostics (resolution dependent)
FFT144                              ! utilities
ORES_2Hx2 OFFT144E                  ! ocean horiz res 2x2.5deg

IO_DRV                              ! new i/o 

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV TQUS_DRV                    ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

#include "latlon_source_files"
#include "modelE4_source_files"
#include "dynamic_ocn_source_files"
#include "tracer_shared_source_files"
OCN_Int_LATLON                      ! atm-ocn regrid routines

#include "ocarbon_cycle_oR_files" /* both gas exch and ocean tracer oR model */

Components:
tracers Ent shared MPI_Support solvers giss_LSM dd2d

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB

Data input files:
AIC=AIC.RES_F40.D771201.nc         ! observed init cond (atm. only) ISTART=2
!!! AIC=1JAN1931.rsfE8F40o32        ! Larissa's restart              ISTART=8
!!! AIC=1JAN1918.rsfE31F40o32       ! XXXXXXXXX dont use this AIC, just for test purposess
GIC=GIC.144X90.DEC01.1.ext_1.nc   ! initial ground conditions      ISTART=2

OIC=OIC.E2HX2.L32.D1201.nc         ! Levitus ocean intial conditions
TOPO=Z144X90N_nocasp.1.nc          ! surface fractions and topography
TOPO_OC=OZ144X90N_nocasp.1.nc      ! ocean fraction and topography

!OIC=OIC_compatible_with_TOPO_OC ! Levitus ocean intial conditions
!TOPO=Z2HX2fromZ1QX1N            ! surface fractions and topography
!TOPO_OC=OZ2HX2fromZ1QX1N.nc     ! ocean fraction and topography

OFTAB=OFTABLE_NEW               ! ocean function table
KBASIN=KB144X90.modelE.nc       ! ocean basin designations
OSTRAITS=OSTRAITS_144x90.nml    ! parameterized straits info
CDN=CD144X90.ext.nc             ! neutral drag coefficient
VEG=V144X90_no_crops.ext.nc     ! vegatation file
CROPS=CROPS_144X90N_nocasp.ext.nc ! crops
SOIL=S144X900098M.ext.nc        ! soil properties
REG=REG2X2.5                    ! special regions-diag
RVR=RD_modelE_F.nc             ! river direction file
NAMERVR=RD_modelE_F.names.txt  ! named river outlets

#include "rad_input_files"
#include "rad_144x90_input_files"

TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc             ! topographic variation for gwdrag
MSU_wts=MSU_SSU_RSS_weights.txt
GLMELT=GLMELT_144X90_gas.OCN.nc   ! glacial melt distribution
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
!!!!!!!!!!!!!!!!!!! obio  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cfle1=abw25b.dat                         ! seawater spectral absorp.
                                         ! and scatt. coefs
cfle2=acbc25b.dat                        ! phytoplankton spectrl absorp.
                                         ! and scatt. coefs
!!!!pco2table=pco2.tbl.asc               ! table to compute pco2 values
                                         ! from sst,sss,dic,alk
                                         ! if not defined pCO2_ONLINE
nitrates_inicond=no3_nodc_annmean_90x144.nc    ! initial cond for nitrates (NODC)
silicate_inicond=sio2_nodc_annmean_90x144.nc   ! initial cond for silicate (NODC)
dic_inicond=dic_glodap_annmean_90x144.nc       ! initial cond for dic (GLODAP)
alk_inicond=alk_glodap_annmean_90x144.nc       ! initial cond/forc for alk(GLODAP)
!!!oasimdirect=oasimdirect_20w_new       ! spectral light components
                                         ! if not def OBIO_RAD_coupling
atmFe_inicond=iron_gocart_1x1mon_90x144.nc     ! GOCART iron flux
atmFedirect1=iron_ron_195x180_20w.asc    ! Ron Miller's dust fluxes
facirr=facirr.asc                        ! factors for mean irrad w/in water
eda_esa_ratios=eda_esa_ratios.asc        ! ratios of rad spectrl components
!!!!!!!!!!!!!!!!!!! obio_rad  input data   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHL_DATA=CHL_WG_2x2.5zavg                !CHL_WG_4x5 or CHL_WG_2x2.5
!! CHL_DATA=CHL_WG_2x2.5                    !CHL_WG_4x5 or CHL_WG_2x2.5
                                         !in GISS ocean grid
                                         !to be used with CHL_from_SeaWIFs


Label and Namelist:
E4arobio_g6c (The latest version of modelE)

DTFIX=180
&&PARAMETERS
! parameters set for coupled ocean runs:
KOCEAN=1        ! ocn is prognostic
variable_lk=1
init_flake=1

! parameters usually not changed when switching to coupled ocean:

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_SDRAG=1         ! conserve ang. mom.
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38  39  40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000045 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.2       ! default is 0.5
CDEF=1.95       ! deformation drag coefficient
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF

OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)
OTIDE = 0           !  Ocean tides are not used

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)


U00a=0.74    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.65    ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.60       ! tune this first to get: glob. ann. mean plan.alb=30%   (U00ice up=>albedo down)
U00wtrX=1.47     ! this to get: glob. ann. mean net heat at surf. = 0   (U00wtrX+.01=>NetHtSrf+.7)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
madaer=3    ! updated aerosols
aer_rad_forc=0
cloud_rad_forc=1

#include "atmCompos_1850_params"
variable_orb_par=-2

! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

! parameters that may have to be changed in emergencies:
DTsrc=1800.
DT=225.
DTO=225.
NIsurf=2        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=192       ! use =48 except on halem
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time if isccp-diags are not essential
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
nssw=48         ! obio needs that in order to always restart from hour 0
                ! then we need to do setups for a whole day

!parameters that affect CO2 gas exchange
!!! atmCO2=368.6      !uatm for year 2000
!!! atmCO2=289.9      !uatm for preindustrial runs
atmCO2=0.             !prognostic atmCO2
to_volume_MixRat=1    ! for tracer printout
!!!solFe=0.02            ! default iron solubility
solFe=0.05            ! enhanced iron solubility

&&END_PARAMETERS

 &INPUTZ
!  YEARI=1900,MONTHI=12,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
!  YEARE=2001,MONTHE=12,DATEE=1,HOURE=0, KDIAG=13*0,
!  ISTART=2,IRANDI=0, YEARE=1900,MONTHE=12,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,

   YEARI=1931,MONTHI=1,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=1931,MONTHE=1,DATEE=2,HOURE=0, KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1931,MONTHE=1,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
/
