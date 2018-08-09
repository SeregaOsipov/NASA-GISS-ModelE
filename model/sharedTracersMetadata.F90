#include "rundeck_opts.h"
!------------------------------------------------------------------------------
module sharedTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  sharedTracersMetadata_mod encapsulates the metadata shared among various
!@+    tracers.
!@auth NCCS ASTG
  use OldTracer_mod, only: nPart
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_HSTAR
  use OldTracer_mod, only: set_F0
  use OldTracer_mod, only: set_tr_RKD
  use OldTracer_mod, only: set_tr_DHD
  use OldTracer_mod, only: set_trdecay
  use OldTracer_mod, only: tr_RKD 
  use OldTracer_mod, only: set_needtrs
#ifdef TRACERS_SPECIAL_Lerner
  use TRACERS_MPchem_COM, only: nMPtable
  use OldTracer_mod, only: set_iMPtable
  use OldTracer_mod, only: set_tcscale
#endif  /* TRACERS_SPECIAL_Lerner */
  use OldTracer_mod, only: dodrydep
  use OldTracer_mod, only: F0
  use OldTracer_mod, only: HSTAR
  use OldTracer_mod, only: ngas, nPART
  use OldTracer_mod, only: set_emisPerFireByVegType
  use OldTracer_mod, only: set_pm2p5fact
  use OldTracer_mod, only: set_pm10fact
  use OldTracer_mod, only: set_has_chemistry
  use TRACER_COM, only : whichEPFCs, seasonalNH3src
  use TRACER_COM, only: n_H2O2, n_NH3,  n_NH4, n_DMS, n_SO2, n_H2O2_s, &
    n_CH4, n_N2O, n_Rn222
  use Dictionary_mod, only: sync_param
  use RunTimeControls_mod, only: tracers_drydep
  use RunTimeControls_mod, only: tracers_special_lerner
  use RunTimeControls_mod, only: dynamic_biomass_burning  
  implicit none
  private 

  public DMS_setSpec
  public SO2_setSpec
  public H2O2_setSpec
  public NH3_setSpec
  public NH4_setSpec
  public H2O2_s_setSpec
  public CH4_setSpec
  public N2O_setSpec
  public Rn222_setSpec

!@param convert_HSTAR converts from mole/Joule to mole/(L*atm)
  real(8), parameter :: convert_HSTAR = 1.01325d2
  public convert_HSTAR

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  subroutine DMS_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_DMS = n
    call set_ntm_power(n, -12)
    call set_tr_mm(n, 62.d+0)
    call set_needtrs(n, .true.)
    call set_has_chemistry(n, .true.)
  end subroutine DMS_setSpec

  subroutine SO2_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_SO2 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 64.d+0)
    call set_tr_RKD(n, 0.0118d0 ) !mole/J or  1.2  M/atm
    call set_tr_DHD(n, -2.62d4) ! in J/mole= -6.27 kcal/mol
    call set_tr_wd_type(n, ngas)
    if (tracers_drydep) CALL SET_HSTAR(N, 1.D5)
#ifdef DYNAMIC_BIOMASS_BURNING
    if (dynamic_biomass_burning) then
      ! 12 below are the 12 VDATA veg types or Ent remapped to them,
      ! from Olga Pechony's EPFC.xlsx e-mailed to Greg 1/13/2013
      ! Note that she also provided numbers for SO4, but we don't use
      ! those.
      call sync_param("whichEPFCs",whichEPFCs)
      select case(whichEPFCs)
      case(1) ! AR5
        call set_emisPerFireByVegType(n, [0.d0,1.11d-7,4.96d-8,3.22d-8, &
        & 7.63d-8,7.69d-8,7.32d-8,1.25d-7,0.d0,0.d0,0.d0,0.d0] )
      case(2) ! GFED3
        call set_emisPerFireByVegType(n, [0.d0,9.95d-8,7.46d-8,1.16d-8, &
        & 3.14d-8,4.98d-8,9.98d-8,7.35d-8,0.d0,0.d0,0.d0,0.d0] )
      case(3) ! GFED2
        call set_emisPerFireByVegType(n, [0.d0,6.84d-8,2.46d-8,3.13d-8, &
        & 3.68d-8,4.23d-8,7.69d-8,7.19d-8,0.d0,0.d0,0.d0,0.d0] )
      case(4) ! MOPITT
        call set_emisPerFireByVegType(n, [0.d0,2.75d-8,8.83d-8,1.85d-8, &
        & 5.94d-8,8.19d-8,9.29d-9,2.24d-8,0.d0,0.d0,0.d0,0.d0] )
      case default
        call stop_model('whichEPFCs unknown',255)
      end select
    end if
#endif
    call set_has_chemistry(n, .true.)
  end subroutine SO2_setSpec

  subroutine H2O2_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_H2O2 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 34.016d0)
    call set_tr_RKD(n, 9.869d2    ) ! in mole/J = 1.d5 mole/(L atm)
    call set_tr_DHD(n, -5.52288d4 ) ! in J/mole = -13.2 kcal/mole.
    if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    call set_F0(n,  1.d0)
    call set_has_chemistry(n, .true.)
  end subroutine H2O2_setSpec

  subroutine NH3_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_NH3 = n
    call set_ntm_power(n, -10)
    call set_tr_mm(n, 17.d0)
    call set_tr_RKD(n, 100.d0) ! higher than nominal; effective Henry
!    call set_tr_RKD(n, 0.7303d0   ) !tr_RKD=74 M/atm
    call set_tr_DHD(n, -2.84d4  ) !tr_DHD=-6.80 kcal/mole
    call set_tr_wd_type(n, ngas)
    if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    call sync_param("seasonalNH3src", seasonalNH3src)
#ifdef DYNAMIC_BIOMASS_BURNING
    if (dynamic_biomass_burning) then
      ! 12 below are the 12 VDATA veg types or Ent remapped to them,
      ! from Olga Pechony's EPFC.xlsx e-mailed to Greg 1/13/2013
      call sync_param("whichEPFCs",whichEPFCs)
      select case(whichEPFCs)
      case(1) ! AR5
        call set_emisPerFireByVegType(n, [0.d0,4.65d-7,1.43d-7,1.06d-7, &
        & 2.06d-7,1.94d-7,2.07d-7,3.95d-7,0.d0,0.d0,0.d0,0.d0] )
      case(2) ! GFED3
        call set_emisPerFireByVegType(n, [0.d0,4.06d-7,2.91d-7,1.30d-7, &
        & 6.24d-8,1.56d-7,2.11d-7,1.44d-7,0.d0,0.d0,0.d0,0.d0] )
      case(3) ! GFED2
        call set_emisPerFireByVegType(n, [0.d0,1.70d-7,1.23d-7,1.46d-7, &
        & 1.40d-7,1.41d-7,1.42d-7,1.89d-7,0.d0,0.d0,0.d0,0.d0] )
      case(4) ! MOPITT
        call set_emisPerFireByVegType(n, [0.d0,1.60d-7,3.07d-7,8.51d-8, &
        & 2.16d-7,2.25d-7,8.28d-8,7.31d-8,0.d0,0.d0,0.d0,0.d0] )
      case default
        call stop_model('whichEPFCs unknown',255)
      end select
    end if
#endif
    call set_has_chemistry(n, .true.)
  end subroutine NH3_setSpec

  subroutine H2O2_s_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_H2O2_s = n
    call set_ntm_power(n, -10)
    call set_tr_mm(n, 34.016d0)
    call set_tr_RKD(n, 986.9d0)
    call set_tr_DHD(n, -5.52288d4 ) ! in J/mole = -13.2 kcal/mole.
    call set_tr_wd_type(n, ngas)
    if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    call set_F0(n,  1.d0)
    call set_has_chemistry(n, .true.)
  end subroutine H2O2_s_setSpec

    subroutine CH4_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_CH4 = n
      call set_tr_mm(n, 16.d0)
#ifdef TRACERS_SPECIAL_Lerner
      if (tracers_special_lerner) then
        call set_ntm_power(n, -9)
        nMPtable=nMPtable+1
        call set_iMPtable(n, nMPtable)
        call set_tcscale(n, 1.d0)
      end if
#endif
      call set_ntm_power(n, -8)

#ifdef DYNAMIC_BIOMASS_BURNING
    if (dynamic_biomass_burning) then
      call sync_param("whichEPFCs",whichEPFCs)
      ! 12 below are the 12 VDATA veg types or Ent remapped to them,
      ! from Olga Pechony's EPFC.xlsx e-mailed to Greg 1/13/2013
      select case(whichEPFCs)
      case(1) ! AR5
        call set_emisPerFireByVegType(n, [0.d0,5.06d-7,3.74d-7,1.90d-7, &
        & 5.62d-7,5.11d-7,4.29d-7,7.98d-7,0.d0,0.d0,0.d0,0.d0] )
      case(2) ! GFED3
        call set_emisPerFireByVegType(n, [0.d0,1.26d-6,4.55d-7,1.37d-7, &
        & 2.20d-7,3.17d-7,5.32d-7,8.21d-7,0.d0,0.d0,0.d0,0.d0] )
      case(3) ! GFED2
        call set_emisPerFireByVegType(n, [0.d0,3.92d-7,1.74d-7,1.68d-7, &
        & 2.69d-7,2.58d-7,3.62d-7,8.59d-7,0.d0,0.d0,0.d0,0.d0] )
      case(4) ! MOPITT
        call set_emisPerFireByVegType(n, [0.d0,3.25d-7,2.38d-7,1.29d-7, &
        & 3.77d-7,3.32d-7,2.63d-7,4.19d-7,0.d0,0.d0,0.d0,0.d0] )
      case default
        call stop_model('whichEPFCs unknown',255)
      end select
    end if
#endif
    call set_has_chemistry(n, .true.)
    end subroutine CH4_setSpec

    subroutine N2O_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_N2O = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, 44.d0)
#ifdef TRACERS_SPECIAL_Lerner
      if (tracers_special_lerner) then
        nMPtable=nMPtable+1
        call set_iMPtable(n, nMPtable)
        call set_tcscale(n, 1.d0)
      end if
#endif
      call set_has_chemistry(n, .true.)
    end subroutine N2O_setSpec

    subroutine Rn222_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Rn222 = n 
      call set_ntm_power(n, -21)
      call set_tr_mm(n, 222.d0)
      call set_trdecay(n,  2.1d-6)
      call set_has_chemistry(n, .true.)
    end subroutine Rn222_setSpec

    subroutine NH4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_NH4 = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 18.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine NH4_setSpec


  ! TOMAS duplicates with nitrate: NH3_setSpec, NH4_setSpec

end module sharedTracersMetadata_mod
