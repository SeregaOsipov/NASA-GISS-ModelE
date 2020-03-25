#!/bin/bash -l

#SBATCH --job-name=scaleacc
#######SBATCH --output=output.txt
#######SBATCH --error=error.txt
#SBATCH --partition=workq
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --time=23:59:00
####SBATCH --time=00:59:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Sergey.Osipov@kaust.edu.sa

echo "do_pp.sh START"

# THIS WILL FIX THE ncgen: NetCDF: NC_MAX_VARS exceeded
# THIS restriction is remove in netcdf library after v4.5.0
#module unloadn cray-netcdf/4.4.1.1.6
#module load cray-netcdf/4.6.1.3

#run=TobaAO_${scale}_${setup} #where the experiment is
#dest=run_ao_${scale}_${setup} #where to move the output

run=$1
dest=$run
mask=195*
setup=1111

echo "inside do_pp.sh"
echo "run is " $run ", mask is " $mask ", destination is " $dest

rootFP="/project/k1090/osipovs/modelE/ModelE_Support/huge_space"
folderFP=$rootFP/$run
storageRootFP="/project/k1090/osipovs/Data/PinatuboInitialStage/modelE"
exec="srun -n 1 /project/k1090/osipovs/modelE/modelE2/model/mk_diags/scaleacc"

#######################
#process the ACC files
########################

fileMask=$folderFP/*$mask.acc$run.nc
echo "acc, will process these files:"
ls $fileMask

storageFP=$storageRootFP/$dest/acc
echo "will store acc files in: $storageFP"

mkdir -p $storageFP

mv $fileMask $storageFP
fileMask=$storageFP/*$mask.acc$run.nc
echo "acc, calling the scaleacc"

cd $storageFP
for filename in $fileMask; do
  [ -e "$filename" ] || continue
  echo "scaleacc is processing " $filename
  $exec $filename all
done

echo "acc section is done"

# merge the acc output into single netcdf file
/project/k1090/osipovs/Data/PinatuboInitialStage/modelE/do_concatenating_netcdf.sh aijToba __all__
/project/k1090/osipovs/Data/PinatuboInitialStage/modelE/do_concatenating_netcdf.sh aijlToba __all__

taij_diag_vars=ext_band6_N_AKK_1,ext_band6_N_AKK_1_hemis,ext_band6_N_ACC_1,ext_band6_N_ACC_1_hemis,ext_band6_N_DD1_1,ext_band6_N_DD1_1_hemis,ext_band6_N_DS1_1,ext_band6_N_DS1_1_hemis,ext_band6_N_DD2_1,ext_band6_N_DD2_1_hemis,ext_band6_N_DS2_1,ext_band6_N_DS2_1_hemis,ext_band6_N_SSA_1,ext_band6_N_SSA_1_hemis,ext_band6_N_SSC_1,ext_band6_N_SSC_1_hemis,ext_band6_N_OCC_1,ext_band6_N_OCC_1_hemis,ext_band6_N_BC1_1,ext_band6_N_BC1_1_hemis,ext_band6_N_BC2_1,ext_band6_N_BC2_1_hemis,ext_band6_N_BC3_1,ext_band6_N_BC3_1_hemis,ext_band6_N_DBC_1,ext_band6_N_DBC_1_hemis,ext_band6_N_BOC_1,ext_band6_N_BOC_1_hemis,ext_band6_N_BCS_1,ext_band6_N_BCS_1_hemis,ext_band6_N_MXX_1,ext_band6_N_MXX_1_hemis,M_AKK_SU_cond_ls,M_AKK_SU_clwevap_ls,M_AKK_SU_reevap_ls,M_AKK_SU_conclw_ls,M_AKK_SU_precip_ls,M_AKK_SU_washout_ls,M_AKK_SU_washout_ls_hemis,M_AKK_SU_washout_mc,M_AKK_SU_washout_mc_hemis,M_AKK_SU_reevap_mc,M_AKK_SU_precip_mc,M_AKK_SU_conclw_mc,M_AKK_SU_downeva_mc,M_AKK_SU_cond_mc,M_AKK_SU_gs_dep,M_AKK_SU_gs_dep_hemis,M_AKK_SU_dry_dep,M_AKK_SU_dry_dep_hemis,M_AKK_SU_wet_dep,M_AKK_SU_wet_dep_hemis,M_AKK_SU_StratTropflux,M_AKK_SU_StratTropflux_hemis,M_AKK_SU_Total_Mass,M_AKK_SU_Total_Mass_hemis,M_OCC_SU_cond_ls,M_OCC_SU_clwevap_ls,M_OCC_SU_reevap_ls,M_OCC_SU_conclw_ls,M_OCC_SU_precip_ls,M_OCC_SU_washout_ls,M_OCC_SU_washout_ls_hemis,M_OCC_SU_washout_mc,M_OCC_SU_washout_mc_hemis,M_OCC_SU_reevap_mc,M_OCC_SU_precip_mc,M_OCC_SU_conclw_mc,M_OCC_SU_downeva_mc,M_OCC_SU_cond_mc,M_OCC_SU_gs_dep,M_OCC_SU_gs_dep_hemis,M_OCC_SU_dry_dep,M_OCC_SU_dry_dep_hemis,M_OCC_SU_wet_dep,M_OCC_SU_wet_dep_hemis,M_OCC_SU_StratTropflux,M_OCC_SU_StratTropflux_hemis,M_OCC_SU_Total_Mass,M_OCC_SU_Total_Mass_hemis,M_ACC_SU_cond_ls,M_ACC_SU_clwevap_ls,M_ACC_SU_reevap_ls,M_ACC_SU_conclw_ls,M_ACC_SU_precip_ls,M_ACC_SU_washout_ls,M_ACC_SU_washout_ls_hemis,M_ACC_SU_washout_mc,M_ACC_SU_washout_mc_hemis,M_ACC_SU_reevap_mc,M_ACC_SU_precip_mc,M_ACC_SU_conclw_mc,M_ACC_SU_downeva_mc,M_ACC_SU_cond_mc,M_ACC_SU_gs_dep,M_ACC_SU_gs_dep_hemis,M_ACC_SU_dry_dep,M_ACC_SU_dry_dep_hemis,M_ACC_SU_wet_dep,M_ACC_SU_wet_dep_hemis,M_ACC_SU_StratTropflux,M_ACC_SU_StratTropflux_hemis,M_ACC_SU_Total_Mass,M_ACC_SU_Total_Mass_hemis,SO2_Total_Mass,SO2_Total_Mass_hemis,N2O5_Total_Mass,N2O5_Total_Mass_hemis,O3_Total_Mass,O3_Total_Mass_hemis,HNO3_Total_Mass,HNO3_Total_Mass_hemis

# the 0111 experiments has different populations output
if [ $setup == 0111 ] || [ $setup == 0111_P ]
then
	taij_diag_vars=ext_band6_NO3p,ext_band6_Silt4,ext_band6_Silt3,ext_band6_Silt2,ext_band6_Silt1,ext_band6_Clay4,ext_band6_Clay3,ext_band6_Clay2,ext_band6_Clay1,ext_band6_OCB,ext_band6_OCIA,ext_band6_BCB,ext_band6_BCIA,ext_band6_SO4,ext_band6_seasalt2,ext_band6_seasalt1,ext_band6_isopp1a,ext_band6_NO3p_hemis,ext_band6_Silt4_hemis,ext_band6_Silt3_hemis,ext_band6_Silt2_hemis,ext_band6_Silt1_hemis,ext_band6_Clay4_hemis,ext_band6_Clay3_hemis,ext_band6_Clay2_hemis,ext_band6_Clay1_hemis,ext_band6_OCB_hemis,ext_band6_OCIA_hemis,ext_band6_BCB_hemis,ext_band6_BCIA_hemis,ext_band6_SO4_hemis,ext_band6_seasalt2_hemis,ext_band6_seasalt1_hemis,ext_band6_isopp1a_hemis
fi
/project/k1090/osipovs/Data/PinatuboInitialStage/modelE/do_concatenating_netcdf.sh taijToba $taij_diag_vars

/project/k1090/osipovs/Data/PinatuboInitialStage/modelE/do_concatenating_netcdf.sh taijlToba __all__
/project/k1090/osipovs/Data/PinatuboInitialStage/modelE/do_concatenating_netcdf.sh tajlToba __all__

echo "concatenating is done"

# echo "COMMENT EXIT to bring back the SUBDD processing"
# exit 0

#######################
#process the SUBDD files
########################

fileMask=$folderFP/*$mask.subdd$run.nc
echo "subdd, ill process these files:"
ls $fileMask

storageFP=$storageRootFP/$dest/subdd
echo "will store SUBDD files in:"
echo $storageFP

mkdir -p $storageFP

mv $fileMask $storageFP
fileMask=$storageFP/*$mask.subdd$run.nc

echo "subdd, calling the scaleacc"

cd $storageFP
for filename in $fileMask; do
  [ -e "$filename" ] || continue
  echo "scaleacc is processing " $filename
  $exec $filename all
done

echo "subdd section is done"

echo "do_pp.sh script is DONE"
#srun /project/k1090/osipovs/modelE/modelE2/model/mk_diags/pdE_nc Toba DEC1952.accToba*.nc
#srun /project/k1090/osipovs/modelE/modelE2/model/mk_diags/scaleacc JAN1950.subddTobaAO_1110.nc all
