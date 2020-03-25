#!/bin/bash -l

#SBATCH --job-name=do_ensemble_pp
#######SBATCH --output=output.txt
#######SBATCH --error=error.txt
#SBATCH --partition=workq
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --time=23:59:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Sergey.Osipov@kaust.edu.sa

echo "do_ensemble_pp.sh STARTS"

# temp fix while nco is not built for intel enviroment on shaheen
module swap  PrgEnv-intel PrgEnv-gnu
module load nco

mask=195* #basically a year, 1950

# REMEMBER, that C run does not scale (no 10x prefix)
scale=0x
setups=C

# 1x Pinatubo
setups=1111
scale=1x
scale=3x

# 10x Pinatubo
#setups=1111,1110,1101,1011,0111
#scale=10x

# 100x Pinatubo, Toba
#scale=100x
#setups=1111,1110,1101,1011,0111,0111_P

for setup in ${setups//,/ }
do
	run=TobaAO_${scale}_${setup} #where the experiment is
	dest=run_ao_${scale}_${setup} #where to move the output
	echo "run is " $run ", scale is" $scale ", mask is " $mask ", destination is " $dest

	storageRootFP="/project/k1090/osipovs/Data/PinatuboInitialStage/modelE"
	storageFP=$storageRootFP/$dest/acc
	#cd $storageFP
	echo "current folder is: $PWD"

	/project/k1090/osipovs/Data/PinatuboInitialStage/modelE/do_pp.sh $setup $scale $mask >& logs/log.${run}

	#run the diagnostics of your choice here
	#/project/k1090/osipovs/Data/PinatuboInitialStage/modelE/do_concatenating_netcdf.sh aijToba __all__
	#/project/k1090/osipovs/Data/PinatuboInitialStage/modelE/do_concatenating_netcdf.sh aijlToba __all__

	#taij_diag_vars=ext_band6_N_AKK_1,ext_band6_N_AKK_1_hemis,ext_band6_N_ACC_1,ext_band6_N_ACC_1_hemis,ext_band6_N_DD1_1,ext_band6_N_DD1_1_hemis,ext_band6_N_DS1_1,ext_band6_N_DS1_1_hemis,ext_band6_N_DD2_1,ext_band6_N_DD2_1_hemis,ext_band6_N_DS2_1,ext_band6_N_DS2_1_hemis,ext_band6_N_SSA_1,ext_band6_N_SSA_1_hemis,ext_band6_N_SSC_1,ext_band6_N_SSC_1_hemis,ext_band6_N_OCC_1,ext_band6_N_OCC_1_hemis,ext_band6_N_BC1_1,ext_band6_N_BC1_1_hemis,ext_band6_N_BC2_1,ext_band6_N_BC2_1_hemis,ext_band6_N_BC3_1,ext_band6_N_BC3_1_hemis,ext_band6_N_DBC_1,ext_band6_N_DBC_1_hemis,ext_band6_N_BOC_1,ext_band6_N_BOC_1_hemis,ext_band6_N_BCS_1,ext_band6_N_BCS_1_hemis,ext_band6_N_MXX_1,ext_band6_N_MXX_1_hemis
	#if [ $setup == 0111 ] || [ $setup == 0111_P ]
	#then
	#	echo "special set of variables for 0111*"
	#	taij_diag_vars=ext_band6_NO3p,ext_band6_Silt4,ext_band6_Silt3,ext_band6_Silt2,ext_band6_Silt1,ext_band6_Clay4,ext_band6_Clay3,ext_band6_Clay2,ext_band6_Clay1,ext_band6_OCB,ext_band6_OCIA,ext_band6_BCB,ext_band6_BCIA,ext_band6_SO4,ext_band6_seasalt2,ext_band6_seasalt1,ext_band6_isopp1a,ext_band6_NO3p_hemis,ext_band6_Silt4_hemis,ext_band6_Silt3_hemis,ext_band6_Silt2_hemis,ext_band6_Silt1_hemis,ext_band6_Clay4_hemis,ext_band6_Clay3_hemis,ext_band6_Clay2_hemis,ext_band6_Clay1_hemis,ext_band6_OCB_hemis,ext_band6_OCIA_hemis,ext_band6_BCB_hemis,ext_band6_BCIA_hemis,ext_band6_SO4_hemis,ext_band6_seasalt2_hemis,ext_band6_seasalt1_hemis,ext_band6_isopp1a_hemis
	#fi
	#/project/k1090/osipovs/Data/PinatuboInitialStage/modelE/do_concatenating_netcdf.sh taijToba $taij_diag_vars
done

echo "waiting for child processes to finish"

wait

echo "DONE do_ensemble_pp"
exit 0
