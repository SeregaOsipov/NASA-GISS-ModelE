#!/bin/bash -l

#SBATCH --job-name=merge_netcdf
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --partition=workq
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --time=60:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Sergey.Osipov@kaust.edu.sa

mask=$1 #such as aijToba or taijlToba
vars=$2 #list of variables
echo "the file mask is: $mask"
echo "the list of selected variables is: $vars"
echo "current folder is: $PWD"

files=$(ls *.${mask}*.nc | sort -k1.4,1.7 -k1.1,1.3M)
echo -e "script will concatenate the following files:\n $files"
mkdir -p merged

#remove the previous temp file if it exists
rm -f merged/merged.nc

if [ "$vars" == "__all__" ]
then
	echo "running ncecat WITHOUT vars selection"
	ncecat $files merged/merged.nc
else
	echo "running ncecat WITH vars selection"
	ncecat -v $vars $files merged/merged.nc
fi
echo "ncecat done"
mv merged/merged.nc merged/${mask}_merged.nc
#ncecat $files merged/${mask}_merged.nc

#rename newly created record dimension into time
ncrename -d record,time merged/${mask}_merged.nc

echo "concatenate script is DONE"
