#!/bin/bash
## This script will check the dependencies required for Armadillo ##

echo Checking system dependencies...
for program in samtools bwa gfClient gfServer python3
do
	echo -ne "\t - ${program}..."
	[[ $(which ${program}) ]] && echo -e OK || { echo -e "Error: $program needs to be installed. Use 'sudo apt-get install $program'"; }
done

echo Checking python dependencies...
for module in argparse gzip multiprocessing os pyfaidx re sys statistics subprocess
do 
	echo -ne "\t - ${module}... "
	[[ $(python3 -c 'import pkgutil; print(1 if pkgutil.find_loader('"'$module'"') else 0)') -eq 1 ]] && echo -e OK || pip3 install ${module}
done

##Set the installation path##
installation_path="installation_path=$(dirname $(realpath $0))/"
sed -i -e "4s:.*:${installation_path}:g" $(dirname $0)/scripts/default_config.txt

