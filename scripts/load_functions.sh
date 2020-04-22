#!/bin/sh
#Define functions
generalusage(){
 echo " Program: Armadillo
 Version: 1.0
 Contact: Pablo Bousquets (bousquetspablo@uniovi.es)

 Usage:
   - armadillo run \t\t\t Run armadillo
   - armadillo data-prep \t\t Create necessary files to run armadillo for a set of ROIs
   - armadillo config-file \t\t Create a configuration file to avoid using arguments in run mode.

 XA Lab - 2019"
}

usage(){ #Create a function to display the help message
    echo "
    ### ARMADILLO ###
    To run armadillo, two genomes (tumor and control) are required.
    Also, remember that this program was writen for hg19 aligned genomes, so the coordinates provided by default files belong to hg19 genome. \n
    Usage: armadillo run -i ID -C control.bam -T tumor.bam [options] || or || armadillo run configuration_file.txt \n
    Input arguments:
    \t -n,  --name \t \t \t Case or sample name to analyze
    \t -rd, --root_dir \t \t Root directory where the genomes are stored
    \t -C,  --control_genome \t \t Control sample genome
    \t -T,  --tumor_genome \t \t Tumour sample genome
    \t -r,  --rois_list \t \t List of ROIs to analyze  \n
    Databases options:
    \t -ad, --armadillo_data \t \t Path to armadillo data-prep command output
    \t -r,  --repeatsDB \t \t Database of genome repeats
    \t -s,  --scripts_dir \t \t Directory where this pipeline's scripts are stored \n
    Analysis parameters:
    \t -cc, --control_coverage \t Coverage of control genome [30]
    \t -tc, --tumor_coverage \t \t Coverage of tumor genome [30]
    \t -cm, --control_max \t \t Maximum variant coverage allowed in the control. [3]
    \t -tt, --tumor_threshold \t Minimum coverage required for a variant to believe it's a good candidate  [6]
    \t -q,  --base_quality \t \t Minimum base quality required to the tumour genome  [30]
    \t -Q,  --control_qual \t \t Minimum base quality required to the control genome [0]
    \t -m,  --mapq \t \t \t Minimum MapQ for reads after being collapsed (most of them should be ~60) [40]
    \t -rl, --read_length \t \t Reads length [150]
    \t -gc, --GCcutoff \t \t Maximum GC% allowed in the reads  [80] \n
    Other:
    \t -t,  --threads \t \t Threads running in parallel [3]
    \t -S,  --skip \t \t \t Skip bam alignment. Useful to reanalyse a case with other parameters [FALSE]
    \t -p,  --port \t \t \t Port used to perform blat analysis [9001]
    \t -P,  --print \t \t \t Print the variants lost in each step [FALSE]
    \t -h,  --help \t \t \t Display help message
"
}

parse_arguments(){
	while [ $# != 0 ]; do
	    PARAM=`echo $1 | awk -F= '{print $1}'`
	    VALUE=`echo $2 | awk -F= '{print $1}'`
	    case $PARAM in
			-n | --name)
			name=$VALUE
			;;
			-s | --scripts_dir)
			scripts_dir=$VALUE
			;;
			-ad | --armadillo_data)
			armadillo_data=$VALUE
			;;
			-R | --ref_genome)
			ref_genome=$VALUE
			;;
			-r | --repeatsDB)
			repeatsDB=$VALUE
			;;
			-l | --rois_list)
			rois_list=$VALUE
			;;
			-rl | --read_length)
			read_length=$VALUE
			;;
			-rd | --root_dir)
			root_dir=$VALUE
			;;
			-cc | --control_coverage)
			control_coverage=$VALUE
			;;
			-T | --tumor_genome)
			tumor_genome=$VALUE
			;;
			-C | --control_genome)
			control_genome=$VALUE
			;;
			-tc | --tumor_coverage)
			tumor_coverage=$VALUE
			;;
			-cm | --control_max)
			control_max=$VALUE
			;;
			-tt | --tumor_threshold)
			tumor_threshold=$VALUE
			;;
			-m | --mapq)
			mapq=$VALUE
			;;
			-p | --port)
			port=$VALUE
			;;
			-P | --print)
			VALUE=$(echo $VALUE | tr '[:upper:]' '[:lower:]')
			case "$VALUE" in
			    "true" | "false")
			    print=$VALUE
			    ;;
			    *)
			    echo "ERROR: --print argument possible values are: TRUE or FALSE\n"
			    usage
			    exit
			    ;;
			esac
			;;
			-t | --threads)
			threads=$VALUE
			;;
			-gc | --GCcutoff)
			GCcutoff=$(echo $VALUE | tr '[:upper:]' '[:lower:]')
			;;
			-q | --base_quality)
			base_quality=$(echo $VALUE | tr '[:upper:]' '[:lower:]')
			;;
			-Q | --control_qual)
			control_qual=$(echo $VALUE | tr '[:upper:]' '[:lower:]')
			;;
			-tc | --tumor_cutoff)
			tumor_cutoff=$VALUE
			;;
			-S | --skip)
			VALUE=$(echo $VALUE | tr '[:upper:]' '[:lower:]')
			case "$VALUE" in
			    "true" | "false")
			    skip=$VALUE
			    ;;
			    *)
			    echo "ERROR: --skip argument possible values are: TRUE or FALSE\n"
			    usage
			    exit
			    ;;
			esac
			;;
			-h | --help)
			usage
			exit
			;;
			*)
			echo "ERROR: unknown parameter \"$PARAM\" \n"
			usage
			exit 1
			;;
		esac
		shift
		shift
	done
}

evaluate(){
	file=$1
	tag="$2 $3"
	if [ ! -f $file ]
	then
		echo "${tag} ERROR:\nThe file couldn't be found. Please, check if the path was correctly introduced: $file \n"
		exit 0
	fi
}

queue_threads(){
	subp=$(($(ps --no-headers -o pid --ppid=$1 |wc -w)-1))
	if [ "${subp}" -gt "${2}" ]
	then
		wait -n
	fi
}
