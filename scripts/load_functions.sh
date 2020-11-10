#!/bin/sh
#Define functions

generalusage(){
 echo -e " Software: Armadillo
 Version: 1.0
 Contact: Pablo Bousquets (bousquetspablo@uniovi.es)

 Usage:
   - armadillo run \t\t\t Run armadillo
   - armadillo data-prep \t\t Create necessary files to run armadillo for a set of ROIs
   - armadillo config-file \t\t Create a configuration file to avoid using arguments in run mode.

 XA Lab - 2020"
}

usage(){ #Create a function to display the help message
    echo -e "
    ### ARMADILLO ###

    Usage: armadillo run [config_file] [options] -N NAME -C control.bam -T tumor.bam \n
    Input arguments:
    \t -N,  --name \t \t \t Case or sample name to analyze
    \t -bd, --bam_dir \t \t Root directory where the genomes are stored
    \t -C,  --control_genome \t \t Control sample genome
    \t -T,  --tumor_genome \t \t Tumour sample genome
    \t -ad, --armadillo_data \t \t Path to armadillo data-prep command output 
    \t -r,  --rois_list \t \t List of ROIs to analyze [all armadillo_data rois] \n
    Analysis arguments:
    \t -cc, --control_coverage \t Coverage of control genome [30]
    \t -tc, --tumor_coverage \t \t Coverage of tumor genome [30]
    \t -cm, --control_threshold \t Maximum variant coverage allowed in the control. [3]
    \t -tt, --tumor_threshold \t Minimum coverage required for a variant to believe it's a good candidate  [6]
    \t -Q,  --base_quality \t \t Minimum base quality required for the tumour genome  [30]
    \t -q,  --mapq \t \t \t Minimum MapQ for reads after being collapsed (most of them should be ~60) [30]
    \t -gc, --GCcutoff \t \t Maximum GC% allowed in the reads  [80] \n
    Other:
    \t -t,  --threads \t \t Threads running in parallel [3]
    \t -S,  --skip \t \t \t Skip bam alignment. Useful to reanalyse a case with other parameters [FALSE]
    \t -p,  --port \t \t \t Port used to perform blat analysis [9006]
    \t -P,  --print \t \t \t Print the variants lost in each step [FALSE]
    \t -h,  --help \t \t \t Display help message
"
}

parse_arguments(){
	while [ $# != 0 ]; do
	    PARAM=`echo $1 | awk -F= '{print $1}'`
	    VALUE=`echo $2 | awk -F= '{print $1}'`
	    case $PARAM in
			-N | --name)
			name=$VALUE
			;;
			-ad | --armadillo_data)
			armadillo_data=$VALUE
			;;
			-R | --ref_genome)
			ref_genome=$VALUE
			;;
			-r | --rois_list)
			rois_list=$VALUE
			;;
			-bd | --bam_dir)
			bam_dir=$VALUE
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
			-cm | --control_threshold)
			control_threshold=$VALUE
			;;
			-tt | --tumor_threshold)
			tumor_threshold=$VALUE
			;;
			-q | --mapq)
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
			-Q | --base_quality)
			base_quality=$(echo $VALUE | tr '[:upper:]' '[:lower:]')
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
		echo "FILE NOT FOUND ERROR:\nThe $2 $3 couldn't be found. Please, check if the path was correctly introduced: $file \n"
		exit 0
	fi
}

check_chr(){
	bam=$1
	chr_names="$(samtools view -H $bam | grep SN:chr | wc -l)"
    if [ "$chr_names" -gt 0 ]
	then 
		echo True
	else
		echo False
	fi

}