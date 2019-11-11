#!/bin/sh
#Define functions
generalusage(){
 echo "Usage:
 - armadillo run \t\t\t Run armadillo
 - armadillo data-prep \t\t\t Create necessary files to run armadillo for a set of ROIs
 - armadillo config-file \t\t Create a configuration file to avoid using arguments in run mode.

Pablo Bousquets - XA Lab
(2019)"
}
usage(){ #Create a function to display the help message
    echo "
    ### ARMADILLO ###
    To run armadillo, two genomes (case and control) are required.
    Also, remind this program was writen for hg19 aligned genomes, so the coordinates provided by default files belong to hg19 genome. \n
    Usage: armadillo run -i ID -C control.bam -T tumor.bam [options] || or || armadillo run configuration_file.txt \n
    Input options:
    \t -i,  --input \t \t \t Case or sample name to analyze
    \t -b,  --bamDir \t \t \t Root directory where the genomes are stored
    \t -C,  --control_genome \t \t Control sample genome
    \t -T,  --tumor_genome \t \t Tumour sample genome
    \t -l,  --list \t \t \t List of ROIs to analyze  \n
    Databases options:
    \t -B,  --blat_coords \t \t Directory where the blat coords of each ROI are stored
    \t -f,  --miniFasta \t \t Directory where the fasta of each exon is stored
    \t -r,  --repeatsDB \t \t Database of genome repeats
    \t -R,  --ref_genome \t \t Repeats reference genome
    \t -s,  --scriptsDir \t \t Directory where this pipeline's scripts are stored \n
    Cutoff options and parameters:
    \t -c,  --control_cov \t \t Minimum coverage with regard to the tumor sample [80%]
    \t -cc, --control_cutoff \t \t Maximum variant coverage allowed in the control. [3]
    \t -ct, --control_contamination \t % tumor cellularity in the control sample [15]
    \t -tc, --tumor_cutoff \t \t Minimum coverage required for a variant to believe it's a good candidate  [6]
    \t -q,  --tum_qual \t \t Minimum base quality required to the tumour genome  [30]
    \t -Q,  --control_qual \t \t Minimum base quality required to the control genome [0]
    \t -m,  --map_qual \t \t Minimum MapQ for reads after being collapsed (note that most of them should be ~60) [40]
    \t -se, --seq_error \t \t Estimation of sequencing error rate [0.00035]
    \t -e,  --max_errors \t \t Maximum sequencing errors in reads [2%]
    \t -rl, --read_length \t \t Reads length [150 bp]
    \t -g,  --gc_content \t \t Maximum GC% allowed in the reads  [80] \n
    Other:
    \t -t,  --threads \t \t Threads running in parallel (only applies if parallelize mode on) [3]
    \t -S,  --skip \t \t \t Skip bam alignment. Useful to reanalyse a case with other parameters [FALSE]
    \t -p,  --port \t \t \t Port used to perform blat analysis [9001]
    \t -P,  --print \t \t \t Print the variants lost step by step [FALSE]
    \t -h,  --help \t \t \t Display help message
"
}

parse_arguments(){
	while [ $# != 0 ]; do
	    PARAM=`echo $1 | awk -F= '{print $1}'`
	    VALUE=`echo $2 | awk -F= '{print $1}'`
	    case $PARAM in
	    	-i | --input)
				case=$VALUE
				;;
	        -s | --scriptsDir)
	            scripts_dir=$VALUE
	            ;;
	        -f | --miniFasta)
				miniFasta_dir=$VALUE
				;;
	        -B | --blat_coords)
	            blat_coords=$VALUE
	            ;;
	      	-R | --ref_genome)
				ref_genome=$VALUE
				;;
	        -r | --repeatsDB)
	            repeatsDB=$VALUE
	            ;;
	        -l | --list)
	            rois_list=$VALUE
	            ;;
	        -e | --max_errors)
	            max_errors=$VALUE
	            ;;
	        -rl | --read_length)
	            read_length=$VALUE
	            ;;
	        -b | --bamDir)
				root_dir=$VALUE
				;;
			-c | --control_cov)
				mincov=$VALUE
				;;
	        -T | --tumor_genome)
				tumor_genome=$VALUE
				;;
	        -C | --control_genome)
				control_genome=$VALUE
				;;
			-cc | --control_cutoff)
				control_cutoff=$VALUE
				;;
               -ct | --control_contamination)
				control_contamination=$VALUE
				;;
			-tc | --tumor_cutoff)
				tumor_cutoff=$VALUE
				;;
			-m | --map_qual)
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
					echo "ERROR: Skip possible arguments are: TRUE or FALSE\n"
					usage
					exit
					;;
				esac
				;;
			-t | --threads)
				threads=$VALUE
				;;
			-g | --gc_content)
				gc_content=$(echo $VALUE | tr '[:upper:]' '[:lower:]')
				;;
			-q | --tum_qual)
				tum_qual=$(echo $VALUE | tr '[:upper:]' '[:lower:]')
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
					echo "ERROR: Skip possible arguments are: TRUE or FALSE\n"
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

extract_minibam(){
	time=$(date +%x%t%X)
	echo ${time}: Generating ${2} minibams... | tee -a pipeline.log
	case=$1
	sample=$2
	type=$3
	blat_coords=$4 #Dir of blat coords of each exon
	miniFasta_dir=$5
	rois_list=$6 #List of regions of interest
	threads=$7
	#list of regions
	lines=$(cat $rois_list | wc -l )

	cat $rois_list | tr ':-' '\t' | while read -r chr st end #Create a tmp minibam for each region
	do
		file=$(echo ${chr}:${st}-${end}) #By reading this way the rois list file, we can use both chr:st-end and bed format
		if [ -f ${blat_coords}/${file} ]
		then
			lines=$((lines - 1))
			#Header
			SM=${sample}
			PL=illumina
			LB=WGS #Or WES
			PU=${sample}
			RG="@RG\\tID:${sample}\\tSM:${SM}\\tPL:${PL}\\tLB:${LB}\\tPU:${PU}"
			allCoords=$(cat ${blat_coords}/${file}) #Save its coords
			samtools view -u -G 0x400 -G 0x4 ${sample} $allCoords | samtools fastq -N - 2>/dev/null|  bwa mem -t ${threads} -R ${RG} ${miniFasta_dir}/${file}.fa - 2>/dev/null | samtools view -uS - > ${type}_tmp_files/${file}.bam
			printf "\r Initial time: ${time}. Files left: $lines        "
		else
			echo "~" ${blat_coords}/${file} does not exist >> pipeline.log
		fi
	done
	time=$(date +%x%t%X)
	echo ${time}: Merging ${case} minibams. This may take a while | tee -a pipeline.log

	iter=0
	cd ${type}_tmp_files
	while [ $(ls | grep -v tmp | wc -l) -ne 0 ]
	do
		iter=$((iter+1))
		bams=$(ls | grep -v tmp | head -n 1000)
		samtools merge tmp_${iter}.bam ${bams}
		rm ${bams}
	done

	if [ ${iter} -gt 1 ]
	then
		samtools merge - tmp_*.bam | samtools sort -@ threads -m 6000000000 -o ../${case}${type}_merged.bam -
		rm tmp_*.bam
	else
		samtools sort -@ ${threads} -m 6000000000 -o ../${case}${type}_merged.bam tmp_1.bam
		rm tmp_1.bam
	fi
	cd ../
	samtools index ${case}${type}_merged.bam


	time=$(date +%x%t%X)
	echo ${time}: Bam generation successful "\n"| tee -a pipeline.log
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
