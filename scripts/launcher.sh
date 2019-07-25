#!/bin/sh

##Load main variables
STARTTIME=$(date +%s)
case='NA'
root_dir=''
scripts_dir='/home/pbousquets/scripts/RR_VaC/scripts'
blat_coords='/home/pbousquets/scripts/RR_VaC/lib/blatCoords'
miniFasta_dir='/home/pbousquets/scripts/RR_VaC/lib/miniFASTA'
ref_genome='/home/pbousquets/scripts/RR_VaC/lib/repetitive_regions.fa'
rois_list="/home/pbousquets/scripts/RR_VaC/lib/ROIs"
repeatsDB='/home/pbousquets/scripts/RR_VaC/lib/hg19_repeatmasker_DB_sorted.bed.gz'
tumor_genome='NA'
control_genome='NA'
tumor_cutoff=6
control_cutoff=3
gc_content=80
tum_qual=30
control_qual=0
mincov=80
mapq=40
threads=3
port='9001'
skip="false"
print='false'
printopt=''
##Load scripts
repeatmasker_candidates_filter=${scripts_dir}'/repeatmasker_candidates_filter.py'
mq2vcf_TDvsND=${scripts_dir}'/mq2vcf_TDvsND.py'

#Define functions
usage(){ #Create a function to display the help message
    echo "
    To run the program, two genomes (case and control) are required. Please, name them as shown below.
    Also, remind this program was writen for hg19 aligned genomes, so the coordinates provided by default files belong to hg19 genome. \n
    Usage: sh thisscript.sh -i case -C control.bam -T tumor.bam [options] \n
    Input options: 
    \t -i,  --input \t \t Case or sample name to analyze  
    \t -b,  --bamDir \t \t Root directory where the genomes are stored 
    \t -C,  --control_genome \t Control sample genome 
    \t -T,  --tumor_genome \t Tumour sample genome
    \t -l,  --list \t \t List of ROIs to analyze  \n
    Databases options:
    \t -B,  --blat_coords \t Directory where the blat coords of each ROI are stored 
    \t -f,  --miniFasta \t Directory where the fasta of each exon is stored  
    \t -r,  --repeatsDB \t Database of genome repeats  
    \t -R,  --ref_genome \t Repeats reference genome 
    \t -s,  --scriptsDir \t Directory where this pipeline's scripts are stored \n
    Cutoff options:
    \t -c,  --control_coverage \t Minimum proportional coverage with regard to the tumor sample [80 (%)]
    \t -cc, --control_cutoff \t Maximum coverage allowed for a variant in the control. Warning: for blood samples, a 0 cutoff may lead to loss of candidates [3] 
    \t -tc, --tumor_cutoff \t Minimum coverage required for a variant to believe it's a good candidate  [6]
    \t -q,  --tum_qual \t Minimum base quality required to the tumour genome  [30]
    \t -Q,  --control_qual \t Minimum base quality required to the control genome [0] 
    \t -m,  --map_qual \t Minimum MapQ for reads after being collapsed (note that most of them should be ~60) [40] 
    \t -g,  --gc_content \t Maximum GC% allowed in the reads  [80] \n
    Other:
    \t -t,  --threads \t Threads running in parallel (only applies if parallelize mode on) [3]    
    \t -S,  --skip \t \t Skip bam alignment. Useful to reanalyse a case with other parameters [FALSE]
    \t -p,  --port \t \t Port used to perform blat analysis [9001]
    \t -P,  --print \t \t Print the variants lost step by step [FALSE]
    \t -h,  --help \t \t Display help message
"
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
	
	cat $rois_list | while read -r file #Create a tmp minibam for each region
	do
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
	        samtools view -u -G 0x400 ${sample} $allCoords | samtools view -u -G 0x4 - | (samtools fastq -N - 2> /dev/null) |  (bwa mem -t ${threads} -R ${RG} ${miniFasta_dir}/${file}.fa - 2> /dev/null) | samtools view -uS - > tmp_files/${file}.bam 
	        printf "\r Initial time: ${time}. Files left: $lines        "
	    else
	        echo "~" ${blat_coords}/${file} does not exist >> pipeline.log
	    fi
	done
	time=$(date +%x%t%X)
	echo ${time}: Merging ${case} minibams. This may take a while | tee -a pipeline.log

	iter=0
	cd tmp_files
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
		echo "$tag ERROR:\nThe file couldn't be found. Please, check if the path was correctly introduced: $file \n"
		exit 0
	fi
}

##Display help
if [ $# -eq 0 ] #If no arguments provided, display help
then 
	usage
	exit 1
fi

##Parse arguments
while [ "$1" != "" ]; do
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
        -b | --bamDir)
			root_dir=$VALUE
			;;
		-c | --control_coverage)
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

if [ "$print" = 'true' ]
then
	printopt='-f'
fi 

##Check required variables
if [ "$case" = 'NA' ]
then
	echo "ERROR: No input provided."
	usage
	exit
fi

##Create the directories and start the case
if [ -d ${case} ]
then 
	cd ${case}
else
	mkdir ${case} 
	cd ${case}
	mkdir tmp_files
fi	

TD=${root_dir}${tumor_genome}
ND=${root_dir}${control_genome}
NDpileup="${case}control_merged.bam"			
TDpileup="${case}tumor_merged.bam"

##Check if files exist
evaluate $ref_genome Reference genome
evaluate $repeatsDB repeatsDB

##Run pipeline
echo "\nYour options: \n  Input: $case \n  Root dir: $root_dir \n  Control genome: $control_genome \n  Tumor genome: $tumor_genome \n  List of ROIs: $rois_list \n  Blat coordinates directory: $blat_coords \n  miniFasta directory: $miniFasta_dir \n  RepeatsDB: $repeatsDB \n  Reference genome: $ref_genome \n  Scripts directory: $scripts_dir \n  Min. coverage: $mincov \n  Control cutoff: $control_cutoff \n  Tumor cutoff: $tumor_cutoff \n  Control quality: $control_qual \n  Tumor quality: $tum_qual \n  Mapping quality cutoff: $mapq \n  GC content cutoff: $gc_content \n  Threads: $threads \n  Print: $print \n  Skip: $skip \n  Port: $port \n" | tee -a pipeline.log

if [ ${skip} = 'false' ]
then
	##Check if needed files exist
	evaluate $rois_list ROIs genome
	evaluate $TD Tumor genome
	evaluate $ND Control_genome

	{ extract_minibam $case ${TD} tumor ${blat_coords} ${miniFasta_dir} ${rois_list} ${threads} ; } 2>> pipeline.log
	{ extract_minibam $case ${ND} control ${blat_coords} ${miniFasta_dir} ${rois_list} ${threads} ; } 2>> pipeline.log
	rm -rf tmp_files

else
	echo "Skipped minibam extraction." | tee -a pipeline.log
fi 
	
time=$(date +%x%t%X)
echo ${time}: Finding candidates...

if [ -f ${case}_candidates.vcf ]
then
	rm ${case}_candidates.vcf 
fi 

samtools mpileup --output-QNAME -Q ${control_qual} -q ${mapq} -R -f ${ref_genome} ${TDpileup} ${NDpileup} | python3 ${mq2vcf_TDvsND} -i - -b ${TDpileup} -n ${case} -td ${tumor_cutoff} -nd $control_cutoff -gc ${gc_content} -q ${tum_qual} -r ${ref_genome} -c ${mincov} -t ${threads} -p ${port} ${printopt} | python3 ${repeatmasker_candidates_filter} $repeatsDB $tumor_cutoff 95 ${printopt}> ${case}_candidates.vcf

lines=$(wc -l ${case}_candidates.vcf | awk '{print $1}')
if [ $lines -eq 3 ]
then
	echo "### No changes found ###" >> ${case}_candidates.vcf
	echo "No changes found" | tee -a pipeline.log
fi

END=$(date +%s)
duration=$((END-STARTTIME))
time=$(date +%x%t%X)
echo "${time}: Done! \n$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." | tee -a pipeline.log


