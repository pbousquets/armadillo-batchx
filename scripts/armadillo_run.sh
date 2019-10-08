#!/bin/sh
## RUN ARMADILLO ##
STARTTIME=$(date +%s)

##Display help if no arguments provided
if [ $# -eq 0 ]
then
	usage
	exit 1
fi

##Load parameters##
if [ -f $1 ] #Arguments or configuration file
then
	. $(dirname $1)/$(basename $1) # Read a config file
	shift
else
	parse_arguments $@ #Parse the arguments pass through the command line
fi

##Load data##
ref_genome=${armadillo_data}'/armadillo_reference_genome.fa'
blat_coords=${armadillo_data}'/rois_copies_coords'
miniFasta_dir=${armadillo_data}'/miniFASTA'
repeatmasker_candidates_filter=${scripts_dir}'/repeatmasker_candidates_filter.py'
mq2vcf_TDvsND=${scripts_dir}'/mq2vcf_TDvsND.py'

if [ "$print" = 'true' ]
then
	printopt='-f'
fi

##Check required variables
if [ "$case" = '' ] || [ "$tumor_genome" = '' ] || [ "$control_genome" = '' ]
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
ND_minibam="${case}control_merged.bam"
TD_minibam="${case}tumor_merged.bam"

##Check if files exist
evaluate $ref_genome Reference genome
evaluate $repeatsDB repeatsDB

##Run pipeline
echo "Final options:   \n  Input: $case   \n  Root dir: $root_dir   \n  Control genome: $control_genome   \n  Tumor genome: $tumor_genome   \n  List of ROIs: $rois_list   \n  Blat coordinates directory: $blat_coords   \n  miniFasta directory: $miniFasta_dir   \n  RepeatsDB: $repeatsDB   \n  Reference genome: $ref_genome   \n  Scripts directory: $scripts_dir   \n  Min. coverage: $mincov   \n  Control max. mut reads: $control_cutoff   \n  Sequencing error rate: $seq_error   \n  Control contamination (%): $control_contamination    \n  Tumor cutoff: $tumor_cutoff   \n  Control quality: $control_qual   \n  Base quality: $tum_qual   \n  Mapping quality cutoff: $mapq   \n  Read length: $read_length   \n  Max errors: $max_errors   \n  GC content cutoff: $gc_content   \n  Threads: $threads   \n  Print: $print   \n  Skip: $skip   \n  Port: $port \n" | tee -a pipeline.log

if [ ${skip} = 'false' ]
then
	##Check if needed files exist##
	evaluate $rois_list ROIs genome
	evaluate $TD Tumor genome
	evaluate $ND Control_genome

	##Check if it was already analysed##
	if [ -f $TD_minibam ] || [ -f $ND_minibam ]
	then
		echo "FileExist ERROR:\nMake sure that $TD_minibam or $ND_minibam don't exist already. Else, run armadillo with --skip true."
		exit 0
	fi

	##Minibam extraction step##
	extract_minibam $case ${TD} tumor ${blat_coords} ${miniFasta_dir} ${rois_list} ${threads} | tee -a pipeline.log
	extract_minibam $case ${ND} control ${blat_coords} ${miniFasta_dir} ${rois_list} ${threads} | tee -a pipeline.log
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

if [ -f discarded_variants.log ]
then
	rm discarded_variants.log
fi

##Filter step##
samtools mpileup --output-QNAME -Q ${control_qual} -q ${mapq} -R -f ${ref_genome} ${TD_minibam} ${ND_minibam} | python3 ${mq2vcf_TDvsND} -i - -tb ${TD_minibam} -cb ${ND_minibam} -n ${case} -r ${ref_genome} -tt ${tumor_cutoff} -nt ${control_contamination} -se ${seq_error} -nm ${control_cutoff} -rl ${read_length} -e ${max_errors} -gc ${gc_content} -q ${tum_qual} -c ${mincov} -t ${threads} -p ${port} ${printopt} | python3 ${repeatmasker_candidates_filter} $repeatsDB 20 100 ${printopt}  > ${case}_candidates.vcf #The 20 specifies the max percentage of reads of a mutation that can appear in more mutations. The 100 is the length of the flanking regions added during the data preparation. By default is 100.

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
