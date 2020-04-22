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
    parse_arguments $@
else
	parse_arguments $@ #Parse the arguments pass through the command line
fi

##Load data##
ref_genome=${armadillo_data}'/armadillo_reference_genome.fa'
blat_coords=${armadillo_data}'/rois_copies_coords'
miniFasta_dir=${armadillo_data}'/miniFASTA'
remove_dups=${scripts_dir}'/remove_dups.py'
mq2vcf_TDvsND=${scripts_dir}'/mq2vcf_TDvsND.py'
extract_minibam=${scripts_dir}'/extract_minibam.sh'

if [ "$print" = 'true' ]
then
	printopt='-f'
fi

##Check required variables
if [ "$name" = '' ] || [ "$tumor_genome" = '' ] || [ "$control_genome" = '' ]
then
	echo "ERROR: No input provided."
	usage
	exit
fi

TD=${root_dir}${tumor_genome}
ND=${root_dir}${control_genome}
ND_minibam="${name}control_merged.bam"
TD_minibam="${name}tumor_merged.bam"

##Check if files exist
evaluate $ref_genome Reference genome
evaluate $repeatsDB repeatsDB

if [ ${skip} = 'false' ]
then
	if [ ! -d ${name} ]
	then
		mkdir -p ${name}/tumor_tmp_files
		mkdir ${name}/control_tmp_files
		cd ${name}

		##Run pipeline##
        echo " Final command:\n armadillo run --name ${name} --root_dir $root_dir --control_genome $control_genome --tumor_genome $tumor_genome --rois_list $rois_list --armadillo_data $armadillo_data --repeatsDB $repeatsDB --scripts_dir $scripts_dir --control_coverage $control_coverage --tumor_coverage $tumor_coverage --control_max $control_max --tumor_threshold $tumor_threshold --base_quality $base_quality --control_qual $control_qual --mapq $mapq --read_length $read_length --GCcutoff $GCcutoff --threads $threads --skip $skip --port $port --print $print \n" | tee -a pipeline.log
        echo " Final options:   \n  Sample name: ${name}   \n  Genomes dir: $root_dir   \n  Control genome: $control_genome   \n  Tumor genome: $tumor_genome   \n  ROIs list: $rois_list   \n  Armadillo data path: $armadillo_data   \n  Repeats database: $repeatsDB   \n  Scripts directory: $scripts_dir   \n  Control coverage: $control_coverage   \n  Tumor coverage: $tumor_coverage    \n  Control maximum mutant reads: $control_max   \n  Tumor minimum mutant reads: $tumor_threshold   \n  Tumor base quality threshold: $base_quality   \n  Control base quality threshold: $control_qual   \n  Mapping quality threshold: $mapq   \n  Read length: $read_length   \n  GC content maximum: $GCcutoff   \n  Threads: $threads   \n  Skip: $skip   \n  Port: $port   \n  Print: $print" | tee -a pipeline.log

	else
		echo "That case already exists. If you want to reanalyse it, please, use '--skip true'"
		exit 0
	fi

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

	bash $extract_minibam $name ${TD} tumor ${blat_coords} ${miniFasta_dir} ${rois_list} ${threads}
	bash $extract_minibam $name ${ND} control ${blat_coords} ${miniFasta_dir} ${rois_list} ${threads}
	
else
	if [ -d ${name} ]
	then
		cd $name
        echo " Final command:\n armadillo run --name ${name} --root_dir $root_dir --control_genome $control_genome --tumor_genome $tumor_genome --rois_list $rois_list --armadillo_data $armadillo_data --repeatsDB $repeatsDB --scripts_dir $scripts_dir --control_coverage $control_coverage --tumor_coverage $tumor_coverage --control_max $control_max --tumor_threshold $tumor_threshold --base_quality $base_quality --control_qual $control_qual --mapq $mapq --read_length $read_length --GCcutoff $GCcutoff --threads $threads --skip $skip --port $port --print $print \n" | tee -a pipeline.log
        echo " Final options:   \n  Sample name: ${name}   \n  Genomes dir: $root_dir   \n  Control genome: $control_genome   \n  Tumor genome: $tumor_genome   \n  ROIs list: $rois_list   \n  Armadillo data path: $armadillo_data   \n  Repeats database: $repeatsDB   \n  Scripts directory: $scripts_dir   \n  Control coverage: $control_coverage   \n  Tumor coverage: $tumor_coverage    \n  Control maximum mutant reads: $control_max   \n  Tumor minimum mutant reads: $tumor_threshold   \n  Tumor base quality threshold: $base_quality   \n  Control base quality threshold: $control_qual   \n  Mapping quality threshold: $mapq   \n  Read length: $read_length   \n  GC content maximum: $GCcutoff   \n  Threads: $threads   \n  Skip: $skip   \n  Port: $port   \n  Print: $print" | tee -a pipeline.log
		echo "Skipped minibam extraction." | tee -a pipeline.log
	else
		echo "$name doesn't seem to exist. Please, verify the it exists in your current directory or use '--skip false'"
		exit 0
	fi
fi

time=$(date +%x%t%X)
echo ${time}: Finding candidates...

if [ -f ${name}_candidates.vcf ]
then
	rm ${name}_candidates.vcf
fi

if [ -f discarded_variants.log ]
then
	rm discarded_variants.log
fi

##Filter step##
samtools mpileup --output-QNAME -Q ${base_quality} -q ${mapq} -R -f ${ref_genome} ${TD_minibam} ${ND_minibam} | python3 ${mq2vcf_TDvsND} -i - -tb ${TD_minibam} -cb ${ND_minibam} -tc ${tumor_coverage} -n ${name} -r ${ref_genome} -tt ${tumor_threshold} -cm ${control_max} -rl ${read_length} -gc ${GCcutoff} -q ${control_qual} -cc ${control_coverage} -t ${threads} -p ${port} ${printopt} > ${name}_candidates.vcf #The 20 specifies the max percentage of reads of a mutation that can appear in more mutations. The 100 is the length of the flanking regions added during the data preparation. By default is 100.
cat ${name}_candidates.vcf | python3 ${remove_dups} 100 > ${name}_nodupscandidates.vcf
lines=$(wc -l ${name}_candidates.vcf | awk '{print $1}')
if [ $lines -eq 3 ]
then
	echo "### No changes found ###" >> ${name}_candidates.vcf
	echo "No changes found" | tee -a pipeline.log
fi

END=$(date +%s)
duration=$((END-STARTTIME))
time=$(date +%x%t%X)
echo "${time}: Done! \n$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." | tee -a pipeline.log
