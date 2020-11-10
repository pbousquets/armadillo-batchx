#!/bin/bash
time=$(date +%x%t%X)
echo ${time}: Generating ${2} minibams... | tee -a pipeline.log
name=$1
sample=$2
type=$3
blat_coords=$4 #Dir of blat coords of each exon
miniFasta_dir=$5
rois_list=$6 #List of regions of interest
threads=$7
chr_bam=$8
maxRam=$9

lines=$(cat $rois_list | wc -l )
cat $rois_list | tr ':-' '\t' | while read -r chr st end #Create a tmp minibam for each region
do
    file=$(echo ${chr}:${st}-${end})

    if [ -f ${blat_coords}/${file} ]
    then
        lines=$((lines - 1))
        echo -en "\r Files left: $lines                                                "
        #Header
        SM=${sample}
        PL=illumina
        LB=WGS #Or WES
        PU=${sample}
        RG="@RG\\tID:${sample}\\tSM:${SM}\\tPL:${PL}\\tLB:${LB}\\tPU:${PU}"

        if [ "$chr_bam" = "False" ]
        then
            allCoords=$(sed -e 's/chr//g' ${blat_coords}/${file}) 
        else
            allCoords=$(sed -e 's/chr//g' ${blat_coords}/${file} | awk '{print "chr"$0}')
        fi
        
        samtools view -@ 3 -u -f 1 -F 3072 -G 0x400 -G 0x4 ${sample} $allCoords 2>> pipeline.log | samtools fastq -@ 2 -N - 2>> pipeline.log |  bwa mem -t 3 -R ${RG} ${miniFasta_dir}/${file}.fa - 2>/dev/null | samtools view -uS - > ${type}_tmp_files/${file}.bam  &
        
        subp=$(($(jobs -r -p | wc -l)+1))
        if [ "$subp" -gt "$((threads - 1))" ]
        then
            wait -n
        fi

    else
        echo "~ ${blat_coords}/${file} does not exist. Make sure it was in the armadillo data-prep list" >> pipeline.log
    fi

    if [ "$lines" -eq 0 ]
    then
        wait
    fi
done

time=$(date +%x%t%X)
echo -e "\n"${time}: Merging ${name} minibams. This may take a while | tee -a pipeline.log

iter=0
cd ${type}_tmp_files

while [ $(ls | grep -v tmp | wc -l) -ne 0 ]
do
    iter=$((iter+1))
    bams=$(ls | grep -v tmp | head -n 1000)
    samtools merge -@ ${threads} tmp_${iter}.bam ${bams}
    rm ${bams}
done

if [ ${iter} -gt 1 ]
then
    samtools merge -@ ${threads} - tmp_*.bam | samtools sort -@ ${threads} -m ${maxRam}00000000 -o ../${name}${type}_merged.bam -
    rm tmp_*.bam
else
    samtools sort -@ ${threads} -m ${maxRam}00000000 -o ../${name}${type}_merged.bam tmp_${iter}.bam
    rm tmp_${iter}.bam
fi
cd ../
rm -d ${type}_tmp_files
samtools index ${name}${type}_merged.bam

time=$(date +%x%t%X)
echo -e ${time}: Bam generation successful "\n"| tee -a pipeline.log
