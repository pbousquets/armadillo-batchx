#!/usr/bin/env python3
import argparse
import sys
import gzip
import os
from pyfaidx import Fasta
from subprocess import check_output

def parse_args():
        parser = argparse.ArgumentParser(description = 'Prepares the data needed by Armadillo by providing just the coordinates of the regions of interest in BED format.')
        parser.add_argument(
        '-g', '--genome_ref', type = str, required = True, metavar = 'FASTA',
        help = 'Reference genome')
        parser.add_argument(
        '-i', '--rois', type = str, required = True, metavar = 'LIST',
        help = 'Input file with regions of interest (BED-formatted)')
        parser.add_argument(
        '-p', '--port', type = str, required = True, metavar = 'INT',
        help = 'Port where gfServer was loaded (BLAT)')
        parser.add_argument(
        '-o', '--output', type = str, required = False, metavar = 'STR', default = 'armadillo',
        help = 'Set name. It will be used for output dir. (default: %(default)s)')
        return parser.parse_args()

def blat_parser(blat, filename):
    for line in blat:
        column=line.strip().split()
        chr_stend=column[0].split(":")
        st_end=chr_stend[1].split("-")
        end=st_end[1].split("_")
        querylength=int(end[0])-int(st_end[0])
        if float(column[2])>90 and float(column[3])/querylength*100 >= 85 and float(column[7])/querylength*100 <=115 and float(column[5])<=1: #90% identity, alignment length +/- 2%, gaps 1
            log=open("rois_copies_coords/"+filename, "a+") #We'll have a log file where we'll append any useful coordinate, so we don't use it twice
            if int(column[8])<int(column[9]): #We consider the start as the shortest coord, so if the it's in the negative strand, flip the coords
                log.write(column[1]+":"+column[8]+"-"+column[9]+"\n")
            else:
                log.write(column[1]+":"+column[9]+"-"+column[8]+"\n")
            log.close()

args = parse_args()
href = Fasta(args.genome_ref, rebuild=False) #Open the reference genome

##PREPARE THE DIRECTORIES##
os.mkdir
refFASTA=open("armadillo_reference_genome.fa", "w+") #We'll write a new reference genome

try:
    os.mkdir(args.output + "_data")
except FileExistsError:
    print("Warning: armadillo_data directory already exists in this path.")

os.chdir(os.getcwd()+"/armadillo_data")

try:
    os.mkdir("rois_copies_coords")
except FileExistsError:
    pass

try:
    os.mkdir("miniFASTA")
except FileExistsError:
    pass

##READ THE INPUT ##
if ".gz" in args.rois:
        rois = gzip.open(args.rois, "rt")
else:
        rois = open(args.rois)

for line in rois:
        if line.startswith("#"): #Skip the header
                continue
        col = line.strip().split()
        try:
                chrom, start, end = col[0:3]
        except ValueError:
                sys.exit("Input error: expected a three column file (chrosome, start, end). Input has " + str(len(col))+ " columns")

        coord = str(chrom) + ":" + str(int(start) - 100) + "-" + str(int(end) + 100) #We add +100 bp to be sure that all reads align, even in those in the ends of the exons
        orig_coord = str(chrom) + ":" + str(start) + "-" + str(end)

        miniFASTA = open("miniFASTA/" + coord + ".fa", "w+") #Create a FASTA for each region
        miniFASTA.write(">" + coord + "\n" + href[chrom][int(start)-100:int(end)+100].seq+"\n") #We'll align the reads against a quite larger region so that the ones that overlap only in the flanks can align too
        refFASTA.write(">" + coord + "\n" + href[chrom][int(start)-100:int(end)+100].seq+"\n")
        blat_input = ">"+orig_coord     + "\n" + href[chrom][int(start):int(end)].seq #Get the sequence of the region to run blat
        blat_command = ['gfClient', '-out=blast8','localhost', args.port , '', 'stdin', 'stdout']
        blat_result = check_output(blat_command, input=blat_input.encode()).decode().strip().split("\n")
        if len(blat_result) > 1:
                blat_parser(blat_result, coord)
        miniFASTA.close()

##Index the fasta files##
fastas = os.listdir("miniFASTA")
for fasta in fastas:
	if any(substring in fasta for substring in [".sa", ".amb", ".ann", ".pac", ".bwt", "fai"]): #Don't try to index the indexes if exist
		pass
	else:
		os.system("bwa index miniFASTA/" + fasta)
		os.system("samtools faidx miniFASTA/" + fasta)
os.system("samtools faidx " + args.output)
