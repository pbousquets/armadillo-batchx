#!/usr/bin/env python3
import argparse
import sys
import gzip
import os
from pyfaidx import Fasta
from subprocess import check_output


def parse_args():
	parser = argparse.ArgumentParser(description = 'Prepares the data needed by Armadillo by providing just the coordinates of the regions of interest.')
	parser.add_argument(
	'-g', '--genome_ref', type = str, required = True,
	help = 'Reference genome')
	parser.add_argument(
	'-i', '--rois', type = str, required = True,
	help = 'Input file with regions of interest (BED or chr:st-end formatted)')
	parser.add_argument(
	'-I', '--identity', type = int, default = 90, required = False,
	help = 'Minimum identity for a hit to be considered a copy of the ROI (default: %(default)s).')
	parser.add_argument(
	'-L', '--lendiff', type = int, required = False, default = 15,
	help = 'Maximum length difference allowed between each hit and the input sequence (default: %(default)s%%).')
	parser.add_argument(
	'-m', '--mlen', type = int, default = 100, required = False,
	help = 'Minimum length (bp) allowed to each gene (default: %(default)s).')
	parser.add_argument(
	'-o', '--output', type = str, required = False, default = 'armadillo_data',
	help = 'Set name for the output directory (default: %(default)s).')
	parser.add_argument(
	'-p', '--port', type = str, required = True,
	help = 'Port where gfServer was loaded')
	return parser.parse_args()

def blat_parser(blat, filename):
	print_fasta = False
	printed = 0
	for line in blat:
		column=line.strip().split()
		chr_stend=column[0].split(":")
		st_end=chr_stend[1].split("-")
		end=st_end[1].split("_")
		querylength=int(end[0])-int(st_end[0])
		if float(column[2])>args.identity and float(column[3])/querylength*100 >= 100-args.lendiff and float(column[7])/querylength*100 <=100+args.lendiff and float(column[5])<=1: #90% identity, alignment length +/- 15%, gaps 1
			print_fasta = True
			printed += 1
			log=open("rois_copies_coords/"+filename, "a+") #We'll have a log file where we'll append any useful coordinate, so we don't use it twice
			if int(column[8])<int(column[9]): #We consider the start as the shortest coord, so if the it's in the negative strand, flip the coords
				log.write(column[1]+":"+column[8]+"-"+column[9]+"\n")
			else:
				log.write(column[1]+":"+column[9]+"-"+column[8]+"\n")
			log.close()
	if printed == 1:
		os.remove("rois_copies_coords/"+filename)
		print_fasta = False
	return(print_fasta)

args = parse_args()
href = Fasta(args.genome_ref, rebuild=False) #Open the reference genome

##PREPARE THE DIRECTORIES##
input_file = os.path.abspath(args.rois)
try:
	os.mkdir(args.output)
except FileExistsError:
	answer = input("Warning: The output directory already exists. Do you want to overwrite (y/n)? \n")
	if answer.lower() == "yes" or answer.lower() == "y":
		pass
	else:
		print("Exiting...")
		sys.exit(1)

os.chdir(os.getcwd()+"/"+args.output)

try:
	os.mkdir("rois_copies_coords")
except FileExistsError:
	pass

try:
	os.mkdir("miniFASTA")
except FileExistsError:
	pass

##READ THE INPUT ##
print("Performing blat queries to find repetitive regions of interest", file=sys.stderr)
if ".gz" in input_file:
	rois = gzip.open(input_file, "rt")
else:
	rois = open(input_file)

for line in rois:
	if line.startswith("#"): #Skip the header
		continue

	col = line.strip().split()

	if len(col)==1:
		col = line.replace(":", "\t").replace("-", "\t").strip().split() #convert chr:st-end format to BED

	try:
		chrom, start, end = col[0:3]
	except ValueError:
		sys.exit("Input error: expected a BED file or a chr:start-end formatted file")

	coord = str(chrom) + ":" + str(int(start) - 100) + "-" + str(int(end) + 100) #We add +100 bp to be sure that all reads align, even in those in the ends of the exons
	
	if int(end) - int(start) < args.mlen: #Remove too short genes
		continue

	if "miniFASTA/"+coord+".fa" not in os.listdir("miniFASTA"): #Do not analyse those that already exist
		orig_coord = str(chrom) + ":" + str(start) + "-" + str(end)
		blat_input = ">"+orig_coord + "\n" + href[chrom][int(start)-1:int(end)].seq #Get the sequence of the region to run blat. We substract 1 because 0 is the start reference for pyfasta
		blat_command = ['gfClient', '-out=blast8','localhost', args.port , '', 'stdin', 'stdout']
		blat_result = check_output(blat_command, input=blat_input.encode()).decode().strip().split("\n")
		if len(blat_result) > 1:
			print_fasta = blat_parser(blat_result, orig_coord)
			if print_fasta:
				miniFASTA = open("miniFASTA/" + orig_coord + ".fa", "w+") #Create a FASTA for each region
				miniFASTA.write(">" + coord + "\n" + href[chrom][int(start)-101:int(end)+100].seq+"\n") #We'll align the reads against a quite larger region so that the ones that overlap only in the flanks can align too
				miniFASTA.close()
		else:
			pass
	else:
		continue

## Remove duplicates
files = [file for file in os.listdir("rois_copies_coords")]
for file in files:
	for line in open("rois_copies_coords/"+file):
		line = line.strip()
		if line != file and line in files and file in files:
			os.remove("rois_copies_coords/"+line)
			os.remove("miniFASTA/"+line+".fa")
			files.remove(line)
			print("Removed ", line, ". Duplicate of: ", file, sep = "", end = "\n")
		else:
			continue

## Index the fasta files ##
fastas = [file for file in os.listdir("miniFASTA") if file.endswith(".fa")]
concat = ''.join([open("miniFASTA/"+f).read() for f in fastas])
refFASTA = open("armadillo_reference_genome.fa", "w+") #We'll write a merged reference genome
refFASTA.write(concat)
refFASTA.close()
os.system("samtools faidx armadillo_reference_genome.fa")

for fasta in fastas:
	if os.path.isfile("miniFASTA/"+fasta+".fai"): #Don't reindex if already done
		continue
	else:
		os.system("bwa index miniFASTA/" + fasta)
		os.system("samtools faidx miniFASTA/" + fasta)

## Create a list of all final rois
roisfile = open("rois", "w+")
for roi in [file for file in os.listdir("rois_copies_coords")]:
	roisfile.write(roi+"\n")
roisfile.close()
