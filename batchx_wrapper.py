#!/usr/bin/python3

"""
Parse batchx user json for launching armadillo
"""
import json
import os
import subprocess
import tarfile
import time
from sys import exit


## Set environment
bxMemory = int(int(os.environ['BX_MEMORY'])/1000) # To gb
bxVcpus = int(os.environ['BX_VCPUS'])
port = 9007
tmpdir = "/tmp/"
defaults = {
        "controlCoverage": 30,
        "tumorCoverage": 30,
        "controlThreshold": 3,
        "tumorThreshold": 6,
        "baseQuality": 30,
        "mapq": 30,
        "GCcutoff": 80,
        "returnBams": False,
        "returnDiscarded": False}
        
## Parse json
with open("/batchx/input/input.json", "r") as inputFile:
    inputJson = inputFile.read()
parsedJson = json.loads(inputJson)
parsedJson = {**defaults, **parsedJson}
project_name = parsedJson["name"]

assert (parsedJson["armadilloData"].endswith(".gz")), "Armadillo data-prep must be compressed as a gz file."

## Prepare bams
bamdir = '/tmp/bams/'
parsedJson["bamDir"] = bamdir
os.mkdir(bamdir)

os.symlink(parsedJson["controlGenome"], bamdir+os.path.basename(parsedJson["controlGenome"]))
os.symlink(parsedJson["tumorGenome"], bamdir+os.path.basename(parsedJson["tumorGenome"]))

parsedJson["controlGenome"] = os.path.basename(parsedJson["controlGenome"])
parsedJson["tumorGenome"] = os.path.basename(parsedJson["tumorGenome"])

assert (parsedJson["controlGenome"].lower().endswith(".bam")), "Control genome is expected to be a BAM file. Are you sure it is indeed BAM formated/named?"
assert (parsedJson["tumorGenome"].lower().endswith(".bam")), "Tumor genome is expected to be a BAM file. Are you sure it is indeed BAM formated/named?"

if "tumorIndex" not in parsedJson:
    print("Indexing tumor bam", flush = True)
    p1 = subprocess.Popen(f'samtools index -@ {bxVcpus/2} {bamdir+os.path.basename(parsedJson["tumorGenome"])}', shell = True)
else:
    os.symlink(parsedJson["tumorIndex"], bamdir + os.path.basename(parsedJson["tumorIndex"]))

if "controlIndex" not in parsedJson:
    print("Indexing control bam", flush = True)
    p2 = subprocess.Popen(f'samtools index -@ {bxVcpus/2} {bamdir+os.path.basename(parsedJson["controlGenome"])}', shell = True)
else:
    os.symlink(parsedJson["controlIndex"], bamdir + os.path.basename(parsedJson["controlIndex"]))


## Prepare ref
custom_rois = True if "roisList" in parsedJson else False

print("Extracting armadillo_data gz file", flush = True)
tar = tarfile.open(parsedJson['armadilloData'])
tar.extractall('armadilloData')
tar.close()
subdir = os.listdir("armadilloData")

if len(subdir) == 1: #Check if armadillodata prep was compressed as one single folder or not
    parsedJson['armadilloData'] = os.path.abspath(f'armadilloData/{subdir[0]}')
else:
    parsedJson['armadilloData'] = os.path.abspath('armadilloData')

parsedJson['roisList'] = parsedJson['roisList'] if custom_rois else parsedJson['armadilloData']+'/rois'

run_args = ["name", "controlGenome", "tumorGenome", "roisList", "bamDir", "model", "armadilloData"]
cmd_run = ''
for attribute, value in parsedJson.items():
    cmd_run += f'--{attribute} {value} ' if attribute in run_args else ''

params = {key: parsedJson[key] for key in ["controlCoverage", "tumorCoverage", "controlThreshold", "tumorThreshold", "baseQuality", "mapq", "GCcutoff"]}

for attribute, value in params.items():
    cmd_run += f'--{attribute} {value} '

cmd_final = f'armadillo run {cmd_run} --threads {bxVcpus} --maxRam {bxMemory} --skip false --print false'

try: # Wait for bam indexes if needed
    p1.wait()
    p2.wait()
except NameError:
    pass

## Prepare output
outputDir = "/batchx/output/" 
outputJson = {}

try:
    r = subprocess.run(cmd_final, shell = True, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    if r.returncode == 0:
        print("The run was successful! Closing the environment", flush = True)
    else:
        print("Error running Armadillo. Exiting now" , flush = True)
        print(r.stdout)
        print(r.stderr.decode())
        exit(1)

    outputJson["vcf"] = f'{outputDir}/{project_name}/{project_name}_final.vcf' 
    
    if parsedJson["returnBams"]:
        outputJson["tumorMinibam"] = f'{outputDir}/{project_name}/{project_name}tumor_merged.bam' 
        outputJson["tumorMinibamIdx"] = f'{outputDir}/{project_name}/{project_name}tumor_merged.bam.bai' 
        outputJson["controlMinibam"] = f'{outputDir}/{project_name}/{project_name}control_merged.bam' 
        outputJson["controlMinibamIdx"] = f'{outputDir}/{project_name}/{project_name}control_merged.bam.bai' 

    if parsedJson["returnDiscarded"]:
        outputJson["discardedVariants"] = f'{outputDir}/{project_name}/discarded_variants.log' 

except subprocess.CalledProcessError as e:
        print(e)
        exit(e.returncode)

with open('/batchx/output/output.json', 'w+') as json_file:
    json.dump(outputJson, json_file)
