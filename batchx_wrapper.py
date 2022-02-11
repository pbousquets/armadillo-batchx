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
        "identity": 95,
        "lendiff": 15,
        "mlen": 100,
        "returnBams": False,
        "returnDiscarded": False,
        "returnDataprep": False}

def compress_dir(source_dir, name):
    """
    Compress directory into a tar.gz
    """
    with tarfile.open(name, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

def validateBlatPorts(check_port = port, timeout = 120):
    """
    Run gfServer status in order to know if the started port is actually ready. Else, keep trying until timeout (seconds).
    """
    cmd_status = f'gfServer status localhost {check_port}'
    start = time.time()
    print("Checking ports. This may take a minute.", flush = True)
    wait = True
    while wait:
        wait = True if subprocess.run(cmd_status, shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode != 0 else False 
        if time.time() - start > timeout:
            print(f'Blat server setting up timeout. Server could not be initialized after {timeout} seconds', flush = True)
            exit(1)

## Parse json
with open("/batchx/input/input.json", "r") as inputFile:
    inputJson = inputFile.read()
parsedJson = json.loads(inputJson)
parsedJson = {**defaults, **parsedJson}
project_name = parsedJson["name"]

## Start BLAT
twobit = parsedJson["referenceBlat"]
os.symlink(parsedJson["referenceBlat"], "/" + os.path.basename(parsedJson["referenceBlat"]))
subprocess.Popen(f'gfServer start localhost {port} {twobit}', shell=True) # Use Popen so it runs on background

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


## Check data-prep 
runDataPrep = False if "armadilloData" in parsedJson else True
custom_rois = True if "armadilloData" in parsedJson and "roisList" in parsedJson else False

if runDataPrep:
    assert ("genomeRef" in parsedJson), "Missing genome reference fasta. Needed for running data-prep module. Please, include this file or a dataprep.gz file"
    assert ("roisList" in parsedJson), "Missing rois list. Needed for running data-prep module. Please, include this file or a dataprep.gz file with dataprep"
    fa = parsedJson["genomeRef"]
    destfa =  tmpdir + os.path.basename(fa)
    os.symlink(fa, destfa)

    if "genomeRefIndex" in parsedJson:
        fai = parsedJson["genomeRefIndex"]
        os.symlink(fai, tmpdir + os.path.basename(fai))
    else:
        print("Indexing reference genome", flush = True)
        subprocess.check_call(f'samtools faidx {destfa}', shell=True)

    parsedJson["rois"] = parsedJson["roisList"]

    dataprep_params = {key: parsedJson[key] for key in ["identity", "lendiff", "mlen"]}
    parsedJson["genomeRef"] = destfa

    dataprep_args = ["genomeRef", "rois", "identity", "lendiff", "mlen"]

    cmd_dataprep = f'armadillo data-prep --port {port} --output armadillo_data --rois {parsedJson["rois"]} '
    for attribute, value in dataprep_params.items():
        cmd_dataprep += f'--{attribute} {value} ' if attribute in dataprep_args else ''

    validateBlatPorts(port)
    print("Running armadillo-data prep", "Command:", cmd_dataprep, sep = "\n", flush = True)
    subprocess.run(cmd_dataprep, shell = True)
else:
    print("Extracting armadillo_data gz file", flush = True)
    tar = tarfile.open(parsedJson['armadilloData'])
    tar.extractall('.')
    tar.close()

parsedJson['armadilloData'] = os.path.abspath("armadillo_data")


## Select rois
parsedJson['roisList'] = parsedJson['roisList'] if custom_rois else os.path.abspath("armadillo_data") + "/rois"

run_args = ["name", "controlGenome", "tumorGenome", "roisList", "bamDir", "model", "armadilloData"]
cmd_run = ''
for attribute, value in parsedJson.items():
    cmd_run += f'--{attribute} {value} ' if attribute in run_args else ''

params = {key: parsedJson[key] for key in ["controlCoverage", "tumorCoverage", "controlThreshold", "tumorThreshold", "baseQuality", "mapq", "GCcutoff"]}

for attribute, value in params.items():
    cmd_run += f'--{attribute} {value} '

cmd_final = f'armadillo run {cmd_run} --threads {bxVcpus} --maxRam {bxMemory}  --port {port} --skip false --print false'

validateBlatPorts(port)

try: # Wait for bam indexes if needed
    p1.wait()
    p2.wait()
except NameError:
    pass

print("Blat ports are ready. Environment is also ready.", "Launching Armadillo with command:", cmd_final, sep = "\n", flush = True)

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

    if parsedJson["returnDataprep"]:
        compress_dir(os.path.abspath("armadilloData"), "/batchx/output/armadillo_data.tar.gz")
        outputJson["dataprep"] = "/batchx/output/armadillo_data.tar.gz"

    if parsedJson["returnDiscarded"]:
        outputJson["discardedVariants"] = f'{outputDir}/{project_name}/discarded_variants.log' 

except subprocess.CalledProcessError as e:
        print(e)
        exit(e.returncode)

with open('/batchx/output/output.json', 'w+') as json_file:
    json.dump(outputJson, json_file)
