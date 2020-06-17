# Armadillo

Mutations in repetitive regions are usually lost by standard variant callers as low mapping quality leads to non-reliable variant calls. We introduce Armadillo, a pipeline that, by accepting not to know the exact affected copy, allows us to find mutations in repetitive genes.

## Getting Started

### Prerequisites

* [BWA](http://bio-bwa.sourceforge.net/) (v.0.7.17)
* [Samtools](http://www.htslib.org/doc/samtools.html) (v.1.9)
* [gfClient](https://genome.ucsc.edu/goldenPath/help/blatSpec.html#gfClientUsage)  (v.35)
* [Python3](https://www.python.org) (v.3.7)



### Installing

Use git to download Armadillo:

```
git clone https://github.com/pbousquets/armadillo
```

Check that the dependencies are installed:

```
cd armadillo

./configure #run it as sudo to install the python packages and creating an alias globally 
```

## Getting armadillo ready

### 1. Open a port for blat 
Before we can run armadillo, we need to get a port prepared to run gfClient. In order to do that, we just run:

```
gfServer start localhost $PORT /path/to/reference_genome.2bit &
```
The port will stay opened unless we kill the task or shut the computer down.

### 2. Get the custom armadillo_data for your ROIs
The regions of interest must be analysed before running armadillo to keep just those which are repetitive and get the coords of their copies. Armadillo can do that just by providing a reference genome, a BED-formatted list of regions of interest and the port previously opened for gfClient:

```
armadillo data-prep -i /path/to/rois.bed -g /path/to/reference_genome -p port [ -m min_len  -o output_dir ]
```

### Running armadillo

It's possible to use a configuration file to run armadillo. Just by using config-file option a configuration file will be copied at current working directory. Then, just change any parameter you want and run armadillo.

```
armadillo config-file
armadillo run configuration_file.txt
```
Options can be also passed directly through the command line:

```
armadillo config-file
armadillo run -n CASE -C control.bam -T tumor.bam --armadillo_data /path/to/armadillo_data [options]
```
__Important consideration before running armadillo:__

When analysing the ROIs during the "data-prep" step, the gfServer must used a port with the **same reference genome** used to align the genomes that will be provided later to armadillo. However, when using "armadillo run", gfServer is used just to check if the reads align perfectly anywhere in the genome so we can discard these reads. Thus, **we can use the latest reference genome** in that step, even though if the alignment was performed with a previous reference genome. It will allow us remove false positives associated to regions that actually don't exist in previous versions.

### Output

The program will print multiple files. We'll get two minibams that are generated in the first step of the pipeline. 

Also, we'll find two VCF files. CASE_candidates.vcf is the first one to be printed. It's an intermediate with mutations' readnames, which are then used by remove_dups.py to remove duplications (same mutation appear in multiple similar regions of interest). This will script print the final CASE_nodupscandidates.vcf file.  

## Authors

* **Pablo Bousquets-Muñoz** - bousquetspablo@uniovi.es
* **Ander Díaz-Navarro**
* **Xose Antón Suárez-Puente**
