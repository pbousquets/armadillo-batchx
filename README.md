# Armadillo

Mutations in repetitive regions are usually lost when analysing NGS data, as the alignment of repetitive regions are seldom trustful. Therefore, Armadillo is a pipeline designed to realign specifically those regions and recover mutations that were invisible for regular variant callers.


## Getting Started

### Prerequisites

```
bwa (>= v.0.7.17)
gfClient (>= v.35)
gfServer (>= v.35)
python3 (>= v.3.7)
samtools (>= v.1.9)

Python will also use: argparse, gzip, multiprocessing, os, pyfaidx, re, sys, statistics and subprocess packages.
```

### Installing

Use git to download Armadillo:

```
git clone ssh://git@156.35.56.116:10022/pbousquets/armadillo.git
```

Check that the dependencies are installed:

```
cd armadillo

bash check_dependencies.bash #run it as sudo to install the python packages globally
```

## Getting armadillo ready

Before we can run armadillo, we need to get a port prepared to run gfClient. In order to do that, we just run:

```
gfServer start localhost 9008 /path/to/reference_genome.2bit &
```
The port will stay opened unless we kill the task or shut the computer down.

Also, the regions of interest must be analysed before running armadillo to keep just those which are repetitive and get the coords of their copies. Armadillo can do that just by providing a reference genome, a BED-formatted list of regions of interest and the port previously opened for gfClient:

```
./armadillo data-prep -i /path/to/rois.bed -g /path/to/reference_genome -p port -o output_dir
```

### Running armadillo

It's possible to use a configuration file to run armadillo. Just by using config-file option a configuration file will be copied at current working directory. Then, just change any parameter you want and run armadillo.

```
./armadillo config-file
./armadillo run configuration_file.txt
```
Options can be also passed directly through the command line:

```
./armadillo config-file
./armadillo run -i ID -C control.bam -T tumor.bam [options]
```
__Important consideration before running armadillo:__
When analysing the ROIs during the "data-prep" step, the gfServer must used a port with the **same reference genome** used to align the genomes that will be provided later to armadillo. However, when using "armadillo run", gfServer is used just to check if the read aligns perfectly anywhere in the genome so we can discard these reads. Thus, **we can use the latest reference genome** in that step, even though the alignment was performed on a previous version, as it's more complete. It will allow us remove more false positives associated to regions that actually don't exist in the reference genome used to align the genomes.
 
## Built with:

* [BWA](http://bio-bwa.sourceforge.net/) - Genome aligner
* [Samtools](http://www.htslib.org/doc/samtools.html) - Bam files' management
* [Python3](https://www.python.org) - Data filtering
* [gfClient](https://genome.ucsc.edu/goldenPath/help/blatSpec.html#gfClientUsage) - Blat

## Authors

* **Pablo Bousquets Muñoz**
* **Ander Díaz-Navarro**
* **Xose Antón Suárez Puente**
