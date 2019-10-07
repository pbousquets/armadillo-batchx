# Armadillo

Armadillo is a pipeline design to recover mutations that belong to repetitive regions, which usually are discarded due to low confidence on the alignments.

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

bash check_dependencies.bash (run it as sudo to install the python packages globally)
```

## Getting armadillo ready

Before we can run armadillo, we need to get a port prepared to run gfClient. In order to do that, we just run:

```
gfServer start localhost 9008 /path/to/reference_genome.2bit &
```
The port will stay opened unless we kill the task or shut down the computer.

Also, the regions of interest must be analysed before running armadillo, to keep just those which are repetitive and get the coords of their copies. Armadillo can do that just by providing a reference genome, a BED-formatted list of regions of interest and the port previously opened for gfClient:

```

./armadillo data-prep -i /path/to/rois.bed -g /path/to/reference_genome -p port -o output_dir
```

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

## Built With

* [BWA](http://bio-bwa.sourceforge.net/) - Genome aligner
* [Samtools](http://www.htslib.org/doc/samtools.html) - Bam files' management
* [Python3](https://www.python.org) - Data filtering
* [gfClient](https://genome.ucsc.edu/goldenPath/help/blatSpec.html#gfClientUsage) - Blat

## Authors

* **Pablo Bousquets Muñoz**
* **Ander Díaz-Navarro**
* **Xose Antón Suárez Puente**
