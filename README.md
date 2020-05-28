# Armadillo

<style>
body {
text-align: justify}
</style>


Mutations in repetitive regions are usually lost by standard variant callers as low mapping quality leads to non-reliable variant calls. We introduce Armadillo, a pipeline that, by accepting to lose the exact position of mutations in the genome, allows us to find mutations in repetitive genes.

## Getting Started

### Prerequisites

```
bwa (>= v.0.7.17)
gfClient (>= v.35)
gfServer (>= v.35)
python3 (>= v.3.7)
samtools (>= v.1.9)

The next python packages are also required: 
argparse, multiprocessing, numpy, pandas, pyfaidx, pysam, re, scipy, statistics, subprocess, sys, tabix and os. If any of them is missing, it'll be automatically installed when configuring armadillo.
```

### Installing

Use git to download Armadillo:

```
git clone ssh://git@156.35.56.116:10022/pbousquets/armadillo.git
```

Check that the dependencies are installed:

```
cd armadillo

./configure #run it as sudo to install the python packages and creating an alias globally 
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

When analysing the ROIs during the "data-prep" step, the gfServer must used a port with the **same reference genome** used to align the genomes that will be provided later to armadillo. However, when using "armadillo run", gfServer is used just to check if the reads align perfectly anywhere in the genome so we can discard these reads. Thus, **we can use the latest reference genome** in that step, even though if the alignment was performed with a previous reference genome. It will allow us remove false positives associated to regions that actually don't exist in the reference genome used to align the genomes.

### Output

The program will output multiple files. We'll have two minibams that are generated in the first step of the pipeline. **Important note on these bams**: in order to make bwa align the reads in the same region, the reference given to bwa was a small fasta for the region of interest (for example, a gene). All reads from any copy of that gene were aligned against that reference. Finally, all minibams for each region of interest were merged into one single bam.

Also, we'll find two VCF files. The main one is \*\_nodupscandidates.vcf as here we remove duplicated hits. If gene A and B are identical and we analyse both, \*\_candidates.vcf should report the mutation in both regions. In \*\_nodupscandidates.vcf only one of them is reported.

The characteristics field in the VCF prints the next variables:
- Somatic posterior probability. Bayesian probability of the tumor's mutation VAF being greater than control's. A 0.95 cutoff is set by default to consider a variant as a somatic mutation.
- Number of mutant reads in the tumor sample
- Tumor's coverage of mutation-phased WT reads
- Total coverage in tumor
- Number of mutant reads in the control sample
- Control's coverage of mutation-phased WT reads
- Total coverage in control
- Reverse strand posterior probability. Bayesian probability of mutant reverse reads being greater than forward.
- Forward strand posterior probability. Bayesian probability of mutant forward reads being greater than reverse.
- Discarded reads. Amount of mutant reads discarded along the analysis. If too high, maybe the variant is not reliable.
- Mutant reads that didn't match the phasing due to sequence errors.
- Average noise. Average variants per position along mutant region. If high, the region may not be reliable.
- Noise standard deviation.

## Built with:

* [BWA](http://bio-bwa.sourceforge.net/)
* [Samtools](http://www.htslib.org/doc/samtools.html)
* [gfClient](https://genome.ucsc.edu/goldenPath/help/blatSpec.html#gfClientUsage)
* [Python3](https://www.python.org) 

## Authors

* **Pablo Bousquets-Muñoz** - bousquetspablo@uniovi.es
* **Ander Díaz-Navarro**
* **Xose Antón Suárez-Puente**
