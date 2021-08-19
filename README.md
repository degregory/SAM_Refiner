[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/samrefiner/README.html)
![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)

# SAM_Refiner
A program for gathering variant information from a SAM formated files.  SAM_Refiner requires a python interpreter to run.


### BioConda Install and Use

```bash
$ conda install samrefiner
$ SAM_Refiner -h
```

### Default Usage

```bash
$ python SAM_Refiner -r reference.fasta
```

## Introduction

SAM Refiner processes SAM formated files generated from sequencing mapping programs such as [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [MiniMap2](https://github.com/lh3/minimap2) to collect variant information relative to a reference sequence and remove chimeric sequences.  For details on the processing and outputs, please refer to [the manuscript][link].  If you use SAM Refiner for any published work, please cite [Gregory, DA. et al.][link].

## Updates
2021-07-11
Separated --min_abundance1 into --min_count and --min_samp_abund, --min_abundance2 now --min_col_abund

2021-06-28
Added --mp to set number of processes, fixed some memory/open file issues.

2021-06-18
Tweaked Covariant Deconvolution algorithm to handle auto-passed covars by least to greatest covar abundance.

2021-06-17
Added --read function to output the read ID and variants for each SAM read line

2021-06-07
Added AA reporting to indel in seq and covar outputs.  Reports AA for inframe indels based on the read for seq and covar outputs.  May be inconsistent w/ mismatch reporting if mismatches flank the indel.  


2021-06-06 
Added handling of non canonical NT calls.  Anything other than ATCG will not be recorded.
Changed default --max_dist to 40




## In Progress

Currently working to have a mode for whole genome sequencing.  Enabled with --wgs 1, but not yet fully implimented and tested.

Currently working on a nt call output for only variant calls.  Enabled with --ntvar 1, but not yet fully implimented and tested.


[doi]: https://doi.org/10.1101/2021.06.24.21259469
[link]: https://www.mdpi.com/1999-4915/13/8/1647



