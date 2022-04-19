[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/samrefiner/README.html)
![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)

# SAM_Refiner
A program for gathering variant information from a SAM formatted files.  SAM_Refiner requires a python interpreter to run.  Please report any errors or innacurate outputs so the program may be improved.


### BioConda Install and Useage

```bash
$ conda install samrefiner
$ SAM_Refiner -h
```
or for the current, possibly unstable version
```bash
$ wget https://raw.githubusercontent.com/degregory/SAM_Refiner/main/SAM_Refiner.py
$ python3 /path/to/SAM_Refiner.py -h
```


### Default Usage

```bash
$ python SAM_Refiner -r reference.fasta
```

## Introduction

SAM Refiner processes SAM formatted files generated from sequencing mapping programs such as [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [MiniMap2](https://github.com/lh3/minimap2) to collect variant information relative to a reference sequence and remove chimeric sequences.  For details on the processing and outputs, please refer to [the manuscript][link].  If you use SAM Refiner for any published work, please cite [Gregory, DA. et al.][link].

## Updates
2022-04-19
Overhaul of some of the code logic.  There should be little effect on usage or output.  Usage of .gb ref and specifically --AAcentered is cleaner.  Processing for some SAM files will be improved.

2022-04-03
Added insertion reporting to nt call output.  Added amino acid centered outputs for unique seqs and covars when using .gb files for reference (--AAcentered 1).  Added a tile based method for calculating abundance for covariants in --wgs mode (more computaitonally intensive, more accurate).

2022-02-12
Improved MNP processing, added --ntcover to place a minumum coverage count for nt call output

2022-01-12
Fixed errors, improved MNP processing, added preliminary function to use gb formatted file as reference (-r)

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
Added handling of non canonical NT calls.  Anything other than ATCG- will not be recorded.
Changed default --max_dist to 40




## In Progress
Currently working on option to use gb format file for reference.  Implemented, but not extensively tested.  Simply use -r to point to a gb file.  Likely won't work on files without CDS entries.

Mode for whole genome sequencing operable.  Enabled with --wgs 1, but not yet extensively tested.

A nt call output for only variant calls operable.  Enabled with --ntvar 1, but not yet extensively tested.


[doi]: https://doi.org/10.1101/2021.06.24.21259469
[link]: https://www.mdpi.com/1999-4915/13/8/1647



