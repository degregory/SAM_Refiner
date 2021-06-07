# SAM_Refiner
A program for gathering variant information from a SAM formated files.  SAM_Refiner requires a python iterator to run.

## Default Usage
	python SAM_Refiner -r reference.fasta

## Introduction

SAM Refiner processes SAM formated files generated from sequencing mapping programs such as Bowtie2 or BWA to collect variant information relative to a reference sequence and remove chimeric sequences.  For details on the processing and outputs, please refer to .  If you use SAM Refiner for any published work, please cite .

## Updates
2021-06-06 
Added handling of non canonical NT calls.  Anything other than ATCG will not be recorded.
Changed default --max_dist to 40
2021-06-06
Added AA reporting to indel in seq and covar outputs.  Reports AA for inframe indels based on the read for seq and covar outputs.  May be inconsistent w/ mismatch reporting if mismatches flank the indel.  



## In Progress

Currently working to have a mode for whole genome sequencing.  Enabled with --wgs 1, but not yet fully implimented and tested.

Currently working on a nt call output for only variant calls.  Enabled with --ntvar 1, but not yet fully implimented and tested.






