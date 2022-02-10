#!/bin/env python3
# Writen by Devon Gregory with assistance from Christopher Bottoms
# University of Missouri
# Distributed under GNU GENERAL PUBLIC LICENSE v3
import os
import sys
import argparse
import itertools
import time
from multiprocessing import Process, Pool
from pathlib import Path
"""
To Do:
Improve MNP processing for gb references
Improve nt call and nt var processing and file writing
Maybe add ins to nt call output
Clean up / polish / add comments
add --verbose --quiet options
"""



# process for collecing command line arguments
def arg_parser():

    parser = argparse.ArgumentParser(
        description='process Sam files for variant information'
    )

    parser.add_argument(
        '-r', '-reference',
        type=argparse.FileType('r'),
        dest='ref',
        help='reference fasta'
    )
    parser.add_argument(
        '-S', '--Sam_files',
        nargs='*',
        dest='Sam_files',
        action='append',
        help='optional .sam files, can use multiple files i.e. "-S Sample1.sam -S Sample2.sam" or "-S Sample1.sam Sample2.sam"'
    )
    parser.add_argument(
        '--use_count',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) use of counts in sequence IDs, default enabled (--use_count 1)'
    )
    parser.add_argument(
        '--min_count',
        type=int,
        default=10,
        help='Minimum observations required to be included in sample reports; >= 1 occurance count; < 1 %% observed (.1 = 10%%), (default: .001)'
    )
    parser.add_argument(
        '--min_samp_abund',
        type=float,
        default=0.001,
        help='Minimum abundance required for inclusion in sample reports; %% observed (.1 = 10%%), (default: .001)'
    )
    parser.add_argument(
        '--min_col_abund',
        type=float,
        default=.01,
        help='Minimum abundance required for variants to be included in collection reports; must be non-negative and  < 1, %% observed (.1 = 10%%), (default: .01)'
    )
    parser.add_argument(
        '--ntabund',
        type=float,
        default=.001,
        help='Minimum abundance relative to total reads required for a position to be reported in the nt call output; must be non-negative and < 1, %% observed (.1 = 10%%), (default: .001)'
    )
    parser.add_argument(
        '--ntcover',
        type=int,
        default=5,
        help='Minimum coverage at a position to be reported in the nt call output. (default: 5)'
    )
    parser.add_argument(
        '--max_dist',
        type=int,
        default=40,
        help='Maximum number of variances from the reference a sequence can have to be consider in covars processing (default: 40)'
    )
    parser.add_argument(
        '--max_covar',
        type=int,
        default=8,
        help='Maximum number of variances from the reference to be reported in covars (default: 8)'
    )
    parser.add_argument(
        '--AAreport',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) amino acid reporting, default enabled (--AAreport 1)'
    )
    parser.add_argument(
        '--AAcodonasMNP',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) reporting multiple nt changes in a single codon as one polymorphism, default enabled (--AAcodonasMNP 1), requires AAreport enabled'
    )
    parser.add_argument(
        '--chim_in_abund',
        type=float,
        default=.001,
        help='Minimum abundance a unique sequence must have to be considered in chimera removal / deconvolution (default: .001)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=1.2,
        help='Modifier for chim_rm chimera checking, default 1.2.  Higher = more sensitive, more false chimeras removed; lower = less sensitive, fewer chimeras removed'
    )
    parser.add_argument(
        '--foldab',
        type=float,
        default=1.8,
        help='Threshold for potential parent / chimera abundance ratio for chim_rm; default is 1.8'
    )
    parser.add_argument(
        '--redist',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) redistribution of chimera counts for chim_rm, default enabled (--redist 1)'
    )
    parser.add_argument(
        '--max_cycles',
        type=int,
        default=100,
        help='Max number of times chimera removal will be performed for chim_rm; default is 100'
    )
    parser.add_argument(
        '--beta',
        type=float,
        default=1,
        help='Modifier for covar pass checking, default 1.  Higher = more sensitive, more failed checks; lower = less sensitive, fewer failed checks'
    )
    parser.add_argument(
        '--autopass',
        type=float,
        default=.3,
        help='threshold for a sequence to automatically pass the covar pass checking'
    )
    parser.add_argument(
        '--colID',
        type=str,
        default='',
        help='ID to prepend collections'
    )
    parser.add_argument(
        '--collect',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) collection step, default enabled (--collect 1)'
    )
    parser.add_argument(
        '--read',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) reads output, default disabled (--nt_call 0)'
    )
    parser.add_argument(
        '--nt_call',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) nt_call output, default enabled (--nt_call 1)'
    )
    parser.add_argument(
        '--ntvar',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) nt_call output, default enabled (--nt_call 1)'
    )
    parser.add_argument(
        '--indel',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) indel output, default enabled (--indel 1)'
        )
    parser.add_argument(
        '--seq',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) unique seq output, default enabled (--seq 1)'
    )
    parser.add_argument(
        '--covar', type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) covar output, default enabled (--covar 1)'
    )
    parser.add_argument(
        '--pass_out',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) covar_pass output, default disabled (--pass_out 0)'
    )
    parser.add_argument(
        '--chim_rm',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) chim_rm output, default enabled (--chim_rm 1)'
    )
    parser.add_argument(
        '--deconv',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) covar deconv, default enabled (--deconv 1)'
    )
    parser.add_argument(
        '--wgs',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) covar deconv, default enabled (--deconv 1)'
    )
    parser.add_argument(
        '--mp',
        type=int,
        default=4,
        choices=range(1,21),
        help='set number of processes SAM Refiner will run in parallel, default = 4 (--mp 4)'
    )

    args = parser.parse_args()

# checking for proper range of some parameters and consistency/compatibility
    if args.wgs == 1:
        if args.deconv == 1 or args.chim_rm == 1:
            args.deconv = 0
            args.chim_rm = 0
            print('WGS mode enabled, disabling chimera removal methods')

    if args.min_count < 1:
        print(f"--min_count must be 1 or greter, setting to 1")
        args.min_count=1

    if args.min_samp_abund < 0:
        print(f"--min_samp_abund must be non-negative and < 1, defaulting to .001")
        args.min_samp_abund=0.001
    elif args.min_samp_abund >= 1:
        print(f"--min_samp_abund must be non-negative and < 1, defaulting to .001")
        args.min_samp_abund=0.001

    if args.min_col_abund < 0:
        print(f"--min_col_abund must be non-negative and < 1, defaulting to .01")
        args.min_col_abund=0.01
    elif args.min_col_abund >= 1:
        print(f"--min_col_abund must be non-negative and < 1, defaulting to .01")
        args.min_col_abund=0.01

    if args.ntabund < 0:
        print(f"--ntabund must be non-negative and < 1, defaulting to .001")
        args.ntabund=0.001
    elif args.ntabund >= 1:
        print(f"--ntabund must be non-negative and < 1, defaulting to .001")
        args.ntabund=0.001

    if args.ntcover < 1:
        print(f"--ntcover must be positive, defaulting to 5")
        args.ntcover=5

    if args.max_dist < 0:
        print(f"--max_dist must be non-negative, defaulting to 40")
        args.max_dist=40

    if args.max_covar < 0:
        print(f"--max_covar must be non-negative, defaulting to 8")
        args.max_covar=8

    if args.chim_in_abund < 0:
        print(f"--chim_in_abund must be non-negative, defaulting to 0.001")
        args.chim_in_abund=0.001
    elif args.chim_in_abund >= 1:
        print(f"--in_abund must be non-negative and < 1, defaulting to 001")
        args.chim_in_abund=0.001

    if args.alpha <= 0:
        print("--alpha must be positive, defaulting to 1.2")
        args.alpha = 1.2

    if args.foldab <= 0:
        print("--foldab must be positive, defaulting to 1.8")
        args.foldab = 1.8

    if args.max_cycles <= 0:
        print("--max_cycles must be positive, defaulting to 100")
        args.max_cycles = 100

    if args.beta <= 0:
        print("--beta must be positive, defaulting to 1")
        args.beta = 1

    if args.autopass <= 0:
        print("--autopass must be positive, defaulting to .3")
        args.autopass = .3

    return(args)

def get_ref(args): # get the reference ID and sequence from the FASTA file.  Will only get the first.
"""
From the reference provide, attempts to obtain reference ID and NT sequence, and AA sequence(s) if AAreport is enabled.
For fasta formatted files, only the first entry is parsed.
For genbank formatted files, CDS AA sequences will be pulled along with their gene ID
"""
    n=0
    refID = ''
    reftype = ''
    refseq = ''
    refORFS = {}
    if args.ref:
        ref = args.ref
        firstline = ref.readline()
        if firstline.startswith('>'):
            reftype = 'fasta'
            n+=1
            refID = firstline[1:].strip("\n\r")
            for line in ref:
                if line.startswith('>'):
                    n+=1
                    if n > 1:
                        break
                    refID = line[1:].strip("\n\r")
                elif n == 1:
                    refseq = refseq + line.strip("\n\r")
            refseq = refseq.upper()
            refprot = ''
            if args.AAreport == 1:
                for x in range(0, (len(refseq))//3):
                    AA = AAcall(refseq[x*3]+refseq[x*3+1]+refseq[x*3+2])
                    refprot = refprot + AA
                if (len(refseq))%3 != 0:
                    refprot = refprot + '?'
                refORFS = [refID, refprot]
        elif firstline.upper().startswith("LOCUS"):
            reftype = 'gb'
            collect = "Null"
            ORFS = {}
            rfs = []
            trans = 0
            gene = ""
            nts = ""
            for line in ref:
                if collect == "Null":
                    split_line = line.strip("\n\r").split(" ")

                    if split_line[0].upper() == "VERSION":
                        refID = split_line[-1]
                    elif "CDS" in line:
                        collect = "CDS"
                        if "join" in line:
                            startstops = split_line[-1].strip("join()").split(",")
                            for startstop in startstops:
                                rfs.append([int(startstop.split(".")[0]) , int(startstop.split(".")[-1])])

                        else:
                            startstop = split_line[-1].split(".")
                            rfs.append([int(startstop[0]) , int(startstop[-1])])

                    if split_line[0].upper() == "ORIGIN":
                        collect = "SEQ"
                elif collect == "CDS":
                    if "gene=" in line:
                        gene = line.split("=")[1].strip('"\n\r')
                        geneID = gene
                        n = 1
                        if gene in ORFS:
                            newgene = gene + '.' + str(n)
                            while newgene in ORFS:
                                n = n + 1
                                newgene = gene + '.' + str(n)
                            geneID = newgene

                        ORFS[geneID] = { "reading frames" : rfs

                                        }
                        rfs = []
                    elif "translation" in line:
                        try:
                            ORFS[geneID]["AAs"] = line.strip("\r\n").split('"')[1]
                        except:
                            gene = 'gene'
                            geneID = gene
                            n = 1
                            if gene in ORFS:
                                newgene = gene + '.' + str(n)
                                while newgene in ORFS:
                                    n = n + 1
                                    newgene = gene + '.' + str(n)
                                geneID = newgene
                            ORFS[geneID] = { "reading frames" : rfs
                                        }
                            rfs = []
                            ORFS[geneID]["AAs"] = line.strip("\r\n").split('"')[1]

                        if not line.strip("\r\n")[-1] == '"':
                            trans = 1
                        else:
                            ORFS[geneID]["AAs"] += "*"
                            gene = ""
                            geneID = ""
                            collect = "Null"

                    elif trans == 1:
                        ORFS[geneID]["AAs"] = ORFS[geneID]["AAs"] + line.strip(' "\n\r')
                        if line.strip("\r\n")[-1] == '"':
                            ORFS[geneID]["AAs"] += "*"
                            trans = 0
                            gene = ""
                            geneID = ""
                            collect = "Null"


                elif collect == "SEQ":
                    if not "//" in line:
                        ntsline = ""
                        for c in line:
                            if c.isalpha():
                                ntsline += c
                        nts += ntsline
            if args.AAreport == 1:
                for gene in ORFS:
                    orfnts = ''
                    for rf in ORFS[gene]['reading frames']:
                        orfnts += nts[rf[0]-1:rf[1]]
                    ORFS[gene]["nts"] = orfnts.upper()
            refORFS = ORFS
            refseq = nts.upper()

    # print(refID)


    return(refID, refseq, reftype, refORFS)

def AAcall(codon): # amino acid / codon dictionary to return encoded AAs
    AAdict = {
        'TTT' : 'F',
        'TTC' : 'F',
        'TTA' : 'L',
        'TTG' : 'L',
        'TCT' : 'S',
        'TCC' : 'S',
        'TCA' : 'S',
        'TCG' : 'S',
        'TAT' : 'Y',
        'TAC' : 'Y',
        'TAA' : '*',
        'TAG' : '*',
        'TGT' : 'C',
        'TGC' : 'C',
        'TGA' : '*',
        'TGG' : 'W',

        'CTT' : 'L',
        'CTC' : 'L',
        'CTA' : 'L',
        'CTG' : 'L',
        'CCT' : 'P',
        'CCC' : 'P',
        'CCA' : 'P',
        'CCG' : 'P',
        'CAT' : 'H',
        'CAC' : 'H',
        'CAA' : 'Q',
        'CAG' : 'Q',
        'CGT' : 'R',
        'CGC' : 'R',
        'CGA' : 'R',
        'CGG' : 'R',

        'ATT' : 'I',
        'ATC' : 'I',
        'ATA' : 'I',
        'ATG' : 'M',
        'ACT' : 'T',
        'ACC' : 'T',
        'ACA' : 'T',
        'ACG' : 'T',
        'AAT' : 'N',
        'AAC' : 'N',
        'AAA' : 'K',
        'AAG' : 'K',
        'AGT' : 'S',
        'AGC' : 'S',
        'AGA' : 'R',
        'AGG' : 'R',

        'GTT' : 'V',
        'GTC' : 'V',
        'GTA' : 'V',
        'GTG' : 'V',
        'GCT' : 'A',
        'GCC' : 'A',
        'GCA' : 'A',
        'GCG' : 'A',
        'GAT' : 'D',
        'GAC' : 'D',
        'GAA' : 'E',
        'GAG' : 'E',
        'GGT' : 'G',
        'GGC' : 'G',
        'GGA' : 'G',
        'GGG' : 'G',

        '---' : '-'
    }
    AA = '?'
    if codon in AAdict:
        AA = AAdict[codon]

    return(AA)

def singletCodon(ntPOS, nt, ref): # process to return the AA and protein seq. position based on the reference and provided nt seq position and nt
    AAPOS = (ntPOS-1)//3
    AAmod = (ntPOS-1)%3
    codon = ""
    try:
        if AAmod == 0:
            codon = nt+ref[ntPOS]+ref[ntPOS+1]
        elif AAmod == 1:
            codon = ref[ntPOS-2]+nt+ref[ntPOS]
        elif AAmod == 2:
            codon = ref[ntPOS-3]+ref[ntPOS-2]+nt
    except:
        codon = "XXX"

    return(AAPOS+1, AAcall(codon))

def getCombos(qlist, clen): # returns combinations of single polymorphisms in a sequence
    combos = []
    if (clen == 0 or clen > len(qlist)):
        clen = len(qlist)
    for N in range(1, clen+1):
        for comb in itertools.combinations(qlist, N):
            combos.append(' '.join(comb))
    return(combos)

def faSAMparse(args, ref, file): # process SAM files

    samp=file[0: -4]
    print(f"Starting {samp} processing")
    # print(ref[1])
    nt_call_dict_dict = {}
    indel_dict = {}
    seq_species = {}
    sam_read_count = 0
    sam_line_count = 0
    coverage = {}
    if args.read == 1:
        readID = ''
        reads_fh = open(samp+'_reads.tsv', "w")

    sam_fh = open(file, "r")
    for line in sam_fh:
        if not line.startswith('@'): # ignore header lines
            splitline = line.split("\t")
            if ref[0].upper().startswith(splitline[2].upper()): # check map ID matches referecne ID
                if int(splitline[4]) > 0:  # Check mapping score is positive

                    reads_count=1
                    if args.use_count == 1: # get the unique sequence counts
                        if '-' in splitline[0] and '=' in splitline[0]:
                            try:
                                eq_split = splitline[0].split('=')
                                dash_split = splitline[0].split('-')
                                if len(eq_split[-1]) > len(dash_split[-1]):
                                    reads_count=int(dash_split[-1])
                                else:
                                    reads_count=int(eq_split[-1])
                            except:
                                pass

                        elif '-' in splitline[0]:
                            try:
                                reads_count=int(splitline[0].split('-')[-1])
                            except:
                                # print(splitline[0])
                                pass

                        elif '=' in splitline[0]:
                            try:
                                reads_count=int(splitline[0].split('=')[-1])
                            except:
                                pass

                    sam_read_count += reads_count
                    sam_line_count += 1

                    # if sam_read_count % 5000 == 0:
                        # print(f"At read {sam_read_count} of {samp}")
                    # if sam_line_count % 5000 == 0:
                        # print(f"At line {sam_line_count} of {samp} SAM")

                    CIGAR = splitline[5]
                    POS = int(splitline[3])

                    readID = splitline[0]
                    query_seq = splitline[9].upper()
                    run_length = 0
                    query_seq_parsed = ''
                    query_pos = 0
                    q_pars_pos = 0
                    mutations = []
                    offsets = {}
                    running_offset = 0

                    for C in CIGAR: # process sequence based on standard CIGAR line
                        if C == 'M' or C == 'I' or C == 'D' or C == 'S' or C == 'H':
                            if C == 'S':
                                query_pos = query_pos + run_length
                            # if C == 'H':

                            if C == 'I':
                                if query_pos > 0:
                                    # add insertion to dict
                                    iPOS = q_pars_pos+POS

                                    running_offset += run_length

                                    iSeq = query_seq[query_pos: query_pos+run_length]
                                    istring = str(iPOS)+'-insert'+iSeq

                                    if args.AAreport == 1 and (run_length % 3 == 0):

                                        iProt = ''
                                        if iPOS % 3 == 1:
                                            for x in range(0, (run_length//3)):
                                                AA = AAcall(iSeq[x*3]+iSeq[x*3+1]+iSeq[x*3+2])
                                                iProt = iProt + AA
                                            istring = istring + '(' + str(((iPOS-1)//3)+1) + iProt + ')'
                                        elif iPOS % 3 == 2:
                                            if query_pos > 0:
                                                ipSeq = query_seq[query_pos-1:query_pos+run_length+5]
                                            else:
                                                ipSeq = "XXX"+query_seq[query_pos+2:query_pos+run_length+5]
                                            for x in range(0, (run_length//3)+2):
                                                AA = AAcall(ipSeq[x*3]+ipSeq[x*3+1]+ipSeq[x*3+2])
                                                iProt = iProt + AA
                                            istring = istring + '(' + ref[3][1][(iPOS-1)//3:((iPOS-1)//3)+2] + str(((iPOS-1)//3)+1) + '-' + str(((iPOS-1)//3)+2) + iProt + ')'
                                        else:
                                            if query_pos > 1:
                                                ipSeq = query_seq[query_pos-2:query_pos+run_length+4]
                                            else:
                                                ipSeq = "XXX"+query_seq[query_pos+1:query_pos+run_length+4]
                                            for x in range(0, (run_length//3)+2):
                                                AA = AAcall(ipSeq[x*3]+ipSeq[x*3+1]+ipSeq[x*3+2])
                                                iProt = iProt + AA
                                            istring = istring + '(' + ref[3][1][(iPOS-1)//3:((iPOS-1)//3)+2] + str(((iPOS-1)//3)+1) + '-' + str(((iPOS-1)//3)+2) + iProt + ')'
                                    elif args.AAreport == 1:
                                        istring = istring+'(fs)'

                                    mutations.append(istring)

                                    if args.indel ==  1:
                                        try:
                                            indel_dict[istring]
                                        except:
                                            indel_dict[istring] = reads_count
                                        else:
                                            indel_dict[istring] += reads_count

                                    query_pos = query_pos + run_length

                            elif C == 'D':
                                for X in range(0, run_length):
                                    query_seq_parsed += '-'

                                delstring = str(q_pars_pos+POS)+'-'+str(q_pars_pos+POS+run_length-1)+'Del'

                                running_offset -= run_length
                                for i in range(q_pars_pos+POS+1, q_pars_pos+POS+1+run_length):
                                    offsets[i] = running_offset


                                if args.AAreport == 1 and (run_length % 3 == 0) and not ((q_pars_pos+POS) % 3 == 1 ):
                                    if (q_pars_pos+POS) % 3 == 2:
                                        newcodon = query_seq[query_pos-1:query_pos+2]
                                        newAArefpos = (q_pars_pos+POS-1) // 3
                                        delstring = delstring + '(' + ref[3][1][newAArefpos:newAArefpos+(run_length//3)+1] + str(newAArefpos+1) + '-' + str(newAArefpos+1+run_length//3) + AAcall(newcodon)
                                        for i in range(0, run_length//3):
                                            delstring = delstring + '-'
                                        delstring = delstring + ')'
                                    else:
                                        newcodon = query_seq[query_pos-2:query_pos+1]
                                        newAArefpos = (q_pars_pos+POS-1) // 3
                                        delstring = delstring + '(' + ref[3][1][newAArefpos:newAArefpos+(run_length//3)+1] + str(newAArefpos+1) + '-' + str(newAArefpos+1+run_length//3) + AAcall(newcodon)
                                        for i in range(0, run_length//3):
                                            delstring = delstring + '-'
                                        delstring = delstring + ')'
                                elif args.AAreport == 1 and (run_length % 3 == 0):
                                    newAArefpos = (q_pars_pos+POS-1) // 3
                                    if run_length // 3 == 1:

                                        delstring = delstring + '(' + ref[3][1][newAArefpos:newAArefpos+(run_length//3)] + str(newAArefpos+1) + '-' #  + str(newAArefpos+run_length//3)

                                    else:
                                        delstring = delstring + '(' + ref[3][1][newAArefpos:newAArefpos+(run_length//3)] + str(newAArefpos+1) + '-' + str(newAArefpos+run_length//3)
                                        for i in range(0, run_length//3):
                                            delstring = delstring + '-'
                                    delstring = delstring + ')'
                                elif args.AAreport == 1:
                                    delstring = delstring + '(fs)'

                                mutations.append(delstring)

                                if args.nt_call == 1:
                                    for N in range(q_pars_pos+POS, q_pars_pos+POS+int(run_length)):
                                        try:
                                            nt_call_dict_dict[N]
                                        except:
                                            nt_call_dict_dict[N] = {'A' : 0,
                                                                    'T' : 0,
                                                                    'C' : 0,
                                                                    'G' : 0,
                                                                    '-' : 0}
                                            nt_call_dict_dict[N]['-'] = reads_count
                                        else:
                                            nt_call_dict_dict[N]['-'] += reads_count

                                if args.indel ==  1:
                                    try:
                                        indel_dict[delstring]
                                    except:
                                        indel_dict[delstring] = int(reads_count)
                                    else:
                                        indel_dict[delstring] += int(reads_count)

                                q_pars_pos = q_pars_pos + run_length

                            elif C == 'M':
                                offset = q_pars_pos-query_pos
                                refPOS = POS+offset

                                for ntPOS in range(query_pos, query_pos+run_length):
                                    offsets[refPOS+ntPOS+1] = running_offset
                                    if query_seq[ntPOS] == 'N':
                                        if args.AAreport == 0 or args.AAcodonasMNP == 0:
                                            mutations.append(ref[1][refPOS+ntPOS-1]+str(refPOS+ntPOS)+query_seq[ntPOS])
                                    if query_seq[ntPOS] == 'A' or query_seq[ntPOS] == 'T' or query_seq[ntPOS] == 'C' or query_seq[ntPOS] == 'G' or query_seq[ntPOS] == '-':
                                        if query_seq[ntPOS] != ref[1][refPOS+ntPOS-1]:
                                            if args.AAreport == 1 and args.AAcodonasMNP == 0:
                                                AAinfo = singletCodon(refPOS+ntPOS, query_seq[ntPOS], ref[1])
                                                mutations.append(ref[1][refPOS+ntPOS-1]+str(refPOS+ntPOS)+query_seq[ntPOS]+'('+ref[3][1][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+')')
                                            else:
                                                mutations.append(ref[1][refPOS+ntPOS-1]+str(refPOS+ntPOS)+query_seq[ntPOS])
                                        if args.nt_call == 1:
                                            try:
                                                nt_call_dict_dict[refPOS+ntPOS]
                                            except:
                                                nt_call_dict_dict[refPOS+ntPOS] = {'A' : 0,
                                                                                   'T' : 0,
                                                                                   'C' : 0,
                                                                                   'G' : 0,
                                                                                   '-' : 0}
                                                nt_call_dict_dict[refPOS+ntPOS][query_seq[ntPOS]] = reads_count
                                            else:
                                                nt_call_dict_dict[refPOS+ntPOS][query_seq[ntPOS]] += reads_count


                                q_pars_pos = q_pars_pos + run_length
                                query_pos = query_pos + run_length

                            run_length = 0


                        else:
                            run_length = (10 * run_length) + int(C)
                    # END CIGAR PARSE



                    if len(mutations) == 0: # record reference counts
                        if args.wgs == 0:
                            try:
                                seq_species['Reference'] += reads_count
                            except:
                                seq_species['Reference'] = reads_count
                        else:
                            try:
                                seq_species[str(POS)+' Ref '+str(POS+q_pars_pos)] += reads_count
                            except:
                                seq_species[str(POS)+' Ref '+str(POS+q_pars_pos)] = reads_count
                        if args.read == 1:
                            reads_fh.write(f"{readID}\tReference\n")

                    else: # record variants and counts
                        if args.AAreport == 1 and args.AAcodonasMNP == 1:
                            last_codon = -1
                            codonchecked = []
                            MNP = []
                            for mut in mutations:
                                if '(fs)' in mut:
                                    if MNP:
                                        mutstart = 0
                                        if 'Del' in mut:
                                            mutstart = int(mut.split('Del')[0].split('-')[0])
                                        elif 'insert' in mut:
                                            mutstart = int((mut.split('(')[0].split('-')[0]))

                                        MNPend = 0
                                        if 'Del' in MNP[-1]:
                                            MNPend = int(MNP[-1].split('Del')[0].split('-')[1])
                                        elif 'insert' in MNP[-1]:
                                            MNPend = int((MNP[-1].split('(')[0].split('-')[0]))
                                        else:
                                            MNPend = int(MNP[-1][1:-1])

                                        if '(fs)' in MNP[-1]:
                                            if mutstart - MNPend < 10:
                                                ntshift = 0
                                                for single in [MNP[-1], mut]:
                                                    if 'Del' in single:
                                                        ntshift -= int(single.split('Del')[0].split('-')[1]) - int(single.split('Del')[0].split('-')[0]) + 1
                                                    elif 'insert' in single:
                                                        ntshift += int(len(single.split('(')[0].split('insert')[1]))
                                                if ntshift % 3 == 0:
                                                    MNP.append(mut)
                                                else:
                                                    codonchecked.append(":".join(MNP))
                                                    MNP = [mut]
                                            else:
                                                codonchecked.append(":".join(MNP))
                                                MNP = [mut]
                                        else:
                                            mutstartcodon = ((mutstart-1)//3)+1
                                            MNPendcodon = ((MNPend-1)//3)+1
                                            if 'insert' in MNP[-1]:
                                                if '-' in MNP[-1].split('insert')[1]:
                                                    MNPendcodon += 1
                                            if mutstartcodon == MNPendcodon:
                                                MNP.append(mut)
                                            else:
                                                codonchecked.append(":".join(MNP))
                                                MNP = [mut]

                                    else:
                                        MNP.append(mut)
                                else:
                                    cur_start_codon = 0
                                    cur_end_codon = 0
                                    if 'Del' in mut:
                                        cur_start_codon = ((int(mut.split('Del')[0].split('-')[0])-1)//3)+1
                                        cur_end_codon = ((int(mut.split('Del')[0].split('-')[1])-1)//3)+1
                                    elif 'insert' in mut:
                                        cur_start_codon = ((int(mut.split('-')[0])-1)//3)+1
                                        if '-' in mut.split('insert')[1]:
                                            cur_end_codon = ((int(mut.split('-')[0])-1)//3)+2
                                        else:
                                            cur_end_codon = ((int(mut.split('-')[0])-1)//3)+1
                                    else:
                                        cur_start_codon = ((int(mut[1:-1])-1)//3)+1
                                        cur_end_codon = ((int(mut[1:-1])-1)//3)+1
                                    if MNP:
                                        if cur_start_codon == last_codon:
                                            MNP.append(mut)
                                        else:
                                            codonchecked.append(":".join(MNP))
                                            MNP = [mut]
                                    else:
                                        MNP.append(mut)
                                    last_codon = cur_end_codon
                            if MNP:
                                codonchecked.append(":".join(MNP))
                            MNP_muts = []
                            for mut in codonchecked:
                                if ":" in mut:
                                    if 'Del' in mut or 'insert' in mut:
                                        fshift = 0
                                        fs_count = 0
                                        if '(fs)' in mut:
                                            for single in mut.split(':'):
                                                if '(fs)' in single:
                                                    fs_count += 1
                                                    if 'Del' in single:
                                                        fshift -= int(single.split('Del')[0].split('-')[1]) - int(single.split('Del')[0].split('-')[0]) + 1
                                                    elif 'insert' in single:
                                                        fshift += int(len(single.split('(')[0].split('insert')[1]))
                                        if ((fshift % 3) == 0):
                                            start = 0
                                            end = 0
                                            ntshift = 0
                                            split_MNP = mut.split(":")

                                            for x in range(0, len(split_MNP)):
                                                if x == 0:
                                                    if 'Del' in split_MNP[x]:
                                                        start = ((int(split_MNP[x].split('Del')[0].split('-')[0])-1)//3)+1
                                                        ntshift -= int(split_MNP[x].split('Del')[0].split('-')[1]) - int(split_MNP[x].split('Del')[0].split('-')[0]) + 1
                                                    elif 'insert' in split_MNP[x]:
                                                        start = ((int(split_MNP[x].split('-')[0])-1)//3)+1
                                                        ntshift += int(len(split_MNP[x].split('(')[0].split('insert')[1]))
                                                    else:
                                                        start = ((int(split_MNP[x][1:-1])-1)//3)+1
                                                else:
                                                    if 'Del' in split_MNP[x]:
                                                        end = ((int(split_MNP[x].split('Del')[0].split('-')[1])-1)//3)+1
                                                        ntshift -= int(split_MNP[x].split('Del')[0].split('-')[1]) - int(split_MNP[x].split('Del')[0].split('-')[0]) + 1
                                                    elif 'insert' in split_MNP[x]:
                                                        ntshift += int(len(split_MNP[x].split('(')[0].split('insert')[1]))
                                                        if '-' in split_MNP[x].split('insert')[1]:
                                                            end = ((int(split_MNP[x].split('-')[0])-1)//3)+2
                                                        else:
                                                            end = ((int(split_MNP[x].split('-')[0])-1)//3)+1
                                                    else:
                                                        end = ((int(split_MNP[x][1:-1])-1)//3)+1
                                            start_nt = ((start-1)*3)+1
                                            end_nt = ((end-1)*3)+3
                                            new_nt_seq = query_seq[start_nt-POS-offsets[start_nt-1]:end_nt-POS+1-offsets[start_nt-1]+ntshift]


                                            if ntshift < 0:
                                                while ntshift < 0:
                                                    new_nt_seq += '-'
                                                    ntshift += 1
                                            new_aa_seq = ""
                                            for i in range(0, (len(new_nt_seq)//3)):
                                                new_aa_seq += AAcall(new_nt_seq[i*3:(i*3)+3])

                                            MNP_muts.append(ref[1][start_nt-1:end_nt]+str(((start-1)*3)+1)+'-'+str(((end-1)*3)+3)+new_nt_seq+ '(' + ref[3][1][start-1:end] + str(start) + '-' + str(end) + new_aa_seq +')')
                                        else:
                                            MNP_muts.append(mut)

                                    else:
                                        split_MNP = mut.split(":")
                                        codon = ((int(split_MNP[0][1:-1])-1)//3)
                                        old_triplet = ref[1][(codon*3):(codon*3)+3]
                                        new_triplet = [old_triplet[0], old_triplet[1], old_triplet[2]] # ['r', 'r', 'r']
                                        for PM in split_MNP:
                                            new_triplet[((int(PM[1:-1])-1)%3)] = PM[-1]
                                        MNP =  old_triplet + str((codon*3)+1) + "".join(new_triplet) + "(" + ref[3][1][codon] + str(codon+1) + AAcall("".join(new_triplet)) + ")"
                                        MNP_muts.append(MNP)
                                elif not 'Del' in mut and not 'insert' in mut:
                                    AAinfo = singletCodon(int(mut[1:-1]), mut[-1], ref[1])
                                    MNP_muts.append(mut+'('+ref[3][1][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+')')
                                else:
                                    MNP_muts.append(mut)
                            mutations = MNP_muts

                        mutations = " ".join(mutations)

                        if args.wgs == 0:
                            try:
                                seq_species[mutations] += reads_count
                            except:
                                seq_species[mutations] = reads_count
                            if args.read == 1:
                                reads_fh.write(f"{readID}\t{mutations}\n")
                        else:
                            try:
                                seq_species[str(POS)+' '+mutations+' '+str(POS+q_pars_pos)] += reads_count
                            except:
                                seq_species[str(POS)+' '+mutations+' '+str(POS+q_pars_pos)] = reads_count
                            if args.read == 1:
                                reads_fh.write(f"{readID}\t{str(POS)} {mutations} {str(POS+q_pars_pos)}\n")

                    for i in range(POS, POS+q_pars_pos): # update coverage
                        try:
                            coverage[i] += reads_count
                        except:
                            coverage[i] = reads_count
    if args.read == 1:
        reads_fh.close()
    sam_fh.close()
    # END SAM LINES
    print(f"End SAM parse for {samp}")
    # print(coverage)

    if sam_read_count == 0:
        print(f"No Reads for {samp}")

    else:

        if args.seq == 1: # output the sequence
            seq_fh = open(samp+'_unique_seqs.tsv', "w")
            seq_fh.write(samp+"("+str(sam_read_count)+")\n")
            seq_fh.write("Unique Sequence\tCount\tAbundance\n")

            sorted_seq = sorted(seq_species, key=seq_species.__getitem__, reverse=True)
            for key in sorted_seq:
                if seq_species[key] >= args.min_count:
                    if (seq_species[key] / sam_read_count >= args.min_samp_abund) and args.wgs == 0:
                        seq_fh.write(f"{key}\t{seq_species[key]}\t{(seq_species[key]/sam_read_count):.3f}\n")
                    elif args.wgs == 1:
                        splitseqs = key.split()
                        cov = []
                        for x in range(int(splitseqs[0]), int(splitseqs[-1])):
                            cov.append(coverage[x])
                        min_cov = min(cov)
                        if (seq_species[key]/min_cov >= args.min_samp_abund):
                            seq_fh.write(f"{key}\t{seq_species[key]}\t{(seq_species[key]/min_cov):.3f}\n")
                else:
                    break

            seq_fh.close()
            # END SEQ OUT
            # print(f"End unqiue seq out for {samp}")

        if args.indel == 1 and len(indel_dict) > 0: # output indels, if there are any
            sorted_indels = sorted(indel_dict, key=indel_dict.__getitem__, reverse=True)
            indels_to_write = []
            for key in sorted_indels:
                if indel_dict[key] >= args.min_count:
                    if indel_dict[key] / sam_read_count >= args.min_samp_abund and args.wgs == 0:
                        indels_to_write.append(f"{key}\t{indel_dict[key]}\t{(indel_dict[key]/sam_read_count):.3f}\n")
                    elif args.wgs == 1:
                        indelPOS = ''
                        for c in key:
                            if c.isdigit():
                                indelPOS += c
                            else:
                                break
                        indelPOS = int(indelPOS)
                        if indel_dict[key] / coverage[indelPOS] >= args.min_samp_abund:
                            indels_to_write.append(f"{key}\t{indel_dict[key]}\t{(indel_dict[key] / coverage[indelPOS]):.3f}\n")
                else:
                    break
            if len(indels_to_write) > 0:
                indel_fh = open(samp+'_indels.tsv', "w")
                indel_fh.write(samp+"("+str(sam_read_count)+")\n")
                indel_fh.write("Indel\tCount\tAbundance\n")
                for indel_entry in indels_to_write:
                    indel_fh.write(indel_entry)
                indel_fh.close()
            # END INDEL OUT
            # print(f"End indel out for {samp}")

        if args.nt_call == 1: # out put nt calls
            ntcall_lines = {'line' : {},
                            'variant' : {}
                            }
            ntcall_fh = open(samp+'_nt_calls.tsv', "w")
            ntcall_fh.write(samp+"("+str(sam_read_count)+")\n")
            if args.ntvar == 1:
                ntcallv_fh = open(samp+'_nt_calls_varonly.tsv', "w")
                ntcallv_fh.write(samp+"("+str(sam_read_count)+")\n")
            sorted_POS = sorted(nt_call_dict_dict)
            if args.AAreport == 1:

                ntcall_fh.write("Position\tref NT\tAA POS\tref AA\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                if args.ntvar == 1:
                    ntcallv_fh.write("Position\tref NT\tAA POS\tref AA\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                consensus = {}
                for POS in sorted_POS:
                    try:
                        total = coverage[POS]
                    except:
                        total = 0
                    if total >= (sam_read_count * args.ntabund) and total >= args.ntcover:
                        # AAinfo = singletCodon(POS, ref[1][POS-1], ref)
                        POS_calls = {}
                        for key in nt_call_dict_dict[POS]:
                            POS_calls[key] = nt_call_dict_dict[POS][key]
                        sorted_calls = sorted(POS_calls, key=POS_calls.__getitem__, reverse=True)

                        ntcall_lines['line'][POS] = (str(POS)+"\t"+ref[1][POS-1]+"\t"+str(((POS-1)//3)+1)+"\t"+ref[3][1][((POS-1)//3)])
                        ntcall_lines['line'][POS] +=("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                        ntcall_lines['line'][POS] +=("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[0]]))
                        ntcall_lines['line'][POS] +=(f"\t{(nt_call_dict_dict[POS][sorted_calls[0]]/total):.3f}")

                        consensus[POS] = sorted_calls


                for POS in sorted_POS:
                    try:
                        ntcall_lines['line'][POS]
                    except:
                        pass
                    else:
                        if consensus[POS][0] != ref[1][POS-1]:
                            mod = (POS)%3

                            if mod == 0:
                                try:
                                    codon = consensus[POS-2][0]+consensus[POS-1][0]+ consensus[POS][0]
                                except:
                                    codon = 'NNN'
                            elif mod == 2:
                                try:
                                    codon = consensus[POS-1][0]+consensus[POS][0]+ consensus[POS+1][0]
                                except:
                                    codon = 'NNN'
                            elif mod == 1:
                                try:
                                    codon = consensus[POS][0]+consensus[POS+1][0]+ consensus[POS+2][0]
                                except:
                                    codon = 'NNN'
                            ntcall_lines['line'][POS] +=("\t"+AAcall(codon)+"\t"+singletCodon(POS, consensus[POS][0], ref[1])[1])
                        else:
                            ntcall_lines['line'][POS] +=("\t\t")

                        if (nt_call_dict_dict[POS][consensus[POS][1]] >= args.min_count) and ((nt_call_dict_dict[POS][consensus[POS][1]] / total) >= args.min_samp_abund):
                            ntcall_lines['line'][POS] +=(f"\t{consensus[POS][1]}\t{nt_call_dict_dict[POS][consensus[POS][1]]}\t{(nt_call_dict_dict[POS][consensus[POS][1]]/total):.3f}"+"\t"+singletCodon(POS, consensus[POS][1], ref[1])[1])

                            if (nt_call_dict_dict[POS][consensus[POS][2]] >= args.min_count) and (nt_call_dict_dict[POS][consensus[POS][2]] / total >= args.min_samp_abund):
                                ntcall_lines['line'][POS] +=(f"\t{consensus[POS][2]}\t{nt_call_dict_dict[POS][consensus[POS][2]]}\t{(nt_call_dict_dict[POS][consensus[POS][2]]/total):.3f}\t{singletCodon(POS, consensus[POS][2], ref[1])[1]}")

                for POS in ntcall_lines['line']:
                    ntcall_fh.write(ntcall_lines['line'][POS])
                    ntcall_fh.write("\n")
                    if args.ntvar == 1:
                        try:
                            ntcall_lines['variant'][POS]
                        except:
                            pass
                        else:
                            ntcallv_fh.write(ntcall_lines['line'][POS])
                            ntcallv_fh.write("\n")
                if args.ntvar == 1:
                        ntcallv_fh.close()

            else:
                ntcall_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")
                if args.ntvar == 1:
                    ntcallv_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")

                for POS in sorted_POS:
                    try:
                        total = coverage[POS] # sum(nt_call_dict_dict[POS].values())
                    except:
                        total = 0
                    if total >= (sam_read_count * args.ntabund) and total >= args.ntcover:
                        POS_calls = {}
                        for key in nt_call_dict_dict[POS]:
                            POS_calls[key] = nt_call_dict_dict[POS][key]
                        sorted_calls = sorted(POS_calls, key=POS_calls.__getitem__, reverse=True)

                        ntcall_fh.write(str(POS)+"\t"+ref[1][POS-1])
                        ntcall_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                        ntcall_fh.write("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[0]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[0]]/total):.3f}")

                        if (nt_call_dict_dict[POS][sorted_calls[1]] >= args.min_count) and (nt_call_dict_dict[POS][sorted_calls[1]] / total) >= args.min_samp_abund:
                            if args.ntvar == 1:
                                ntcallv_fh.write(str(POS)+"\t"+ref[1][POS-1])
                                ntcallv_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                                ntcallv_fh.write("\t"+str(total)+"\t\t")
                                ntcallv_fh.write(f"\t")
                                ntcallv_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                            ntcall_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                            if (nt_call_dict_dict[POS][sorted_calls[2]] > args.min_count) and (nt_call_dict_dict[POS][sorted_calls[2]] /total > args.min_samp_abund):
                                ntcall_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                                if args.ntvar == 1:
                                    ntcallv_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                            if args.ntvar == 1:
                                ntcallv_fh.write("\n")

                        ntcall_fh.write("\n")

            ntcall_fh.close()
            if args.ntvar == 1:
                ntcallv_fh.close()
            # END NT CALL OUT
            # print(f"End nt call out for {samp}")
        if args.covar == 1: # output covariants
            testtrack = 0
            combinations = {}
            for sequence in seq_species:
                if args.wgs == 0:
                    singles = sequence.split()
                else:
                    singles = (sequence.split())[1:-1]
                if len(singles) <= args.max_dist and singles[0] != 'Ref':
                    for combo in getCombos(singles, args.max_covar):
                        if not combo in combinations:
                            combinations[combo] = seq_species[sequence]
                        else:
                            combinations[combo] += seq_species[sequence]

            covar_fh = open(samp+'_covars.tsv', "w")
            covar_fh.write(samp+"("+str(sam_read_count)+")\n")
            covar_fh.write("Co-Variants\tCount\tAbundance\n")
            sortedcombos = sorted(combinations, key=combinations.__getitem__, reverse=True)
            # print(sortedcombos)
            for key in sortedcombos:
                if (combinations[key] >= args.min_count):
                    if (combinations[key] / sam_read_count >= args.min_samp_abund) and args.wgs == 0:
                        covar_fh.write(key+"\t"+str(combinations[key])+"\t"+f"{(combinations[key]/sam_read_count):.3f}\n")
                    elif args.wgs == 1:
                        coveragepercent = 0
                            # print(key)
                        splitcombos = key.split()
                        if len(splitcombos) == 1:
                            coveragePOS = ''
                            for c in key.strip('ATGC'):
                                if c.isdigit():
                                    coveragePOS += c
                                else:
                                    break
                            # print(key)
                            # print(coveragePOS)
                            coveragepercent = combinations[key] / coverage[int(coveragePOS)]
                        else:
                            startcovPOS = ''
                            for c in splitcombos[0].strip('ATGC'):
                                if c.isdigit():
                                    startcovPOS += c
                                else:
                                    break
                            endcovPOS = ''
                            for c in splitcombos[-1].strip('ATGC'):
                                if c.isdigit():
                                    endcovPOS += c
                                else:
                                    break
                            coveragevals = []
                            for i in range(int(startcovPOS), int(endcovPOS)+1):
                                coveragevals.append(coverage[i])
                            mincov = min(coverval for coverval in coveragevals)
                            coveragepercent = combinations[key] / mincov
                        if coveragepercent >= args.min_samp_abund:
                            covar_fh.write(f"{key}\t{combinations[key]}\t{coveragepercent:.3f}\n")

                    # \t{coveragepercent:.3f}
            covar_fh.close()
            # END COVAR OUT
            # print(f"End covar out for {samp}")

def gbSAMparse(args, ref, file): # process SAM files

    samp=file[0: -4]
    print(f"Starting {samp} processing")
    nt_call_dict_dict = {}
    indel_dict = {}
    seq_species = {}
    sam_read_count = 0
    sam_line_count = 0
    coverage = {}
    if args.read == 1:
        readID = ''
        reads_fh = open(samp+'_reads.tsv', "w")

    sam_fh = open(file, "r")
    for line in sam_fh:
        if not line.startswith('@'): # ignore header lines
            splitline = line.split("\t")

            if ref[0].upper().startswith(splitline[2].upper()): # check map ID matches referecne ID
                if int(splitline[4]) > 0:  # Check mapping score is positive

                    reads_count=1
                    if args.use_count == 1: # get the unique sequence counts
                        if '-' in splitline[0] and '=' in splitline[0]:
                            eq_split = splitline[0].split('=')
                            dash_split = splitline[0].split('-')
                            if len(eq_split[-1]) > len(dash_split[-1]):
                                reads_count=int(dash_split[-1])
                            else:
                                reads_count=int(eq_split[-1])

                        elif '-' in splitline[0]:
                            try:
                                reads_count=int(splitline[0].split('-')[-1])
                            except:
                                # print(splitline[0])
                                pass

                        elif '=' in splitline[0]:
                            try:
                                reads_count=int(splitline[0].split('=')[-1])
                            except:
                                pass


                    sam_read_count += reads_count
                    sam_line_count += 1

                    # if sam_read_count % 5000 == 0:
                        # print(f"At read {sam_read_count} of {samp}")
                    # if sam_line_count % 5000 == 0:
                        # print(f"At line {sam_line_count} of {samp} SAM")

                    CIGAR = splitline[5]
                    POS = int(splitline[3])

                    readID = splitline[0]
                    query_seq = splitline[9].upper()
                    run_length = 0
                    query_seq_parsed = ''
                    query_pos = 0
                    q_pars_pos = 0
                    mutations = []

                    for C in CIGAR: # process sequence based on standard CIGAR line
                        if C == 'M' or C == 'I' or C == 'D' or C == 'S' or C == 'H':
                            if C == 'S':
                                query_pos = query_pos + run_length
                            # if C == 'H':


                            if C == 'I':
                                if query_pos > 0:
                                    # add insertion to dict
                                    iPOS = q_pars_pos+POS

                                    iSeq = query_seq[query_pos: query_pos+run_length]
                                    istring = str(iPOS)+'-insert'+iSeq


                                    if args.AAreport == 1:

                                        iProt = ''
                                        for orf in ref[3]:
                                            orflength = 0
                                            for rf in ref[3][orf]['reading frames']:
                                                if iPOS >= rf[0] and iPOS <= rf[1]:
                                                    iProt += "(" + orf + "_"
                                                    orfPOS = 1 + iPOS - rf[0] + orflength
                                                    iProt += "nt:" + str(orfPOS)
                                                    if (run_length % 3 == 0):
                                                        if orfPOS % 3 == 1:
                                                            for x in range(0, (run_length//3)):
                                                                AA = AAcall(iSeq[x*3]+iSeq[x*3+1]+iSeq[x*3+2])
                                                        elif orfPOS % 3 == 2:
                                                            if query_pos > 0:
                                                                ipSeq = query_seq[query_pos-1:query_pos+run_length+2]
                                                            else:
                                                                ipSeq = "XXX"+query_seq[query_pos+2:query_pos+run_length+2]
                                                            for x in range(0, (run_length//3)+1):
                                                                AA = AAcall(ipSeq[x*3]+ipSeq[x*3+1]+ipSeq[x*3+2])
                                                        else:
                                                            if query_pos > 1:
                                                                ipSeq = query_seq[query_pos-2:query_pos+run_length+1]
                                                            else:
                                                                ipSeq = "XXX"+query_seq[query_pos+1:query_pos+run_length+1]

                                                            for x in range(0, (run_length//3)+1):
                                                                AA = AAcall(ipSeq[x*3]+ipSeq[x*3+1]+ipSeq[x*3+2])
                                                        iProt = iProt + "_AA:" + str((orfPOS//3)+1) + AA
                                                    else:
                                                        iProt += "_AA:" + str((orfPOS//3)+1) + "fs"
                                                    iProt += ")"
                                                orflength += rf[1] - rf[0] + 1

                                        istring += iProt




                                    mutations.append(istring)

                                    if args.indel ==  1:
                                        try:
                                            indel_dict[istring]
                                        except:
                                            indel_dict[istring] = reads_count
                                        else:
                                            indel_dict[istring] += reads_count

                                    query_pos = query_pos + run_length

                            elif C == 'D':
                                for X in range(0, run_length):
                                    query_seq_parsed += '-'

                                delPOS = q_pars_pos+POS
                                delstring = str(delPOS)+'-'+str(delPOS+run_length-1)+'Del'
                                if args.AAreport == 1:
                                    delProt = ""
                                    for orf in ref[3]:
                                        orflength = 0
                                        for rf in ref[3][orf]['reading frames']:
                                            if delPOS >= rf[0] and delPOS <= rf[1]:
                                                delProt += "(" + orf + "_"
                                                orfPOS = 1 + delPOS - rf[0] + orflength
                                                delProt += "nt:" + str(orfPOS) + "-" + str(orfPOS+run_length-1)
                                                if (delPOS+run_length-1) <= rf[1]:
                                                    if (run_length % 3 == 0):
                                                        delProt += "_AA:"
                                                        if ((orfPOS) % 3 == 1 ):
                                                            delProt +=  ref[3][orf]["AAs"][orfPOS//3:orfPOS//3+(run_length//3)] +str((orfPOS//3)+1) + "-" + str((orfPOS//3)+(int(run_length/3))) + 'del'
                                                        else:
                                                            newAArefpos = ((orfPOS) // 3)
                                                            if (orfPOS) % 3 == 2:
                                                                newcodon = query_seq[query_pos-1:query_pos+2]
                                                            else:
                                                                newcodon = query_seq[query_pos-2:query_pos+1]
                                                            delProt += ref[3][orf]["AAs"][newAArefpos] + str(newAArefpos+1) + AAcall(newcodon)
                                                            if run_length > 3:
                                                                delProt += ":" + ref[3][orf]["AAs"][(orfPOS//3)+1:orfPOS//3+1+(run_length//3)] + str((orfPOS//3)+2) + "-" + str((orfPOS//3)+1+(int(run_length/3))) + 'del'

                                                    else:
                                                        delProt += "_AA:" + str((orfPOS//3)+1) + "fs"
                                                    delProt += ")"

                                                else:
                                                    delProt += "Frame shift / Splicing / Terminating codon disrupted"
                                            elif (delPOS+run_length-1) >= rf[0] and (delPOS+run_length-1) <= rf[1]:
                                                delProt += "Frame shift / Splicing / Start codon disrupted"
                                            orflength += rf[1] - rf[0] + 1
                                    delstring += delProt


                                mutations.append(delstring)

                                if args.nt_call == 1:
                                    for N in range(q_pars_pos+POS, q_pars_pos+POS+int(run_length)):
                                        try:
                                            nt_call_dict_dict[N]
                                        except:
                                            nt_call_dict_dict[N] = {'A' : 0,
                                                                    'T' : 0,
                                                                    'C' : 0,
                                                                    'G' : 0,
                                                                    '-' : 0}
                                            nt_call_dict_dict[N]['-'] = reads_count
                                        else:
                                            nt_call_dict_dict[N]['-'] += reads_count

                                if args.indel ==  1:
                                    try:
                                        indel_dict[delstring]
                                    except:
                                        indel_dict[delstring] = int(reads_count)
                                    else:
                                        indel_dict[delstring] += int(reads_count)

                                q_pars_pos = q_pars_pos + run_length

                            elif C == 'M':
                                offset = q_pars_pos-query_pos
                                refPOS = POS+offset

                                for ntPOS in range(query_pos, query_pos+run_length):
                                    if query_seq[ntPOS] == 'A' or query_seq[ntPOS] == 'T' or query_seq[ntPOS] == 'C' or query_seq[ntPOS] == 'G' or query_seq[ntPOS] == '-':
                                        PMPOS = refPOS+ntPOS
                                        if query_seq[ntPOS] != ref[1][PMPOS-1]:
                                            PM = ref[1][PMPOS-1] +str(PMPOS)+query_seq[ntPOS]
                                            if args.AAreport == 1 and args.AAcodonasMNP == 0:
                                                for orf in ref[3]:
                                                    orflength = 0
                                                    for rf in ref[3][orf]['reading frames']:
                                                        if PMPOS >= rf[0] and PMPOS <= rf[1]:
                                                            ORFPOS = PMPOS-rf[0]+1+orflength
                                                            AAinfo = singletCodon(ORFPOS, query_seq[ntPOS], (ref[3][orf]['nts']))
                                                            PM += '(' + orf + "_nts:" + str(ORFPOS) + "_AA:" + ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+')'
                                                        orflength += rf[1] - rf[0] + 1
                                            mutations.append(PM)
                                        if args.nt_call == 1:
                                            try:
                                                nt_call_dict_dict[PMPOS]
                                            except:
                                                nt_call_dict_dict[PMPOS] = {'A' : 0,
                                                                                   'T' : 0,
                                                                                   'C' : 0,
                                                                                   'G' : 0,
                                                                                   '-' : 0}
                                                nt_call_dict_dict[PMPOS][query_seq[ntPOS]] = reads_count
                                            else:
                                                nt_call_dict_dict[PMPOS][query_seq[ntPOS]] += reads_count


                                q_pars_pos = q_pars_pos + run_length
                                query_pos = query_pos + run_length

                            run_length = 0


                        else:
                            run_length = (10 * run_length) + int(C)
                    # END CIGAR PARSE



                    if len(mutations) == 0: # record reference counts
                        if args.wgs == 0:
                            try:
                                seq_species['Reference'] += reads_count
                            except:
                                seq_species['Reference'] = reads_count
                        else:
                            try:
                                seq_species[str(POS)+' Ref '+str(POS+q_pars_pos)] += reads_count
                            except:
                                seq_species[str(POS)+' Ref '+str(POS+q_pars_pos)] = reads_count
                        if args.read == 1:
                            reads_fh.write(f"{readID}\tReference\n")

                    else: # record variants and counts
                        if args.AAreport == 1 and args.AAcodonasMNP == 1:
                            mutorfs = {}
                            for orf in ref[3]:
                                for mut in mutations:
                                    if not 'Del' in mut and not 'insert' in mut:
                                        mutPOS = int(''.join([c for c in mut if c.isdigit()]))
                                        frame = 0
                                        orflength = 0
                                        for rf in ref[3][orf]['reading frames']:
                                            if mutPOS >= rf[0] and mutPOS <= rf[1]:
                                                orfPOS = mutPOS-ref[3][orf]['reading frames'][frame][0]+1+orflength
                                                try:
                                                    mutorfs[orf]['muts'].append([frame, mut, orfPOS])
                                                except:
                                                    mutorfs[orf] = {'muts' : [[frame, mut, orfPOS]],
                                                                    'mutstrings' : {}
                                                                    }
                                            frame += 1
                                            orflength += rf[1] - rf[0] + 1
                            for orf in mutorfs:
                                codon = ''
                                skip = 0
                                MNP = ''
                                for i in range(0, len(mutorfs[orf]['muts'])): # checking for MNP

                                    if skip > 0:
                                        skip -= 1
                                    else:
                                        mut1POS = mutorfs[orf]['muts'][i][2]
                                        try:
                                            mutorfs[orf]['muts'][i+1]
                                        except:
                                            AAinfo = singletCodon(mut1POS, mutorfs[orf]['muts'][i][1][-1], ref[3][orf]['nts'])
                                            try:
                                                mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] += ", nts:" + str(mutorfs[orf]['muts'][i][2])+'_AAs:'+ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]
                                            except:
                                                mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] = "nts:" + str(mutorfs[orf]['muts'][i][2])+'_AAs:'+ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]
                                        else:
                                            if mut1POS % 3 == 0:
                                                AAinfo = singletCodon(mut1POS, mutorfs[orf]['muts'][i][1][-1], ref[3][orf]['nts'])
                                                try:
                                                    mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] += ", nts:" + str(mutorfs[orf]['muts'][i][2])+'_AAs:'+ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]
                                                except:
                                                    mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] = "nts:" + str(mutorfs[orf]['muts'][i][2])+'_AAs:'+ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]
                                            else:
                                                mut2POS = mutorfs[orf]['muts'][i+1][2]
                                                if mut2POS - mut1POS < 3:
                                                    AAPOS = mut1POS // 3
                                                    if mut1POS % 3 == 1:
                                                        if mut2POS % 3 == 0:
                                                            codon = mutorfs[orf]['muts'][i][1][-1]+ref[3][orf]['nts'][mut1POS]+mutorfs[orf]['muts'][i+1][1][-1]
                                                            MNP = mutorfs[orf]['muts'][i][1][-1]+'r'+mutorfs[orf]['muts'][i+1][1][-1]
                                                            skip = 1
                                                        else:
                                                            try:
                                                                mutorfs[orf]['muts'][i+2]
                                                            except:
                                                                codon = mutorfs[orf]['muts'][i][1][-1]+mutorfs[orf]['muts'][i+1][1][-1]+ref[3][orf]['nts'][mut2POS]
                                                                MNP = mutorfs[orf]['muts'][i][1][-1]+mutorfs[orf]['muts'][i+1][1][-1]+'r'
                                                                skip = 1
                                                            else:
                                                                mut3POS = mutorfs[orf]['muts'][i+2][2]
                                                                if mut2POS == mut3POS -1:
                                                                    codon = mutorfs[orf]['muts'][i][1][-1]+mutorfs[orf]['muts'][i+1][1][-1]+mutorfs[orf]['muts'][i+2][1][-1]
                                                                    MNP = mutorfs[orf]['muts'][i][1][-1]+mutorfs[orf]['muts'][i+1][1][-1]+mutorfs[orf]['muts'][i+2][1][-1]
                                                                    skip = 2
                                                                else:
                                                                    codon = mutorfs[orf]['muts'][i][1][-1]+mutorfs[orf]['muts'][i+1][1][-1]+ref[3][orf]['nts'][mut2POS]
                                                                    MNP = mutorfs[orf]['muts'][i][1][-1]+mutorfs[orf]['muts'][i+1][1][-1]+'r'
                                                                    skip = 1
                                                        try:
                                                            mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] += ", nts:"+ ref[3][orf]['nts'][mut1POS-1:mut1POS+2] + str(mut1POS)+MNP+'_AAs:'+ref[3][orf]["AAs"][AAPOS]+str(AAPOS+1)+AAcall(codon)
                                                        except:
                                                            mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] = "nts:"+ ref[3][orf]['nts'][mut1POS-1:mut1POS+2] + str(mut1POS)+MNP+'_AAs:'+ref[3][orf]["AAs"][AAPOS]+str(AAPOS+1)+AAcall(codon)
                                                    elif mut2POS - mut1POS == 1:
                                                        codon = ref[3][orf]['nts'][mut1POS-2]+mutorfs[orf]['muts'][i][1][-1]+mutorfs[orf]['muts'][i+1][1][-1]
                                                        MNP = "r" + mutorfs[orf]['muts'][i][1][-1]+mutorfs[orf]['muts'][i+1][1][-1]
                                                        skip = 1
                                                        try:
                                                            mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] += ", nts:"+ ref[3][orf]['nts'][mut1POS-2:mut1POS+1] + str(mut1POS)+MNP+'_AAs:'+ref[3][orf]["AAs"][AAPOS]+str(AAPOS+1)+AAcall(codon)
                                                        except:
                                                            mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] = "nts:"+ ref[3][orf]['nts'][mut1POS-2:mut1POS+1] + str(mut1POS)+MNP+'_AAs:'+ref[3][orf]["AAs"][AAPOS]+str(AAPOS+1)+AAcall(codon)
                                                    else:
                                                        AAinfo = singletCodon(mut1POS, mutorfs[orf]['muts'][i][1][-1], ref[3][orf]['nts'])
                                                        try:
                                                            mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] += ", nts:"+ str(mutorfs[orf]['muts'][i][2]) +'_AAs:'+ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]
                                                        except:
                                                            mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] = "nts:"+ str(mutorfs[orf]['muts'][i][2]) +'_AAs:'+ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]
                                                else:
                                                    AAinfo = singletCodon(mut1POS, mutorfs[orf]['muts'][i][1][-1], ref[3][orf]['nts'])
                                                    try:
                                                        mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] += ", nts:"+ str(mutorfs[orf]['muts'][i][2]) +'_AAs:'+ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]
                                                    except:
                                                        mutorfs[orf]['mutstrings'][mutorfs[orf]['muts'][i][1]] = "nts:"+ str(mutorfs[orf]['muts'][i][2]) +'_AAs:'+ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]



                            codonchecked = []
                            for mut in mutations:
                                newmutstring = mut
                                for orf in mutorfs:
                                    try:
                                        newmutstring += "(" + orf + "_" + mutorfs[orf]['mutstrings'][mut] + ")"
                                    except:
                                        pass
                                codonchecked.append(newmutstring)

                            mutations = codonchecked

                        mutations = " ".join(mutations)

                        if args.wgs == 0:
                            try:
                                seq_species[mutations] += reads_count
                            except:
                                seq_species[mutations] = reads_count
                        else:
                            try:
                                seq_species[str(POS)+' '+mutations+' '+str(POS+q_pars_pos)] += reads_count
                            except:
                                seq_species[str(POS)+' '+mutations+' '+str(POS+q_pars_pos)] = reads_count
                        if args.read == 1:
                            reads_fh.write(f"{readID}\t{mutations}\n")

                    for i in range(POS, POS+q_pars_pos): # update coverage
                        try:
                            coverage[i] += reads_count
                        except:
                            coverage[i] = reads_count
    if args.read == 1:
        reads_fh.close()
    sam_fh.close()
    # END SAM LINES
    print(f"End SAM parse for {samp}")

    if sam_read_count == 0:
        print(f"No Reads for {samp}")

    else:
        if args.seq == 1: # output the sequence
            seq_fh = open(samp+'_unique_seqs.tsv', "w")
            seq_fh.write(samp+"("+str(sam_read_count)+")\n")
            seq_fh.write("Unique Sequence\tCount\tAbundance\n")

            sorted_seq = sorted(seq_species, key=seq_species.__getitem__, reverse=True)
            for key in sorted_seq:
                if seq_species[key] >= args.min_count:
                    if (seq_species[key] / sam_read_count >= args.min_samp_abund) and args.wgs == 0:
                        seq_fh.write(f"{key}\t{seq_species[key]}\t{(seq_species[key]/sam_read_count):.3f}\n")
                    elif args.wgs == 1:
                        splitseqs = key.split()
                        cov = []
                        for x in range(int(splitseqs[0]), int(splitseqs[-1])):
                            cov.append(coverage[x])
                        min_cov = min(cov)
                        if (seq_species[key]/min_cov >= args.min_samp_abund):
                            seq_fh.write(f"{key}\t{seq_species[key]}\t{(seq_species[key]/min_cov):.3f}\n")
                else:
                    break

            seq_fh.close()
            # END SEQ OUT
            # print(f"End unqiue seq out for {samp}")

        if args.indel == 1 and len(indel_dict) > 0: # output indels, if there are any
            sorted_indels = sorted(indel_dict, key=indel_dict.__getitem__, reverse=True)
            indels_to_write = []
            for key in sorted_indels:
                if indel_dict[key] >= args.min_count:
                    if indel_dict[key] / sam_read_count >= args.min_samp_abund and args.wgs == 0:
                        indels_to_write.append(f"{key}\t{indel_dict[key]}\t{(indel_dict[key]/sam_read_count):.3f}\n")
                    elif args.wgs == 1:
                        indelPOS = ''
                        for c in key:
                            if c.isdigit():
                                indelPOS += c
                            else:
                                break
                        indelPOS = int(indelPOS)
                        if indel_dict[key] / coverage[indelPOS] >= args.min_samp_abund:
                            indels_to_write.append(f"{key}\t{indel_dict[key]}\t{(indel_dict[key] / coverage[indelPOS]):.3f}\n")
                else:
                    break
            if len(indels_to_write) > 0:
                indel_fh = open(samp+'_indels.tsv', "w")
                indel_fh.write(samp+"("+str(sam_read_count)+")\n")
                indel_fh.write("Indel\tCount\tAbundance\n")
                for indel_entry in indels_to_write:
                    indel_fh.write(indel_entry)
                indel_fh.close()
            # END INDEL OUT
            # print(f"End indel out for {samp}")

        if args.nt_call == 1: # out put nt calls
            ntcall_lines = {'line' : {},
                            'variant' : {}
                            }
            ntcall_fh = open(samp+'_nt_calls.tsv', "w")
            ntcall_fh.write(samp+"("+str(sam_read_count)+")\n")
            if args.ntvar == 1:
                ntcallv_fh = open(samp+'_nt_calls_varonly.tsv', "w")
                ntcallv_fh.write(samp+"("+str(sam_read_count)+")\n")
            sorted_POS = sorted(nt_call_dict_dict)
            if args.AAreport == 1:
                OrfPosDict = {}
                for orf in ref[3]:
                    orflength = 0
                    for rf in ref[3][orf]['reading frames']:
                        for i in range(rf[0],rf[1]+1):
                            orfPOS = i-rf[0]+1+orflength
                            try:
                                OrfPosDict[i].append([orf , orfPOS, ref[3][orf]['AAs'][(orfPOS-1)//3] ,(((orfPOS-1)//3)+1)])
                            except:
                                OrfPosDict[i] = [[orf , orfPOS, ref[3][orf]['AAs'][(orfPOS-1)//3] ,(((orfPOS-1)//3)+1)]]
                        orflength += rf[1] - rf[0] + 1

                ntcall_fh.write("Position\tref NT\tAAs\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                if args.ntvar == 1:
                    ntcallv_fh.write("Position\tref NT\tAAs\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                consensus = {}
                ORFmismatch = {}
                for POS in sorted_POS:
                    try:
                        total = coverage[POS]
                    except:
                        total = 0
                    if total >= (sam_read_count * args.ntabund) and total >= args.ntcover:
                        # AAinfo = singletCodon(POS, ref[1][POS-1], ref)
                        POS_calls = {}
                        for key in nt_call_dict_dict[POS]:
                            POS_calls[key] = nt_call_dict_dict[POS][key]
                        sorted_calls = sorted(POS_calls, key=POS_calls.__getitem__, reverse=True)

                        ntcall_lines['line'][POS] =(str(POS)+"\t"+ref[1][POS-1]+"\t")
                        try:
                            OrfPosDict[POS]
                        except:
                            pass # ntcall_lines['line'][POS] +=("\t")
                        else:
                            orfAAreports = []
                            for entry in OrfPosDict[POS]:
                                orfAAreports.append(entry[0] +"_nt:"+ str(entry [1]) +"_AA:"+ entry[2] + str(entry[3]))
                            ntcall_lines['line'][POS] +=(", ".join(orfAAreports))

                        ntcall_lines['line'][POS] +=("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                        ntcall_lines['line'][POS] +=("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[0]]))
                        ntcall_lines['line'][POS] +=(f"\t{(nt_call_dict_dict[POS][sorted_calls[0]]/total):.3f}")

                        consensus[POS] = sorted_calls
                        if consensus[POS][0] != ref[1][POS-1]:
                            try:
                                for entry in OrfPosDict[POS]:
                                    try:
                                        ORFmismatch[entry[0]]
                                    except:
                                        ORFmismatch[entry[0]] = {entry[1] : sorted_calls[0]}
                                    else:
                                        ORFmismatch[entry[0]][entry[1]] = sorted_calls[0]
                            except:
                                pass

                for POS in sorted_POS:
                    try:
                        total = coverage[POS]
                    except:
                        total = 0

                    try:
                        ntcall_lines['line'][POS]
                    except:
                        pass
                    else:
                        if consensus[POS][0] != ref[1][POS-1]:
                            ntcall_lines['variant'][POS] = 1
                            try:
                                OrfPosDict[POS]
                            except:
                                if (nt_call_dict_dict[POS][consensus[POS][1]] > args.min_count) and (nt_call_dict_dict[POS][consensus[POS][1]] /total > args.min_samp_abund):
                                    ntcall_lines['line'][POS] +=(f"\t\t\t{consensus[POS][1]}\t{nt_call_dict_dict[POS][consensus[POS][1]]}\t{(nt_call_dict_dict[POS][consensus[POS][1]]/total):.3f}")
                                    if nt_call_dict_dict[POS][consensus[POS][2]] > args.min_count and nt_call_dict_dict[POS][consensus[POS][2]] /total  > args.min_samp_abund:
                                        ntcall_lines['line'][POS] +=(f"\t\t{consensus[POS][2]}\t{nt_call_dict_dict[POS][consensus[POS][2]]}\t{(nt_call_dict_dict[POS][consensus[POS][2]]/total):.3f}")
                            else:
                                orfAAreports = [[],[]]
                                for ORFentry in OrfPosDict[POS]:
                                    OrfPOS = ORFentry[1]
                                    orf = ORFentry[0]
                                    mod = (OrfPOS)%3
                                    codon = ['n','n','n']

                                    if mod == 0:
                                        codon[2] = consensus[POS][0]
                                        try:
                                            codon[0] = ORFmismatch[orf][OrfPOS-2]
                                        except:
                                            codon[0] = ref[3][orf]['nts'][OrfPOS-3]
                                        try:
                                            codon[1] = ORFmismatch[orf][OrfPOS-1]
                                        except:
                                            codon[1] = ref[3][orf]['nts'][OrfPOS-2]

                                    elif mod == 2:
                                        codon[1] = consensus[POS][0]
                                        try:
                                            codon[2] = ORFmismatch[orf][OrfPOS+1]
                                        except:
                                            codon[2] = ref[3][orf]['nts'][OrfPOS]
                                        try:
                                            codon[0] = ORFmismatch[orf][OrfPOS-1]
                                        except:
                                            codon[0] = ref[3][orf]['nts'][OrfPOS-2]

                                    elif mod == 1:
                                        codon[0] = consensus[POS][0]
                                        try:
                                            codon[2] = ORFmismatch[orf][OrfPOS+2]
                                        except:
                                            codon[2] = ref[3][orf]['nts'][OrfPOS+1]
                                        try:
                                            codon[1] = ORFmismatch[orf][OrfPOS+1]
                                        except:
                                            codon[1] = ref[3][orf]['nts'][OrfPOS]
                                    orfAAreports[0].append(orf+"_"+AAcall("".join(codon)))
                                    orfAAreports[1].append(orf+"_"+singletCodon(OrfPOS, consensus[POS][0], ref[3][orf]['nts'])[1])
                                ntcall_lines['line'][POS] +=("\t"+", ".join(orfAAreports[0])+"\t"+", ".join(orfAAreports[1]))
                                if (nt_call_dict_dict[POS][consensus[POS][1]] > args.min_count) and (nt_call_dict_dict[POS][consensus[POS][1]] /total > args.min_samp_abund):
                                        ntcall_lines['line'][POS] +=(f"\t{consensus[POS][1]}\t{nt_call_dict_dict[POS][consensus[POS][1]]}\t{(nt_call_dict_dict[POS][consensus[POS][1]]/total):.3f}")
                                        orfAAreports = []
                                        for ORFentry in OrfPosDict[POS]:
                                            OrfPOS = ORFentry[1]
                                            orf = ORFentry[0]
                                            orfAAreports.append(orf+"_"+singletCodon(OrfPOS, consensus[POS][1], ref[3][orf]['nts'])[1])
                                        ntcall_lines['line'][POS] += ("\t"+", ".join(orfAAreports))
                                        if nt_call_dict_dict[POS][consensus[POS][2]] > args.min_count:
                                            if nt_call_dict_dict[POS][consensus[POS][2]] /total  > args.min_samp_abund:
                                                ntcall_lines['line'][POS] +=(f"\t{consensus[POS][2]}\t{nt_call_dict_dict[POS][consensus[POS][2]]}\t{(nt_call_dict_dict[POS][consensus[POS][2]]/total):.3f}")
                                                orfAAreports = []
                                                for ORFentry in OrfPosDict[POS]:
                                                    OrfPOS = ORFentry[1]
                                                    orf = ORFentry[0]
                                                    orfAAreports.append(orf+"_"+singletCodon(OrfPOS, consensus[POS][2], ref[3][orf]['nts'])[1])
                                                ntcall_lines['line'][POS] += ("\t"+", ".join(orfAAreports))


                        elif (nt_call_dict_dict[POS][consensus[POS][1]] >= args.min_count) and ((nt_call_dict_dict[POS][consensus[POS][1]] / total) >= args.min_samp_abund):
                            ntcall_lines['variant'][POS] = 1

                            ntcall_lines['line'][POS] +=("\t\t")
                            ntcall_lines['line'][POS] +=(f"\t{consensus[POS][1]}\t{nt_call_dict_dict[POS][consensus[POS][1]]}\t{(nt_call_dict_dict[POS][consensus[POS][1]]/total):.3f}")

                            try:
                                OrfPosDict[POS]
                            except:
                                if nt_call_dict_dict[POS][consensus[POS][2]] /total  > args.min_samp_abund:
                                    ntcall_lines['line'][POS] +=(f"\t\t{consensus[POS][2]}\t{nt_call_dict_dict[POS][consensus[POS][2]]}\t{(nt_call_dict_dict[POS][consensus[POS][2]]/total):.3f}")
                            else:
                                orfAAreports = []
                                for ORFentry in OrfPosDict[POS]:
                                    OrfPOS = ORFentry[1]
                                    orf = ORFentry[0]
                                    orfAAreports.append(orf+"_"+singletCodon(OrfPOS, consensus[POS][1], ref[3][orf]['nts'])[1])
                                ntcall_lines['line'][POS] += ("\t"+", ".join(orfAAreports))
                                if nt_call_dict_dict[POS][consensus[POS][2]] > args.min_count:
                                    if nt_call_dict_dict[POS][consensus[POS][2]] /total  > args.min_samp_abund:
                                        ntcall_lines['line'][POS] +=(f"\t{consensus[POS][2]}\t{nt_call_dict_dict[POS][consensus[POS][2]]}\t{(nt_call_dict_dict[POS][consensus[POS][2]]/total):.3f}")
                                        orfAAreports = []
                                        for ORFentry in OrfPosDict[POS]:
                                            OrfPOS = ORFentry[1]
                                            orf = ORFentry[0]
                                            orfAAreports.append(orf+"_"+singletCodon(OrfPOS, consensus[POS][2], ref[3][orf]['nts'])[1])
                                        ntcall_lines['line'][POS] += ("\t"+", ".join(orfAAreports))

                        ntcall_lines['line'][POS] +=("\n")
                for POS in ntcall_lines['line']:
                    ntcall_fh.write(ntcall_lines['line'][POS])
                    if args.ntvar == 1:
                        try:
                            ntcall_lines['variant'][POS]
                        except:
                            pass
                        else:
                            ntcallv_fh.write(ntcall_lines['line'][POS])
                if args.ntvar == 1:
                        ntcallv_fh.close()
            else:
                ntcall_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")
                if args.ntvar == 1:
                    ntcallv_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")

                for POS in sorted_POS:
                    try:
                        total = coverage[POS] # sum(nt_call_dict_dict[POS].values())
                    except:
                        total = 0
                    if total >= (sam_read_count * args.ntabund) and total >= args.ntcover:
                        POS_calls = {}
                        for key in nt_call_dict_dict[POS]:
                            POS_calls[key] = nt_call_dict_dict[POS][key]
                        sorted_calls = sorted(POS_calls, key=POS_calls.__getitem__, reverse=True)

                        ntcall_fh.write(str(POS)+"\t"+ref[1][POS-1])
                        ntcall_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                        ntcall_fh.write("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[0]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[0]]/total):.3f}")

                        if sorted_calls[0] != ref[1][POS-1]:
                            if args.ntvar == 1:
                                ntcallv_fh.write(str(POS)+"\t"+ref[1][POS-1])
                                ntcallv_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                                ntcallv_fh.write("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[0]]))
                                ntcallv_fh.write(f"\t{(nt_call_dict_dict[POS][sorted_calls[0]]/total):.3f}")
                            if (nt_call_dict_dict[POS][sorted_calls[1]] > args.min_count) and (nt_call_dict_dict[POS][sorted_calls[1]] / total > args.min_samp_abund):
                                ntcall_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                                if sorted_calls[1] != ref[1][POS-1] and args.ntvar == 1:
                                    ntcallv_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                                if (nt_call_dict_dict[POS][sorted_calls[2]] > args.min_count) and (nt_call_dict_dict[POS][sorted_calls[2]] /total > args.min_samp_abund):
                                    ntcall_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                                    if sorted_calls[2] != ref[1][POS-1] and args.ntvar == 1:
                                        ntcallv_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                            if args.ntvar == 1:
                                ntcallv_fh.write("\n")

                        elif (nt_call_dict_dict[POS][sorted_calls[1]] > args.min_count) and (nt_call_dict_dict[POS][sorted_calls[1]] /total > args.min_samp_abund):
                            if args.ntvar == 1:
                                ntcallv_fh.write(str(POS)+"\t"+ref[1][POS-1])
                                ntcallv_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                                ntcallv_fh.write("\t"+str(total)+"\t\t")
                                ntcallv_fh.write(f"\t")
                                ntcallv_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                            ntcall_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                            if (nt_call_dict_dict[POS][sorted_calls[2]] > args.min_count) and (nt_call_dict_dict[POS][sorted_calls[2]] /total > args.min_samp_abund):
                                ntcall_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                                if sorted_calls[2] != ref[1][POS-1] and args.ntvar == 1:
                                    ntcallv_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                            if args.ntvar == 1:
                                ntcallv_fh.write("\n")

                        ntcall_fh.write("\n")

            ntcall_fh.close()
            if args.ntvar == 1:
                ntcallv_fh.close()
            # END NT CALL OUT
            # print(f"End nt call out for {samp}")
        if args.covar == 1: # output covariants
            testtrack = 0
            combinations = {}
            for sequence in seq_species:
                if args.wgs == 0:
                    singles = sequence.split()
                else:
                    singles = (sequence.split())[1:-1]
                if len(singles) <= args.max_dist and singles[0] != 'Ref':
                    for combo in getCombos(singles, args.max_covar):
                        if not combo in combinations:
                            combinations[combo] = seq_species[sequence]
                        else:
                            combinations[combo] += seq_species[sequence]

            covar_fh = open(samp+'_covars.tsv', "w")
            covar_fh.write(samp+"("+str(sam_read_count)+")\n")
            covar_fh.write("Co-Variants\tCount\tAbundance\n")
            sortedcombos = sorted(combinations, key=combinations.__getitem__, reverse=True)
            for key in sortedcombos:
                if (combinations[key] >= args.min_count):
                    if (combinations[key] / sam_read_count >= args.min_samp_abund) and args.wgs == 0:
                        covar_fh.write(key+"\t"+str(combinations[key])+"\t"+f"{(combinations[key]/sam_read_count):.3f}\n")
                    elif args.wgs == 1:
                        coveragepercent = 0
                            # print(key)
                        splitcombos = key.split()
                        if len(splitcombos) == 1:
                            coveragePOS = ''
                            for c in key:
                                if c.isdigit():
                                    coveragePOS += c
                                else:
                                    break
                            coveragepercent = combinations[key] / coverage[int(coveragePOS)]
                        else:
                            startcovPOS = ''
                            for c in splitcombos[0]:
                                if c.isdigit():
                                    startcovPOS += c
                                else:
                                    break
                            endcovPOS = ''
                            for c in splitcombos[-1]:
                                if c.isdigit():
                                    endcovPOS += c
                                else:
                                    break
                            coveragevals = []
                            for i in range(int(startcovPOS), int(endcovPOS)+1):
                                coveragevals.append(coverage[i])
                            mincov = min(coverval for coverval in coveragevals)
                            coveragepercent = combinations[key] / mincov
                        if coveragepercent >= args.min_samp_abund:
                            covar_fh.write(f"{key}\t{combinations[key]}\t{coveragepercent:.3f}\n")

                    # \t{coveragepercent:.3f}
            covar_fh.close()


            # END COVAR OUT
            # print(f"End covar out for {samp}")

def cvdeconv(args, samp, covardict, seqdict): # covar deconvolution process

    passedseqs = {}
    preservedseqs = {}
    for seq in seqdict: # pass check for actual : expected abundance
        if seq != 'total' and seq != 'singles':

            splitseq = seq.split(' ')
            abund = 1
            for sing in splitseq:
                try:
                    abund = abund * (covardict[sing] / covardict['total'])
                except:
                    abund = abund * (seqdict[seq] / seqdict['total'])

            try:
                covarabund = covardict[seq]/covardict['total']
            except:
                covarabund = seqdict[seq]/seqdict['total']
                covardict[seq] = seqdict[seq]

            if covarabund >= args.autopass:
                preservedseqs[seq] = max(1, args.beta, (covarabund / abund))
            elif covarabund >= abund * args.beta:
                passedseqs[seq] = covarabund / abund
            elif len(seq.split(' ')) == 1:
                passedseqs[seq] = max(1, args.beta)

    if args.min_samp_abund < 1:
        min_count = args.min_samp_abund * covardict['total']
    else:
        min_count = args.min_samp_abund

    if args.pass_out == 1: # write passed covars to file if enabled
        fh_pass = open(samp+"_covar_pass.tsv", "w")
        fh_pass.write(f"{samp}({covardict['total']})\tCount\tAbundance\tPass Ratio\n")
        for seq in preservedseqs:
            if covardict[seq] >= min_count:
                fh_pass.write(f"{seq}\t{covardict[seq]}\t{(covardict[seq]/covardict['total']):.3f}\t{preservedseqs[seq]}*\n")
        for seq in passedseqs:
            if covardict[seq] >= min_count:
                fh_pass.write(f"{seq}\t{covardict[seq]}\t{(covardict[seq]/covardict['total']):.3f}\t{passedseqs[seq]}\n")
        fh_pass.close()

    # sort passed covars
    lensortedpassed = sorted(passedseqs, key = lambda key : len(key.split(' ')), reverse=True)
    ratiolensortedpassed = sorted(lensortedpassed, key = lambda key : passedseqs[key], reverse = True)
    sortedsingles = sorted(covardict['singles'], key = covardict['singles'].__getitem__)
    # print(lensortedpassed)
    # print(ratiolensortedpassed)
    deconved = {}
    for seq in ratiolensortedpassed: # reassign counts
        # print(seq)
        singles = seq.split(' ')
        first = 0
        rem_count = 0
        for sing in sortedsingles:
            # print(sing+' '+str(covardict['singles'][sing]))
            if sing in singles:
                if covardict['singles'][sing] > 0:
                    if first == 0:
                        first = 1
                        rem_count = covardict['singles'][sing]
                        covardict['singles'][sing] = 0
                        deconved[seq] = rem_count
                    else:
                        covardict['singles'][sing] = covardict['singles'][sing] - rem_count
                else:
                    break
        sortedsingles = sorted(covardict['singles'], key = covardict['singles'].__getitem__)

    sortedpreserved = sorted(preservedseqs, key = lambda key : covardict[key])

    for seq in sortedpreserved:
        # print(seq)
        singles = seq.split(' ')
        first = 0
        rem_count = 0
        for sing in sortedsingles:
            # print(sing+' '+str(covardict['singles'][sing]))
            if sing in singles:
                if covardict['singles'][sing] > 0:
                    if first == 0:
                        first = 1
                        rem_count = covardict['singles'][sing]
                        covardict['singles'][sing] = 0
                        deconved[seq] = rem_count
                    else:
                        covardict['singles'][sing] = covardict['singles'][sing] - rem_count
                else:
                    break
        sortedsingles = sorted(covardict['singles'], key = covardict['singles'].__getitem__)


    newtotal = sum(deconved.values())
    fh_deconv = open(samp+"_covar_deconv.tsv", "w")
    fh_deconv.write(f"{samp}({covardict['total']}) | ({newtotal})\tCount\tAbundance\n")
    sorted_deconved = sorted(deconved, key = deconved.__getitem__, reverse = True)
    for seq in sorted_deconved: # write deconv
        if deconved[seq] >= min_count:
            fh_deconv.write(f"{seq}\t{deconved[seq]}\t{(deconved[seq]/newtotal):.3f}\n")
    fh_deconv.close()

    print(f"End covar deconv out for {samp}") # END COVAR DECONV OUT

    return()

def dechim(args, seqs): # processing sequence dictionary to remove chimeras


    total = seqs['total']
    del seqs['total']
    sorted_seqs = sorted(seqs, key=seqs.__getitem__) # sort sequences by abundance, least to greatest
    chimeras = []
    for seq in sorted_seqs:
        pot_chim = ['Reference']+seq.split()+['Reference']
        chim_halves = []
        for i in range(0, len(pot_chim)-1): # create dimera halves
            chim_halves.append([pot_chim[:i+1], pot_chim[i+1:]])
        parent_pairs = []

        for dimera in chim_halves: # check for potential parents

            pot_left = []
            pot_rt = []
            lft_len = len(dimera[0])
            rt_len = len(dimera[1])
            for pseq in sorted_seqs:
                if not seq == pseq:
                    if (seqs[pseq] >= (seqs[seq] * args.foldab)):
                        pot_par = ['Reference']+pseq.split()+['Reference']
                        if dimera[0] == pot_par[:lft_len]:
                            pot_left.append(pot_par)
                        if ((len(pot_par) >= rt_len) and (dimera[1] == pot_par[(len(pot_par)-rt_len):])):
                            pot_rt.append(pot_par)

            if (len(pot_left) > 0 and len(pot_rt) > 0 ):
                for left_par in pot_left: # check potential parents' pairing
                    for rt_par in pot_rt:
                        if left_par != rt_par:
                            left_break = left_par[lft_len]
                            rt_break = rt_par[(len(rt_par)-rt_len)-1]
                            if left_break == 'Reference' or rt_break == 'Reference':
                                parent_pairs.append([' '.join(left_par[1:-1]), ' '.join(rt_par[1:-1])])
                            else:
                                left_break_POS = ''
                                for c in left_break:
                                    if c.isdigit():
                                        left_break_POS += c
                                    else:
                                        if left_break_POS:
                                            break

                                rt_break_POS = ''
                                for c in rt_break:
                                    if c.isdigit():
                                        rt_break_POS += c
                                    else:
                                        if rt_break_POS:
                                            break

                                if int(left_break_POS) > int(rt_break_POS):
                                    parent_pairs.append([' '.join(left_par[1:-1]), ' '.join(rt_par[1:-1])])

        par_tot_abund = 0
        pair_probs = []
        for parents in parent_pairs: # calc 'expected' abundance
            pair_prob = (seqs[parents[0]] / total) * (seqs[parents[1]] / total)
            par_tot_abund += pair_prob
            pair_probs.append(pair_prob)

        recomb_count = par_tot_abund * total

        if not seqs[seq] >= recomb_count * args.alpha: # chimera check
            redist_count = float(seqs[seq])
            chimeras.append(seq)
            seqs[seq] = 0
            if args.redist == 1: # redist counts of chimera
                toadd = {}
                for i in range(0, len(parent_pairs)):
                    counts_to_redist = (redist_count * (pair_probs[i]/par_tot_abund))/2
                    seqs[parent_pairs[i][0]] += counts_to_redist
                    seqs[parent_pairs[i][1]] += counts_to_redist



    for chim in chimeras: # remove chimeras
        del seqs[chim]

    total = sum(seqs.values())


    # total = sum(seqs.values)
    seqs['total'] = total

    return(seqs)

def chimrm(args, samp, seqs): # chimera removed process

    pre_len = len(seqs)
    inf_loop_shield = 0
    while True: # send sequences for chimera removal while chimeras are still found
        dechim(args, seqs)
        post_len = len(seqs)
        inf_loop_shield += 1
        # print(f"{inf_loop_shield} {pre_len} {post_len}")
        if post_len >= pre_len:
            break
        if inf_loop_shield > args.max_cycles:
            break
        pre_len = len(seqs)

    total = seqs['total']
    del seqs['total']
    if args.min_samp_abund < 1:
        min_count = args.min_samp_abund * total
    else:
        min_count = args.min_samp_abund
    sorted_seqs = sorted(seqs, key=seqs.__getitem__, reverse=True)
    fh_dechime = open(f"{samp}_a{args.alpha}f{args.foldab}rd{args.redist}_chim_rm.tsv",'w')
    fh_dechime.write(f"{samp}({int(total)})\n")
    fh_dechime.write("Sequences\tCount\tAbundance\n")
    for seq in seqs: # write chim_rm seqs
        abund = seqs[seq]/total
        if seqs[seq] >= min_count:
            fh_dechime.write(f"{seq}\t{round(seqs[seq])}\t{abund:.3f}\n")

    fh_dechime.close()
    print(f"End chim_rm out for {samp}") # END CHIM RM DECONV OUT
    return()

def chimproc(args, samp):
    if args.deconv == 1:
        in_covars = {}
        in_seqs = {}
        try:
            seqin_file = open(samp+'_unique_seqs.tsv', 'r')
            for line in seqin_file:
                lineparts = line.strip("\n\r").split("\t")
                try:
                    lineparts[1]
                except:
                    in_seqs = {'total' : int(lineparts[0].split("(")[1][0:-1])}
                else:
                    if lineparts[1] != 'Count':
                        if float(lineparts[2]) >= args.chim_in_abund:
                            in_seqs[lineparts[0]] = float(lineparts[1])
            seqin_file.close()
        except:
            print(f"Reading of {samp}_unique_seqs.tsv failed")

        try:
            covin_file = open(samp+'_covars.tsv', 'r')
            for line in covin_file:
                lineparts = line.strip("\n\r").split("\t")
                try:
                    lineparts[1]
                except:
                    in_covars = {'total' : int(lineparts[0].split("(")[1][0:-1]),
                                              'singles' : {}
                                              }
                else:
                    if lineparts[1] != 'Count':
                        if float(lineparts[2]) >= args.chim_in_abund:
                            in_covars[lineparts[0]] = int(lineparts[1])
                            if len(lineparts[0].split(' ')) == 1:
                                in_covars['singles'][lineparts[0]] = int(lineparts[1])
            covin_file.close()
            if in_covars and in_seqs:
                cvdeconv(args, samp, in_covars, in_seqs)
        except:
            print(f"Reading of {samp}_covars.tsv failed")

        if args.chim_rm == 1:
            chimrm(args, samp, in_seqs)

    elif args.chim_rm == 1:
            in_covars = {}
            in_seqs = {}
            try:
                seqin_file = open(samp+'_unique_seqs.tsv', 'r')
                for line in seqin_file:
                    lineparts = line.strip("\n\r").split("\t")
                    try:
                        lineparts[1]
                    except:
                        in_seqs = {'total' : int(lineparts[0].split("(")[1][0:-1])}
                    else:
                        if lineparts[1] != 'Count':
                            if float(lineparts[2]) >= args.chim_in_abund:
                                in_seqs[lineparts[0]] = float(lineparts[1])
                seqin_file.close()
                chimrm(args, samp, in_seqs)
            except:
                print(f"Failed to Process {samp}_unique_seqs.tsv")

def main():

    args = arg_parser() # getting command line arguments

    #if args.sams == 1:

    if args.ref:
        ref = get_ref(args) # get the reference ID and sequence from the FASTA file
        # args.ref.close()
        if ref[1] == '':
            print('Reference not recognized as a Fasta or Genebank format, skipping SAM parsing')
        else:

            # print(ref[3])
            # collect SAM files to process, either from the command line or the working directory
            SAMs = []
            try:
                args.Sam_files[0]
            except:
                for file in os.listdir(os.getcwd()):
                    if (file.lower()).endswith('.sam'):
                        SAMs.append(file)
            else:
                for files in args.Sam_files:
                    for file in files:
                        if os.path.isfile(file):
                            SAMs.append(file)
                        else:
                            print(f"Can't find {file}, skipping")

            args.ref = ''
            if ref[2] == 'fasta':
                with Pool(processes=args.mp) as pool:
                    pool.starmap(faSAMparse, zip(itertools.repeat(args), itertools.repeat(ref), SAMs))
            elif ref[2] == 'gb':
                with Pool(processes=args.mp) as pool:
                    pool.starmap(gbSAMparse, zip(itertools.repeat(args), itertools.repeat(ref), SAMs))
            print(f"End Sam Parsing Output")
    else:
        print('No reference provided, skipping SAM parsing')


    seq_files = []
    # Begin chimera removal if enabled
    if args.chim_rm == 1 or args.deconv == 1:
        for file in os.listdir(os.getcwd()):
            if file.endswith('_seqs.tsv'): # get unique sequences for chimera removal
                seq_files.append(file[0:-16])
        with Pool(processes=args.mp) as pool:
            pool.starmap(chimproc, zip(itertools.repeat(args), seq_files))

# begin collection of sample outputs
    if args.collect == 1:
        covar_dict_dict = {}
        seq_dict_dict = {}
        deconv_dict_dict = {}
        pass_dict_dict = {}
        cr_dict_dict = {}

        all_covars = {}
        all_seqs = {}
        all_deconv = {}
        all_pass = {}
        all_cr = {}

        sample_line = ''

        for file in os.listdir(os.getcwd()):
            if file.endswith('_covars.tsv'):
                try:
                    samp=open(file, "r")
                except:
                    print("can't open "+file)
                else:
                    for line in samp:
                        splitline = line.strip("\n\r").split("\t")
                        try:
                            splitline[1]
                        except:
                            sample_line = splitline[0]
                            covar_dict_dict[sample_line] = {}

                        else:
                            if not splitline[1] == 'Count':
                                if float(splitline[2]) >= args.min_col_abund:
                                    covar_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_covars[splitline[0]] = 1
                    samp.close()


            if file.endswith('_seqs.tsv'):
                try:
                    samp=open(file, "r")
                except:
                    print("can't open "+file)
                else:
                    for line in samp:
                        splitline = line.strip("\n\r").split("\t")
                        try:
                            splitline[1]
                        except:
                            sample_line = splitline[0]
                            seq_dict_dict[sample_line] = {}

                        else:
                            if not splitline[1] == 'Count':
                                if float(splitline[2]) >= args.min_col_abund:
                                    seq_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_seqs[splitline[0]] = 1
                    samp.close()

            if file.endswith('_deconv.tsv'):
                try:
                    samp=open(file, "r")
                except:
                    print("can't open "+file)
                else:
                    for line in samp:
                        splitline = line.strip("\n\r").split("\t")
                        try:
                            splitline[1]
                        except:
                            sample_line = splitline[0]
                            deconv_dict_dict[sample_line] = {}

                        else:
                            if splitline[1] == 'Count':
                                sample_line = splitline[0]
                                deconv_dict_dict[sample_line] = {}
                            else:
                                if float(splitline[2]) >= args.min_col_abund:
                                    deconv_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_deconv[splitline[0]] = 1
                    samp.close()

            if file.endswith('_pass.tsv'):
                try:
                    samp=open(file, "r")
                except:
                    print("can't open "+file)
                else:
                    for line in samp:
                        splitline = line.strip("\n\r").split("\t")
                        try:
                            splitline[1]
                        except:
                            sample_line = splitline[0]
                            pass_dict_dict[sample_line] = {}

                        else:
                            if splitline[1] == 'Count':
                                sample_line = splitline[0]
                                pass_dict_dict[sample_line] = {}
                            else:
                                if float(splitline[2]) >= args.min_col_abund:
                                    pass_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_pass[splitline[0]] = 1
                    samp.close()

            if file.endswith('_chim_rm.tsv'):
                try:
                    samp=open(file, "r")
                except:
                    print("can't open "+file)
                else:
                    for line in samp:
                        splitline = line.strip("\n\r").split("\t")
                        try:
                            splitline[1]
                        except:
                            sample_line = splitline[0]
                            cr_dict_dict[sample_line] = {}

                        else:
                            if not splitline[1] == 'Count':
                                if float(splitline[2]) >= args.min_col_abund:
                                    cr_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_cr[splitline[0]] = 1
                    samp.close()

        if len(covar_dict_dict) > 0:
            if args.colID == '':
                Abund_Poly_fh = open('Collected_Covariances.tsv',"w")
            else:
                Abund_Poly_fh = open(args.colID+'_Collected_Covriances.tsv',"w")
            sorted_covars = sorted(all_covars)
            Abund_Poly_fh.write("\t")
            for sampline in covar_dict_dict:
                Abund_Poly_fh.write(sampline+"\t\t")
            Abund_Poly_fh.write("\nVariants\t")
            for sampline in covar_dict_dict:
                Abund_Poly_fh.write("Count\tAbundance\t")
            Abund_Poly_fh.write("\n")
            for covar in sorted_covars:
                Abund_Poly_fh.write(covar+"\t")
                for sample in covar_dict_dict:
                    try:
                        Abund_Poly_fh.write(covar_dict_dict[sample][covar][0]+"\t"+covar_dict_dict[sample][covar][1]+"\t")
                    except:
                        Abund_Poly_fh.write("\t\t")

                Abund_Poly_fh.write("\n")
            Abund_Poly_fh.close()
        else:
            print('No covar files found')

        if len(seq_dict_dict) > 0:
            if args.colID == '':
                Col_Seq_fh = open('Collected_Unique_Seqs.tsv',"w")
            else:
                Col_Seq_fh = open(args.colID+'_Collected_Unique_Seqs.tsv',"w")
            sorted_seqs = sorted(all_seqs)
            Col_Seq_fh.write("\t")
            for sampline in seq_dict_dict:
                Col_Seq_fh.write(sampline+"\t\t")
            Col_Seq_fh.write("\nUnique Sequences\t")
            for sampline in seq_dict_dict:
                Col_Seq_fh.write("Count\tAbundance\t")
            Col_Seq_fh.write("\n")
            for seq in sorted_seqs:
                Col_Seq_fh.write(seq+"\t")
                for sample in seq_dict_dict:
                    try:
                        Col_Seq_fh.write(seq_dict_dict[sample][seq][0]+"\t"+seq_dict_dict[sample][seq][1]+"\t")
                    except:
                        Col_Seq_fh.write("\t\t")

                Col_Seq_fh.write("\n")
            Col_Seq_fh.close()
        else:
            print('No seq files found')

        if len(deconv_dict_dict) > 0:
            if args.colID == '':
                Col_deconv_fh = open('Collected_Covar_Deconv.tsv',"w")
            else:
                Col_deconv_fh = open(args.colID+'_Collected_Covar_Deconv.tsv',"w")
            sorted_deconvs = sorted(all_deconv)
            Col_deconv_fh.write("\t")
            for sampline in deconv_dict_dict:
                Col_deconv_fh.write(sampline+"\t\t")
            Col_deconv_fh.write("\nCovariant\t")
            for sampline in deconv_dict_dict:
                Col_deconv_fh.write("Count\tAbundance\t")
            Col_deconv_fh.write("\n")
            for deconv in sorted_deconvs:
                Col_deconv_fh.write(deconv+"\t")
                for sample in deconv_dict_dict:
                    try:
                        Col_deconv_fh.write(deconv_dict_dict[sample][deconv][0]+"\t"+deconv_dict_dict[sample][deconv][1]+"\t")
                    except:
                        Col_deconv_fh.write("\t\t")

                Col_deconv_fh.write("\n")
            Col_deconv_fh.close()
        else:
            print('No deconv files found')

        if len(pass_dict_dict) > 0:
            if args.colID == '':
                Col_pass_fh = open('Collected_Covar_Pass.tsv',"w")
            else:
                Col_pass_fh = open(args.colID+'_Collected_Covar_Pass.tsv',"w")
            sorted_pass = sorted(all_pass)
            Col_pass_fh.write("\t")
            for sampline in pass_dict_dict:
                Col_pass_fh.write(sampline+"\t\t")
            Col_pass_fh.write("\nCovariant\t")
            for sampline in pass_dict_dict:
                Col_pass_fh.write("Count\tAbundance\t")
            Col_pass_fh.write("\n")
            for passed in sorted_pass:
                Col_pass_fh.write(passed+"\t")
                for sample in pass_dict_dict:
                    try:
                        Col_pass_fh.write(pass_dict_dict[sample][passed][0]+"\t"+pass_dict_dict[sample][passed][1]+"\t")
                    except:
                        Col_pass_fh.write("\t\t")

                Col_pass_fh.write("\n")
            Col_pass_fh.close()
        else:
            print('No pass files found')

        if len(cr_dict_dict) > 0:
            if args.colID == '':
                Col_cr_fh = open('Collected_Chimeras_Removed.tsv',"w")
            else:
                Col_cr_fh = open(args.colID+'_Collected_Chimeras_Removed.tsv',"w")
            sorted_crs = sorted(all_cr)
            Col_cr_fh.write("\t")
            for sampline in cr_dict_dict:
                Col_cr_fh.write(sampline+"\t\t")
            Col_cr_fh.write("\nUnique Sequence\t")
            for sampline in cr_dict_dict:
                Col_cr_fh.write("Count\tAbundance\t")
            Col_cr_fh.write("\n")
            for cr in sorted_crs:
                Col_cr_fh.write(cr+"\t")
                for sample in cr_dict_dict:
                    try:
                        Col_cr_fh.write(cr_dict_dict[sample][cr][0]+"\t"+cr_dict_dict[sample][cr][1]+"\t")
                    except:
                        Col_cr_fh.write("\t\t")

                Col_cr_fh.write("\n")
            Col_cr_fh.close()
        else:
            print('No chim_rm files found')


if __name__ == '__main__':
    main()