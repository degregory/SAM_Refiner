#!/bin/env python3

# Writen by Devon Gregory with assistance from Christopher Bottoms
# University of Missouri
# Distributed under GNU GENERAL PUBLIC LICENSE v3
# last editted on 20220501

import os
import sys
import argparse
import itertools
from multiprocessing import Process, Pool

bc = False
try:
    import pysam
    pysam.set_verbosity(0)
    bc = True
except:
    print("Failed to import pysam.  Processing of bams/crams skipped")


"""
To Do:
replace none atcgn- nt calls function
handle snp '-' deletion in MNP processing
Clean up / polish / add comments
add --verbose --quiet options
"""

def arg_parser():
    """
    Called to get the arguments passed by the command line and process them
    Functionality:
    uses argparse module to set up argument values based on the command line then does some conflict checking
    to make sure incompatible values aren't assigned or values aren't out of usable bounds
    Returns stored argument values
    """

    parser = argparse.ArgumentParser(
        description='process Sam files for variant information'
    )

    parser.add_argument(
        '-r', '-reference',
        type=argparse.FileType('r'),
        dest='ref',
        help='reference fasta or genbank file.  Only chimera removal and collections will be performed if omitted.'
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
        help='Enable/Disable (1/0) use of counts in sequence IDs, default enabled (--use_count 1) (default: 1)'
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
        '--covar_tile_coverage',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) using tiles covering positions instead of minimum nt coverage to calculate abundance of covariants (--covar_tile_coverage), require --wgs 1 (default: 0)'
    )
    parser.add_argument(
        '--AAreport',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) amino acid reporting, default enabled (--AAreport 1) (default: 1)'
    )
    parser.add_argument(
        '--AAcodonasMNP',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) reporting multiple nt changes in a single codon as one polymorphism, default enabled (--AAcodonasMNP 1), requires AAreport enabled (default: 1)'
    )
    parser.add_argument(
        '--AAcentered',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) amino acid centered seq and covar outputs for .gb processing (--AAcentered 0), requires AAreport enabled (default: 0)'
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
        help='threshold for a sequence to automatically pass the covar pass checking (default: 0.3)'
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
        help='Enable/Disable (1/0) collection step, default enabled (--collect 1) (default: 1)'
    )
    parser.add_argument(
        '--read',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) reads output, default disabled (--read 0) (default: 0)'
    )
    parser.add_argument(
        '--nt_call',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) nt_call output, default enabled (--nt_call 1) (default: 1)'
    )
    parser.add_argument(
        '--ntvar',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) nt_call output, default enabled (--ntvar 1) (default: 0)'
    )
    parser.add_argument(
        '--indel',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) indel output, default enabled (--indel 1) (default: 1)'
        )
    parser.add_argument(
        '--seq',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) unique seq output, default enabled (--seq 1) (default: 1)'
    )
    parser.add_argument(
        '--covar', type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) covar output, default enabled (--covar 1) (default: 1)'
    )
    parser.add_argument(
        '--pass_out',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) covar_pass output, default disabled (--pass_out 0) (default: 0)'
    )
    parser.add_argument(
        '--chim_rm',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) chim_rm output, default enabled (--chim_rm 1) (default: 1)'
    )
    parser.add_argument(
        '--deconv',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) covar deconv, default enabled (--deconv 1) (default: 1)'
    )
    parser.add_argument(
        '--wgs',
        type=int,
        default=0,
        choices=[0, 1],
        help='Enable/Disable (1/0) covar deconv, default enabled (--wgs 1)(default: 0)'
    )
    parser.add_argument(
        '--mp',
        type=int,
        default=4,
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

    if args.min_samp_abund < 0 or args.min_samp_abund >= 1:
        print(f"--min_samp_abund must be non-negative and < 1, defaulting to .001")
        args.min_samp_abund=0.001

    if args.min_col_abund < 0 or args.min_col_abund >= 1:
        print(f"--min_col_abund must be non-negative and < 1, defaulting to .01")
        args.min_col_abund=0.01

    if args.ntabund < 0 or args.ntabund >= 1:
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

    if args.chim_in_abund < 0 or args.chim_in_abund >= 1:
        print(f"--chim_in_abund must be non-negative, defaulting to 0.001")
        args.chim_in_abund = 0.001

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

def get_ref(args):
    """
    Called to get reference data from file provided in the args
    Parameters:
    args - argument values
    Functionality:
    From the reference provide, attempts to obtain reference ID and NT sequence, and AA sequence(s) if AAreport is enabled.
    For fasta formatted files, only the first entry is parsed.
    For genbank formatted files, CDS AA sequences will be pulled along with their gene ID

    Returns reference ID, nt sequence, file type and amino acid sequences encoded by the nt seq for fasta files or the CDS reported by genbank files
    """
    n=0
    ref_id = ''
    ref_type = ''
    ref_seq = ''
    ref_orfs = {}
    if args.ref:
        ref = args.ref
        first_line = ref.readline()
        if first_line.startswith('>'):
            ref_type = 'fasta'
            n+=1
            ref_id = first_line[1:].strip("\n\r").split(" ")[0]
            for line in ref:
                if line.startswith('>'):
                    n+=1
                    if n > 1:
                        break
                    ref_id = line[1:].strip("\n\r")
                elif n == 1:
                    ref_seq = ref_seq + line.strip("\n\r")
            ref_seq = ref_seq.upper()
            ref_prot = ''
            if args.AAreport == 1:
                for x in range(0, (len(ref_seq))//3):
                    amino_acid = aa_call(ref_seq[x*3]+ref_seq[x*3+1]+ref_seq[x*3+2])
                    ref_prot = ref_prot + amino_acid
                if (len(ref_seq))%3 != 0:
                    ref_prot = ref_prot + '?'
                ref_orfs = [ref_id, ref_prot]
        elif first_line.upper().startswith("LOCUS"):
            ref_type = 'gb'
            collect = "Null"
            orfs = {}
            reading_frames = []
            trans = 0
            nts = ""
            for line in ref:
                if collect == "Null":
                    split_line = line.strip("\n\r").split(" ")

                    if split_line[0].upper() == "VERSION":
                        ref_id = split_line[-1]
                    elif "CDS" in line:
                        collect = "CDS"
                        if "join" in line:
                            startstops = split_line[-1].strip("join()").split(",")
                            for startstop in startstops:
                                reading_frames.append([int(startstop.split(".")[0]) , int(startstop.split(".")[-1])])

                        else:
                            startstop = split_line[-1].split(".")
                            reading_frames.append([int(startstop[0]) , int(startstop[-1])])

                    if split_line[0].upper() == "ORIGIN":
                        collect = "SEQ"
                elif collect == "CDS":
                    if "gene=" in line:
                        gene = line.split("=")[1].strip('"\n\r')
                    elif "product=" in line:
                        product = line.split("=")[1].strip('"\n\r')

                    elif "/translation=" in line:
                        try:
                            gene
                        except:
                            gene = 'unlabeled'
                        try:
                            product
                        except:
                            product = gene

                        orf_id = product.replace(' ', '_')
                        n = 1
                        if orf_id in orfs:
                            new_id = orf_id + '.' + str(n)
                            while new_id in orfs:
                                n += 1
                                new_id = orf_id + '.' + str(n)
                            orf_id = new_id
                        orfs[orf_id] = { "reading frames" : reading_frames
                                    }
                        reading_frames = []
                        orfs[orf_id]["AAs"] = line.strip("\r\n").split('"')[1]

                        if not line.strip("\r\n")[-1] == '"':
                            trans = 1
                        else:
                            orfs[orf_id]["AAs"] += "*"
                            gene = ""
                            orf_id = ""
                            collect = "Null"

                    elif trans == 1:
                        orfs[orf_id]["AAs"] = orfs[orf_id]["AAs"] + line.strip(' "\n\r')
                        if line.strip("\r\n")[-1] == '"':
                            orfs[orf_id]["AAs"] += "*"
                            trans = 0
                            gene = ""
                            orf_id = ""
                            collect = "Null"

                elif collect == "SEQ":
                    if not "//" in line:
                        ntsline = ""
                        for c in line:
                            if c.isalpha():
                                ntsline += c
                        nts += ntsline
            if args.AAreport == 1:
                for gene in orfs:
                    orfnts = ''
                    for rf in orfs[gene]['reading frames']:
                        orfnts += nts[rf[0]-1:rf[1]]
                    orfs[gene]["nts"] = orfnts.upper()
            ref_orfs = orfs
            ref_seq = nts.upper()

    return(ref_id, ref_seq, ref_type, ref_orfs)

def aa_call(codon):
    """
    Called to return an amino acid enoded by the passed codon
    Parameters:
    codon - 3 nt sequence
    Functionality:
    looks up the codon in the dict and returns its value
    Returns amino acid or '?' if the codon isn't a valid codon
    """
    AADict = {
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
    if codon in AADict:
        AA = AADict[codon]

    return(AA)

def singlet_codon_call(nt_pos, nt, ref_nts):
    """
    Called to determine the amino acid and its position based on a single nt change relative to a referernce sequence
    Parameters:
        nt_pos - position in the reference sequence of the changed nt (1 indexed)
        nt - the new nucleotide
        ref_nts - string of the reference nts (0 indexed) encoding a peptide
    Functionality:
        calculates the position of the AA coded from the codon that the new nt is a part of and the new codon and encoded amino acid
    Returns the amino acid position (1 indexed) and the aa_call of the new codon
    """
    aa_pos = (nt_pos-1)//3
    aa_mod = (nt_pos-1)%3
    codon = ""
    try:
        if aa_mod == 0:
            codon = nt+ref_nts[nt_pos]+ref_nts[nt_pos+1]
        elif aa_mod == 1:
            codon = ref_nts[nt_pos-2]+nt+ref_nts[nt_pos]
        elif aa_mod == 2:
            codon = ref_nts[nt_pos-3]+ref_nts[nt_pos-2]+nt
    except:
        codon = "XXX"

    return(aa_pos+1, aa_call(codon))

def sam_line_parser(args, ref, file):
    """
    Called to pasrse through a sam file line by line and collect information for further processing
    Parameters:
    args - argument values
    ref - tuple of the information pulled from the reference file
    file - name of the sam file to be processed
    Functionality:
    opens the indicated file and runs through each line that isn't a header (starts with '@')
    and where the sequence mapped is mapped to the reference
    based on the cigar string in the line, the sequence is used to collect the nt called at each position, the indels,
    the snp variations, total number of reads mapped and the coverage
    information from each line is collected and returned
    Returns a dictionary of dictionaries for each position's nt calls, a dictionary of all the insertions,
    a list of each read id and its variant nt sequence, a dictionary of each unique variant nt sequence with total counts as values,
    the total count of reads mapped to the reference, and a dictionary for coverage of each position
    """

    samp=file.replace(file.split(".")[-1], '')[:-1]
    nt_call_dict_dict = {}
    ins_nt_dict = {}
    reads_list = []
    col_reads = {}
    sam_read_count = 0
    sam_line_count = 0
    coverage = {}

    comped = False
    if file.lower().endswith(".sam"):
        sam_fh = open(file, "r")
    elif file.lower().endswith(".bam"):
        sam_fh = pysam.AlignmentFile(file, "rb")
        comped = True
    elif file.lower().endswith(".cram"):
        sam_fh = pysam.AlignmentFile(file, "rc")
        comped = True
    if comped and not bc:
        print("bam/cram parsing failed.  need to install pysam")
        return
    for line in sam_fh:
        if comped:
            line = str(line)
        if not line.startswith('@'): # ignore header lines
            split_line = line.split("\t")
            if ref[0].upper() == split_line[2].upper() or comped: # check map ID matches referecne ID
                if int(split_line[4]) > 0:  # Check mapping score is positive

                    sam_line_count += 1
                    query_seq = split_line[9].upper()
                    if query_seq.strip('ATCGN-'):
                        print(f"Invalid character(s) {query_seq.strip('ATCGN-')} in sequence for line {sam_line_count}, ID {split_line[0]}")
                        print("skipping")
                        continue

                    CIGAR = split_line[5]
                    cigar_num = CIGAR
                    for c in "MIDSH":
                        cigar_num = cigar_num.replace(c,"")

                    if not cigar_num.isnumeric():
                        print(f"Non-standard CIGAR string {CIGAR} {cigar_num} for line {sam_line_count}, ID {split_line[0]}. Skipping")
                        continue

                    reads_count=1
                    if args.use_count == 1: # get the unique sequence counts
                        if '-' in split_line[0] and '=' in split_line[0]:
                            try:
                                eq_split = split_line[0].split('=')
                                dash_split = split_line[0].split('-')
                                if len(eq_split[-1]) > len(dash_split[-1]):
                                    reads_count=int(dash_split[-1])
                                else:
                                    reads_count=int(eq_split[-1])
                            except:
                                pass
                        elif '-' in split_line[0]:
                            try:
                                reads_count=int(split_line[0].split('-')[-1])
                            except:
                                pass
                        elif '=' in split_line[0]:
                            try:
                                reads_count=int(split_line[0].split('=')[-1])
                            except:
                                pass
                    sam_read_count += reads_count

                    read_start_pos = int(split_line[3])
                    readID = split_line[0]
                    run_length = 0
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
                                    insert_position = q_pars_pos+read_start_pos
                                    iSeq = query_seq[query_pos: query_pos+run_length]
                                    istring = str(insert_position)+'-insert'+iSeq
                                    mutations.append(istring)
                                    if args.nt_call == 1:
                                        # add insertion to dict
                                        try:
                                            ins_nt_dict[insert_position]
                                        except:
                                            ins_nt_dict[insert_position] = {istring : reads_count}
                                        else:
                                            try:
                                                ins_nt_dict[insert_position][istring] += reads_count
                                            except:
                                                ins_nt_dict[insert_position][istring] = reads_count

                                    query_pos = query_pos + run_length

                            elif C == 'D':
                                delPos = q_pars_pos+read_start_pos
                                delnts = ref[1][delPos-1:delPos+run_length-1]
                                delstring = delnts+str(delPos)+'-'+str(delPos+run_length-1)+'del'
                                mutations.append(delstring)
                                if args.nt_call == 1:
                                    for N in range(delPos, delPos+int(run_length)):
                                        try:
                                            nt_call_dict_dict[N]
                                        except:
                                            nt_call_dict_dict[N] = {'A' : 0,
                                                                    'T' : 0,
                                                                    'C' : 0,
                                                                    'G' : 0,
                                                                    'N' : 0,
                                                                    '-' : 0}
                                            nt_call_dict_dict[N]['-'] = reads_count
                                        else:
                                            nt_call_dict_dict[N]['-'] += reads_count
                                q_pars_pos = q_pars_pos + run_length

                            elif C == 'M':
                                offset = q_pars_pos-query_pos
                                refPos = read_start_pos+offset
                                for ntPos in range(query_pos, query_pos+run_length):
                                    if query_seq[ntPos] != ref[1][refPos+ntPos-1]:
                                        mutations.append(ref[1][refPos+ntPos-1]+str(refPos+ntPos)+query_seq[ntPos])
                                    if args.nt_call == 1:
                                        if not query_seq[ntPos] in ['A', 'T' ,'C', 'G', 'N', '-']:
                                            mut_nt = 'N'
                                        else:
                                            mut_nt = query_seq[ntPos]
                                        try:
                                            nt_call_dict_dict[refPos+ntPos]
                                        except:
                                            nt_call_dict_dict[refPos+ntPos] = {'A' : 0,
                                                                               'T' : 0,
                                                                               'C' : 0,
                                                                               'G' : 0,
                                                                               'N' : 0,
                                                                               '-' : 0}
                                            nt_call_dict_dict[refPos+ntPos][mut_nt] = reads_count
                                        else:
                                            nt_call_dict_dict[refPos+ntPos][mut_nt] += reads_count
                                q_pars_pos = q_pars_pos + run_length
                                query_pos = query_pos + run_length

                            run_length = 0

                        else:
                            run_length = (10 * run_length) + int(C)
                    # END CIGAR PARSE

                    seq_end_pos = read_start_pos+q_pars_pos-1
                    if args.wgs == 1 or args.nt_call == 1:
                        for i in range(read_start_pos, seq_end_pos+1): # update coverage
                            try:
                                coverage[i] += reads_count
                            except:
                                coverage[i] = reads_count

                    if len(mutations) == 0: # record reference counts
                        if args.read == 1:
                            reads_list.append([readID, 'Reference', read_start_pos, seq_end_pos])
                        if args.seq == 1 or args.covar == 1 or args.indel == 1:
                            try:
                                col_reads['Reference']
                            except:
                                col_reads['Reference'] = { str(read_start_pos)+'x'+str(seq_end_pos): reads_count}
                            else:
                                try:
                                    col_reads['Reference'][str(read_start_pos)+'x'+str(seq_end_pos)] += reads_count
                                except:
                                    col_reads['Reference'][str(read_start_pos)+'x'+str(seq_end_pos)] = reads_count
                    else:
                        mutations = " ".join(mutations)
                        if args.read == 1:
                            reads_list.append([readID, mutations, read_start_pos, seq_end_pos])
                        if args.seq == 1 or args.covar == 1 or args.indel == 1:
                            try:
                                col_reads[mutations]
                            except:
                                col_reads[mutations] = { str(read_start_pos)+'x'+str(seq_end_pos): reads_count}
                            else:
                                try:
                                    col_reads[mutations][str(read_start_pos)+'x'+str(seq_end_pos)] += reads_count
                                except:
                                    col_reads[mutations][str(read_start_pos)+'x'+str(seq_end_pos)] = reads_count

    sam_fh.close()
    # END SAM LINES
    print(f"End SAM line parsing for {samp}")

    return(nt_call_dict_dict, ins_nt_dict, reads_list, col_reads, sam_read_count, coverage)

def fasta_snp_call(mut, ref):
    """
    Called to determine the changed amino encoding from a single mutation event 1 indel or 1 snp based on a fasta reference
    Parameters:
    mut - string of the snp event
    ref - tuple of the information pulled from the reference file
    Functionality: for indels, determines if the change effects multiple reference codons and how. for snp, singlet_codon_call function is called to get changes
    Returns original mutation string appended with amino acid change string
    """

    if 'ins' in mut:
        istring = ''
        iSeq = mut.split('insert')[1]
        insert_position = int(mut.split('-')[0])
        run_length = len(iSeq)
        if (run_length % 3 == 0):
            iProt = ''
            if insert_position % 3 == 1:
                for x in range(0, (run_length//3)):
                    AA = aa_call(iSeq[x*3]+iSeq[x*3+1]+iSeq[x*3+2])
                    iProt += AA
                istring += mut+'(' + str(((insert_position-1)//3)+1) + iProt + ')'
            elif insert_position % 3 == 2:
                ipSeq = ref[1][insert_position-2]+iSeq+ref[1][insert_position-1:insert_position+1]
                for x in range(0, (run_length//3)+1):
                    AA = aa_call(ipSeq[x*3]+ipSeq[x*3+1]+ipSeq[x*3+2])
                    iProt += AA
                istring += ref[1][insert_position-2:insert_position+1]+str(insert_position-1)+'-'+str(insert_position+1)+ipSeq+'insert(' + ref[3][1][(insert_position-1)//3] + str(((insert_position-1)//3)+1) + iProt + ')'
            else:
                ipSeq = ref[1][insert_position-3:insert_position-1]+iSeq+ref[1][insert_position-1]
                for x in range(0, (run_length//3)+1):
                    AA = aa_call(ipSeq[x*3]+ipSeq[x*3+1]+ipSeq[x*3+2])
                    iProt += AA
                istring += ref[1][insert_position-3:insert_position]+str(insert_position-2)+'-'+str(insert_position)+ipSeq+'insert(' + ref[3][1][(insert_position-1)//3] + str(((insert_position-1)//3)+1) + iProt + ')'
        else:
            istring += mut+'('+ str(((insert_position-1)//3)+1) +'fs)'
        return(istring)

    elif 'del' in mut:
        delPos = ''
        delnts = ''
        run_length = 0
        for c in mut.split('-')[0]:
            if c.isdigit():
                delPos += c
            else:
                run_length += 1
        delPos = int(delPos)
        delstring = ''
        delendPos = delPos+run_length-1
        newAArefpos = (delPos-1) // 3
        if (run_length % 3 == 0):
            if (delPos) % 3 == 1:
                if run_length // 3 == 1:
                    delstring += mut + '(' + ref[3][1][newAArefpos:newAArefpos+(run_length//3)] + str(newAArefpos+1) + 'del)'
                else:
                    delstring += mut + '(' + ref[3][1][newAArefpos:newAArefpos+(run_length//3)] + str(newAArefpos+1) + '-' + str(newAArefpos+run_length//3) + 'del' + ')'
            else:
                if (delPos) % 3 == 2:
                    oldnts = ref[1][delPos-2:delendPos+2]
                    newcodon = ref[1][delPos-2]+ref[1][delendPos:delendPos+2]
                    delstring += oldnts+str(delPos-1)+'-'+str(delendPos+2)+newcodon
                elif (delPos) % 3 == 0:
                    oldnts = ref[1][delPos-3:delendPos+1]
                    newcodon = ref[1][delPos-3:delPos-1]+ref[1][delendPos]
                    delstring += oldnts+str(delPos-2)+'-'+str(delendPos+1)+newcodon
                delstring += 'del(' + ref[3][1][newAArefpos:newAArefpos+(run_length//3)+1] + str(newAArefpos+1) + '-' + str(newAArefpos+1+run_length//3) + aa_call(newcodon) + 'del)'
        else:
            delstring += mut + '(' + str(newAArefpos+1) + 'fs)'
        return(delstring)
    else:
        AAinfo = singlet_codon_call(int(mut[1:-1]), mut[-1], ref[1])
        AAstring = '('+ref[3][1][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+')'
        return(mut+AAstring)

def gb_snp_call(mut, ref):
    """
    Called to determine the changed amino encoding from a single mutation event 1 indel or 1 snp based on a gb reference
    Parameters:
    mut - string of the snp event
    ref - tuple of the information pulled from the reference file
    Functionality: Determines amino acid changes for each orf effected.  For indels, determines if the change effects multiple reference codons and how.
    for snp, singlet_codon_call function is called to get changes
    Returns original mutation string appended with amino acid change string
    """

    if 'ins' in mut:
        iSeq = mut.split('insert')[1]
        insert_position = int(mut.split('-')[0])
        run_length = len(iSeq)
        iProt = ''
        for orf in ref[3]:
            orflength = 0
            for rf in ref[3][orf]['reading frames']:
                if insert_position >= rf[0] and insert_position <= rf[1]:
                    iProt += "|(" + orf + ":"
                    orfPos = 1 + insert_position - rf[0] + orflength
                    AA = ''
                    if (run_length % 3 == 0):
                        if orfPos % 3 == 1:
                            for x in range(0, (run_length//3)):
                                AA += aa_call(iSeq[x*3]+iSeq[x*3+1]+iSeq[x*3+2])
                            iProt += str(orfPos)+'insert'+iSeq+'(' + str(((orfPos-1)//3)+1) + AA
                        else:
                            if orfPos % 3 == 2:
                                ipSeq = ref[1][insert_position-2]+iSeq+ref[1][insert_position-1:insert_position+1]
                                for x in range(0, (run_length//3)+1):
                                    AA += aa_call(ipSeq[x*3]+ipSeq[x*3+1]+ipSeq[x*3+2])
                                iProt += ref[1][insert_position-2:insert_position+1]+str(orfPos-1)+'-'+str(orfPos+1)+ipSeq
                            else:
                                ipSeq = ref[1][insert_position-3:insert_position-1]+iSeq+ref[1][insert_position-1]
                                for x in range(0, (run_length//3)+1):
                                    AA += aa_call(ipSeq[x*3]+ipSeq[x*3+1]+ipSeq[x*3+2])
                                iProt += ref[1][insert_position-3:insert_position]+str(orfPos-2)+'-'+str(orfPos)+ipSeq
                            iProt += "insert(" + ref[3][orf]["AAs"][(orfPos-1)//3] + str(((orfPos-1)//3)+1) + AA
                    else:
                        iProt += str(orfPos)+'insert'+iSeq+'(' + str(((orfPos-1)//3)+1) + 'fs'
                    iProt += "))"
                orflength += rf[1] - rf[0] + 1
        return(mut+iProt)

    elif 'del' in mut:
        delPos = ''
        delnts = ''
        run_length = 0
        for c in mut.split('-')[0]:
            if c.isdigit():
                delPos += c
            else:
                run_length += 1
        delPos = int(delPos)
        delProt = ''
        delendPos = delPos+run_length-1
        for orf in ref[3]:
            orflength = 0
            for rf in ref[3][orf]['reading frames']:
                if delPos > rf[1] or delendPos < rf[0]:
                    pass
                elif delPos >= rf[0] and delendPos <= rf[1]:
                    delProt += "|(" + orf + ":"
                    orfPos = 1 + delPos - rf[0] + orflength
                    orfendPos = orfPos+run_length-1
                    startcodon = ((orfPos-1)//3)+1
                    if (run_length % 3 == 0):
                        endcodon = ((orfendPos-1)//3)+1
                        if ((orfPos) % 3 == 1 ):
                            delProt += delnts + str(orfPos) +'-'+ str(orfendPos) + 'del(' + ref[3][orf]["AAs"][startcodon-1:endcodon] + str(startcodon)
                            if run_length // 3 == 1:
                                delProt += 'del'
                            else:
                                delProt += '-' + str(endcodon) + 'del'
                        else:
                            delAA = ''
                            if (orfPos) % 3 == 2:
                                oldnts = ref[1][delPos-2:delendPos+2]
                                newcodon = ref[1][delPos-2]+ref[1][delendPos:delendPos+2]
                                delProt += oldnts + str(orfPos-1) +'-'+ str(orfendPos+2) + newcodon
                            else:
                                oldnts = ref[1][delPos-3:delendPos+1]
                                newcodon = ref[1][delPos-3:delPos-1]+ref[1][delendPos]
                                delProt += oldnts + str(orfPos-2) +'-'+ str(orfendPos+1) + newcodon
                            delProt += 'del(' + ref[3][orf]["AAs"][startcodon-1:endcodon] + str(startcodon) \
                                    + '-' + str(endcodon) + aa_call(newcodon)+ 'del'
                    else:
                        delProt += delnts + str(orfPos) +'-'+ str(orfendPos) + 'del(' + str(startcodon) + 'fs'
                    delProt += "))"
                elif delendPos >= rf[0] and delPos <= rf[0]:
                    delProt += "|(" + orf
                    if orflength == 0:
                        delProt += ":Start_disrupted)"
                    else:
                        delProt += ":fs/Splicing_disrupted)"
                else:
                    delProt += "|(" + orf + ":fs/Splicing/Termination_disrupted)"
                orflength += rf[1] - rf[0] + 1
        return(mut+delProt)
    else:
        PM = ''
        PMPos = int(mut[1:-1])
        for orf in ref[3]:
            orflength = 0
            for rf in ref[3][orf]['reading frames']:
                if PMPos >= rf[0] and PMPos <= rf[1]:
                    ORFPos = PMPos-rf[0]+1+orflength
                    AAinfo = singlet_codon_call(ORFPos, mut[-1], (ref[3][orf]['nts']))
                    PM += '|(' + orf + ":" + mut[0] + str(ORFPos) + mut[-1] + "(" + ref[3][orf]["AAs"][AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+'))'
                orflength += rf[1] - rf[0] + 1
        return(mut+PM)

def get_combos(qlist, clen):
    """
    Called to get combinations of cosegregating polymorphisms
    Parameters:
    qlist - list of polymorphisms
    clen - number of polymorphisms to report together

    Functionality:
        uses itertools.combinations to get sets of combinations for each length
    Returns list of combinations
    """
    combos = []
    if (clen == 0 or clen > len(qlist)):
        clen = len(qlist)
    for N in range(1, clen+1):
        for comb in itertools.combinations(qlist, N):
            combos.append(' '.join(comb))
    return(combos)

def print_reads(samp, reads_list, col_reads, args):
    """
    Called to print the reads output
    Parameters:
    samp - name of the sam being processed
    reads_list - list of read lines
    col_reads - dict of collapsed sequences with amino acid seq if applicable
    args - arguement values

    Functionality:
        opens output file and writes lines
    Returns nothing
    """
    reads_fh = open(samp+'_reads.tsv', "w")
    for line in reads_list:
        if args.AAreport == 1:
            try:
                reads_fh.write(f"{line[0]}\t{line[2]} {col_reads[line[1]]['AA_sequence']} {line[3]}\n")
            except:
                reads_fh.write(f"{line[0]}\t{line[2]} {line[1]} {line[3]}\n")
        else:
            reads_fh.write(f"{line[0]}\t{line[2]} {line[1]} {line[3]}\n")
    reads_fh.close()

def print_unique_seq(samp, sam_read_count, col_reads, coverage, args):
    """
    Called to print the unique sequences output
    Parameters:
    samp - name of the sam being processed
    sam_read_count - number of total reads mapped in the sam file
    col_reads - dict of collapsed sequences with amino acid seq if applicable
    coverage - dict of coverage at reference positions
    args - arguement values

    Functionality:
        opens output file(s), collects and sorts unique sequences for wgs mode or not, and writes lines
    Returns nothing
    """
    seq_fh = open(samp+'_unique_seqs.tsv', "w")
    seq_fh.write(samp+"("+str(sam_read_count)+")\n")
    seq_fh.write("Unique Sequence\tCount\tAbundance\n")
    if args.AAcentered == 1 and args.AAreport == 1:
        AAseq_fh = open(samp+'_AA_unique_seqs.tsv', "w")
        AAseq_fh.write(samp+"("+str(sam_read_count)+")\n")
        AAseq_fh.write("Unique Sequence\tCount\tAbundance\n")

    seq_species = {}

    for SNP_sequence in col_reads:
        for start_end in col_reads[SNP_sequence]:
            if 'x' in start_end:
                if args.wgs == 1:
                    seq_species[start_end+'x'+SNP_sequence] = col_reads[SNP_sequence][start_end]
                else:
                    try:
                        seq_species[SNP_sequence] += col_reads[SNP_sequence][start_end]
                    except:
                        seq_species[SNP_sequence] = col_reads[SNP_sequence][start_end]


    sorted_seq = sorted(seq_species, key=seq_species.__getitem__, reverse=True)
    for key in sorted_seq:
        if seq_species[key] >= args.min_count:
            splitseqs = key.split('x')
            if args.wgs == 0:
                abund = seq_species[key] / sam_read_count
            elif args.wgs == 1:
                cov = []
                for x in range(int(splitseqs[0]), int(splitseqs[1])+1):
                    cov.append(coverage[x])
                abund = seq_species[key] / min(cov)
            if (abund >= args.min_samp_abund):
                try:
                    splitseqs[2]
                except:
                    if args.AAreport == 1:
                        try:
                            seq_fh.write(f"{col_reads[key]['AA_sequence']}\t{seq_species[key]}\t{abund:.3f}\n")
                        except:
                            seq_fh.write(f"{key}\t{seq_species[key]}\t{abund:.3f}\n")
                    else:
                        seq_fh.write(f"{key}\t{seq_species[key]}\t{abund:.3f}\n")
                    if args.AAcentered == 1 and args.AAreport == 1:
                        try:
                            seq_fh.write(f"{col_reads[key]['AA_centered']}\t{seq_species[key]}\t{abund:.3f}\n")
                        except:
                            if not key == 'Reference':
                                print('Non Reference seq without AA centereed output : '+splitseqs[2])
                else:
                    if args.AAreport == 1:
                        try:
                            seq_fh.write(f"{splitseqs[0]} {col_reads[splitseqs[2]]['AA_sequence']} {splitseqs[1]}\t{seq_species[key]}\t{abund:.3f}\n")
                        except:
                            seq_fh.write(f"{splitseqs[0]} {splitseqs[2]} {splitseqs[1]}\t{seq_species[key]}\t{abund:.3f}\n")
                    else:
                        seq_fh.write(f"{splitseqs[0]} {splitseqs[2]} {splitseqs[1]}\t{seq_species[key]}\t{abund:.3f}\n")
                    if args.AAcentered == 1 and args.AAreport == 1:
                        try:
                            AAseq_fh.write(f"{splitseqs[0]} {col_reads[splitseqs[2]]['AA_centered']} {splitseqs[1]}\t{seq_species[key]}\t{abund:.3f}\n")
                        except:
                            if not splitseqs[2] == 'Reference':
                                print('Non Reference seq without AA centereed output : '+splitseqs[2])
        else:
            break

    seq_fh.close()
    if args.AAcentered == 1 and args.AAreport == 1:
        AAseq_fh.close()

def print_indels(samp, sam_read_count, indel_dict, coverage, args):
    """
    Called to print the indel output
    Parameters:
    samp - name of the sam being processed
    sam_read_count - number of total reads mapped in the sam file
    indel_reads - dict of indels
    coverage - dict of coverage at reference positions
    args - arguement values

    Functionality:
        sorts indels and gets abundance for wgs mode or not, if any indels pass abundance check, opens file and writes lines
    Returns nothing
    """
    sorted_indels = sorted(indel_dict, key=indel_dict.__getitem__, reverse=True)
    indels_to_write = []
    for key in sorted_indels:
        if indel_dict[key] >= args.min_count:
            if args.wgs == 0:
                abund = indel_dict[key] / sam_read_count
            elif args.wgs == 1:
                indelPos = ''
                for c in key.strip('ATCGN'):
                    if c.isdigit():
                        indelPos += c
                    else:
                        break
                abund = indel_dict[key] / coverage[int(indelPos)]
            if abund >= args.min_samp_abund:
                indels_to_write.append(f"{key}\t{indel_dict[key]}\t{abund:.3f}\n")
        else:
            break
    if len(indels_to_write) > 0:
        indel_fh = open(samp+'_indels.tsv', "w")
        indel_fh.write(samp+"("+str(sam_read_count)+")\n")
        indel_fh.write("Indel\tCount\tAbundance\n")
        for indel_entry in indels_to_write:
            indel_fh.write(indel_entry)
        indel_fh.close()

def print_covars(samp, sam_read_count, col_reads, coverage, args, aa_centered):
    """
    Called to print the covars output
    Parameters:
    samp - name of the sam being processed
    sam_read_count - number of total reads mapped in the sam file
    col_reads - dict of collapsed sequences with amino acid seq if applicable
    coverage - dict of coverage at reference positions
    args - arguement values
    aa_centered- dict of coverage at posistion of polymoprhisms in the aa centered form

    Functionality:
        parses each unique variant sequence, normal and aa centered, if applicable, to get the cosegregating variations and tiled coverage if enabled,
        if there are combos that meet the count threshold output files are opened and lines written
    Returns nothing
    """
    ntcombinations = {}
    AAcombinations = {}
    tiles = {}
    for sequence in col_reads:
        ntsingles = []
        AAsingles = []

        if args.wgs == 1 and sequence == "Reference":
            continue
        try:
            ntsingles = col_reads[sequence]['AA_sequence'].split(' ')
        except:
            ntsingles = sequence.split(' ')
        if aa_centered:
            AAsingles = col_reads[sequence]['AA_centered'].split(' ')
        sequence_count = 0
        for value in col_reads[sequence].values():
            if isinstance(value, int):
                sequence_count += value

        if args.wgs == 1 and args.covar_tile_coverage == 1:
            for start_end in col_reads[SNP_sequence]:
                if 'x' in start_end:
                    start_pos, end_pos = start_end.split('x')
                    start_pos = int(start_pos)
                    end_pos = int(end_pos)
                    try:
                        tiles[start_pos]
                    except:
                        tiles[start_pos] = {end_pos : col_reads[SNP_sequence][start_end]}
                    else:
                        try:
                            tiles[start_pos][end_pos] += col_reads[SNP_sequence][start_end]
                        except:
                            tiles[start_pos][end_pos] = col_reads[SNP_sequence][start_end]
        if ntsingles and len(ntsingles) <= args.max_dist:
            for combo in get_combos(ntsingles, args.max_covar):
                if not combo in ntcombinations:
                    ntcombinations[combo] = sequence_count
                else:
                    ntcombinations[combo] += sequence_count

        if AAsingles and len(AAsingles) <= args.max_dist:
            for combo in get_combos(AAsingles, args.max_covar):
                if not combo in AAcombinations:
                    AAcombinations[combo] = sequence_count
                else:
                    AAcombinations[combo] += sequence_count
    for combinations, samp_tag in ((ntcombinations, ''), (AAcombinations, '_AA')):
        try:
            max_count = max(combinations.values())
        except:
            pass
        else:
            if max_count >= args.min_count:
                covar_fh = open(samp+samp_tag+'_covars.tsv', "w")
                covar_fh.write(samp+"("+str(sam_read_count)+")\n")
                covar_fh.write("Co-Variants\tCount\tAbundance\n")
                sortedcombos = sorted(combinations, key=combinations.__getitem__, reverse=True)
                for key in sortedcombos:
                    if combinations[key] >= args.min_count:
                        if args.wgs == 0:
                            abund = combinations[key] / sam_read_count
                        elif args.wgs == 1:
                            splitcombos = key.split()
                            if not samp_tag:
                                startcovPos = ''
                                for c in splitcombos[0].strip('ATGC'):
                                    if c.isdigit():
                                        startcovPos += c
                                    else:
                                        break
                                endcovPos = ''
                                split_mut = splitcombos[-1].split('(')[0]
                                if 'del' in split_mut:
                                    pos_string = split_mut.strip('ATGCN').split('-')[1]
                                else:
                                    pos_string = split_mut.strip('ATGCN')
                                for c in pos_string:
                                    if c.isdigit():
                                        endcovPos += c
                                    else:
                                        break
                                if not endcovPos:
                                    print(key)
                                    print(split_mut)
                                    print(pos_string)
                            else:
                                startcovPos = aa_centered[splitcombos[0]][0]
                                endcovPos = aa_centered[splitcombos[-1]][-1]
                            if startcovPos == endcovPos:
                                abund = combinations[key] / coverage[int(startcovPos)]

                            elif args.covar_tile_coverage == 0:
                                coveragevals = []
                                for i in range(int(startcovPos), int(endcovPos)+1):
                                    coveragevals.append(coverage[i])
                                    abund = combinations[key] / min(coveragevals)
                            elif args.covar_tile_coverage == 1:
                                coverageval = 0
                                for tile_start in tiles:
                                    if int(startcovPos) >= int(tile_start):
                                        for tile_end in tiles[tile_start]:
                                            if int(endcovPos) <= int(tile_end):
                                                coverageval += tiles[tile_start][tile_end]
                                abund = combinations[key] / coverageval
                        if abund >= args.min_samp_abund:
                            covar_fh.write(f"{key}\t{combinations[key]}\t{abund:.3f}\n")
                    else:
                        break

                covar_fh.close()

def get_nt_indels(col_reads):
    """
    Called to collect indels if amino acid reporting is disabled
    Parameters:
    col_reads - dict of unique variant sequences
    Functionality: goes through unique sequences looking for indels and adds them to a new dict
    Returns the new indel dict
    """
    indel_dict = {}
    for SNP_sequence in col_reads:
        if 'del' in SNP_sequence or 'insert' in SNP_sequence:
            for mut in SNP_sequence.split(' '):
                if 'del' in mut or 'ins' in mut:
                    try:
                        indel_dict[mut] += sum(col_reads[SNP_sequence].values())
                    except:
                        indel_dict[mut] = sum(col_reads[SNP_sequence].values())
    return(indel_dict)

def fa_sam_parse(args, ref, file):
    """
    Called to handle main sam info parsing logic when using a fasta reference
    Parameters:
    args - argument values
    ref - tuple of referecne information
    file - name of sam file
    Functionality:
    calls sam_line_parser to get info from sam, then processes the information to get amino acid change information if applicable.
    prints outputs, mainly by function calls.  nt call output is the exception
    Returns nothing
    """

    samp=file.replace(file.split(".")[-1], '')[:-1]
    print(f"Starting {samp} processing")
    nt_call_dict_dict, ins_nt_dict, reads_list, col_reads, sam_read_count, coverage = sam_line_parser(args, ref, file)

    if sam_read_count == 0:
        print(f"No Reads for {samp}")
    else:

        indel_dict = {}
        nt_to_AA_dict = {}

        if args.AAreport == 1 and (args.seq == 1 or args.read == 1 or args.covar == 1 or args.indel == 1):
            for SNP_sequence in col_reads:
                if 'Reference' == SNP_sequence:
                    pass
                elif args.AAcodonasMNP == 1:
                    try:
                        col_reads[SNP_sequence]['AA_sequence']
                    except:
                        MNPs = []
                        curMNP = ''
                        last_codon = -1
                        fshift = 0
                        for mut in SNP_sequence.split(' '): # collect together mutations that affect the same codon(s)
                            startPos = int(mut.split('-')[0].strip('ATCGN'))
                            endPos = startPos
                            if 'del' in mut:
                                endPos = int(mut.split('-')[1].strip('del'))

                            startcodon = (((startPos)-1)//3)+1
                            endcodon = ((endPos-1)//3)+1

                            mutshift = 0
                            if 'insert' in mut:
                                mutshift += len(mut.split('(')[0].split('insert')[1])
                            elif 'del' in mut:
                                mutshift -= len(mut.split('-')[0].strip('0123456789'))

                            if curMNP:
                                if not mutshift % 3 == 0:
                                    if (fshift % 3) == 0:
                                        if startcodon == last_codon:
                                            curMNP.append([mut, startPos])
                                            fshift += mutshift
                                        else:
                                            MNPs.append(curMNP)
                                            curMNP = [[mut, startPos]]
                                            fshift = mutshift

                                    else:
                                        if startcodon < last_codon + 5:
                                            curMNP.append([mut, startPos])
                                            fshift += mutshift
                                        else:

                                            ## todo splitfsMNP() make process
                                            fsfreeMNP = []
                                            for SNP in curMNP:
                                                if 'fs' in SNP:
                                                    if fsfreeMNP:
                                                        MNPs.append(fsfreeMNP)
                                                        fsfreeMNP = []
                                                    MNPs.append([SNP])
                                                else:
                                                    fsfreeMNP.append(SNP)
                                            if fsfreeMNP:
                                                MNPs.append(fsfreeMNP)
                                            curMNP = [[mut, startPos]]
                                            fshift = mutshift
                                else:
                                    if startcodon == last_codon:
                                        curMNP.append([mut, startPos])
                                    else:
                                        if (fshift % 3) == 0:
                                            MNPs.append(curMNP)
                                        else:
                                            ## splitfsMNP() make process
                                            fsfreeMNP = []
                                            for SNP in curMNP:
                                                if 'fs' in SNP:
                                                    if fsfreeMNP:
                                                        MNPs.append(fsfreeMNP)
                                                        fsfreeMNP = []
                                                    MNPs.append([SNP])
                                                else:
                                                    fsfreeMNP.append(SNP)
                                            if fsfreeMNP:
                                                MNPs.append(fsfreeMNP)
                                        curMNP = [[mut, startPos]]
                                        fshift = 0
                            else:
                                curMNP = [[mut, startPos]]
                                if 'insert' in mut:
                                    fshift = len(mut.split('(')[0].split('insert')[1])
                                elif 'del' in mut:
                                    fshift = len(mut.split('-')[0].strip('0123456789'))
                            last_codon = endcodon
                        if curMNP:
                            MNPs.append(curMNP)

                        if MNPs:
                            combinedmut = []
                            for entry in MNPs:
                                if len(entry) > 1: # re-calcs sequence for codons affected by multiple changes
                                    nt_to_AA_key = []
                                    for snpmut in entry:
                                        nt_to_AA_key.append(snpmut[0])
                                    nt_to_AA_key = ' '.join(nt_to_AA_key)
                                    try:
                                        newmut = nt_to_AA_dict[nt_to_AA_key]
                                    except:
                                        MNPreport = ''
                                        mutntseq = ''
                                        orfstartpos = entry[0][1]
                                        orfendpos = entry[-1][1]

                                        if 'del' in entry[-1][0]:
                                            orfendpos += len(entry[-1][0].split('-')[0].strip('0123456789')) - 1

                                        cur_nts = {}
                                        for i in range(orfstartpos, orfendpos+1):
                                            cur_nts[i] = ref[1][i-1]
                                        for curPM in entry:
                                            if 'insert' in curPM[0]:
                                                cur_nts[curPM[1]] = curPM[0].split('insert')[-1] + cur_nts[curPM[1]]
                                            elif 'del' in curPM[0]:
                                                for i in range(curPM[1], curPM[1]+len(curPM[0].split('-')[0].strip('0123456789'))):
                                                    cur_nts[i] = ''
                                            else:
                                                cur_nts[curPM[1]] = curPM[0][-1]
                                        for pos in cur_nts:
                                            mutntseq += cur_nts[pos]

                                        startcodonpos = ((orfstartpos-1)//3)+1
                                        endcodonpos = ((orfendpos-1)//3)+1
                                        startmod = orfstartpos % 3
                                        endmod = orfendpos % 3
                                        wtorfstartpos = orfstartpos
                                        wtorfendpos = orfendpos
                                        if startmod == 2:
                                            wtorfstartpos = orfstartpos-1
                                            mutntseq = ref[1][orfstartpos-2] + mutntseq
                                        elif startmod == 0:
                                            wtorfstartpos = orfstartpos-2
                                            mutntseq = ref[1][orfstartpos-3:orfstartpos-1] + mutntseq
                                        if endmod == 2:
                                            wtorfendpos = orfendpos+1
                                            mutntseq += ref[1][orfendpos]
                                        elif endmod == 1:
                                            wtorfendpos = orfendpos+2
                                            mutntseq += ref[1][orfendpos:orfendpos+2]
                                        wtntseq = ref[1][wtorfstartpos-1:wtorfendpos]
                                        wtAAseq = ref[3][1][startcodonpos-1:endcodonpos]
                                        mutAAseq = ''
                                        for i in range(0, (len(mutntseq)//3)):
                                            mutAAseq += aa_call(mutntseq[i*3:(i*3)+3])
                                        if len(mutntseq) < len(wtntseq):
                                            mutntseq += 'del'
                                        elif len(mutntseq) > len(wtntseq):
                                            mutntseq += 'insert'
                                        if len(mutntseq) % 3 != 0:
                                            mutAAseq += 'fs'
                                        if startcodonpos == endcodonpos:
                                            newmutstring = f"{wtntseq}{wtorfstartpos}-{wtorfendpos}{mutntseq}({wtAAseq}{startcodonpos}{mutAAseq})"
                                        else:
                                            newmutstring = f"{wtntseq}{wtorfstartpos}-{wtorfendpos}{mutntseq}({wtAAseq}{startcodonpos}-{endcodonpos}{mutAAseq})"

                                        nt_to_AA_dict[nt_to_AA_key] = newmutstring
                                        newmut = newmutstring
                                else:
                                    if entry[0][0] in nt_to_AA_dict:
                                        newmut = nt_to_AA_dict[entry[0][0]]
                                    else:
                                        mutstring = fasta_snp_call(entry[0][0], ref)
                                        nt_to_AA_dict[entry[0][0]] = mutstring
                                        newmut = mutstring
                                combinedmut.append(newmut)
                                if args.indel == 1:
                                    if 'del' in newmut or 'insert' in newmut:
                                        try:
                                            indel_dict[newmut] += sum(col_reads[SNP_sequence].values())
                                        except:
                                            indel_dict[newmut] = sum(col_reads[SNP_sequence].values())
                        col_reads[SNP_sequence]['AA_sequence'] = " ".join(combinedmut)
                else:
                    try:
                        col_reads[SNP_sequence]['AA_sequence']
                    except:
                        mutations =[]
                        for mut in SNP_sequence.split(' '):
                            try:
                                newmut = nt_to_AA_dict[mut]
                            except:
                                mutstring = fasta_snp_call(mut, ref)
                                nt_to_AA_dict[mut] = mutstring
                                newmut = mutstring
                            mutations.append(newmut)
                            if args.indel == 1:
                                if 'del' in mut or 'insert' in mut:
                                    try:
                                        indel_dict[newmut] += sum(col_reads[SNP_sequence].values())
                                    except:
                                        indel_dict[newmut] = sum(col_reads[SNP_sequence].values())
                        col_reads[SNP_sequence]['AA_sequence'] = " ".join(mutations)

        elif args.indel == 1:
            indel_dict = get_nt_indels(col_reads)

        if args.read == 1:
            print_reads(samp, reads_list, col_reads, args)

        if args.seq == 1:
            print_unique_seq(samp, sam_read_count, col_reads, coverage, args)

        if args.indel == 1 and len(indel_dict) > 0:
            print_indels(samp, sam_read_count, indel_dict, coverage, args)

        if args.covar == 1:
            print_covars(samp, sam_read_count, col_reads, coverage, args, {})

        if args.nt_call == 1:
            ntcall_lines = {'line' : {},
                            'variant' : {}
                            }
            ntcall_fh = open(samp+'_nt_calls.tsv', "w")
            ntcall_fh.write(samp+"("+str(sam_read_count)+")\n")
            if args.ntvar == 1:
                ntcallv_fh = open(samp+'_nt_calls_varonly.tsv', "w")
                ntcallv_fh.write(samp+"("+str(sam_read_count)+")\n")
            sorted_Pos = sorted(nt_call_dict_dict)
            if args.AAreport == 1:

                ntcall_fh.write("Position\tref NT\tAA Pos\tref AA\tA\tT\tC\tG\t-\tN\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                if args.ntvar == 1:
                    ntcallv_fh.write("Position\tref NT\tAA Pos\tref AA\tA\tT\tC\tG\t-\tN\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                consensus = {}
                for Pos in sorted_Pos:

                    try:
                        total = coverage[Pos]
                    except:
                        total = 0
                        print(f"coverage of position {Pos} not found")

                    if (total >= (sam_read_count * args.ntabund) or args.wgs == 1) and total >= args.ntcover:

                        try:
                            ins_nt_dict[Pos]
                        except:
                            pass
                        else:
                            i = 1
                            for insertion in ins_nt_dict[Pos]:
                                insert_position = Pos+(i/1000)
                                i_nts = insertion.split('insert')[1]
                                try:
                                    i_AAs = nt_to_AA_dict[insertion].split('(')[1].strip(')')
                                except:
                                    i_AAs = fasta_snp_call(insertion, ref).split('(')[1].strip(')')
                                AA_pos = ((Pos-1)//3)+1
                                iabund = ins_nt_dict[Pos][insertion]/total
                                ntcall_lines['line'][insert_position] = f"{Pos}\t-\t{AA_pos}\t-\t\t\t\t\t\t\t{total}\t{i_nts}\t{ins_nt_dict[Pos][insertion]}"
                                ntcall_lines['line'][insert_position] += f"\t{iabund:.3f}"
                                if 'fs' in i_AAs:
                                    ntcall_lines['line'][insert_position] += f"\t\tfs"
                                else:
                                    split_AAs = i_AAs.split(str(AA_pos))
                                    if split_AAs[0]:
                                        ntcall_lines['line'][insert_position] += f"\t\t{split_AAs[0]}->{split_AAs[1][:-1]}"
                                    else:
                                        ntcall_lines['line'][insert_position] += f"\t\t{split_AAs[1]}"
                                if ( ins_nt_dict[Pos][insertion] >= args.min_count) and (iabund >= args.min_samp_abund):
                                    ntcall_lines['variant'][insert_position] = 1
                                i += 1

                        Pos_calls = {}
                        for key in nt_call_dict_dict[Pos]:
                            Pos_calls[key] = nt_call_dict_dict[Pos][key]
                        sorted_calls = sorted(Pos_calls, key=Pos_calls.__getitem__, reverse=True)

                        ntcall_lines['line'][Pos] = (str(Pos)+"\t"+ref[1][Pos-1]+"\t"+str(((Pos-1)//3)+1)+"\t"+ref[3][1][((Pos-1)//3)])
                        ntcall_lines['line'][Pos] +=("\t"+str(nt_call_dict_dict[Pos]['A'])+"\t"+str(nt_call_dict_dict[Pos]['T'])+"\t"+str(nt_call_dict_dict[Pos]['C']))
                        ntcall_lines['line'][Pos] +=("\t"+str(nt_call_dict_dict[Pos]['G'])+"\t"+str(nt_call_dict_dict[Pos]['-'])+"\t"+str(nt_call_dict_dict[Pos]['N']))
                        ntcall_lines['line'][Pos] +=("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[Pos][sorted_calls[0]]))
                        ntcall_lines['line'][Pos] +=(f"\t{(nt_call_dict_dict[Pos][sorted_calls[0]]/total):.3f}")

                        consensus[Pos] = sorted_calls


                for Pos in sorted_Pos:
                    try:
                        ntcall_lines['line'][Pos]
                    except:
                        pass
                    else:
                        if consensus[Pos][0] != ref[1][Pos-1]:
                            ntcall_lines['variant'][Pos] = 1
                            mod = (Pos)%3

                            if mod == 0:
                                try:
                                    codon = consensus[Pos-2][0]+consensus[Pos-1][0]+ consensus[Pos][0]
                                except:
                                    codon = 'NNN'
                            elif mod == 2:
                                try:
                                    codon = consensus[Pos-1][0]+consensus[Pos][0]+ consensus[Pos+1][0]
                                except:
                                    codon = 'NNN'
                            elif mod == 1:
                                try:
                                    codon = consensus[Pos][0]+consensus[Pos+1][0]+ consensus[Pos+2][0]
                                except:
                                    codon = 'NNN'
                            ntcall_lines['line'][Pos] +=("\t"+aa_call(codon)+"\t"+singlet_codon_call(Pos, consensus[Pos][0], ref[1])[1])
                        else:
                            ntcall_lines['line'][Pos] +=("\t\t")

                        if (nt_call_dict_dict[Pos][consensus[Pos][1]] >= args.min_count) and ((nt_call_dict_dict[Pos][consensus[Pos][1]] / total) >= args.min_samp_abund):
                            ntcall_lines['variant'][Pos] = 1
                            ntcall_lines['line'][Pos] +=(f"\t{consensus[Pos][1]}\t{nt_call_dict_dict[Pos][consensus[Pos][1]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][1]]/total):.3f}"+"\t"+singlet_codon_call(Pos, consensus[Pos][1], ref[1])[1])

                            if (nt_call_dict_dict[Pos][consensus[Pos][2]] >= args.min_count) and (nt_call_dict_dict[Pos][consensus[Pos][2]] / total >= args.min_samp_abund):
                                ntcall_lines['variant'][Pos] = 1
                                ntcall_lines['line'][Pos] +=(f"\t{consensus[Pos][2]}\t{nt_call_dict_dict[Pos][consensus[Pos][2]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][2]]/total):.3f}\t{singlet_codon_call(Pos, consensus[Pos][2], ref[1])[1]}")

                for Pos in ntcall_lines['line']:
                    ntcall_fh.write(ntcall_lines['line'][Pos])
                    ntcall_fh.write("\n")
                    if args.ntvar == 1:
                        try:
                            ntcall_lines['variant'][Pos]
                        except:
                            pass
                        else:
                            ntcallv_fh.write(ntcall_lines['line'][Pos])
                            ntcallv_fh.write("\n")
                if args.ntvar == 1:
                        ntcallv_fh.close()

            else:
                ntcall_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tN\t\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")
                if args.ntvar == 1:
                    ntcallv_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tN\t\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")

                for Pos in sorted_Pos:
                    try:
                        total = coverage[Pos]
                    except:
                        total = 0
                        print(f"Coverage of position {Pos} not found")

                    if (total >= (sam_read_count * args.ntabund) or args.wgs == 1) and total >= args.ntcover:

                        try:
                            ins_nt_dict[Pos]
                        except:
                            pass
                        else:
                            i = 1
                            for insertion in ins_nt_dict[Pos]:
                                insert_position = Pos+(i/1000)
                                i_nts = insertion.split('insert')[1]
                                AA_pos = ((Pos-1)//3)+1
                                iabund = ins_nt_dict[Pos][insertion]/total
                                ntcall_lines['line'][insert_position] = f"{Pos}\t-\t\t\t\t\t\t\t{total}\t{i_nts}\t{ins_nt_dict[Pos][insertion]}"
                                ntcall_lines['line'][insert_position] += f"\t{iabund:.3f}"

                                if ( ins_nt_dict[Pos][insertion] >= args.min_count) and (iabund >= args.min_samp_abund):
                                    ntcall_lines['variant'][insert_position] = 1
                                i += 1

                        Pos_calls = {}
                        for key in nt_call_dict_dict[Pos]:
                            Pos_calls[key] = nt_call_dict_dict[Pos][key]
                        sorted_calls = sorted(Pos_calls, key=Pos_calls.__getitem__, reverse=True)

                        ntcall_lines['line'][Pos] = str(Pos)+"\t"+ref[1][Pos-1]
                        ntcall_lines['line'][Pos] += "\t"+str(nt_call_dict_dict[Pos]['A'])+"\t"+str(nt_call_dict_dict[Pos]['T'])+"\t"+str(nt_call_dict_dict[Pos]['C'])\
                            +"\t"+str(nt_call_dict_dict[Pos]['G'])+"\t"+str(nt_call_dict_dict[Pos]['-'])+"\t"+str(nt_call_dict_dict[Pos]['N'])
                        ntcall_lines['line'][Pos] += "\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[Pos][sorted_calls[0]])+"\t"+f"{(nt_call_dict_dict[Pos][sorted_calls[0]]/total):.3f}"
                        if not ref[1][Pos-1] == sorted_calls[0]:
                            ntcall_lines['variant'][Pos] = 1
                        if (nt_call_dict_dict[Pos][sorted_calls[1]] >= args.min_count) and (nt_call_dict_dict[Pos][sorted_calls[1]] / total) >= args.min_samp_abund:
                            ntcall_lines['variant'][Pos] = 1
                            ntcall_lines['line'][Pos] += "\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[Pos][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[Pos][sorted_calls[1]]/total):.3f}"
                            if (nt_call_dict_dict[Pos][sorted_calls[2]] > args.min_count) and (nt_call_dict_dict[Pos][sorted_calls[2]] /total > args.min_samp_abund):
                                ntcall_lines['variant'][Pos] = 1
                                ntcall_lines['line'][Pos] += "\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[Pos][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[Pos][sorted_calls[2]]/total):.3f}"

                for Pos in ntcall_lines['line']:
                    ntcall_fh.write(ntcall_lines['line'][Pos])
                    ntcall_fh.write("\n")
                    if args.ntvar == 1:
                        try:
                            ntcall_lines['variant'][Pos]
                        except:
                            pass
                        else:
                            ntcallv_fh.write(ntcall_lines['line'][Pos])
                            ntcallv_fh.write("\n")

            ntcall_fh.close()
            if args.ntvar == 1:
                ntcallv_fh.close()

    print(f'End {file} main output')

def gb_sam_parse(args, ref, file):
    """
    Called to handle main sam info parsing logic when using a genbank reference
    Parameters:
    args - argument values
    ref - tuple of referecne information
    file - name of sam file
    Functionality:
    calls sam_line_parser to get info from sam, then processes the information to get amino acid change information if applicable.
    prints outputs, mainly by function calls.  nt call output is the exception
    Returns nothing
    """

    samp=file.replace(file.split(".")[-1], '')[:-1]
    print(f"Starting {samp} processing")
    nt_call_dict_dict, ins_nt_dict, reads_list, col_reads, sam_read_count, coverage = sam_line_parser(args, ref, file)

    if sam_read_count == 0:
        print(f"No Reads for {samp}")
    else:

        indel_dict = {}
        aa_genome_pos_dict = {}
        nt_to_AA_dict = {}

        if args.AAreport == 1 and (args.seq == 1 or args.covar == 1 or args.indel == 1 or args.read == 1):
            for SNP_sequence in col_reads:
                if 'Reference' == SNP_sequence:
                    pass
                elif args.AAcodonasMNP == 1: # similar to fasta ref processing, but performs aa checks for each reading frame
                    try:
                        col_reads[SNP_sequence]['AA_sequence']
                    except:
                        mutorfs = {}
                        for orf in ref[3]:
                            try:
                                nt_to_AA_dict[orf]
                            except:
                                nt_to_AA_dict[orf] = {}
                            MNPs = []
                            curMNP = ''
                            last_codon = -1
                            fshift = 0
                            orflength = 0
                            for rf in ref[3][orf]['reading frames']:
                                for mut in SNP_sequence.split(' '):
                                    startPos = int(mut.split('-')[0].strip('ATCGN'))
                                    endPos = startPos
                                    if 'del' in mut:
                                        endPos = int(mut.split('-')[1].strip('del'))
                                    if startPos >= rf[0] and startPos <= rf[1]:
                                        orfstartpos = 1 + startPos - rf[0] + orflength
                                        orfendpos = 1 + endPos - rf[0] + orflength
                                        startcodon = (((orfstartpos)-1)//3)+1
                                        endcodon = ((orfendpos-1)//3)+1
                                        if curMNP:
                                            mutshift = 0
                                            if 'insert' in mut:
                                                mutshift += len(mut.split('|')[0].split('insert')[1])
                                            elif 'del' in mut:
                                                mutshift -= len(mut.split('-')[0].strip('0123456789'))
                                            if not mutshift % 3 == 0:
                                                if (fshift % 3) == 0:
                                                    if startcodon == last_codon:
                                                        curMNP.append([mut, orfstartpos])
                                                        fshift += mutshift
                                                    else:
                                                        MNPs.append(curMNP)
                                                        curMNP = [[mut, orfstartpos]]
                                                        fshift = mutshift

                                                else:
                                                    if startcodon < last_codon + 5:
                                                        curMNP.append([mut, orfstartpos])
                                                        fshift += mutshift
                                                    else:
                                                        ## splitfsMNP() make process
                                                        fsfreeMNP = []
                                                        for SNP in curMNP:
                                                            if 'fs' in SNP:
                                                                if fsfreeMNP:
                                                                    MNPs.append(fsfreeMNP)
                                                                    fsfreeMNP = []
                                                                MNPs.append([SNP])
                                                            else:
                                                                fsfreeMNP.append(SNP)
                                                        if fsfreeMNP:
                                                            MNPs.append(fsfreeMNP)
                                                        curMNP = [[mut, orfstartpos]]
                                                        fshift = mutshift
                                            else:
                                                if startcodon == last_codon:
                                                    curMNP.append([mut, orfstartpos])
                                                else:
                                                    if (fshift % 3) == 0:
                                                        MNPs.append(curMNP)
                                                    else:
                                                        ## splitfsMNP() make process
                                                        fsfreeMNP = []
                                                        for SNP in curMNP:
                                                            if 'fs' in SNP:
                                                                if fsfreeMNP:
                                                                    MNPs.append(fsfreeMNP)
                                                                    fsfreeMNP = []
                                                                MNPs.append([SNP])
                                                            else:
                                                                fsfreeMNP.append(SNP)
                                                        if fsfreeMNP:
                                                            MNPs.append(fsfreeMNP)
                                                    curMNP = [[mut, orfstartpos]]
                                                    fshift = 0
                                        else:
                                            curMNP = [[mut, orfstartpos]]
                                            if 'insert' in mut:
                                                fshift = len(mut.split('|')[0].split('insert')[1])
                                            elif 'del' in mut:
                                                fshift = len(mut.split('-')[0].strip('0123456789'))
                                        last_codon = endcodon
                                    elif endPos >= rf[0] and startPos <= rf[0]:
                                        if curMNP:
                                            MNPs.append(curMNP)
                                        MNPs.append([[mut, -1]])
                                        curMNP = ''
                                orflength += rf[1] - rf[0] + 1
                            if curMNP:
                                MNPs.append(curMNP)

                            if MNPs:
                                mutorfs[orf] = {}
                                for entry in MNPs:
                                    if len(entry) > 1:
                                        nt_to_AA_key = []
                                        for snpmut in entry:
                                            nt_to_AA_key.append(snpmut[0])
                                        nt_to_AA_key = ' '.join(nt_to_AA_key)
                                        if nt_to_AA_key in nt_to_AA_dict[orf]:
                                            for SNP in nt_to_AA_dict[orf][nt_to_AA_key]:
                                                mutorfs[orf][SNP] = nt_to_AA_dict[orf][nt_to_AA_key][SNP]
                                        else:
                                            MNPreport = ''
                                            mutntseq = ''
                                            orfstartpos = entry[0][1]
                                            orfendpos = entry[-1][1]
                                            if 'del' in entry[-1][0]:
                                                orfendpos += len(entry[-1][0].split('-')[0].strip('0123456789')) - 1

                                            cur_nts = {}

                                            for i in range(orfstartpos, orfendpos+1):
                                                cur_nts[i] = ref[3][orf]['nts'][i-1]
                                            for curPM in entry:
                                                if 'insert' in curPM[0]:
                                                    cur_nts[curPM[1]] = curPM[0].split('insert')[-1] + cur_nts[curPM[1]]
                                                elif 'del' in curPM[0]:
                                                    for i in range(curPM[1], curPM[1]+len(curPM[0].split('-')[0].strip('0123456789'))):
                                                        cur_nts[i] = ''
                                                else:
                                                    cur_nts[curPM[1]] = curPM[0][-1]
                                            for pos in cur_nts:
                                                mutntseq += cur_nts[pos]

                                            startcodonpos = ((orfstartpos-1)//3)+1
                                            endcodonpos = ((orfendpos-1)//3)+1
                                            startmod = orfstartpos % 3
                                            endmod = orfendpos % 3
                                            wtorfstartpos = orfstartpos
                                            wtorfendpos = orfendpos
                                            rf_end_del = 0
                                            if startmod == 2:
                                                wtorfstartpos = orfstartpos-1
                                                mutntseq = ref[3][orf]['nts'][orfstartpos-2] + mutntseq
                                            elif startmod == 0:
                                                wtorfstartpos = orfstartpos-2
                                                mutntseq = ref[3][orf]['nts'][orfstartpos-3:orfstartpos-1] + mutntseq
                                            if endmod == 2:
                                                wtorfendpos = orfendpos+1
                                                try:
                                                    mutntseq += ref[3][orf]['nts'][orfendpos]
                                                except:
                                                    rf_end_del = 1
                                            elif endmod == 1:
                                                wtorfendpos = orfendpos+2
                                                try:
                                                    mutntseq += ref[3][orf]['nts'][orfendpos:orfendpos+2]
                                                except:
                                                    rf_end_del = 1
                                            wtntseq = ref[3][orf]['nts'][wtorfstartpos-1:wtorfendpos]
                                            wtAAseq = ref[3][orf]['AAs'][startcodonpos-1:endcodonpos]
                                            mutAAseq = ''
                                            for i in range(0, (len(mutntseq)//3)):
                                                mutAAseq += aa_call(mutntseq[i*3:(i*3)+3])
                                            if len(mutntseq) < len(wtntseq):
                                                mutAAseq += 'del'
                                                mutntseq += 'del'
                                            elif len(mutntseq) > len(wtntseq):
                                                mutAAseq += 'insert'
                                                mutntseq += 'insert'
                                            if len(mutntseq) % 3 != 0:
                                                mutAAseq += 'fs'
                                            if rf_end_del == 1:
                                                mutAAseq += 'term/splice_distrupted'
                                                mutntseq += 'term/splice_distrupted'

                                            for curPM in entry:
                                                if startcodonpos == endcodonpos:
                                                    newmut = f"({orf}:{wtntseq}{wtorfstartpos}-{wtorfendpos}{mutntseq}({wtAAseq}{startcodonpos}{mutAAseq}))"
                                                else:
                                                    newmut = f"({orf}:{wtntseq}{wtorfstartpos}-{wtorfendpos}{mutntseq}({wtAAseq}{startcodonpos}-{endcodonpos}{mutAAseq}))"
                                                mutorfs[orf][curPM[0]] = newmut
                                                try:
                                                    nt_to_AA_dict[orf][nt_to_AA_key][curPM[0]] = newmut
                                                except:
                                                    nt_to_AA_dict[orf][nt_to_AA_key]= {curPM[0] : newmut}

                                    else:
                                        curPM = entry[0][0]
                                        if curPM in nt_to_AA_dict[orf]:
                                            newmut = nt_to_AA_dict[orf][curPM]
                                        else:
                                            mutstring = gb_snp_call(curPM, ref)
                                            newmut = []
                                            for submut in mutstring.split('|')[1:]:
                                                if orf in submut:
                                                    newmut.append(submut)
                                            newmut = '|'.join(newmut)
                                            nt_to_AA_dict[orf][curPM] = newmut
                                        mutorfs[orf][curPM] = newmut
                        MNPchecked = []
                        for mut in SNP_sequence.split(' '):
                            checkedmut = mut
                            for orf in mutorfs:
                                try:
                                    checkedmut += '|'+mutorfs[orf][mut]
                                except:
                                    pass
                            MNPchecked.append(checkedmut)
                            if args.indel == 1:
                                if 'del' in checkedmut or 'insert' in checkedmut:
                                    try:
                                        indel_dict[checkedmut] += sum(col_reads[SNP_sequence].values())
                                    except:
                                        indel_dict[checkedmut] = sum(col_reads[SNP_sequence].values())
                        mutations = " ".join(MNPchecked)
                        col_reads[SNP_sequence]['AA_sequence'] = mutations
                else:
                    try:
                        col_reads[SNP_sequence]['AA_sequence']
                    except:
                        mutations =[]
                        for mut in SNP_sequence.split(' '):
                            try:
                                newmut = nt_to_AA_dict[mut]
                            except:
                                mutstring = gb_snp_call(mut, ref)
                                nt_to_AA_dict[mut] = mutstring
                                newmut = mutstring
                            mutations.append(newmut)
                            if args.indel == 1:
                                if 'del' in mut or 'insert' in mut:
                                    try:
                                        indel_dict[newmut] += sum(col_reads[SNP_sequence].values())
                                    except:
                                        indel_dict[newmut] = sum(col_reads[SNP_sequence].values())
                        col_reads[SNP_sequence]['AA_sequence'] = " ".join(mutations)

        elif args.indel == 1:
            indel_dict = get_nt_indels(col_reads)

        # converts nt based variant lines into AA based lines
        if args.AAcentered == 1 and args.AAreport == 1:
            for SNP_sequence in col_reads:
                aa_seq = {}
                try:
                    col_reads[SNP_sequence]['AA_sequence']
                except:
                    pass
                else:
                    for mut in col_reads[SNP_sequence]['AA_sequence'].split(' '):
                        if '|' in mut:
                            split_mut = mut.split('|')

                            for orfmut in split_mut[1:]:
                                orfmut_split = orfmut[1:-1].split(':')
                                orf = orfmut_split[0]
                                split_orfmut_split = orfmut_split[1].split('(')
                                try:
                                    AAchange = split_orfmut_split[1][:-1]
                                    orfnts = split_orfmut_split[0]
                                except:
                                    AAchange = split_orfmut_split[0]
                                    orfnts = ''

                                entry = f"{orf}:{AAchange}({orfnts})"

                                startPos = int(split_mut[0].split('-')[0].strip('ATCGN'))
                                endPos = startPos
                                if 'del' in split_mut[0]:
                                    endPos = int(split_mut[0].split('-')[1].strip('del'))

                                try:
                                    aa_seq[entry]
                                except:
                                    aa_seq[entry] = [startPos, endPos]
                                else:
                                    aa_seq[entry].append(startPos)
                                    aa_seq[entry].append(endPos)

                        else:
                            startPos = int(mut.split('-')[0].strip('ATCGN'))
                            endPos = startPos
                            if 'del' in mut:
                                endPos = int(mut.split('-')[1].strip('del'))
                            aa_seq['Non-Coding:'+mut] = [startPos, endPos]


                    aa_sequence = []
                    for entry in aa_seq:
                        aa_sequence.append(entry)
                        aa_genome_pos_dict[entry] = [min(aa_seq[entry]), max(aa_seq[entry])]
                    aa_sequence = " ".join(aa_sequence)

                    col_reads[SNP_sequence]['AA_centered'] = aa_sequence

        if args.read == 1:
            print_reads(samp, reads_list, col_reads, args)

        if args.seq == 1:
            print_unique_seq(samp, sam_read_count, col_reads, coverage, args)

        if args.indel == 1 and len(indel_dict) > 0:
            print_indels(samp, sam_read_count, indel_dict, coverage, args)

        if args.covar == 1:
            print_covars(samp, sam_read_count, col_reads, coverage, args, aa_genome_pos_dict)

        if args.nt_call == 1:
            ntcall_lines = {'line' : {},
                            'variant' : {}
                            }
            ntcall_fh = open(samp+'_nt_calls.tsv', "w")
            ntcall_fh.write(samp+"("+str(sam_read_count)+")\n")
            if args.ntvar == 1:
                ntcallv_fh = open(samp+'_nt_calls_varonly.tsv', "w")
                ntcallv_fh.write(samp+"("+str(sam_read_count)+")\n")
            sorted_Pos = sorted(nt_call_dict_dict)
            if args.AAreport == 1:

                ntcall_fh.write("Position\tref NT\tAAs\tA\tT\tC\tG\t-\tN\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                if args.ntvar == 1:
                    ntcallv_fh.write("Position\tref NT\tAAs\tA\tT\tC\tG\t-\tN\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                consensus = {}
                OrfPosDict = {}
                for orf in ref[3]:
                    orflength = 0
                    for rf in ref[3][orf]['reading frames']:
                        for i in range(rf[0],rf[1]+1):
                            orfPos = i-rf[0]+1+orflength
                            try:
                                OrfPosDict[i].append([orf , orfPos, ref[3][orf]['AAs'][(orfPos-1)//3] ,(((orfPos-1)//3)+1)])
                            except:
                                OrfPosDict[i] = [[orf , orfPos, ref[3][orf]['AAs'][(orfPos-1)//3] ,(((orfPos-1)//3)+1)]]
                        orflength += rf[1] - rf[0] + 1

                ORFmismatch = {}
                for Pos in sorted_Pos:
                    try:
                        total = coverage[Pos]
                    except:
                        total = 0
                        print(f"coverage of position {Pos} not found")

                    if (total >= (sam_read_count * args.ntabund) or args.wgs == 1) and total >= args.ntcover:


                        orfAAreports = []
                        try:
                            OrfPosDict[Pos]
                        except:
                            pass
                        else:
                            for entry in OrfPosDict[Pos]:
                                orfAAreports.append(entry[0] +"_nt:"+ str(entry [1]) +"_AA:"+ entry[2] + str(entry[3]))

                        try:
                            ins_nt_dict[Pos]
                        except:
                            pass
                        else:
                            i = 1
                            for insertion in ins_nt_dict[Pos]:
                                insert_position = Pos+(i/1000)
                                i_nts = insertion.split('insert')[1]
                                try:
                                    i_AAs = '|'.join(gb_snp_call(insertion, ref).split('|')[1:])
                                except:
                                    i_AAs = ''
                                AA_pos = ", ".join(orfAAreports)
                                iabund = ins_nt_dict[Pos][insertion]/total
                                ntcall_lines['line'][insert_position] = f"{Pos}\t-\t{AA_pos}\t\t\t\t\t\t\t{total}\t{i_nts}\t{ins_nt_dict[Pos][insertion]}"
                                ntcall_lines['line'][insert_position] += f"\t{iabund:.3f}\t\t{i_AAs}"
                                ntcall_lines['line'][insert_position] += "\n"
                                i += 1

                        Pos_calls = {}
                        for key in nt_call_dict_dict[Pos]:
                            Pos_calls[key] = nt_call_dict_dict[Pos][key]
                        sorted_calls = sorted(Pos_calls, key=Pos_calls.__getitem__, reverse=True)

                        ntcall_lines['line'][Pos] =(str(Pos)+"\t"+ref[1][Pos-1]+"\t")
                        ntcall_lines['line'][Pos] +=(", ".join(orfAAreports))

                        ntcall_lines['line'][Pos] += "\t"+str(nt_call_dict_dict[Pos]['A'])+"\t"+str(nt_call_dict_dict[Pos]['T'])+"\t"+str(nt_call_dict_dict[Pos]['C'])
                        ntcall_lines['line'][Pos] += "\t"+str(nt_call_dict_dict[Pos]['G'])+"\t"+str(nt_call_dict_dict[Pos]['-'])+"\t"+str(nt_call_dict_dict[Pos]['N'])
                        ntcall_lines['line'][Pos] +=("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[Pos][sorted_calls[0]]))
                        ntcall_lines['line'][Pos] +=(f"\t{(nt_call_dict_dict[Pos][sorted_calls[0]]/total):.3f}")

                        consensus[Pos] = sorted_calls
                        if consensus[Pos][0] != ref[1][Pos-1]:
                            try:
                                for entry in OrfPosDict[Pos]:
                                    try:
                                        ORFmismatch[entry[0]]
                                    except:
                                        ORFmismatch[entry[0]] = {entry[1] : sorted_calls[0]}
                                    else:
                                        ORFmismatch[entry[0]][entry[1]] = sorted_calls[0]
                            except:
                                pass

                for Pos in sorted_Pos:
                    try:
                        total = coverage[Pos]
                    except:
                        total = 0

                    try:
                        ntcall_lines['line'][Pos]
                    except:
                        pass
                    else:
                        if consensus[Pos][0] != ref[1][Pos-1]:
                            ntcall_lines['variant'][Pos] = 1
                            try:
                                OrfPosDict[Pos]
                            except:
                                if (nt_call_dict_dict[Pos][consensus[Pos][1]] > args.min_count) and (nt_call_dict_dict[Pos][consensus[Pos][1]] /total > args.min_samp_abund):
                                    ntcall_lines['line'][Pos] +=(f"\t\t\t{consensus[Pos][1]}\t{nt_call_dict_dict[Pos][consensus[Pos][1]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][1]]/total):.3f}")
                                    if nt_call_dict_dict[Pos][consensus[Pos][2]] > args.min_count and nt_call_dict_dict[Pos][consensus[Pos][2]] /total  > args.min_samp_abund:
                                        ntcall_lines['line'][Pos] +=(f"\t\t{consensus[Pos][2]}\t{nt_call_dict_dict[Pos][consensus[Pos][2]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][2]]/total):.3f}")
                            else:
                                orfAAreports = [[],[]]
                                for ORFentry in OrfPosDict[Pos]:
                                    OrfPos = ORFentry[1]
                                    orf = ORFentry[0]
                                    mod = (OrfPos)%3
                                    codon = ['n','n','n']

                                    if mod == 0:
                                        codon[2] = consensus[Pos][0]
                                        try:
                                            codon[0] = ORFmismatch[orf][OrfPos-2]
                                        except:
                                            codon[0] = ref[3][orf]['nts'][OrfPos-3]
                                        try:
                                            codon[1] = ORFmismatch[orf][OrfPos-1]
                                        except:
                                            codon[1] = ref[3][orf]['nts'][OrfPos-2]

                                    elif mod == 2:
                                        codon[1] = consensus[Pos][0]
                                        try:
                                            codon[2] = ORFmismatch[orf][OrfPos+1]
                                        except:
                                            codon[2] = ref[3][orf]['nts'][OrfPos]
                                        try:
                                            codon[0] = ORFmismatch[orf][OrfPos-1]
                                        except:
                                            codon[0] = ref[3][orf]['nts'][OrfPos-2]

                                    elif mod == 1:
                                        codon[0] = consensus[Pos][0]
                                        try:
                                            codon[2] = ORFmismatch[orf][OrfPos+2]
                                        except:
                                            codon[2] = ref[3][orf]['nts'][OrfPos+1]
                                        try:
                                            codon[1] = ORFmismatch[orf][OrfPos+1]
                                        except:
                                            codon[1] = ref[3][orf]['nts'][OrfPos]
                                    orfAAreports[0].append(orf+"_"+aa_call("".join(codon)))
                                    orfAAreports[1].append(orf+"_"+singlet_codon_call(OrfPos, consensus[Pos][0], ref[3][orf]['nts'])[1])
                                ntcall_lines['line'][Pos] +=("\t"+", ".join(orfAAreports[0])+"\t"+", ".join(orfAAreports[1]))
                                if (nt_call_dict_dict[Pos][consensus[Pos][1]] > args.min_count) and (nt_call_dict_dict[Pos][consensus[Pos][1]] /total > args.min_samp_abund):
                                        ntcall_lines['line'][Pos] +=(f"\t{consensus[Pos][1]}\t{nt_call_dict_dict[Pos][consensus[Pos][1]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][1]]/total):.3f}")
                                        orfAAreports = []
                                        for ORFentry in OrfPosDict[Pos]:
                                            OrfPos = ORFentry[1]
                                            orf = ORFentry[0]
                                            orfAAreports.append(orf+"_"+singlet_codon_call(OrfPos, consensus[Pos][1], ref[3][orf]['nts'])[1])
                                        ntcall_lines['line'][Pos] += ("\t"+", ".join(orfAAreports))
                                        if nt_call_dict_dict[Pos][consensus[Pos][2]] > args.min_count:
                                            if nt_call_dict_dict[Pos][consensus[Pos][2]] /total  > args.min_samp_abund:
                                                ntcall_lines['line'][Pos] +=(f"\t{consensus[Pos][2]}\t{nt_call_dict_dict[Pos][consensus[Pos][2]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][2]]/total):.3f}")
                                                orfAAreports = []
                                                for ORFentry in OrfPosDict[Pos]:
                                                    OrfPos = ORFentry[1]
                                                    orf = ORFentry[0]
                                                    orfAAreports.append(orf+"_"+singlet_codon_call(OrfPos, consensus[Pos][2], ref[3][orf]['nts'])[1])
                                                ntcall_lines['line'][Pos] += ("\t"+", ".join(orfAAreports))


                        elif (nt_call_dict_dict[Pos][consensus[Pos][1]] >= args.min_count) and ((nt_call_dict_dict[Pos][consensus[Pos][1]] / total) >= args.min_samp_abund):
                            ntcall_lines['variant'][Pos] = 1

                            ntcall_lines['line'][Pos] +=("\t\t")
                            ntcall_lines['line'][Pos] +=(f"\t{consensus[Pos][1]}\t{nt_call_dict_dict[Pos][consensus[Pos][1]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][1]]/total):.3f}")

                            try:
                                OrfPosDict[Pos]
                            except:
                                if nt_call_dict_dict[Pos][consensus[Pos][2]] /total  > args.min_samp_abund:
                                    ntcall_lines['line'][Pos] +=(f"\t\t{consensus[Pos][2]}\t{nt_call_dict_dict[Pos][consensus[Pos][2]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][2]]/total):.3f}")
                            else:
                                orfAAreports = []
                                for ORFentry in OrfPosDict[Pos]:
                                    OrfPos = ORFentry[1]
                                    orf = ORFentry[0]
                                    orfAAreports.append(orf+"_"+singlet_codon_call(OrfPos, consensus[Pos][1], ref[3][orf]['nts'])[1])
                                ntcall_lines['line'][Pos] += ("\t"+", ".join(orfAAreports))
                                if nt_call_dict_dict[Pos][consensus[Pos][2]] > args.min_count:
                                    if nt_call_dict_dict[Pos][consensus[Pos][2]] /total  > args.min_samp_abund:
                                        ntcall_lines['line'][Pos] +=(f"\t{consensus[Pos][2]}\t{nt_call_dict_dict[Pos][consensus[Pos][2]]}\t{(nt_call_dict_dict[Pos][consensus[Pos][2]]/total):.3f}")
                                        orfAAreports = []
                                        for ORFentry in OrfPosDict[Pos]:
                                            OrfPos = ORFentry[1]
                                            orf = ORFentry[0]
                                            orfAAreports.append(orf+"_"+singlet_codon_call(OrfPos, consensus[Pos][2], ref[3][orf]['nts'])[1])
                                        ntcall_lines['line'][Pos] += ("\t"+", ".join(orfAAreports))

                        ntcall_lines['line'][Pos] +=("\n")
                for Pos in ntcall_lines['line']:
                    ntcall_fh.write(ntcall_lines['line'][Pos])
                    if args.ntvar == 1:
                        try:
                            ntcall_lines['variant'][Pos]
                        except:
                            pass
                        else:
                            ntcallv_fh.write(ntcall_lines['line'][Pos])
                if args.ntvar == 1:
                        ntcallv_fh.close()
            else:
                ntcall_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tN\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")
                if args.ntvar == 1:
                    ntcallv_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tN\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")

                for Pos in sorted_Pos:
                    try:
                        total = coverage[Pos]
                    except:
                        total = 0
                        print(f"coverage of position {Pos} not found")

                    if (total >= (sam_read_count * args.ntabund) or args.wgs == 1) and total >= args.ntcover:

                        try:
                            ins_nt_dict[Pos]
                        except:
                            pass
                        else:
                            i = 1
                            for insertion in ins_nt_dict[Pos]:
                                insert_position = Pos+(i/1000)
                                i_nts = insertion.split('insert')[1]
                                iabund = ins_nt_dict[Pos][insertion]/total
                                ntcall_lines['line'][insert_position] = f"{Pos}\t-\t\t\t\t\t\t\t{total}\t{i_nts}\t{ins_nt_dict[Pos][insertion]}"
                                ntcall_lines['line'][insert_position] += f"\t{iabund:.3f}"
                                ntcall_lines['line'][insert_position] += "\n"
                                i += 1

                        Pos_calls = {}
                        for key in nt_call_dict_dict[Pos]:
                            Pos_calls[key] = nt_call_dict_dict[Pos][key]
                        sorted_calls = sorted(Pos_calls, key=Pos_calls.__getitem__, reverse=True)

                        ntcall_lines['line'][Pos] = str(Pos)+"\t"+ref[1][Pos-1]
                        ntcall_lines['line'][Pos] += "\t"+str(nt_call_dict_dict[Pos]['A'])+"\t"+str(nt_call_dict_dict[Pos]['T'])+"\t"+str(nt_call_dict_dict[Pos]['C'])+"\t"+str(nt_call_dict_dict[Pos]['G'])
                        ntcall_lines['line'][Pos] += "\t"+str(nt_call_dict_dict[Pos]['-'])+"\t"+str(nt_call_dict_dict[Pos]['N'])
                        ntcall_lines['line'][Pos] += "\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[Pos][sorted_calls[0]])+"\t"+f"{(nt_call_dict_dict[Pos][sorted_calls[0]]/total):.3f}"

                        if (nt_call_dict_dict[Pos][sorted_calls[1]] > args.min_count) and (nt_call_dict_dict[Pos][sorted_calls[1]] /total > args.min_samp_abund):
                            ntcall_lines['line'][Pos] += "\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[Pos][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[Pos][sorted_calls[1]]/total):.3f}"
                            if (nt_call_dict_dict[Pos][sorted_calls[2]] > args.min_count) and (nt_call_dict_dict[Pos][sorted_calls[2]] /total > args.min_samp_abund):
                                ntcall_lines['line'][Pos] += "\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[Pos][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[Pos][sorted_calls[2]]/total):.3f}"

                        ntcall_lines['line'][Pos] += "\n"
                for Pos in ntcall_lines['line']:
                    ntcall_fh.write(ntcall_lines['line'][Pos])
                    if args.ntvar == 1:
                        try:
                            ntcall_lines['variant'][Pos]
                        except:
                            pass
                        else:
                            ntcallv_fh.write(ntcall_lines['line'][Pos])
            ntcall_fh.close()
            if args.ntvar == 1:
                ntcallv_fh.close()

    nt_call_dict_dict, ins_nt_dict, reads_list, col_reads, sam_read_count, coverage = (None, None, None, None, None, None)
    print(f'End {file} main output')

def covar_deconv(args, samp, covariance_dict, sequence_dict):
    """
    Called to perform chimera removal based on the covars output and print the covar deconv output
    Parameters:
    args - argument values
    samp - sam sample name
    covariance_dict - dict of covar lines
    sequence_dict - dict of unique seq lines

    Functionality:
    A sequence is first determined to likely be a true or chimeric sequence.  The ratio of the frequency of a given covariant sequence
    relative to the product of the abundances of each individual polymorphism that is present in that covariant sequence is calculated.
    If the ratio of the sequence abundance to the product is equal to or greater than 1 (beta), that covariant passes the check.
    Any sequence that has an abundance of 0.3 or greater is automatically passed.  Then the passed sequences, in order of greatest
    observed/expected ratio to least, are assigned a new occurrence count based on their constituent individual polymorphisms. The count
    of the least abundant individual polymorphism is assigned to the sequence and constituent polymorphisms making up the sequence have
    their count reduced by the amount of the least abundant polymorphism. Any sequence not yet processed in which that polymorphism is
    present is removed. This process is repeated until all sequences have been reassessed or removed.

    Returns nothing
    """

    passedseqs = {}
    preservedseqs = {}
    covar_dict = covariance_dict
    for seq in sequence_dict: # pass check for actual : expected abundance
        if seq != 'total' and seq != 'singles':

            splitseq = seq.split(' ')
            abund = 1
            for sing in splitseq:
                try:
                    abund = abund * (covar_dict[sing] / covar_dict['total'])
                except:
                    abund = abund * (sequence_dict[seq] / sequence_dict['total'])

            try:
                covarabund = covar_dict[seq]/covar_dict['total']
            except:
                covarabund = sequence_dict[seq]/sequence_dict['total']
                covar_dict[seq] = sequence_dict[seq]

            if covarabund >= args.autopass:
                preservedseqs[seq] = max(1, args.beta, (covarabund / abund))
            elif covarabund >= abund * args.beta:
                passedseqs[seq] = covarabund / abund
            elif len(seq.split(' ')) == 1:
                passedseqs[seq] = max(1, args.beta)

    if args.min_samp_abund < 1:
        min_count = args.min_samp_abund * covar_dict['total']
    else:
        min_count = args.min_samp_abund

    if args.pass_out == 1: # write passed covars to file if enabled
        fh_pass = open(samp+"_covar_pass.tsv", "w")
        fh_pass.write(f"{samp}({covar_dict['total']})\nSequences\tCount\tAbundance\tPass Ratio\n")
        for seq in preservedseqs:
            if covar_dict[seq] >= min_count:
                fh_pass.write(f"{seq}\t{covar_dict[seq]}\t{(covar_dict[seq]/covar_dict['total']):.3f}\t{preservedseqs[seq]}*\n")
        for seq in passedseqs:
            if covar_dict[seq] >= min_count:
                fh_pass.write(f"{seq}\t{covar_dict[seq]}\t{(covar_dict[seq]/covar_dict['total']):.3f}\t{passedseqs[seq]}\n")
        fh_pass.close()

    # sort passed covars
    lensortedpassed = sorted(passedseqs, key = lambda key : len(key.split(' ')), reverse=True)
    ratiolensortedpassed = sorted(lensortedpassed, key = lambda key : passedseqs[key], reverse = True)
    sortedsingles = sorted(covar_dict['singles'], key = covar_dict['singles'].__getitem__)
    deconved = {}
    for seq in ratiolensortedpassed: # reassign counts
        singles = seq.split(' ')
        first = 0
        rem_count = 0
        for sing in sortedsingles:
            if sing in singles:
                if covar_dict['singles'][sing] > 0:
                    if first == 0:
                        first = 1
                        rem_count = covar_dict['singles'][sing]
                        covar_dict['singles'][sing] = 0
                        deconved[seq] = rem_count
                    else:
                        covar_dict['singles'][sing] = covar_dict['singles'][sing] - rem_count
                else:
                    break
        sortedsingles = sorted(covar_dict['singles'], key = covar_dict['singles'].__getitem__)

    sortedpreserved = sorted(preservedseqs, key = lambda key : covar_dict[key])

    for seq in sortedpreserved:
        singles = seq.split(' ')
        first = 0
        rem_count = 0
        for sing in sortedsingles:
            if sing in singles:
                if covar_dict['singles'][sing] > 0:
                    if first == 0:
                        first = 1
                        rem_count = covar_dict['singles'][sing]
                        covar_dict['singles'][sing] = 0
                        deconved[seq] = rem_count
                    else:
                        covar_dict['singles'][sing] = covar_dict['singles'][sing] - rem_count
                else:
                    break
        sortedsingles = sorted(covar_dict['singles'], key = covar_dict['singles'].__getitem__)


    newtotal = sum(deconved.values())
    fh_deconv = open(samp+"_covar_deconv.tsv", "w")
    fh_deconv.write(f"{samp}({covar_dict['total']}) | ({newtotal})\nSequences\tCount\tAbundance\n")
    sorted_deconved = sorted(deconved, key = deconved.__getitem__, reverse = True)
    for seq in sorted_deconved: # write deconv
        if deconved[seq] >= min_count:
            fh_deconv.write(f"{seq}\t{deconved[seq]}\t{(deconved[seq]/newtotal):.3f}\n")
    fh_deconv.close()

    print(f"End covar deconv out for {samp}") # END COVAR DECONV OUT

    return()

def dechim(args, seqs):
    """
    Called to perform chimera removal based on the unique seqs output
    Parameters:
    args - argument values
    seqs - dict of unique seq lines
    Functionality:
    The individual unique sequences, starting with the lowest abundance, are broken up into all possible dimeric halves. Each pair
    is then compared to all the other sequences to detect potential parents. A sequence is flagged as a potential parent if its
    abundance is greater than or equal to the abundance of the potential chimera multiplied by foldab and there is at least one
    other sequence that would be a matched parent to the complimentary dimeric half. When a pair of dimeric halves have potential
    parents, the abundances of parent pairs are multiplied. The products from each potential parent pairings are summed as an
    expected abundance value and compared to the observed abundance of the potential chimera. If the abundance of the potential
    chimera is less than that of the expected value multiplied by alpha, that sequence is flagged as a chimera and removed. The
    counts attributed to the flagged chimeric sequence are then redistributed to the parent sequences based on the relative expected
    contribution to recombination.
    Returns new sequence dict without found chimeras
    """


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
                                left_break_Pos = ''
                                for c in left_break:
                                    if c.isdigit():
                                        left_break_Pos += c
                                    else:
                                        if left_break_Pos:
                                            break

                                rt_break_Pos = ''
                                for c in rt_break:
                                    if c.isdigit():
                                        rt_break_Pos += c
                                    else:
                                        if rt_break_Pos:
                                            break

                                if int(left_break_Pos) > int(rt_break_Pos):
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

    seqs['total'] = total

    return(seqs)

def chim_rm(args, samp, sequences):
    """
    Called to reiterate the chim_rm chimera removal process based on unique seqs and print the chim_rm output
    Parameters:
    args - argument values
    samp - sample name from sam file
    sequencess - dict of unique sequences
    Functionality:
    Looped calling of chimera removal until no more chimeras are found or limit hit.  Results are then printed
    Returns nothing
    """

    seqs = sequences
    pre_len = len(seqs)
    inf_loop_shield = 0
    while True: # send sequences for chimera removal while chimeras are still found
        dechim(args, seqs)
        post_len = len(seqs)
        inf_loop_shield += 1
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
    print(f"End chim_rm out for {samp}")
    return()

def chim_process(args, samp):
    """
    Called to collect sequences from covar and unique seq files, then pass those to the chimera removal functions
    Parameters:
    args - argument values
    samp - samples name from sam file
    Functionality:
    Tries to open files based on based sample name.  If successful, calls the specific chimera removal
    functions
    Returns nothing
    """
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

    if args.deconv == 1:
        in_covars = {}
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
                covar_deconv(args, samp, in_covars, in_seqs)
        except:
            print(f"Reading of {samp}_covars.tsv failed")

    if args.chim_rm == 1:
        chim_rm(args, samp, in_seqs)

def collect_lines(dict_dict, all_dict, file, args):
    """
    Called to collect the lines of a file
    Parameters:
    dict_dict - dict of dict of lines from files
    all_dict - dict of all unique sequences from files
    file - name of the file to collect from
    args - argument values
    Functionality:
    Parses through the file to collect sequences and info
    Returns updated dicts
    """

    sample_line = ''
    try:
        samp=open(file, "r")
    except:
        print("can't open "+file)
    else:
        for line in samp:
            split_line = line.strip("\n\r").split("\t")
            try:
                split_line[1]
            except:
                sample_line = split_line[0]
                dict_dict[sample_line] = {}

            else:
                if not split_line[1] == 'Count':
                    if float(split_line[2]) >= args.min_col_abund:
                        dict_dict[sample_line][split_line[0]] = [split_line[1], split_line[2]]
                        all_dict[split_line[0]] = 1
        samp.close()
        return(dict_dict, all_dict)

def collection(args):
    """
    Called to perform collection logic and printing
    Parameters:
    args - argument vlaues
    Functionality:
    Looks through directory for sample files to collect sequence info from and calls function to collect info.
    Then prints collected info to new files.
    Returns
    """

    covar_dict_dict = {}
    seq_dict_dict = {}
    deconv_dict_dict = {}
    pass_dict_dict = {}
    cr_dict_dict = {}

    all_covar = {}
    all_seq = {}
    all_deconv = {}
    all_pass = {}
    all_cr = {}

    for file in os.listdir(os.getcwd()):
        if file.endswith('_covars.tsv'):
            covar_dict_dict, all_covar = collect_lines(covar_dict_dict, all_covar, file, args)

        if file.endswith('_seqs.tsv'):
            seq_dict_dict, all_seq = collect_lines(seq_dict_dict, all_seq, file, args)

        if file.endswith('_deconv.tsv'):
            deconv_dict_dict, all_deconv = collect_lines(deconv_dict_dict, all_deconv, file, args)

        if file.endswith('_pass.tsv'):
            pass_dict_dict, all_pass = collect_lines(pass_dict_dict, all_pass, file, args)

        if file.endswith('_chim_rm.tsv'):
            cr_dict_dict, all_cr = collect_lines(cr_dict_dict, all_cr, file, args)

    dict_collection = (
    ("Covariances", covar_dict_dict, all_covar),
    ("Unique_Seqs", seq_dict_dict, all_seq),
    ("Covar_Deconv", deconv_dict_dict, all_deconv),
    ("Covar_Pass", pass_dict_dict, all_pass),
    ("Chimeras_Removed", cr_dict_dict, all_cr)
    )

    for dicts in dict_collection:
        if dicts[1]:
            if args.colID == '':
                collection_fh = open(f"Collected_{dicts[0]}.tsv","w")
            else:
                collection_fh = open(args.colID+f"_Collected_{dicts[0]}.tsv","w")
            sorted_seqs = sorted(dicts[2])
            collection_fh.write("\t")
            for sampline in dicts[1]:
                collection_fh.write(sampline+"\t\t")
            collection_fh.write("\nSequences\t")
            for sampline in dicts[1]:
                collection_fh.write("Count\tAbundance\t")
            collection_fh.write("\n")
            for seq in sorted_seqs:
                collection_fh.write(seq+"\t")
                for sample in dicts[1]:
                    try:
                        collection_fh.write(dicts[1][sample][seq][0]+"\t"+dicts[1][sample][seq][1]+"\t")
                    except:
                        collection_fh.write("\t\t")

                collection_fh.write("\n")
            collection_fh.close()
        else:
            print(f"No {dicts[0]} files found")

def main():
    """
    Main functional logic
    Calls functions to get argument values, then as relevant, get reference, process sam files, perform chimera removal and collect sample's
    """
    args = arg_parser() # getting command line arguments

    if args.ref:
        ref = get_ref(args) # get the reference ID and sequence from the FASTA file
        if ref[1] == '':
            print('Reference not recognized as a Fasta or Genebank format, skipping SAM parsing')
        else:

            # collect SAM files to process, either from the arguments or the working directory
            SAMs = []
            try:
                args.Sam_files[0]
            except:
                for file in os.listdir(os.getcwd()):
                    if file.lower().endswith('.sam'):
                        SAMs.append(file)
                    elif bc and (file.lower().endswith('.bam') or file.lower().endswith('.cram')):
                        SAMs.append(file)
            else:
                for files in args.Sam_files:
                    for file in files:
                        if os.path.isfile(file):
                            if file.lower().endswith(".sam") or (bc and (file.lower().endswith('.bam') or file.lower().endswith('.cram'))):
                                SAMs.append(file)
                        else:
                            print(f"Can't find {file}, skipping")

            args.ref = ''
            if ref[2] == 'fasta':
                with Pool(processes=args.mp) as pool:
                    pool.starmap(fa_sam_parse, zip(itertools.repeat(args), itertools.repeat(ref), SAMs))
            elif ref[2] == 'gb':
                if len(SAMs) == 1:
                    gb_sam_parse(args, ref, SAMs[0])
                else:
                    with Pool(processes=args.mp) as pool:
                        pool.starmap(gb_sam_parse, zip(itertools.repeat(args), itertools.repeat(ref), SAMs))
            print(f"End Sam Parsing Output")
    else:
        print('No reference provided, skipping SAM parsing')


    seq_files = []
    # Begin chimera removal if enabled
    if args.chim_rm == 1 or args.deconv == 1:
        for file in os.listdir(os.getcwd()):
            if file.endswith('_seqs.tsv'): # get unique sequence files for chimera removal
                seq_files.append(file[0:-16])
        with Pool(processes=args.mp) as pool:
            pool.starmap(chim_process, zip(itertools.repeat(args), seq_files))

# begin collection of sample outputs
    if args.collect == 1:
        collection(args)
    exit()

if __name__ == '__main__':
    main()