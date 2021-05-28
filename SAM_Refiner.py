#!/bin/env python3
# Writen by Devon Gregory with assistance from Christopher Bottoms
# University of Missouri
# Licence TBD
import os
import sys
import argparse
import json
import itertools
from multiprocessing import Process
from pathlib import Path

DEBUG = False


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
        type=argparse.FileType('r'),
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
        '--min_abundance1',
        type=float,
        default=0.001,
        help='Minimum observations required to be included in sample reports; >= 1 occurance count; < 1 %% observed (.1 = 10%%), (default: .001)'
    )
    parser.add_argument(
        '--min_abundance2',
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
        '--max_dist',
        type=int,
        default=50,
        help='Maximum number of variances from the reference a sequence can have to be consider in covars processing (default: 50)'
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
        '--sams',
        type=int,
        default=1,
        choices=[0, 1],
        help='Enable/Disable (1/0) sam processing, default enabled (--sams 1)'
    )
    parser.add_argument(
        '--nt_call',
        type=int,
        default=1,
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

    args = parser.parse_args()


    if args.wgs == 1:
        if args.deconv == 1 or args.chim_rm == 1:
            args.deconv = 0
            args.chim_rm = 0
            print('WGS mode enabled, disabling chimera removal methods')

    if args.sams == 1:
        if not args.ref:
            print('Must have a reference(-r, --referecne) to parse sams, skipping sam parsing')
            args.sams = 0

    if args.min_abundance1 < 0:
        print(f"--min_abundance1 must be non-negative, defaulting to 0.001")
        args.min_abundance1=0.001

    if args.min_abundance2 < 0:
        print(f"--min_abundance2 must be non-negative and < 1, defaulting to 01")
        args.min_abundance2=0.01
    elif args.min_abundance2 >= 1:
        print(f"--min_abundance2 must be non-negative and < 1, defaulting to 01")
        args.min_abundance2=0.01

    if args.ntabund < 0:
        print(f"--ntabund must be non-negative and < 1, defaulting to .001")
        args.ntabund=0.001
    elif args.ntabund >= 1:
        print(f"--ntabund must be non-negative and < 1, defaulting to .001")
        args.ntabund=0.001

    if args.max_dist < 0:
        print(f"--max_dist must be non-negative, defaulting to 50")
        args.max_dist=50

    if args.max_covar < 0:
        print(f"--max_covar must be non-negative, defaulting to 8")
        args.max_dist=8

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

def get_ref():

    n=0
    refID = ''
    refseq = ''
    for line in args.ref:
        if line.startswith('>'):
            n+=1
            if n > 1:
                break
            refID = line[1:].strip("\n\r")
        elif n == 1:
            refseq = refseq + line.strip("\n\r")


    return(refID, refseq)

def AAcall(codon):
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
        'GGG' : 'G'
    }
    AA = '?'
    if codon in AAdict:
        AA = AAdict[codon]

    return(AA)

def singletCodon(ntPOS, nt):
    AAPOS = (ntPOS-1)//3
    AAmod = (ntPOS-1)%3
    codon = ""

    if AAmod == 0:
        codon = nt+ref[1][ntPOS]+ref[1][ntPOS+1]
    elif AAmod == 1:
        codon = ref[1][ntPOS-2]+nt+ref[1][ntPOS]
    elif AAmod == 2:
        codon = ref[1][ntPOS-3]+ref[1][ntPOS-2]+nt

    return(AAPOS+1, AAcall(codon))

def getCombos(qlist, clen):
    combos = []
    if (clen == 0 or clen > len(qlist)):
        clen = len(qlist)
    for N in range(1, clen+1):
        for comb in itertools.combinations(qlist, N):
            combos.append(' '.join(comb))
    return(combos)

def SAMparse(file):
    samp=file.name[0: -4]
    print(f"Starting {samp} processing")
    nt_call_dict_dict = {}
    indel_dict = {}
    seq_species = {}
    sam_read_count = 0
    sam_line_count = 0
    firstPOS = 0
    lastPOS = 0
    coverage = {}


    for line in file:
        if not line.startswith('@'):
            splitline = line.split("\t")
            if ref[0].startswith(splitline[2]):
                if int(splitline[4]) > 0:

                    abund_count=1
                    if args.use_count == 1:
                        if '-' in splitline[0] and '=' in splitline[0]:
                            eq_split = splitline[0].split('=')
                            dash_split = splitline[0].split('-')
                            if len(eq_split[-1]) > len(dash_split[-1]):
                                abund_count=int(dash_split[-1])
                            else:
                                abund_count=int(eq_split[-1])

                        elif '-' in splitline[0]:
                            try:
                                abund_count=int(splitline[0].split('-')[-1])
                            except:
                                # print(splitline[0])
                                pass

                        elif '=' in splitline[0]:
                            try:
                                abund_count=int(splitline[0].split('=')[-1])
                            except:
                                pass


                    sam_read_count += abund_count
                    sam_line_count += 1

                    if sam_read_count % 5000 == 0:
                        print(f"At read {sam_read_count} of {samp}")
                    if sam_line_count % 5000 == 0:
                        print(f"At line {sam_line_count} of {samp} SAM")

                    CIGAR = splitline[5]
                    POS = int(splitline[3])
                    if firstPOS == 0:
                        firstPOS = POS
                    elif POS < firstPOS:
                        firstPOS = POS

                    query_seq = splitline[9]
                    run_length = 0
                    query_seq_parsed = ''
                    query_pos = 0
                    q_pars_pos = 0
                    mutations = []

                    for C in CIGAR:
                        if C == 'M' or C == 'I' or C == 'D' or C == 'S' or C == 'H':
                            if C == 'S':
                                query_pos = query_pos + run_length
                            if C == 'H':
                                query_pos = query_pos + run_length
                                q_pars_pos = q_pars_pos + run_length

                            if C == 'I':
                                if query_pos > 0:
                                    # add insertion to array
                                    iPOS = q_pars_pos+POS

                                    iSeq = query_seq[query_pos: query_pos+run_length]

                                    try:
                                        indel_dict[str(iPOS)+'-insert'+iSeq]
                                    except:
                                        indel_dict[str(iPOS)+'-insert'+iSeq] = abund_count
                                    else:
                                        indel_dict[str(iPOS)+'-insert'+iSeq] += abund_count

                                    mutations.append(str(iPOS)+'-insert'+iSeq)

                                    query_pos = query_pos + run_length

                            elif C == 'D':
                                for X in range(0, run_length):
                                    query_seq_parsed += '-'

                                mutations.append(str((q_pars_pos+POS))+'-'+str((q_pars_pos+POS+int(run_length)-1))+'Del')

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
                                            nt_call_dict_dict[N]['-'] = abund_count
                                        else:
                                            nt_call_dict_dict[N]['-'] += abund_count

                                try:
                                    indel_dict[str(q_pars_pos+POS)+'-'+str(q_pars_pos+POS+run_length-1)+'Del']
                                except:
                                    indel_dict[str(q_pars_pos+POS)+'-'+str(q_pars_pos+POS+run_length-1)+'Del'] = int(abund_count)
                                else:
                                    indel_dict[str(q_pars_pos+POS)+'-'+str(q_pars_pos+POS+run_length-1)+'Del'] += int(abund_count)

                                q_pars_pos = q_pars_pos + run_length

                            elif C == 'M':
                                offset = q_pars_pos-query_pos
                                refPOS = POS+offset
                                for ntPOS in range(query_pos, query_pos+run_length):
                                    if query_seq[ntPOS] != ref[1][refPOS+ntPOS-1]:
                                        if args.AAreport == 1 and args.AAcodonasMNP == 0:
                                            AAinfo = singletCodon(refPOS+ntPOS, query_seq[ntPOS])
                                            mutations.append(str(refPOS+ntPOS)+query_seq[ntPOS]+'('+refprot[AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+')')
                                        else:
                                            mutations.append(str(refPOS+ntPOS)+query_seq[ntPOS])
                                    if args.nt_call == 1:
                                        try:
                                            nt_call_dict_dict[refPOS+ntPOS]
                                        except:
                                            nt_call_dict_dict[refPOS+ntPOS] = {'A' : 0,
                                                                               'T' : 0,
                                                                               'C' : 0,
                                                                               'G' : 0,
                                                                               '-' : 0}
                                            nt_call_dict_dict[refPOS+ntPOS][query_seq[ntPOS]] = abund_count
                                        else:
                                            nt_call_dict_dict[refPOS+ntPOS][query_seq[ntPOS]] += abund_count


                                q_pars_pos = q_pars_pos + run_length
                                query_pos = query_pos + run_length

                            run_length = 0


                        else:
                            run_length = (10 * run_length) + int(C)
                    # END CIGAR PARSE

                    if len(mutations) == 0:
                        if args.wgs == 0:
                            try:
                                seq_species['Reference'] += abund_count
                            except:
                                seq_species['Reference'] = abund_count
                        else:
                            try:
                                seq_species[str(POS)+' Ref '+str(POS+q_pars_pos)] += abund_count
                            except:
                                seq_species[str(POS)+' Ref '+str(POS+q_pars_pos)] = abund_count

                    else:
                        if args.AAreport == 1 and args.AAcodonasMNP == 1:
                            codonchecked = []
                            codon = ''
                            skip = 0
                            MNP = ''
                            for i in range(0, len(mutations)):
                                if 'Del' in mutations[i] or 'insert' in mutations[i]:
                                    codonchecked.append(mutations[i])
                                elif skip > 0:
                                    skip -= 1
                                else:
                                    mut1POS = int(''.join([c for c in mutations[i] if c.isdigit()]))
                                    try:
                                        mutations[i+1]
                                    except:
                                        AAinfo = singletCodon(mut1POS, mutations[i][-1])
                                        codonchecked.append(mutations[i]+'('+refprot[AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+')')
                                    else:

                                        if mut1POS % 3 == 0:
                                            AAinfo = singletCodon(mut1POS, mutations[i][-1])
                                            codonchecked.append(mutations[i]+'('+refprot[AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+')')
                                        else:
                                            mut2POS = int(''.join([c for c in mutations[i+1] if c.isdigit()]))
                                            if mut2POS - mut1POS < 3:
                                                AAPOS = mut1POS // 3
                                                if mut1POS % 3 == 1:
                                                    if mut2POS % 3 == 0:
                                                        codon = mutations[i][-1]+ref[1][mut1POS]+mutations[i+1][-1]
                                                        MNP = mutations[i][-1]+'r'+mutations[i+1][-1]
                                                        skip = 1
                                                    else:
                                                        try:
                                                            mutations[i+2]
                                                        except:
                                                            codon = mutations[i][-1]+mutations[i+1][-1]+ref[1][mut2POS]
                                                            MNP = mutations[i][-1]+mutations[i+1][-1]+'r'
                                                            skip = 1
                                                        else:
                                                            mut3POS = int(''.join([c for c in mutations[i+2] if c.isdigit()]))
                                                            if mut2POS == mut3POS -1:
                                                                codon = mutations[i][-1]+mutations[i+1][-1]+mutations[i+2][-1]
                                                                MNP = mutations[i][-1]+mutations[i+1][-1]+mutations[i+2][-1]
                                                                skip = 2
                                                            else:
                                                                codon = mutations[i][-1]+mutations[i+1][-1]+ref[1][mut2POS]
                                                                MNP = mutations[i][-1]+mutations[i+1][-1]+'r'
                                                                skip = 1
                                                    codonchecked.append(str(mut1POS)+MNP+'('+refprot[AAPOS]+str(AAPOS+1)+AAcall(codon)+')')
                                                elif mut2POS - mut1POS == 1:
                                                    codon = ref[1][mut1POS-2]+mutations[i][-1]+mutations[i+1][-1]
                                                    MNP = mutations[i][-1]+mutations[i+1][-1]
                                                    skip = 1
                                                    codonchecked.append(str(mut1POS)+MNP+'('+refprot[AAPOS]+str(AAPOS+1)+AAcall(codon)+')')
                                                else:
                                                    AAinfo = singletCodon(mut1POS, mutations[i][-1])
                                                    codonchecked.append(mutations[i]+'('+refprot[AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+')')


                                            else:
                                                AAinfo = singletCodon(mut1POS, mutations[i][-1])
                                                codonchecked.append(mutations[i]+'('+refprot[AAinfo[0]-1]+str(AAinfo[0])+AAinfo[1]+')')
                            mutations = codonchecked

                        if args.wgs == 0:
                            try:
                                seq_species[" ".join(mutations)] += abund_count
                            except:
                                seq_species[" ".join(mutations)] = abund_count
                        else:
                            try:
                                seq_species[str(POS)+' '+" ".join(mutations)+' '+str(POS+q_pars_pos)] += abund_count
                            except:
                                seq_species[str(POS)+' '+" ".join(mutations)+' '+str(POS+q_pars_pos)] = abund_count

                    if lastPOS < POS+q_pars_pos:
                        lastPOS = POS+q_pars_pos
                    for i in range(POS, POS+q_pars_pos):
                        try:
                            coverage[i] += abund_count
                        except:
                            coverage[i] = abund_count

    # END SAM LINES
    print(f"End sam parse for {samp}")
    # print(coverage)

    if sam_read_count == 0:
        print(f"No Reads for {samp}")

    else:

        if args.min_abundance1 < 1:
            min_abund = args.min_abundance1
        else:
            min_abund = args.min_abundance1 / sam_read_count

        if args.seq == 1:
            seq_fh = open(samp+'_unique_seqs.tsv', "w")
            seq_fh.write(samp+"("+str(sam_read_count)+")\n")
            seq_fh.write("Unique Sequence\tCount\tAbundance\n")

            sorted_seq = sorted(seq_species, key=seq_species.__getitem__, reverse=True)
            for key in sorted_seq:
                if (seq_species[key] / sam_read_count >= min_abund) and args.wgs == 0:
                    seq_fh.write(f"{key}\t{seq_species[key]}\t{(seq_species[key]/sam_read_count):.3f}\n")
                elif args.wgs == 1:
                    splitseqs = key.split()
                    cov = []
                    for x in range(int(splitseqs[0]), int(splitseqs[-1])):
                        cov.append(coverage[x])
                    min_cov = min(cov)
                    if seq_species[key]/min_cov >= min_abund:
                        seq_fh.write(f"{key}\t{seq_species[key]}\t{(seq_species[key]/min_cov):.3f}\n")

            seq_fh.close()
            # END SEQ OUT
            print(f"End unqiue seq out for {samp}")

        if args.indel == 1 and len(indel_dict) > 0:
            indel_fh = open(samp+'_indels.tsv', "w")
            indel_fh.write(samp+"("+str(sam_read_count)+")\n")
            indel_fh.write("Indel\tCount\tAbundance\n")
            sorted_indels = sorted(indel_dict, key=indel_dict.__getitem__, reverse=True)
            for key in sorted_indels:
                if indel_dict[key] / sam_read_count >= min_abund and args.wgs == 0:
                    indel_fh.write(f"{key}\t{indel_dict[key]}\t{(indel_dict[key]/sam_read_count):.3f}\n")
                elif args.wgs == 1:
                    indelPOS = ''
                    for c in key:
                        if c.isdigit():
                            indelPOS += c
                        else:
                            break
                    indelPOS = int(indelPOS)
                    if indel_dict[key] / coverage[indelPOS] >= min_abund:
                        indel_fh.write(f"{key}\t{indel_dict[key]}\t{(indel_dict[key] / coverage[indelPOS]):.3f}\n")
                    
            indel_fh.close()
            # END INDEL OUT
            print(f"End indel out for {samp}")

        if args.nt_call == 1:
            ntcall_fh = open(samp+'_nt_calls.tsv', "w")
            ntcall_fh.write(samp+"("+str(sam_read_count)+")\n")
            ntcallv_fh = open(samp+'_nt_calls_varonly.tsv', "w")
            ntcallv_fh.write(samp+"("+str(sam_read_count)+")\n")
            sorted_POS = sorted(nt_call_dict_dict)
            if args.AAreport == 1:
                ntcall_fh.write("Position\tref NT\tAA POS\tref AA\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                ntcallv_fh.write("Position\tref NT\tAA POS\tref AA\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tPrimary Seq AA\tsingle nt AA\tSecondary NT\tCounts\tAbundance\tAA\tTertiary NT\tCounts\tAbundance\tAA\n")
                for POS in sorted_POS:
                    try:
                        total = coverage[POS]
                    except:
                        total = 0
                    if total >= (sam_read_count * args.ntabund):
                        AAinfo = singletCodon(POS, ref[1][POS-1])
                        POS_calls = {}
                        for key in nt_call_dict_dict[POS]:
                            POS_calls[key] = nt_call_dict_dict[POS][key]
                        sorted_calls = sorted(POS_calls, key=POS_calls.__getitem__, reverse=True)

                        ntcall_fh.write(str(POS)+"\t"+ref[1][POS-1]+"\t"+str(AAinfo[0])+"\t"+AAinfo[1])
                        ntcall_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                        ntcall_fh.write("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[0]]))
                        ntcall_fh.write(f"\t{(nt_call_dict_dict[POS][sorted_calls[0]]/total):.3f}")
                        if sorted_calls[0] != ref[1][POS-1]:
                            ntcallv_fh.write(str(POS)+"\t"+ref[1][POS-1]+"\t"+str(AAinfo[0])+"\t"+AAinfo[1])
                            ntcallv_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                            ntcallv_fh.write("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[0]]))
                            ntcallv_fh.write(f"\t{(nt_call_dict_dict[POS][sorted_calls[0]]/total):.3f}")

                            mod = (POS)%3

                            if mod == 0:
                                try:
                                    codon = ref[1][POS-3]+ref[1][POS-2]+sorted_calls[0]
                                except:
                                    codon = 'NNN'
                            elif mod == 2:
                                try:
                                    codon = ref[1][POS-2]+sorted_calls[0]+ref[1][POS]
                                except:
                                    codon = 'NNN'
                            elif mod == 1:
                                try:
                                    codon = sorted_calls[0]+ref[1][POS]+ref[1][POS+1]
                                except:
                                    codon = 'NNN'
                            ntcall_fh.write("\t"+AAcall(codon)+"\t"+singletCodon(POS, sorted_calls[0])[1])
                            ntcallv_fh.write("\t"+AAcall(codon)+"\t"+singletCodon(POS, sorted_calls[0])[1])
                            if nt_call_dict_dict[POS][sorted_calls[1]] /total > min_abund:
                                if sorted_calls[1] != ref[1][POS-1]:
                                    ntcallv_fh.write(f"\t{sorted_calls[1]}\t{nt_call_dict_dict[POS][sorted_calls[1]]}\t{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}"+"\t"+singletCodon(POS, sorted_calls[1])[1])
                                ntcall_fh.write(f"\t{sorted_calls[1]}\t{nt_call_dict_dict[POS][sorted_calls[1]]}\t{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}"+"\t"+singletCodon(POS, sorted_calls[1])[1])
                                if nt_call_dict_dict[POS][sorted_calls[2]] /total  > min_abund:
                                    if sorted_calls[2] != ref[1][POS-1]:
                                        ntcallv_fh.write(f"\t{sorted_calls[2]}\t{nt_call_dict_dict[POS][sorted_calls[2]]}\t{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}\t{singletCodon(POS, sorted_calls[2])[1]}")
                                    ntcall_fh.write(f"\t{sorted_calls[2]}\t{nt_call_dict_dict[POS][sorted_calls[2]]}\t{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}\t{singletCodon(POS, sorted_calls[2])[1]}")

                            ntcallv_fh.write("\n")
                        elif nt_call_dict_dict[POS][sorted_calls[1]]  /total > min_abund:
                            ntcallv_fh.write(str(POS)+"\t"+ref[1][POS-1]+"\t"+str(AAinfo[0])+"\t"+AAinfo[1])
                            ntcallv_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                            ntcallv_fh.write("\t"+str(total)+"\t\t")
                            ntcallv_fh.write(f"\t")
                            ntcall_fh.write("\t\t")
                            ntcallv_fh.write("\t\t")
                            ntcall_fh.write(f"\t{sorted_calls[1]}\t{nt_call_dict_dict[POS][sorted_calls[1]]}\t{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}"+"\t"+singletCodon(POS, sorted_calls[1])[1])
                            ntcallv_fh.write(f"\t{sorted_calls[1]}\t{nt_call_dict_dict[POS][sorted_calls[1]]}\t{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}"+"\t"+singletCodon(POS, sorted_calls[1])[1])
                            if nt_call_dict_dict[POS][sorted_calls[2]] /total  > min_abund:
                                ntcall_fh.write(f"\t{sorted_calls[2]}\t{nt_call_dict_dict[POS][sorted_calls[2]]}\t{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}\t{singletCodon(POS, sorted_calls[2])[1]}")
                                if sorted_calls[2] != ref[1][POS-1]:
                                    ntcallv_fh.write(f"\t{sorted_calls[2]}\t{nt_call_dict_dict[POS][sorted_calls[2]]}\t{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}\t{singletCodon(POS, sorted_calls[2])[1]}")
                            ntcallv_fh.write("\n")

                        ntcall_fh.write("\n")
            else:
                ntcall_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")
                ntcallv_fh.write("Position\tref NT\tA\tT\tC\tG\t-\tTotal\tPrimary NT\tCounts\tAbundance\tSecondary NT\tCounts\tAbundance\tTertiary NT\tCounts\tAbundance\n")

                for POS in sorted_POS:
                    try:
                        total = coverage[POS] # sum(nt_call_dict_dict[POS].values())
                    except:
                        total = 0
                    if total >= (sam_read_count * args.ntabund):
                        POS_calls = {}
                        for key in nt_call_dict_dict[POS]:
                            POS_calls[key] = nt_call_dict_dict[POS][key]
                        sorted_calls = sorted(POS_calls, key=POS_calls.__getitem__, reverse=True)

                        ntcall_fh.write(str(POS)+"\t"+ref[1][POS-1])
                        ntcall_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                        ntcall_fh.write("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[0]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[0]]/total):.3f}")

                        if sorted_calls[0] != ref[1][POS-1]:
                            ntcallv_fh.write(str(POS)+"\t"+ref[1][POS-1])
                            ntcallv_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                            ntcallv_fh.write("\t"+str(total)+"\t"+sorted_calls[0]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[0]]))
                            ntcallv_fh.write(f"\t{(nt_call_dict_dict[POS][sorted_calls[0]]/total):.3f}")
                            if nt_call_dict_dict[POS][sorted_calls[1]] /total  > min_abund:
                                ntcall_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                                if sorted_calls[1] != ref[1][POS-1]:
                                    ntcallv_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                                if nt_call_dict_dict[POS][sorted_calls[2]] /total  > min_abund:
                                    ntcall_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                                    if sorted_calls[2] != ref[1][POS-1]:
                                        ntcallv_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                            ntcallv_fh.write("\n")

                        elif nt_call_dict_dict[POS][sorted_calls[1]] /total  > min_abund:
                            ntcallv_fh.write(str(POS)+"\t"+ref[1][POS-1])
                            ntcallv_fh.write("\t"+str(nt_call_dict_dict[POS]['A'])+"\t"+str(nt_call_dict_dict[POS]['T'])+"\t"+str(nt_call_dict_dict[POS]['C'])+"\t"+str(nt_call_dict_dict[POS]['G'])+"\t"+str(nt_call_dict_dict[POS]['-']))
                            ntcallv_fh.write("\t"+str(total)+"\t\t")
                            ntcallv_fh.write(f"\t")
                            ntcall_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                            if sorted_calls[1] != ref[1][POS-1]:
                                ntcallv_fh.write("\t"+sorted_calls[1]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[1]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[1]]/total):.3f}")
                            if nt_call_dict_dict[POS][sorted_calls[2]] /total  > min_abund:
                                ntcall_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                                if sorted_calls[2] != ref[1][POS-1]:
                                    ntcallv_fh.write("\t"+sorted_calls[2]+"\t"+str(nt_call_dict_dict[POS][sorted_calls[2]])+"\t"+f"{(nt_call_dict_dict[POS][sorted_calls[2]]/total):.3f}")
                            ntcallv_fh.write("\n")

                        ntcall_fh.write("\n")

            ntcall_fh.close()
            ntcallv_fh.close()
            # END NT CALL OUT
            print(f"End nt call out for {samp}")

        if args.covar ==1:
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
                if combinations[key] / sam_read_count >= min_abund and args.wgs == 0:
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
                    if coveragepercent >= min_abund:
                        covar_fh.write(f"{key}\t{combinations[key]}\t{coveragepercent:.3f}\n")
                
                    # \t{coveragepercent:.3f}
            covar_fh.close()


            # END COVAR OUT
            print(f"End covar out for {samp}")

def cvdeconv(samp, covardict, seqdict):

    passedseqs = {}
    for seq in seqdict:
        if seq != 'total' and seq != 'singles':
            if len(seq.split(' ')) == 1:
                passedseqs[seq] = max(1, args.beta)
            else:
                splitseq = seq.split(' ')
                abund = 1
                for sing in splitseq:
                    abund = abund * (covardict[sing] / covardict['total'])

                try:
                    covarabund = covardict[seq]/covardict['total']
                except:
                    covarabund = seqdict[seq]/seqdict['total']
                    covardict[seq] = seqdict[seq]

                if covarabund >= abund * args.beta:
                    passedseqs[seq] = covarabund / abund
                if covarabund >= args.autopass:
                    passedseqs[seq] = max(1, args.beta, (covarabund / abund))

    if args.min_abundance1 < 1:
        min_abund = args.min_abundance1 * covardict['total']
    else:
        min_abund = args.min_abundance1

    if args.pass_out == 1:
        fh_pass = open(samp+"_covar_pass.tsv", "w")
        fh_pass.write(f"{samp}({covardict['total']})\tCount\tAbundance\tPass Ratio\n")
        for seq in passedseqs:
            if covardict[seq] >= min_abund:
                fh_pass.write(f"{seq}\t{covardict[seq]}\t{(covardict[seq]/covardict['total']):.3f}\t{passedseqs[seq]:.3f}\n")
        fh_pass.close

    lensortedpassed = sorted(passedseqs, key = lambda key : len(key.split(' ')), reverse=True)
    ratiolensortedpassed = sorted(lensortedpassed, key = lambda key : passedseqs[key], reverse = True)
    sortedsingles = sorted(covardict['singles'], key = covardict['singles'].__getitem__)

    deconved = {}
    for seq in ratiolensortedpassed:
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
    for seq in sorted_deconved:
        if deconved[seq] >= min_abund:
            fh_deconv.write(f"{seq}\t{deconved[seq]}\t{(deconved[seq]/newtotal):.3f}\n")
    fh_deconv.close

    print(f"End covar deconv out for {samp}") # END COVAR DECONV OUT

    return()

def dechim(seqs):
    total = seqs['total']
    del seqs['total']
    sorted_seqs = sorted(seqs, key=seqs.__getitem__)
    chimeras = []
    for seq in sorted_seqs:
        pot_chim = ['Reference']+seq.split()+['Reference']
        chim_halves = []
        for i in range(0, len(pot_chim)-1):
            chim_halves.append([pot_chim[:i+1], pot_chim[i+1:]])
        parent_pairs = []

        for dimera in chim_halves:

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
                for left_par in pot_left:
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
                                        break

                                rt_break_POS = ''
                                for c in rt_break:
                                    if c.isdigit():
                                        rt_break_POS += c
                                    else:
                                        break
                                # if 'Del' in rt_break:
                                    # rt_break_POS = rt_break.split('-')[1].split('D')[0]
                                # elif 'insert' in rt_break:
                                    # rt_break_POS = rt_break.split('-')[0]
                                # elif 'r' in
                                # else:
                                    # rt_break_POS = rt_break.split('(')[0][:-1]

                                if int(left_break_POS) > int(rt_break_POS):
                                    parent_pairs.append([' '.join(left_par[1:-1]), ' '.join(rt_par[1:-1])])

        par_tot_abund = 0
        pair_probs = []
        for parents in parent_pairs:
            pair_prob = (seqs[parents[0]] / total) * (seqs[parents[1]] / total)
            par_tot_abund += pair_prob
            pair_probs.append(pair_prob)

        recomb_count = par_tot_abund * total

        if not seqs[seq] >= recomb_count * args.alpha:
            redist_count = float(seqs[seq])
            chimeras.append(seq)
            seqs[seq] = 0
            if args.redist == 1:
                toadd = {}
                for i in range(0, len(parent_pairs)):
                    # abundx = (( seqs[parent_pairs[i][0]] / total ) * ( seqs[parent_pairs[i][1]] / total ))
                    # print(f"{seq} {seqs[seq]} {redist_count} {parent_pairs[i]} {abundx} {pair_probs[i]} {par_tot_abund}")
                    counts_to_redist = (redist_count * (pair_probs[i]/par_tot_abund))/2
                    seqs[parent_pairs[i][0]] += counts_to_redist
                    seqs[parent_pairs[i][1]] += counts_to_redist



    for chim in chimeras:
        del seqs[chim]

    total = sum(seqs.values())


    # total = sum(seqs.values)
    seqs['total'] = total

    return(seqs)

def chimrm(samp, seqs):

    pre_len = len(seqs)
    inf_loop_shield = 0
    while True:
        dechim(seqs)
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
    if args.min_abundance1 < 1:
        min_abund = args.min_abundance1 * total
    else:
        min_abund = args.min_abundance1
    sorted_seqs = sorted(seqs, key=seqs.__getitem__, reverse=True)
    fh_dechime = open(f"{samp}_a{args.alpha}f{args.foldab}rd{args.redist}_chim_rm.tsv",'w')
    fh_dechime.write(f"{samp}({int(total)})\n")
    fh_dechime.write("Sequences\tCount\tAbundance\n")
    for seq in seqs:
        abund = seqs[seq]/total
        if seqs[seq] >= min_abund:
            fh_dechime.write(f"{seq}\t{round(seqs[seq])}\t{abund:.3f}\n")

    fh_dechime.close
    print(f"End chim_rm out for {samp}") # END CHIM RM DECONV OUT
    return()

if __name__ == '__main__':

    args = arg_parser()

    if args.sams == 1:

        ref = get_ref()
        if ref[1] == '':
            print('Reference not recognized as a Fasta format, skipping SAM parsing')
        else:
            SAMs = []
            try:
                args.Sam_files[0]
            except:
                for file in os.listdir(os.getcwd()):
                    if (file.lower()).endswith('.sam'):
                        SAMs.append(open(file, "r"))
            else:
                for files in args.Sam_files:
                    for file in files:
                        SAMs.append(file)
            refprot = ''
            if args.AAreport == 1:

                for x in range(0, (len(ref[1])-1)//3):
                    AA = AAcall(ref[1][x*3]+ref[1][x*3+1]+ref[1][x*3+2])
                    refprot = refprot + AA
                if (len(ref[1])-1)%3 != 0:
                    refprot = refprot + '?'


            processes = []
            for file in SAMs:
                p = Process(target=SAMparse, args=(file,))
                p.start()
                processes.append(p)
            # # END SAM FILES
            for p in processes:
                p.join()
            print(f"End Sam Parsing Output")

    in_covars = {}
    in_seqs = {}

    if args.chim_rm == 1 or args.deconv == 1:
        for file in os.listdir(os.getcwd()):
            sampline = []
            if file.endswith('_covars.tsv'):
                if args.deconv == 1:
                    try:
                        open(file, 'r')
                    except:
                        print(f"Can't open {file}, skipping")
                    else:
                        covar_fh = open(file, 'r')
                        for line in covar_fh:
                            lineparts = line.strip("\n\r").split("\t")
                            try:
                                lineparts[1]
                            except:
                                sampline = lineparts[0].split("(")
                                in_covars[sampline[0]] = {'total' : int(sampline[1][0:-1]),
                                                          'singles' : {}
                                                          }
                            else:
                                if lineparts[1] != 'Count':
                                    if float(lineparts[2]) >= args.chim_in_abund:
                                        in_covars[sampline[0]][lineparts[0]] = int(lineparts[1])
                                        if len(lineparts[0].split(' ')) == 1:
                                            in_covars[sampline[0]]['singles'][lineparts[0]] = int(lineparts[1])

                        covar_fh.close


            if file.endswith('_seqs.tsv'):
                try:
                    open(file, 'r')
                except:
                    print(f"Can't open {file}, skipping")
                else:
                    seq_fh = open(file, 'r')
                    for line in seq_fh:
                        lineparts = line.strip("\n\r").split("\t")
                        try:
                            lineparts[1]
                        except:
                            sampline = lineparts[0].split("(")
                            in_seqs[sampline[0]] = {'total' : int(sampline[1][0:-1])}
                        else:
                            if lineparts[1] != 'Count':
                                if float(lineparts[2]) >= args.chim_in_abund:
                                    in_seqs[sampline[0]][lineparts[0]] = float(lineparts[1])
                    seq_fh.close

    if args.deconv == 1:
        deconv_procs = []
        for samp in in_covars:
            # cvdeconv(samp, in_covars[samp], in_seqs[samp])
            deconv_p = Process(target=cvdeconv, args=(samp, in_covars[samp], in_seqs[samp],))
            deconv_p.start()
            deconv_procs.append(deconv_p)

        for deconv_p in deconv_procs:
            deconv_p.join()




    if args.chim_rm == 1:
        cr_procs = []
        for samp in in_seqs:
            # chimrm(samp, in_seqs[samp])
            chimrm_p = Process(target=chimrm, args=(samp, in_seqs[samp],))
            chimrm_p.start()
            cr_procs.append(chimrm_p)

        for chimrm_p in cr_procs:
            chimrm_p.join()


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
                                if float(splitline[2]) >= args.min_abundance2:
                                    covar_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_covars[splitline[0]] = 1


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
                                if float(splitline[2]) >= args.min_abundance2:
                                    seq_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_seqs[splitline[0]] = 1

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
                                if float(splitline[2]) >= args.min_abundance2:
                                    deconv_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_deconv[splitline[0]] = 1

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
                                if float(splitline[2]) >= args.min_abundance2:
                                    pass_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_pass[splitline[0]] = 1

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
                                if float(splitline[2]) >= args.min_abundance2:
                                    cr_dict_dict[sample_line][splitline[0]] = [splitline[1], splitline[2]]
                                    all_cr[splitline[0]] = 1

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