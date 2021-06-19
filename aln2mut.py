#!/usr/bin/env python
"""
This script identifies and reports mutations from a FASTA-format sequence alignment.

Example command:
    python aln2mut.py -i samples.aln -o samples -d var -r reference

Dependencies: Python 3, pandas

Copyright (C) 2021 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 19 June 2021; the latest update: 19 June 2021
"""

import os
import sys
import pandas
from argparse import ArgumentParser


def parse_argument():
    parser = ArgumentParser(description = "Identify mutations from a FASTA-format sequence alignment")
    parser.add_argument('-i', '--input', dest = 'input', type = str, required = True, help = "Input alignment")
    parser.add_argument('-o', '--output', dest = 'output', type = str, required = True, help = "Prefix for output files")
    parser.add_argument('-d', '--dir', dest = 'dir', type = str, required = False, default = '.', help = "Output directory")
    parser.add_argument('-r', '--ref', dest = 'ref', type = str, required = True, help = "Name of the reference sequence in the alignment")
    return parser.parse_args()


def main():
    args = parse_argument()
    aln = import_alignment(args.input)  # aln = {sequence ID : sequence}
    ref_name = args.ref
    validate_params(ref_name, aln, args.dir)  # The script exists when this validation fails.
    ref = aln[ref_name]  # The reference sequence (may contain '-' when insertions are present in sample sequences)
    print("Length of the reference sequence: %i" % (len(ref) - ref.count('-')), file = sys.stderr)  # The sequence length does not count the size of gaps.
    aln = {key : value for key, value in aln.items() if key != ref_name}  # Take a subset of the dictionary aln
    vcf = aln2vcf(ref, aln)  # Create a VCF-like table (pandas data frame) of six columns: Sample, Pos, Ref, Alt, Type (of mutation), Aux_pos
    vcf.to_csv(os.path.join(args.dir, args.output + ".vcf"), index = False, sep = '\t')
    mat = vcf2mat(vcf = vcf, ref = ref, all_samples = list(aln.keys()))  # Convert a VCF data frame into a matrix of alterations (sample x variant positions) that can be aligned to a phylogenetic tree
    mat.to_csv(os.path.join(args.dir, args.output + ".tsv"), index = False, sep = '\t')
    return


def import_alignment(fasta):
    """ Import sequences from the input FASTA file in which sequences may contain '-' characters. """
    aln = dict()
    with open(fasta, 'r') as f:
        lines = f.read().splitlines()
    s = ''  # Sequence cache
    i = None  # ID of the sequence in the cache
    for line in lines:
        if line.startswith('>'):  # A new sequence record is encountered
            if i != None:
                aln[i] = s
                s = ''  # Reset the sequence cache
            i = line.split(' ')[0]  # Drop the sequence annotation from the header line
            i = i[1 : ]  # Drop '>' from the ID
        else:
            s += line
    aln[i] = s  # Save the last record
    return aln


def validate_params(ref, aln, outdir):
    """
    Check whether the reference sequence is in the alignment file and whether the output directory is accessible
    """
    msg = None
    if len(aln) == 1:
        msg = "Error: Input alignment must consist of at least two sequences."
    elif not (ref in aln.keys()):
        msg = "Error: Reference sequence " + ref + " is not found in the alignment file."
    elif not os.path.exists(outdir):
        msg = "Error: output directory " + outdir + " does not exist."
    if msg != None:
        print(msg, file = sys.stderr)
        sys.exit(1)
    return


def aln2vcf(ref, aln):
    """
    Mutation calling: to create a VCF-like table of five columns: Sample, Pos, Ref, Alt, Type (of mutation)
    Code for alteration types: S (substitution), I (insertion to the reference), D (deletion from the reference)
    This function implements the core algorithm of this script.
    """
    n = len(ref)  # In the alignment file, very sequence (including '-') must have the same length.
    samples = list()
    coords = list()  # Positions in the reference sequence
    aux_coords = list()  # Auxiliary positions (real numbers) for sorting the output VCF file based on positions
    refs = list()  # Bases/amino acids (AAs) in the reference sequence
    alts = list()  # Alternative bases/AAs in the sample sequence
    types = list()  # Types of alterations
    for sam, seq in aln.items():  # Of note, a sample does not appear in the VCF if it is identical to the reference.
        p = 0  # A pointer for the current character in the reference sequence, excluding '-' characters. (real position - 1)
        ins = False  # A flag for whether the current position is in an insertion ('-' characters in the reference sequence)
        ins_up = 0  # Immediately upstream position of an insertion; Variables ins_up and p mark the flanking positions of the current insertion.
        ins_seq = ''  # The inserted sequence
        for i in range(0, n):  # Python character indexes across sequences in the alignment, including '-' characters.
            r = ref[i]  # Reference base/AA
            s = seq[i]  # Sample base/AA
            if r == s:
                if s != '-':  # s = r = '-' when there is a larger insertion in another sequence overlapping the current sample and reference sequences.
                    p += 1
                    if ins:  # The pointer reaches the end of an insertion.
                        samples.append(sam)
                        coords.append('%i^%i' % (ins_up, p))  # E.g., '25^36' means an insertion between positions 25 and 26. ins_up equals zero when the insertion happens before the start codon of the reference sequence.
                        aux_coords.append(ins_up + 0.5)  # The extra 0.5 marks the insertion site between two consecutive positions.
                        refs.append('-')
                        alts.append(ins_seq)
                        types.append('I')
                        ins = False
                        ins_seq = ''
            else:  # Record an alteration. Note that in a multisequence alignment, r and s may both be '-'.
                if r != '-':  # The alteration is a deletion, a substitution, or the end of an insertion.
                    samples.append(sam)
                    p += 1
                    if ins:  # Now the insertion region ends
                        coords.append('%i^%i' % (ins_up, p))
                        aux_coords.append(ins_up + 0.5)
                        refs.append('-')
                        alts.append(ins_seq)
                        types.append('I')
                        ins = False
                        ins_seq = ''
                    else:  # Substitution or deletion
                        coords.append(p)
                        aux_coords.append(p)
                        refs.append(r)
                        alts.append(s)
                        types.append('D' if s == '-' else 'S')
                else:  # The alteration is an insertion to the reference. The variable p does not increase in this situation.
                    if ins:  # The current position remains in an insertion
                        ins_seq += s
                    else:  # At the start of an insertion: start to record the current insertion.
                        ins = True
                        ins_up = p  # Record where the insertion starts
                        ins_seq = s  # Start to record the current insertion region
        if ins:  # Finish the insertion that happens at the end of the reference sequence
            samples.append(sam)
            coords.append(str(ins_up) + '^')
            aux_coords.append(ins_up + 0.5)
            refs.append('-')
            alts.append(ins_seq)
            types.append('I')
            ins = False
            ins_seq = ''
    return pandas.DataFrame({'Sample' : samples, 'Pos' : coords, 'Ref' : refs, 'Alt' : alts, 'Type' : types, 'Aux_pos' : aux_coords})


def vcf2mat(vcf, ref, all_samples):
    """ Convert a VCF data frame into a matrix of alterations (sample x variant positions) """
    mat = None
    return mat

if __name__ == '__main__':
    main()