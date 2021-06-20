#!/usr/bin/env python
"""
This script identifies and reports mutations from a FASTA-format sequence alignment.

Example command:
    python aln2mut.py -i samples.aln -o samples -d var -r reference -l

Output files:
    [output prefix]_var.tsv: a VCF-like tab-delimited file listing alterations
    [output prefix]_mat.tsv: a matrix of alterations (sample-by-variant positions)
    [output prefix]_lst.tsv (optional): a list of alteration in a conventional format (e.g. W25N, 36-38del)

Dependencies: Python 3, pandas

Copyright (C) 2021 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 19 June 2021; the latest update: 20 June 2021
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
    parser.add_argument('-l', '--list', dest = 'list', action = 'store_true', help = "Create a list of alterations in a conventional format (e.g., W25N)")
    return parser.parse_args()


def main():
    args = parse_argument()
    aln = import_alignment(args.input)  # aln = {sequence ID : sequence}
    ref_name = args.ref
    ref, ref_gap_free = validate_params(ref_name, aln, args.dir)  # The script exists when this validation fails.
    print("Info: Length of the reference sequence: %i" % (len(ref) - ref.count('-')), file = sys.stderr)  # The sequence length does not count the size of gaps.
    aln = {key : value for key, value in aln.items() if key != ref_name}  # Take a subset of the dictionary aln
    vcf = aln2vcf(ref, aln)  # Create a VCF-like table (pandas data frame) of six columns: Sample, Pos, Ref, Alt, Type (of mutation), Aux_pos
    vcf.to_csv(os.path.join(args.dir, args.output + '_var.tsv'), index = False, sep = '\t')  # This table can be read into Python using pandas.read_csv('XXXX.vcf', sep = '\t').
    mat, vcf_samples = vcf2mat(vcf, ref_gap_free, ref_name, list(aln.keys()))  # Convert a VCF data frame into a matrix of alterations (sample x variant positions) that can be aligned to a phylogenetic tree
    mat.to_csv(os.path.join(args.dir, args.output + '_mat.tsv'), index = False, sep = '\t')
    if args.list:
        lst = vcf2lst(vcf, vcf_samples)
        lst.to_csv(os.path.join(args.dir, args.output + '_lst.tsv'), index = False, sep = '\t')
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


def validate_params(ref_name, aln, outdir):
    """
    Check whether the reference sequence is in the alignment file and whether the output directory is accessible
    """
    msg = None
    if len(aln) == 1:
        msg = "Error: Input alignment must consist of at least two sequences."
    elif not (ref_name in aln.keys()):
        msg = "Error: Reference sequence " + ref_name + " is not found in the alignment file."
    elif not os.path.exists(outdir):
        msg = "Error: output directory " + outdir + " does not exist."
    else:
        ref_seq = aln[ref_name]  # The reference sequence (may contain '-' when insertions are present in sample sequences).
        ref_gap_free = ref_seq.replace('-', '')  # Otherwise, the coordinates in function get_ref_row do not match those in the reference sequence.
        if ref_gap_free == '':
            msg = "Error: the reference sequence cannot be a gap (namely, consisting of only dash characters)."
    if msg != None:
        print(msg, file = sys.stderr)
        sys.exit(1)
    return ref_seq, ref_gap_free


def aln2vcf(ref, aln):
    """
    Mutation calling: to create a VCF-like table of five columns: Sample, Pos, Ref, Alt, Type (of mutation)
    Code for alteration types: S (substitution), I (insertion to the reference), D (deletion from the reference)
    This function implements the core algorithm of this script.
    """
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
        for i in range(0, len(ref)):  # Python character indexes across sequences in the alignment, including '-' characters. In the alignment file, every sequence (including '-') must have the same length.
            r = ref[i]  # Reference base/AA
            s = seq[i]  # Sample base/AA
            if r == s:
                if r != '-':  # s = r = '-' when there is a larger insertion in another sequence overlapping the current sample and reference sequences.
                    p += 1
                    if ins:  # The pointer reaches the end of an insertion.
                        samples.append(sam)
                        if p == 1:
                            coords.append('^1')  # The insertion happens before the first character of the reference. Do not need to write '0^1' here.
                        else:
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
                        coords.append('^1' if p == 1 else '%i^%i' % (ins_up, p))
                        aux_coords.append(ins_up + 0.5)
                        refs.append('-')
                        alts.append(ins_seq)
                        types.append('I')
                        ins = False
                        ins_seq = ''
                    else:  # Substitution or deletion
                        coords.append(str(p))
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


def vcf2mat(vcf, ref_gap_free, ref_name, all_samples):
    """ Convert a VCF data frame into a matrix of alterations (sample x variant positions) """
    var_pos = vcf[['Aux_pos', 'Pos']].drop_duplicates()  # Select two columns from the data frame VCF and drop duplicated rows (There is no need to use pos.Aux_pos.unique() later as Aux_pos and Pos are linked)
    var_pos = var_pos.sort_values(by = ['Aux_pos'], ascending = True)  # Sort the data frame by auxiliary positions
    var_pos = var_pos.reset_index()
    n = len(var_pos.index)  # Row count of var_pos
    mat = pandas.DataFrame(columns = ['Sample'] + var_pos.Pos.tolist())
    mat = mat.append(pandas.Series([ref_name] + get_ref_row(ref_gap_free = ref_gap_free, coords = var_pos.Aux_pos.tolist()),\
          index = mat.columns), ignore_index = True)
    vcf_samples = vcf.Sample.unique()  # Returns an array object
    for sam in all_samples:  # all_samples does not include the reference sequence.
        new_row = [sam]  # One row per sample
        if sam in vcf_samples:
            vcf_sam = vcf.loc[vcf['Sample'] == sam]  # Select rows corresponding to the current sample
            vcf_sam_auxpos = vcf_sam.Aux_pos.tolist()
            for _, row in var_pos.iterrows():  # Iterate through positions of variants (namely, non-ID columns of the output matrix)
                p = row['Aux_pos']  # A numeric value.
                if p in vcf_sam_auxpos:  # The current sample has an alteration at position p.
                    new_row.append(vcf_sam.loc[vcf_sam['Aux_pos'] == p, 'Alt'].iloc[0])  # A single row will be extracted from vcf_sam. Then get the first (also the only one) value from the column Alt.
                else:
                    new_row.append('.')  # The current character is identical to that in the reference.
        else:
            print("Info: sequence " + sam + " is identical to the reference sequence.", file = sys.stderr)
            new_row += ['.'] * n
        mat = mat.append(pandas.Series(new_row, index = mat.columns), ignore_index = True)  # Append a row to mat
    return mat, vcf_samples


def get_ref_row(ref_gap_free, coords):
    """
    A subordinate function of vcf2mat for making a list of reference characters corresponding to alteration positions
    """
    ref_chars = list()
    for p in coords:
        if int(p) == p:  # Insertion: int(p) < p.
            try:
                i = int(p) - 1
                ref_chars += [ref_gap_free[i]]
            except IndexError:
                print("Runtime error: index " + str(i) + " exceeds the range of the gap-free reference sequence.", file = sys.stderr)
                sys.exit(1)
        else:
            ref_chars += ['-']
    return(ref_chars)


def vcf2lst(vcf, vcf_samples):
    """ Convert the VCF into a list of alterations """
    lst = pandas.DataFrame(columns = ['Sample', 'Alteration'])
    for sam in vcf_samples:
        vcf_sam = vcf.loc[vcf['Sample'] == sam]
        alt = list()
        for _, row in vcf_sam.iterrows():
            t = row['Type']
            if t == 'S':
                alt.append(row['Ref'] + row['Pos'] + row['Alt'])  # E.g. V80F
            elif t == 'I':
                alt.append(row['Pos'].replace('^', 'ins' + row['Alt']))  # E.g. 10insATRQ11 (The prefix 'ins' seems redundant here, but it can be used as a keyword by users for quickly identifying insertions)
            else:  # t == 'D'
                alt.append(row['Ref'] + row['Pos'] + 'del')  # E.g. R80del. Users may want to manually merge consecutive deletions into a single one.
        lst = lst.append(pandas.Series([sam, ','.join(alt)], index = lst.columns), ignore_index = True)  # Sample name'\t'A comma-delimited list of alterations
    return lst


if __name__ == '__main__':
    main()