# Identifying mutations from FASTA-format sequence alignments

Development environment: Python v3.8.5, pandas v1.2.4.



## Demonstration

```bash
python aln2mut.py -i test/demo.aln -o demo -d test -r reference -l
```

Length of the reference sequence: 21 aa.



## Output

- \[output prefix\]\_var.tsv: a VCF-like tab-delimited file listing alterations

- \[output prefix\]\_mat.tsv: a matrix of alterations (sample-by-variant positions)

- \[output prefix\]\_lst.tsv (optional): a list of alteration in a conventional format (e.g. W25N, 36-38del) for amino acid alterations. Users may want to skip the generation of this file when they are processing nucleotide alignments and create their own alteration lists from the VCF-like outputs.