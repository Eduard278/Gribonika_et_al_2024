# 2404-iR-0005 Additional B-cell analysis

Prepared on 2024-04-16.

## Preprocessing

Sequences from the _pep.csv were realigned with IgBlast for CDR/FR annotation. This was to ensure a common
starting point (CDR2) for fair SHM calculation.

## cdr2_fr4_airr

Derived from IgBlast output. Includes the most important columns in AIRR format:
V, J, C, CDR3 junction aa, CDR3 junction nt, CDR1 -> FR4 peptide, CDR1 -> FR4 nucleotide,
read/UMI count, and SHM calculations. These files were used for all downstream analysis.

## cdr3_pep_hd1

This graph shows a network of CDR3 peptide sequences, where each sequence is 
treated as a node, and edges are defined by having a hamming distance of 1
between a pair of nodes. Node size is based off relative expression within 
the sample. Node color is based off isotype, as shown in the legend (the 
majority of CDR3 peptide sequences only have 1 isotype, but in cases where 
there are multiple, the highest frequency isotype is chosen).

## individual treemaps

Each rectangle represents a unique clone (CDR3 peptide + V + J combination), with its size
proportional to its expression level in the sample. In these locally generated treemaps,
the same CDR3 always has the same color, regardless of sample.

## data points

Definitions for each index are included in `repertoire_metric_definitions.txt`.

WARNING: Please be wary of the the class-switching indexes like CSR0/CSR1 which
are based off of low-frequency isotypes. For example, if IgG3 and IgG2 are both
rare in the samples, we wouldn't recommend using the index `IGHG2-IGHG3_CSR1`.

## isotype upset plots

UpSet plot showing combinations of isotypes on individual CDR3 peptide sequences.
