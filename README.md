# liftover
This script assigns SNP positions in different editions of the same reference genome.

It was used in the Good Wines project to assign SNP positions from Cabezas et al. (2011) and
Laucou et al. (2018) papers into the reference Pinot Noir PN40024 v4 edition in 2023.
The papers used an older edition of Pinot Noir from 2012, the PN40024 12x V0. To compare the recommended positions
from the two biologists, we needed to find the corresponding positions in our reference and
check their positions with our analysis.

Usage:
    python3 liftOver.py rf1.fasta rf2.fasta snp.txt

The algorithm takes three inputs:
1. rf1.fasta    | a fasta file with the old edition of the reference
2. rf2.fasta    | a fasta file with the new edition of the reference
3. snp.txt       | a txt file with all the SNPs positions from old edition to be aligned into new edition
					and their corresponding chromosome

The algorithm outputs:
1. results.txt   | a txt file with chromosome ID, old SNP position, new SNP position
