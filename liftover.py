"""
This script assigns SNP positions in different editions of the same reference genome.

It was used in the Good Wines project to assign SNP positions from Cabezas et al. (2011) and
Laucou et al. (2018) papers into the reference Pinot Noir PN40024 v4 edition in 2023.
The papers used an older edition of Pinot Noir from 2012, the PN40024 12x V0. To compare the recommended positions
from the two biologists, we needed to find the corresponding positions in our reference and
check their positions with our analysis.

Usage:
    python3 liftOver.py ref1.fasta ref2.fasta snp.txt

The algorithm takes three inputs:
1. ref1.fasta    | a fasta file with the old edition of the reference
2. ref2.fasta    | a fasta file with the new edition of the reference
3. snp.txt       | a txt file with all the SNPs positions from old edition to be aligned into new edition
					and their corresponding chromosome

The algorithm outputs:
1. results.txt   | a txt file with chromosome ID, old SNP position, new SNP position
"""

import os
import sys
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pandas as pd

def read_snp_file(snp_file):
	"""Reads the SNP file and returns a DataFrame."""
	return pd.read_csv(snp_file, header=None, names=['chr', 'snp_pos_ref1'])

def split_fasta(fa_file):
	"""	Split the fasta file into multiple fasta"""
	with open(fa_file,'r') as handle:
		for seq in SeqIO.parse(handle,'fasta'):
			with open(f"ref2_{seq.id}.fasta","w") as sm_fa:
				sm_fa.write(f">{seq.id},ref2\n")
				sm_fa.write(f"{seq.seq}")

	return 0

def make_fasta_snp(v1_fasta, chr_id, snp_list):
	"""Creates fasta files for each SNP. 
	Takes the SNP positions from snp.txt and creates a small fasta file
	with 200 bases and the SNP placed in the middle i.e at 100th position.
	The fasta id will be > chrN, snp_SNPOS,ref1
	"""
	for snp in snp_list:
		print(f"chromosome: {chr_id} ---> SNP pos: {snp}")
		with open(v1_fasta, 'r') as v1_file:
			queries = SeqIO.parse(v1_file, 'fasta')
			for q in queries:
				if q.id == chr_id:
					break

			sequence = q.seq[snp - 100:snp + 100].upper()

			output_file = f'ref1_{chr_id}_snp_{snp}.fasta'
			with open(output_file, 'w') as handle:
				handle.write(f">{q.id},snp_{snp},ref1\n")
				handle.write(f"{sequence}\n")

			print(f"Fasta file {chr_id} - {snp} written")

def match_snp_ref1_to_ref2(chr_id, ref1_snp_pos):
	"""Matches SNP positions from ref1 to ref2.
	Blasts the small fasta files from make_fasta_snp with their corresponding chromosome.
	The first alignment 100% is supposed to be the correct one aligned.
	"""
	query_file = f"ref1_{chr_id}_snp_{ref1_snp_pos}.fasta"
	subject_file = f"ref2_{chr_id}.fasta"
	cmd = f"blastn -query {query_file} -subject {subject_file} -outfmt 5 -out blast_output_{chr_id}_{ref1_snp_pos}.xml"

	print(f"BLAST query: {query_file} --> subject: {subject_file}")
	subprocess.run(cmd, shell=True)

	blast_output = f'blast_output_{chr_id}_{ref1_snp_pos}.xml'
	ref2_snp_pos = -1

	with open(blast_output) as blast_file:
		blast_records = NCBIXML.parse(blast_file)
		for blast_record in blast_records:
			print(f'Query : {blast_record.query}')
			if blast_record.alignments:
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						subject_start = hsp.sbjct_start
						percentage = (hsp.identities / hsp.align_length)*100
						print(f"Subject: {alignment.title}\nStart Position: {subject_start}\nIdentity % : {percentage:.2f}")
						ref2_snp_pos = subject_start + 100
						print(f"{hsp.query[50:150]}")
						print(f"{hsp.match[50:150]}")
						print(f"{hsp.sbjct[50:150]}")
			else:
				print("No alignments found.")
	return ref2_snp_pos

def pipeline(ref1_fasta, ref2_fasta, snp_file, output_file):
	"""Main pipeline function to process SNPs."""
	split_fasta(ref2_fasta)
	df_snp_file = read_snp_file(snp_file)

	results = []

	for chr_id in df_snp_file['chr'].unique():
		snp_list = df_snp_file[df_snp_file['chr'] == chr_id]['snp_pos_ref1'].tolist()
		print(f"Processing {chr_id}")
		make_fasta_snp(ref1_fasta, chr_id, snp_list)

		for snp_pos in snp_list:
			snp_to_ref2 = match_snp_ref1_to_ref2(chr_id, snp_pos)
			results.append((chr_id, snp_pos, snp_to_ref2))
			print(f"{chr_id} --> {snp_pos} --> {snp_to_ref2}")

	with open(output_file, 'w') as f:
		for result in results:
			f.write(f"{result[0]}\t{result[1]}\t{result[2]}\n")

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("Usage: python3 liftOver.py ref1.fasta ref2.fasta snp.txt")
		sys.exit(1)

	ref1_fasta = sys.argv[1]
	ref2_fasta = sys.argv[2]
	snp_file = sys.argv[3]
	output_file = "results.txt"

	pipeline(ref1_fasta, ref2_fasta, snp_file, output_file)
	print(f"Results written to {output_file}")

