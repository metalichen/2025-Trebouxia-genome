import sys
from Bio import SeqIO
import gzip
from os.path import join
import pandas as pd
import re

inputvcf = sys.argv[1]
inputref = sys.argv[2]
outvcf = sys.argv[3]
minlen = int(sys.argv[4])

bigscf = []
totalseqs, filtercount = [0, 0]


#Get list of samples: only one
samples = ['GTX0536']

# Read fasta file
for seq_record in SeqIO.parse(inputref, "fasta"):
	totalseqs += 1
	if len(seq_record) >= minlen:
		bigscf.append(seq_record.id)
		filtercount += 1

print("Number of sequences in input fasta file: " + str(totalseqs))
print("Number of sequences larger than " + str(minlen) + " bps: " + str(filtercount))

vcfopen = gzip.open(inputvcf, 'rt') # 't' is text mode, to interpret the tabs and new lines
outputvcf = open(outvcf, "w")

for line in vcfopen:
	if "Chrom" in line:
		outputvcf.write("##fileformat=VCFv4.X\n")

		# Add explanations
		outputvcf.write(f"##FORMAT=<ID=Cons,Type=String,Description='Consensus in IUPAC confidence'>\n")
		outputvcf.write(f"##FORMAT=<ID=Cov,Number=R,Type=Integer,Description='Coverage'>\n")
		outputvcf.write(f"##FORMAT=<ID=Reads1,Number=R,Type=Integer,Description='Allele 1'>\n")
		outputvcf.write(f"##FORMAT=<ID=Reads2,Number=R,Type=Integer,Description='Allele 2'>\n")
		outputvcf.write(f"##FORMAT=<ID=Freq,Type=String,Description='Frequency of Allele 2'>\n")
		outputvcf.write(f"##FORMAT=<ID=P-value,Number=R,Type=Float,Description='P-value'>\n")

		# Get the names of the samples
		individuals = '\t'.join(samples)

		# Print a new header
		outputvcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{individuals}\n")

	else:
		line = re.sub(r"[\s]", '\t', line) # Replace white spaces from VarScan with tabs (wtf, why are they there anyway?)
		columns = line.rstrip("\n").split('\t')
		CHROM, POS, REF, ALT, QUAL, FILTER, SamplesRef, SamplesHet, SamplesHom, SamplesNC = columns[0:10]

		FORMAT = "Cons:Cov:Reads1:Reads2:Freq:P-value"
		indivs = '\t'.join(columns[10:])
		indivs = indivs.rstrip("\t")

		if CHROM in bigscf:
			newline = f"{CHROM}\t{POS}\t.\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t.\t{FORMAT}\t{indivs}\n"
			outputvcf.write(newline)

vcfopen.close()
outputvcf.close()
