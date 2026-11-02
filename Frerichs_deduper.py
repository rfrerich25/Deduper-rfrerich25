#!/usr/bin/env python

#Import modules
import argparse
import gzip
import re

#define argparse
def get_args():
	parser = argparse.ArgumentParser(description="Program to remove PCR duplicates from a sorted sam file. Also reports number of headers, duplicates, unique reads, wrong UMIs, and reads per chromosome. It is recomended to use samtools sort prior to running. If duplicates are found the first read encountered is output. This program cannot handle hardclipping.")
	parser.add_argument("-f", "--file", help="The file path to sorted sam file", required=True, type=str)
	parser.add_argument("-o", "--outfile", help="The file path to output sam file", required=True, type=str)
	parser.add_argument("-u", "--umi", help="The file containing a list of UMIs", required=True, type=str)
	return parser.parse_args()

args = get_args()

#Functions
def convert_strand(flag:int) -> bool:
	'''Converts the bitwise flag into a bool where true means it is the reverse strand'''
	if ((flag & 16) == 16): #see if 16th bit is flipped
		reverse = True
		return reverse
	else:
		reverse = False
		return reverse

def convert_position(Pos:int, reverse:bool, CIGAR:str) -> int:
	'''Converts the Pos to 5' pos and accounts for softclipping'''
	if reverse == True:
		length = 0 
		matches = re.findall(r"([0-9]+)[MDN=X]", CIGAR) #finds all digits that come before MDN=X
		for match in matches:
			length += int(match) #sums all digits to find length
		endclip = re.search(r"([0-9]+)S$", CIGAR) #Check for end softclipping
		if endclip == None:
			Position = int(Pos) + int(length) - 1
		else:
			Position = int(Pos) + int(length) + int(endclip.group(1)) - 1
	else:
		startclip = re.search(r"^([0-9]+)S", CIGAR) #Check for beginning softclipping
		if startclip == None:
			Position = Pos
		else:
			Position = Pos - int(startclip.group(1))
	return Position

# Initialize dictionaries/variables
umi_dict = {}
seen_dict = {}
stats_dict = {}
stats_dict["header"] = 0
stats_dict["wrong_umi"] = 0
stats_dict["duplicate"] = 0
stats_dict["unique"] = 0
chrom_reads = {}
current_chrom = 0

#Create umi_dict from file
with open(args.umi) as um:
	for line in um:
		line = line.strip("\n")
		umi_dict[line] = 0

#Parse file to deduplicate
with open(args.file, "r") as fh:
	with open(args.outfile, "w") as fw:
		for full_line in fh:
			if full_line.startswith("@"): #write headers
				fw.write(full_line)
				stats_dict["header"] += 1
			else:
				line = full_line.strip().split("\t")
				#Check UMIs
				almost_umi = line[0]
				umi = almost_umi.split(":")[-1]
				if umi not in umi_dict:
					stats_dict["wrong_umi"] += 1
					continue
				#Assign variables
				chrom = line[2]
				strand = convert_strand(int(line[1]))
				pos = convert_position(int(line[3]), strand, line[5])
				#Compare to previous reads
				if chrom != current_chrom:
					seen_dict.clear()
				if (pos, strand, umi) in seen_dict:
					stats_dict["duplicate"] += 1
					continue
				else:
					fw.write(full_line)
					stats_dict["unique"] += 1
					if chrom in chrom_reads:
						chrom_reads[chrom] += 1
					else:
						chrom_reads[chrom] = 1
					current_chrom = chrom
					seen_dict[pos, strand, umi] = 1

print(stats_dict)
print(chrom_reads)