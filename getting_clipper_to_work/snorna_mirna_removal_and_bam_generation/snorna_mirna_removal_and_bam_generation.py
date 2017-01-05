import os, cmath, math, sys, glob, subprocess, re, argparse, shutil, datetime, csv, commands
import numpy as np
import pandas as pd
from collections import defaultdict
from operator import itemgetter
import matplotlib as mpl
import matplotlib.pyplot as plt
from optparse import OptionParser
mpl.rcParams['savefig.dpi'] = 2 * mpl.rcParams['savefig.dpi']
mpl.rcParams['path.simplify'] = True
csv.register_dialect("textdialect",delimiter='\t')

snoRNAmasker = 'sno_coordinates_hg19_formatted.bed'
miRNAmasker = 'miR_sort_clean.bed'
genome = '-s hg19'
genomeFile = 'hg19.sizes'

def filter_snoRNAs(negAndPosMerged, snoRNAmasker, miRNAmasker):
	# Usage: Filter snoRNA and miRNAs from protein coding reads.
	# Input: .bed file with protein coding reads.
	# Output: snoRNA and miR filtered .bed file.
	program='intersectBed'
	proteinWithoutmiRNAs = negAndPosMerged.replace('.bed','_snoRNAremoved_miRNAremoved.bed')

	cmd1 = "bedtools intersect -a {} -b {} -wa -v -s | sort -k1,1 -k2,2n".format(negAndPosMerged, snoRNAmasker)
	cmd2_1 = "bedtools intersect -a - -b {} -wa -v -s -sorted".format(miRNAmasker)
	cmd2_2 = "awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {print $1,$2,$3,$4 \"_\" NR,$5,$6}'"
	cmd2_3 = proteinWithoutmiRNAs
	cmd = cmd1 + ' | ' + cmd2_1 + ' | ' + cmd2_2 + ' > ' + cmd2_3
	os.system(cmd)
	
	return proteinWithoutmiRNAs

def makeBamAndRunClipper(RTclusterfile,genome,genomeFile):
	# Usage: Process the mergedRT file and pass through CLIPper FDR script.
	# Input: Merged RT file.
	# Output: CLIPper input (.bed) file and output file.
	
	bamfile_sorted = RTclusterfile.replace('.bed','.srt')  
	cmd1 = "bedToBam -i {} -g {} | samtools sort - {}".format(RTclusterfile, genomeFile, bamfile_sorted)
	os.system(cmd1)
	
	bamfile_sorted += ".bam"
	mapStats=bamfile_sorted.replace('.srt.bam','.mapStats.txt') 
	cmd2 = "samtools flagstat {} > {}".format(bamfile_sorted, mapStats)
	cmd3 = "samtools index {}".format(bamfile_sorted)
	os.system(cmd2)
	os.system(cmd3)

	CLIPPERout_dup = RTclusterfile.replace('.bed','_CLIP_clusters_dupl') 
	cmd4 = "clipper --bam {} {} --outfile={} > /dev/null 2>&1".format(bamfile_sorted, genome, CLIPPERout_dup)
	os.system(cmd4)

	# added by BD 4/12/15 to merge adjacent clip clusters and remove duplicates
	CLIPPERout = CLIPPERout_dup.replace('_CLIP_clusters_dupl','_CLIP_clusters') 
	with open(CLIPPERout_dup,'r') as ifile, open(CLIPPERout,'w') as ofile:
		reader = csv.reader(ifile, 'textdialect')
		writer = csv.writer(ofile, 'textdialect')
		currRow = ['chr1',0,0,0,0,'+']
		for row in reader:
			currStart = int(currRow[1])
			currEnd = int(currRow[2])
			newStart = int(row[1])
			newEnd = int(row[2])
			if currStart==newStart and currEnd==newEnd: continue #duplicates
			if math.fabs(newStart-currEnd) <= 15 and currRow[5]==row[5]: #overlap and same strand
				if int(currRow[1]) != 0: #not the first one
					currRow[2]=newEnd #merge the two adjacent clusters
			else: #not overlap
				if int(currRow[1]) != 0:
					writer.writerow(currRow)
				currRow = row #cycle continues
		writer.writerow(currRow) #fencepost
	
	return CLIPPERout
	

