#!/bin/bash

#previous analyses I have done:
#cmd = hnM_TAM_CLIP_clusters_cleaned.bed hg19 homer_TAM_ryan_clusters_exons_background_default_size -rna -len 4,6,8 -basic -bg /Users/Sam/bin/FAST-iCLIP/docs/hg19_transcriptome_collapse_exon.bed

#What input to use for homer peak analysis? 

#Huelga et al., 2012 paper (findMotifs.pl with parameters –len 5 –homer1 –chopify –norevopp –rna –fasta)

#I will use the individual RT stops with 30 basepairs flanking each RT-stop
#For background, I will use introns bound by hnRNPM
#For size, use size given

#DMSO - error coming up

findMotifsGenome.pl ../../hnM_DMSO_clipper/hnM_DMSO_RT_stops_clipper_sig_extension_with_id hg19 dmso_homer_results -size given -rna -len 6 -basic -bg ../../hnM_DMSO_clipper/hnM_DMSO_RT_stops_clipper_sig_introns_intersect 