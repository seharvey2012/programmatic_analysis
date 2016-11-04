#!/bin/bash

#run dCLIP using the iCLIP parameter.
#Note - dCLIP is designed to use the unprocessed RT-stops, but after removing the header and such.
#This means taking the mappedToGenome_sorted.bam files for each sample and concatenating them
#Note that dCLIP doesn't take into account differential expression...

#Concatenate the original files and place them in the dCLIP analysis folder

#samtools merge -@ 8 -O SAM dmso_r1_r2_sam_merged.sam /media/sam/Data1/hnRNPM_iCLIP/new_ryan_iclip_download_2016_03_12/HS_hnRNPM_DMSO_T4_clipper/rawdata_and_stats/HS_hnRMPM_DMSO_4nt_R1_trimmed_mappedToGenome_sorted.bam /media/sam/Data1/hnRNPM_iCLIP/new_ryan_iclip_download_2016_03_12/HS_hnRNPM_DMSO_T4_clipper/rawdata_and_stats/HS_hnRMPM_DMSO_4nt_R2_trimmed_mappedToGenome_sorted.bam

#samtools merge -@ 8 -O SAM tam_r1_r2_sam_merged.sam /media/sam/Data1/hnRNPM_iCLIP/new_ryan_iclip_download_2016_03_12/HS_hnRNPM_TAM_T4_clipper/rawdata_and_stats/HS_hnRMPM_TAM_4nt_R1_trimmed_mappedToGenome_sorted.bam /media/sam/Data1/hnRNPM_iCLIP/new_ryan_iclip_download_2016_03_12/HS_hnRNPM_TAM_T4_clipper/rawdata_and_stats/HS_hnRMPM_TAM_4nt_R2_trimmed_mappedToGenome_sorted.bam

#Run dCLIP in -iCLIP mode. Note the important parameters
#-filter - there is no control so they recommend > 10 reads for a >40 Million mapped read dataset 
#DMSO = 63.2 million TAM = 47.9 million

dCLIP.pl -iCLIP 7 -f1 dmso_r1_r2_sam_merged.sam -f2 tam_r1_r2_sam_merged.sam -temp temp -dir dclip_output

#I can also run dCLIP on a set of already preprocessed reads
#-process filter = upload a user-processed tag intensity count file. What does this file need to look like? The readme details it, it would take a little work to make this