#!/home/sam/anaconda/bin/python

import pandas as pd
import numpy as np
import os
import csv
import argparse

#add arguments
parser = argparse.ArgumentParser(description='Process rMATS output using an FDR cutoff, PSI cutoff, and read cutoff')
parser.add_argument('-i', metavar='INPUT', help="rMATS skipped exon output file pathname", required=True)
parser.add_argument('-fdr', metavar='FDR_cutoff', type=float, help="Splice events must have < FDR. Default is 0.05.", default=0.05)
parser.add_argument('-psi', metavar='PSI_cutoff', type=float, help="Differential events must differ by at >= PSI. Default is 0.2.", default=0.2)
parser.add_argument('-read', metavar='read_cutoff', type=int, help="For a given splicing event, each sample must have at least this many reads when adding included + excluded reads. Default is 20.", default=20)

#define parameters
args = parser.parse_args()
filename = args.i
FDR_cutoff = args.fdr
PSI_cutoff = args.psi
read_cutoff = args.read

def parse_rMATS_SE(filename,FDR_cutoff,PSI_cutoff,read_cutoff):
    '''This function will parse rMATS output in filename to output only events that meet an
    FDR cutoff and PSI cutoff. The output include both positive and negative delta PSI as well
    as each separated
    Input: filename, FDR_cutoff which is a float 0-1, and PSI cutoff which is a float > 1, and read_cutoff which is and int'''
    rMATS_file = open(filename).readlines()
    rMATS_list = []
    for line in rMATS_file:
        line = line.strip().split('\t')
        fixed_line = []
        for x in line:
            x = x.strip('"')
            fixed_line.append(x)
        rMATS_list.append(fixed_line)
    #add a new column to the list withe the splice ID
    #Note - I would like to figure out how to do this with an existing dataframe
    for item in rMATS_list[0:1]:
        item.insert(0,"splice_id")
    for item in rMATS_list[1:]:
        splice_id = str(item[3]+':'+item[7]+'-'+item[8]+':'+item[5]+'-'+item[6]+':'+item[9]+'-'+item[10]+':'+item[4]+":"+item[2])
        item.insert(0,splice_id)
    #add new columns for sample_1_read_number and sample_2_read_number
    for item in rMATS_list[0:1]:
        item.insert(15,"SAMPLE_1_AVERAGE_READ")
    for item in rMATS_list[1:]:
        sample_1_ic_list = item[13].strip().split(',')
        sample_1_ic_numbers = [int(x) for x in sample_1_ic_list]
        sample_1_sc_list = item[14].strip().split(',')
        sample_1_sc_numbers = [int(x) for x in sample_1_sc_list]
        sample_1_average_read = np.mean(sample_1_ic_numbers) + np.mean(sample_1_sc_numbers)
        item.insert(15,sample_1_average_read)
    for item in rMATS_list[0:1]:
        item.insert(18,"SAMPLE_2_AVERAGE_READ")
    for item in rMATS_list[1:]:
        sample_2_ic_list = item[16].strip().split(',')
        sample_2_ic_numbers = [int(x) for x in sample_2_ic_list]
        sample_2_sc_list = item[17].strip().split(',')
        sample_2_sc_numbers = [int(x) for x in sample_2_sc_list]
        sample_2_average_read = np.mean(sample_2_ic_numbers) + np.mean(sample_2_sc_numbers)
        item.insert(18,sample_2_average_read)    
    #add new column for enhance or silence effect
    for item in rMATS_list[0:1]:
        item.append("splicing_factor_effect")
    for item in rMATS_list[1:]:
        if float(item[-1]) < 0:
            item.append("silence")
        else:
            item.append("enhance")
    #add new columns to have the mean of the PSI values since it is a comma separate list
    for item in rMATS_list[0:1]:
        item.insert(24,"inc_level_1_mean")
    for item in rMATS_list[1:]:
        mean = item[23].strip().split(',')
        mean_no_na = filter(lambda x: x != "NA",mean)
        mean_numbers = [float(i) for i in mean_no_na]
        mean_value = np.mean(mean_numbers)
        item.insert(24,mean_value)
    for item in rMATS_list[0:1]:
        item.insert(26,"inc_level_2_mean")
    for item in rMATS_list[1:]:
        mean = item[25].strip().split(',')
        mean_no_na = filter(lambda x: x != "NA",mean)
        mean_numbers = [float(i) for i in mean_no_na]
        mean_value = np.mean(mean_numbers)
        item.insert(26,mean_value)
    #write this file as a csv for easier import into pandas. Then import it
    with open('tempfile','wb') as temp:
        tempwriter = csv.writer(temp,delimiter='\t')
        for item in rMATS_list:
            tempwriter.writerow(item)
    with open('tempfile') as temp:
        rMATS_df = pd.read_csv(temp,sep='\t')
    os.remove('tempfile')
    #remove samples less than read cutoff. This means that each sample must have inc + exc >= read cutoff
    rMATS_df = rMATS_df[rMATS_df.SAMPLE_1_AVERAGE_READ >= read_cutoff]
    rMATS_df = rMATS_df[rMATS_df.SAMPLE_2_AVERAGE_READ >= read_cutoff]
    #remove FDR > FDR_cutoff
    rMATS_df = rMATS_df[rMATS_df.FDR < FDR_cutoff]
    #remove those beneath delta PSI cutoff
    rMATS_df = rMATS_df[rMATS_df.IncLevelDifference.abs() >= PSI_cutoff]
    #now write this to a text file and an Excel
    with open(filename.split('.')[0]+'_FDR_'+str(FDR_cutoff)+'_dPSI_'+str(PSI_cutoff)+'_read_cutoff_'+str(read_cutoff)+'.txt','w') as txt_file:
        rMATS_df.to_csv(txt_file,sep='\t',index=False)
    #writer = pd.ExcelWriter(filename.split('.')[0]+'_FDR_'+str(FDR_cutoff)+'_dPSI_'+str(PSI_cutoff)+'_read_cutoff_'+str(read_cutoff)+'.xlsx')
    #rMATS_df.to_excel(writer,index=False)
    #writer.save()
    print('Done')

def main():
    parse_rMATS_SE(filename,FDR_cutoff,PSI_cutoff,read_cutoff)

if __name__ == "__main__":
    main()
