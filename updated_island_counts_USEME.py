#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 17:25:41 2024

@author: robinch
"""

import pandas as pd
import os


names = pd.read_csv('/home/whirling-in-rags/projects/phase_hic/28_hiczin_2.0/hiczin_output/taxa/taxa_list_inputs_FIXED_230922_USEME.csv')



#get list of all MGI from viral genomes for each genome from csv
def load_and_concatenate_csvs(directory_path):
    csv_files = [f for f in os.listdir(directory_path) if f.endswith('.csv')]
    
    #create empty df list
    df_list = []
    
    #iterate through csv files, load into df_list
    for i, file in enumerate(csv_files):
        file_path = os.path.join(directory_path,file)
        if i ==0:
            df = pd.read_csv(file_path)
        else:
            df = pd.read_csv(file_path,header=0) #skips headers, they duplicate
        df_list.append(df)
        
    #concatenate into single df
    combined_df = pd.concat(df_list,ignore_index=True)
    
    return combined_df

directory_path = '/home/whirling-in-rags/projects/phase_hic/38_viral_coverage_MGI/contig_groupby'

combined_df = load_and_concatenate_csvs(directory_path)

#checking presesence or absence of vMAG in metagenome based on condition:
    #if number of MGI in group_contig,metagenome pair > 50% of number of windows in vMAG
    #then it is absent from metagenome
    
#if window mean < thresh_hold mean, assigned 1, otherwise 0
combined_df['low_cov'] = combined_df.apply(lambda row: 1 if row['window_mean'] < row['threshold_mean'] else 0, axis=1)

#group by [group_contig, metagenome] and calculate max window and sum of low_cov

grouped_df = combined_df.groupby(['group_contig','metagenome']).agg(
    max_window_id=('window_id','max'),
    sum_low_cov=('low_cov','sum'),
    avg_threshold_mean = ('threshold_mean','mean'),
    avg_window_mean=('window_mean','mean')
).reset_index()

grouped_df['pres_abs'] = grouped_df.apply(
    lambda row: 0 if row['sum_low_cov'] > (row['max_window_id'] * 0.25) else 1, axis=1
)


#now going to create a more descriptive table from grouped_df

#first, group by 'group_contig' and calculate max_window_id, highest avg_window_mean, and total pres_abs
total_df = grouped_df.groupby('group_contig').agg(
    max_window_id=('max_window_id','max'),
    highest_avg_window_mean=('avg_window_mean','max'),
    total_pres_abs=('pres_abs','sum')
).reset_index()

metagenomes = grouped_df['metagenome'].unique()

for metagenome in metagenomes:
    # Create a new column for each metagenome and initialize it with 0
    total_df[metagenome] = 0

    # Iterate through each group_contig
    for i, row in total_df.iterrows():
        group_contig = row['group_contig']
        
        # Find the corresponding 'pres_abs' value for this group_contig and metagenome
        match = grouped_df[(grouped_df['group_contig'] == group_contig) & (grouped_df['metagenome'] == metagenome)]
        
        # If there's a match, assign the pres_abs value, otherwise keep it as 0
        if not match.empty:
            total_df.at[i, metagenome] = match['pres_abs'].values[0]


total_df['votu'] = None
for i, row in total_df.iterrows():
    group_contig = row['group_contig']
    match = names[names.iloc[:,1] == group_contig]
    if not match.empty:
        total_df.at[i,'votu'] = match.iloc[0,-2]
        
#total_df.to_csv('WIP_ImportSheet_vMAG_metagenomes.csv')

