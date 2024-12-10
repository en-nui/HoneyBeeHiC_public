#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 15:48:23 2024

@author: robinch
"""

import pandas as pd

#This purpose of this code is to import a bed file via bedtools bga and calculate the mean coverage for each contig
#Then, it will identify regions (consecutive bases >= 100bp) where the mean coverage is 0.25*mean_coverage of contig

g86_bed = pd.read_csv('/home/whirling-in-rags/projects/phase_hic/38_viral_coverage_MGI/bed/depth/g86_viruses_depth.tsv',sep="\t")
o91_bed = pd.read_csv('/home/whirling-in-rags/projects/phase_hic/38_viral_coverage_MGI/bed/depth/o91_viruses_depth.tsv',sep="\t")
w76_bed = pd.read_csv('/home/whirling-in-rags/projects/phase_hic/38_viral_coverage_MGI/bed/depth/w76_viruses_depth.tsv',sep="\t")


#fluff below
g86_bed = g86_bed.set_axis(["contig", "position","coverage"], axis="columns")
o91_bed = o91_bed.set_axis(["contig", "position","coverage"], axis="columns")
w76_bed = w76_bed.set_axis(["contig", "position","coverage"], axis="columns")

#below code groups by 'contig' column and gets mean of column 'coverage' for each group
g86_bed_groups = g86_bed.groupby('contig')['coverage'].mean()
g86_bed_groups_df = g86_bed_groups.to_frame().reset_index()
g86_bed_groups_df['cutoff_25'] = g86_bed_groups_df['coverage']*0.25
#Now, I want to identify groups contigs that are >5x and contigs that are between 2x - 5x
g86_bed_zero = g86_bed_groups_df[(g86_bed_groups_df['coverage']<2)]
g86_bed_minimum = g86_bed_groups_df[(g86_bed_groups_df['coverage'] > 2) & (g86_bed_groups_df['coverage'] < 5)]
g86_bed_good = g86_bed_groups_df[g86_bed_groups_df['coverage']>5]

#duplicate check
duplicate_check = pd.concat([g86_bed_zero, g86_bed_minimum, g86_bed_good])
duplicates = duplicate_check[duplicate_check.duplicated()]



def process_contig_groups(df, window_size=100):
    filtered_windows_df = pd.DataFrame()
    unique_contig_window_ids = []
    all_window_info = []  # List to store information for all windows
    window_id = 0
    for contig, group_df in df.groupby('contig'):
        group_mean = group_df['coverage'].mean()
        threshold_mean = group_mean * 0.25
        
        num_windows = (len(group_df) // window_size) + 1  # Calculate number of windows
        
        for i in range(num_windows):
            window_start = i * window_size
            window_end = min((i + 1) * window_size, len(group_df))
            window_df = group_df.iloc[window_start:window_end]

            # Add group_contig, window_id, group_mean, threshold_mean, and window_mean columns
            window_df['group_contig'] = contig
            window_df['window_id'] = window_id
            window_df['group_mean'] = group_mean
            window_df['threshold_mean'] = threshold_mean
            window_df['window_mean'] = window_df['coverage'].mean()  # Add window_mean column
            
            #filter
            if (window_df['group_mean'] > 5).any() and (window_df['window_mean'] < threshold_mean).any():
                filtered_windows_df = pd.concat([filtered_windows_df, window_df])

                # Add values to unique_contig_window_ids, including window_mean and threshold_mean
                unique_contig_window_ids.append((contig, window_id, window_df['window_mean'].iloc[0], threshold_mean))

            # Add information for all windows to all_window_info list
            all_window_info.append((contig, window_id, window_df['window_mean'].iloc[0], threshold_mean))

            window_id += 1
        window_id=0
    # Create unique_df DataFrame with additional columns
    unique_df = pd.DataFrame(unique_contig_window_ids, columns=['group_contig', 'window_id', 'window_mean', 'threshold_mean'])

    # Create all_windows_df DataFrame with the same structure as unique_df
    all_windows_df = pd.DataFrame(all_window_info, columns=['group_contig', 'window_id', 'window_mean', 'threshold_mean'])

    return filtered_windows_df, unique_df, all_windows_df  # Return all three DataFrames


g86_filtered_windows_df, g86_unique_df, g86_all_windows_df = process_contig_groups(g86_bed)

o91_filtered_windows_df, o91_unique_df, o91_all_windows_df = process_contig_groups(o91_bed)

w76_filtered_windows_df, w76_unique_df, w76_all_windows_df = process_contig_groups(w76_bed)



