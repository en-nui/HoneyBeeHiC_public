#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:20:07 2024

@author: whirling-in-rags
"""

import pandas as pd


#Small function to count consecutive integer sequences in windows to determine length of each 100bp window
#It also calculates the window mean from the means of 100 bp windows
#Its input is a 4 column dataframe with columns ['group_contig','window_id','window_id_mean','threshold_mean']. 

g86 = pd.read_csv('/home/whirling-in-rags/projects/phase_hic/38_viral_coverage_MGI/filtered_window/g86_filtered_window.csv')
o91 = pd.read_csv('/home/whirling-in-rags/projects/phase_hic/38_viral_coverage_MGI/filtered_window/o91_filtered_window.csv')
w76 = pd.read_csv('/home/whirling-in-rags/projects/phase_hic/38_viral_coverage_MGI/filtered_window/w76_filtered_window.csv')


def consequence_sequence_count(df):
    #list to store length of windows
    window_lengths = []
    #list to store means of windows
    window_means = []
    #list to hold current series of integers
    current_window = []
    #list to hold means from current series of integers
    current_mean = 0
    #function to loop through rows. Checks to see if current row is empty OR if window_id value is equal to current window_id value  + 1
    #adds current window to current_window list, adds current window mean to current_mean list
    #otherwise it's a new window and appends previous current_window to window_lengths
    
    for i in range(len(df)):
        if not current_window or df['window_id'].iloc[i] == current_window[-1] + 1:
            current_window.append(df['window_id'].iloc[i])
            current_mean += df['window_mean'].iloc[i]
        else:
            window_lengths.append(len(current_window))
            window_means.append(current_mean / len(current_window))
            current_window = [df['window_id'].iloc[i]]
            current_mean = df['window_mean'].iloc[i]
    if current_window:
        window_lengths.append(len(current_window))
        window_means.append(current_mean / len(current_window))
    return pd.DataFrame({'window_length': window_lengths,'window_mean':window_means})

g86_result = consequence_sequence_count(g86.copy())
o91_result = consequence_sequence_count(o91.copy())
w76_result = consequence_sequence_count(w76.copy())
all_result_tmp = [g86_result,o91_result,w76_result]
all_result = pd.concat(all_result_tmp)
g86_counts = g86_result['window_length'].value_counts()
o91_counts = o91_result['window_length'].value_counts()
w76_counts = w76_result['window_length'].value_counts()
g86_counts_summed = pd.DataFrame({'n_windows': g86_counts.index, 'count': g86_counts.values})
o91_counts_summed = pd.DataFrame({'n_windows': o91_counts.index, 'count': o91_counts.values})
w76_counts_summed = pd.DataFrame({'n_windows': w76_counts.index, 'count': w76_counts.values})

unique_list_tmp_counts = all_result['window_length'].value_counts()
all_result_uniques_summed = pd.DataFrame({'n_windows': unique_list_tmp_counts.index, 'count': unique_list_tmp_counts.values})

g86_counts_summed.to_csv('g86_MGI_counts_summed.csv')
o91_counts_summed.to_csv('o91_MGI_counts_summed.csv')
w76_counts_summed.to_csv('w76_MGI_counts_summed.csv')
all_result_uniques_summed.to_csv('concat_MGI_counts_summed.csv')