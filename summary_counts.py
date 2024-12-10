#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 17:25:41 2024

@author: robinch
"""
import pandas as pd

df = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/38_viral_coverage_MGI/annotation_popgen_statistics/all_metagenome_annotation_coord_summary_statistics_presabs.csv')
df['total_metagenome_sum'] = df[['g86','o91','w76']].sum(axis=1)
sum_input = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/38_viral_coverage_MGI/annotation_popgen_statistics/summary_figre/nucDiv_theta_summaries_input.csv')
df_2_meta = df[df['total_metagenome_sum'] == 2]
df_3_meta = df[df['total_metagenome_sum']==3]


#this is to get cumulative number of islands across and mean # islands within metagenome
g86_windows = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/38_viral_coverage_MGI/filtered_window/g86_filtered_window.csv')
o91_windows = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/38_viral_coverage_MGI/filtered_window/o91_filtered_window.csv')
w76_windows = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/38_viral_coverage_MGI/filtered_window/w76_filtered_window.csv')
g86_windows['metagenome']="g86"
o91_windows['metagenome']="o91"
w76_windows['metagenome']="w76"

joined_windows = pd.concat([g86_windows,o91_windows,w76_windows])
grouped_df = joined_windows.groupby('group_contig')['window_id'].count().reset_index()
grouped_df_mean = joined_windows.groupby('group_contig')['window_mean'].mean().reset_index()

#coverage data
g86_cov = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/37_annotation_coverage/viruses/g86_viruses_depth.tsv',sep="\t")
o91_cov = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/37_annotation_coverage/viruses/o91_viruses_depth.tsv',sep="\t")
w76_cov = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/37_annotation_coverage/viruses/w76_viruses_depth.tsv',sep="\t")
g86_cov.columns=['contig','pos','cov']
g86_cov['metagenome']="g86"
o91_cov.columns=['contig','pos','cov']
o91_cov['metagenome']="o91"
w76_cov.columns=['contig','pos','cov']
w76_cov['metagenome']="w76"



g86_average_cov = g86_cov.groupby('contig')['cov'].mean().reset_index()
g86_average_cov['metagenome']="g86"
o91_average_cov = o91_cov.groupby('contig')['cov'].mean().reset_index()
o91_average_cov['metagenome']="o91"
w76_average_cov = w76_cov.groupby('contig')['cov'].mean().reset_index()
w76_average_cov['metagenome']="w76"



joined_cov = pd.concat([g86_cov,o91_cov,w76_cov])

joined_cov_avg = joined_cov.groupby('contig')['cov'].mean().reset_index()
joined_cov_max = pd.concat([g86_average_cov,o91_average_cov,w76_average_cov])

filtered_cov_max = (
    joined_cov_max.groupby(['contig'])['cov']
    .max()
    .reset_index()
)




annotation_groupby = df.groupby('contig').count().reset_index()

joined_cov_filtered = joined_cov_avg[joined_cov_avg["contig"].isin(annotation_groupby["contig"])]
filtered_cov_max2 = filtered_cov_max[filtered_cov_max["contig"].isin(annotation_groupby["contig"])]

g86_cov_filtered = g86_average_cov[g86_average_cov["contig"].isin(annotation_groupby["contig"])]
o91_cov_filtered = o91_average_cov[o91_average_cov["contig"].isin(annotation_groupby["contig"])]
w76_cov_filtered = w76_average_cov[w76_average_cov["contig"].isin(annotation_groupby["contig"])]
window_df_filtered = grouped_df[grouped_df["group_contig"].isin(annotation_groupby["contig"])]
window_df_mean_filtered = grouped_df_mean[grouped_df_mean["group_contig"].isin(annotation_groupby["contig"])]



'''
joined_cov_filtered.to_csv('joined_cov_avg.csv')
filtered_cov_max2.to_csv('filtered_cov_max.csv')
window_df_filtered.to_csv('window_df_filtered.csv')
window_df_mean_filtered.to_csv('window_df_mean.csv')
'''


sum_stats = sum_input.groupby('contig')[['pNpS_gene_reference']].mean().reset_index()
sum_stats.to_csv('sum_pNpS.csv')