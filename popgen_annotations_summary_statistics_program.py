#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 15:40:30 2024

@author: robinch
"""

import pandas as pd
import math

def select_matching_rows(df1, df2):
    """
    Selects rows from df2 that overlap with regions defined in df1,
    calculates various diversity statistics, and returns a combined DataFrame.

    Args:
        df1: A DataFrame containing columns 'scaffold', 'start_position', 'end_position', and 'gene_position'.
        df2: A DataFrame containing columns 'contig', 'pos', and other relevant data.

    Returns:
        A DataFrame containing merged results with calculated diversity statistics.
    """

    output_list = []
    for scaffold in df1['scaffold'].unique():
        matching_df1 = df1[df1['scaffold'] == scaffold]
        matching_df2 = pd.DataFrame()

        for i, row in matching_df1.iterrows():
            start_position = row['start_position']
            end_position = row['end_position']
            gene_position = row['gene_position']

            temp_df = df2[(df2['contig'] == scaffold) & (df2['pos'] >= start_position) & (df2['pos'] <= end_position)]
            temp_df['start_position'] = start_position
            temp_df['end_position'] = end_position
            temp_df['gene_position'] = gene_position

            matching_df2 = pd.concat([matching_df2, temp_df], ignore_index=True)

        output_list.append(matching_df2.copy())

    gene_position_list = []
    for df in output_list:
        for gene_position in df['gene_position'].unique():
            gene_df = df[df['gene_position'] == gene_position]
            gene_position_list.append(gene_df.copy())

    for df in gene_position_list:
        df['freq_sum'] = df.apply(
            lambda row: ((((row['A'] / row['coverage'])**2 +
                            (row['T'] / row['coverage'])**2 +
                            (row['C'] / row['coverage'])**2 +
                            (row['G'] / row['coverage'])**2))),
            axis=1
        )

        df['normalization'] = df.apply(
            lambda row: (row['coverage']/(row['coverage'] - 1)),
            axis=1
        )

        df['heterozygosity'] = df.apply(
            lambda row: ((1- row['freq_sum']) * row['normalization']),
            axis=1
        )

    unique_contig_list = []
    for df in gene_position_list:
        contig = df['contig'].iloc[0]
        start_position = df['start_position'].iloc[0]
        end_position = df['end_position'].iloc[0]
        summed_heterozygosity = df['heterozygosity'].sum()
        average_coverage = df['coverage'].mean()
        snpLength = len(df)

        unique_df = pd.DataFrame({
            'contig': contig,
            'gene_position': df['gene_position'].unique(),
            'start_position': start_position,
            'end_position': end_position,
            'nucleotide_diversity': summed_heterozygosity,
            'segSite_mean_coverage': average_coverage,
            'snpLength': snpLength
        })

        unique_contig_list.append(unique_df)

    merged_dataframes = []
    unique_contigs = set(df['contig'].iloc[0] for df in unique_contig_list)

    for contig in unique_contigs:
        contig_dataframes = [df for df in unique_contig_list if df['contig'].iloc[0] == contig]
        merged_df = pd.concat(contig_dataframes, ignore_index=True) 
        merged_dataframes.append(merged_df)

    merged_dataframes2 = []    
    for dataframe in merged_dataframes:
        for index, row in dataframe.iterrows():
            n_segSites = row['snpLength']
            length = (row['end_position']-row['start_position'])
            div = row['nucleotide_diversity']
            div_by_length = div / length
            theta = row['nucleotide_diversity']
            n = row['segSite_mean_coverage']
    
            varPi = (((n + 1) / (3 * (n - 1))) * theta) + ((theta**2) * (2 * ((n**2) + n + 3)) / (9 * n * (n - 1)))
    
            Wtheta_denom = sum(1/i for i in range(1,int(n)))
            Wtheta = n_segSites / Wtheta_denom
            Wtheta_by_length = Wtheta/length
            tajima_d_var_a1 = sum(1/i for i in range(1,int(n)))
            tajima_d_var_a2 = sum(1/i**2 for i in range(1,int(n)))
            
            tajima_d_var_b1 = (n + 1) / (3 * (n - 1))
            tajima_d_var_b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
            tajima_d_var_c1 = tajima_d_var_b1 - 1 / tajima_d_var_a1
            tajima_d_var_c2 = tajima_d_var_b2 - (n + 2) / (tajima_d_var_a1 * n) + tajima_d_var_a2 / tajima_d_var_a1**2
            tajima_d_var_e1 = tajima_d_var_c1/tajima_d_var_a1
            tajima_d_var_e2 = tajima_d_var_c2 / (tajima_d_var_a1**2 + tajima_d_var_a2)
            
            tajima_d_var = math.sqrt((tajima_d_var_e1 * n_segSites) + tajima_d_var_e2 * n_segSites * (n_segSites-1))
            
            tajima_D = (div - Wtheta)/tajima_d_var
            
            # create new dataframe with calculated statistics
            new_df = pd.DataFrame({
                'contig': row['contig'],
                'gene_position': row['gene_position'],
                'start_position': row['start_position'],
                'end_position': row['end_position'],
                'gene_length': length,
                'nucleotide_diversity': row['nucleotide_diversity'],
                'segSite_mean_coverage': row['segSite_mean_coverage'],
                'snpLength': row['snpLength'],
                'nucDiv_by_snp': div_by_length,
                'varPi': varPi,
                'Wtheta': Wtheta,
                'Wtheta_by_snp': Wtheta_by_length,
                'varD': tajima_d_var,
                'tajimaD': tajima_D
            }, index=[index])
            
            merged_dataframes2.append(new_df)
            
            # final operations
    merged_dataframes2 = pd.concat(merged_dataframes2)
    return merged_dataframes2


df1 = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/32_vmag_correlations/annotations.tsv',sep='\t')

df2 = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/38_viral_coverage_MGI/annotation_popgen_statistics/g86_annotation_extracted_snp.csv')
df3 = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/38_viral_coverage_MGI/annotation_popgen_statistics/o91_annotation_extracted_snp.csv')
df4 = pd.read_csv('/home/robinch/projects/phase_hic/HoneyBeeHiC/38_viral_coverage_MGI/annotation_popgen_statistics/w76_annotation_extracted_snp.csv')

result_g86 = select_matching_rows(df1,df2)
result_o91 = select_matching_rows(df1,df3)
result_w76 = select_matching_rows(df1,df4)

result_g86.to_csv('g86_annotation_coord_summary_statistics.csv')
result_o91.to_csv('o91_annotation_coord_summary_statistics.csv')
result_w76.to_csv('w76_annotation_coord_summary_statistics.csv')


