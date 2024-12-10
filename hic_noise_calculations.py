#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 15:17:41 2023

@author: whirling-in-rags
"""

import pandas as pd
import numpy as np
import os

os.chdir("/home/whirling-in-rags/projects/phase_hic/28_hiczin_2.0/hiczin_output/normalized_outputs")

g86_df = pd.read_csv("/home/whirling-in-rags/projects/phase_hic/28_hiczin_2.0/hiczin_output/normalized_outputs/g86_HiCZin_normalized.csv")
o91_df = pd.read_csv("/home/whirling-in-rags/projects/phase_hic/28_hiczin_2.0/hiczin_output/normalized_outputs/o91_HiCZin_normalized.csv")
w76_df = pd.read_csv("/home/whirling-in-rags/projects/phase_hic/28_hiczin_2.0/hiczin_output/normalized_outputs/w76_HiCZin_normalized.csv")

#calculation from hwang et al 2023
#calculating raw noise and relaxed noise 
#raw noise:
    # num inter-mMAG contacts/num intra-mMAG contacts
#relaxed noise:
    #num inter-order contacts/num intra-mMAG contacts
    
##ref code
#come back to this, bascailly remove all rows where mMAG is linking to another mMAG 
#removes mMAG x mMAG contacts
#raw_contacts_cleaned = (raw_contacts.loc[~((raw_contacts[6] == 'mMAG') & (raw_contacts[9]=='mMAG'))])
#same as above, also removes all vMAG x vMAG and pMAG x pMAG contacts
#raw_contacts_cleaner = (raw_contacts_cleaned.loc[~((raw_contacts_cleaned[4] == 'vMAG') & raw_contacts_cleaned[6]=='vMAG')])
#raw_contacts_bleached = (raw_contacts_cleaner.loc[~((raw_contacts_cleaner[4] == 'pMAG') & raw_contacts_cleaner[6]=='pMAG')])
 
   
raw_g86_intra_mMAG_contacts = g86_df.loc[~((g86_df['source_genus'] == g86_df['target_genus']))]
raw_g86_inter_mMAG_contacts = g86_df.loc[~((g86_df['source_genus'] != g86_df['target_genus']))]

relaxed_g86_intra_mMAG_contacts = g86_df.loc[~((g86_df['source_phylotype'] == g86_df['target_phylotype']))]
relaxed_g86_inter_mMAG_contacts = g86_df.loc[~((g86_df['source_phylotype'] != g86_df['target_phylotype']))]

raw_g86_noise_signal_ratio = (len(raw_g86_inter_mMAG_contacts)/len(raw_g86_intra_mMAG_contacts))
relaxed_g86_noise_signal_ratio = (len(relaxed_g86_inter_mMAG_contacts)/len(relaxed_g86_intra_mMAG_contacts))

raw_o91_intra_mMAG_contacts = o91_df.loc[~((o91_df['source_genus'] == o91_df['target_genus']))]
raw_o91_inter_mMAG_contacts = o91_df.loc[~((o91_df['source_genus'] != o91_df['target_genus']))]

relaxed_o91_intra_mMAG_contacts = o91_df.loc[~((o91_df['source_phylotype'] == o91_df['target_phylotype']))]
relaxed_o91_inter_mMAG_contacts = o91_df.loc[~((o91_df['source_phylotype'] != o91_df['target_phylotype']))]

raw_o91_noise_signal_ratio = (len(raw_o91_inter_mMAG_contacts)/len(raw_o91_intra_mMAG_contacts))
relaxed_o91_noise_signal_ratio = (len(relaxed_o91_inter_mMAG_contacts)/len(relaxed_o91_intra_mMAG_contacts))

raw_w76_intra_mMAG_contacts = w76_df.loc[~((w76_df['source_genus'] == w76_df['target_genus']))]
raw_w76_inter_mMAG_contacts = w76_df.loc[~((w76_df['source_genus'] != w76_df['target_genus']))]

relaxed_w76_intra_mMAG_contacts = w76_df.loc[~((w76_df['source_phylotype'] == w76_df['target_phylotype']))]
relaxed_w76_inter_mMAG_contacts = w76_df.loc[~((w76_df['source_phylotype'] != w76_df['target_phylotype']))]

raw_w76_noise_signal_ratio = (len(raw_w76_inter_mMAG_contacts)/len(raw_w76_intra_mMAG_contacts))
relaxed_w76_noise_signal_ratio = (len(relaxed_w76_inter_mMAG_contacts)/len(relaxed_w76_intra_mMAG_contacts))