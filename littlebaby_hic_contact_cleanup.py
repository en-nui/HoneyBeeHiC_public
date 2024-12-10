import pandas as pd
import numpy as np
import os
#import scipy
#import pickletools
#import pickle
#import rpy2
#import networkx as nx
os.chdir("/home/robinch/projects/phase_hic/HoneyBeeHiC/28_hiczin_2.0/hiczin_output/normalized_outputs")

#raw_contacts = pd.read_csv("g86_HiCZin_normalized.csv",header=None)

#come back to this, bascailly remove all rows where mMAG is linking to another mMAG 
#removes mMAG x mMAG contacts
#raw_contacts_cleaned = (raw_contacts.loc[~((raw_contacts[6] == 'mMAG') & (raw_contacts[9]=='mMAG'))])
#same as above, also removes all vMAG x vMAG and pMAG x pMAG contacts
#raw_contacts_cleaner = (raw_contacts_cleaned.loc[~((raw_contacts_cleaned[4] == 'vMAG') & raw_contacts_cleaned[6]=='vMAG')])
#raw_contacts_bleached = (raw_contacts_cleaner.loc[~((raw_contacts_cleaner[4] == 'pMAG') & raw_contacts_cleaner[6]=='pMAG')])

#load clean normalized contacts
clean_normalized_contacts = pd.read_csv('o91_raw_contacts_cleaned_updated_pMAG.csv')




clean_contacts_derep = clean_normalized_contacts.groupby(['source_genus','target_genus'],as_index=False)['value'].sum()
clean_contacts_derep2 =  clean_contacts_derep[clean_contacts_derep['source_genus'] != clean_contacts_derep['target_genus']]

clean_contacts_counts = clean_normalized_contacts[['source_genus','target_genus']].value_counts()
clean_counts_counts_df = clean_contacts_counts.to_frame().reset_index()
clean_res = pd.merge(clean_contacts_derep2,clean_counts_counts_df)



#raw_contacts_derep = raw_contacts_cleaned.groupby([4, 7], as_index=False)[2].sum()
#raw_contacts_counts = raw_contacts_cleaned[[4,7]].value_counts()
#raw_contacts_counts_df = raw_contacts_counts.to_frame().reset_index()

#res = pd.merge(raw_contacts_derep,raw_contacts_counts_df)
#pd.DataFrame(res).to_csv("g86_raw_contacts_cleaned.csv")

#pd.DataFrame(raw_contacts_cleaned).to_csv("w76_raw_contacts_cleaned.csv")
pd.DataFrame(clean_res).to_csv('o91_normalized_df_FINAL_updated_20230922.csv')