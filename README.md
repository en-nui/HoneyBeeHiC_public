# Honey bee HiC python scripts

Author: Chris Robinson
# Bioinformatic pipeline used for data analysis of HiC data. Includes a summary of the pipeline and some of the important Python scripts. Please contact me for any requests or questions regarding how something was done.
# Please refer to manuscript methods section first regarding software used and flags.

1) After QC, trimming, assembly, and binning, vMAGs were called for SNVs using Anvi'o (variability-profile). SNV and allele frequencies were used as input for the [population genetic scripts] (https://github.com/en-nui/HoneyBeeHiC_public/blob/main/popgen_annotations_summary_statistics_program.py)

2) vMAGs were used to recruit freads from each genome to identify low coverage windows from [this script] (https://github.com/en-nui/HoneyBeeHiC_public/blob/main/sliding_window_genome_island_viruses_WIP.py).

3) HiC noise calculations were calculating using this [script] (https://github.com/en-nui/HoneyBeeHiC_public/blob/main/hic_noise_calculations.py).

4) Low coverage vMAG windows were combined with annotations using this [script] (https://github.com/en-nui/HoneyBeeHiC_public/blob/main/updated_island_counts_USEME.py). 
