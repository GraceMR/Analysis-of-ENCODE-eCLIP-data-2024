####Combining outputs of intersects.py and biomart_script.r####
#The purpose of this script is to exclusively produce a list of transcripts bound by DDX6 in their 3'UTRs
#with associated gene names/descriptions and original experiment information, derived from the POSTAR3 study
import pandas as pd


DDX6_overlapping_intervals = pd.read_csv('DDX6_overlapping_intervals.csv', sep = ',')
print(DDX6_overlapping_intervals.head())
print(len(DDX6_overlapping_intervals.index))

gene_info = pd.read_csv('R_analysis/output/overlapping_intervals_gene_names.csv', sep = ',')
print(gene_info.head())
print(len(gene_info.index))

print(set(DDX6_overlapping_intervals['start']) - set(gene_info['start']))
print(set(DDX6_overlapping_intervals['end']) - set(gene_info['end']))

#merged_df = pd.merge(gene_info, DDX6_overlapping_intervals, on=['start', 'end', 'chrom'], how='inner')
#print(merged_df.head())
#print(len(merged_df.index))