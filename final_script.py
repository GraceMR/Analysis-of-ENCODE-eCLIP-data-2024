####Combining outputs of intersects.py and biomart_script.r####
#The purpose of this script is to produce a final list of transcripts bound by DDX6 in their 3'UTRs
#with associated gene names/descriptions and original experiment information, derived from the POSTAR3 study
import pandas as pd

final_table = pd.read_csv('final_table.csv', sep = ',')
print(final_table.head())
print(len(final_table.index))

gene_info = pd.read_csv('R_analysis/output/overlapping_intervals_gene_names.csv', sep = ',')
print(gene_info.head())
print(len(gene_info.index))

merged_df = pd.merge(gene_info, final_table, on='X', how='inner')
print(merged_df.head())
print(len(merged_df.index))
#consistent length- structure looks good. Can now remove 'X' column and rename the chrom, start and end columns. Can also
#possibly re-order the columns

merged_df = merged_df.drop('X', axis=1)

#merged_df.rename(columns={'chrom_x': '3UTR_chrom', 'start_x': '3UTR_start', 'end_x': '3UTR_end', 'chrom_y': 'DDX6_chrom', 'start_y': 'DDX6_start', 'end_y': 'DDX6_end'}, inplace=True)
#merged_df = merged_df[['DDX6_chrom', 'DDX6_start', 'DDX6_end', '3UTR_chrom', '3UTR_start',
    #'3UTR_end', 'ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'description', 'strand', 'confidence_score', 'sample_or_tissue_used']]

print(merged_df.head())

#export as .csv
#merged_df.to_csv('DDX6_3UTR_binding_sites.csv', index=False)
#DONE
