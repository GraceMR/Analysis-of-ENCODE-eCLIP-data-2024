###Producing a list of DDX6 binding sites in 3'UTRs across the genome###
#Need to compare DDX6 list of binding coordinates with list of 3UTR coordinates. Ensure both starting lists are in the same format.
import pandas as pd
import bioframe as bf
import numpy as np

UTR_coords = pd.read_csv('3UTR.csv', header = None, sep = '\t')
#Now give columns names. First three columns comprise coords- give same header names as for DDX6 binding site table.
UTR_coords.columns = ["chrom", "start", "end", "identifier", "unknown", "strand"]
#create a new identifier column with useful values (i.e. only the GENCODE transcript identifier)
UTR_coords['GENCODE_ID'] = UTR_coords['identifier'].str.split('_').str[0]
print(UTR_coords.head())
#We only want the first three columns.
UTR_coords_subset = UTR_coords[["chrom", "start", "end"]]
#The first three columns are named as such so that the dframe is compatible with bioframe
print(UTR_coords_subset.head())

#Now want the same info from the list of DDX6 binding sites.
DDX6_coords = pd.read_csv('output/DDX6_binding_coords.csv', header = 'infer' , sep = ',')
#Again, we only want the first three columns.
DDX6_coords_subset = DDX6_coords[["chrom", "start", "end"]]
print(DDX6_coords_subset.head())
#Now that we have our two lists of genomic coordinates, we need to identify areas of overlap.

#ADDITIONAL PRE-PROCESSING
#need to order coords in both dframes in ascending order
UTR_coords_subset = UTR_coords_subset.copy()
DDX6_coords_subset = DDX6_coords_subset.copy()
UTR_coords_subset.sort_values(by='start', ascending=True, inplace=True)
DDX6_coords_subset.sort_values(by='start', ascending=True, inplace=True)
print(UTR_coords_subset.head())
print(len(UTR_coords_subset.index))
print(DDX6_coords_subset.head())
print(len(DDX6_coords_subset.index))

overlapping_intervals = bf.overlap(UTR_coords_subset, DDX6_coords_subset, how='inner', suffixes=('_1','_2'))
overlapping_intervals = overlapping_intervals.copy()
overlapping_intervals.sort_values(by='start_1', ascending=True, inplace=True)
#print the length of one of the columns- make into numpy array
print(overlapping_intervals.head())
print(len(overlapping_intervals.index))
#_1 = UTR coords; _2 = DDX6_coords. Should only return coords in both dframes where there is overlap.
#Contents of overlapping_intervals seems consistent. However, it looks like the row number changes every time. Not sure what this means.

#Now need to append additional important info to the overlapping_intervals dframe, including 'sample_or_tissue_used',
#'confidence_score', 'strand' (dependent on _2 values); and gene name.
#May be easiest to create a custom function for this.

#May need to drop first three columns from overlapping_intervals for the function to work; the first three
#columns correspond to 3'UTR coordinates from across the genome. Still need to retain first three columns for additional
#data filtering.

overlapping_intervals = overlapping_intervals.drop_duplicates(subset=['chrom_2', 'start_2', 'end_2'], keep=False)
DDX6_overlapping_intervals = overlapping_intervals.filter(regex='_2')  
UTR_overlapping_intervals = overlapping_intervals.filter(regex='_1')

#try re-naming the columns
DDX6_overlapping_intervals = DDX6_overlapping_intervals.copy()
UTR_overlapping_intervals = UTR_overlapping_intervals.copy()
DDX6_overlapping_intervals.rename({'chrom_2': 'chrom', 'start_2': 'start', 'end_2': 'end'}, axis=1, inplace=True)
UTR_overlapping_intervals.rename({'chrom_1': 'chrom', 'start_1': 'start', 'end_1': 'end'}, axis=1, inplace=True)
print(DDX6_overlapping_intervals.head())
print(len(DDX6_overlapping_intervals.index))
print(UTR_overlapping_intervals.head())
print(len(UTR_overlapping_intervals.index))
#done! Seems to be 3201 3'UTR binding sites

def get_additional_columns(row):
    coords = (row['chrom'], row['start'], row['end'])
    matching_row = DDX6_coords[
        (DDX6_coords['chrom'] == coords[0]) &
        (DDX6_coords['start'] == coords[1]) &
        (DDX6_coords['end'] == coords[2])
    ]
    if not matching_row.empty:
        strand = matching_row['strand'].values[0]
        confidence_score = matching_row['confidence_score'].values[0]
        sample_or_tissue_used = matching_row['sample_or_tissue_used'].values[0]
        additional_values = {
            'strand': strand,
            'confidence_score': confidence_score,
            'sample_or_tissue_used': sample_or_tissue_used
        }
        return additional_values
    else:
        return None

# Apply get_additional_columns to overlapping_intervals
additional_values_df = DDX6_overlapping_intervals.apply(get_additional_columns, axis=1)
DDX6_overlapping_intervals[['strand', 'confidence_score', 'sample_or_tissue_used']] = additional_values_df.apply(pd.Series)
print(DDX6_overlapping_intervals.head())
print(len(DDX6_overlapping_intervals.index))

#it works! Now need to figure out how to get gene/transcript names from the coordinates
#The easiest way to do this might be to take the identifier value in the original UTR database.
#need to remove any value after '_' in the identifier- anything after this is superfluous- DONE
#now need to make a new function that does the same as get_additional_columns, but with UTR_coords

def get_GENCODE_ID (row):
    coords = (row['chrom'], row['start'], row['end'])
    matching_row = UTR_coords[
        (UTR_coords['chrom'] == coords[0]) &
        (UTR_coords['start'] == coords[1]) &
        (UTR_coords['end'] == coords[2])
    ]
    if not matching_row.empty:
        GENCODE_ID = matching_row['GENCODE_ID'].values[0]
        return GENCODE_ID
    else:
        return None

gencode_values_df = UTR_overlapping_intervals.apply(get_GENCODE_ID, axis=1)
UTR_overlapping_intervals['GENCODE_ID'] = gencode_values_df.apply(pd.Series)
print(UTR_overlapping_intervals.head())
DDX6_overlapping_intervals.to_csv('DDX6_overlapping_intervals.csv', index=True)
UTR_overlapping_intervals.to_csv('UTR_overlapping_intervals.csv', index=True)

####USING BIOMART TO CONVERT GENCODE/ENSEMBL IDS TO GENE NAMES
#Easiest to accomplish this in R using the Bioconductor package Biomart; used UTR_overlapping_intervals.
#UTR_overlapping_intervals.to_csv('UTR_overlapping_intervals.csv', index=False)
#Output of biomart_script.r = overlapping_intervals_gene_names

####COMBINING GENE INFO WITH ORIGINAL POSTAR3 EXPERIMENTAL INFO####
#need to merge overlapping_intervals_gene_names with DDX6_overlapping_intervals. This should possibly be done in a separate script.
#Produce a .csv file of DDX6_overlapping_intervals
#IMPORTANT- need to initally add 3'utr binding coordinates to DDX6_overlapping intervals; should be two sets of coordinates.
#OR retain index number when exporting

