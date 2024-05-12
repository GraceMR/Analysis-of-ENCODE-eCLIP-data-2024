###Producing a list of DDX6 binding sites in 3'UTRs across the genome###
#Need to compare DDX6 list of binding coordinates with list of 3UTR coordinates. Ensure both starting lists are in the same format.
import pandas as pd
UTR_coords = pd.read_csv('3UTR.csv', header = None, sep = '\t')
#Now give columns names. First three columns comprise coords- give same header names as for DDX6 binding site table.
UTR_coords.columns = ["chrom", "start", "end", "identifier", "unknown", "strand"]
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
UTR_coords_subset.sort_values(by='start', ascending=True, inplace=True)
DDX6_coords_subset.sort_values(by='start', ascending=True, inplace=True)
print(UTR_coords_subset.head())
print(DDX6_coords_subset.head())

import bioframe as bf
overlapping_intervals = bf.overlap(UTR_coords_subset, DDX6_coords_subset, how='inner', suffixes=('_1','_2'))
print(overlapping_intervals.head())
#_1 = UTR coords; _2 = DDX6_coords. Should only return coords in both dframes where there is overlap.

#Now need to append additional important info to the overlapping_intervals dframe, including 'sample_or_tissue_used',
#'confidence_score', 'strand' (dependent on _2 values); and gene name.
#May be easiest to create a custom function for this.

def get_additional_columns(row):
    coords = row[['chrom', 'start', 'end']]
    matching_row = DDX6_coords_subset[DDX6_coords_subset[['chrom', 'start', 'end']] == coords]
    if not matching_row.empty:
        return matching_row[['strand', 'confidence_score', 'sample_or_tissue_used']].values[0]
    else:
        return None

# Apply the custom function to each row in dframeB
#overlapping_intervals[['strand', 'confidence_score', 'sample_or_tissue_used']] = overlapping_intervals.apply(get_additional_columns, axis=1)
