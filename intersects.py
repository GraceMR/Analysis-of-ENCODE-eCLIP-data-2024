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
UTR_coords_subset = UTR_coords_subset.copy()
DDX6_coords_subset = DDX6_coords_subset.copy()
UTR_coords_subset.sort_values(by='start', ascending=True, inplace=True)
DDX6_coords_subset.sort_values(by='start', ascending=True, inplace=True)
print(UTR_coords_subset.head())
print(len(UTR_coords_subset.index))
print(DDX6_coords_subset.head())
print(len(DDX6_coords_subset.index))

import bioframe as bf

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
#columns correspond to 3'UTR coordinates from across the genome. 

overlapping_intervals.drop(columns=['chrom_1', 'start_1', 'end_1'], inplace = True)
#try re-naming the columns
overlapping_intervals.rename({'chrom_2': 'chrom', 'start_2': 'start', 'end_2': 'end'}, axis=1, inplace=True)
overlapping_intervals = overlapping_intervals.drop_duplicates(subset=['chrom', 'start', 'end'], keep=False)
print(overlapping_intervals.head())
print(len(overlapping_intervals.index))
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
additional_values_df = overlapping_intervals.apply(get_additional_columns, axis=1)
overlapping_intervals[['strand', 'confidence_score', 'sample_or_tissue_used']] = additional_values_df.apply(pd.Series)
print(overlapping_intervals.head())

#it works! Now need to figure out how to get gene/transcript names from the coordinates