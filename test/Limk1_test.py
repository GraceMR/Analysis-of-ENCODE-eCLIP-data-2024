#Checking if bioframe is working- looking for Limk1 binding site/3'UTR overlap
import pandas as pd
import bioframe as bf
import numpy as np

Limk1_UTR_coords = pd.read_csv('Limk1_3UTR.csv', header='infer', sep=',')
#Now give columns names. First three columns comprise coords- give same header names as for DDX6 binding site table.
##Limk1_UTR_coords.columns = ["chrom", "start", "end", "identifier", "unknown", "strand"]
#create a new identifier column with useful values (i.e. only the GENCODE transcript identifier)
Limk1_UTR_coords['GENCODE_ID'] = Limk1_UTR_coords['identifier'].str.split('_').str[0]
print(Limk1_UTR_coords.head())
#We only want the first three columns.
Limk1_UTR_coords_subset = Limk1_UTR_coords[["chrom", "start", "end"]]
#The first three columns are named as such so that the dframe is compatible with bioframe
print(Limk1_UTR_coords_subset.head())

#Now want the same info from the list of DDX6 binding sites.
Limk1_DDX6_coords = pd.read_csv('Limk1_DDX6_binding_coords.csv', header = 'infer' , sep = ',')
Limk1_DDX6_coords.reset_index(drop=True, inplace=True)
#The only two TRUE binding sites are chr7	74122298	74122353
#and chr7	74122308	74122367. These fall within the 3UTR coords. The otherwise
#binding sites are partial overlaps with the 3UTR, or completely unrelated to the 3UTR.
#Again, we only want the first three columns.
Limk1_DDX6_coords_subset = Limk1_DDX6_coords[["chrom", "start", "end"]]
print(Limk1_DDX6_coords_subset.head())
#Now that we have our two lists of genomic coordinates, we need to identify areas of overlap.

#ADDITIONAL PRE-PROCESSING
#need to order coords in both dframes in ascending order
Limk1_UTR_coords_subset = Limk1_UTR_coords_subset.copy()
Limk1_DDX6_coords_subset = Limk1_DDX6_coords_subset.copy()
Limk1_UTR_coords_subset.sort_values(by='start', ascending=True, inplace=True)
Limk1_UTR_coords_subset.reset_index(drop=True, inplace=True)
Limk1_DDX6_coords_subset.sort_values(by='start', ascending=True, inplace=True)
Limk1_DDX6_coords_subset.reset_index(drop=True, inplace=True)
print(Limk1_UTR_coords_subset.head())
print(len(Limk1_UTR_coords_subset.index))
print(Limk1_DDX6_coords_subset.head())
print(len(Limk1_DDX6_coords_subset.index))

Limk1_overlapping_intervals = bf.overlap(Limk1_UTR_coords_subset, Limk1_DDX6_coords_subset, how='inner', suffixes=('_1','_2'))
Limk1_overlapping_intervals = Limk1_overlapping_intervals.copy()
Limk1_overlapping_intervals.sort_values(by='start_1', ascending=True, inplace=True)
Limk1_overlapping_intervals.reset_index(drop=True, inplace=True)
#print the length of one of the columns- make into numpy array
print(Limk1_overlapping_intervals.head())
print(len(Limk1_overlapping_intervals.index))

#so it is able to identify that the coordinates are either entirely overlapping, partially overlapping, or not overlapping at all.
#Only returned the coordinates that are full or partial overlaps.
#The only issue came with loading the data in for the 3UTR table- it refused to recognise how the data were separated. Had
#to manually add header names.
