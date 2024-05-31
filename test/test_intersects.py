#Need to create functions to test elements of intersets.py; this is normally to test functions, but it is not necessarily functions we want to test.
#most important parts to test where things might go wrong:

overlapping_intervals = bf.overlap(UTR_coords_subset, DDX6_coords_subset, how='inner', suffixes=('_1','_2'))

overlapping_intervals.sort_values(by='start_1', ascending=True, inplace=True)

#Above: want to check that bf.overlap yields appropriate overlaps. Also want to check that result is sorted in ascending order.
#there are a few ways it could return data.
#1: full overlap- both start and end of binding sites within UTR.
#2: full overlap 2- neither start or end of binding site within UTR, but UTR is within the binding site.
#3: partial overlap: only start coord in the UTR.
#4: partial overlap: only end coord in the UTR.
#5: no overlap: neither start nor end of binding site within UTR- UTR not within binding site.
#Start with the above.

#TEST_INPUTS
#3UTR_test = input test 3UTR coordinates
#DDX6_test = input test DDX6 binding site coordinates (reflects two examples of each of the above 1-5)
#TEST_OUTPUTS
#looks like the bioframe.overlap() function returns all full and partial overlaps.
#all_full_partial_overlaps = first possible output (all full + partial overlaps)
#full_overlaps_only = second possible output (all full overlaps only)
#partial_overlaps_only = third possible output (all partial overlaps)

overlapping_intervals = overlapping_intervals.drop_duplicates(subset=['chrom_2', 'start_2', 'end_2'], keep=False)
#not sure if maybe this line is introducing error? Why use it?
