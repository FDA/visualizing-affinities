###############################
#### Parameter descriptions ###
###############################
#   
#'perc_rank'
#    1 to use percent rank, 0 to use affinity(nM)

perc_rank=1

#'seq_file_path':
#    Path to FASTA file containing the (segmented) sequence
#    Every sequence within the FASTA file must have a unique name
#

seq_file_path = "./test.fasta"

#'HLA_file_path':
#    Path to a text file where each line contains an allele name, followed by a space, followed by a percent frequency
#    For example, one row could be "DRB1*0401 6.64"
#    

HLA_file_path = "./World.txt"

#'allele_csv'
#	path to a *.csv file containing Column headings with different population subgroups. each line after the heading contains an allele name in the format DRB1_0101 
#   followed by comma seperated values for each population subtype.

allele_csv= "./Nworld.csv"

#'columns'
#	Choose which columns of the csv above to use in your plot. 
#	6 - World
#	0 - African
#	1 - European
# 	2 - Native American
#	3 - North American
#	4 - Northeast Asian
#	5 - ASountheast Asian
# 
#	This variable must still be a list even for one population. ex. 'columns=[6]' Will produce a chart identical to the world only population percentages. 
columns=[6,0,1,2,3,4,5]

#'output_folder_path':
#    Either '-1', or the path to the directory where all outputs will be stored (the directory will be created if it does not already exist)
#    If '-1', then a directory called "Output" will be created in the same location as the sequence file
#    

output_folder_path = "./output"

#'mutation_inds':
#    Either '-1', or a list of lists giving the 1-indexed mutation locations of each mutation ex. [[15,17]]
#    

mutation_inds = -1

#    Note that this variable has recently been extended to allow for named colors to be passed to the plotting function in the form of dict
#    Set each desired color to be a key in a dict, with it's corresponding value being the list-of-lists containing the matching indices 
#


#'cutoffs':
#	Either '-1', or a list of threshold values to use when constructing heatmaps and promiscuity plots
#	If '-1' cutoffs will either be 50nM and 500nM or 2.0% and 10.0% depending on the value of 'perc_rank' above
#	If additional values are desired, set 'cutoffs' to a list of values not including the defaults. i.e if percent ranks of 0.5%,2.0% and 10.0%
#	are desired, use 'cutoffs=[0.5]'

cutoffs = -1

#'my_title':
#    The full title that you wish to appear on the plots
#    
my_title = "Title"

#'seq_name':
#    A short identifier to be used when saving the output; should be related to the sequence
#
seq_name = "seq_name"




