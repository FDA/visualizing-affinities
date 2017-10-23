#!/usr/bin/env python2.7

from Affinity_Visualization import Affinity_Visualization
import sys
    
if __name__ == "__main__": 
    
    param_file_path = './parameters.py' # expects a path to a .py file

    # Loads parameters from parameters.py
    
    with open(param_file_path) as infile:
        for line in infile.readlines():
            exec line

    print "\nParameter file has been loaded"

    My_Affinity_Visualization = Affinity_Visualization(perc_rank,seq_file_path, HLA_file_path, allele_csv, columns, mutation_inds, cutoffs, my_title, output_folder_path, seq_name)
    print "Heatmap object has been created"
    
    My_Affinity_Visualization.compute_affinities()
    print "Affinities have been computed and saved"
    
    My_Affinity_Visualization.generate_image_mat()
    print "Affinities have been reshaped and saved"
    
    My_Affinity_Visualization.get_heatmaps()
    print "Heatmaps have been generated and saved"
    
    My_Affinity_Visualization.get_promiscuity_plots()
    print "Promiscuity plots have been generated and saved"

    My_Affinity_Visualization.get_multi_promiscuity_plots()
    print "Promiscuity plots have been generated and saved"

    print "\nAll Done!\n"
    
