import subprocess
import os
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb

class Affinity_Visualization:
    
    def __init__(self,perc_rank,seq_file, allele_file,alleles_file,col_select, mutation_inds, cutoffs, my_title, output_folder, short_name):
        '''
        
        Initialize an object of class "Affinity_Visualization"
        
        
        Class members: 
            
        'seq_file' - path to the file containing the amino acids sequence(s) (string)
        'seq_df' - contains the sequence data 'self.seq_file' (pandas.DataFrame)
        'allele_freq_series' - contains the frequencies corresponding to each allele a single population (pandas.Series)
        'allele_freq_matrix' - contains the frequencies corresponfing to each allele for each population (np.Array)
        'allele_list' - contains the cleaned names of the alleles to be used (list(string))
        'mutation_inds' - positions of the mutated amino acids (either -1, or list(list(int)), OR dict(string: list(list(int)))
        'cutoffs' - nanomolar affinity thresholds to use when generating plots (list(int))
        'my_title' - title to be used on the generated plots (string)
        'output_folder' - directory where all output will be written (string)
        'output_prefix' - short name to be used as a prefix on any output files (string)
        'affinity_df' - contains all processed outputs from running netMHCIIpan (pandas.DataFrame)
        'image_mat' - reshaped version of self.affinity_df to be used in visualization (pandas.DataFrame)
        
        
        Class methods:
            
        '_get_seq_df' - initialize self.seq_df
        '_get_cutoffs' - initialize self.cutoffs
        '_get_output_folder' - initialize self.output_folder
        '_get_allele_freq_series' - initialize self.allele_freq_series
        '_get_allele_freq_matrix' - initialize self.allele_freq_matrix
        '_transform_mut_inds' - transform the self.mutation_inds from sequence format to image_mat format
	    '_transform_mut_inds_dict' - transform the self.mutation_inds if mut_inds is a dict rather than a list
        'compute_affinities' - compute the binding affinities using netMHCIIpan and store them in self.affinity_df
        'generate_image_mat' - transform self.affinity_df into self.image_mat
        'get_heatmaps' - generate and save the heatmaps
        'get_promiscuity_plots' - generate and save the promiscuity plots
        
        '''
        self.perc_rank=perc_rank
        self.seq_file = seq_file
        self.seq_df = self._get_seq_df()
        self.allele_freq_series = self._get_allele_freq_series(allele_file)
        self.col_select= col_select
        self.allele_freq_matrix = self._get_allele_freq_matrix(alleles_file,col_select)
        self.allele_list = self.allele_freq_series.index
        self.mutation_inds = mutation_inds
        self.cutoffs = self._get_cutoffs(cutoffs)
        self.my_title = my_title
        self.output_folder = self._get_output_folder(output_folder)
        self.output_prefix = short_name
        
        self.affinity_df = None
        self.image_mat = None
        self.allele_numbering = {}
        self.identity_list = []
        

    def _get_seq_df(self):
        with open(self.seq_file) as infile:
            fasta_sequences = SeqIO.parse(infile,'fasta')
            row_list = [(this_seq.id,str(this_seq.seq).upper()) for this_seq in fasta_sequences]
            seq_df = pd.DataFrame(data=row_list, columns=['Name','Sequence'])
            return seq_df
    
    
    def _get_cutoffs(self, cutoffs):
        if cutoffs == -1:
            if self.perc_rank==0:
                return [50,500]
            if self.perc_rank==1:
                return [2.0,10.0]

        else:
            if self.perc_rank==0:
                return [50,500]+cutoffs #assume that 'cuttoffs' is in the proper form if it's not -1
            if self.perc_rank==1:
                return [2.0,10.0]+cutoffs
    
    # Make the directory if it doesn't already exist, and return the path
    def _get_output_folder(self, output_folder):
        if output_folder == -1:
            new_folder = '/'.join(self.seq_file.split('/')[:-1]) + '/Output'
        else:
            new_folder = output_folder
        if not os.path.exists(new_folder):
            os.makedirs(new_folder)
        return new_folder
    
    
    def _get_allele_freq_series(self, allele_file):
        allele_freq_df = pd.read_csv(allele_file, sep=' ', header=None, names=["Allele","Freq"])
        allele_freq_series = allele_freq_df['Freq']
        #The allele format depends on the type of allele. Refer to netMHCIIpan-3.1/data/convert_pseudo.dat
        allele_freq_series.index = [str(allele).replace(':','').replace('*','_') for allele in allele_freq_df['Allele']]
        return allele_freq_series
    
    def _get_allele_freq_matrix(self, alleles_file, col_select):
        allele_freq_matrix_df = pd.read_csv(alleles_file, sep=',', header=0, index_col=0)
        if col_select!=-1:
            allele_freq_matrix = allele_freq_matrix_df.ix[:,col_select]
        else:
            allele_freq_matrix = allele_freq_matrix_df
        #The allele format depends on the type of allele. Refer to netMHCIIpan-3.1/data/convert_pseudo.dat
        # It is assumed that the allele names are in the correct format in the file allele_csv
        return allele_freq_matrix   
                
    # Transform the mutation inds so that they match with the column labels
    # Remember that we're traversing a list of lists, here
    def _transform_mut_inds(self):
        transformed_mut_inds = []
        seq_shift = 0
        for i,inds in enumerate(self.mutation_inds):
            transformed_mut_inds += [seq_shift+this_ind-8 for this_ind in inds if this_ind>7 and this_ind<len(self.seq_df['Sequence'][i])-6]
            seq_shift += (len(self.seq_df['Sequence'][i]) - 14 + 2) # remember to account for the 2 blank columns
        return transformed_mut_inds
    
    # Same as previous function, but return a list of 2-tuples, with the first entry containing the color, and
    # the second entry containing the index
    def _transform_mut_inds_dict(self):
        transformed_mut_inds = []
        for color,inds_list in self.mutation_inds.iteritems():
            seq_shift = 0
            for i,inds in enumerate(inds_list):
                transformed_mut_inds += [(color,seq_shift+this_ind-8) for this_ind in inds if this_ind>7 and this_ind<len(self.seq_df['Sequence'][i])-6]
                seq_shift += (len(self.seq_df['Sequence'][i]) - 14 + 2) # remember to account for the 2 blank columns
        return transformed_mut_inds

    # If there's already a computed affinity file, just load it rather than recomputing
    # NOTE: THIS WILL PRODUCE INCORRECT OUTPUT IF THE CONTENTS OF 'self.seq_file'
    # HAVE LINES ENDING WITH '\r\n' RATHER THAN SIMPLY '\n'
    def compute_affinities(self):
        full_outfile_path = self.output_folder+'/'+self.output_prefix+'_netMHC_results.csv'
        if os.path.isfile(full_outfile_path):
            self.affinity_df = pd.read_csv(full_outfile_path)
            self.affinity_df["Identity"] = map(str, self.affinity_df["Identity"])

        else:
            bash_command = 'for allele in ' + ' '.join(self.allele_list) + '; do netMHCIIpan -fast 1 -f \"' + self.seq_file + '\" -a $allele; done'
            process = subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE)
            line_list = process.communicate()[0].split('\n')

            
            n = 0
            for line in line_list:
                split_line = line.split()
                if len(split_line)>0 and split_line[0].isdigit():
                    if str(split_line[3]) not in self.identity_list:
                        self.identity_list.append(str(split_line[3]))

            
            trimmed_lines = []
            for line in line_list:
                split_line = line.split()
                if len(split_line)>0 and split_line[0].isdigit():
                    if (split_line[0] == '0') and (str(split_line[3]) == self.identity_list[0]):  # going off the idea that each new allele will start with window 0
                        n += 1
                    split_line.insert(1, n)
                    #if len(split_line) == 9: # add binding level to the rows lacking it
                    if len(split_line) == 10: # adjusted to 10 to account for the allele # column I inserted
                        split_line.append('')
                    trimmed_lines.append(split_line)

            
            headers = ['Pos_nmer','Allele #','Allele','Peptide','Identity','Pos_core','Core','CoreRel','Affinity_log50k','Affinity_nM','PercentRank','ExpBind','BindingLevel']
            self.affinity_df = pd.DataFrame(trimmed_lines, columns=headers)
            self.affinity_df.to_csv(full_outfile_path)
            
            # MANUALLY SET 'Pos_nmer' COLUMN TO BE INTEGER, OTHERWISE THE PIVOT STEP IN 'generate_heatmaps' GETS CONFUSED
            self.affinity_df['Pos_nmer'] = map(int, self.affinity_df['Pos_nmer'])
            # MANUALLY SET 'Affinity_nM' COLUMN TO BE FLOAT, FOR SIMILAR REASONS
            self.affinity_df['Affinity_nM'] = map(float, self.affinity_df['Affinity_nM'])
            # MANUALLY SET 'PercentRank' COLUMN TO BE FLOAT, FOR SIMILAR REASONS
            self.affinity_df['PercentRank'] = map(float, self.affinity_df['PercentRank'])
            # MANUALLY SET 'IDENTITY' COLUMN TO BE STRING
            self.affinity_df["Identity"] = map(str, self.affinity_df["Identity"])
            
    
    # Reshape computed affinities into the final heatmap
    def generate_image_mat(self):
        mat_list = []
        self.get_allele_numbering()
        self.allele_numbering["Allele #"] = "Allele name"
        name_list = []
        for key in sorted(self.allele_numbering.keys()):
            name_list.append(self.allele_numbering[key])
        for i,row in self.seq_df.iterrows():
            if self.perc_rank==0:
                this_mat = self.affinity_df.loc[self.affinity_df['Identity']==row['Name']].pivot(index='Allele #', columns='Pos_nmer', values='Affinity_nM')
            #Index contains duplicate entries, cannot reshape
            if self.perc_rank==1:
                this_mat = self.affinity_df.loc[self.affinity_df['Identity']==row['Name']].pivot(index='Allele #', columns='Pos_nmer', values='PercentRank')
            
            ### JA's code here ###
            
            headers = {"Allele #":"Allele"}
            for position in this_mat.columns:
                headers[position] = row["Sequence"][position+7]

            this_mat.rename(columns = headers, inplace = True)

            blank_mat = pd.DataFrame(data=5e4*np.ones((this_mat.shape[0],2)), columns=['...','...'], index=this_mat.index)
            if len(mat_list) == 0:
                mat_list.append(this_mat)
            else:
                mat_list.append(blank_mat)
                mat_list.append(this_mat)    
        self.image_mat = pd.concat(mat_list, axis=1)    
        self.image_mat.to_csv(self.output_folder+'/'+self.output_prefix+'_reshaped_affinity_df.csv')

    def get_allele_numbering(self):
        self.allele_numbering = {}
        numbers = self.affinity_df["Allele #"].unique()
        for n in numbers:
            allele_name = self.affinity_df.loc[self.affinity_df["Allele #"] == n,"Allele"].unique()[0]
            self.allele_numbering[n] = allele_name
       
    # Generate and save the len(self.cutoffs) many heatmaps
    def get_heatmaps(self):

        for this_cutoff in self.cutoffs:
            # Make the heatmap
            fig, ax = plt.subplots()
            this_heatmap = ax.pcolor(self.image_mat, cmap=plt.cm.hot, vmin=0, vmax=int(1.1*this_cutoff))
            
            
            # Construct colorbar legend
            if self.perc_rank==0:
                these_ticks = range(0, int(1.1*this_cutoff), int(this_cutoff/5)+1 )
                cbar = fig.colorbar(this_heatmap, ticks=these_ticks, orientation='vertical', pad=0.05)#KK:pad controls hm-colorbar distance
                cbar.ax.set_yticklabels([str(this_tick)+'nM' for this_tick in these_ticks])
            else:
                these_ticks = range(0, int(this_cutoff+1), int(this_cutoff/5)+1)
                cbar = fig.colorbar(this_heatmap, ticks=these_ticks, orientation='vertical', pad=0.05)#KK:pad controls hm-colorbar distance
                cbar.ax.set_yticklabels([str(this_tick)+'%' for this_tick in these_ticks])
        

            # remove the unnecessary extra space at the edges
            plt.xlim(0, self.image_mat.shape[1])
            plt.ylim(0, self.image_mat.shape[0])
            
            # put the major ticks at the middle of each cell
            ax.set_xticks(np.arange(self.image_mat.shape[1]) + 0.5, minor=False)
            ax.set_yticks(np.arange(self.image_mat.shape[0]) + 0.5, minor=False)
            ax.invert_yaxis()
            ax.xaxis.tick_top()
            
            # tick labels
            ax.set_xticklabels(self.image_mat.columns, minor=False) # JA: testing quick fix
            ax.set_yticklabels([self.allele_numbering[num] for num in self.image_mat.index], minor=False, size=11) # JA: testing quick fix
            
            # hide tick marks
            for t in ax.xaxis.get_ticklines(): t.set_visible(False) 
            for t in ax.yaxis.get_ticklines(): t.set_visible(False) 
            
            # Add a line dividing each of the segmented chunks (between the two ellipses)
            split_inds = np.where(self.image_mat.columns == '...')[0]
            split_inds = [ind for i,ind in enumerate(split_inds) if i%2==1]# we only want every-other entry
            if len(split_inds)>0: # i.e. if #sequences > 1
                plt.vlines(x=split_inds, ymin=0, ymax=self.image_mat.shape[0])
            
            # Transform mutation inds, and color the tick labels as needed
            # Take the proper steps depending upon whether "self.mutation_inds" is a list or a dict (or neither)
            if self.mutation_inds != -1:
                if type(self.mutation_inds) == list:
                    for tick_ind in self._transform_mut_inds():
                        ax.get_xticklabels()[tick_ind].set_color('red')
                elif type(self.mutation_inds) == dict:
                    for color,tick_ind in self._transform_mut_inds_dict():
                        ax.get_xticklabels()[tick_ind].set_color(color)
            
            # Generate the plot's title
            plt.suptitle(self.my_title, fontsize=16)
            
            # Size and save the plot; we want the plot dimensions 
            # to be roughly proportional to self.image_mat, with .25 inches per row, and .20 inches per column
            my_dpi = 150
            if (self.image_mat.shape[1]*.20 + 5.5)*my_dpi > 32000: # official limit of 32768 pixels, computed as width*dpi
                fig.set_size_inches(32000/my_dpi, self.image_mat.shape[0]*.25 + 1.5) # cap the image size
            else:
                fig.set_size_inches(self.image_mat.shape[1]*.20 + 5.5, self.image_mat.shape[0]*.25 + 1.5)

            # Save and close the plot
            fig.set_tight_layout({'rect': (0, 0, 1, 0.95)})
            if self.perc_rank==0:
                fig.savefig(self.output_folder+'/'+self.output_prefix+'_Heatmap_'+str(this_cutoff)+'nM.png', dpi=my_dpi)
            if self.perc_rank==1:
                fig.savefig(self.output_folder+'/'+self.output_prefix+'_Heatmap_'+str(this_cutoff)+'%.png', dpi=my_dpi)
            plt.close()
            
    
    # Generate and save the len(self.cutoffs) many promiscuity plots
    def get_promiscuity_plots(self):
        
        total_freq = np.sum(self.allele_freq_series)
        for this_cutoff in self.cutoffs:
            # Get y-values for promiscuity plot, and mask the elipses
            binary_mat = self.image_mat <= this_cutoff

            freq_array = np.array(self.allele_freq_series).dot(binary_mat)
            
            freq_array = np.array(freq_array) / total_freq
            freq_array = 1-(1-freq_array)**2
            
            freq_array=100*freq_array
                      
            masked_freq_array = np.ma.masked_array(freq_array, self.image_mat.columns=='...')

            # Make the plot
            fig, ax = plt.subplots()
            ax.plot(masked_freq_array, lw=3)
            
            # Scale and label the axes
            plt.xlim(0, len(freq_array))
            plt.ylim(0, 100)
            plt.ylabel('Promiscuity Score', size=18)

            # Place the tick marks as needed
            ax.set_xticks(np.arange(self.image_mat.shape[1]), minor=False)
            ax.set_yticks(range(0,110,10), minor=False)
            
            # tick labels
            ax.set_xticklabels(self.image_mat.columns, minor=False)
            ax.set_yticklabels(map(str,range(0,110,10)), minor=False, size=20)
            
            # hide tick marks on x-axis
            for t in ax.xaxis.get_ticklines(): t.set_visible(False) 
            
            # Add a line dividing each of the segmented chunks (between the two ellipses)
            split_inds = np.where(self.image_mat.columns == '...')[0]
            split_inds = [ind for i,ind in enumerate(split_inds) if i%2==1]# we only want every-other entry
            if len(split_inds)>0: # i.e. if #sequences > 1
                plt.vlines(x=split_inds, ymin=0, ymax=100)
            
            # Transform mutation inds, and color the tick labels as needed
            # Take the proper steps depending upon whether "self.mutation_inds" is a list or a dict (or neither)
            if self.mutation_inds != -1:
                if type(self.mutation_inds) == list:
                    for tick_ind in self._transform_mut_inds():
                        ax.get_xticklabels()[tick_ind].set_color('red')
                elif type(self.mutation_inds) == dict:
                    for color,tick_ind in self._transform_mut_inds_dict():
                        ax.get_xticklabels()[tick_ind].set_color(color)
                                    
            # Generate the plot's title
            plt.title(self.my_title+' Promiscuity Plot ('+str(this_cutoff)+'nM cutoff)', fontsize=18)
            plt.title(self.my_title, fontsize=16)

            # Size and save the plot; we want the plot to be twice as wide as it is tall,
            # with roughly .20 inches per letter
            my_dpi = 150
            if (self.image_mat.shape[1]*.20 + 2)*my_dpi > 32000: # official limit of 32768 pixels, computed as width*dpi
                fig.set_size_inches(32000/my_dpi, 6400/my_dpi) # cap the image size; height=width/5
            else:
                fig.set_size_inches(self.image_mat.shape[1]*.20 + 2, self.image_mat.shape[1]*.20 / 2 + 2)
            
            # Save and close the plot
            fig.set_tight_layout(True)
            if self.perc_rank==0:
                fig.savefig(self.output_folder+'/'+self.output_prefix+'_Promiscuity_'+str(this_cutoff)+'nM.png', dpi=my_dpi)
            if self.perc_rank==1: 
                fig.savefig(self.output_folder+'/'+self.output_prefix+'_Promiscuity_'+str(this_cutoff)+'%.png', dpi=my_dpi)   
            plt.close()

    def get_multi_promiscuity_plots(self):
        
        total_freq = np.sum(self.allele_freq_matrix, axis=0)
        for this_cutoff in self.cutoffs:
            # Get y-values for promiscuity plot, and mask the elipses
            binary_mat = self.image_mat <= this_cutoff
            #freq_matrix_array = self.allele_freq_series.dot(binary_mat)
            freq_matrix_array = np.transpose(np.array(self.allele_freq_matrix)).dot(binary_mat)
            #freq_matrix_array = np.divide(total_freq,np.array(freq_matrix_array)) # / total_freq
            freq_matrix_array = 100 * freq_matrix_array / total_freq[:, np.newaxis]

            plt.plot(freq_matrix_array)
            freq_matrix_array=np.transpose(freq_matrix_array)
            
            # Make the plot
            fig, ax = plt.subplots()
            ax.plot(freq_matrix_array, lw=1.5)
            
            # Scale and label the axes
            plt.xlim(0, len(freq_matrix_array))
            plt.ylim(0, 100)
            plt.ylabel('Promiscuity Score', size=18)

            # Place the tick marks as needed
            ax.set_xticks(np.arange(self.image_mat.shape[1]), minor=False)
            ax.set_yticks(range(0,110,10), minor=False)
            
            # tick labels
            ax.set_xticklabels(self.image_mat.columns, minor=False)
            ax.set_yticklabels(map(str,range(0,110,10)), minor=False, size=20)
            
            # hide tick marks on x-axis
            for t in ax.xaxis.get_ticklines(): t.set_visible(False) 
            
            # Add a line dividing each of the segmented chunks (between the two ellipses)
            split_inds = np.where(self.image_mat.columns == '...')[0]
            split_inds = [ind for i,ind in enumerate(split_inds) if i%2==1]# we only want every-other entry
            if len(split_inds)>0: # i.e. if #sequences > 1
                plt.vlines(x=split_inds, ymin=0, ymax=100)
            
            # Transform mutation inds, and color the tick labels as needed
            # Take the proper steps depending upon whether "self.mutation_inds" is a list or a dict (or neither)
            if self.mutation_inds != -1:
                if type(self.mutation_inds) == list:
                    for tick_ind in self._transform_mut_inds():
                        ax.get_xticklabels()[tick_ind].set_color('red')
                elif type(self.mutation_inds) == dict:
                    for color,tick_ind in self._transform_mut_inds_dict():
                        ax.get_xticklabels()[tick_ind].set_color(color)
            
            plt.legend(total_freq.index, prop={'size':10})

            # Generate the plot's title
            if self.perc_rank==0:
                plt.suptitle(self.my_title+' Promiscuity Plot ('+str(this_cutoff)+'nM cutoff)', fontsize=18)
            if self.perc_rank==1:
                plt.suptitle(self.my_title+' Promiscuity Plot ('+str(this_cutoff)+'% cutoff)', fontsize=18)
            plt.suptitle(self.my_title)

            #Size and save the plot; we want the plot to be twice as wide as it is tall,
            # with roughly .20 inches per letter
            my_dpi = 150
            if (self.image_mat.shape[1]*.20 + 2)*my_dpi > 32000: # official limit of 32768 pixels, computed as width*dpi
                fig.set_size_inches(32000/my_dpi, 6400/my_dpi) # cap the image size; height=width/5
            else:
                fig.set_size_inches(self.image_mat.shape[1]*.20 + 2, self.image_mat.shape[1]*.20 / 2 + 2)
            
            # Save and close the plot
            fig.set_tight_layout(True)
            if self.perc_rank==0:
                fig.savefig(self.output_folder+'/'+self.output_prefix+'_Multi_Promiscuity_'+str(this_cutoff)+'nM.png', dpi=my_dpi)
            if self.perc_rank==1:
                fig.savefig(self.output_folder+'/'+self.output_prefix+'_Multi_Promiscuity_'+str(this_cutoff)+'%.png', dpi=my_dpi)

            plt.close()
            