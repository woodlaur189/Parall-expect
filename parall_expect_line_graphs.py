#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 14:44:39 2021

@author: lwoo0005
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 12:50:18 2021

@author: lwoo0005
"""

 
import pandas
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from glob import glob

#Input folder containing all Parall-expect runs of interest
for folder in glob('/Users/lwoo0005/Documents/Laura_stuff/Thesis/Jake_YPE_parallel_expect/yeast/*/'):
    property_file=folder+"Properties_of_all_genes.csv"
    gene_tally_file=str(glob(folder+"All*tallies.csv")[0])
    print(property_file)
    print(gene_tally_file)
    properties=pandas.read_csv(property_file)
    gene_tallies=pandas.read_csv(gene_tally_file)
    out_folder=folder
    data_for_corr=properties.merge(gene_tallies[['Gene','Maximum null value','Experimental hits']], on="Gene")
    
    by_max=data_for_corr.groupby('Maximum null value',as_index='False')
    by_experimental=data_for_corr.groupby('Experimental hits',as_index='False')

    fig, ax = plt.subplots(1,1, figsize=(10, 6))
    gene_prop = properties.columns[1]
    ax.errorbar(by_max.mean().index.to_list(),by_max[gene_prop].mean(),yerr=[i*1.96 for i in by_max[gene_prop].sem().fillna(0).tolist()],label="Mean for null maximum hit +/- 95% CI",color='grey')
    ax.errorbar(by_experimental.mean().index.to_list(),by_experimental[gene_prop].mean(),yerr=[i*1.96 for i in by_experimental[gene_prop].sem().fillna(0).tolist()],label="Mean for experimental hit +/- 95% CI",color='red')
    ax.set_xlim(left=1)
    ax.set_ylabel(str(gene_prop), size=30)
    ax.set_xlabel('Number of populations with gene hit', size=30)
    ax.xaxis.label.set_size(30)
    ax.yaxis.label.set_size(30)
    ax.tick_params(axis='both', which='major', labelsize=30)
    max_by_max=max(by_max[gene_prop].mean().tolist())
    max_by_experimental=max(by_experimental[gene_prop].mean().tolist())
    max_y=max([max_by_max,max_by_experimental])             
    ax.legend(fontsize=25)
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Avenir']
    plt.tight_layout()
    plt.savefig(out_folder+'/gene_length_line_plots.png',dpi=300)
    
    plt.show()
 