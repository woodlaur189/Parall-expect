#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 23:08:24 2021

@author: lwoo0005
"""


import pandas 
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from Bio import SeqIO
import numpy as np
from glob import glob

#Make a loop to do each chromosome
genes=[]
av_rep_times=[]
for ref_seq in glob("/Users/lwoo0005/Documents/Laura_stuff/Thesis/Yeast_s288c/chrm*.gb"):
    chrm=ref_seq.split("/")[-1].split(".gb")[0]
    print(chrm)
    rep_file="/Users/lwoo0005/Documents/Laura_stuff/Thesis/Yeast_s288c/"+chrm+"_yeast_rep_timing.csv"
    with open(rep_file, 'r') as csv_file:
        rep_data=(pandas.read_csv(csv_file, delimiter='\t', encoding='utf-8',header =0, names=['Kb','Time'],usecols=[1,0])).dropna()  
    csv_file.close()
    pos = [x*1000 for x in rep_data['Kb'].tolist()]
    times = [x for x in rep_data['Time'].tolist()]
    #f = interp1d(pos, times, kind='cubic')
    f=InterpolatedUnivariateSpline(pos, times, k=3)
    records = [rec for rec in SeqIO.parse(ref_seq, "genbank")]
    starts=[]
    ends=[]
    #genes=[]
    #av_rep_times=[]
    for record in records:
        for feature in record.features:
            if feature.type=="CDS":
                try:
                    gene_name=str(feature.qualifiers['gene'][0])
                except KeyError:
                    gene_name=str(feature.qualifiers['locus_tag'][0])
                starts.append(int(feature.location.start))
                ends.append(int(feature.location.end))
                genes.append(gene_name)
                underscore_gene_name=gene_name.replace('-','_')
                av_rep_times.append(np.mean(f(range(int(feature.location.start),int(feature.location.end),1))))
 
data_tuples = list(zip(genes,av_rep_times))
out_data = pandas.DataFrame(data_tuples, columns=['Gene','Average_replication_time'])

#Making simple gene,rep_time output
out_data.to_csv('/Users/lwoo0005/Documents/Laura_stuff/Thesis/Yeast_s288c/Gene_avg_rep_time.csv',index=False)

"""  
import matplotlib.pyplot as plt
plt.plot(range(0,max(pos),1), f(range(0,max(pos),1)), '-')
plt.legend(['cubic'], loc='best')
plt.show()
"""