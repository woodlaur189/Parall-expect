
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 17:25:30 2020
@author: lwoo0005
"""

#Import relevant modules
import csv
from collections import Counter
from numpy.random import choice
import pandas
import os, errno
import scipy.stats as st
from numpy import std, average, unique
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
import matplotlib
import matplotlib.pyplot as plt


ref_seq = str(raw_input("Input the path to the genome reference file.\n"))
num_muts = int(raw_input("How many mutations?\n"))
num_reps = float(raw_input("How many simulations should be run?\n"))
observed1 = str(raw_input("Please input the path to the first observation CSV file in which genes are in row 1 and frequencies in row 2.\n"))
out_folder = str(raw_input("Where do you want results to go? Enter a path.\n"))

#Beta: synonymous mutations
syn_data='No'

rep_timings='No'
rep_times='/Users/lwoo0005/Documents/Laura_stuff/Thesis/Yeast_s288c/Gene_avg_rep_time.csv'

#Overrides other options
emp_mut_vars='No'
empirical_weight_data=''

use_only_mutable='Yes'
use_length='Yes'
#Modified by gene length:
#Replication times
#Fraction NS sites

#Not modified by gene length:
#Gene-specific mutation rates

if rep_timings=='Yes':
    with open(rep_times, 'r') as csv_file:
        rep_time_data=(pandas.read_csv(csv_file, delimiter=',', encoding='utf-8',header=0, names=['Gene','Average_replication_time'],usecols=[1,0])).dropna() 
    csv_file.close()
    in_mutvar_genes = [x for x in rep_time_data['Gene'].tolist()]
    in_mutvar_rates = [x for x in rep_time_data['Average_replication_time'].tolist()]
elif rep_timings=='No':
    in_mutvar_genes=[]
    in_mutvar_rates=[]

#Empirical gene mutation rate overrides all other rate calculations at end.
#Updates, so don't need all genes
if emp_mut_vars=='Yes':
    in_mutvar_genes = [x for x in empirical_weight_data['Gene'].tolist()]
    in_mutvar_rates = [x for x in empirical_weight_data['Average_replication_time'].tolist()]
#in_mutvar_genes=[]
#in_mutvar_rates=[]


if not os.path.exists(out_folder):
    try:
        os.makedirs(out_folder)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
          
#Collect identities and information for genes of interest

with open(observed1, 'r') as obs:
    ob_data1=(pandas.read_csv(obs, encoding='utf-8',header = 0, names=['gene','freq'],usecols=[1,0])).dropna()
    #ob_data1=(pandas.read_csv(obs, header = 0, names=['gene','freq'],usecols=[1,0])).dropna()
    obs.close()              
#ob_genes1 = [((x.encode('utf-8')).replace('\xe2\x80\x91', '-')) for x in ob_data1['gene'].tolist()]
ob_genes1 = [x for x in ob_data1['gene'].tolist()]
exp_freqs=ob_data1['freq'].tolist()
num_muts=sum(exp_freqs)
ob_gene_freq_dic1 = dict(zip(ob_genes1,exp_freqs))

#ob_genes=[]

records = [rec for rec in SeqIO.parse(ref_seq, "genbank")]

#aa redundancy table (generic)

aas='FLIMVSPTAYHQNKDECWRG'
NSvals='84796466688888888946'
Svals='15203533311111111053'

lengths=[]
true_lengths=[]
genes=[]
genes_for_display=[]
gois=[]
lois=[]
true_lois=[]
ordered_exp_freqs=[]
mods=[]
products=[]
ordered_rep_times=[]
ob_rep_times=[]

i=0
for record in records:
    for feature in record.features:
        if feature.type=="CDS":
            try:
                gene_name=str(feature.qualifiers['gene'][0])
                product=str(feature.qualifiers['product'][0])
            except KeyError:
                gene_name=str(feature.qualifiers['locus_tag'][0])
            if feature.qualifiers.get('translation') is not None:  
                start=int(feature.location.start)
                end=int(feature.location.end)
                length= abs(start-end)
                if use_length=='No':
                    length=1
                modifier=1
                if use_only_mutable=='Yes':
                    aa_seq=str(feature.qualifiers['translation'][0])
                    if syn_data=='No':
                        table = aa_seq.maketrans(aas,NSvals)
                    if syn_data=='Yes':
                        table = aa_seq.maketrans(aas,Svals)
                    trans_aa_seq=aa_seq.translate(table)
                    NS_total=0
                    for aa in trans_aa_seq:
                        NS_total+=int(aa)
                    all_muts=float(3)*float(length)
                    modifier=float(NS_total)/float(all_muts)
                #print(modifier)
                genes.append(gene_name)
                underscore_gene_name=gene_name.replace('-','_')
                #print(gene_name)
                rep_time=""
                for in_mod_gene, in_mod_rate in zip(in_mutvar_genes, in_mutvar_rates):
                    if gene_name == in_mod_gene or underscore_gene_name == in_mod_gene:
                        modifier=float(modifier)*float(in_mod_rate)
                        rep_time=in_mod_rate
                ordered_rep_times.append(rep_time)
                mod_length=float(modifier)*float(length)
                #print(length)
                #print(mod_length)
                lengths.append(mod_length)
                true_lengths.append(length)
                mods.append(modifier)
            else:
                print("CDS not translated")
            for ob_gene, exp_freq in zip(ob_genes1, exp_freqs):
                if gene_name == ob_gene or underscore_gene_name == ob_gene:
                    ordered_exp_freqs.append(exp_freq)
                    gois.append(gene_name)
                    true_lois.append(length)
                    lois.append(mod_length)
                    ob_rep_times.append(rep_time)
            gene_name=""
            underscore_gene_name=""
            start=""
            end=""
            length=""
            mod_length=""
            all_muts=""
                  


with open(str(out_folder)+"/Properties_of_observed_genes.csv","w") as results:
    wtr = csv.writer(results)
    table_header=["Gene","Gene length"]
    row_data=[]
    for g, l, m, r in zip(gois,true_lois,lois,ob_rep_times):
        row=[]
        row.extend([g,l])
        if emp_mut_vars=='No':
            if use_only_mutable=='Yes':
                row.append(m*3)
            if rep_timings=='Yes':
                row.append(r)
        row_data.append(row)
    if emp_mut_vars=='No':
        if use_only_mutable=='Yes':
            table_header.append("#NS bp changes possible")
        if rep_timings=='Yes':
            table_header.append("Replication timing")
    wtr.writerow(table_header)
    for row in row_data:
        wtr.writerow(row)
    results.close()

with open(str(out_folder)+"/Properties_of_all_genes.csv","w") as results:
    wtr = csv.writer(results)
    table_header=["Gene","Gene length"]
    row_data=[]
    for g, l, m, r in zip(genes,true_lengths,lengths,ordered_rep_times):
        row=[]
        row.extend([g,l])
        if emp_mut_vars=='No':
            if use_only_mutable=='Yes':
                row.append(m*3)
            if rep_timings=='Yes':
                row.append(r)
        row_data.append(row)
    if emp_mut_vars=='No':
        if use_only_mutable=='Yes':
            table_header.append("#NS bp changes possible")
        if rep_timings=='Yes':
            table_header.append("Replication timing")
    wtr.writerow(table_header)
    for row in row_data:
        wtr.writerow(row)
    results.close()
             

#Output file to show overall frequency distribution of gene lengths     
gene_freq_dic = {}
for k,v in Counter(true_lengths).items():
                gene_freq_dic[k]=v
with open(out_folder+"/gene_length_distribution.csv","w") as results2:
        wtr = csv.writer(results2)
        wtr.writerow(("Gene Size", "Frequency"))
        for key in sorted(gene_freq_dic.keys()):
                wtr.writerow((key,gene_freq_dic[key]))

if emp_mut_vars=='No':
    total_genes = len(genes)
    weights=[float(l)/float(sum(lengths)) for l in lengths]
elif emp_mut_vars=='Yes':
    genes=in_mutvar_genes
    weights=in_mutvar_rates     
else:
    print("Error: Indicate whether empirical mutation rates are being used: Yes or No.")

#Choosing given number of target genes by probability density fucntion defined by gene size
#This is repeated for the selected number of simulations

def mut_simulator(genes, num_muts, weights):
    end_muts = {}
    sim_muts=choice(genes, size=int(num_muts), replace =True, p=weights)
    for mut in sim_muts:
        if mut in end_muts.keys():
            end_muts[(mut)]+=1
        else:
            end_muts[str(mut)]=1
    return end_muts
sims=[mut_simulator(genes, num_muts, weights) for i in range(int(num_reps))]


#Update headers for null data
current_max=0
for sim in sims:
    for k,v in sim.items():
        if v>current_max:
            current_max=v
null_headers=["Gene", "Times never hit"]
for i in range(1,current_max+1):
    if i==1:
        a=str("Times hit in 1 population")
    else:
        a=str("Times hit in " +str(i) +" populationss")
    null_headers.append(a)
null_headers.append('Maximum null value')
null_headers.append('Experimental hits')
with open(str(out_folder)+"/All-Genes_simulation_genes-"+str(num_muts)+"_reps-"+str(int(num_reps))+"_tallies.csv", 'w') as results3:
    with open(str(out_folder)+"/Observed-Genes_simulation_genes-"+str(num_muts)+"_reps-"+str(int(num_reps))+"_tallies.csv", 'w') as results4:
        wtr=csv.writer(results3)
        wtr.writerow(null_headers)
        wtr2=csv.writer(results4)
        wtr2.writerow(null_headers)
        for gene in genes:
            row_data=[gene]
            sub_zero=0
            n_collect=[]
            mut_sum_dic = {}
            mut_sum = []
            mut_counter=[]
            gene_max=0
            for i in range(1,(current_max+1)):
                n=0
                for dic in sims:
                    if gene in dic:
                       if dic[gene]==i:
                           n+=1
                n_collect.append(n)
                if n>0:
                    gene_max=i
                sub_zero+=n
            row_data.append(int(num_reps-sub_zero))
            row_data.extend(n_collect)
            row_data.append(gene_max)
            if gene not in gois:
                row_data.append(0)
            else:
                for goi,freq in zip(gois,ordered_exp_freqs):
                    if gene == goi:
                        row_data.append(freq)
                        wtr2.writerow(row_data)
            wtr.writerow(row_data)
    results4.close()
results3.close()
with open(str(out_folder)+"/All-Genes_simulation_genes-"+str(num_muts)+"_reps-"+str(int(num_reps))+"_multis_by_run.csv", 'w') as results5:
    wtr=csv.writer(results5)
    wtr.writerow(("Run_number", "Gene", "Times_hit"))
    i=0
    for sim in sims:
        i+=1
        for k, v in sim.items():
            if v >=2:
                wtr.writerow((i,k,v))
    results5.close()
#return 
#Collecting count stats. Add more count categories if desired
#Here, multihit genes are not expected to exceed 5 hits
CI_z=st.norm.ppf((1-(float(1)/2)/(num_reps)))
#Count_0 is a list of number of genes of total not picked in a dic
count_0 = [(len(genes)-len(dic)) for dic in sims]
counter_list=[count_0]
for i in range(1,(current_max+1)):
    i_counter=[list(dic.values()).count(i) for dic in sims]
    counter_list.append(i_counter)
max_picks = []
min_picks = []
averages = []
pos_CIs = []
neg_CIs = []
y_err_bars = []
for count_list in counter_list:
    max_picks.append(max(count_list))
    min_picks.append(min(count_list))
    averages.append(average(count_list))
    #y_err_bar = (st.norm.ppf(1-float(1)/num_reps))*std(count_list)/(float(num_reps))**0.5
    y_err_bar = CI_z*(std(count_list,ddof=num_reps-1))/(float(num_reps)**0.5)
    y_err_bars.append(y_err_bar)
    pos_CIs.append(average(count_list)+y_err_bar)
    neg_CIs.append(average(count_list)-y_err_bar)
with open(str(out_folder)+"/STATSb_"+str(num_muts)+"_reps-"+str(int(num_reps))+".csv", 'w') as out:
    stats = csv.writer(out)
    stats.writerow(("Times picked","Max","Min","Average","CI"))
    for i in range(len(counter_list)):
        stats.writerow((i,max_picks[i],min_picks[i],averages[i],str(neg_CIs[i])+", "+str(pos_CIs[i])))
    out.close()
 