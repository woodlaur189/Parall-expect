#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 20:41:26 2019

@author: lwoo0005
"""

#Import relevant modules
import csv
import re
from collections import Counter
from numpy.random import choice
import pandas
import os, errno
import scipy.stats as st
from numpy import std, average, unique

#Input and output definitions
print "/Users/lwoo0005/Documents/Laura_stuff/Misc Bioinformatics/w303_ref2.gff"
ref_seq = str(raw_input("Input the path to the genome reference file.\n"))
num_muts = int(raw_input("How many mutations?\n"))
num_reps = float(raw_input("How many simulations should be run?\n"))
print "\n/Users/lwoo0005/Documents/Laura_stuff/Misc Bioinformatics/final_para_list_nodoubles_LW_wt.csv"
observed1 = str(raw_input("Please input the path to the first observation CSV file in which genes are in row 1 and frequencies in row 2.\n"))
observed2 = str(raw_input("Please input the path to the second observation CSV file in which genes are in row 1 and frequencies in row 2.\n"))
print "\n/Users/lwoo0005/Documents/Laura_stuff/ExpEvol_Program_Tests/16-11-17"
out_folder = str(raw_input("Where do you want results to go? Enter a path.\n"))

#Collect identities and information for genes of interest
with open(observed1, 'r') as obs:
    ob_data1=(pandas.read_csv(obs, sep=',',header = 0, names=['freq','gene'],usecols=[0,1])).dropna()
    obs.close()       
with open(observed2, 'r') as obs2:
    ob_data2=(pandas.read_csv(obs2, sep=',',header = 0, names=['freq','gene'],usecols=[0,1])).dropna()
    obs2.close()        
ob_genes1 = ob_data1['gene'].tolist()
ob_genes2 = ob_data2['gene'].tolist()
ob_genes = unique(ob_genes1+ob_genes2)
ob_gene_freq_dic1, ob_gene_freq_dic2 = dict(zip(ob_genes1,ob_data1['freq'].tolist())), dict(zip(ob_genes2,ob_data2['freq'].tolist()))

#Open GFF annotation file for species of interest and collect gene name and length information
lengths=[]
genes=[]
with open(ref_seq ,"r") as gff_file:
        reader = csv.reader(gff_file, delimiter="\t")
        #Skip header
        next(reader, None)
        #Make new directory if required
        if not os.path.exists(out_folder):
            try:
                os.makedirs(out_folder)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise
        with open(str(out_folder)+"/lengths_of_observed_genes.csv","w") as results:
                wtr = csv.writer(results)
                wtr.writerow(("Gene","Start Position", "End Position", "Gene Length"))
                for row in reader:
                        #limit type of genes here if required
                        if str(row[2]) == "gene":
                                start_pos = int(row[3])
                                end_pos = int(row[4])
                                length = abs(start_pos-end_pos)
                                lengths.append(length)
                                #Specific to GFF format; modify for alternatives
                                gene_name = re.sub(r'%..', r' ', row[8])
                                genes.append("".join(gene_name.split("Note=")))
                                for ob_gene in ob_genes:
                                        if "gene="+str(ob_gene)+";" in gene_name:
                                            wtr.writerow((ob_gene,start_pos,end_pos,length))
        gff_file.close()
        
#Output file to show overall frequency distribution of gene lengths     
gene_freq_dic = {}
for k,v in Counter(lengths).iteritems():
                gene_freq_dic[k]=v
with open(out_folder+"/gene_length_distribution.csv","w") as results2:
        wtr = csv.writer(results2)
        wtr.writerow(("Gene Size", "Frequency"))
        for key in sorted(gene_freq_dic.iterkeys()):
                wtr.writerow((key,gene_freq_dic[key]))

total_genes = len(genes)
weights=[float(l)/float(sum(lengths)) for l in lengths]

#Choosing given number of target genes by probability density fucntion defined by gene size
#This is repeated for the selected number of simulations
def many_sims(num_reps,sort_gene_names,num_muts,weights):
    def mut_simulator(sort_gene_names, num_muts, weights):
        end_muts = {}
        sim_muts=choice(sort_gene_names, size=num_muts, replace =True, p=weights)
        for mut in sim_muts:
            if str(mut) in end_muts.keys():
                end_muts[str(mut)]+=1
            else:
                end_muts[str(mut)]=1
        return end_muts
    all_reps=[mut_simulator(sort_gene_names, num_muts, weights) for i in range(int(num_reps))]
    return all_reps

#Applying function with selected arguments
sims=many_sims(num_reps,genes,num_muts,weights)

#Writing simulation results to CSV
def sim_writes(sim_dic_list,gene_names,genes_of_interest,file_name,num_muts,num_reps) :
    with open(str(out_folder)+"/"+file_name+"_simulation_genes-"+str(num_muts)+"_reps-"+str(int(num_reps))+".csv", 'w') as results3:
        wtr=csv.writer(results3)
        wtr.writerow(("Gene","Expected value for one simulation given "+str(num_reps)+" trials", "Never picked in "+str(num_reps)+" trials", "Picked once", "Picked twice","Picked three times", "Picked four times", "Picked five times", "Picked six times or more", "CI_lower", "CI_upper"))
        for gene in gene_names:
            if gene in genes_of_interest:
                mut_sum_dic = {}
                mut_sum = []
                not_pick = 0
                pick_one = 0
                pick_two = 0
                pick_three = 0
                pick_four = 0
                pick_five = 0
                pick_six_plus = 0
                for dic in sim_dic_list:
                    if gene in dic:
                        mut_sum.append(float(dic[gene]))
                        if dic[gene] ==1:
                            pick_one += 1
                        if dic[gene] ==2:
                            pick_two += 1
                        if dic[gene] ==3:
                            pick_three += 1
                        if dic[gene] == 4:
                            pick_four +=1
                        if dic[gene] == 5:
                            pick_five += 1
                        if dic[gene] >= 6:
                            pick_six_plus += 1
                    else:
                        not_pick += 1
                        mut_sum.append(0)
                #Alter confidence interval as required
                #Take CI as input in future versions
                CI_int = 1.96*std(mut_sum)/((float(num_reps))**0.5)
                avg = average(mut_sum)
                mut_sum_dic[gene]=[avg,avg+CI_int,avg-CI_int]
                for k,v in mut_sum_dic.items():
                    wtr.writerow((k,v[0],not_pick,pick_one,pick_two,pick_three,pick_four,pick_five, pick_six_plus, v[2], v[1]))
        results3.close()
    return 

#Output files are produced for simulations
#The first file displays the results for all genes
all_counts=sim_writes(sims,genes, genes, "ALL-GENES",num_muts, num_reps)
#The second file displays the results only for genes of interest
#Note that this modulates information written only so results are the same in both files.
null_counts=sim_writes(sims,genes, ob_genes, "OBSERVED-GENES",num_muts, num_reps)

#Collecting count stats. Add more count categories if desires
#Here, gmultihit genes are not expected to exceed 5 hits
CI_z=st.norm.ppf((1-(float(1)/2)/(num_reps)))
#Count_0 is a list of number of genes of total not picked in a dic
count_0 = [(len(genes)-len(dic)) for dic in sims]
count_1=[dic.values().count(1) for dic in sims]
count_2=[dic.values().count(2) for dic in sims]
count_3=[dic.values().count(3) for dic in sims]
count_4=[dic.values().count(4) for dic in sims]
count_5=[dic.values().count(5) for dic in sims]
max_picks = []
min_picks = []
averages = []
pos_CIs = []
neg_CIs = []
y_err_bars = []
for count_list in [count_0, count_1,count_2,count_3,count_4,count_5]:
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
    for i in range(0,6):
        stats.writerow((i,max_picks[i],min_picks[i],averages[i],str(neg_CIs[i])+", "+str(pos_CIs[i])))
    out.close()