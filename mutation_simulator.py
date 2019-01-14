import csv
import re
from collections import Counter
from numpy.random import choice
import pandas
#import sympy
import scipy.stats as st
from numpy import std, average, unique
import matplotlib.pyplot as plt
import matplotlib as mpl

print "/Users/lwoo0005/Documents/Laura_stuff/Misc Bioinformatics/w303_ref2.gff"
ref_seq = str(raw_input("Input the path to the genome reference file.\n"))
num_muts = int(raw_input("How many mutations?\n"))
num_reps = float(raw_input("How many simulations should be run?\n"))
print "\n/Users/lwoo0005/Documents/Laura_stuff/Misc Bioinformatics/final_para_list_nodoubles_LW_wt.csv"
observed1 = str(raw_input("Please input the path to the first observation CSV file in which genes are in row 1 and frequencies in row 2.\n"))
observed2 = str(raw_input("Please input the path to the second observation CSV file in which genes are in row 1 and frequencies in row 2.\n"))
print "\n/Users/lwoo0005/Documents/Laura_stuff/ExpEvol_Program_Tests/16-11-17"
out_folder = str(raw_input("Where do you want results to go? Enter a path.\n"))

with open(observed1, 'r') as obs:
        ob_data1=(pandas.read_csv(obs, sep=',',header = 0, names=['freq','gene'],usecols=[0,1])).dropna()
        obs.close()
with open(observed2, 'r') as obs2:
        ob_data2=(pandas.read_csv(obs2, sep=',',header = 0, names=['freq','gene'],usecols=[0,1])).dropna()
        obs2.close()        
ob_genes1 = ob_data1['gene'].tolist()
ob_genes2 = ob_data2['gene'].tolist()
ob_genes_all = ob_genes1+ob_genes2
ob_genes = unique(ob_genes_all).tolist()
ob_nums1 = ob_data1['freq'].tolist()
ob_nums2 = ob_data2['freq'].tolist()
ob_freqs1, ob_freqs2 = [int(i) for i in ob_nums1], [int(i) for i in ob_nums2]

ob_gene_freq_dic1, ob_gene_freq_dic2 = dict(zip(ob_genes1,ob_freqs1)), dict(zip(ob_genes2,ob_freqs2))

lengths=[]
gene_len_dic = {}
gene_freq_dic = {}
test = []

with open(ref_seq ,"r") as gff_file:
        reader = csv.reader(gff_file, delimiter="\t")
        next(reader, None)
        with open(str(out_folder)+"/gene_lengths.csv","w") as results:
                wtr = csv.writer(results)
                wtr.writerow(("Gene","Start Position", "End Position", "Gene Length"))
                for row in reader:
                        #limit type of genes here if required
                        if str(row[2]) == "gene":
                                start_pos = int(row[3])
                                end_pos = int(row[4])
                                length = abs(start_pos-end_pos)
                                lengths.append(length)
                                gene_name = re.sub(r'%..', r' ', row[8])
                                gene_name = "".join(gene_name.split("Note="))
                                for ob_gene in ob_genes:
                                        search_gene = "gene="+str(ob_gene)+";"
                                        if search_gene in gene_name:
                                                test.append(ob_gene)
                                                gene_name=ob_gene
                                wtr.writerow((gene_name,start_pos,end_pos,length))
                                gene_len_dic[gene_name] = float(length)

for k,v in Counter(lengths).iteritems():
                gene_freq_dic[k]=v
                
with open(out_folder+"/gene_length_distribution.csv","w") as results2:
        wtr = csv.writer(results2)
        wtr.writerow(("Gene Size", "Frequency"))
        for key in sorted(gene_freq_dic.iterkeys()):
                wtr.writerow((key,gene_freq_dic[key]))

weights = []
#weights_normal = []
sort_gene_names = []
sort_lengths = []

total_genes = len(gene_len_dic)
for key, value in sorted(gene_len_dic.iteritems(), key=lambda (k,v): (v,k)):
                sort_gene_names.append(key)
                sort_lengths.append(value)
                
for length in sort_lengths:
    weight = float(length)/float(sum(sort_lengths))
    weights.append(weight)

def many_sims(num_reps,sort_gene_names,num_muts,weights):
    all_reps = []
    def mut_simulator(sort_gene_names, num_muts, weights):
        end_muts = {}
        sim_muts=choice(sort_gene_names, size=num_muts, replace =True, p=weights)
        for mut in sim_muts:
            if str(mut) in end_muts.keys():
                end_muts[str(mut)]+=1
            else:
                end_muts[str(mut)]=1
        return end_muts
    for i in range(int(num_reps)):
        sim=mut_simulator(sort_gene_names, num_muts, weights)
        all_reps.append(sim)
    return all_reps

sims=many_sims(num_reps,sort_gene_names,num_muts,weights)
"""
#old code
    def all_mut_sim(sort_gene_names, genes_of_interest, num_muts, weights, num_reps, file_name):

    all_reps = []
    def mut_simulator(sort_gene_names, num_muts, weights_normal):
        end_muts = {}
        sim_muts=choice(sort_gene_names, size=num_muts, replace =True, p=weights_normal)
        for mut in sim_muts:
            if str(mut) in end_muts.keys():
                end_muts[str(mut)]+=1
            else:
                end_muts[str(mut)]=1
        return end_muts
    for i in range(int(num_reps)):
        sim=mut_simulator(sort_gene_names, num_muts, weights)
        all_reps.append(sim)
"""
def sim_writes(sim_dic_list,sort_gene_names,genes_of_interest,file_name,num_muts,num_reps) :
    with open(str(out_folder)+"/"+file_name+"_simulation_genes-"+str(num_muts)+"_reps-"+str(int(num_reps))+".csv", 'w') as results3:
        wtr=csv.writer(results3)
        wtr.writerow(("Gene","Expected value for one simulation given "+str(num_reps)+" trials", "Never picked in "+str(num_reps)+" trials", "Picked once", "Picked twice","Picked three times", "Picked four times", "Picked five times", "Picked six times or more", "CI_lower", "CI_upper"))
        for gene in sort_gene_names:
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
                CI_int = 1.96*std(mut_sum)/((float(num_reps))**0.5)
                avg = average(mut_sum)
                mut_sum_dic[gene]=[avg,avg+CI_int,avg-CI_int]
                for k,v in mut_sum_dic.items():
                    wtr.writerow((k,v[0],not_pick,pick_one,pick_two,pick_three,pick_four,pick_five, pick_six_plus, v[2], v[1]))
        results3.close()
    return 

all_counts=sim_writes(sims,sort_gene_names, sort_gene_names, "ALL-GENES",num_muts, num_reps)
null_counts=sim_writes(sims,sort_gene_names, ob_genes, "OBSERVED-GENES",num_muts, num_reps)


#Collecting count stats
CI_z=st.norm.ppf((1-(float(1)/2)/(num_reps)))
count_0 = [(len(sort_gene_names)-len(dic)) for dic in sims]
count_1=[dic.values().count(1) for dic in sims]
count_2=[dic.values().count(2) for dic in sims]
count_3=[dic.values().count(3) for dic in sims]
count_4=[dic.values().count(4) for dic in sims]
count_5=[dic.values().count(5) for dic in sims]
averages = []
pos_CIs = []
neg_CIs = []
y_err_bars = []
for count_list in [count_0, count_1,count_2,count_3,count_4,count_5]:
    averages.append(average(count_list))
    #y_err_bar = (st.norm.ppf(1-float(1)/num_reps))*std(count_list)/(float(num_reps))**0.5
    y_err_bar = CI_z*(std(count_list,ddof=num_reps-1))/(float(num_reps)**0.5)
    y_err_bars.append(y_err_bar)
    pos_CIs.append(average(count_list)+y_err_bar)
    neg_CIs.append(average(count_list)-y_err_bar)
with open(str(out_folder)+"/STATS_"+str(num_muts)+"_reps-"+str(int(num_reps))+".csv", 'w') as out:
    stats = csv.writer(out)
    stats.writerow(("Times picked","Average","CI"))
    for i in range(0,6):
        stats.writerow((i,averages[i],str(neg_CIs[i])+","+str(pos_CIs[i])))
    out.close()
    
  
#Plotting simulations by count
cmap = mpl.cm.YlOrBr
mpl_color = plt.figure()
c=0

for dic in sims:
    c+=1
    #x=["0","1","2","3","4","5"]
    x=["1","2","3","4","5"]
    y=[]
    y.append(len(sort_gene_names)-len(dic))
    for i in range(1,6):
        if i in dic.values():
            y_values=(dic.values().count(i))
        else:
            y_values=0
        y.append(y_values)
    #plt.plot(x,y,lw=0.1,color=cmap(float(c)/2500))
    plt.plot(x,y,lw=0.1,color="orange")

#Plotting average of simulations by count
#plt.plot(["0","1","2","3","4","5"],averages, lw=0.5, color="red")
#plt.fill_between(["0","1","2","3","4","5"], pos_CIs, neg_CIs, facecolor='yellow', alpha=0.5)
#plt.errorbar(["0","1","2","3","4","5"],averages,yerr=y_err_bars,color="red")    

plt.plot(["1","2","3","4","5"],averages, lw=0.5, color="red")
plt.fill_between(["1","2","3","4","5"], pos_CIs, neg_CIs, facecolor='yellow', alpha=0.5)
plt.errorbar(["1","2","3","4","5"],averages,yerr=y_err_bars,color="red")    
 
  
#Plotting observed data by count
ob_colours = ["blue","cyan"]
c=0
for dic in [ob_gene_freq_dic1,ob_gene_freq_dic2]:
    #x=["0","1","2","3","4","5"]
    x=["1","2","3","4","5"]
    y=[]
    #y.append(len(sort_gene_names)-len(dic))
    for i in range(1,6):
        if i in dic.values():
            y_values=(dic.values().count(i))
        else:
            y_values=0
        y.append(y_values)
    plt.plot(x,y,lw=0.8,color=ob_colours[c])
    c+=1
plt.show()



#plotting sims by gene

c=0
cmap = mpl.cm.YlOrBr
mpl_color = plt.figure()

for dic in sims:
    c+=1
    x=[]
    y=[]
    for gene in ob_genes:
        for k,v in dic.items():
            if k == gene:
                x.append(k)
                y.append(v)
        if gene not in dic.keys():
            x.append(gene)
            y.append(0)
    plot_dic=dict(zip(x,y))
    x=[]
    y=[]
    for key in sorted(plot_dic.iterkeys()):
        x.append(key)
        y.append(plot_dic[key])
    plt.plot(x,y,color=cmap(float(c)/2000),lw=0.2)
    
#Plotting observed genes by gene
ob_colours = ["blue","cyan"]
c=0
for dic in [ob_gene_freq_dic1,ob_gene_freq_dic2]:
    x=[]
    y=[]
    for gene in ob_genes:
        for k,v in dic.items():
            if k == gene:
                x.append(k)
                y.append(v)
        if gene not in dic.keys():
            x.append(gene)
            y.append(0)
    plot_dic=dict(zip(x,y))
    x=[]
    y=[]
    for key in sorted(plot_dic.iterkeys()):
        x.append(key)
        y.append(plot_dic[key])
    plt.plot(x,y,color=ob_colours[c],lw=0.5)
    c+=1
plt.ylabel('Counts')
plt.xlabel('Genes')
plt.xticks(x, rotation='vertical')
plt.ylim(ymin=0)
plt.tight_layout()
plt.show()  

#14-11-18 I don't think this bit was used so I've triple-quoted it out.
"""


    #max likelihood methods
    
    print freq
    print gene
    pos_int=float(freq)+(1.96*(float(freq)))**0.5
    print pos_int
    neg_int=float(freq)-(1.96*(float(freq)))**0.5
    print neg_int

    CI = (1.96/num_muts)*((freq*(num_muts-freq))/num_muts)**0.5
    pos_int = ob_prob+CI
    neg_int = ob_prob-CI
    print gene
    print freq
    print ob_prob
    print pos_int
    print neg_int

#re-weighting including unobserved genes and re-run of sim
pos_indices = []
neg_indices = []
for sorted_gene in sort_gene_names:
    if sorted_gene in ob_genes:
        pos_indices.append(int(sort_gene_names.index(sorted_gene)))
    else:
        neg_indices.append(int(sort_gene_names.index(sorted_gene)))
     
i = 0
factors = []
for index in pos_indices:
    start_weight=float(weights[index])
    weights[index]=float(ob_freqs[i])/sum(ob_freqs)
    factors.append(weights[index]/start_weight)
    i+=1
    
divisor=float(sum(factors))
for index in neg_indices:
    start_weight=float(weights[index])
    weights[index]=float(start_weight/divisor)

for weight in weights:
    weight = weight/float(sum(weights))
    weights_normal.append(weight)

ob_counts=all_mut_sim(sort_gene_names, ob_genes, num_muts, weights_normal, num_reps, "OBSERVED")

#reweighting with only observed genes with max likelihood-derived independent event probabilities
new_weights = [i/num_muts for i in ob_freqs]
null_for_ob_genes = all_mut_sim(ob_genes, ob_genes, num_muts, new_weights, num_reps, "RE-WEIGHTED_OBSERVED")

#solving for 0.95 = (nCx)(p^x)((1-p)^n-x)    
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction

def nCk(n,k): 
  return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )

for freq in ob_freqs:
    n_C_k = nCk(int(num_muts),int(freq)
    y = 0.95/n_C_k
    solve(x**)

"""

print "Analysis complete and located in output folder--Thanks!"
         

        
