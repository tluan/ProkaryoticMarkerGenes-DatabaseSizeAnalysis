import sys
import os
import random
import numpy as np
from Bio import SeqIO
def cluster_analysis(clstr,cog_,pid_,replicate_): # parses the cd-hit *.clstr file to determine number of multispecies clusters
	r = {}
	geneLens = []
	species = {}
	single_species = 0
	multispecies = 0
	
	for i in open(clstr):
		if i[0] == ">":
			cluster = i.strip().replace(">","")
			r[cluster] = []
		else:
			tmp = i.strip().split('\t')[1]
			geneLens.append(int(tmp.split(',')[0].replace('aa','')))
			taxa = tmp.split(">")[1].split('...')[0]
			#print(taxa)
			species[taxa] = 0
			r[cluster].append(taxa)
	
	for k,v in r.items():
		num_species = len(set(v))
		if num_species == 1: single_species += 1
		elif num_species > 1: multispecies += 1
	
	total_species = len(species)
	num_genes = len(geneLens)
	gene_mean_len = np.mean(geneLens)
	gene_std_len = np.std(geneLens)
	num_clusters = len(r)
	print("total_species")
	print(total_species)
	print(str(replicate)+'\t'+str(num_genes)+'\t'+str(total_species)+'\t'+str(num_clusters)+'\t'+str(single_species)+'\t'+str(multispecies)+'\t'+str(gene_mean_len)+'\t'+str(gene_std_len)+'\n')
	with open(cog_+'.'+pid_+'_results.txt','a') as out3:
		out3.write(str(replicate_)+'\t'+str(num_genes)+'\t'+str(total_species)+'\t'+str(num_clusters)+'\t'+str(single_species)+'\t'+str(multispecies)+'\t'+str(gene_mean_len)+'\t'+str(gene_std_len)+'\n')

fasta = sys.argv[1]
cog = fasta.replace('.fna.taxa','')
pid = sys.argv[2] # for cd-hit, 0.99 is equal to 99% identity

seqs = [">%s\n%s\n" %(str(i.description),str(i.seq)) for i in SeqIO.parse(fasta,'fasta')]
L = len(seqs)

with open(cog+'.'+pid+'_results.txt','w') as out:
	out.write('Replicate\tNum_genes\tNum_species\tNum_clusters\tSingle_species\tMultispecies\tGene_mean_len\tGene_std_len\n')

with open(cog+'.'+pid+'_results.txt','w') as out2:
	for replicate in range(30): # this loop is number of bootstraps
		random.shuffle(seqs) # for each replicate of experiment, order of sequences is shuffled
		database = []
		for h,i in enumerate(seqs):
			if h == 0:
				database.append(i)
			elif h % 10000 == 0: # clustering and data collection will occur for every 100 sequences sampled
				database.append(i)
				with open(cog+'_tmp_'+pid,'w') as out:
					for seq in database:
						out.write(seq)
				os.system('/cbcbhomes/tluan/cdhit/cd-hit -i %s_tmp_%s -o %s_tmp.p%sq100 -c %s -g 1 -G 0 -aS 1.0 -d 50 -T 16' % (cog,pid,cog,pid,pid))
				cluster_analysis(cog+'_tmp.p%sq100.clstr' % (pid),cog,pid,replicate)
			else:
				database.append(i)
