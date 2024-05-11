import sys
import os
from collections import defaultdict
from collections import Counter
import random
import numpy as np
from Bio import SeqIO

def cluster_analysis(clstr): # parses the cd-hit *.clstr file to determine number of multispecies clusters
	r = defaultdict(list)
	
	for i in open(clstr):
		if i[0] == ">":
			cluster = int(i.strip().split(' ')[1])
		else:
			tmp = i.strip().split('\t')[1]
			taxa = '_'.join(tmp.split(">")[1].split('...')[0].split('_')[:-3])
			r[cluster].append(taxa)
	
	with open(cog+'.'+pid+'_results.txt','a') as out:
		for k,v in sorted(r.items()):
			species,counts = [],[]
			for z,w in sorted(Counter(v).items(),reverse=True,key=lambda x:x[1]):
				species.append(z)
				counts.append(str(w))
			out.write('%s\t%s\t%s\t%s\t%s\n' % (replicate,str(h),str(k),','.join(species),','.join(counts)))

fasta = sys.argv[1]
cog = fasta.replace('.fna.taxa','')
pid = sys.argv[2] # for cd-hit, 0.99 is equal to 99% identity

seqs = [">%s\n%s\n" %(str(i.description),str(i.seq)) for i in SeqIO.parse(fasta,'fasta')]
L = len(seqs)

with open(cog+'.'+pid+'_results.txt','w') as out:
	out.write('Replicate\tNum_genes\tCluster_number\tSpecies\tSpecies_counts\n')

for replicate in range(100): # this loop is number of bootstraps
	random.shuffle(seqs) # for each replicate of experiment, order of sequences is shuffled
	database = []
	for h,i in enumerate(seqs):
		if h == 0:
			database.append(i)
		elif h % 1000 == 0: # clustering and data collection will occur for every 1000 sequences sampled
			database.append(i)
			with open(cog+'_tmp_'+pid,'w') as out:
				for seq in database:
					out.write(seq)
			os.system('cd-hit -T 12 -i %s_tmp_%s -o %s_tmp.p%sq100 -d 1000 -c %s -g 1 -G 0 -aS 1.0' % (cog,pid,cog,pid,pid))
			cluster_analysis(cog+'_tmp.p%sq100.clstr' % (pid))
		else:
			database.append(i)
