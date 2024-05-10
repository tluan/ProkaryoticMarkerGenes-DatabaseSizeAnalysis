from Bio import SeqIO

# nodes.dmp and fullnamelineage.dmp are files downloaded from the FTP site for the NCBI Taxonomy database

species_nodes = {i.split('|')[0].strip() for i in open('nodes.dmp') if i.split('|')[2].strip() == 'species'}

for i in open('fullnamelineage.dmp'):
	tmp = i.split('|')
	id,lineage = tmp[0].strip(),tmp[2].strip().lower()
	if id in species_nodes:
		if 'unclassified' in lineage:
			species_nodes.remove(id)

names = {i.strip().lower().split('|')[1].strip():i.split('|')[0].strip() for i in open('names.dmp') if i.split('|')[0].strip() in species_nodes}

print('Step A')

with open('SILVA_138.1_SSURef_tax_silva_DNA.fasta','w') as out:
	for i in SeqIO.parse('SILVA_138.1_SSURef_tax_silva.fasta','fasta'):
		d,s = str(i.description),str(i.seq)
		lineage = d.lower()
		# filter out sequences with the following labels
		for j in ['uncultured','eukaryota','mitochondria','chloroplast','unidentified','symbiont','unclassified','unknown','metagenome']:
			if j in lineage: break
		else:
			tmp = lineage.split(';')[-1].strip().split(' ')
			species = ' '.join(tmp[:2])
			# filter out sequences that don't have a species label
			if species in names:
				speciesID = names[species]
				s = s.upper().replace('U','T')
				out.write(">%s\n%s\n" % (speciesID,s))
			else:
				for j in tmp[2:]:
					species += ' '+j
					if species.endswith('sp.'): continue
					elif species in names:
						speciesID = names[species]
						s = s.upper().replace('U','T')
						out.write(">%s\n%s\n" % (speciesID,s))
						break
