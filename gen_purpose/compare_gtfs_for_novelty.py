#!/usr/bin/env python2.7

def main():
	args = getArgs()
	ref, new = args.ref, args.new
	### store both gtfs
	### dict in form of gene, transcript, list of coords for the transcript
	### dict[gene][transcript] = [1:4, 6:8, 10:14]
	refGtf, newGtf = storeGtfs(ref = ref, new = new)
	#refSJ, newSJ = makeSJdict(ref = refGtf, new = newGtf)
	genesNovelIsos = compareGtfs(ref = refGtf, new = newGtf)
	#compareJuncs(ref = refGtf, new = newGtf)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = "compare transcripts between gtfs to count up number of novel and known isoforms. Known isoforms are those that are present in the reference gtf.")
	parser.add_argument('-ref', action = 'store', type = str, required = True, help = "reference gtf file")
	parser.add_argument('-new', action = 'store', type = str, required = True, help = "new gtf file")
	args = parser.parse_args()
	return(args)

def makeSJdict(ref, new):
	from collections import defaultdict
	refSJ, newSJ = defaultdict(dict), defaultdict(dict)
	fh = open(ref, 'r')
	for line in fh:
		line = line.strip()
		geneArr = arr[8].split(';')[0].split('=')[1].split('.')
		gene = geneArr[0] + '.' + geneArr[1]
		#print(line)
		tidArr = arr[8].split(';')[0].split('=')[1].split('.')
		tid = geneArr[0] + '.' + geneArr[1] + '.' + geneArr[2]
		if(gene in refSJ):
			if(left in refSJ[gene]):
				if(right not in refSJ[gene][left]):
					refSJ[gene][left].append(right)
			else:
				refSJ[gene][left] = [right]
		else:
			refSJ[gene][left] = [right]

	fh = open(new, 'r')
	for line in fh:
		line = line.strip()
		arr = line.split('\t')
		if(arr[2] == "transcript"):
			continue
		else:
			gene = line.split('\"')[3]
			left = arr[3]
			right = arr[4]
			tid = line.split('\"')[1]
			chromo = arr[0]
			if(gene in newSJ):
				if(left in newSJ[gene]):
					if(right not in newSJ[gene][left]):
						newSJ[gene][left].append(right)
				else:
					newSJ[gene][left] = [right]
			else:
				newSJ[gene][left] = [right]
	return(refSJ, newSJ)

def storeGtfs(ref, new):
	import io,re
	from collections import defaultdict
	refGtf, newGtf = defaultdict(dict), defaultdict(dict)
	ref, new = open(ref, 'r'), open(new, 'r')
	for line in ref:
		line = line.strip()
		if(re.search('exon', line)):
			arr = line.split()
			chromo = arr[0]
			### change these depending on the format of the gtf file
			#  Chr01   phytozomev10    exon    1660    2502    .       -       .       ID=Potri.001G000100.1.v3.0.exon.1;Parent=Potri.001G000100.1.v3.0;pacid=27043735
			geneArr = arr[8].split(';')[0].split('=')[1].split('.')
			gene = geneArr[0] + '.' + geneArr[1]
			#print(line)
			tidArr = arr[8].split(';')[0].split('=')[1].split('.')
			tid = geneArr[0] + '.' + geneArr[1] + '.' + geneArr[2]
			########
			five, three = arr[3], arr[4]
			coords = five + ':' + three
			if(gene in refGtf):
				if(tid in refGtf[gene]):
					refGtf[gene][tid].append(coords)
				else:
					refGtf[gene][tid] = [chromo, coords]
			else:
				refGtf[gene][tid] = [chromo, coords]

	for line in new:
		line = line.strip()
		if(re.search('exon', line)):
			arr = line.split('\t')
			### change these depending on the format of the gtf file
			chromo = arr[0]
			gene = arr[8].split(';')[1].split('\"')[1]
			tid = arr[8].split(';')[0].split('\"')[1]
			########
			five, three = arr[3], arr[4]
			coords = five + ':' + three
			if(gene in newGtf):
				if(tid in newGtf[gene]):
					newGtf[gene][tid].append(coords)
				else:
					newGtf[gene][tid] = [chromo, coords]
			else:
				newGtf[gene][tid] = [chromo, coords]
	return(refGtf, newGtf)

def compareGtfs(ref, new):
	novelIso, knownIso, detIso, notDetIso = 0,0,0,0
	genesNovelIsos = []
	novelIsoList = []
	for gene in new:
		#print(gene)
		#print(str(ref[gene]))
		for tid in new[gene]:
			check = 0
			
			for oldTid in ref[gene]:
				#print(gene + '\t' + tid + '\t' + str(new[gene][tid]))
				if(new[gene][tid] == ref[gene][oldTid]):
					check = 1
			if(check == 1):
				knownIso += 1
				#print(gene + '\t' + tid + '\t' + "known")
			else:
				novelIso += 1
				novelIsoList.append(tid)
				if gene not in genesNovelIsos:
					genesNovelIsos.append(gene)
				#print(gene + '\t' + tid + '\t' + "novel")
	print("master xome transcripts in reference:\t" + str(knownIso))
	print("novel transcripts not in reference:\t" + str(novelIso))
	#print('\n'.join(novelIsoList))	
	for gene in ref:
		for tid in ref[gene]:
			check = 0
			#if gene in new:
			for newTid in new[gene]:
				if(ref[gene][tid] == new[gene][newTid]):
					check = 1
			if(check == 1):
				detIso += 1
				#print(gene + '\t' + tid + '\tdetected' )
			else:
				notDetIso +=1
				#print(gene + '\t' + tid + '\t not detected' )
	#print("reference transcripts in new xome:\t" + str(detIso))
	#print("reference transcripts not in new xome:\t" + str(notDetIso))
	return(genesNovelIsos)

def compareJuncs(ref,new):
	from collections import defaultdict
	refSJ, newSJ = defaultdict(dict), defaultdict(dict)
	novelSJgenes = []
	for gene in ref:
		for iso in ref[gene]:
			arr = ref[gene][iso]
			chromo = arr[0]
			#print(str(arr))
			for i in range(1,(len(arr) - 1)):
				sj = arr[i].split(':')[1] + ":" + arr[i+1].split(':')[0]
				if chromo in refSJ:
					refSJ[chromo].append(sj)
				else:
					refSJ[chromo] = [sj]
	for gene in new:
		for iso in new[gene]:
			arr = new[gene][iso]
			chromo = arr[0]
			for i in range(1,(len(arr) -1)):
				sj = arr[i].split(':')[1] + ":" + arr[i+1].split(':')[0]
				if chromo in newSJ:
					newSJ[chromo].append(sj)
				else:
					newSJ[chromo] = [sj]
	knownSJ, novelSJ, detSJ, notDetSJ = 0,0,0,0
	for chromo in newSJ:
		for sj in newSJ[chromo]:
			if sj in refSJ[chromo]:
				knownSJ +=1 
			else:
				novelSJ += 1
				for gene in new:
					for iso in new[gene]:
						arr = new[gene][iso]
						for i in range(1,len(arr)-1):
							sjCheck = arr[i].split(':')[1] + ":" + arr[i+1].split(':')[0]
							if sj == sjCheck:
								#print(gene)
								if gene not in novelSJgenes:
									novelSJgenes.append(gene)
								
	#for g in novelSJgenes:
	#	print(g)
	#print("SJs present in ref:\t" + str(knownSJ))
	#print("SJs novel:\t" + str(novelSJ))

	for chromo in refSJ:
		for sj in refSJ[chromo]:
			if sj in newSJ[chromo]:
				detSJ += 1
			else:
				notDetSJ += 1
	#print("SJs present in ref and xome:\t" + str(detSJ))
	#print("SJs in ref not detected:\t" + str(notDetSJ))			









if __name__ == "__main__":
	main()
