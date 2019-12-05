#!/usr/bin/env python2.7

def main():
	import argparse, io
	parser = argparse.ArgumentParser(description = 'find premature stop codons using transdecoder .pep and .fa file made from gffread')
	parser.add_argument('-fa', action = 'store', type = str, required = True, help = '.fa file')
	parser.add_argument('-pep', action = 'store', type = str, required = True, help = '.pep file')
	args = parser.parse_args()
	nt = args.fa
	pep = args.pep
	ntFh = open(nt, 'r')
	pepFh = open(pep, 'r')
	process_files(nt = ntFh, pep = pepFh)

def process_files(nt, pep):
	juncDict, isoDict, longIsoDict, strandDict = makeNtDict(nt = nt)
	ptcDict, TSSdict = findPTC(juncs = juncDict, pep = pep)	
	sharedTSSisos = findSharedTSS(TSSdict = TSSdict, isoDict = isoDict, juncDict = juncDict, strandDict = strandDict)
	validatePTC(sharedTSS = sharedTSSisos, ptcDict = ptcDict, longIsos = longIsoDict, isoDict = isoDict, TSSdict = TSSdict, strandDict = strandDict, juncDict = juncDict)

def makeNtDict(nt):
	from collections import defaultdict
	import re
	juncDict = defaultdict(dict)
	isoDict = defaultdict(dict)
	longIsoDict = defaultdict(dict)
	strandDict = defaultdict(dict)
	##  >TCONS_00014356 gene=Potri.001G000400 loc:Chr01|8278-11349|- exons:8278-8666,9503-9619,10900-10972,11082-11177,11253-11349 segs:1-97,98-193,194-266,267-383,384-772
	##  AGCAAACCTTTAAAACTTACACCAGACAAGGCCTGCAAAATCACAAGACCAATCTCTTACCTAGAGACTC
	##  TAAATCAACTATTTCAAAATTAGCAATGCTCAAGAGCAATGAACTGGGTTCAACGTAAGATCTACCTCTA

	for line in nt:	
		line = line.strip()
		if re.search(r'^>', line):
			#print(line)
			if( len(line.split(' ')) > 1 ):
				tid = line.split(' ')[0].strip('>')
				gene = line.split(' ')[1].split('=')[1]
				strand = line.split(' ')[2].split('|')[2]
				strandDict[tid] = strand
				juncs = line.split(' ')[4].split(':')[1].split(',')
				exons = line.split(' ')[3].split(':')[1].split(',')
				juncDict[tid]["exons"] = exons
				juncDict[tid]["juncs"] = juncs
				longIsoDict = longCheck(longDict = longIsoDict, defline = line, juncs = juncs)
				
				isoDict[tid] = gene
				
	nt.close()
	return(juncDict, isoDict, longIsoDict, strandDict)		

def longCheck(longDict, defline, juncs):
	gene = defline.split(' ')[1].split('=')[1]
	tid = defline.split(' ')[0].strip('>')
	length =  int(juncs[(len(juncs)-1)].split('-')[1]) - int(juncs[0].split('-')[0])
	if(gene in longDict):
		#print(longDict[gene][tid])
		for iso in longDict[gene]:
			if(length > longDict[gene][iso]):
				del(longDict[gene][iso])
				longDict[gene][tid] = length
				break
				#longDict[gene] = iso 
	else:
		longDict[gene][tid] = length
		#longDict[gene] = tid
	return(longDict)

def findPTC(juncs, pep):
	import re
	from collections import defaultdict
	ptcDict = defaultdict(dict)
	TSSdict = defaultdict(dict)
	## >TCONS_00000011|m.28 TCONS_00000011|g.28  ORF TCONS_00000011|g.28 TCONS_00000011|m.28 type:5prime_partial len:248 (+) TCONS_00000011:1-744(+)

	for line in pep:
		line = line.strip()
		if re.search(r'^>', line):
			tid = line.split(' ')[0].split(':')[2].strip('>')
			transRange = line.split(' ')[7].split(':')[1].split('(')[0].split('-')  ###trying to get the locations for translation, sorry this is so nasty
			#print(str(transRange))
			#length = int(transRange[0]) - int(transRange[1])
			
			### edited to enable multiple TSS sites in sequence

			#if(tid in TSSdict):
			#	TSSdict[tid].append(transRange)
			#elif(tid not in TSSdict):
			#	TSSdict[tid] = [transRange]
			TSSdict[tid] = transRange


			stop = int(transRange[1])
			juncArr = juncs[tid]["juncs"]
			#print(str(juncArr))
			lastExon = juncArr[(len(juncArr) -1)]
			last3primeExon = int(lastExon.split('-')[0])
			last5primeExon = int(juncArr[len(juncArr) -2].split('-')[1])
			if(stop  < last5primeExon):
				if(stop <= (last5primeExon - 50)):
					ptcDict[tid] = stop
	return(ptcDict, TSSdict)

def findSharedTSS(TSSdict, isoDict, juncDict, strandDict):
	from collections import defaultdict
	sharedTSS = defaultdict(dict)
	for tid in TSSdict:
		gene = isoDict[tid]
		if(strandDict[tid] == "-"):
			tss = (int(TSSdict[tid][0]) -1)
			start = (int(juncDict[tid]["exons"][(len(juncDict[tid]["exons"]) -1)].split('-')[1]) - tss)

		elif(strandDict[tid] == "+"):
			tss = (int(TSSdict[tid][0]) -1)
			start = (int(juncDict[tid]["exons"][0].split('-')[0]) + tss)
			
		#start = TSSdict[tid][0]
		if start not in sharedTSS[gene]:
			sharedTSS[gene][start] = [tid]
		else:
			sharedTSS[gene][start].append(tid)
	return(sharedTSS)

def validatePTC(sharedTSS, ptcDict, longIsos, isoDict, strandDict, TSSdict, juncDict):
	for tid in ptcDict:
		gene = isoDict[tid]
		if(strandDict[tid] == "-"):
			tss = (int(TSSdict[tid][0]) -1)
			start = (int(juncDict[tid]["exons"][(len(juncDict[tid]["exons"]) -1)].split('-')[1]) - tss)
		elif(strandDict[tid] == "+"):
			tss = (int(TSSdict[tid][0]) -1)
			start = (int(juncDict[tid]["exons"][0].split('-')[0]) + tss)
		for longIso in longIsos[gene]:
			if(tid in sharedTSS[gene][start] and longIso in sharedTSS[gene][start] and len(sharedTSS[gene][start]) > 1):
				stop = getStopLoc(loc = ptcDict[tid], juncDict = juncDict, tid = tid)
				print(gene + '\t' + tid + "\t" + "genomic_pos:" + str(stop) + "\ttranscript_pos:" + str(ptcDict[tid]))

def getStopLoc(loc, juncDict, tid):
	for i in range(0, len(juncDict[tid]["juncs"])):
		left = int(juncDict[tid]["juncs"][i].split('-')[0]) 
		right = int(juncDict[tid]["juncs"][i].split('-')[1])
		if(int(loc) >= left and int(loc) <= right):
			diff = right - int(loc)
			exon = juncDict[tid]["exons"][i].split('-')[1]
			exonicPos = int(exon) - diff
	return(exonicPos)


		
if __name__ == "__main__":
        main()
