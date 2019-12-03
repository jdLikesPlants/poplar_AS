#!/usr/bin/env python2.7
import sys,re,io
from collections import defaultdict


def main():
	args = getArgs()
	tmap, cis, trans, geneLocs = args.tmap, args.cis, args.trans, args.geneloc
	geneLocs = storeGeneLocs(geneLocs = geneLocs)
	tmap, geneMap = storeTmap(tmap = tmap)
	snpDict = defaultdict(dict)
	snpDict = storeCisEqtl(eqtl = cis, snpDict = snpDict)
	snpDict = storeTransEqtl(eqtl = trans, snpDict = snpDict) 
	findMult(snpDict = snpDict, tmap = tmap, geneMap = geneMap, geneLocs = geneLocs)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = 'get number of heritable isoforms for regulated/regulator genes')
	parser.add_argument('-tmap', action = 'store', type = str, required = True, help = 'transmap file')
	parser.add_argument('-cis', action = 'store', type = str, required = True, help = 'matrix cis eqtl output file')
	parser.add_argument('-trans', action = 'store', type = str, required = True, help = 'matrix trans eqtl output file')
	parser.add_argument('-geneloc', action = 'store', type = str, required = True, help = 'genomic coordinates of genes')
	
	args = parser.parse_args()
	return(args)


def findMult(tmap, geneMap, snpDict, geneLocs):
	from collections import defaultdict
	snpTrack = defaultdict(dict)
	for gene in tmap:	
		if len(tmap[gene]) > 1:
			for tid in tmap[gene]:
				if tid in snpDict:
					snpTrack[gene][tid] = snpDict[tid][0]
					#snps = snpDict[tid]
					#for snp in snps:
					#	if "snp" in snpTrack and snp not in snpTrack[gene]["snp"]:
					#		snpTrack[gene]["snp"].append(snp)
					#	elif "snp" not in snpTrack[gene]:
					#		snpTrack[gene]["snp"] = [snp]
	#print(str(snpTrack))
	diffSnps = defaultdict(dict)
	for gene in snpTrack:
		diffCheck = 0
		for tid in snpTrack[gene]:
			for secTid in snpTrack[gene]:
				if tid != secTid:
					if snpTrack[gene][tid] != snpTrack[gene][secTid]:
						diffSnps[gene][tid] = snpTrack[gene][tid]
						diffSnps[gene][secTid] = snpTrack[gene][secTid]
						diffCheck = 1
						break
		if diffCheck == 1:
			sys.stdout.write(gene + '\t' + geneLocs[gene])
			for tid in diffSnps[gene]:
				sys.stdout.write('\t' + tid + '\t' + diffSnps[gene][tid])
			sys.stdout.write('\n')
					#if len(diffs) > 0:
					#	comps = tid + "_" + secTid
					#	revComp = secTid + "_" + tid
					#	for diff in diffs:
					#		if gene in diffSnps:
					#			if diff not in diffSnps[gene]:
					#				diffSnps[gene].append(diff)
					#		else:
					#			diffSnps[gene] = [diff]

						#if gene in diffSnps:
						#	if comps not in diffSnps[gene] and revComp not in diffSnps[gene]:
						#		diffSnps[gene][comps] = diffs
						#else:
						#	diffSnps[gene][comps] = diffs	
	#for gene in diffSnps:
		#loc = geneLocs[gene]
		#for comp in diffSnps[gene]:
		#print(gene + '\t' + loc + '\t' + ','.join(diffSnps[gene]))

			

def storeCisEqtl(eqtl, snpDict):
	import io, re, sys
	from collections import defaultdict
	eqtl = open(eqtl, 'r')
	for line in eqtl:
		line = line.strip()
		if re.search('snp', line):
			continue
		arr = line.split()
		fdr = float(arr[4])	
		snp, tid = arr[0], arr[1]
		snp = snp + "_cis"
		if tid in snpDict:
			if fdr < snpDict[tid][1]:
				snpDict[tid] = [snp, fdr]
		else:
			snpDict[tid] = [snp, fdr]
	return(snpDict)

def storeTransEqtl(eqtl, snpDict):
	import io, re, sys
	from collections import defaultdict
	eqtl = open(eqtl, 'r')
	for line in eqtl:
		line = line.strip()
		if re.search('snp', line):
			continue
		arr = line.split()	
		fdr = float(arr[4])
		snp, tid = arr[0], arr[1]
		snp = snp + "_trans"
		if tid in snpDict:
			if fdr < snpDict[tid][1]:
				snpDict[tid] = [snp, fdr]
		else:
			snpDict[tid] = [snp, fdr]
	return(snpDict)

		
def storeTmap(tmap):
	import io,re,sys
	from collections import defaultdict
	tmap = open(tmap, 'r')
	transMap = defaultdict(dict)
	geneMap = defaultdict(dict)
	for line in tmap:
		line = line.strip()
		arr = line.split()
		gene, tid = arr[0], arr[1]
		if gene in transMap:
			transMap[gene].append(tid)
		else:
			transMap[gene] = [tid]
		geneMap[tid] = gene
	return(transMap, geneMap)

def storeGeneLocs(geneLocs):
	locs = defaultdict(dict)
	geneLocs = open(geneLocs, 'r')
	for line in geneLocs:
		line = line.strip()
		if re.search('gene', line):
			continue
		arr = line.split()
		gene = arr[0]
		chromo = arr[1]
		loc = arr[2] + '-' + arr[3]
		loc = chromo + '_' + loc
		locs[gene] = loc
	return(locs)

if __name__ == "__main__":
	main()
