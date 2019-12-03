#!/usr/bin/env python2.7
import sys,re,io
from collections import defaultdict


def main():
	args = getArgs()
	tmap, eqtl = args.tmap, args.eqtl
	tmap, geneMap = storeTmap(tmap = tmap)
	snpDict = storeEqtl(eqtl = eqtl)
	findMult(snpDict = snpDict, tmap = tmap, geneMap = geneMap)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = 'get number of heritable isoforms for regulated/regulator genes')
	parser.add_argument('-tmap', action = 'store', type = str, required = True, help = 'transmap file')
	parser.add_argument('-eqtl', action = 'store', type = str, required = True, help = 'matrix eqtl output file')
	args = parser.parse_args()
	return(args)


def findMult(tmap, geneMap, snpDict):
	from collections import defaultdict
	for snp in snpDict:	
		if len(snpDict[snp]) > 1:
			tidCounts = defaultdict(dict)
			tidTrack = defaultdict(dict)
			for tid in snpDict[snp]:
				gene = geneMap[tid]
				if gene in tidTrack:
					tidTrack[gene].append(tid)
				else:
					tidTrack[gene] = [tid]

				if gene in tidCounts:
					tidCounts[gene] += 1
				else:
					tidCounts[gene] = 1
			for gene in tidCounts:
				if tidCounts[gene] > 1:
					print(snp + '\t' + gene + '\t' + ','.join(tidTrack[gene]) + '\t' + str(tidCounts[gene]))
			

def storeEqtl(eqtl):
	import io, re, sys
	from collections import defaultdict
	snpDict = defaultdict(dict)
	eqtl = open(eqtl, 'r')
	for line in eqtl:
		line = line.strip()
		if re.search('snp', line):
			continue
		arr = line.split()	
		snp, tid = arr[0], arr[1]
		if snp in snpDict and tid not in snpDict[snp]:
			snpDict[snp].append(tid)
		else:
			snpDict[snp] = [tid]
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



if __name__ == "__main__":
	main()
