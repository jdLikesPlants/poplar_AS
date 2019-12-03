#!/usr/bin/env python2.7
import sys,re,io
from collections import defaultdict


def main():
	args = getArgs()
	tmap, herit = args.tmap, args.herit
	tmap, geneMap = storeTmap(tmap = tmap)
	herit = storeHerit(herit = herit)
	getHeritIsosPerGene(herit = herit, tmap = tmap)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = 'get number of heritable isoforms for regulated/regulator genes')
	parser.add_argument('-tmap', action = 'store', type = str, required = True, help = 'transmap file')
	parser.add_argument('-herit', action = 'store', type = str, required = True, help = 'two column file of h2 values for isoforms')
	args = parser.parse_args()
	return(args)

def getHeritIsosPerGene(herit, tmap):
	for gene in tmap:
		if len(tmap[gene]) > 1:	
			heritIso = 0
			for iso in tmap[gene]:
				if iso in herit:
					heritIso += 1
			print(gene + '\t' + str(heritIso))

def storeHerit(herit):
	import io,re,sys
	from collections import defaultdict
	fh = open(herit, 'r')
	herit = defaultdict(dict)
	for line in fh:
		line = line.strip()
		arr = line.split()
		herit[arr[0]] = arr[1]
	return(herit)

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
