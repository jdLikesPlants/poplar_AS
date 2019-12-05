#!/usr/bin/env python2.7
import io,sys,re
from collections import defaultdict

medTids = sys.argv[1]
domGainLoss = sys.argv[2]
tids = []
fh = open(medTids, 'r')
for line in fh:
	line = line.strip()
	tids.append(line)
fh.close()

domDict = defaultdict(dict)

fh = open(domGainLoss, 'r')
for line in fh:
	line = line.strip()
	if re.search('longest', line):
		tid = line.split('\t')[2]
		gene = line.split('\t')[0]
		if tid in tids:
			domDict[gene]["longest"] = line
	else:
		tid = line.split('\t')[1]
		gene = line.split('\t')[0]
		if tid in tids:
			domDict[gene][tid] = line

for gene in domDict:
	if "longest" in domDict[gene]:
		nonLong = 0
		for i in domDict[gene]:
			if i != "longest":
				nonLong += 1
		if nonLong != 0:
			#print(str(domDict[gene]))
			print(domDict[gene]["longest"])
			for i in domDict[gene]:
				if i != "longest":
					print(domDict[gene][i])


