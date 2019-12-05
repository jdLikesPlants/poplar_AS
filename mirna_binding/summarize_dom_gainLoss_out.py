#!/usr/bin/env python2.7
import io,re,sys
from collections import defaultdict
import numpy

### get number of gains/losses per gene
occDict = defaultdict(dict)

fh = open(sys.argv[1], 'r')
for line in fh:
	line = line.strip()
	if re.search("longest", line):
		continue
	arr = line.split('\t')
	gene, tid = arr[0], arr[1]
	occ = arr[2]
	if gene in occDict:
		if occ in occDict[gene]:
			occDict[gene][occ] += 1
		else:
			occDict[gene][occ] = 1
	else:
		occDict[gene][occ] = 1
#print(str(occDict))
### summary stats of number of gains/losses
### how many genes have both

numLosses = []
numGains = []
numBoth = 0

for gene in occDict:
	if "gain" in occDict[gene] and "loss" in occDict[gene]:
		numBoth += 1
	elif "gain" in occDict[gene]:
		numGains.append(occDict[gene]["gain"])
	elif "loss" in occDict[gene]:
		numLosses.append(occDict[gene]["loss"])

print("number of genes with both gain and loss isoforms" + '\t' + str(numBoth))

print("number of genes with loss isoforms\t" + str(len(numLosses)))
print("avg/median num loss isoforms\t" + str(numpy.mean(numLosses)) + '\t' + str(numpy.median(numLosses)))
print("number of genes with gain isoforms\t" + str(len(numGains)))
print("avg/median num loss isoforms\t" + str(numpy.mean(numGains)) + '\t' + str(numpy.median(numGains)))
