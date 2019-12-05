#!/usr/bin/env python2.7
import io,re,sys
from collections import defaultdict

geneCount = defaultdict(dict)

fh = open(sys.argv[1], 'r')
for line in fh:
	line = line.strip()
	arr = line.split()
	gene = arr[0]
	if gene in geneCount:
		geneCount[gene] += 1
	else:
		geneCount[gene] = 1
fh.close()
fh = open(sys.argv[1], 'r')
for line in fh:
	line = line.strip()
	arr = line.split()
	gene = arr[0]
	if geneCount[gene] > 1:
		print(line)
