#!/usr/bin/env python2.7
import io,re,sys
from collections import defaultdict

tmap = open(sys.argv[1], 'r')
inFile = open(sys.argv[2], 'r')
tDict = defaultdict(dict)
for line in tmap:
	arr = line.strip().split()
	gene, tid = arr[0], arr[1]
	tDict[tid] = gene

for line in inFile:
	line = line.strip()
	print(tDict[line])
