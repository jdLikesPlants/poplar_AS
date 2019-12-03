#!/usr/bin/env python2.7
import io,re,sys
from collections import defaultdict

snpDict = defaultdict(dict)
snps = open(sys.argv[1], 'r')
for line in snps:
	line = line.strip()
	if re.search("snp", line):
		continue
	arr = line.split('\t')
	snp = arr[0]
	junc = arr[1]
	if junc in snpDict:
		snpDict[junc].append(snp)
	else:
		snpDict[junc] = [snp]

regFile = open(sys.argv[2])
for line in regFile:
	line = line.strip()
	arr = line.split('\t')
	junc = arr[0] + ':' + arr[4]
	#print(junc)
	if junc in snpDict:
		print(line + '\t' + ','.join(snpDict[junc]))
