#!/usr/bin/env python2.7
import io,re,sys
from collections import defaultdict

expDict = defaultdict(dict)
allXp = open("juncs_norm_asfilt_lowcountsFilt.csv", 'r')
for line in allXp:
	line = line.strip()
	if re.search('UF', line):
		continue
	arr = line.split(',', 1)
	expDict[arr[0]] = arr[1]

query = open(sys.argv[1], 'r')
for line in query:
	line = line.strip()
	print(line + ',' + expDict[line])
