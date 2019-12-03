#!/usr/bin/env python2.7
import io,re,sys
from collections import defaultdict

mapFile = open(sys.argv[1], 'r')
sqtlOut = open(sys.argv[2], 'r')

mapDict = defaultdict(dict)
for line in mapFile:
	line = line.strip()
	arr = line.split()
	feat = arr[3] + ":" + arr[4]
	mapDict[arr[0]][arr[1]] = feat

#print(str(mapDict["Chr18"]))

for line in sqtlOut:
	line = line.strip()
	arr = line.split()
	t = arr[0].split('_')
	if re.search("scaffold", line):
		chromo = t[0] + '_' + t[1]
		loc = t[2]
	else:
		chromo, loc = t[0], t[1]
	#print(line)
	if chromo in mapDict:
		feat = mapDict[chromo][loc]
		feat = feat.split(":")
		print(line + '\t' + '\t'.join(feat))
