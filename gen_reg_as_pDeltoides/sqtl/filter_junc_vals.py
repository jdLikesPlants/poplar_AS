#!/usr/bin/env python2.7
import io,re,sys

fh = open(sys.argv[1], 'r')
for line in fh:
	line = line.strip()
	if re.search('UF',line):
		print(line)
		continue
	arr = line.split(',')
	junc = arr[0]
	count = 0
	thresh = (len(arr)-1)*.05
	for i in range(1, len(arr)):
		if float(arr[i]) > 0:
			count += 1
	if count >= thresh:
		print(line)
