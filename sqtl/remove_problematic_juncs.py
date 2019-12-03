#!/usr/bin/env python2.7
import io,re,sys

fh = open(sys.argv[1], 'r')
for line in fh:
	line = line.strip()
	if re.search("UF",line):
		print(line)
		continue
	check = 0
	arr = line.split(',')
	for i in range(1,len(arr)):
		if float(arr[i]) > 1:
			check = 1
			break
	if check != 1:
		print(line)
