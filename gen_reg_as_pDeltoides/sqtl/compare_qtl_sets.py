#!/usr/bin/env python2.7
import io,re,sys
from collections import defaultdict

sqtl = open(sys.argv[1], 'r')
eqtl = open(sys.argv[2], 'r')

sqtlList = []
eqtlList = []
shared = []

for line in sqtl:
	line = line.strip()
	if line not in sqtlList:
		sqtlList.append(line)
for line in eqtl:
	line = line.strip()
	if line not in eqtlList:
		eqtlList.append(line)

sqtlU = []
for sqtl in sqtlList:
	if sqtl in eqtlList and sqtl not in shared:
		shared.append(sqtl)
	elif sqtl not in eqtlList:
		sqtlU.append(sqtl)	

eqtlU = []
for eqtl in eqtlList:
	if eqtl in sqtlList and eqtl not in shared:
		shared.append(eqtl)
	elif eqtl not in sqtlList:
		eqtlU.append(eqtl)

print("num shared\t" + str(len(shared)))
print("eqtl unique\t" + str(len(eqtlU)))
print("sqtl unique\t" + str(len(sqtlU)))
