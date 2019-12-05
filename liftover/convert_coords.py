#!/usr/bin/env python2.7
import io
import sys
import re
import argparse
from collections import defaultdict

def main():
	parser = argparse.ArgumentParser(description = 'updates old bed file with new coordinates')	
	parser.add_argument('-inbed', action = 'store', type = str, required = True, help = 'old bed file')
	parser.add_argument('-conversions', action = 'store', type = str, required = True, help = 'crossmap output file')
	args = parser.parse_args()
	bedFile = args.inbed
	conversions = args.conversions
	process_files(bedFile, conversions)	

def process_files(bed, conv):
	inBed = open(bed, "r")
	inConv = open(conv, "r")
	
	convDict = getConvInfo(inConv)
	updateCoords(inBed, convDict)

def getConvInfo(conv):
	convDict = {}
	for line in conv:
		line = line.strip()
		if re.search('Fail', line):
			continue
		else:
			linearr = line.split('->')
			for i in range(0, len(linearr)):
				linearr[i] = linearr[i].strip()
			oldPos = linearr[0].split('\t')[2]
			oldChrom = linearr[0].split('\t')[0].strip()
			oldPos = oldChrom + "\t" + oldPos
			newPos = linearr[1].split('\t')[2]
			chromo = linearr[1].split('\t')[0].strip()
			newPos = chromo + "\t" + newPos
			#print(oldPos + '\t' + newPos)
			convDict[oldPos] = newPos
	return convDict

def updateCoords(bed, convs):
	#outfile = open('updated_coords.bed', "w")
	for line in bed:
		line = line.strip()
		linearr = line.split('.', 1)
		chromo = linearr[0].split('\t')[0]
		old = linearr[0].split('\t')[2]
		old = chromo + "\t" + old
		
		#print(str(old))
		#print(str(convs[old]))
		#print(linearr[0] + linearr[1])
		if old in convs:
			print(convs[old] + "\t" + "." + linearr[1])
			#outfile.write(convs[old] + "\t" + "." + linearr[1] + "\n")

if __name__ == "__main__":
	main()
