#!/usr/bin/env python2.7

def main():
	import argparse, io
	parser = argparse.ArgumentParser(description = 'find genomic locations of hmmscan domains')
	parser.add_argument('-fa', action = 'store', type = str, required = True, help = '.fa file')
	parser.add_argument('-mirna', action = 'store', type = str, required = True, help = 'miRNA binding sites from psRNAtarget')
	args = parser.parse_args()
	nt = args.fa
	nt = open(nt, "r")
	rna = open(args.mirna, "r")
	processFiles(nt = nt, rna = rna)

def processFiles(nt, rna):
	from find_PTC import makeNtDict, longCheck, findPTC
	juncDict, isoDict, longIsoDict, strandDict = makeNtDict(nt = nt)
	#ptcDict, TSSdict = findPTC(juncs = juncDict, pep = pep)
	convertCoords(juncDict = juncDict, rna = rna, strandDict = strandDict)

def convertCoords(juncDict, strandDict, rna):
	import re
	print("isoform\tgene\tmirna\teft\tright\ttype\tmultiplicity")
	for line in rna:
		line = line.strip()
		if re.search('#', line):
			continue
		arr = line.split("\t")
		### change array indices here depending on input file format
		tid = arr[1]
		mirna = arr[0]
		#clan = arr[1]
		#print(line)
		desc = arr[10]
		gene = arr[11].split('=')[1]
		ntLeft = int(arr[6])
		ntRight = int(arr[7])
		genomicLeft, genomicRight = getGenomicCoords(juncDict = juncDict, left = ntLeft, right = ntRight, \
			tid = tid, strandDict = strandDict)
		if genomicLeft == 0 or genomicRight == 0:
			continue
		### change to output pfam id instead of clan
		#if pf != "NULL":
		#print(pf + "\t" + tid + "\t" + str(genomicLeft) + "\t" + str(genomicRight) + "\t" + "0" + "\t" + arr[1])
		print(tid + "\t" + gene + '\t' + mirna + '\t' + str(genomicLeft) + "\t" + str(genomicRight) + "\t" + desc + '\t' + arr[-1])
		#else:
		#	print(clan + "\t" + tid + "\t" + str(genomicLeft) + "\t" + str(genomicRight) + "\t" + "0" + "\t" + arr[1])

def getGenomicCoords(juncDict, left, right, tid, strandDict):
	import sys
	#print(tid)
	tidName = tid
	tidArr = juncDict[tidName]["juncs"]
	exons = juncDict[tid]["exons"]
	nleft, nright = left, right
	if strandDict[tid] == "-":
			tss = int(exons[-1].split('-')[1])
	else:
			tss = int(exons[0].split('-')[0])
	newLeft, newRight = 0,0
	if 0 < 1:	
		for i in range (0, len(tidArr)):
			juncLeft = int(tidArr[i].split('-')[0])
			juncRight = int(tidArr[i].split('-')[1])  
			#print("junc left\t" + str(juncLeft) + '\tjunc right\t' + str(juncRight))
			### get the tidArr index of ntLeft and ntRight and then get values from exon arr	
			### diff of left and right exonarr values with nt values to get genomic pos
			#print(str(strandDict[tidName]))
			if(nleft >= juncLeft and nleft <= juncRight and newLeft == 0):
				#print(str(nleft))
				diff = nleft - juncLeft
				#print(str(diff))
				#exon = int(juncDict[tid]["exons"][i].split('-')[0])
				if(strandDict[tidName] == "-"):
					
					exon = int(juncDict[tidName]["exons"][((len(tidArr) -1) - i)].split('-')[1])
					newLeft = exon - diff
					#print(str(newLeft))
				else:
					exon = int(juncDict[tidName]["exons"][i].split('-')[0])
					#print(str(exon))
					newLeft = exon + diff
					#print(str(newLeft))			
		
			if(nright >= juncLeft and nright <= juncRight and newRight == 0):
				diff = juncRight - nright
				#print(str(diff))
				if(strandDict[tidName] == "-"):
					exon = int(juncDict[tidName]["exons"][((len(tidArr) -1) - i)].split('-')[0])
					newRight = exon + diff
				else:
					exon = int(juncDict[tidName]["exons"][i].split('-')[1])
					newRight = exon - diff
				#print(str(newRight))		
			#if newLeft != 0 and newRight != 0:
			#	break
			if newLeft != 0 and newRight != 0:
				break
		#if(tidName == "TCONS_00000112"):
		#	print(str(newLeft) + '\t' + str(newRight) + '\t' + str(juncLeft) + '\t' + str(juncRight))
		if(strandDict[tidName] == "-"):
			nL = newRight
			nR = newLeft
			newLeft = nL
			newRight = nR
		
	#print(str(newLeft) + '\t' + str(newRight))
	return(newLeft, newRight)
	

if __name__ == "__main__":	
	main()
