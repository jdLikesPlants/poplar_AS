#!/usr/bin/env python2.7

def main():
	args = getArgs()
	tmap, dat, dom, gtf, longs = args.tmap, args.dat, args.mirna, args.gtf, args.longs
	longs = getLongTids(longs)
	transMap, geneDict  = storeTmap(tmap = tmap)
	exonDict, chromoDict = parseGtf(gtf = gtf)  
	asDict, asCoordDict = parseDat(dat = dat, geneDict = geneDict, transMap = transMap, exonDict = exonDict)
	domainsDict, geneDomains = parseDomains(dom = dom, geneDict = geneDict)
	findChange(transMap = transMap, asDict = asDict, domainsDict = domainsDict, asCoordDict = asCoordDict, geneDomains = geneDomains, exonDict = exonDict, longTids = longs)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = 'identify gain/loss of protein domains as a result of alternative splicing')
	parser.add_argument('-tmap', action = 'store', type = str, required = True, help = 'transmap file')
	parser.add_argument('-dat', action = 'store', type = str, required = True, help = 'pasa .dat file')
	parser.add_argument('-mirna', action = 'store', type = str, required = True, help = 'mirna genomic domains file')
	parser.add_argument('-gtf', action = 'store', type = str, required = True, help = 'filtered gtf file')
	parser.add_argument('-longs', action = 'store', type = str, required = True, help = '2 col file of prot coding longest tids for each gene')
	args = parser.parse_args()
	return(args)

def getLongTids(longs):
	import sys,io,re
	from collections import defaultdict
	longTids = defaultdict(dict)
	longs = open(longs, 'r')
	for line in longs:
		arr = line.strip().split()
		gene, tid = arr[0], arr[1]
		longTids[gene] = tid
	return(longTids)

def storeTmap(tmap):
	import re,io
	from collections import defaultdict
	transToGene = defaultdict(dict)
	geneToTrans = defaultdict(dict)
	geneList = []
	fh = open(tmap, "r")
	for line in fh:
		line = line.strip()
		arr = line.split('\t')
		gene, tid = arr[0], arr[1]
		if gene in transToGene:
			transToGene[gene].append(tid)
		else:
			transToGene[gene] = [tid]
		if gene not in geneList:
			geneList.append(gene)
		geneToTrans[tid] = gene
	return(transToGene, geneToTrans)

def parseGtf(gtf):
	import re,io
	from collections import defaultdict
	exonDict = defaultdict(dict)
	chromoDict = defaultdict(dict)
	fh = open(gtf, "r")
	for line in fh:
		line = line.strip()
		if(re.search('exon', line)):
			arr = line.split('\t')
			chromo = arr[0]
			tid = arr[8].split(';')[0].split('\"')[1]
			if not tid in chromoDict:
				chromoDict[tid] = chromo
			left = arr[3]
			right = arr[4]
			exons = left + ":" + right
			if tid not in exonDict:
				exonDict[tid] = [exons]
			elif tid in exonDict:
				exonDict[tid].append(exons)
	return(exonDict, chromoDict)

def parseDat(dat, geneDict, transMap, exonDict):
	import io,re
	from collections import defaultdict
	asTidDict = defaultdict(dict)
	asCoordDict = defaultdict(dict)
	typeDict = {"retained_intron" : "RI",
				"spliced_intron" : "SI" ,
				"retained_exon" : "RE" ,
				"skipped_exon" : "SE" , 
				"alternate_exon" : "AE"}
	fh = open(dat, "r")
	for line in fh:
		line = line.strip()
		arr = line.split('\t')
		
		if(re.search('alt_acceptor', line) and re.search('TCONS', line)):
			arr = line.split()
			asmblArr = line.split()[7].split(',')
			tid = ""
			strand = arr[6]
			### get more accurate locus here with getA3Coords
			loc = arr[0] + ":" + arr[4] + ':' + arr[5] + ':' + arr[6]
			if(arr[6] == "+"):
				left, right = int(arr[4]), int(arr[5])
			elif(arr[6] == "-"):
				left, right = int(arr[5]), int(arr[4])

			for i in range(0, len(asmblArr)):
				if(re.search("TCONS", asmblArr[i])):
					tid = asmblArr[i]
					#print(tid)
					gene = geneDict[tid]
					if(len(transMap[gene]) < 2):
						continue
					#left ,right = getA3Coords(left = left, right = right, exonDict = exonDict, tid = tid, loci = loc)
					loc = arr[0] + ":" + str(left) + ":" + str(right) + ":" + arr[6]
					if gene in asCoordDict:
						if "A3" in asCoordDict[gene]:
							if loc not in asCoordDict[gene]["A3"]:
								asCoordDict[gene]["A3"].append(loc)	
						else:
							asCoordDict[gene]["A3"] = [loc]
					if gene not in asCoordDict:
						asCoordDict[gene]["A3"] = [loc]
					
					if(tid in asTidDict):
						if("A3" in asTidDict[tid]):
							asTidDict[tid]["A3"].append(loc) 
						else:
							asTidDict[tid]["A3"] = [loc]
					elif(tid not in asTidDict):
						asTidDict[tid]["A3"] = [loc]
			#asDict[ssType][loc] = incl
			#asAllDict[ssType][loc] = allTid
		
		elif(re.search('alt_donor', line) and re.search('TCONS', line)):
			arr = line.split()
			asmblArr = arr[7].split(',')
			tid = ""
			loc = arr[0] + ":" + arr[4] + ':' + arr[5] + ':' + arr[6]
			if(arr[6] == "+"):
				left , right = int(arr[4]), int(arr[5])
			elif(arr[6] == "-"):
				left, right = int(arr[5]), int(arr[4])

			for i in range(0, len(asmblArr)):
				if(re.search("TCONS", asmblArr[i])):
					tid = asmblArr[i]
					gene = geneDict[tid]
					if(len(transMap[gene]) < 2):
						continue
					#left ,right = getA5Coords(left = left, right = right, exonDict = exonDict, tid = tid, loci = loc)
					loc = arr[0] + ":" + str(left) + ":" + str(right) + ":" + arr[6]
					
					gene = geneDict[tid]
					if(gene in asCoordDict):
						if "A5" in asCoordDict[gene]:
							if loc not in asCoordDict[gene]["A5"]:
								asCoordDict[gene]["A5"].append(loc)
						else:
							asCoordDict[gene]["A5"] = [loc]
					elif(gene not in asCoordDict):
						asCoordDict[gene]["A5"] = [loc]

					if(tid in asTidDict):
						if("A5" in asTidDict[tid]):
							asTidDict[tid]["A5"].append(loc)
						else:
							asTidDict[tid]["A5"] = [loc]
					elif(tid not in asTidDict):
						asTidDict[tid]["A5"] = [loc]
		elif(re.search('alternate_exon', line) and re.search('TCONS', line)):
			check = 0
			arr = line.split()
			asmblArr = arr[7].split(',')
			tid = arr[7]
			loc = arr[0] + ":" + arr[4] + ':' + arr[5] + ':' + arr[6]
			if(arr[6] == "+"):
				left , right = int(arr[4]), int(arr[5])
			elif(arr[6] == "-"):
				left, right = int(arr[5]), int(arr[4])
			for exonCoords in exonDict[tid]:
				coords = exonCoords.split(':')
				coords = [int(x) for x in coords]
				#print(tid + '\t' + str(left) + '\t' + str(right) + '\t' + str(coords))
				if left in coords and right in coords:
					check = 1
			if check == 0:
				continue
				
			for i in range(0, len(asmblArr)):
				if(re.search("TCONS", asmblArr[i])):
					tid = asmblArr[i]
					gene = geneDict[tid]
					if(len(transMap[gene]) < 2):
						continue
					#left ,right = getA5Coords(left = left, right = right, exonDict = exonDict, tid = tid, loci = loc)
					loc = arr[0] + ":" + str(left) + ":" + str(right) + ":" + arr[6]
					gene = geneDict[tid]
					if(gene in asCoordDict):
						if "AE" in asCoordDict[gene]:
							if loc not in asCoordDict[gene]["AE"]:
								asCoordDict[gene]["AE"].append(loc)
						else:
							asCoordDict[gene]["AE"] = [loc]
					elif(gene not in asCoordDict):
						asCoordDict[gene]["AE"] = [loc]

					if(tid in asTidDict):
						if("AE" in asTidDict[tid]):
							asTidDict[tid]["AE"].append(loc)
						else:
							asTidDict[tid]["AE"] = [loc]
					elif(tid not in asTidDict):
						asTidDict[tid]["AE"] = [loc]
		




		elif(re.search('retained', line) or re.search('spliced', line) or re.search('skipped',line) and re.search('TCONS', line)):
 
	
			arr = line.split()
			asmblArr = arr[7].split(',')
			tid = ""
			asType = typeDict[arr[3]]
			for i in range(0, len(asmblArr)):
				if(re.search("TCONS", asmblArr[i])):
					tid = asmblArr[i]
					gene = geneDict[tid]
					loc = arr[0] + ":" + arr[4] + ":" + arr[5] + ':' + arr[6]
					if(gene in asCoordDict):
						if(asType in asCoordDict[gene]):
							if(loc not in asCoordDict[gene][asType]):
								asCoordDict[gene][asType].append(loc)
						else:
							asCoordDict[gene][asType] = [loc]
					if(gene not in asCoordDict):
						asCoordDict[gene][asType] = [loc]

					if(tid in asTidDict):
						if(asType in asTidDict[tid]):
							asTidDict[tid][asType].append(loc)
						else:
							asTidDict[tid][asType] = [loc]
					elif(tid not in asTidDict):
						asTidDict[tid][asType] = [loc]

	return(asTidDict, asCoordDict)


def parseDomains(dom, geneDict):
	import io,re
	from collections import defaultdict
	domainsDict, geneDomains = defaultdict(dict), defaultdict(dict)
	fh = open(dom, "r")
	for line in fh:
		line = line.strip()
		arr = line.split('\t')
		if re.search("isoform", line):
			continue
		#print(line)
		### change here to use clan instead of domain
		domain = arr[2]
		#tid = arr[1].split(':')[2]
		tid = arr[0]
		if tid in geneDict:
			gene = geneDict[tid]
		else:
			continue
		loc = arr[3] + ":" + arr[4]
		### domainsDict[tid][domain] = [loci]		
		#print(str(gene))

		if tid in domainsDict:
			if domain in domainsDict[tid]:
				if loc not in domainsDict[tid][domain]:
					domainsDict[tid][domain].append(loc)
			elif domain not in domainsDict[tid]:
				domainsDict[tid][domain] = [loc]
		else:
			domainsDict[tid][domain] = [loc]

		if gene in geneDomains:
			if domain not in geneDomains[gene]:
				geneDomains[gene].append(domain)
		elif gene not in geneDomains:
			geneDomains[gene] = [domain]
	return(domainsDict, geneDomains)

def findChange(transMap, longTids, asDict, domainsDict, asCoordDict, geneDomains, exonDict):
	from collections import defaultdict
	import re, sys
	from check_AS import checkAS
	#print(str(transMap["Potri.008G065800"]))
	#print(str(asDict["TCONS_00231440"]))
	#print(str(domainsDict["tRNA_int_endo_N"]))  ##gain tids in the dict
	#print(str(domainsDict))  ##gain tids in the dict
	#print(str(asCoordDict["Potri.008G044800"]))
	#print(str(geneDomains["Potri.011G058400"]))
	#print(str(exonDict["TCONS_00116758"]))	

	for gene in transMap:
		if len(transMap[gene]) > 1 and gene in longTids:
		### which tids gain/lose relative to the intercept...
			longTid = longTids[gene]
			
			if longTid in domainsDict:
				domains = defaultdict(dict)
				domains = domainsDict[longTid]
				#for domain in domainsDict[longTid]:
				#	domains.append(domain)
			else:
				domains = "none"
			if domains != "none":
				tidDict = defaultdict(dict)
				for tid in transMap[gene]:
					if tid != longTid:
						tidGains, tidLosses, maintains = [],[], []
						if tid in domainsDict and tid != longTid:
							tidDoms = domainsDict[tid]
							for iDoms in domains:
								if iDoms not in tidDoms:
									tidLosses.append(iDoms)
								else:
									if iDoms not in maintains:
										maintains.append(iDoms)
							for doms in tidDoms:
								if doms not in domains:
									tidGains.append(doms)
						elif tid not in domainsDict and tid != longTid:
							#print(gene + '\t' + tid + str(domains))
							for iDoms in domains:
								tidLosses.append(iDoms)
						if len(tidGains) > 0:
							tidDict[tid]["gains"] = tidGains
						if len(tidLosses) > 0:
							tidDict[tid]["losses"] = tidLosses
						#tidDict[tid]["maintains"] = maintains
				
			


						#if len(tidGains) > 0 or len(tidLosses) > 0 or len(maintains) > 0:
							#tidDict[tid] = str('\t'.join([gene, tid, "gains", ','.join(tidGains), "losses", ','.join(tidLosses), "maintains", ','.join(maintains) ]))
						#if len(tidGains) > 0 or len(tidLosses) > 0:
							
							#tidDict[tid] = str('\t'.join([gene, tid, "gains", ','.join(tidGains), "losses", ','.join(tidLosses) ]))
				if len(tidDict) > 0:
					#print(gene + '\t' + "longest iso\t" + longTid + '\t' + ','.join(domains))
					sys.stdout.write(gene + '\t' + "longest iso\t" + longTid)
					dL = []
					for d in domains:
						dL.append(str(d + ':' + ','.join(domains[d])))
						#sys.stdout.write('\t' + d + ':' + str(domains[d]))
					sys.stdout.write('\t' + '; '.join(dL))
					sys.stdout.write('\n')
 
					checkAS(longTid, asDict, asCoordDict, tidDict, domainsDict, geneDomains, exonDict, gene)		

			#	if len(tidDict) > 0:
			#		print(gene + '\t' + "longest iso\t" + longTid + '\t' + ','.join(domains))
			#		for tid in tidDict:
			#			print(tidDict[tid])
									
			if domains == "none":
				tidDict = defaultdict(dict)
				for tid in transMap[gene]:
					if tid in domainsDict and tid != longTid:
						doms = domainsDict[tid]
						#tidDict[tid] = str(gene + '\t' + tid + '\t' + "gains\t" + ','.join(doms) + "\tlosses\t\tmaintains\t")
						
						tidDict[tid]["gains"] = doms
						tidDict[tid]["losses"] = []

					#elif tid not in domainsDict and tid != longTid:
					#	tidDict[tid] = str(gene + '\t' + tid + '\t' + "gains\t\t" + "\tlosses\t\tmaintains\t")
						
				if len(tidDict) > 0:
					print(gene + '\t' + "longest iso\t" + longTid + '\tnone')
					checkAS(longTid, asDict, asCoordDict, tidDict, domainsDict, geneDomains, exonDict, gene)		

												
def domCheck(asType, lossType, hasLeft, hasRight, lossLeft, lossRight, domLeft, domRight, loci, asDict, lossTid):
	hasCheck, lossCheck = 0, 0
	#print(str(domLeft) + '-'+ str(domRight) + '\t'  + lossType +  '\t' + str(hasLeft) +'\t' + str(hasRight) + '\t' + str(lossLeft)  + '-' + str(lossRight))
	
####### AE ##########
	if(asType == "AE" and hasCheck == 0 and lossType == "AE"):
		if(hasLeft >= domRight and hasRight <= domLeft):								
			hasCheck = 1
		if(hasLeft <= domLeft and hasRight >= domRight):					## 5' or 3' AE, whole domain in exon
			hasCheck = 1		
		elif(hasLeft >= domLeft and hasRight >= domRight and hasLeft <= domRight):	## 3' AE
			hasCheck = 1
		elif(hasLeft <= domLeft and hasRight <= domRight and hasRight >= domLeft):		## 5' AE
			hasCheck = 1
		if(hasCheck == 1):
			lossCheck = 1
			
			
			
####### AE ##########

#### A3 A5 ##########
	elif(asType == "A5" and lossType == "A5"):
		#print("has left\t" + str(hasLeft)+ "\thas Right\t" + str(hasRight) + "\tdom Left\t" + str(domLeft) + "\tdom right\t" + str(domRight))
		if(hasLeft <= domRight and hasLeft >= domLeft):								
			hasCheck = 1
		elif(hasRight >= domLeft and hasRight <= domRight):
			hasCheck = 1
	
		if(lossLeft == hasLeft or lossRight == hasRight):
			if(hasCheck == 0):
				if(lossLeft >= domLeft and lossLeft <= domRight):
					lossCheck, hasCheck = 1, 1
				elif(lossRight <= domRight and lossRight >= domLeft):
					lossCheck, hasCheck = 1, 1
			else:
				lossCheck = 1
				
	elif(asType == "A3" and lossType == "A3"):
		if(hasLeft <= domRight and hasLeft >= domLeft):								
			hasCheck = 1
		elif(hasRight >= domLeft and hasRight <= domRight):
			hasCheck = 1
	
		if(lossLeft == hasLeft or lossRight == hasRight):
			if(hasCheck == 0):
				if(lossLeft >= domLeft and lossLeft <= domRight):
					lossCheck, hasCheck = 1, 1
				elif(lossRight <= domRight and lossRight >= domLeft):
					lossCheck, hasCheck = 1, 1
			else:
				lossCheck = 1		
				
	
#### A3 A5 ##########			
			
#### SI RI ##########

### needs to confirm if in domain
	elif(asType == "SI" and lossType == "RI"):
		hasCheck, lossCheck = checkSI(asDict = asDict, leftCoord = hasLeft, rightCoord = hasRight, loci = loci, riLeft = lossLeft, \
									  riRight = lossRight, domLeft = domLeft, domRight = domRight)
	elif(asType == "RI" and lossType == "SI"):
		hasCheck, lossCheck = checkRI(asDict = asDict, leftCoord = hasLeft, rightCoord = hasRight, loci = loci, siLeft = lossLeft, \
									  siRight = lossRight, domLeft = domLeft, domRight = domRight)			
#### SI RI ##########

#### SE RE ##########
	elif(asType == "RE" and lossType == "SE"):
		hasCheck, lossCheck = checkRE(asDict = asDict, leftCoord = hasLeft, rightCoord = hasRight, domLeft = domLeft, \
									  domRight = domRight, seLeft = lossLeft, seRight = lossRight)
	elif(asType == "SE" and lossType == "RE"):
		hasCheck, lossCheck = checkSE(asDict = asDict, leftCoord = hasLeft, rightCoord = hasRight, domLeft = domLeft, \
									  domRight = domRight, reLeft = lossLeft, reRight = lossRight)
#### SE RE ##########
	return(hasCheck, lossCheck)
################################################			
			
					
										
			

def getA3Coords(left, right, exonDict, tid, loci):
	import re
	exons = exonDict[tid]
	#print(str(exons))
	if(re.search('-', loci)):
		right -= 1
		left -= 1
		for i in range(0, len(exonDict[tid])):
			exonArr = exons[i].split(':')
			#print(tid + '\t' + str(left) + '\t' + exonArr[1])
			if(right == int(exonArr[1])):
				exRight = exons[i + 1].split(':')[0]
				right = int(exRight)
				#print(tid + '\t' + str(left) + '\t' + str(right))
						
	else:
		right += 1 
		for i in range(0, len(exonDict[tid])):
						
			exonArr = exons[i].split(':')
			if(right == int(exonArr[0])):
				exLeft = exons[i -1].split(':')[1]
				left = int(exLeft)

	return(left, right)

def getA5Coords(left, right, exonDict, tid, loci):
	import re
	exons = exonDict[tid]
	if(re.search('-', loci)):
		left += 1 
		for i in range(0, len(exonDict[tid])):
			exonArr = exons[i].split(':')
			if(left == int(exonArr[0])):
				exLeft = exons[i -1].split(':')[1]
				left = int(exLeft)
	else:
		left -= 1
		for i in range(0, len(exonDict[tid])):
			exonArr = exons[i].split(':')
			#print(tid + '\t' + str(left))
			#print(str(exonArr))
			if(left == int(exonArr[1])):
				exRight = exons[i + 1].split(':')[0]
				right = int(exRight)
	return(left, right)

def checkSI(asDict, leftCoord, rightCoord, riLeft, riRight, domLeft, domRight, loci):	
	import re
	domCheck, lossCheck = 0, 0
	#if("RI" in asTypes):
	#for lossTid in lacksDom:
	if(re.search('-', loci)):

		if(leftCoord <= domRight and leftCoord >= domLeft):                                                             
				domCheck = 1
		elif(rightCoord >= domLeft and rightCoord <= domRight):
				domCheck = 1
	else:
		if(leftCoord >= domLeft and leftCoord <= domRight):                                                             
				domCheck = 1
		elif(rightCoord <= domRight and rightCoord >= domLeft):
				domCheck = 1

	if(riLeft == leftCoord and riRight == rightCoord):
		lossCheck = 1
		#if gene == "Potri.001G008400":
		#	print(gene + '\t' + loc + '\t' + str(asTypes) + '\t' + str(lacksAS))
		
	return(domCheck, lossCheck)

def checkRI(asDict, leftCoord, rightCoord, siLeft, siRight, domLeft, domRight, loci):
	import re
	domCheck, lossCheck = 0, 0
	#for lossTid in lacksDom:			
	if(re.search('-', loci)):
		if(leftCoord <= domRight and leftCoord >= domLeft):                                                             
				domCheck = 1
		elif(rightCoord >= domLeft and rightCoord <= domRight):
				domCheck = 1
	else:
		if(leftCoord >= domLeft and leftCoord <= domRight):                                                             
				domCheck = 1
		elif(rightCoord <= domRight and rightCoord >= domLeft):
				domCheck = 1

	if(siLeft == leftCoord and siRight == rightCoord):
		lossCheck = 1
		#if gene == "Potri.001G008400":
		#	print(gene + '\t' + loc + '\t' + str(asTypes) + '\t' + str(lacksAS))
	return(domCheck, lossCheck)

def checkSE(asDict, leftCoord, rightCoord, reLeft, reRight, domLeft, domRight):
	domCheck, lossCheck = 0, 0
	if(leftCoord > domLeft and rightCoord < domRight and domCheck == 0):		
		domCheck = reCheck(seLeft = leftCoord, seRight = rightCoord, reLeft = reLeft, reRight = reRight)	
	elif(leftCoord < domLeft and rightCoord > domRight and domCheck == 0):  ### dom contined in the RE
		domCheck = reCheck(seLeft = leftCoord, seRight = rightCoord, reLeft = reLeft, reRight = reRight)	
	elif(leftCoord > domLeft and rightCoord > domRight and leftCoord < domRight and domCheck == 0):
		domCheck = reCheck(seLeft = leftCoord, seRight = rightCoord, reLeft = reLeft, reRight = reRight)	
	elif(leftCoord < domLeft and rightCoord < domRight and rightCoord > domLeft and domCheck == 0):
		domCheck = reCheck(seLeft = leftCoord, seRight = rightCoord, reLeft = reLeft, reRight = reRight)	
	if(domCheck == 1):
		lossCheck = 1
	return(domCheck, lossCheck)
	
def checkRE(asDict, leftCoord, rightCoord, seLeft, seRight, domLeft, domRight):
	domCheck, lossCheck = 0, 0
	if(leftCoord > domLeft and rightCoord < domRight and domCheck == 0):		
		domCheck = reCheck(reLeft = leftCoord, reRight = rightCoord, seLeft = seLeft, seRight = seRight)	
	elif(leftCoord < domLeft and rightCoord > domRight and domCheck == 0):  ### dom contined in the RE
		domCheck = reCheck(reLeft = leftCoord, reRight = rightCoord, seLeft = seLeft, seRight = seRight)	
	elif(leftCoord > domLeft and rightCoord > domRight and leftCoord < domRight and domCheck == 0):
		domCheck = reCheck(reLeft = leftCoord, reRight = rightCoord, seLeft = seLeft, seRight = seRight)	
	elif(leftCoord < domLeft and rightCoord < domRight and rightCoord > domLeft and domCheck == 0):
		domCheck = reCheck(reLeft = leftCoord, reRight = rightCoord, seLeft = seLeft, seRight = seRight)	
	if(domCheck == 1):
		lossCheck = 1
	return(domCheck, lossCheck)	

def reCheck(reLeft, reRight, seLeft, seRight):
	domCheck = 0
	if(reLeft > seLeft and reRight < seRight):		
		domCheck = 1
	return(domCheck)



	

if __name__ == "__main__":
	main()
