#!/usr/bin/env python2.7

def main():
	args = getArgs()
	dat, junc, tmap = args.dat, args.junc, args.tmap
	dat  = storeDat(dat = dat)
	tmap, geneDict = storeTmap(tmap = tmap)
	getEvents(dat = dat, junc = junc, geneDict = geneDict, tmap = tmap)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = 'clean up pasa indiv_splice_labels_and_coords file such that A3 and A5 have genomic coords of splice junctions')
	parser.add_argument('-dat', action = 'store', type = str, required = True, help = 'pasa .dat file')
	parser.add_argument('-junc', action = 'store', type = str, required = True, help = 'matrix of junction usage in population')
	parser.add_argument('-tmap', action = 'store', type = str, required = True, help = 'transmap')
	args = parser.parse_args()
	return(args)

def storeTmap(tmap):
	import io,re,sys
	from collections import defaultdict
	fh = open(tmap, 'r')
	geneDict = defaultdict(dict)
	tmap = defaultdict(dict)
	for line in fh:
		line = line.strip()
		arr = line.split()
		gene, tid = arr[0], arr[1]
		geneDict[tid] = gene
		if gene in tmap:
			tmap[gene].append(tid)
		else:
			tmap[gene] = [tid]
	return(tmap, geneDict)


def storeDat(dat):
	import re,io
	from collections import defaultdict
	datDict = defaultdict(dict)
	fh = open(dat, "r")
	for line in fh:
		line = line.strip()
		arr = line.split('\t')
		chromo = arr[0]
		asType = arr[3]
		tid = arr[7]
		if asType == "retained_exon" or asType == "spliced_intron":
			continue
		l,r = arr[4], arr[5]
		coords = l +':' + r
		if tid in datDict:
			if coords in datDict[tid]:
				if asType not in datDict[tid][coords]:
					datDict[tid][coords].append(asType)
			else:
				datDict[tid][coords] = [asType]
		else:
			datDict[tid][coords] = [asType]
	return(datDict)

def getEvents(dat, junc, geneDict, tmap):
        import re,io
        from collections import defaultdict
        fh = open(junc, "r")
        for line in fh:
			line = line.strip()
			if re.search('UF',line):
				print(line)
				continue
			arr = line.split(',')
			gene = arr[0].split(':')[0]
			junc = arr[0].split(':')
			junc = junc[1] + ':' + junc[2]
			tids = tmap[gene]
			for tid in tids:
				if tid in dat:
					#print(str(dat[tid]))
					if junc in dat[tid]:
						print(line)
						break




if __name__ == "__main__":
	main()
