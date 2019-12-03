#!/usr/bin/env python2.7

def main():
	args = getArgs()
	dat, sqtl = args.dat, args.sqtl
	dat  = storeDat(dat = dat)
	getEvents(dat = dat, sqtl = sqtl)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = 'clean up pasa indiv_splice_labels_and_coords file such that A3 and A5 have genomic coords of splice junctions')
	parser.add_argument('-dat', action = 'store', type = str, required = True, help = 'pasa .dat file')
	parser.add_argument('-sqtl', action = 'store', type = str, required = True, help = 'sqtl output from matrix etql')
	args = parser.parse_args()
	return(args)

def storeDat(dat):
	import re,io
	from collections import defaultdict
	datDict = defaultdict(dict)
	fh = open(dat, "r")
	for line in fh:
		line = line.strip()
		arr = line.split('\t')
		chromo = arr[0]
		asType = arr[1]
		if asType == "retained_exon" or asType == "spliced_intron":
			continue
		l,r = arr[2], arr[3]
		coords = l +':' + r
		if chromo in datDict:
			if coords in datDict[chromo]:
				if asType not in datDict[chromo][coords]:
					datDict[chromo][coords].append(asType)
			else:
				datDict[chromo][coords] = [asType]
		else:
			datDict[chromo][coords] = [asType]
	return(datDict)

def getEvents(dat, sqtl):
        import re,io
        from collections import defaultdict
        fh = open(sqtl, "r")
        for line in fh:
			line = line.strip()
			if re.search('snps',line):
				continue
			arr = line.split('\t')
			chromo = arr[0].split('_')[0]
			junc = arr[1].split(':')
			junc = junc[1] + ':' + junc[2]
			if junc in dat[chromo]:
				print(arr[0] + '\t' + arr[1] + '\t' + ','.join(dat[chromo][junc]))




if __name__ == "__main__":
	main()
