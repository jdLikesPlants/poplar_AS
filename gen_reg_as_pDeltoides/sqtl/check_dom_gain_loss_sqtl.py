#!/usr/bin/env python2.7

def main():
	args = getArgs()
	dom, sqtl = args.dom, args.sqtl
	sqtl  = storeSqtl(sqtl = sqtl)
	checkEvents(dom = dom, sqtl = sqtl)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = 'clean up pasa indiv_splice_labels_and_coords file such that A3 and A5 have genomic coords of splice junctions')
	parser.add_argument('-dom', action = 'store', type = str, required = True, help = 'domain gain/loss file')
	parser.add_argument('-sqtl', action = 'store', type = str, required = True, help = 'sqtl output from matrix etql')
	args = parser.parse_args()
	return(args)

def storeSqtl(sqtl):
	import re,io
	from collections import defaultdict
	sqtlDict = defaultdict(dict)
	fh = open(sqtl, "r")
	for line in fh:
		line = line.strip()
		arr = line.split('\t')
		if re.search("snps", line):
			continue
		arr = arr[1].split(':')
		gene = arr[0]
		coords = arr[1] + ':' + arr[2]
		if gene in sqtlDict:
			if coords not in sqtlDict[gene]:
				sqtlDict[gene].append(coords)
		else:
			sqtlDict[gene] = [coords]
	return(sqtlDict)

def checkEvents(sqtl, dom):
        import re,io
        from collections import defaultdict
        fh = open(dom, "r")
        for line in fh:
			line = line.strip()
			if re.search('dominant_iso',line):
				continue
			gene = line.split()[0]
			gl = line.split()[2]
			domain = line.split('\t')[3]
			arr = line.split('\t')[-1].split(',')
			for i in arr:
				j = i.split(':')			
				asType = j[0]
				junc = j[2] + ':' + j[3]
				if gene in sqtl:
					if junc in sqtl[gene]:
						print(gene + '\t' + asType + '\t' + gl + '\t' + domain + '\t' +junc)
					



if __name__ == "__main__":
	main()
