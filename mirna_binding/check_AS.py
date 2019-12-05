#!/usr/bin/env python2.7
def main():
	import io

def checkAS(longTid, asDict, asCoordDict, tidDict, domainsDict, geneDomains, exonDict, gene):
	from collections import defaultdict
	events = defaultdict(dict)
		
	for tid in tidDict:
		events[tid]["gains"] = []
		if "gains" in tidDict[tid]:
			gains = tidDict[tid]["gains"]
			gCoords = []
			for g in gains:
				for gLoc in domainsDict[tid][g]:
					#print(str(domainsDict[tid][g]))
					gCoords = gLoc.split(':')
					gL, gR = int(gCoords[0]), int(gCoords[1])
					for e in asDict[tid]:
						for loc in asDict[tid][e]:
							eCoords = loc.split(':')
							eL, eR = int(eCoords[1]), int(eCoords[2])
							lCheck, rCheck = 0,0
							if eL <= gR:
								lCheck = 1
							if eR >= gL:
								rCheck = 1
							if lCheck == 1 and rCheck == 1:
									desc = e + ':' + loc
									if desc not in events[tid]["gains"]:
										events[tid]["gains"].append(desc)
					if len(events[tid]["gains"]) > 0:
						print(gene + '\t' + tid +'\tgain\t' + g +'\t' + gLoc + '\t' + ','.join(events[tid]["gains"]))
		events[tid]["losses"] = []
		if "losses" in tidDict[tid]:
			losses = tidDict[tid]["losses"]
			lCoords = []
			if longTid in domainsDict:
				longDoms = domainsDict[longTid]
			for l in losses:
				for lLoc in domainsDict[longTid][l]:
					lCoords = lLoc.split(':')
					lL, lR = int(lCoords[0]), int(lCoords[1])
					for e in asDict[tid]:
						for loc in asDict[tid][e]:
							eCoords = loc.split(':')
							eL, eR = int(eCoords[1]), int(eCoords[2])
							lCheck, rCheck = 0,0
							if eL <= lR:
								lCheck = 1
							if eR >= lL:
								rCheck = 1
							if lCheck == 1 and rCheck == 1:
								desc = e + ':' + loc
								if desc not in events[tid]["losses"]:
									events[tid]["losses"].append(desc)
					if len(events[tid]["losses"]) > 0:
						print(gene + '\t' + tid +'\tloss\t' + l +'\t' + lLoc + '\t' + ','.join(events[tid]["losses"]))
								

if __name__ == "__main__":
	main()
