#!/usr/bin/env python

import argparse
from functools import partial
from operator import itemgetter
import six


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__file__, description = "Taxonomic Ranks calculation")
    parser.add_argument('-quantFile', '-qf', help = 'Salmon quant.sf file/cedar quant file/ or kallisto', default = None)
    parser.add_argument('-lineageFile', '-lf', help = 'Custom taxonomy file directory', default = None)
    parser.add_argument('-outDir', '-od', help = 'Output directory', default = None)
    parser.add_argument('-dataType', '-dt', help = 'Quantification file type: salmon, kallisto, cedar', default = 'salmon')

args = parser.parse_args()


lineage = list()
salmonCOUNTS = list()
salmonTAXID = list()
results = dict()


def binary_search(L, target):
    start = 0
    end = len(L) - 1

    while start <= end:
    	middle = (start + end) // 2
    	midpoint = L[middle]
    	if midpoint > target:
    		end = middle - 1
    	elif midpoint < target:
    		start = middle + 1
    	else:
    		return middle


with open(args.lineageFile) as f:
    for line in f:
        el = line.split('|')
        if el[9].strip() == 'Eukaryota':
            continue
        else:
            if el[2].strip() == '':
                lineage.append([el[0].strip(), el[9].strip(), el[7].strip(), el[6].strip(), el[5].strip(), \
                        el[4].strip(), el[3].strip(), el[1].strip(), el[1].strip()])
            else:
                lineage.append([el[0].strip(), el[9].strip(), el[7].strip(), el[6].strip(), el[5].strip(), \
                        el[4].strip(), el[3].strip(), el[2].strip(), el[1].strip()])


with open(args.quantFile) as f:
    next(f)
    for line in f:
        elements = line.split('\t')
        salmonTAXID.append(elements[0].strip().split('|')[2])
        if args.dataType == 'cedar':
            salmonCOUNTS.append(float(elements[2].strip()))
        elif args.dataType == 'kallisto':
            salmonCOUNTS.append(float(elements[3].strip()))
        else:
            salmonCOUNTS.append(float(elements[4].strip()))


noMatch = open(args.outDir + '/noMatch.tab', 'w')


sorted_lineage = sorted(lineage, key = itemgetter(0))
tmp_list = [item[0] for item in sorted_lineage]

sorted_newl = [[el or '---' for el in sublist] for sublist in sorted_lineage]

for i in range(0, len(salmonTAXID)):
    index = binary_search(tmp_list, salmonTAXID[i])
    if index is not None:
        if sorted_newl[index][0] in results.keys():
            results[sorted_newl[index][0]][8] += salmonCOUNTS[i]
        else:
            results[sorted_newl[index][0]] = [sorted_newl[index][1], sorted_newl[index][2], \
                                             sorted_newl[index][3], sorted_newl[index][4], sorted_newl[index][5], \
                                            sorted_newl[index][6], sorted_newl[index][7], sorted_newl[index][8], salmonCOUNTS[i]]
    else:
        noMatch.write('\t'.join(['%s', '%f']) % (salmonTAXID[i], salmonCOUNTS[i]))



match = open(args.outDir + '/out.tab', 'w')
match.write('\t'.join(['Superkigndom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Scietific_Name', 'TaxID', 'NumCounts']))
match.write('\n\n')

for k, v in six.iteritems(results):
    match.write('\t'.join(['%s','%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%f']) % \
    			(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], k, v[8]))
    match.write('\n')


