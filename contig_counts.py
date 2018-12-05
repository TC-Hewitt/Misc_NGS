#!/usr/bin/env python

import argparse
import numpy as np
import re

def main():
    
    # Parse arguments.
    parser = argparse.ArgumentParser(description='get contig counts and stats from a multi-fasta file')
    parser.add_argument('-i', '--input', help='indicate input fasta', required=True)
    parser.add_argument('-s', '--show', help='print id of each contig above this length. Indicate min length', type=int, required=False)
    parser.add_argument('-b', '--bins', help='indicate True to print binned contig lengths', required=False)
    args = parser.parse_args()
    
    # establish params
    ctgLens = []
    p = re.compile('(>.+\n)|\W')
    N = re.compile('[nN]')
    Ns = 0
    if args.show:
        show = True
    else:
        show = False
    if args.bins in ['T', 't', 'True', 'true', 'TRUE']:
        bins = True
        binLims = (0, 500, 1000, 2000, 5000, 8000, 12000, 20000, 40000, 60000, 100000, 200000, float('inf'))
        binCounts = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0}
        binDex = ('0-0.5kb', '0.5-1kb', '1-2kb', '2-5kb', '5-8kb', '8-12kb', '12-20kb', '20-40kb', '40-60kb', '60-100kb', '100-200kb', '200kb+')
    else:
        bins = False
    # Open FASTA.
    with open(args.input, 'rU') as fastaIn:
    	tempSeq = ''
        ctgid = 'null'
	for line in fastaIn:
    	    if line.startswith('>'):
                seqLen = len(tempSeq)
                if show and seqLen >= args.show:
                    print(ctgid.strip('>') + '\t' + str(seqLen))
                ctgid = line.strip('\n')
                if bins:
                    i = 0
                    while i < 12:
                        if binLims[i] < seqLen <= binLims[i+1]:
                            binCounts[i] += 1
                            break
                        i += 1
                ctgLens.append(seqLen)
    		Ns += len(N.findall(tempSeq))
                tempSeq = ''
            else:
                tempSeq += p.sub('', line)
    	ctgLens.append(len(tempSeq))
        Ns += len(N.findall(tempSeq))
        if show and seqLen >= args.show:
            print(ctgid.strip('>') + '\t' + str(len(tempSeq)) + '\n')
        if bins:
            i = 0
            while i < 12:
                if binLims[i] < len(tempSeq) <= binLims[i+1]:
                    binCounts[i] += 1
                    break
                i += 1
    	del(tempSeq)
    
    ctgLens.pop(0)
    ctgLens.sort(reverse=True)
    ctgCount = len(ctgLens)
    totalLen = sum(ctgLens)
    ctgMax = max(ctgLens)
    ctgMin = min(ctgLens)
    ctgAvg = int(round(np.mean(ctgLens)))
    ctgMed = int(round(np.median(ctgLens)))
    ctgStd = round(np.std(ctgLens), 2)
    len25 = totalLen/4
    check25 = False
    N25 = ''
    len50 = len25*2
    check50 = False
    N50 = ''
    len75 = len25*3
    N75 = ''
    rollSum = 0
    for num in ctgLens:
        rollSum += num
        if check25 == False and rollSum >= len25:
            N25 = str(num)
            check25 = True
        if check50 == False and rollSum >= len50:
            N50 = str(num)
            check50 = True
        if rollSum >= len75:
            N75 = str(num)
            break


    print('for ' + args.input + ':\nnumber of contigs = ' + str(ctgCount) + '\ncombined length = ' + str(totalLen) + 'bp\nmax length = ' + str(ctgMax) + 'bp\nmin length = ' + str(ctgMin) + 'bp\naverage length = ' + str(ctgAvg) + 'bp\nmedian length = ' + str(ctgMed) + 'bp\nSD = ' + str(ctgStd) + '\nNs = ' + str(Ns) + '\nN25 = ' + N25 + '\nN50 = ' + N50 + '\nN75 = ' + N75 + '\n')
    if bins:
        print('length ranges:')
        for i in range(12):
            print(str(binDex[i]) + ' = ' + str(binCounts[i]))
        print('\n')

if __name__ == '__main__':
    main()
