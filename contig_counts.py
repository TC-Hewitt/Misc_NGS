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
        binLims = (0, 500, 1000, 2000, 5000, 8000, 12000, 20000, 40000, 60000, 100000, 200000, 500000, 1000000, 3000000, 6000000, 10000000, 20000000, float('inf'))
        binCounts = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0, 17:0}
        binDex = ('0-0.5Kb', '0.5-1Kb', '1-2Kb', '2-5Kb', '5-8Kb', '8-12Kb', '12-20Kb', '20-40Kb', '40-60Kb', '60-100Kb', '100-200Kb', '200-500Kb', '0.5-1Mb', '1-3Mb', '3-6Mb', '6-10Mb', '10-20Mb', '20Mb+')
    else:
        bins = False
    # Open FASTA.
    with open(args.input, 'r') as fastaIn:
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
                    while i < 18:
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
            while i < 18:
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
    countL = 0
    L25 = ''
    L50 = ''
    L75 = ''
    for num in ctgLens:
        rollSum += num
        countL += 1
        if check25 == False and rollSum >= len25:
            N25 = str(num)
            L25 = str(countL)
            check25 = True
        if check50 == False and rollSum >= len50:
            N50 = str(num)
            L50 = str(countL)
            check50 = True
        if rollSum >= len75:
            N75 = str(num)
            L75 = str(countL)
            break


    print('# for ' + args.input + ':\n# number of contigs = ' + str(ctgCount) + '\n# combined length = ' + str(totalLen) + 'bp\n# max length = ' + str(ctgMax) + 'bp\n# min length = ' + str(ctgMin) + 'bp\n# average length = ' + str(ctgAvg) + 'bp\n# median length = ' + str(ctgMed) + 'bp\n# SD = ' + str(ctgStd) + '\n# Ns = ' + str(Ns) + '\n# N25 = ' + N25 + '\n# N50 = ' + N50 + '\n# N75 = ' + N75 + '\n# L25 = ' + L25 + '\n# L50 = ' + L50 + '\n# L75 = ' + L75 + '\n#')
    if bins:
        print('# length ranges:')
        for i in range(18):
            print('# ' + str(binDex[i]) + ' = ' + str(binCounts[i]))
        print('\n')

if __name__ == '__main__':
    main()
