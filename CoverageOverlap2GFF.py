#!/usr/bin/env python

import argparse, sys, csv, re
import numpy as np
from os.path import exists

def main():

    ## Parse arguments ##
    parser = argparse.ArgumentParser(description='Compare read coverage from BAM/SAM file (requires samtools depth output) against gff annotation. \
        Generates gff output of coverage intervals both non-overlapping and overlapping with gff features \
        (strand info not retained, attributes: ID=seqid_0start|overlapping gff ID;mean cov;median cov;num gaps;overlap type;overlap len;overlap %;).')
    parser.add_argument('-i', '--indepths', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='samtools depth single sample output file (leave out if using STDIN).')
    parser.add_argument('-g', '--gffref', required=False, help='reference annotation gff file (optional). Only mRNA or gene lines preferred (see option -t). Must contain ID attribute.')
    parser.add_argument('-t', '--typestr', type=str, required=False, help='indicate feature type (as string) to retrieve from reference gff, e.g. mRNA (default=none).')
    parser.add_argument('-m', '--mindep', type=int, required=False, default=2, help='min mean read depth to consider as feature (default=2).')
    parser.add_argument('-d', '--depnuc', type=int, required=False, default=2, help='min per nucleotide depth. Ignore bases below this depth (default=2).')
    parser.add_argument('-l', '--featlen', type=int, required=False, default=100, help='min length to consider as feature (default=100).')
    parser.add_argument('-G', '--gaplen', type=int, required=False, default=20, help='min length of ignored bases to count as gap (default=20).')
    parser.add_argument('-j', '--jumplen', type=int, required=False, default=500, help='max gap length allowed to extend feature, i.e. introns (default=500).')
    parser.add_argument('-o', '--outname', type=str, required=False, default="out", help='specify output prefix (default="out")')
    parser.add_argument('-s', '--showtypes', help='indicate TRUE to show overlap types.')
    args = parser.parse_args()

    ## check input arguments ##
    if args.gffref and not exists(args.gffref):
        sys.exit("Error: " + args.gffref + " not found!\n")

    if args.showtypes in ['T', 't', 'True', 'true', 'TRUE']:
        print('''
    overlap types (1= coverage interval, 2= gff feature):

    null)    1------------------>
    


    tail)    1------------------>
                        2------------------>


    head)               1------------------>
             2------------------>


    inset)           1-------->
                2------------------>


    flank)      1------------------>
                     2-------->


    exact)      1------------------>
                2------------------>
        
    ''')

    fileOut = args.outname + '.gff'

    ## parse depths and record feats to data struct ##
    indepths = csv.reader(args.indepths, delimiter = '\t', quoting=csv.QUOTE_NONE)
    line1 = next(indepths)
    ncol = len(line1) # count num cols
    if ncol != 3:
        sys.exit("Error: depth input has incorrect number of columns!")
    else:
        seqid = line1[0]
        runstart = int(line1[1])
        runend = int(line1[1])
        depths = [int(line1[2])]
        gaps = 0
        covList = [] # coverage feats stored as list of tuples
        for row in indepths:
            if row[2] == '0' or int(row[2]) < args.depnuc: #lines with depth 0 (if samtools depth -a/aa is used) or less than min nt depth ignored
                continue
            jump = int(row[1]) - runend
            if row[0] == seqid and jump <= args.jumplen: #maintain feature interval if seqname same, jump less than cutoff
                seqid = row[0]
                runend = int(row[1])
                depths.append(int(row[2]))
                if jump >= args.gaplen: #count as gap if jump len exceeds gap threshold
                	gaps += 1
            else: #do calcs and store expression feature and start new one - reset params when deviate beyond cutoffs
                runlen = runend - runstart
                depmean = int(round(np.mean(depths)))
                depmedian = int(np.median(depths))
                if runlen >= args.featlen and depmean >= args.mindep:
                    covList.append((seqid, runstart, runend, depmean, depmedian, gaps))
                seqid = row[0]
                runstart = int(row[1])
                runend = int(row[1])
                depths = [int(row[2])]
                gaps = 0
        covList.append((seqid, runstart, runend, depmean, depmedian, gaps))

        print("Coverage features recorded.")

    ## write out coverage feats as is if no ref gff supplied ##
    if not args.gffref:
        print("Reference gff not supplied. Writing coverage features to " + args.outname + ".gff")
        output = open(fileOut, 'w+')
        for feat in covList:
            output.write(feat[0] + "\t.\tcDNA\t" + str(feat[1]) + "\t" + str(feat[2]) + "\t.\t+\t.\tID=" + \
                feat[0] + "_" + str(feat[1]-1) + ";mean=" + str(feat[3]) + ";median=" + str(feat[4]) + ";gaps=" + str(feat[5]) + ";\n")
        print("Finished.\n\nn coverage features: " + str(len(covList)))
        sys.exit()

    ## if ref gff supplied, write gff feats to data struct ##
    gffDict = {}
    gffFile = open(args.gffref, 'r')
    gffIn = csv.reader(gffFile, delimiter = '\t', quoting=csv.QUOTE_NONE)
    for row in gffIn:
        try:
            if args.typestr and args.typestr not in row[2]: #skip row if feature field does not match specified type
                continue
            featid = re.match(r"ID=([^;]+)",row[8]).group(1)
            if row[0] in gffDict:
                gffDict.setdefault(row[0], []).append((int(row[3]),int(row[4]), featid))
            else:
                gffDict[row[0]]=[(int(row[3]),int(row[4]), featid)]
        except (IndexError, ValueError):
            continue
        #generates something like: {'contig_1':[(210,510,'gene_1'),(1215,3211,'gene_2')],'contig_2':[(123,456,'gene_3'),(789,1112,'gene_4'),...}
    gffFile.close()

    ## compare cov feats to gff feats to find overlaps and write out cov feats with modified attributes ##
    print("Comparing coverage intervals to gff features...")
    output = open(fileOut, 'w+')
    #setup overlap types tally
    null = 0
    tail = 0
    head = 0
    inset = 0
    flank = 0
    exact = 0
    for feat in covList: #iterate over cov feats
        #feat=(seqid, runstart, runend, depmean, depmedian)
        if feat[0] in gffDict.keys(): #check if feat id in gff dict, otherwise write out cov feat as is (null type) 
            featLen = float(feat[2] - feat[1])
            attID = feat[0] + "_" + str(feat[1]-1) #unless overlap with gff feat is found
            anymatch = False
        else:
            null += 1
            #attributes: ID=seqid_seqstart;mean cov;median cov;overlap type;
            output.write(feat[0] + "\t.\tcDNA\t" + str(feat[1]) + "\t" + str(feat[2]) + "\t.\t+\t.\tID=" + \
                feat[0] + "_" + str(feat[1]-1) + ";mean=" + str(feat[3]) + ";median=" + str(feat[4]) + ";gaps=" + str(feat[5]) + ";type=null;\n")
            continue
        for (start,end,featid) in gffDict[feat[0]]: #per cov feat, iterate over gff feats with matching seqid comparing coords
                overlapLen = 0
                overlapPcnt = 0.0
                overlapType = ''
                if start < feat[2] <= end and feat[1] < start: #overlap type= tail
                    tail += 1
                    overlapLen = feat[2] - start
                    overlapPcnt = round((overlapLen/featLen)*100, 1)
                    overlapType = 'tail'
                    anymatch = True
                elif start <= feat[1] < end and feat[2] > end: #overlap type= head
                    head += 1
                    overlapLen = end - feat[1]
                    overlapPcnt = round((overlapLen/featLen)*100, 1)
                    overlapType = 'head'
                    anymatch = True
                elif start == feat[1] and feat[2] == end: #overlap type= exact
                    exact += 1
                    overlapLen = int(featLen)
                    overlapPcnt = 100.0
                    overlapType = 'exact'
                    anymatch = True
                elif start <= feat[1] and feat[2] <= end: #overlap type= inset
                    inset += 1
                    overlapLen = feat[2] - feat[1]
                    overlapPcnt = 100.0
                    overlapType = 'inset'
                    anymatch = True
                elif start > feat[1] and feat[2] > end: #overlap type= flank
                    flank += 1
                    overlapLen = end - start
                    overlapPcnt = round((overlapLen/featLen)*100, 1)
                    overlapType = 'flank'
                    anymatch = True
                else:
                    continue
                if anymatch == True: #overlap found, write out feat with attributes: ID=featid;cov parent;mean cov;median cov;n gaps;overlap type;overlap len;overlap %;)
                    attID = featid #from gff line
                    output.write(feat[0] + "\t.\tcDNA\t" + str(feat[1]) + "\t" + str(feat[2]) + "\t.\t+\t.\tID=" + \
                        attID + ";Parent=" + feat[0] + "_" + str(feat[1]-1) + ";mean=" + str(feat[3]) + ";median=" + str(feat[4]) + \
                        ";gaps=" + str(feat[5]) + ";type=" + overlapType + ";len=" + str(overlapLen) + ";pcnt=" + str(overlapPcnt) + ";\n")

        if anymatch == False: #no overlaps with feat found, write out cov feat as is (null type)
            #attributes: ID=seqid_seqstart;mean cov;median cov;overlap type;
            null += 1
            output.write(feat[0] + "\t.\tcDNA\t" + str(feat[1]) + "\t" + str(feat[2]) + "\t.\t+\t.\tID=" + \
                feat[0] + "_" + str(feat[1]-1) + ";mean=" + str(feat[3]) + ";median=" + str(feat[4]) + ";gaps=" + str(feat[5]) + ";type=null;\n")

    ## print final tallies to console and end ##
    print("Finished.\n\nn coverage features: " + str(len(covList)) + "\noverlaps with gff features:\nnull=" + str(null) + "\ntail=" + \
        str(tail) + '\nhead=' + str(head) + '\ninset=' + str(inset) + '\nflank=' + str(flank) + '\nexact=' + str(exact))

if __name__ == '__main__':
    main()
