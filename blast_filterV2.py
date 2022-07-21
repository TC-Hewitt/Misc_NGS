# Copyright (C) 2017 Timothy C. Hewitt - All Rights Reserved
# You may use, distribute and modify this code under the terms of the GNU Public License version 3 (GPLv3)
# You should have recieved a copy of the GPLv3 license with this file. If not, please visit https://github.com/TC-Hewitt/MuTrigo

#!/usr/bin/env python


import argparse, csv, sys
from operator import itemgetter

def main():

    parser = argparse.ArgumentParser(description='filter and sort a blastn output (outfmt 6 or 7) based on the cutoffs you want for each parameter. Assumes standard out format with addition of qcovs and qcovhsp occupying fields 13 and 14, but these can be substituted for other output options in the blastn command.')
    parser.add_argument('-i', '--input', help='indicate input file (as blast tab outfmt 6 or 7)', required=True)
    parser.add_argument('-o', '--output', help='indicate output file', required=True)
    parser.add_argument('-p', '--pcntid', help='set min percent id', default=0, type=float, required=False)
    parser.add_argument('-a', '--alnlen', help='set min alignment length', default=0, type=int, required=False)
    parser.add_argument('-m', '--msmtch', help='set max mismatches', default=float('Inf'), type=int, required=False)
    parser.add_argument('-g', '--gpopen', help='set max gap opens', default=float('Inf'), type=int, required=False)
    parser.add_argument('-qs', '--qstart', help='set min query seq start', default=0, type=int, required=False)
    parser.add_argument('-qe', '--q_end', help='set max query seq end', default=float('Inf'), type=int, required=False)
    parser.add_argument('-ss', '--sstart', help='set min subject seq start', default=0, type=int, required=False)
    parser.add_argument('-se', '--s_end', help='set max subject seq end', default=float('Inf'), type=int, required=False)
    parser.add_argument('-e', '--evalue', help='set max evalue', default=0.01, type=float, required=False)
    parser.add_argument('-b', '--bscore', help='set min bit score', default=0, type=float, required=False)
    parser.add_argument('-qcov', '--qcov', help='set min qcov', default=0, type=int, required=False)
    parser.add_argument('-qcovh', '--qcovhsp', help='set min qcovhsp', default=0, type=int, required=False)
    parser.add_argument('-sort1', '--sort1', help='choose primary paramter to sort rows best to worst by entering string value: <qname|sname|pcntd|alnlen|msmtch|gpopen|qstart|q_end|sstart|s_end|evalue|bscore|qcov|qcovhsp>', type=str, required=False)
    parser.add_argument('-sort2', '--sort2', help='choose secondary paramter to sort rows best to worst by entering string value: <qname|sname|pcntd|alnlen|msmtch|gpopen|qstart|q_end|sstart|s_end|evalue|bscore(default)|qcov|qcovhsp>', type=str, required=False)
    parser.add_argument('-topq', '--gettopq', help='set number of top hits per query to output. Top hits based on preferred sort parameter - must use sort option', type=int, required=False)
    parser.add_argument('-tops', '--gettops', help='set number of top hits per subject to output. Top hits based on preferred sort parameter - must use sort option', type=int, required=False)
    args = parser.parse_args()

    '''based on blast+ outfmt 6 or 7. Each row split into python list (delim = \t). Each index reps following values (default after 'set'):
        0=query id (query name)
        1=subject id (subject name)
        2=% identity (int(row[2])) set min >= 0
        3=align length (int(row[3])) set min >= 0
        4=mismatches (int(row[4])) set max <= float('Inf')
        5=gap opens (int(row[5])) set max <= float('Inf')
        6=q. start (int(row[6])) set min >= 0
        7=q. end (int(row[7])) set max <= float('Inf')
        8=s. start (int(row[8])) set min >= 0
        9=s. end (int(row[9])) set max <= float('Inf')
        10=evalue (float(row[10])) set max <= float(0.01)
        11=bit score (int(row[11])) set min >= 0
        12=qcov (int(row[12])) set min >= 0
        13=qcovhsp (int(row[13])) set min >= 0'''
    
    coldex = {'qname':(0, False), 'sname':(1, False), 'pcntd':(2,True), 'alnlen':(3,True), 'msmtch':(4,False), 'gpopen':(5,False), 'qstart':(6,False), 'q_end':(7,False), 'sstart':(8,False), 's_end':(9,False), 'evalue':(10,False), 'bscore':(11,True), 'qcov':(12,True), 'qcovhsp':(13,True)}
    if args.sort1:
        param1 = coldex[args.sort1] # col to sort is param1[0], ascend or descend is param1[1]
        sortLst = []
    else:
        print('Primary sort not given - output will be unsorted.')

    if args.sort1 and args.sort2:
        param2 = coldex[args.sort2] # col to sort is param2[0], ascend or descend is param2[1]
    else:
        param2 = coldex['bscore']

    if args.gettopq and args.gettops:
        sys.exit('abort: cannot run both -topq and -tops. Choose one.')

    with open(args.input, 'r') as file_in:
        print('Parsing ' + args.input + '...')
        file_out = open(args.output, 'w+')
        file_out.write('#This output based on arguments: \n' + '#input file=' + args.input + '\n#output file=' + args.output + '\n#min percent id=' + str(args.pcntid) + '\n#min alignment length=' + str(args.alnlen) + '\n#max mismatches=' + str(args.msmtch) + '\n#max gap opens=' + str(args.gpopen) + '\n#min query seq start=' + str(args.qstart) + '\n#max query seq end=' + str(args.q_end) + '\n#min subject seq start=' + str(args.sstart) + '\n#max subject seq end=' + str(args.s_end) + '\n#max evalue=' + str(args.evalue) + '\n#min bit score=' + str(args.bscore) + '\n#min pcnt query cov=' + str(args.qcov) + '\n#min pcnt query cov hsp=' + str(args.qcovhsp) + '\n')
        file_out.write('#Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, qcov, qcovhsp \n')
        count = 0
        reader = csv.reader(file_in, delimiter = '\t')
        for row in reader:
            if '#' in row[0]:
                continue
            idxError = 0
            while idxError <= 13:
                try:
                    if all([float(row[2])>=args.pcntid, int(row[3])>=args.alnlen, int(row[4])<=args.msmtch, int(row[5])<=args.gpopen, int(row[6])>=args.qstart, int(row[7])<=args.q_end, int(row[8])>=args.sstart, int(row[9])<=args.s_end, float(row[10])<=args.evalue, float(row[11])>=args.bscore, int(row[12])>=args.qcov, int(row[13])>=args.qcovhsp]):
                        if args.sort1:
                            sortLst.append(row)
                        else:
                            file_out.write(str('\t'.join(row)) + '\n')
                        count += 1
                except IndexError:
                    row.append('0')
                    idxError += 1
                    continue
                break
        
        if args.sort1:
            def force_type(x):
                try:
                    return(int(x))
                except ValueError:
                    try:
                        return(float(x))
                    except ValueError:
                        return(x)

            sortLst = [[force_type(x) for x in row] for row in sortLst]
            sortLst.sort(key=itemgetter(param2[0]), reverse=param2[1])
            sortLst.sort(key=itemgetter(param1[0]), reverse=param1[1])

            if args.gettopq:
                qnames = {}
                tally = 0
                topLst = []
                if args.sort1 == 'qname':
                    for row in sortLst:
                        if row[0] not in qnames.keys():
                            qnames[row[0]] = 1
                            file_out.write(str('\t'.join(str(x) for x in row)) + '\n')
                            tally += 1
                        elif qnames[row[0]] < args.gettopq:
                            qnames[row[0]] += 1
                            file_out.write(str('\t'.join(str(x) for x in row)) + '\n')
                            tally += 1
                        else:
                            continue
                    print('Wrote %d top hits out of %d matches to %s' %(tally, count, args.output))
                else:
                    for row in sortLst:
                        if row[0] not in qnames.keys():
                            qnames[row[0]] = 1
                            topLst.append(row)
                            tally += 1
                        elif qnames[row[0]] < args.gettopq:
                            qnames[row[0]] += 1
                            topLst.append(row)
                            tally += 1
                        else:
                            continue
                    if args.gettopq > 1:
                        topLst.sort(key=itemgetter(0), reverse=False)
                    for row in topLst:
                        file_out.write(str('\t'.join(str(x) for x in row)) + '\n')
                    print('Wrote %d top hits out of %d matches to %s' %(tally, count, args.output))
            elif args.gettops:
                snames = {}
                tally = 0
                topLst = []
                if args.sort1 == 'sname':
                    for row in sortLst:
                        if row[1] not in snames.keys():
                            snames[row[1]] = 1
                            file_out.write(str('\t'.join(str(x) for x in row)) + '\n')
                            tally += 1
                        elif snames[row[1]] < args.gettops:
                            snames[row[1]] += 1
                            file_out.write(str('\t'.join(str(x) for x in row)) + '\n')
                            tally += 1
                        else:
                            continue
                    print('Wrote %d top hits out of %d matches to %s' %(tally, count, args.output))
                else:
                    for row in sortLst:
                        if row[1] not in snames.keys():
                            snames[row[1]] = 1
                            topLst.append(row)
                            tally += 1
                        elif snames[row[1]] < args.gettops:
                            snames[row[1]] += 1
                            topLst.append(row)
                            tally += 1
                        else:
                            continue
                    if args.gettops > 1:
                        topLst.sort(key=itemgetter(1), reverse=False)
                    for row in topLst:
                        file_out.write(str('\t'.join(str(x) for x in row)) + '\n')
                    print('Wrote %d top hits out of %d matches to %s' %(tally, count, args.output))
            else:
                for row in sortLst:
                    file_out.write(str('\t'.join(str(x) for x in row)) + '\n')
                print('Wrote %d matches to %s' %(count, args.output))
        else:
            print('Wrote %d matches to %s' %(count, args.output))
        file_out.close()
        

if __name__ == '__main__':
    main()

