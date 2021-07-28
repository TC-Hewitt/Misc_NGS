#!/usr/bin/env python


import argparse
import re

def main():

    # Parse arguments.
    parser = argparse.ArgumentParser(description='capture all contigs user defined by numbers')
    parser.add_argument('-i', '--input', help='indicate input.fasta', required=True)
    parser.add_argument('-o', '--output', help='indicate output.fasta', required=True)
    parser.add_argument('-l', '--list', help='input space delim list of contig numbers', nargs='*', required=False)
    parser.add_argument('-t', '--table', help='provide txt file or blast tab output containing a list of contigs', required=False)
    parser.add_argument('-s', '--strid', help='identifier string before contig number. Default is <contig_> (only compatible if using -t)', default='contig_', type=str, required=False)
    parser.add_argument('-g', '--getnames', help='enter output file name to get a list of non-redundant contig names from <-t> input', required=False)
    parser.add_argument('-v', '--invert', help='indicate True to retrieve contigs NOT in list provided by --table', default=False, type=bool, required=False)
    args = parser.parse_args()

    # Open FASTA.
    fasta_in = open(args.input, 'rU')
    count = 0
    matched = 0

    if args.list and not args.table:
        fasta_out = open(args.output, 'w+')
        
        # append user defined numbers with regex anchors
        query_str = '(\\D' + '\\s)|(\\D'.join(args.list) + '\\s)'
        print('searching for contigs in %s with exact values inputted: \n %s' % (args.input, query_str))

        # Parse file and write to output.

        for line in fasta_in:
           if line.startswith('>'):
	       if re.search(query_str, line):
                   print('found')
                   fasta_out.write(line)
                   count += 1
                   matched = 1
               else:
                   matched = 0
           elif matched == 1:
               fasta_out.write(line)
        
        # Finish.
        fasta_out.close()
        print('Wrote %d contigs found of %d inputs to %s.' % (count, len(args.list), args.output))

    elif args.table and not args.list:
        table_in = open(args.table, 'rU')
        fasta_out = open(args.output, 'w+')
        print('populating list of queries...')

        #retrieve contigs from fasta_in based on list in args.table
        search_term = '(' + args.strid + '\\w+)[\.\\s,|>:;-]'
        query_get = re.findall(search_term, table_in.read())
	query_raw = [re.sub(args.strid, '', x, count=1) for x in query_get]
	query_set = list(set(query_raw))
        if args.getnames:
            namesOut = open(args.getnames, 'w+')
            [namesOut.write(args.strid + i + '\n') for i in sorted(query_set)]
            namesOut.close()
        num_queries = len(query_set)
	duplicates = len(query_raw) - len(query_set)
        query_hash = set(query_set)
        print(str(duplicates) + ' duplicates removed from list')  
	print(str(num_queries) + ' unique contigs identified from ' + args.table)
        table_in.close()
        if not args.invert:
            print('searching for contigs in %s that match to contents of %s' % (args.input, args.table))
        else:
            print('searching for contigs in %s that do not match to contents of %s' % (args.input, args.table))
	p = re.compile(args.strid + '(\w+)[\.\\s,|]')
        
        for line in fasta_in:
            if line.startswith('>'):
		matched = 0
                try:
                    if not args.invert and p.search(line).group(1) in query_hash:
                        fasta_out.write(line)
                        count += 1
                        matched = 1
                    elif args.invert and p.search(line).group(1) not in query_hash:
                        fasta_out.write(line)
                        count += 1
                        matched = 1
                    else:
                        matched = 0
                except:
                    pass
            elif matched == 1:
                fasta_out.write(line)
        if not args.invert:
            print('Wrote %d contigs of %d listed in %s to %s' % (count, num_queries, args.table, args.output))
        else:
            print('Wrote %d contigs not listed in %s to %s' % (count, args.table, args.output))

    elif not args.list and not args.table:
        print('error: please input either list of contigs (-l) or txt file containing list of contigs (-t)')

    elif args.list and args.table:
        print('error: cannot accept both list (-l) and table (-t). Choose one and try again')


if __name__ == '__main__':
    main()




