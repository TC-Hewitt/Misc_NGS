#!/usr/bin/env python


import argparse


def main():

    # Parse arguments.
    parser = argparse.ArgumentParser(description='Rename contigs in FASTA files.')
    parser.add_argument('-i', '--input', help='indicate input FASTA file', required=True)
    parser.add_argument('-o', '--output', help='indicate output FASTA file', required=True)
    parser.add_argument('-pre', '--prefix', help='indicate prefix name. Default is <contig_>', default=None, type=str, required=False)
    parser.add_argument('-pos', '--suffix', help='indicate suffix name. Default is null', default='', type=str, required=False)
    parser.add_argument('-z', '--zfill', help='include leading zeros. Indicate min digit length', default=0, type=int, required=False)
    parser.add_argument('-s', '--startn', help='indicate starting number. Default = 1', default=1, type=int, required=False)
    parser.add_argument('-k', '--keep', help='indicate True to keep original names but add prefix and/or suffix', default=False, type=bool, required=False)
    args = parser.parse_args()

    # Open FASTA.
    fasta_in = open(args.input, 'rU')

    # Create FASTA output file.
    fasta_out = open(args.output, 'wb')

    # Start counter.
    if args.startn:
        count = args.startn
    else:
        count = 1

    # Parse file and write to output.
    print('Parsing %s...' % args.input)

    if args.keep == True:
        if not args.prefix:
            args.prefix = ''
        for line in fasta_in.readlines():
            if line.startswith('>'):
                nuline = line.strip('>')
                contig_id = '>' + args.prefix + nuline.strip('\n') + args.suffix + '\n'
                fasta_out.write(contig_id)
                count += 1
            else:
                fasta_out.write(line)
    else:
        if not args.prefix:
            args.prefix = 'contig_'
        for line in fasta_in.readlines():
            if line.startswith('>'):
                contig_id = '>' + args.prefix + str(count).zfill(args.zfill) + args.suffix + '\n'
                fasta_out.write(contig_id)
                count += 1
            else:
                fasta_out.write(line)

    # Finish.
    fasta_out.close()
    fasta_in.close()
    print('Wrote %d contigs to %s.' % (count-args.startn, args.output))


if __name__ == '__main__':
    main()


