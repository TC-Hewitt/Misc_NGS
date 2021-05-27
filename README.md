# Misc_NGS
various tools that help with deep sequencing analysis. Use option -h/--help for usage help.

## Tool descriptions
**blast_filterV2.pyc**

filter and sort a blast output (outfmt 6 or 7) based on the cutoffs you want
for each parameter

**contig_counts.py**

get contig counts and stats from a multi-fasta file

**NLRparser_txt2bed.py**

convert NLR-parser (https://github.com/steuernb/NLR-Parser) txt output to bed

**contig_namer.py**

rename and enumerate sequences in a multi-fasta file

**get_contigs.py**

retrieve multiple specific contigs from a multifasta file usinq a list of seqIDs or file containing seqIDs (e.g. BLAST output)

**BWA_SAM_edit_dist_filter.py**

maximum edit distance option (-n) for BWA (last tried with v0.7.17) does not seem to work. Use this tool instead to filter reads by edit distance from SAM output of BWA (may not work for other aligners due to differences in SAM format). BWA output can be piped directly through for samtools processing. For example: 
```
bwa sampe reference.fasta aln1.sai aln2.sai reads1.fq.gz reads2.fq.gz | python BWA_SAM_edit_dist_filter.py -n 2 | samtools view -hub -o output.bam - 
```

**NoiseOutStats.py**

get statistics from output files of Noisefinder.pyc (https://github.com/TC-Hewitt/MuTrigo)
