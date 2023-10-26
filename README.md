# Misc_NGS
various tools that help with seq processing and analysis. Use option -h/--help for usage help.

## Tool descriptions
### blast_filterV2.pyc

filter and sort a blast output (outfmt 6 or 7) based on the cutoffs you want
for each parameter

### contig_counts.py

get contig counts and stats from a multi-fasta file

### NLRparser_txt2bed.py

convert [NLR-parser](https://github.com/steuernb/NLR-Parser) txt output to bed

### contig_namer.py

rename and enumerate sequences in a multi-fasta file

### get_contigs.py

retrieve multiple specific contigs from a multifasta file usinq a list of seqIDs or file containing seqIDs (e.g. BLAST output)

### BWA_SAM_edit_dist_filter.py

maximum edit distance option (-n) for BWA (last tried with v0.7.17) does not seem to work. Use this tool instead to filter reads by edit distance from SAM output of BWA (may not work for other aligners due to differences in SAM format). BWA output can be piped directly through for samtools processing. For example: 
```
bwa sampe reference.fasta aln1.sai aln2.sai reads1.fq.gz reads2.fq.gz | python BWA_SAM_edit_dist_filter.py -n 2 | samtools view -hub -o output.bam - 
```

### CoverageOverlap2GFF.py

Compare read coverage from BAM/SAM file (requires samtools depth output) against gff annotation. Generates gff output of coverage intervals both non-overlapping and overlapping with gff features (strand info not retained, output attributes: <ID=seqid_0start|overlapping gff ID> <mean cov\> <median cov\> <num gaps\> <overlap type\> <overlap len\> <overlap %>). Useful for quick lookup of which features have coverage or not and to assess how well the annotation captures full exome (for RNAseq mapped reads)

**Example: scraping unannotated regions with RNA expression for potential signal peptides**<br />
Sometimes de novo gene prediction fails to capture regions where there is obvious RNA expression, even when guided with RNAseq data. This is a quick and dirty way to identify a class of potential genes that may have failed to be annotated - in this example, ones with a secreted signal.

**1.** run [SAMtools](https://github.com/samtools/samtools) depth on RNAseq mapping to ref and pipe through CoverageOverlap2GFF.py using mRNA lines from gff as reference to find overlaps
```
samtools depth planta_rna.bam | ./CoverageOverlap2GFF.py -g planta_genemodels.gff3 -o rnacov_vs_modelsgff.overlaps
```
**2.** from output, pull out features with "null" overlap type (no existing mRNA annotation)
```
grep null rnacov_vs_modelsgff.overlaps.gff > rnacov_vs_modelsgff.overlaps.null.gff
```
**3.** optional: pull out features satisfying a certain cutoff* (e.g. median coverage 10x or more)
```
awk -F'[;=]' '$6 >= 10' rnacov_vs_modelsgff.overlaps.null.gff > rnacov_vs_modelsgff.overlaps.null.med10X.gff
```
>*fields for awk: $4=mean cov, $6=median cov, $8=n gaps, $10=overlap type<br />

**4.** using [bedtools](https://github.com/arq5x/bedtools2), extract the DNA fasta sequences for these features from the reference fasta
```
bedtools getfasta -fi planta.fasta -bed rnacov_vs_modelsgff.overlaps.null.med10X.gff -fo rnacov_vs_modelsgff.overlaps.null.med10X.fna
```
**5.** using [OrfM](https://github.com/wwood/OrfM), generate translations in all 6 frames for these dna seqs (header line suffix includes _startPosition_frameNumber_orfNumber)
```
orfm -p rnacov_vs_modelsgff.overlaps.null.med10X.fna > rnacov_vs_modelsgff.overlaps.null.med10X.orfm.faa
```
**6.** run [SignalP](https://github.com/fteufel/signalp-6.0) on the OrfM aa seqs
```
signalp -fasta rnacov_vs_modelsgff.overlaps.null.med10X.orfm.faa
```
**7.** pull out IDs with secretory hits from signalp summary file
```
grep 'SP(Sec/SPI)' rnacov_vs_modelsgff.overlaps.null.med10X.orfm_summary.signalp4 > rnacov_vs_modelsgff.overlaps.null.med10X.orfm.signalp4_hits.txt
```
**8.** retrieve the original matching coverage gff lines
```
for i in $(cut -d"-" -f1 rnacov_vs_modelsgff.overlaps.null.med10X.orfm.signalp4_hits.txt | uniq | sed 's/:/_/g')
do grep "$i;" rnacov_vs_modelsgff.overlaps.null.gff >> rnacov_vs_modelsgff.overlaps.null.signalp4.gff
done
```
Note this may not salvage all uncaptured secreted genes, especially if signal peptides span splice junctions - that is why it is quick and dirty. However, it had ~90% recovery in my tests.
