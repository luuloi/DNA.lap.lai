qrsh
module load phuluu/python/2.7.8
module load gi/samtools/1.2
module load gi/bedtools/2.22.0


# download the Tranposable element annatation from ucsc
cd /share/ClusterScratch/phuluu/repeat/data/ucsc
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
gunzip rmsk.txt.gz
# convert bed file to gff3
awk 'BEGIN{OFS="\t"}{print $1, ".", ".", $2, $3, ".", "+" , "." , "Name="$4}' rmsk.txt > rmsk.gff3
awk 'BEGIN{OFS="\t"}{print $6, ".", ".", $7, $8, ".", "+" , "." , "Name="$11}' rmsk.txt > rmsk.gff3

jitterbug=/home/phuluu/bin/jitterbug/jitterbug.py
bam=/share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr22_sorted.bam
# bam=/share/ClusterScratch/phuluu/captureSeq/15092015/bams/PrEC.WG/PrEC.bam
# bam=/share/ClusterScratch/phuluu/captureSeq/15092015/bams/LNCaP.WG/LNCaP.bam
gff3=/share/ClusterScratch/phuluu/repeat/data/ucsc/rmsk.gff3

$jitterbug $bam $gff3


# filter TE in NNNNNN reference for PrEC
# make NNNN regions
cd /share/ClusterScratch/phuluu/repeat/data/ucsc
# wget "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz"
# gunzip cytoBand.txt.gz
wget "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz"
gunzip gap.txt.gz
# convert to gff3
awk 'BEGIN{OFS="\t"}{print $2"\t"$3"\t"$4"\t"$8}' gap.txt > N_annot.gff3

# jitterbug_filter_results_func.py [-h] [-g GFF] [-c CONFIG] [-o OUTPUT]
cd /share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG
jitterbug_filter=/home/phuluu/bin/jitterbug/tools/jitterbug_filter_results_func.py
N_annot=/share/ClusterScratch/phuluu/repeat/data/ucsc/N_annot.gff3
# filter out the low coverage and expect distance
$jitterbug_filter -g jitterbug.TE_insertions_paired_clusters.gff3 -c jitterbug.filter_config.txt -o PrEC.WG.TE_insertions_paired_clusters.filtered.gff3
# filter out the NNNNNN regions
intersectBed -a PrEC.WG.TE_insertions_paired_clusters.filtered.gff3 -b $N_annot -v > PrEC.WG.TE_insertions_paired_clusters.filtered.noNs.gff3


# filter TE in NNNNNN reference for LNCaP
# jitterbug_filter_results_func.py [-h] [-g GFF] [-c CONFIG] [-o OUTPUT]
cd /share/ClusterScratch/phuluu/repeat/data/jitterburg/LNCaP.WG
jitterbug_filter=/home/phuluu/bin/jitterbug/tools/jitterbug_filter_results_func.py
N_annot=/share/ClusterScratch/phuluu/repeat/data/ucsc/N_annot.gff3
# filter out the low coverage and expect distance
$jitterbug_filter -g jitterbug.TE_insertions_paired_clusters.gff3 -c jitterbug.filter_config.txt -o LNCaP.WG.TE_insertions_paired_clusters.filtered.gff3
# filter out the NNNNNN regions
intersectBed -a LNCaP.WG.TE_insertions_paired_clusters.filtered.gff3 -b $N_annot -v > LNCaP.WG.TE_insertions_paired_clusters.filtered.noNs.gff3


# USE CASE 2 : looking for somatic insertions in a tumor/normal pair
# STEP 2.1: run Jitterbug on ND and TD samples separately
# Run Jitterbug and filter results as described above. 
# STEP 2.2: Identify insertions present in TD which are absent in ND
# The goal is to identify somatic insertions: that are present in the tumor sample (ND) and absent from the normal sample (TD). This is done in two steps:
# Step 1: take the high-confidence predictions (ie, filtered) from TD and compare them to the whole set (ie, unfiltered) of ND TEI.
# Step 2: for those that remain, inspect the ND bam to be sure that the absence of that TEI is not due to low coverage, 
#         or to a FP. more than 1% discordant reads in a 200bp window around that locus is taken to be a FP.

normal_bam=/share/ClusterScratch/phuluu/captureSeq/15092015/bams/PrEC.WG/PrEC.bam
tumor_bam=/share/ClusterScratch/phuluu/captureSeq/15092015/bams/LNCaP.WG/LNCaP.bam
normal_gff=/share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG/jitterbug.TE_insertions_paired_clusters.gff3
normal_read_stats=/share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG/jitterbug.read_stats.txt
tumor_gff=/share/ClusterScratch/phuluu/repeat/data/jitterburg/LNCaP.WG/LNCaP.WG.TE_insertions_paired_clusters.filtered.noNs.gff3
jitterbug_diff=/home/phuluu/bin/jitterbug/tools/process_ND_TD.py
cd /share/ClusterScratch/phuluu/repeat/data/jitterburg/DTEI

$jitterbug_diff --TD_F $tumor_gff --ND $normal_gff -b $normal_bam -s $normal_read_stats -o LNCaP.PrEC.somatic.TEI.gff3



normal_bam=/share/ClusterScratch/phuluu/captureSeq/15092015/bams/PrEC.WG/PrEC.bam
tumor_bam=/share/ClusterScratch/phuluu/captureSeq/15092015/bams/LNCaP.WG/LNCaP.bam
normal_gff=/share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG/PrEC.WG.TE_insertions_paired_clusters.filtered.noNs.gff3
normal_read_stats=/share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG/jitterbug.read_stats.txt
tumor_gff=/share/ClusterScratch/phuluu/repeat/data/jitterburg/LNCaP.WG/LNCaP.WG.TE_insertions_paired_clusters.filtered.noNs.gff3
jitterbug_diff=/home/phuluu/bin/jitterbug/tools/process_ND_TD.py
cd /share/ClusterScratch/phuluu/repeat/data/jitterburg/DTEI

$jitterbug_diff --TD_F $tumor_gff --ND $normal_gff -b $normal_bam -s $normal_read_stats -o LNCaP.PrEC.noNs.somatic.TEI.gff3



# download WGS bam file from 1000 genomes
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/
cd /share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam




# Why only chr1 has TEI
# Do with chr5 
# extract chr5
/share/ClusterScratch/phuluu/repeat/data/WGBS/bams/LNCaP$ samtools view -hb LNCaP.bam chr5 \
> /share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr5.bam
# sort chr5
samtools sort -T . -O bam -o /share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr5_sorted.bam \
/share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr5.bam
# index chr5
samtools index /share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr5_sorted.bam
# find TEI
jitterbug=/home/phuluu/bin/jitterbug/jitterbug.py
bam=/share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr5_sorted.bam
# check cut -f1 $gff3| sort| uniq -c
gff3=/share/ClusterScratch/phuluu/repeat/data/ucsc/rmsk.gff3
out=/share/ClusterScratch/phuluu/repeat/data/jitterburg/LNCaP.WG.chr5
mkdir -p $out
cd $out
$jitterbug $bam $gff3
jitterbug_filter=/home/phuluu/bin/jitterbug/tools/jitterbug_filter_results_func.py
N_annot=/share/ClusterScratch/phuluu/repeat/data/ucsc/N_annot.gff3
# filter out the low coverage and expect distance
$jitterbug_filter -g jitterbug.TE_insertions_paired_clusters.gff3 -c jitterbug.filter_config.txt -o LNCaP.WG.chr5.TE_insertions_paired_clusters.filtered.gff3
# filter out the NNNNNN regions
intersectBed -a LNCaP.WG.chr5.TE_insertions_paired_clusters.filtered.gff3 -b $N_annot -v > LNCaP.WG.chr5.TE_insertions_paired_clusters.filtered.noNs.gff3

# merge bams
bam5=/share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr5_sorted.bam
bam22=/share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr22_sorted.bam
bam522=/share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr522_sorted.bam

samtools cat -h $bam22 -o $bam522 $bam5 $bam22
samtools index $bam522
# find TEI
jitterbug=/home/phuluu/bin/jitterbug/jitterbug.py
bam=/share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr522_sorted.bam
# check cut -f1 $gff3| sort| uniq -c
gff3=/share/ClusterScratch/phuluu/repeat/data/ucsc/rmsk.gff3
out=/share/ClusterScratch/phuluu/repeat/data/jitterburg/LNCaP.WG.chr522
mkdir -p $out
cd $out
$jitterbug $bam $gff3



# Make a bam file only not contained
# the -F 3854 flag excludes all the reads with the folowwing bits set: 
# read mapped in proper pair
# read unmapped
# mate unmapped
# not primary alignment
# read fails platform/vendor quality checks
# read is PCR or optical duplicate
# supplementary alignment
# IMPORTANTLY: chimeric and non-primary reads are excluded

# the -F 1550 flag excludes all the reads with the folowwing bits set: 
# read mapped in proper pair
# read unmapped
# mate unmapped
# read fails platform/vendor quality checks
# read is PCR or optical duplicate
# IMPORTANTLY: non-primary and supplementary reads are kept, ie, split reads generated by bwa-mem
bam=/share/ClusterScratch/phuluu/repeat/data/WGBS.chr/bams/LNCaP.WG.chr22_sorted.bam
tumor_bam=/share/ClusterScratch/phuluu/captureSeq/15092015/bams/LNCaP.WG/LNCaP.bam
samtools view -hb -F 3854 $bam > ${bam/.bam/.F3854.bam}

samtools view ${bam/.bam/.F3854.bam}| grep HISEQ:55:C3EWTACXX:2:1110:15150:25679
# HISEQ:55:C3EWTACXX:2:1110:15150:25679	81	chr22	16057494	0	101M	chr14	19785412	0	GAAAAGAGCCACATGAATACAGAAGAGGGGAACCTTTTACACACCAGACCTCATACAGGAAAGGGGGTCGTAGATCAGGGAGTTTGGAGTGGAAGGGACAT	##BBB<70''70FBBBB<0'<<'<<BBBFBB707FFFB<0700''<B7''7'FB7'FFBB00''FFBB0IIFFFB<FIIIIIIFFIIIFFFFFFFFFB<BB	\
# XA:Z:fchr14,+19785412,101M,4;fchr18,+14986851,101M,5;fchr2,+132574793,101M,6;	YC:Z:CT	MD:Z:1C20C39C0C37	YD:Z:r	NM:i:4	AS:i:90	XS:i:90	RG:Z:TKCC20140123_LNCaP_P73	PG:Z:MarkDuplicates
samtools view ${bam}| grep HISEQ:55:C3EWTACXX:2:1110:15150:25679
# HISEQ:55:C3EWTACXX:2:1110:15150:25679	81	chr22	16057494	0	101M	chr14	19785412	0	GAAAAGAGCCACATGAATACAGAAGAGGGGAACCTTTTACACACCAGACCTCATACAGGAAAGGGGGTCGTAGATCAGGGAGTTTGGAGTGGAAGGGACAT	##BBB<70''70FBBBB<0'<<'<<BBBFBB707FFFB<0700''<B7''7'FB7'FFBB00''FFBB0IIFFFB<FIIIIIIFFIIIFFFFFFFFFB<BB	\
# XA:Z:fchr14,+19785412,101M,4;fchr18,+14986851,101M,5;fchr2,+132574793,101M,6;	YC:Z:CT	MD:Z:1C20C39C0C37	YD:Z:r	NM:i:4	AS:i:90	XS:i:90	RG:Z:TKCC20140123_LNCaP_P73	PG:Z:MarkDuplicates



# DTIE and DNA methylation
path=/share/ClusterScratch/phuluu/repeat/data
rna=/share/ClusterScratch/phuluu/repeat/data/RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.tsv
jitterburg=jitterburg/DTEI
fn=LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.gff

cut -f1,4,5 $path/$jitterburg/$fn > $path/$jitterburg/${fn/gff/bed}
readfile=/share/ClusterScratch/phuluu/repeat/data/jitterburg/LNCaP.WG/jitterbug.TE_insertions_paired_clusters.supporting_clusters.table
for i in `cut -f 9 LNCaP.PrEC.somatic.TEI.gff3.Tf.not_ND.context_read_stats.gff| cut -d';' -f4| cut -d'=' -f2`;
do
echo $i
awk '{if($1=="R" && $2=='$i' && $8=="mate") print $18}' $readfile >> /share/ClusterScratch/phuluu/repeat/data/jitterburg/LNCaP.WG/seqs.fasta
done


awk -F 'CG' '{sum = sum + NF -1 }END{print sum}' ../LNCaP.WG/seqs.fasta

# PrEC
readfile=/share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG/jitterbug.TE_insertions_paired_clusters.supporting_clusters.table
for i in `cut -f 9 PrEC.WG.TE_insertions_paired_clusters.filtered.noNs.gff3| cut -d';' -f3| cut -d'=' -f2`;
do
echo $i
awk '{if($1=="R" && $2=='$i' && $8=="mate") print $18}' $readfile >> /share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG/seqs.fasta
done

awk -F 'CG' '{sum = sum + NF -1 }END{print sum}' seqs.fasta

# TE type
cd /share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG
cut -f 9 PrEC.WG.TE_insertions_paired_clusters.filtered.noNs.gff3| cut -d';' -f5| cut -d'=' -f2| xargs echo| tr ',' '\n'| sort| uniq -c| sort -k1,1nr

cd /share/ClusterScratch/phuluu/repeat/data/jitterburg/DTEI 
cut -f 9 LNCaP.PrEC.somatic.TEI.gff3.Tf.not_ND.context_read_stats.gff| cut -d';' -f7| cut -d'=' -f2| xargs echo| tr ',' '\n'| sort| uniq -c| sort -k1,1nr

# TE type superfamily
cd /share/ClusterScratch/phuluu/repeat/data/jitterburg/DTEI
cut -f 9 LNCaP.PrEC.somatic.TEI.gff3.Tf.not_ND.context_read_stats.gff| cut -d';' -f12| cut -d'=' -f2| xargs echo| tr ' ' '\n'| sort| uniq -c| sort -k1,1nr

cd /share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG
cut -f 9 PrEC.WG.TE_insertions_paired_clusters.filtered.noNs.gff3| cut -d';' -f12| cut -d'=' -f2| xargs echo| tr ' ' '\n'| sort| uniq -c| sort -k1,1nr

# TE zygosity
/share/ClusterScratch/phuluu/repeat/data/jitterburg/DTEI$ cut -f 9 LNCaP.PrEC.noNs.somatic.TEI.gff3_Tf.notND.gff| cut -d';' -f12| cut -d'=' -f2\
| xargs echo| tr ' ' '\n'| sort| uniq -c| sort -k1,1nr| Rscript -e 'library(ggplot2);tb=read.table(file("stdin"));ggplot(tb, aes(x=V2, y=V1)) + geom_line() + xlab("Zygosity Score") + ylab("Frequency"); ggsave("zygosity.lncp.png")'

/share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG$ cut -f 9 PrEC.WG.TE_insertions_paired_clusters.filtered.noNs.gff3| cut -d';' -f12| cut -d'=' -f2\
| xargs echo| tr ' ' '\n'| sort| uniq -c| sort -k1,1nr| Rscript -e 'library(ggplot2);tb=read.table(file("stdin"));ggplot(tb, aes(x=V2, y=V1)) + geom_line() + xlab("Zygosity Score") + ylab("Frequency"); ggsave("zygosity.prec.png")'

# How broad the TEI
cat $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.bed| awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
cd /share/ClusterScratch/phuluu/repeat/data/jitterburg/PrEC.WG

# Expression
path=/share/ClusterScratch/phuluu/repeat/data
rna=/share/ClusterScratch/phuluu/repeat/data/RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.tsv
jitterburg=jitterburg/DTEI
fn=LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.gff

cut -f1,4,5 $path/$jitterburg/$fn > $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.bed
awk 'NR>1{print $3"\t"$4"\t"$5"\t"$2"\t"$9}' ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.tsv > ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.bed
grep ^chr ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.bed| sort -k1,1 -k2,2n > ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE_sort.bed
bedtools intersect -a $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.bed \
-b ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE_sort.bed \
-wa -wb > $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.DE.RNAseq.intersect


tb <- matrix(c(4, 4, 62, 21827), 2, 2)
chisq.test(tb)

# meth at repeat
cd /share/ClusterScratch/phuluu/repeat/data/WGBS
awk 'NR>1{if($4>0 && $6>0 && $8>0){print $1"\t"$2"\t"$2+1"\t"$3/$4"\t"$5/$6"\t"$7/$8}}' bigTable.LNCaP/bigTable.tsv > bigTable.LNCaP/bigTable_ratio.tsv

bedtools intersect -a /share/ClusterScratch/phuluu/repeat/data/ucsc/rmsk.CpG.bed \
-b /share/ClusterScratch/phuluu/repeat/data/WGBS/bigTable.LNCaP/bigTable_ratio.tsv
-wa -wb > /share/ClusterScratch/phuluu/repeat/data/WGBS/bigTable.LNCaP/rmsk.CpG.meth.bed

# plot
library(data.table)
library(ggplot2)

# for repeat CpGs
tb <- fread("/share/ClusterScratch/phuluu/repeat/data/WGBS/bigTable.LNCaP/rmsk.CpG.meth.bed")
colnames(tb) <- c("chrom", "S", "E", "chrom1", "S1", "E1", "function", "chrom3", "S3", "E3", "LNCaP", "LNCaP_EPI", "PrEC")
# plot
ggplot(tb) + geom_density(aes(x=LNCaP, colour="LNCaP"), show.legend=TRUE) + geom_density(aes(x=LNCaP_EPI, colour="LNCaP_EPI"), show.legend=TRUE) + 
geom_density(aes(x=PrEC, colour="PrEC"), show.legend=TRUE) + xlab("DNA methylation") + ylab("Frequency") + ggtitle("Repeat Elements") + theme(legend.title=element_blank())     
ggsave("rmsk.CpG.meth.png")

ggplot(tb) + geom_density(aes(x=LNCaP, colour="LNCaP"), show.legend=TRUE) + geom_density(aes(x=PrEC, colour="PrEC"), show.legend=TRUE) + 
xlab("DNA methylation") + ylab("Frequency") + ggtitle("Repeat Elements") + theme(legend.title=element_blank())     
ggsave("rmsk.CpG.meth.png")

ggplot(tb) + geom_boxplot(aes(x=LNCaP, colour="LNCaP"), show.legend=TRUE) + geom_boxplot(aes(x=LNCaP, colour="PrEC"), show.legend=TRUE) +
xlab("DNA methylation") + ylab("Frequency") + ggtitle("Repeat Elements") + theme(legend.title=element_blank())     
ggsave("rmsk.CpG.meth.boxplot.png")

png("rmsk.CpG.meth.boxplot.png")
boxplot(tb$LNCaP, tb$PrEC, notch=TRUE, names=c("LNCaP", "PrEC"),
  col=(c("gold","darkgreen")),
  main="Repeat Elements", xlab="sample")
dev.off()


# for whole genome
tb <- fread("/share/ClusterScratch/phuluu/repeat/data/WGBS/bigTable.LNCaP/bigTable_ratio.tsv")
colnames(tb) <- c("chrom", "S", "E", "LNCaP", "LNCaP_EPI", "PrEC")
# plot
ggplot(tb) + geom_density(aes(x=LNCaP, colour="LNCaP"), show.legend=TRUE) + geom_density(aes(x=PrEC, colour="PrEC"), show.legend=TRUE) + 
xlab("DNA methylation") + ylab("Frequency") + ggtitle("Whole genome") + theme(legend.title=element_blank())     
ggsave("whole.CpG.meth.png")

ggplot(tb) + geom_boxplot(aes(x=LNCaP, colour="LNCaP"), show.legend=TRUE) + geom_boxplot(aes(x=LNCaP, colour="PrEC"), show.legend=TRUE) +
xlab("DNA methylation") + ylab("Frequency") + ggtitle("Repeat Elements") + theme(legend.title=element_blank())     
ggsave("whole.CpG.meth.boxplot.png")

png("whole.CpG.meth.boxplot.png")
boxplot(tb$LNCaP, tb$PrEC, notch=TRUE, names=c("LNCaP", "PrEC"),
  col=(c("gold","darkgreen")),
  main="Whole genome", xlab="sample")
dev.off()


# Expression 2
module load gi/bedtools/2.22.0
hg19=/home/phuluu/methods/darlo/annotations/hg19/chrom.sizes.short
path=/share/ScratchGeneral/phuluu/repeat/data
rna=/share/ScratchGeneral/phuluu/repeat/data/RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.tsv
jitterburg=jitterburg/DTEI
fn=LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.gff

cut -f1,4,5 $path/$jitterburg/$fn > $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.bed
awk 'NR>1{print $3"\t"$4"\t"$5"\t"$2"\t"$9}' ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.tsv > ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.bed
grep ^chr ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.bed| sort -k1,1 -k2,2n > ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE_sort.bed
bedtools intersect -a $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.bed \
-b ../../RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE_sort.bed \
-wa -wb > $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.DE.RNAseq.intersect


tb <- matrix(c(4, 4, 62, 21827), 2, 2)
chisq.test(tb)


# change the DE expression file
module load phuluu/R/3.1.2
R
library(GenomicRanges)
path = "/share/ScratchGeneral/phuluu/repeat/data/RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.tsv"
df <- read.table(path, sep="\t", header=TRUE)
df1 <- df[complete.cases(df),]
df1$chr <- factor(df1$chr)
df1$strand <- factor(df1$strand)
de.gr <- makeGRangesFromDataFrame(df1, keep.extra.columns=TRUE)

path = "/share/ScratchGeneral/phuluu/repeat/data/jitterburg/DTEI/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.slop.500.bed"
df <- read.table(path, sep="\t")
colnames(df) <- c("chr","start","end")
DTEI.gr <- makeGRangesFromDataFrame(df)

# sort
seqlevels(de.gr) <- sort(seqlevels(de.gr))
de.gr <- sort(de.gr)
seqlevels(DTEI.gr) <- sort(seqlevels(DTEI.gr))
DTEI.gr <- sort(DTEI.gr)

ov <- de.gr[de.gr %over% DTEI.gr]
table(ov$DGE.status)
DOWN   NC   UP 
   3   29    2 

table(de.gr$DGE.status)
DOWN    NC    UP 
1415 16499   811 


# make slop for TEI
size=500
bedtools slop -i $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.bed \
-g $hg19 -b $size > $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.slop.$size.bed
# intersect slop TEI with expression
bedtools intersect -a $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.slop.$size.bed \
-b $path/RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE_sort.bed \
-wa -wb > $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.slop.$size.DE.RNAseq.intersect

cut -f8 $path/$jitterburg/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.slop.$size.DE.RNAseq.intersect| sort| uniq -c

# nearest gene with expression
module load phuluu/R/3.1.2
R
library(GenomicRanges)
library(data.table)
library(ggplot2)

path = "/share/ScratchGeneral/phuluu/repeat/data/RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.tsv"
pathout = "/share/ScratchGeneral/phuluu/repeat/data/jitterburg/DTEI/"

df <- read.table(path, sep="\t", header=TRUE)
df1 <- df[complete.cases(df),]
df1$chr <- factor(df1$chr)
df1$strand <- factor(df1$strand)
de.gr <- makeGRangesFromDataFrame(df1, keep.extra.columns=TRUE)

path = "/share/ScratchGeneral/phuluu/repeat/data/jitterburg/DTEI/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.bed"
df <- read.table(path, sep="\t")
colnames(df) <- c("chr","start","end")
DTEI.gr <- makeGRangesFromDataFrame(df)

# sort
seqlevels(de.gr) <- sort(seqlevels(de.gr))
de.gr <- sort(de.gr)
seqlevels(DTEI.gr) <- sort(seqlevels(DTEI.gr))
DTEI.gr <- sort(DTEI.gr)
nearest.gene <- de.gr[nearest(DTEI.gr, de.gr),]
table(nearest.gene$DGE.status)
DOWN   NC   UP 
   4   53    5
table(de.gr$DGE.status)
DOWN    NC    UP 
1415 16499   811 


tb <- matrix(c(4+5, 53, 1415+811, 16499), 2, 2)
chisq.test(tb)


# DNA methylation
DOWN <- DTEI.gr[nearest.gene$DGE.status=="DOWN"]
NC <- DTEI.gr[nearest.gene$DGE.status=="NC"]
UP <- DTEI.gr[nearest.gene$DGE.status=="UP"]
btp <- "/share/ScratchGeneral/phuluu/repeat/data/WGBS/bigTable.LNCaP/bigTable.tsv"
bt <- fread(btp, header=TRUE)
# LNCaP.gr <- GRanges(bt$chr, IRanges(bt$position, (bt$position + 1)), strand="+", score=bt$LNCaP.C/bt$LNCaP.cov)
LNCaP.gr <- GRanges(bt$chr, IRanges(bt$position, (bt$position + 1)), strand="+", score=bt$LNCaP.C/bt$LNCaP.cov, C=bt$LNCaP.C, COV=bt$LNCaP.cov)
PrEC.gr <- GRanges(bt$chr, IRanges(bt$position, (bt$position + 1)), strand="+", score=bt$PrEC.C/bt$PrEC.cov, C=bt$PrEC.C, COV=bt$PrEC.cov)


# get pvalue out of lm ouput
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

# DOWN
DOWNdiff <- lapply(DOWN, function(x) list(LNCaP.gr[LNCaP.gr%over%x],PrEC.gr[PrEC.gr%over%x]))

DOWNp <- lapply(1:length(DOWNdiff), function(x){
	print(x)
	if(length(DOWNdiff[[x]][[1]])!=0){
		temp <- data.frame(pos=paste0(seqnames(DOWNdiff[[x]][[1]]), "_", start(DOWNdiff[[x]][[1]])), LNCaP=DOWNdiff[[x]][[1]]$score, PrEC=DOWNdiff[[x]][[2]]$score)
		temp <- melt(temp, "pos", value.name="DNA.methylation", variable.name="Sample")
		ggplot(aes(x = as.numeric(temp$pos), y = DNA.methylation, colour=Sample), data = temp) + geom_point() + geom_smooth(se=FALSE) +
		xlab("Genome.Coordinate") + ylim(0, 1.0) + scale_x_continuous(breaks=as.numeric(temp$pos), labels=as.character(temp$pos)) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1))
		ggsave(paste0(pathout, "DOWN_", x, ".png"))
		return(lmp(lm(DNA.methylation~Sample, data=temp)))
	}
})
a = 0
for (i in DOWNp){ if(!is.null(i)){if (!is.na(i)){ print(i); if(i<0.01){a=a+1}}}}
a = 2
total = 4

# UP
UPdiff <- lapply(UP, function(x) list(LNCaP.gr[LNCaP.gr%over%x],PrEC.gr[PrEC.gr%over%x]))

UPp <- lapply(1:length(UPdiff), function(x){
	print(x)
	if(length(UPdiff[[x]][[1]])!=0){
		temp <- data.frame(pos=paste0(seqnames(UPdiff[[x]][[1]]), "_", start(UPdiff[[x]][[1]])), LNCaP=UPdiff[[x]][[1]]$score, PrEC=UPdiff[[x]][[2]]$score)
		temp <- melt(temp, "pos", value.name="DNA.methylation", variable.name="Sample")
		ggplot(aes(x = as.numeric(temp$pos), y = DNA.methylation, colour=Sample), data = temp) + geom_point() + geom_smooth(se=FALSE) +
		xlab("Genome.Coordinate") + ylim(0, 1.0) + scale_x_continuous(breaks=as.numeric(temp$pos), labels=as.character(temp$pos)) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1))
		ggsave(paste0(pathout, "UP_", x, ".png"))
		return(lmp(lm(DNA.methylation~Sample, data=temp)))
	}
})
a = 0
for (i in UPp){ if(!is.null(i)){if (!is.na(i)){ print(i); if(i<0.01){a=a+1}}}}
a = 3
total = 5

# NC
NCdiff <- lapply(NC, function(x) list(LNCaP.gr[LNCaP.gr%over%x],PrEC.gr[PrEC.gr%over%x]))

NCp <- lapply(1:length(NCdiff), function(x){
	print(x)
	if(length(NCdiff[[x]][[1]])!=0){
		temp <- data.frame(pos=paste0(seqnames(NCdiff[[x]][[1]]), "_", start(NCdiff[[x]][[1]])), LNCaP=NCdiff[[x]][[1]]$score, PrEC=NCdiff[[x]][[2]]$score)
		temp <- melt(temp, "pos", value.name="DNA.methylation", variable.name="Sample")
		ggplot(aes(x = as.numeric(temp$pos), y = DNA.methylation, colour=Sample), data = temp) + geom_point() + geom_smooth(se=FALSE) +
		xlab("Genome.Coordinate") + ylim(0, 1.0) + scale_x_continuous(breaks=as.numeric(temp$pos), labels=as.character(temp$pos)) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1))
		ggsave(paste0(pathout, "NC_", x, ".png"))
		return(lmp(lm(DNA.methylation~Sample, data=temp)))
	}
})

a = 0
for (i in NCp){ if(!is.null(i)){if (!is.na(i)){ print(i); if(i<0.01){a=a+1}}}}
a = 17
total = 54




# DNA methyloation of LNCaP specific DTEI
module load phuluu/R/3.1.2
R
library(GenomicRanges)
library(data.table)
library(ggplot2)

path = "/share/ScratchGeneral/phuluu/repeat/data/RNAseq/DE/ALL_GENE_LNCaP_PrEC_DGE.tsv"
pathout = "/share/ScratchGeneral/phuluu/repeat/data/jitterburg/DTEI/lethal/"

df <- read.table(path, sep="\t", header=TRUE)
df1 <- df[complete.cases(df),]
df1$chr <- factor(df1$chr)
df1$strand <- factor(df1$strand)
de.gr <- makeGRangesFromDataFrame(df1, keep.extra.columns=TRUE)

path = "/share/ScratchGeneral/phuluu/repeat/data/jitterburg/DTEI/LNCaP.PrEC.noNs.somatic.TEI.gff3.Tf.not_ND.context_read_stats.bed"
df <- read.table(path, sep="\t")
colnames(df) <- c("chr","start","end")
DTEI.gr <- makeGRangesFromDataFrame(df)

# sort
seqlevels(de.gr) <- sort(seqlevels(de.gr))
de.gr <- sort(de.gr)
seqlevels(DTEI.gr) <- sort(seqlevels(DTEI.gr))
DTEI.gr <- sort(DTEI.gr)
nearest.gene <- de.gr[nearest(DTEI.gr, de.gr),]
table(nearest.gene$DGE.status)
DOWN   NC   UP 
   4   53    5
table(de.gr$DGE.status)
DOWN    NC    UP 
1415 16499   811 

tb <- matrix(c(4+5, 53, 1415+811, 16499), 2, 2)
chisq.test(tb)

# DNA methylation
DOWN <- DTEI.gr[nearest.gene$DGE.status=="DOWN"]
NC <- DTEI.gr[nearest.gene$DGE.status=="NC"]
UP <- DTEI.gr[nearest.gene$DGE.status=="UP"]
btp <- "/share/ScratchGeneral/phuluu/repeat/data/WGBS/bigTable.lethal/bigTable.tsv"
bt <- fread(btp, header=TRUE)
# LNCaP.gr <- GRanges(bt$chr, IRanges(bt$position, (bt$position + 1)), strand="+", score=bt$LNCaP.C/bt$LNCaP.cov)
S204C.gr <- GRanges(bt$chr, IRanges(bt$position, (bt$position + 1)), strand="+", score=bt$"204C.C"/bt$"204C.cov", C=bt$"204C.C", COV=bt$"204C.cov")
S1601C.gr <- GRanges(bt$chr, IRanges(bt$position, (bt$position + 1)), strand="+", score=bt$"1601C.C"/bt$"1601C.cov", C=bt$"1601C.C", COV=bt$"1601C.cov")
S1601N.gr <- GRanges(bt$chr, IRanges(bt$position, (bt$position + 1)), strand="+", score=bt$"1601N.C"/bt$"1601N.cov", C=bt$"1601N.C", COV=bt$"1601N.cov")
btp <- "/share/ScratchGeneral/phuluu/repeat/data/WGBS/bigTable.LNCaP/bigTable.tsv"
bt <- fread(btp, header=TRUE)
LNCaP.gr <- GRanges(bt$chr, IRanges(bt$position, (bt$position + 1)), strand="+", score=bt$LNCaP.C/bt$LNCaP.cov, C=bt$LNCaP.C, COV=bt$LNCaP.cov)
PrEC.gr <- GRanges(bt$chr, IRanges(bt$position, (bt$position + 1)), strand="+", score=bt$PrEC.C/bt$PrEC.cov, C=bt$PrEC.C, COV=bt$PrEC.cov)

# get pvalue out of lm ouput
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

# DOWN
DOWNdiff <- lapply(DOWN, function(x) list(S204C.gr[S204C.gr%over%x],S1601C.gr[S1601C.gr%over%x],S1601N.gr[S1601N.gr%over%x],LNCaP.gr[LNCaP.gr%over%x],PrEC.gr[PrEC.gr%over%x]))

DOWNp <- lapply(1:length(DOWNdiff), function(x){
	print(x)
	if(length(DOWNdiff[[x]][[1]])!=0){
		temp <- data.frame(pos=paste0(seqnames(DOWNdiff[[x]][[1]]), "_", start(DOWNdiff[[x]][[1]])), S204C=DOWNdiff[[x]][[1]]$score, S1601C=DOWNdiff[[x]][[2]]$score, S1601N=DOWNdiff[[x]][[3]]$score, LNCaP=DOWNdiff[[x]][[4]]$score, PrEC=DOWNdiff[[x]][[5]]$score)
		temp <- melt(temp, "pos", value.name="DNA.methylation", variable.name="Sample")
		ggplot(aes(x = as.numeric(temp$pos), y = DNA.methylation, colour=Sample), data = temp) + geom_point() + geom_smooth(se=FALSE) +
		xlab("Genome.Coordinate") + ylim(0, 1.0) + scale_x_continuous(breaks=as.numeric(temp$pos), labels=as.character(temp$pos)) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1))
		ggsave(paste0(pathout, "DOWN_", x, ".png"))
		return(lmp(lm(DNA.methylation~Sample, data=temp)))
	}
})
a = 0
for (i in DOWNp){ if(!is.null(i)){if (!is.na(i)){ print(i); if(i<0.01){a=a+1}}}}


# UP
UPdiff <- lapply(UP, function(x) list(S204C.gr[S204C.gr%over%x],S1601C.gr[S1601C.gr%over%x],S1601N.gr[S1601N.gr%over%x],LNCaP.gr[LNCaP.gr%over%x],PrEC.gr[PrEC.gr%over%x]))
UPp <- lapply(1:length(UPdiff), function(x){
	print(x)
	if(length(UPdiff[[x]][[1]])!=0){
		temp <- data.frame(pos=paste0(seqnames(UPdiff[[x]][[1]]), "_", start(UPdiff[[x]][[1]])), S204C=UPdiff[[x]][[1]]$score, S1601C=UPdiff[[x]][[2]]$score, S1601N=UPdiff[[x]][[3]]$score, LNCaP=UPdiff[[x]][[4]]$score, PrEC=UPdiff[[x]][[5]]$score)
		temp <- melt(temp, "pos", value.name="DNA.methylation", variable.name="Sample")
		ggplot(aes(x = as.numeric(temp$pos), y = DNA.methylation, colour=Sample), data = temp) + geom_point() + geom_smooth(se=FALSE) +
		xlab("Genome.Coordinate") + ylim(0, 1.0) + scale_x_continuous(breaks=as.numeric(temp$pos), labels=as.character(temp$pos)) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1))
		ggsave(paste0(pathout, "UP_", x, ".png"))
		return(lmp(lm(DNA.methylation~Sample, data=temp)))
	}
})
a = 0
for (i in UPp){ if(!is.null(i)){if (!is.na(i)){ print(i); if(i<0.01){a=a+1}}}}


# NC
NCdiff <- lapply(NC, function(x) list(S204C.gr[S204C.gr%over%x],S1601C.gr[S1601C.gr%over%x],S1601N.gr[S1601N.gr%over%x],LNCaP.gr[LNCaP.gr%over%x],PrEC.gr[PrEC.gr%over%x]))
NCp <- lapply(1:length(NCdiff), function(x){
	print(x)
	if(length(NCdiff[[x]][[1]])!=0){
		temp <- data.frame(pos=paste0(seqnames(NCdiff[[x]][[1]]), "_", start(NCdiff[[x]][[1]])), S204C=NCdiff[[x]][[1]]$score, S1601C=NCdiff[[x]][[2]]$score, S1601N=NCdiff[[x]][[3]]$score, LNCaP=NCdiff[[x]][[4]]$score, PrEC=NCdiff[[x]][[5]]$score)
		temp <- melt(temp, "pos", value.name="DNA.methylation", variable.name="Sample")
		ggplot(aes(x = as.numeric(temp$pos), y = DNA.methylation, colour=Sample), data = temp) + geom_point() + geom_smooth(se=FALSE) +
		xlab("Genome.Coordinate") + ylim(0, 1.0) + scale_x_continuous(breaks=as.numeric(temp$pos), labels=as.character(temp$pos)) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1))
		ggsave(paste0(pathout, "NC_", x, ".png"))
		return(lmp(lm(DNA.methylation~Sample, data=temp)))
	}
})
a = 0
for (i in NCp){ if(!is.null(i)){if (!is.na(i)){ print(i); if(i<0.01){a=a+1}}}}
