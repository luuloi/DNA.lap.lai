# all types of file format
https://genome.ucsc.edu/FAQ/FAQformat.html
http://www.broadinstitute.org/annotation/argo/help/gff3.html

# direct convert 
# download the tools from here
# example
bedGraphToBigWig
bedToBigBed
bedToExons
bigWigToBedGraph
bigWigToWig
http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v287/
# bedGraph format
# chrom start    end      score
# chr19 49302000 49302300 -1.0
# chr19 49302300 49302600 -0.75
# chr19 49302600 49302900 -0.50

# samtools compute coverage for each position of a bam file
# modules
module load gi/samtools/1.2
bedGraphToBigWig=/home/phuluu/bin/ucsc/bedGraphToBigWig
fetchChromSizes=/home/phuluu/bin/ucsc/fetchChromSizes
# get hg19 genome sizes of each chromosomes
$fetchChromSizes hg19 > hg19.chrom.sizes
# convert bam file to bigwig
samtools depth ERR127306.asd.bam| head
# chrM	1	1
# chrM	2	12
# chrM	3	17
# chrM	4	110
samtools depth ERR127306.asd.bam| awk '{print $1"\t"$2"\t"$2+1"\t"$3}'| head| sort -k1,1 -k2,2n
chrM	1	2	1
chrM	2	3	12
chrM	3	4	17
chrM	4	5	110
# first make a bedGraph file
samtools depth ERR127306.asd.bam| awk '{print $1"\t"$2"\t"$2"\t"$3}'| sort -k1,1 -k2,2n > ERR127306.asd.bedGraph
# convert bedGraph to bigwig file
$bedGraphToBigWig ERR127306.asd.bedGraph hg19.chrom.sizes ERR127306.asd.bw

# convert bam to bed
http://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html
bedtools bamtobed -i reads.bam | head -3
# chr7   118970079   118970129   TUPAC_0001:3:1:0:1452#0/1   37   -
# chr7   118965072   118965122   TUPAC_0001:3:1:0:1452#0/2   37   +
# chr11  46769934    46769984    TUPAC_0001:3:1:0:1472#0/1   37   -

# convert bed file or tab txt file to gff3
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
gunzip rmsk.txt.gz
# convert bed file to gff3
head rmsk.txt| column -t
585  1504  13   4    13  chr1  10000  10468  -249240153  +  (CCCTAA)n  Simple_repeat  Simple_repeat  1      463   0     1
585  3612  114  270  13  chr1  10468  11447  -249239174  -  TAR1       Satellite      telo           -399   1712  483   2
585  437   235  186  35  chr1  11503  11675  -249238946  -  L1MC       LINE           L1             -2236  5646  5449  3
585  239   294  19   10  chr1  11677  11780  -249238841  -  MER5B      DNA            hAT-Charlie    -74    104   1     4
585  318   230  38   0   chr1  15264  15355  -249235266  -  MIR3       SINE           MIR            -119   143   49    5
585  203   162  0    0   chr1  16712  16749  -249233872  +  (TGG)n     Simple_repeat  Simple_repeat  1      37    0     6
585  239   338  148  0   chr1  18906  19048  -249231573  +  L2a        LINE           L2             2942   3104  -322  7
# covert tab txt to gff3
awk 'BEGIN{OFS="\t"}{print $6, ".", ".", $7, $8, ".", "+" , "." , "Name="$11}' rmsk.txt| head| column -t
# chr1  .  .  10000  10468  .  +  .  Name=(CCCTAA)n
# chr1  .  .  10468  11447  .  +  .  Name=TAR1
# chr1  .  .  11503  11675  .  +  .  Name=L1MC
# chr1  .  .  11677  11780  .  +  .  Name=MER5B
# chr1  .  .  15264  15355  .  +  .  Name=MIR3
# chr1  .  .  16712  16749  .  +  .  Name=(TGG)n
# chr1  .  .  18906  19048  .  +  .  Name=L2a
# real convert
awk 'BEGIN{OFS="\t"}{print $6, ".", ".", $7, $8, ".", "+" , "." , "Name="$11}' rmsk.txt > rmsk.gff3
