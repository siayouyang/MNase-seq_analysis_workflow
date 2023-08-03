#author:siayouyang
#input=$1=.sam


samtools view -Sb $1 > "$1".bam;
samtools sort -n -o "$1"name.bam "$1".bam;
samtools fixmate -m "$1"name.bam "$1"fixmate.bam;
samtools sort -o "$1"pos.bam "$1"fixmate.bam;
samtools markdup -r "$1"pos.bam "$1"rd.bam;
samtools index -b "$1"rd.bam
