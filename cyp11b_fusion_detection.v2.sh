#!/usr/bin/bash

## Parameters
fq=$1
sid=$2
outdir=$3
refgenome=/data1/reference/GRCh38.p12/GRCh38.primary_assembly.genome.fa

root=/home/chialanghsu/BINFO0039_GRA
db=$root/db
src=$root/scripts

outdir=$outdir/$sid

if [ ! -d "$outdir" ]; then
    mkdir -p $outdir
    mkdir -p $outdir/aln1
    mkdir -p $outdir/aln2
    mkdir -p $outdir/reads
    mkdir -p $outdir/temp
fi

## Only consider the reads with length larger than 3k
awk '
    BEGIN {FS = "\t" ; OFS = "\n"}
    {header = $0 ;
    getline seq ;
    getline qheader ;
    getline qseq ;
    if (length(seq) >= 2000)
    {print header, seq, qheader, qseq}}
' < $fq > $outdir/${sid}.filtered.fq

## Map reads onto genome using NGMLR and save as BAM file
ngmlr -r $refgenome -t 8 -q $outdir/${sid}.filtered.fq -x ont | \
    samtools sort -O BAM -o $outdir/${sid}.ngmlr.bam
samtools index $outdir/${sid}.ngmlr.bam

## Map reads onto genome using minimap2
#minimap2 -ax map-ont $refgenome $outdir/${sid}.filtered.fq | \
#    samtools sort -O BAM -o $outdir/${sid}.minimap2.bam
#samtools index $outdir/${sid}.minimap2.bam

## Retreive reads mapped to CYP11B1 and CYP11B2
samtools view -f 2048 $outdir/${sid}.ngmlr.bam "chr8:142872206-142918080" > $outdir/${sid}.target.sam
$src/sam2fasta.py $outdir/${sid}.target.sam > $outdir/${sid}.target.fasta
#| samtools fasta -0 ${sid}.target.fasta -

tempDir=$outdir/temp
aln1Dir=$outdir/aln1
aln2Dir=$outdir/aln2
readDir=$outdir/reads

index=1
while IFS= read header
do
    read sequence
    echo $header > $tempDir/read.${index}.fasta
    echo $sequence >> $tempDir/read.${index}.fasta

    ssw_test -c $db/CYP11B1.fasta $tempDir/read.${index}.fasta > $aln1Dir/read.${index}.b1.aln
    ssw_test -c $db/CYP11B2.fasta $tempDir/read.${index}.fasta > $aln1Dir/read.${index}.b2.aln

    $src/edit_read.py $aln1Dir/read.${index}.b1.aln $aln1Dir/read.${index}.b2.aln > $readDir/read.${index}.new.fasta

    ssw_test -c $db/CYP11B1.fasta $readDir/read.${index}.new.fasta > $aln2Dir/read.${index}.b1.aln
    ssw_test -c $db/CYP11B2.fasta $readDir/read.${index}.new.fasta > $aln2Dir/read.${index}.b2.aln

    index=$(expr $index + 1)
done < $outdir/${sid}.target.fasta

$src/read_multialigner.py $aln2Dir > $outdir/${sid}.summary.txt
Rscript $src/typing_plot.R -f $outdir/${sid}.summary.txt -o $outdir -s $sid

#$src/mergeAlnRes_v2.py $index $sid
#rm read.*.aln
