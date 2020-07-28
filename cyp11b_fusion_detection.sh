#!/usr/bin/bash

## Parameters
fq=$1
sid=$2
refgenome=/data1/reference/GRCh38.p12/GRCh38.primary_assembly.genome.fa

root=/home/chialanghsu/nanopore
db=$root/db
src=$root/scripts

## Only consider the reads with length larger than 3k
awk '
    BEGIN {FS = "\t" ; OFS = "\n"}
    {header = $0 ;
    getline seq ;
    getline qheader ;
    getline qseq ;
    if (length(seq) >= 3000)
    {print header, seq, qheader, qseq}}
' < $fq > ${sid}.filtered.fq

## Map reads onto genome using NGMLR and save as BAM file
$CONDA3/ngmlr -r $refgenome -t 8 -q ${sid}.filtered.fq -x ont | \
    samtools sort -O BAM -o ${sid}.ngmlr.bam

## Index BAM file
samtools index ${sid}.ngmlr.bam

## Retreive reads mapped to CYP11B1 and CYP11B2
samtools view -f 2048 ${sid}.ngmlr.bam "chr8:142872206-142918080" > ${sid}.target.sam
$src/sam2fasta.py ${sid}.target.sam > ${sid}.target.fasta
#| samtools fasta -0 ${sid}.target.fasta -

index=1
while IFS= read header
do
    read sequence
    echo $header > read.${index}.fasta
    echo $sequence >> read.${index}.fasta

    $CONDA2/ssw_test -c $db/CYP11B1.fasta read.${index}.fasta > read.${index}.b1.aln
    $CONDA2/ssw_test -c $db/CYP11B2.fasta read.${index}.fasta > read.${index}.b2.aln

    $src/edit_read.py read.${index}.b1.aln read.${index}.b2.aln > read.${index}.new.fasta

    $CONDA2/ssw_test -c $db/CYP11B1.fasta read.${index}.new.fasta > read.${index}.b1.aln
    $CONDA2/ssw_test -c $db/CYP11B2.fasta read.${index}.new.fasta > read.${index}.b2.aln

    $src/pattern_align_v2.py read.${index}.b1.aln read.${index}.b2.aln > read.${index}.aln

    rm read.${index}.fasta read.${index}.b1.aln read.${index}.b2.aln read.${index}.new.fasta
    index=$(expr $index + 1)
done < ${sid}.target.fasta

$src/mergeAlnRes_v2.py $index $sid
rm read.*.aln
