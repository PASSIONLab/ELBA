#!/bin/bash

INPUT_FOLDER=$1
READS=$2
GENOME=$3
OUTPUT=$4

QUAST=$TEST/contigs/quast/quast-5.0.2/quast.py

rm -rf $OUTPUT
mkdir $OUTPUT

cat $INPUT_FOLDER/*contigs_rank_*.fa > $OUTPUT/elba.contigs.fa

#python $QUAST $OUTPUT/elba.contigs.fa -o $OUTPUT/quast_results -r $GENOME --min-alignment 20000 --extensive-mis-size 5000000 --min-identity 90 --threads 32
python $QUAST $OUTPUT/elba.contigs.fa -o $OUTPUT/quast_results -r $GENOME --threads 32

minimap2 -x asm5 -t64 -a $GENOME $OUTPUT/elba.contigs.fa | samtools view -b | samtools sort > $OUTPUT/elba.contigs.bam && samtools index $OUTPUT/elba.contigs.bam
#minimap2 -x map-hifi -a -t64 $GENOME $READS | samtools view -b | samtools sort > $OUTPUT/elba.reads.bam && samtools index $OUTPUT/elba.reads.bam

#cp $GENOME $OUTPUT

#tar cvzf $TRANSFER/$OUTPUT.tar.gz $OUTPUT
