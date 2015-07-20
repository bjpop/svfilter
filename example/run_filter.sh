#!/bin/bash

DIR=/XXX/YYY

SAMPLES="A B C D"

for sample in $SAMPLES; do
	input="$DIR/$sample/${sample}.lumpy.vcf"
	output="lumpy.filtered.vcf"
	#svfilter --coords mmr_genes.bed --type vcf $input >> $output
	svfilter --coords genes.bed --type vcf --custom custom_filter.py $input >> $output
done
