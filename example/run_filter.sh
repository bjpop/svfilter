#!/bin/bash

DIR=/scratch/VR0002/bjpop/colorectal_lynch_syndrome_buchanan/

SAMPLES="0131313009 0151081002 0636307011 0131326001 0151095001 0656045001 0757045010 E60176 0151032037 0350124001 0656050001 9930087001 0151052001 0350355001 0757017010 0151078001 0636307001 0757045001 C4055110001"

for sample in $SAMPLES; do
	input="$DIR/$sample/${sample}.lumpy.vcf"
	output="${sample}.lumpy.filtered.vcf"
	echo $input
	echo $output
	svfilter --coords genes_promoter.bed --type lumpy $input > $output
done
