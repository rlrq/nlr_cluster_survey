#!/bin/bash

gff_in=$1
bed_out=$2
func=$3

if [ -z "${func}" ]; then
    echo "leaving chromosomes as is"
    awk '{print $1 "\t" $4-1 "\t" $5 "\t" $2 " \t" $3 "\t" $6 " \t" $7 "\t" $8 "\t" $9}' $gff_in > $bed_out
elif [ $func == 1 ]; then
    echo "removing chromosome"
    awk '{print $1 "\t" $4-1 "\t" $5 "\t" $2 " \t" $3 "\t" $6 " \t" $7 "\t" $8 "\t" $9}' $gff_in | awk -F "\t" '{gsub("Chr","chr",$1)}1' OFS="\t" | awk -F "\t" '$5 != "chromosome" {print}' OFS="\t" > $bed_out
elif [ $func == 2 ]; then
    echo "retaining chromosome"
    # retain chromosome
    awk '{print $1 "\t" $4-1 "\t" $5 "\t" $2 " \t" $3 "\t" $6 " \t" $7 "\t" $8 "\t" $9}' $gff_in | awk -F "\t" '{gsub("Chr","chr",$1)}1' OFS="\t" > $bed_out
fi
