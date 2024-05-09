#!/bin/bash

#1 working dir
D="$1"
echo "$D"
#2 read1 file
read1="$2"
#3 read2 file
read2="$3"
#4 threshold 
t="$4"
echo "$t"
#5 clones
clones="$5"
echo "$clones"

# #decompress files
# gzip -fd "$D/$read1"
# gzip -fd "$D/$read2"
# 
# wc -l "$D/"*
# cat "$D/$read1" | sed -n '1~4s/^@/>/p;2~4p' >> "$D"/R1.fa
# 
# wc -l "$D"/R1.fa
# 
# grep -v '>' "$D"/R1.fa >> "$D"/inter.txt
# wc -l "$D"/inter.txt
# cat "$D"/inter.txt | sort -T "$D" >> "$D"/inter_sort.txt
# cat "$D"/inter_sort.txt | uniq >> "$D"/inter_uniq.txt
# wc -l "$D"/inter_uniq.txt
# 
# cut -c 1-12 <<< cat "$D"/inter.txt >> "$D"/inter_UMI.txt
# cut -c 30-37 <<< cat "$D"/inter.txt >> "$D"/inter_BC.txt
# cut -c 38-43 <<< cat "$D"/inter.txt >> "$D"/inter_UCI.txt
# cut -c 69-93 <<< cat "$D"/inter.txt >> "$D"/inter_shRNA.txt
# 
# wc -l "$D"/inter*.txt

Rscript /home/plasmid_inter2.R /data/scratch "$t" "$clones"