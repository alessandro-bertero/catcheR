#!/bin/bash

#1 working dir
D="$1"
echo "$D"
#2 read1 file
read1="$2"
#3 threshold percentage
DIs="$3"
echo "$DIs"
#4 clones
clones="$4"
echo "$clones"

#decompress files
gzip -fd "$D/$read1"
#gzip -fd "$D/$read2"

wc -l "$D/"*
#2,880,774
cat "$D/$read1" | sed -n '1~4s/^@/>/p;2~4p' >> "$D"/R1.fa

wc -l "$D"/R1.fa
#14,403,874

grep -v '>' "$D"/R1.fa >> "$D"/final.txt
wc -l "$D"/final.txt
#7,201,937

#cat "$D"/final.txt | sort -T /20tb/ratto >> final_sort.txt ?
cat "$D"/final.txt | sort -T  "$D" >>  "$D"/final_sort.txt
cat "$D"/final_sort.txt | uniq >> "$D"/final_uniq.txt
wc -l "$D"/final_uniq.txt
#7,086,995

cut -c 1-12 <<< cat "$D"/final.txt >> "$D"/final_UMI.txt
cut -c 30-37 <<< cat "$D"/final.txt >> "$D"/final_BC.txt
cut -c 38-43 <<< cat "$D"/final.txt >> "$D"/final_UCI.txt
cut -c 44-66 <<< cat "$D"/final.txt >> "$D"/final_reference.txt
#cut -c 69-93 <<< cat "$D"/final.txt >> "$D"/final_shRNA.txt

wc -l "$D"/final*.txt

 Rscript /home/plasmid_final2.R /data/scratch/ "$DIs" "$clones"