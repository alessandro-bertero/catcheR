#!/bin/bash

#1 working dir
D="$1"
echo "$D"
#2 read1 file
read1="$2"
#3 read2 file
read2="$3"
#4 threshold percentage
t_p="$4"
echo "$t_p"
#5 threshold 
t="$5"
echo "$t"
#5 clones
clones="$6"
echo "$clones"

# #decompress files
# gzip -fd "$D/$read1"
# gzip -fd "$D/$read2"
# 
# wc -l "$D/"*
# #2,880,774
# cat "$D/$read1" | sed -n '1~4s/^@/>/p;2~4p' >> "$D"/R1.fa
# 
# wc -l "$D"/R1.fa
# #14,403,874
# 
# grep -v '>' "$D"/R1.fa >> "$D"/final.txt
# wc -l "$D"/final.txt
# #7,201,937
# 
# #cat "$D"/final.txt | sort -T /20tb/ratto >> final_sort.txt ?
# cat "$D"/final.txt | sort -T  "$D" >>  "$D"/final_sort.txt
# cat "$D"/final_sort.txt | uniq >> "$D"/final_uniq.txt
# wc -l "$D"/final_uniq.txt
# #7,086,995
# 
# cut -c 1-12 <<< cat "$D"/final.txt >> "$D"/final_UMI.txt
# cut -c 30-37 <<< cat "$D"/final.txt >> "$D"/final_BC.txt
# cut -c 38-43 <<< cat "$D"/final.txt >> "$D"/final_UCI.txt
# #cut -c 69-93 <<< cat "$D"/final.txt >> "$D"/final_shRNA.txt
#EMPTY
cut -c 27-49 <<< cat "$D"/final_uniq.txt >> "$D"/final_reference_empty.txt
cut -c 44-66 <<< cat "$D"/final_uniq.txt >> "$D"/final_reference.txt
echo "Unique reads:"
wc -l "$D"/final_uniq.txt
echo "References:"
grep GGCGCGTTCATCTGGGGGAGCCG "$D"/final_reference.txt | wc -l
echo "Empty:"
grep TACGCGTTCATCTGGGGGAGCCG "$D"/final_reference_empty.txt | wc -l

wc -l "$D"/final*.txt

#Rscript /home/plasmid_final2.R /data/scratch/ "$t_p" "$t" "$clones"