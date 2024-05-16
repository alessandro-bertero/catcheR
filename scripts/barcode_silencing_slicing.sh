#!/bin/bash

#1 = folder
D="$1"
echo "$D"
#2 fastq r1
#3 fastq r2
#4 gene exp matrix
#5 reference GGCGCGTTCATCTGGGGGAGCCG
ref=$5
#6 length of UCI
uci_len=$6
#7 threads
thread=$7
UMI_count=$8
percentage=$9
mode=$10


#DECOMPRESS AND GET FASTA 
#gzip -d ./1st2nd_hiPSC_S3_R1_001.fastq.gz
gzip -fd "$D/$2" 
gzip -fd "$D/$3"  
# Check if the filename ends with ".gz"
if [[ $2 == *.gz ]]; then
    # Remove ".gz" from the filename
    filename1="${2%.gz}"
else
    filename1="$2"
fi

echo "$filename1"

# Check if the filename ends with ".gz"
if [[ $3 == *.gz ]]; then
    # Remove ".gz" from the filename
    filename2="${3%.gz}"
else
    filename2="$3"
fi

echo "$filename2"

cat "$D/$filename1" | sed -n '1~4s/^@/>/p;2~4p' >> "$D"/R1.fa
cat "$D/$filename2" | sed -n '1~4s/^@/>/p;2~4p' > "$D"/R2.fa

#gzip -d ./1st2nd_hiPSC_S3_R2_001.fastq.gz
wc -l "$D"/* >> "$D"/log.txt


#GET SEQ WITH NO ID
#paste ./R1.fa ./R2.fa > read1and2_check.fasta
grep -v '>' "$D"/R1.fa >> "$D"/R1.txt
grep -v '>' "$D"/R2.fa >> "$D"/R2.txt


#REVERSE COMPLEMENT R2 AND BARCODES
cat "$D"/R2.txt | tr ACGTacgt TGCAtgca >> "$D"/R2_c.txt
cat "$D"/R2_c.txt | rev >> "$D"/R2_rc.txt
wc -l "$D"/R*.txt >> "$D"/log.txt

#CONCATENATE R1 AND R2
paste -d '' "$D"/R1.txt "$D"/R2_rc.txt > "$D"/read1and2.txt
wc -l "$D"/read1and2.txt >> "$D"/log.txt
sort "$D"/read1and2.txt | uniq -c > "$D"/read1and2_uniq.txt
wc -l "$D"/read1and2_uniq.txt >> "$D"/log.txt


#EXTRACT CELL NAMES FROM 10X MATRIX 
#select only reads with cellID from 10x at the beginning. names.txt contains cellIDs from 10x and preceded by "^"
head -n 1 "$D/$4" | tr ',' '\n' | tr -d '"' | tail -n +2 | tr -d '.' | tr -d '0123456789' | tr -d '-' | sed -e 's/^/^[0-9]{1,3}\\s/' >> "$D"/names.txt
#'s/^/^[0-9]{1,2}\\s/'
wc -l "$D"/names.txt >> "$D"/log.txt
split -l$((`wc -l < "$D"/names.txt`/$thread)) "$D"/names.txt "$D"/names.split.txt
#sort ./names.txt | uniq > names_uniq.txt
#wc -l ./names_uniq.txt >> ./log.txt
#select only reads with cellID from 10x at the beginning. names.txt contains cellIDs from 10x and preceded by "^"
#Select with indeces cell ID and UMI
sed 's/^[ \t]*//' "$D"/read1and2_uniq.txt >> "$D"/read1and2_uniq2.txt
for file in $(find "$D" -type f -name 'names.split*'); do
    /home/grep.sh "$file" "$D"/read1and2_uniq2.txt &
done

wait

rm "$D"/read1and2_uniq.txt
cat "$D"/names.split.txt*.out >> "$D"/with_cells_seq_c.txt
awk '{print length}' "$D"/with_cells_seq_c.txt | sort -n | uniq -c
#tr -s '[:space:]' < "$D"/with_cells_seq_ctmp.txt | sed '/^\s*$/d' > "$D"/with_cells_seq_c.txt
rm "$D"/names.split.txt*
awk '{print $2}' "$D"/with_cells_seq_c.txt >> "$D"/with_cells_seq.txt
awk '{print $1}' "$D"/with_cells_seq_c.txt >> "$D"/with_cells_reads_counts.txt
cut -c 1-16 <<< cat "$D"/with_cells_seq.txt >> "$D"/with_cells_cellID.txt
cut -c 17-28 <<< cat "$D"/with_cells_seq.txt >> "$D"/with_cells_UMI.txt
sort "$D"/with_cells_cellID.txt | uniq > "$D"/with_cells_cellID_uniq.txt

#Reverse and select reference
ref_len="${#ref}"
cat "$D"/with_cells_seq.txt | rev | cut -c 1-"$ref_len" | rev >> "$D"/with_cells_reference.txt
grep -v 'N' "$D"/with_cells_reference.txt > "$D"/with_cells_reference2.txt
grep -v $6 "$D"/with_cells_reference2.txt > "$D"/with_cells_reference3.txt


#UCI
beg=$(($ref_len + 1))
end=$(($beg + $uci_len - 1))
cat "$D"/with_cells_seq.txt | rev | cut -c "$beg"-"$end" | rev >> "$D"/with_cells_UCI.txt

#BARCODES
bar_b=$(($end + 1)) #starts at end of UCI
#extract barcode length from file
IFS=',' read -r string _ < "$D"/rc_barcodes_genes.csv
length="${#string}"
bar_end=$(($bar_b + $length - 1))
cat "$D"/with_cells_seq.txt | rev | cut -c "$bar_b"-"$bar_end" | rev >> "$D"/with_cells_barcode.txt
wc -l "$D"/with_cells*.txt >> "$D"/log.txt
#Extract barcode +-1
#cat ./with_cells_seq.txt | rev | cut -c 29-38 | rev >> ./with_cells_barcode_in_del.txt

#count matches
#grep -B 1 --no-group-separator -f ./rc_barcodes.txt ./with_cells_barcode.txt > ./perfect_match_barcodes.txt
#wc -l ./perfect_match_barcodes.txt
#grep -B 1 --no-group-separator -f ./rc_barcodes.txt ./with_cells_barcode_in_del.txt > ./in_del_barcodes.txt
#wc -l ./in_del_barcodes.txt

#run next step - explorative
Rscript /home/barcode_silencing_explorative_analysis.R "$D" "$ref" "$mode"
wait

#run cell filtering
if [ -z "$UMI_count" ]; then
    echo "\$UMI_count is either unset or null: using suggested threshold."
    filename="/data/scratch/UMI_threshold.txt"
    read -r UMI_count < "$filename"
else
    echo "\$UMI_count is not null."
    echo "$UMI_count"
fi
Rscript /home/barcode_silencing_cell_filtering.R "$D" 15 "$UMI_count" "$4"
wait

#run epmty selection
Rscript /home/barcode_silencing_empty_selection.R "$D" "$4" "$UMI_count"
wait

chmod 777 /data/scratch