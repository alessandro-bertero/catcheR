#!/bin/bash

#1 = folder
D="$1"
echo "$D"
#2 gene exp matrix
matrix=$2
#3 reference GGCGCGTTCATCTGGGGGAGCCG
ref=$3
#4 length of UCI
uci_len=$4
#5 threads
thread=$5
UMI_count=$6
percentage=$7
mode=$8

#DECOMPRESS AND GET FASTA
#gzip -d ./1st2nd_hiPSC_S3_R1_001.fastq.gz
echo "Decompressing"
gzip -fd "$D"/fastq/*
echo "Finished decompressing. Copying read 1 files..."
mkdir "$D"/fastq_R1
cp "$D"/fastq/*R1*.fastq "$D"/fastq_R1
echo "Finished copying read 1 files. Copying read 2 files..."
mkdir "$D"/fastq_R2
cp "$D"/fastq/*R2*.fastq "$D"/fastq_R2
echo "Finished copying read 2 files."
#rm -r "$D"/fasta
mkdir "$D"/fasta

echo "From fastq to fasta read1..."
for file in "$D"/fastq_R1/*; do
    if [ -f "$file" ]; then
        prefix=$(basename "$file" | cut -c1-3)
        output_file="$D/fasta/$prefix.R1.fa"

        # Skip processing if the output file already exists
        if [ -f "$output_file" ]; then
            echo "Output file $output_file already exists. Skipping file: $file"
            continue
        fi

        echo "Processing file: $file, Prefix: $prefix"
        cat "$file" | sed -n '1~4s/^@/>/p;2~4p' > "$output_file"
        echo "Converted file $file"
    else
        echo "Skipping non-regular file: $file"
    fi
done

echo "Done. From fastq to fasta read2..."
mkdir $D/fasta_R2
for file in "$D"/fastq_R2/*; do
    if [ -f "$file" ]; then
        prefix=$(basename "$file" | cut -c1-3)
        output_file="$D/fasta_R2/$prefix.R2.fa"

        # Skip processing if the output file already exists
        if [ -f "$output_file" ]; then
            echo "Output file $output_file already exists. Skipping file: $file"
            continue
        fi

        echo "Processing file: $file, Prefix: $prefix"
        cat "$file" | sed -n '1~4s/^@/>/p;2~4p' > "$output_file"
        echo "Converted file $file"
    else
        echo "Skipping non-regular file: $file"
    fi
done
echo "Done
"
wc -l "$D"/fasta_R2/* >> "$D"/log.txt

rm -r "$D"/fastq_R1
rm -r "$D"/fastq_R2

#GET SEQ WITH NO ID AND ADD DEMULTIPLEXING INFO TO READ ONE
#rm -r "$D"/txt
mkdir "$D"/txt
echo "Adding demultiplexing info to read 1..."

for file in "$D"/fasta/*; do
    if [ -f "$file" ]; then
        prefix=$(basename "$file" | cut -c1-3)
        output_file="$D/txt/$prefix.R1.txt"

        # Skip processing if the output file already exists
        if [ -f "$output_file" ]; then
            echo "Output file $output_file already exists. Skipping file: $file"
            continue
        fi

        echo "Processing file: $file, Prefix: $prefix"
        grep -v '>' "$file" >> "$output_file"
        sed -i "s/^/$prefix/" "$output_file"
        echo "Converted file $file"
    else
        echo "Skipping non-regular file: $file"
    fi
done
echo "Done. Now doing reverse complement of read2..."

#REVERSE COMPLEMENT R2 AND BARCODES
#cat "$D"/R2.txt | tr ACGTacgt TGCAtgca | rev >> "$D"/R2_rc.txt
mkdir "$D"/txt_R2/
for file in "$D"/fasta_R2/*; do
    if [ -f "$file" ]; then
        prefix=$(basename "$file" | cut -c1-3)
        output_file="$D/txt_R2/$prefix.R2.txt"

        # Skip processing if the output file already exists
        if [ -f "$output_file" ]; then
            echo "Output file $output_file already exists. Skipping file: $file"
            continue
        fi

        echo "Processing file: $file, Prefix: $prefix"
        grep -v '>' "$file" | tr ACGTacgt TGCAtgca | rev >> "$output_file"
        echo "Converted file $file"
    else
        echo "Skipping non-regular file: $file"
    fi
done
echo "Done. Now concatenating read1 and read2..."

rm -r "$D"/fasta
rm -r "$D"/fasta_R2

#CONCATENATE R1 AND R2
for file in "$D"/txt/*; do
    if [ -f "$file" ]; then
        prefix=$(basename "$file" | cut -c1-3)
        output_file="$D/$prefix.read1and2.txt"

        # Skip processing if the output file already exists
        if [ -f "$output_file" ]; then
            echo "Output file $output_file already exists. Skipping file: $file"
            continue
        fi

        echo "Processing file: $file, Prefix: $prefix"
        paste -d '' "$file" "$D"/txt_R2/$prefix.R2.txt > "$output_file"
        echo "Converted file $file"
    else
        echo "Skipping non-regular file: $file"
    fi
done
echo "Done! Now concatenating all files..."

#rm "$D"/Und.read1and2.txt
cat "$D"/*.read1and2.txt >> "$D"/read1and2.txt
rm "$D"/*.read1and2.txt
wc -l "$D"/read1and2.txt >> "$D"/log.txt
sort "$D"/read1and2.txt | uniq > "$D"/read1and2_uniq.txt
wc -l "$D"/read1and2_uniq.txt >> "$D"/log.txt

#EXTRACT CELL NAMES FROM 10X MATRIX
#head -n 1 "$D/$4" | tr ',' '\n' | tr -d '"' | tail -n +2 | tr -d '.' | tr -d '0123456789' | tr -d '-' | sed -e 's/^/^/' >> "$D"/names.txt
#wc -l "$D"/names.txt >> "$D"/log.txt
#sort ./names.txt | uniq > names_uniq.txt
#wc -l ./names_uniq.txt >> ./log.txt
#select only reads with cellID from 10x at the beginning. names.txt contains cellIDs from 10x and preceded by "^"
#Select with indeces cell ID and UMI
#grep --no-group-separator -f "$D"/names.txt "$D"/read1and2_uniq.txt > "$D"/with_cells_seq.txt
echo "Done! Now isolating cellID and UMI..."
cut -c 1-3,12-21 <<< cat "$D"/read1and2_uniq.txt >> "$D"/with_cells_cellID.txt
cut -c 4-11 <<< cat "$D"/read1and2_uniq.txt >> "$D"/with_cells_UMI.txt
#cut -c 1-28 <<< cat ./with_cells_seq.txt >> ./with_cells_cellID_UMI_tally.txt
sort "$D"/with_cells_cellID.txt | uniq > "$D"/with_cells_cellID_uniq.txt
echo "Done! Now isolating reference..."

#Reverse and select reference
ref_len="${#ref}"
cat "$D"/read1and2_uniq.txt | rev | cut -c 1-"$ref_len" | rev >> "$D"/with_cells_reference.txt
grep -v 'N' "$D"/with_cells_reference.txt > "$D"/with_cells_reference2.txt
#grep -v 'GGCGCGTTCATCTGGGGGAGCCG' ./with_cells_reference2.txt > ./with_cells_reference3.txt
grep -v $3 "$D"/with_cells_reference2.txt > "$D"/with_cells_reference3.txt

#Extract whole barcode (reference+UCI+barcode), UCI and barcode
#cat ./with_cells_seq.txt | rev | cut -c 1-37 | rev >> ./with_cells_whole_barcode.txt
echo "Done! Now isolating UCI and barcode..."
#UCI
beg=$(($ref_len + 1))
end=$(($beg + $uci_len - 1))
cat "$D"/read1and2_uniq.txt | rev | cut -c "$beg"-"$end" | rev >> "$D"/with_cells_UCI.txt

#barcode
bar_b=$(($end + 1))
IFS=',' read -r string _ < "$D"/rc_barcodes_genes.csv
length="${#string}"
bar_end=$(($bar_b + $length - 1))
cat "$D"/read1and2_uniq.txt | rev | cut -c "$bar_b"-"$bar_end" | rev >> "$D"/with_cells_barcode.txt
wc -l "$D"/with_cells*.txt >> "$D"/log.txt
#Extract barcode +-1
#cat ./with_cells_seq.txt | rev | cut -c 29-38 | rev >> ./with_cells_barcode_in_del.txt

#count matches
#grep -B 1 --no-group-separator -f ./rc_barcodes.txt ./with_cells_barcode.txt > ./perfect_match_barcodes.txt
#wc -l ./perfect_match_barcodes.txt
#grep -B 1 --no-group-separator -f ./rc_barcodes.txt ./with_cells_barcode_in_del.txt > ./in_del_barcodes.txt
#wc -l ./in_del_barcodes.txt

#run next step
echo "All done! Now running barcode_silencing_explorative_analysis.R ..."
#run next step - explorative
Rscript /home/sci_barcode_silencing_explorative_analysis.R "$D" "$ref" "$mode" "$matrix"
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
Rscript /home/barcode_silencing_cell_filtering.R "$D" 15 "$UMI_count" "$matrix"
wait

#run epmty selection
Rscript /home/barcode_silencing_empty_selection.R "$D" "$matrix"
wait

chmod 777 /data/scratch
