#!/bin/bash

#1 working dir
D="$1"
#D="/20tb/ratto/catcheR/oligo/"
echo "$D"
#2 sequences file
filename="$2"
#filename="shRNAs.csv"
#3 gibson 5' AGTTCCCTATCAGTGATAGAGATCCC
gibsonfive="$3"
#gibsonfive="AGTTCCCTATCAGTGATAGAGATCCC"
#4 fixed sequence between poly t and UCI
#GTCGACATTTAAATGGCGCGCC
fixed="$4"
#fixed="GTCGACATTTAAATGGCGCGCC"
#5 gbosn 3 GTAGCTCGCTGATCAGC
gibsonthree="$5"
#gibsonthree="GTAGCTCGCTGATCAGC"
#6 chech for restriction sites
restr="$6"
#restr="restr.txt"

full_path="$D/$filename"
echo $full_path
# Check if the file exists
if [ ! -f "$full_path" ]; then
    echo "File not found: $full_path"
    exit 1
fi

restr_file="$D/restr.txt"


output_file="$D/output.txt"
rm "$D/tmp.txt"

# Read each line from the file
while IFS=',' read -r line barcode gene; do
    #echo "$line"
    #echo "$barcode"
    # Select the first 52 characters of the line
    first_52=$(echo "$line" | cut -c 1-52 | rev | cut -c 1-48 | rev)
    #echo "$first_52"
    if [[ "${first_52:0:1}" != "A" && "${first_52:0:1}" != "G" ]]; then
        # Add "G" at the beginning
        first_52="G$first_52"
    fi
    
    # Add $gibsonfive before $first_52
    result=""
    result="$gibsonfive$first_52"
    
    # Add 7 T characters
    result+="TTTTTTT"
    
    #Add fixed
    result+="$fixed"
    
    # Add 6 random characters choosing between A, C, G, and T
    # Generate a random sequence of 6 numbers from 0 to 3
    #random_sequence=""
    #for ((i = 0; i < 6; i++)); do
    #  random_number=$((RANDOM % 4))  # Generate a random number from 0 to 3
    # case $random_number in
    #    0) random_sequence+="A";;
     #   1) random_sequence+="C";;
    #    2) random_sequence+="G";;
    #    3) random_sequence+="T";;
     # esac
    #done

    #echo "$random_sequence"
    #result+="$random_sequence"
    result+="NNNNNN"
    
    # Add $barcode
    result+="$barcode"
    
    # Add $gibsonthree
    result+="$gibsonthree"
    
    #add barcode and gene
    result+=,"$barcode","$gene"
    
    echo "$result" >> "$D/tmp.txt"
done < "$full_path"


#restriztion
if [ -f "$restr_file" ]; then
  while IFS= read -r subsequence; do
    sed -i "s/$subsequence/$(echo "$subsequence" | tr '[:upper:]' '[:lower:]')/gI" "$output_file"
  done < "$restr_file"
fi

grep -v -f "$restr_file" "$D/tmp.txt" > $output_file
grep -f "$restr_file" "$D/tmp.txt" > "$D/bad_oligos.txt"

chmod 777 /data/scratch/*
