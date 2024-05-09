#!/bin/bash

file=$1
reads=$2

# Your script logic using $file as input
echo "Processing file: $file"
# Add your script logic here
grep --no-group-separator -E -f $file $reads >> "${file}.out"
#grep --no-group-separator -E -f $file "$D"/read1and2_uniq2.txt >> "$D"/with_cells_seq_c.txt