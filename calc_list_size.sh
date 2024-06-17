#!/bin/bash

file_list=$1
if [[ -z "$file_list" ]]; then
	echo "Usage: $0 <file_list>"
	exit 1
fi

total_size=0
flist=$(cat ${file_list})
nfiles=0
for file in $flist
do
  if [[ -f "$file" ]]; then
    file_size=$(du -b "$file" | cut -f1)
    total_size=$((total_size + file_size))
		nfiles=$((nfiles + 1))
	else
		echo "File $file not found"
  fi
done
echo "Number of files: $nfiles"
echo "Total size: $total_size bytes"
total_size_mb=$(echo "scale=2; $total_size / 1048576" | bc)
echo "Total size: $total_size_mb MB"