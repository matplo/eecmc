#!/bin/bash

file_list=$1
if [[ -z "$file_list" ]]; then
	echo "Usage: $0 <file_list>"
	exit 1
fi

total_size=0
flist=$(cat ${file_list})
nfiles=0
files_to_hadd=""
for file in $flist
do
  if [[ -f "$file" ]]; then
    file_size=$(du -b "$file" | cut -f1)
    total_size=$((total_size + file_size))
		nfiles=$((nfiles + 1))

		# no check for now
		files_to_hadd="${files_to_hadd} ${file}"
		# check_root_file.py ${file}
		# if [ $? -eq 0 ]; then
		# 	echo "File ${file} OK"
		# 	files_to_hadd="${files_to_hadd} ${file}"
		# else
		# 	echo "* File ${file} not OK"
		# fi

	else
		echo "File $file not found"
  fi
done
echo "Number of files: $nfiles"
echo "Total size: $total_size bytes"
total_size_mb=$(echo "scale=2; $total_size / 1048576" | bc)
echo "Total size: $total_size_mb MB"

if [ ! -z "${files_to_hadd}" ]; then
	echo "hadd -f $(dirname ${file_list})/hadded_$(basename ${file_list} .root).root ${files_to_hadd}"
	hadd -f $(dirname ${file_list})/hadded_$(basename ${file_list} .root).root ${files_to_hadd}
fi
