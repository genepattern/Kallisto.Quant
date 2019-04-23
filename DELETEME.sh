#!/bin/sh

# Reading each line in the file list
filename=$1
list_of_files=""
while read line; do
  # reading each line
  #echo $line
  list_of_files="$list_of_files $line"
done < $filename
echo "list of files provided:"
echo $list_of_files
