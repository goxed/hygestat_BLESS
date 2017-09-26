#!/usr/bin/sh

file=$1;
chr=`less $file| awk '{print $3}' | sort | uniq`
echo "chromosome" "number_of_reads">$file."chr_dist"
for i in $chr
  do
    reads=`grep -c $i $file`
    echo $i $reads >> $file."chr_dist"
  done
