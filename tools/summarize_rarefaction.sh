#!/usr/bin/bash
inp=$1
file_no=$(echo "$1" | awk -F '_' '{print $3}')

dir=$(ls | grep _$file_no| grep output_close)
dir=(${dir//" "/})
files="${dir[0]}/output.txt"
var=""
tes=8
size=${#dir[@]}
size=$((size-1))
for i in `seq 1 $size`
 do
    files=$files" ${dir[$i]}/output.txt"
    if [ $tes -ne 8 ]
      then
    	var=$var"",$tes,$((tes +1)),$((tes +2)),$((tes +4))""
    else
    	var=$var""$tes,$((tes +1)),$((tes +2)),$((tes +4))""
    fi
    tes=$((tes + 7))
done
paste `  echo $files` | cut -f $var --complement > "/tmp/"$inp"_summary.txt"
