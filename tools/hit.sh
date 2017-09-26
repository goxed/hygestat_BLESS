#!/usr/bin/env bash
EXPECTED_ARGS=2
#EXPECTED_ARGS1=2

if [ $# -ne $EXPECTED_ARGS ]
then
  exit
fi

if [ $# -ne $EXPECTED_ARGS ]
then
  exit
fi

if [ ! -z "$1" ] && echo "Argument #1 = $1" && [ ! -z "$2" ] && echo "Argument #2 = $2"
then
  echo "At least 2 arguments passed to script."
else
  echo "Less than 2 arguments passed to the script."
fi

if [ $2 == "mouse" ]
then

for((x=1;x<=19;x++));
  do  
    cat ${1}  | grep + | grep  -P "chr${x}\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr${x}+.hits.txt &       
  done
    wait
    cat ${1}  | grep + | grep  -P "chrX\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrX+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrY\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrY+.hits.txt        

for((x=1;x<=19;x++));
  do
    echo chr${x}+ >> chr_all_plus.txt
    tail -3 chr${x}+.hits.txt >> chr_all_plus.txt
 done
    echo chrX+ >> chr_all_plus.txt
    tail -3 chrX+.hits.txt >> chr_all_plus.txt
    echo chrY+ >> chr_all_plus.txt
    tail -3 chrY+.hits.txt >> chr_all_plus.txt

for((x=1;x<=19;x++));
  do  
    cat ${1}  | grep - | grep  -P "chr${x}\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr${x}-.hits.txt &       
  done
  wait
  cat ${1}  | grep - | grep  -P "chrX\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrX-.hits.txt        
  cat ${1}  | grep - | grep  -P "chrY\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrY-.hits.txt        


for((x=1;x<=19;x++));
  do
    echo chr${x}- >> chr_all_minus.txt
    tail -3 chr${x}-.hits.txt >> chr_all_minus.txt
 done
    echo chrX- >> chr_all_minus.txt
    tail -3 chrX-.hits.txt >> chr_all_minus.txt
    echo chrY- >> chr_all_minus.txt
    tail -3 chrY-.hits.txt >> chr_all_minus.txt

elif [ $2 == "human" ]
then

for((x=1;x<=22;x++));
  do
    cat ${1}  | grep + | grep  -P "chr${x}\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr${x}+.hits.txt &
  done
    wait
    cat ${1}  | grep + | grep  -P "chrX\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrX+.hits.txt
    cat ${1}  | grep + | grep  -P "chrY\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrY+.hits.txt

for((x=1;x<=22;x++));
  do
    echo chr${x}+ >> chr_all_plus.txt
    tail -3 chr${x}+.hits.txt >> chr_all_plus.txt
 done
    echo chrX+ >> chr_all_plus.txt
    tail -3 chrX+.hits.txt >> chr_all_plus.txt
    echo chrY+ >> chr_all_plus.txt
    tail -3 chrY+.hits.txt >> chr_all_plus.txt

for((x=1;x<=22;x++));
  do
    cat ${1}  | grep - | grep  -P "chr${x}\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr${x}-.hits.txt &
  done
  wait
  cat ${1}  | grep - | grep  -P "chrX\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrX-.hits.txt
  cat ${1}  | grep - | grep  -P "chrY\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrY-.hits.txt


for((x=1;x<=22;x++));
  do
    echo chr${x}- >> chr_all_minus.txt
    tail -3 chr${x}-.hits.txt >> chr_all_minus.txt
 done
    echo chrX- >> chr_all_minus.txt
    tail -3 chrX-.hits.txt >> chr_all_minus.txt
    echo chrY- >> chr_all_minus.txt
    tail -3 chrY-.hits.txt >> chr_all_minus.txt

elif [ $2 == "yeast" ]
then

    cat ${1}  | grep + | grep  -P "chrI\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr1+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr2+.hits.txt  
    cat ${1}  | grep + | grep  -P "chrIII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr3+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrIV\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr4+.hits.txt
    cat ${1}  | grep + | grep  -P "chrV\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr5+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrVI\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr6+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrVII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr7+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrVIII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr8+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrIX\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr9+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrX\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr10+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrXI\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr11+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrXII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr12+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrXIII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr13+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrXIV\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr14+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrXV\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr15+.hits.txt        
    cat ${1}  | grep + | grep  -P "chrXVI\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr16+.hits.txt
    cat ${1}  | grep + | grep  -P "chrM\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrM+.hits.txt                
        
    echo chr1+ >> chr_all_plus.txt
    tail -3 chr1+.hits.txt >> chr_all_plus.txt
    echo chr2+ >> chr_all_plus.txt
    tail -3 chr2+.hits.txt >> chr_all_plus.txt
    echo chr3+ >> chr_all_plus.txt
    tail -3 chr3+.hits.txt >> chr_all_plus.txt
    echo chr4+ >> chr_all_plus.txt
    tail -3 chr4+.hits.txt >> chr_all_plus.txt
    echo chr5+ >> chr_all_plus.txt
    tail -3 chr5+.hits.txt >> chr_all_plus.txt
    echo chr6+ >> chr_all_plus.txt
    tail -3 chr6+.hits.txt >> chr_all_plus.txt
    echo chr7+ >> chr_all_plus.txt
    tail -3 chr7+.hits.txt >> chr_all_plus.txt
    echo chr8+ >> chr_all_plus.txt
    tail -3 chr8+.hits.txt >> chr_all_plus.txt
    echo chr9+ >> chr_all_plus.txt
    tail -3 chr9+.hits.txt >> chr_all_plus.txt
    echo chr10+ >> chr_all_plus.txt
    tail -3 chr10+.hits.txt >> chr_all_plus.txt
    echo chr11+ >> chr_all_plus.txt
    tail -3 chr11+.hits.txt >> chr_all_plus.txt
    echo chr12+ >> chr_all_plus.txt
    tail -3 chr12+.hits.txt >> chr_all_plus.txt
    echo chr13+ >> chr_all_plus.txt
    tail -3 chr13+.hits.txt >> chr_all_plus.txt
    echo chr14+ >> chr_all_plus.txt
    tail -3 chr14+.hits.txt >> chr_all_plus.txt
    echo chr15+ >> chr_all_plus.txt
    tail -3 chr15+.hits.txt >> chr_all_plus.txt
    echo chr16+ >> chr_all_plus.txt
    tail -3 chr16+.hits.txt >> chr_all_plus.txt
    echo chrM+ >> chr_all_plus.txT
    tail -3 chrM+.hits.txt >> chr_all_plus.txt

    cat ${1}  | grep - | grep  -P "chrI\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr1-.hits.txt
    cat ${1}  | grep - | grep  -P "chrII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr2-.hits.txt
    cat ${1}  | grep - | grep  -P "chrIII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr3-.hits.txt
    cat ${1}  | grep - | grep  -P "chrIV\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr4-.hits.txt
    cat ${1}  | grep - | grep  -P "chrV\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr5-.hits.txt
    cat ${1}  | grep - | grep  -P "chrVI\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr6-.hits.txt
    cat ${1}  | grep - | grep  -P "chrVII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr7-.hits.txt
    cat ${1}  | grep - | grep  -P "chrVIII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr8-.hits.txt
    cat ${1}  | grep - | grep  -P "chrIX\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr9-.hits.txt
    cat ${1}  | grep - | grep  -P "chrX\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr10-.hits.txt
    cat ${1}  | grep - | grep  -P "chrXI\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr11-.hits.txt
    cat ${1}  | grep - | grep  -P "chrXII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr12-.hits.txt
    cat ${1}  | grep - | grep  -P "chrXIII\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr13-.hits.txt
    cat ${1}  | grep - | grep  -P "chrXIV\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr14-.hits.txt
    cat ${1}  | grep - | grep  -P "chrXV\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr15-.hits.txt
    cat ${1}  | grep - | grep  -P "chrXVI\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chr16-.hits.txt
    cat ${1}  | grep - | grep  -P "chrM\t" | awk '{print $4}' | sort -g | uniq -c | sort -g -k+1 > chrM-.hits.txt

    echo chr1- >> chr_all_plus.txt
    tail -3 chr1-.hits.txt >> chr_all_plus.txt
    echo chr2- >> chr_all_plus.txt
    tail -3 chr2-.hits.txt >> chr_all_plus.txt
    echo chr3- >> chr_all_plus.txt
    tail -3 chr3-.hits.txt >> chr_all_plus.txt
    echo chr4- >> chr_all_plus.txt
    tail -3 chr4-.hits.txt >> chr_all_plus.txt
    echo chr5- >> chr_all_plus.txt
    tail -3 chr5-.hits.txt >> chr_all_plus.txt
    echo chr6- >> chr_all_plus.txt
    tail -3 chr6-.hits.txt >> chr_all_plus.txt
    echo chr7- >> chr_all_plus.txt
    tail -3 chr7-.hits.txt >> chr_all_plus.txt
    echo chr8- >> chr_all_plus.txt
    tail -3 chr8-.hits.txt >> chr_all_plus.txt
    echo chr9- >> chr_all_plus.txt
    tail -3 chr9-.hits.txt >> chr_all_plus.txt
    echo chr10- >> chr_all_plus.txt
    tail -3 chr10-.hits.txt >> chr_all_plus.txt
    echo chr11- >> chr_all_plus.txt
    tail -3 chr11-.hits.txt >> chr_all_plus.txt
    echo chr12- >> chr_all_plus.txt
    tail -3 chr12-.hits.txt >> chr_all_plus.txt
    echo chr13- >> chr_all_plus.txt
    tail -3 chr13-.hits.txt >> chr_all_plus.txt
    echo chr14- >> chr_all_plus.txt
    tail -3 chr14-.hits.txt >> chr_all_plus.txt
    echo chr15- >> chr_all_plus.txt
    tail -3 chr15-.hits.txt >> chr_all_plus.txt
    echo chr16- >> chr_all_plus.txt
    tail -3 chr16-.hits.txt >> chr_all_plus.txt
    echo chrM- >> chr_all_plus.txT
    tail -3 chrM-.hits.txt >> chr_all_plus.txt

else 
    echo "No matching name"
fi

