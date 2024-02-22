#!/bin/bash
inpath=$1
Chromosome_file=$2
start_seq=$3
end_seq=$4
repeat_length=$5
linecount=$6

pathname=$(pwd)

cd $pathname

line_start=$(echo $((start_seq/linecount)))
line_end=$(echo $((end_seq/linecount)))

# extract this part of the file
sed -n "${line_start},${line_end}p" ${inpath}/${Chromosome_file} > ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.fasta

echo  ${inpath}/${Chromosome_file}

echo $line_start
echo $line_end

# load emboss and search for repeats of that length
#module load nixpkgs/16.09  gcc/7.3.0 emboss/6.6.0
equicktandem -sequence ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.fasta -maxrepeat ${repeat_length} -threshold 500 -outfile out.txt

mv out.txt ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out 

sed -n "23,27p" ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out
#  Start     End  Strand   Score   Size  Count
#  93852   99234       +    1391      4   1345

# line 24 is first repeat
linep=$(grep -n " ${repeat_length} " ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out| cut -d ':' -f 1)
echo $linep

# extract values from output
emb_start=$(sed -n "${linep},${linep}p" ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out  | tr -s " " "_" | cut -d '_' -f 2)
echo $emb_start

emb_end=$(sed -n "${linep},${linep}p" ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out | tr -s " " "_" | cut -d '_' -f 3)
echo $emb_end

#repeatlength=$(sed -n "${linep}p" ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out | tr -s " " "_" | cut -d '_' -f 6)
#echo $repeatlength

# get line numbers associated with character positions
line_emb_start=$(echo $((emb_start/linecount)))
line_emb_end=$(echo $((emb_end/linecount)))

# check fasta at this line
sed -n "${line_emb_start},${line_emb_end}p" ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.fasta > ${pathname}/${Chromosome_file}_${repeat_length}_${emb_start}-${emb_end}.fasta

# print the various values
echo $(echo $(cat ${pathname}/${Chromosome_file}_${repeat_length}_${emb_start}-${emb_end}.fasta))
echo $(sed -n "4,4p" ${pathname}/${Chromosome_file}_${repeat_length}_${emb_start}-${emb_end}.fasta)
echo $(sed -n "4,4p" ${pathname}/${Chromosome_file}_${repeat_length}_${emb_start}-${emb_end}.fasta |cut -c1-${repeat_length})

