#!/bin/bash
inpath=$1
Chromosome_file=$2
start_seq=$3
end_seq=$4
repeat_length=$5
linecount=$6
threshold=$7

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
equicktandem -sequence ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.fasta -maxrepeat ${repeat_length} -threshold ${threshold} -outfile out.txt

mv out.txt ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out 

sed -n "23,$p" ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out
#  Start     End  Strand   Score   Size  Count
#  93852   99234       +    1391      4   1345

# line -4 is last largest repeat - so don't search 
#linep=$(grep -n " ${repeat_length} " ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out| cut -d ':' -f 1)
#echo $linep

# extract values from output
emb_start=$(tail -n4 ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out | head -n1 |  tr -s " " "_" | cut -d '_' -f 2)
echo $emb_start

emb_end=$(tail -n4 ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out | head -n1 |  tr -s " " "_" | cut -d '_' -f 3)
echo $emb_end

repeatlength=$(tail -n4 ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.out | head -n1 |  tr -s " " "_" | cut -d '_' -f 6)
echo $repeatlength

# get line numbers associated with character positions
line_emb_start=$(echo $(((emb_start/linecount)+1)))
line_emb_end=$(echo $((emb_end/linecount)))

# check fasta at this line
sed -n "${line_emb_start},${line_emb_end}p" ${pathname}/${Chromosome_file}_${start_seq}-${end_seq}.fasta > ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}.fasta

# add header
sed -e "1i\>${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}" ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}.fasta > ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}_1.fasta

#module load StdEnv/2020 seqkit/2.3.1
seqkit seq -w ${repeatlength} ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}_1.fasta > ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}_2.fasta

# print the various values
echo $(echo "$(<${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}_2.fasta)")

#############
# future ideas to print the whole repeat if ${repeatlength} > ${linecount}

# if [ ${repeatlength} < ${linecount} ]
# then
# echo $(sed -n "4,4p" ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}.fasta |cut -c1-${repeatlength})
# else
# line2count=$(echo $((2*${linecount})))
# if [ ${repeatlength} < ${line2count}]
# then
# lengthdiff=$(echo $((-(${repeatlength}-${line2count}))))
# echo $(sed -n "4,4p" ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}.fasta)
# echo $(sed -n "5,5p" ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}.fasta |cut -c1-${lengthdiff})
# else
# line3count=$(echo $((3*${linecount})))
# if [ ${repeatlength} < ${line3count}]
# then
# lengthdiff=$(echo $((-(${repeatlength}-${line3count}))))
# echo $(sed -n "4,4p" ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}.fasta)
# echo $(sed -n "5,5p" ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}.fasta)
# echo $(sed -n "6,6p" ${pathname}/${Chromosome_file}_${repeatlength}_${emb_start}-${emb_end}.fasta |cut -c1-${lengthdiff})
# fi
# fi 
# fi
