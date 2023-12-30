#################################
# Notes for running R code on the server
# June 2023 
# Repeat runs for centromere paper
###################################
# to run with default script

cd ~/scratch/repeats/auto_script/COLCEN

# download the Setup_Run_Repeats.sh
wget https://github.com/celphin/RepeatOBserverV1/blob/main/Setup_Run_Repeats.sh
chmod +x Setup_Run_Repeats.sh

# download your genome
wget https://www.arabidopsis.org/download_files/Genes/Col-CEN%20genome%20assembly%20release/ColCEN.fasta


cat << EOF > Arabidopsis_repeats.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=128000M

module load StdEnv/2020
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i COLCEN -f ColCEN.fasta -h H0 -c 15 -mem 128000M

EOF

sbatch Arabidopsis_repeats.sh

#--------------------------
# the following code is how each genome was run for the manuscript
# all the following code can be used as an example if the automatic script is not working properly 

# Arabidopsis
# https://www.biorxiv.org/content/10.1101/2021.05.30.446350v2.full

# 5 chromosomes

cd /home/celphin/scratch/repeats/input_chromosomes/Arab_COLCEN

# download link from: https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FCol-CEN+genome+assembly+release
wget https://www.arabidopsis.org/download_files/Genes/Col-CEN%20genome%20assembly%20release/ColCEN.fasta

module load seqkit/2.3.1
seqkit seq -w 60 ColCEN.fasta > ColCEN_2.fasta

# Move each sequence into its own fasta      
mkdir chromosome_files/ 
cd ./chromosome_files/                                                                                                                                                                
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../ColCEN_2.fasta

rename Chr COLCEN_H0_Chr *


##############################
# Get other genomes from ncbi and other websites

# Orchid
# https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.15306

# https://www.ncbi.nlm.nih.gov/genome/41611
# https://www.ncbi.nlm.nih.gov/genome/17836

# 19 chromosomes, 2 species, 
#Dendrobium huoshanense
#Dendrobium nobile

cd /home/celphin/scratch/repeats/input_chromosomes/Orchid

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/618/105/GCA_016618105.1_ASM1661810v1/GCA_016618105.1_ASM1661810v1_genomic.fna.gz #Dhuo
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/539/455/GCA_022539455.1_Dnobile/GCA_022539455.1_Dnobile_genomic.fna.gz #Dnob

# Unzip
gunzip GCA_016618105.1_ASM1661810v1_genomic.fna.gz
gunzip GCA_022539455.1_Dnobile_genomic.fna.gz  

module load seqkit/2.3.1
seqkit seq -w 60 GCA_016618105.1_ASM1661810v1_genomic.fna > Dhuo_2.fna
seqkit seq -w 60 GCA_022539455.1_Dnobile_genomic.fna > Dnob_2.fna

# Move each sequence into its own fasta         
mkdir chromosome_files
cd chromosome_files      
                                                                                                                                                        
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Dhuo_2.fna
# remove contigs
find -type f -size -10000000c -delete
# rename
rename CM02834 Dhuo_H0_Chr0 *
rename CM02835 Dhuo_H0_Chr1 *
rename CM02836 Dhuo_H0_Chr2 *
rename .1.fasta .fasta *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Dnob_2.fna
# remove contigs
find -type f -size -10000000c -delete
# rename
rename CM03971 Dnob_H0_Chr0 *
rename CM03972 Dnob_H0_Chr1 *
rename CM03973 Dnob_H0_Chr2 *
rename .1.fasta .fasta *

rename Dhuo_H0_Chr Orchid_H0_ChrDhuo *
rename Dnob_H0_Chr Orchid_H0_ChrDnob *

########################
Dhuo
Chr ID
1	42
2	43
3	44
4	45
5	46
6	47
7	48
8	49
9	50
10	51
11	52
12	53
13	54
14	55
15	56
16	57
17	58
18	59
19	60



########################
Dnob
Chr ID
1	18
2	19
3	20
4	21
5	22
6	23
7	24
8	25
9	26
10	27
11	28
12	29
13	30
14	31
15	32
16	33
17	34
18	35
19	36
	


#####################################
# sunflower HA89 and HA 412
#------------------------------
cd /home/celphin/projects/rpp-rieseber/celphin/sunflower/HA412/chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Ha412HOv2.0-20181130.fasta
rm -v Ha412HOChr00c*.fasta
rename Ha412HOChr Ha412_HO_Chr *

#----------------------
cd /home/celphin/projects/rpp-rieseber/celphin/sunflower/HA89/
# need to fold HA89/chromosome_files/
module load seqkit/2.3.1
seqkit seq -w 60 HA89.fasta > HA89_2.fasta

mkdir chromosome_files/
cd chromosome_files
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../HA89_2.fasta
rm -v HanHA89Chr00c*.fasta
rename HanHA89Chr HA89_H0_Chr *


#########################################################
# Fern - 39 chromosomes, 7 gbp, N50 2Mbp
# https://www.ncbi.nlm.nih.gov/assembly/GCA_020310875.1

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Fern

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/020/310/875/GCA_020310875.1_C.richardii_v2/GCA_020310875.1_C.richardii_v2_genomic.fna.gz

# Unzip
gunzip GCA_020310875.1_C.richardii_v2_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_020310875.1_C.richardii_v2_genomic.fna > Fern_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Fern_2.fasta

# remove contigs
find -type f -size -10000000c -delete

rename CM0354 Fern_H0_Chr *
rename .1.fasta .fasta *

# note -6 to get true chromosome names

wc -c *
262 351 372 Fern_H0_Chr06.fasta
 230912062 Fern_H0_Chr07.fasta
 230892866 Fern_H0_Chr08.fasta
 222352590 Fern_H0_Chr09.fasta
 213963964 Fern_H0_Chr10.fasta
 210423836 Fern_H0_Chr11.fasta
 204903031 Fern_H0_Chr12.fasta
 203840995 Fern_H0_Chr13.fasta
 203832987 Fern_H0_Chr14.fasta
 202230780 Fern_H0_Chr15.fasta
 202132940 Fern_H0_Chr16.fasta
 201506194 Fern_H0_Chr17.fasta
 200566167 Fern_H0_Chr18.fasta
 196757974 Fern_H0_Chr19.fasta
 193484094 Fern_H0_Chr20.fasta
 190514130 Fern_H0_Chr21.fasta
 189143611 Fern_H0_Chr22.fasta
 185694338 Fern_H0_Chr23.fasta
 184995749 Fern_H0_Chr24.fasta
 184715226 Fern_H0_Chr25.fasta
 181549257 Fern_H0_Chr26.fasta
 181086115 Fern_H0_Chr27.fasta
 179160593 Fern_H0_Chr28.fasta
 168457847 Fern_H0_Chr29.fasta
 167638652 Fern_H0_Chr30.fasta
 165088357 Fern_H0_Chr31.fasta
 161965825 Fern_H0_Chr32.fasta
 159556217 Fern_H0_Chr33.fasta
 158470585 Fern_H0_Chr34.fasta
 158036985 Fern_H0_Chr35.fasta
 156132220 Fern_H0_Chr36.fasta
 152473430 Fern_H0_Chr37.fasta
 152299297 Fern_H0_Chr38.fasta
 152259759 Fern_H0_Chr39.fasta
 149313362 Fern_H0_Chr40.fasta
 144249974 Fern_H0_Chr41.fasta
 138742306 Fern_H0_Chr42.fasta
 136342097 Fern_H0_Chr43.fasta
 116499002 Fern_H0_Chr44.fasta
7 094 536 786 total

#--------------
# 12 chromsomes and 2.1Gbp per chr, 25Gbp assembly

# divide chromosomes into 100Mbp segments?
cd /home/celphin/scratch/repeats/input_chromosomes/Fern/chromosome_files

# test
#split -l 1666666 Fern_H0_Chr06.fasta Fern_H0_Chr06_part --numeric-suffixes=1
#rm Fern_H0_Chr06.fasta

# loop
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

###########################
Chr ID
1	6
2	7
3	8
4	9
5	10
6	11
7	12
8	13
9	14
10	15
11	16
12	17
13	18
14	19
15	20
16	21
17	22
18	23
19	24
20	25
21	26
22	27
23	28
24	29
25	30
26	31
27	32
28	33
29	34
30	35
31	36
32	37
33	38
34	39
35	40
36	41
37	42
38	43
39	44



#############################################
# Maize
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02029-9
# abnormal chr 10 - meiotic drive - check this??

# 10 chromosomes, 2.6Gbp, N50 200Mbp

# https://www.ncbi.nlm.nih.gov/assembly/GCA_029775835.1

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Maize

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/775/835/GCA_029775835.1_ASM2977583v1/GCA_029775835.1_ASM2977583v1_genomic.fna.gz

# Unzip
gunzip GCA_029775835.1_ASM2977583v1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_029775835.1_ASM2977583v1_genomic.fna > Maize_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Maize_2.fasta

rename CP1167 Maize_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;
#########################
Chr ID
1	39
2	40
3	41
4	42
5	43
6	44
7	45
8	46
9	47
10	48
	
#######################################
# try Maize again with this genome:
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902167145.1/

# name it Corn not Maize to differentiate

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Corn

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz

# Unzip
gunzip GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna > Corn_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Corn_2.fasta

find -type f -size -5000000c -delete
rename  NC_050 Corn_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


##############################################
# Fruit fly: https://elifesciences.org/articles/49002
# https://www.ncbi.nlm.nih.gov/bioproject/PRJNA545704

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_008121275.1/ lowei
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_004329205.1/ pseudoobscura
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_008121225.1/ athabasca
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_008121235.1/ subobscura

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Fly

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/121/275/GCA_008121275.1_UCBerk_Dlow_1.0/GCA_008121275.1_UCBerk_Dlow_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/329/205/GCA_004329205.1_UCBerk_Dpse_1.0/GCA_004329205.1_UCBerk_Dpse_1.0_genomic.fna.gz 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/121/225/GCA_008121225.1_UCBerk_Dath_EA_1.0/GCA_008121225.1_UCBerk_Dath_EA_1.0_genomic.fna.gz 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/121/235/GCF_008121235.1_UCBerk_Dsub_1.0/GCF_008121235.1_UCBerk_Dsub_1.0_genomic.fna.gz

# Unzip
gunzip GCA_008121275.1_UCBerk_Dlow_1.0_genomic.fna.gz
gunzip GCA_004329205.1_UCBerk_Dpse_1.0_genomic.fna.gz 
gunzip GCA_008121225.1_UCBerk_Dath_EA_1.0_genomic.fna.gz 
gunzip GCF_008121235.1_UCBerk_Dsub_1.0_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_008121275.1_UCBerk_Dlow_1.0_genomic.fna > FlyDlow_2.fasta
seqkit seq -w 60 GCA_004329205.1_UCBerk_Dpse_1.0_genomic.fna > FlyDpse_2.fasta
seqkit seq -w 60 GCA_008121225.1_UCBerk_Dath_EA_1.0_genomic.fna > FlyDath_2.fasta
seqkit seq -w 60 GCF_008121235.1_UCBerk_Dsub_1.0_genomic.fna > FlyDsub_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../FlyDlow_2.fasta
find -type f -size -5000000c -delete
rename CM0177 Fly_H0_ChrDlow *
rename .1.fasta .fasta *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../FlyDpse_2.fasta
find -type f -size -5000000c -delete
rename CM0147 Fly_H0_ChrDpse *
rename .1.fasta .fasta *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../FlyDath_2.fasta
find -type f -size -5000000c -delete
rename CM0177 Fly_H0_ChrDath *
rename .1.fasta .fasta *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../FlyDsub_2.fasta
find -type f -size -5000000c -delete
rename NC_0485 Fly_H0_ChrDsub *
rename .1.fasta .fasta *


# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

#############################################
# Cotton: https://pubmed.ncbi.nlm.nih.gov/36971835/ # not clearly long arrays of satellites, only retrotransposon repeats
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_025698475.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Cotton

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/698/475/GCA_025698475.1_ASM2569847v1/GCA_025698475.1_ASM2569847v1_genomic.fna.gz

# Unzip
gunzip GCA_025698475.1_ASM2569847v1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_025698475.1_ASM2569847v1_genomic.fna > Cotton_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Cotton_2.fasta
find -type f -size -5000000c -delete
rename CM0469 Cotton_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


##############################################
# Brassica rapa
# https://onlinelibrary.wiley.com/doi/full/10.1111/pbi.14015 
# http://39.100.233.196:82/download_genome/Brassica_Genome_data/Brara_Chiifu_V4.0/ 

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Brassica

# go to FTP directory for GenBank assembly
wget http://39.100.233.196:82/download_genome/Brassica_Genome_data/Brara_Chiifu_V4.0/Brapa_genome_v4.0_chrom.fasta.gz

# Unzip
gunzip Brapa_genome_v4.0_chrom.fasta.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 Brapa_genome_v4.0_chrom.fasta > Brassica_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Brassica_2.fasta

find -type f -size -5000000c -delete
rename A Brassica_H0_Chr *


# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


###########################################
# ZigZag eel: https://www.ncbi.nlm.nih.gov/datasets/taxonomy/205130/ 
# https://pubmed.ncbi.nlm.nih.gov/34253240/#&gid=article-figures&pid=fig-2-uid-1

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Eel

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/455/535/GCA_019455535.1_Marm_hapX/GCA_019455535.1_Marm_hapX_genomic.fna.gz

# Unzip
gunzip GCA_019455535.1_Marm_hapX_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_019455535.1_Marm_hapX_genomic.fna > Eel_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Eel_2.fasta

find -type f -size -5000000c -delete
rename CM0335 Eel_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

##########################################
# Alligator
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_030867095.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Alligator

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/867/095/GCF_030867095.1_rAllMis1/GCF_030867095.1_rAllMis1_genomic.fna.gz

# Unzip
gunzip GCF_030867095.1_rAllMis1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCF_030867095.1_rAllMis1_genomic.fna > Alligator_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Alligator_2.fasta

find -type f -size -5000000c -delete
rename NC_0818 Alligator_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


##########################################
# Taxus (Yew)
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_018340775.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Yew

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/340/775/GCA_018340775.1_ASM1834077v1/GCA_018340775.1_ASM1834077v1_genomic.fna.gz

# Unzip
gunzip GCA_018340775.1_ASM1834077v1_genomic.fna.gz

# if need to fold 
salloc -c1 --time 2:50:00 --mem 192000M --account def-rieseber
module load seqkit/2.3.1
seqkit seq -w 60 GCA_018340775.1_ASM1834077v1_genomic.fna > Yew_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Yew_2.fasta

find -type f -size -5000000c -delete
rename CM0312 Yew_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

###################################
# Holocentric centromeres
# https://link.springer.com/article/10.1007/s10577-012-9292-1 

#     Chionographis japonica 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_947650365.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Chionographis

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/947/650/365/GCA_947650365.1_Cjaponica_Koshiki_v1.0/GCA_947650365.1_Cjaponica_Koshiki_v1.0_genomic.fna.gz

# Unzip
gunzip GCA_947650365.1_Cjaponica_Koshiki_v1.0_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_947650365.1_Cjaponica_Koshiki_v1.0_genomic.fna > Chionographis_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Chionographis_2.fasta

find -type f -size -10000000c -delete
rename OX3935 Chionographis_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


#--------------------------------------
# Luzula sylvatica
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_946800325.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Luzula

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/946/800/325/GCA_946800325.1_lpLuzSylv1.1/GCA_946800325.1_lpLuzSylv1.1_genomic.fna.gz

# Unzip
gunzip GCA_946800325.1_lpLuzSylv1.1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_946800325.1_lpLuzSylv1.1_genomic.fna > Luzula_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Luzula_2.fasta

find -type f -size -5000000c -delete
rename OX3269 Luzula_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;



#########################################
# Pineapple
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_029339315.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Pineapple

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/339/315/GCA_029339315.1_ASM2933931v1/GCA_029339315.1_ASM2933931v1_genomic.fna.gz

# Unzip
gunzip GCA_029339315.1_ASM2933931v1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_029339315.1_ASM2933931v1_genomic.fna > Pineapple_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Pineapple_2.fasta

find -type f -size -5000000c -delete
rename 	CM0555 Pineapple_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


#######################################
# Ginkgo 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_024626585.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Ginkgo

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/626/585/GCA_024626585.1_ASM2462658v1/GCA_024626585.1_ASM2462658v1_genomic.fna.gz

# Unzip
gunzip GCA_024626585.1_ASM2462658v1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_024626585.1_ASM2462658v1_genomic.fna > Ginkgo_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Ginkgo_2.fasta

find -type f -size -5000000c -delete
rename JANKJI0100000 Ginkgo_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

#######################################
# Polygonum 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_934048045.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Polygonum

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/934/048/045/GCA_934048045.1_dcPolAvic1.1/GCA_934048045.1_dcPolAvic1.1_genomic.fna.gz

# Unzip
gunzip GCA_934048045.1_dcPolAvic1.1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_934048045.1_dcPolAvic1.1_genomic.fna > Polygonum_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Polygonum_2.fasta

find -type f -size -5000000c -delete
rename OW20402 Polygonum_H0_Chr *
rename OW20403 Polygonum_H0_Chr1 *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

#######################################
# Fagopyrum 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_002319775.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Fagopyrum

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/319/775/GCA_002319775.1_Ft1.0/GCA_002319775.1_Ft1.0_genomic.fna.gz

# Unzip
gunzip GCA_002319775.1_Ft1.0_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_002319775.1_Ft1.0_genomic.fna > Fagopyrum_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Fagopyrum_2.fasta

find -type f -size -5000000c -delete
rename CM0082 Fagopyrum_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;
#############################################
# Fabaceae - pea
# https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010633
# https://www.nature.com/articles/s41588-019-0480-1
# chr6 centromere

# 7 chromosomes, 3Gbp

# https://www.ncbi.nlm.nih.gov/assembly/GCF_024323335.1

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Pea

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/323/335/GCA_024323335.2_CAAS_Psat_ZW6_1.0/GCA_024323335.2_CAAS_Psat_ZW6_1.0_genomic.fna.gz

# Unzip
gunzip GCA_024323335.2_CAAS_Psat_ZW6_1.0_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_024323335.2_CAAS_Psat_ZW6_1.0_genomic.fna > Pea_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Pea_2.fasta
find -type f -size -8000000c -delete
rename CM0443 Pea_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

#################################
Chr ID
1	45
2	46
3	47
4	48
5	49
6	50
7	51


#############################################
# Human
# https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.3/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Human

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.3_T2T-CHM13v1.1/GCA_009914755.3_T2T-CHM13v1.1_genomic.fna.gz

# Unzip
gunzip GCA_009914755.3_T2T-CHM13v1.1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_009914755.3_T2T-CHM13v1.1_genomic.fna > Human_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Human_2.fasta

rename CP0682 Human_H0_Chr *
rename .2.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

#---------------------------------
#Y chromosome: 
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/020/881/995/GCA_020881995.2_ASM2088199v2/

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/020/881/995/GCA_020881995.2_ASM2088199v2/GCA_020881995.2_ASM2088199v2_genomic.fna.gz

gunzip GCA_020881995.2_ASM2088199v2_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_020881995.2_ASM2088199v2_genomic.fna > Human_3.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Human_3.fasta

rename CP0682 Human_H0_Chr *
rename .2.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


#######################
Chr	ID
1	77
2	76
3	75
4	74
5	73
6	72
7	71
8	70 # past centromere?
9	69
10	68
11	67
12	66
13	65
14	64
15	63
16	62
17	61
18	60
19	59
20	58
21	57
22	56
X	55
MT	54


############################################
# Mouse
# Older: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27

# Try: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030265425.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Mouse

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/265/425/GCA_030265425.1_NEI_Mmus_1.0/GCA_030265425.1_NEI_Mmus_1.0_genomic.fna.gz

# Unzip
gunzip GCA_030265425.1_NEI_Mmus_1.0_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_030265425.1_NEI_Mmus_1.0_genomic.fna > Mouse_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Mouse_2.fasta
find -type f -size -10000000c -delete

rename CM0582 Mouse_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

##################################
Chr	ID
1   41
2   42
3 43
4 44
5 45
6
7
8
9 49
10 50
11 51
12 52
13 53
16 56
X   60
Y   61 # too short not included

############################################
# Tomato

# Older
# Solanum galapagense 
# Solanum corneliomulleri
# https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA809001
# https://www.nature.com/articles/s41588-023-01340-y
# pangenome - inversions on chromosome 3
# "clade IV-specific inversion on chromosome 3 from 47.5 Mb to 54.6 Mb is shown, as evidenced by abnormally strong interactions around the inversion breakpoints"
# also check chr 6 and 7 for inversions
#--------------------------
# download
# cd /home/celphin/scratch/repeats/input_chromosomes/Tomato

# # go to FTP directory for GenBank assembly
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/704/825/GCA_027704825.1_ASM2770482v1/GCA_027704825.1_ASM2770482v1_genomic.fna.gz #Sgal
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/704/805/GCA_027704805.1_ASM2770480v1/GCA_027704805.1_ASM2770480v1_genomic.fna.gz #Scor

# # Unzip
# gunzip GCA_027704825.1_ASM2770482v1_genomic.fna.gz
# gunzip GCA_027704805.1_ASM2770480v1_genomic.fna.gz 

# # if need to fold 
# module load seqkit/2.3.1
# seqkit seq -w 60 GCA_027704825.1_ASM2770482v1_genomic.fna  > Sgal_2.fasta
# seqkit seq -w 60 GCA_027704805.1_ASM2770480v1_genomic.fna > Scor_2.fasta

# mkdir chromosome_files/
# cd ./chromosome_files/

# awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Sgal_2.fasta
# find -type f -size -10000000c -delete
# rename JAKWIZ0100000 Sgal_H0_Chr *
# rename .1.fasta .fasta *


# awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Scor_2.fasta
# find -type f -size -10000000c -delete
# rename JAKWIY0100030 Scor_H0_Chr *
# rename .1.fasta .fasta *

# rename Scor_H0_Chr Tomato_H0_ChrScor *
# rename Sgal_H0_Chr Tomato_H0_ChrSgal *


#############################################
# Rice (Oryza sativa L.)
# https://www.nature.com/articles/s41559-022-01974-x 
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02861-9 

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001433935.1/
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000231095.2/

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2034625/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Rice

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/231/095/GCF_000231095.2_ObraRS2/GCF_000231095.2_ObraRS2_genomic.fna.gz

# Unzip
gunzip GCF_000231095.2_ObraRS2_genomic.fna.gz 
gunzip GCF_001433935.1_IRGSP-1.0_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCF_001433935.1_IRGSP-1.0_genomic.fna  > Rice_2Sat.fasta
seqkit seq -w 60 GCF_000231095.2_ObraRS2_genomic.fna  > Rice_2Bra.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Rice_2Sat.fasta
find -type f -size -10000000c -delete
rename NC_0292 Rice_H0_ChrSat *
rename .1.fasta .fasta *


awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Rice_2Bra.fasta
find -type f -size -10000000c -delete
rename 	NC_0231 Rice_H0_ChrBra *
rename .1.fasta .fasta *



#############################################
# Chicken
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3913482/ 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027557775.1/


# download
cd /home/celphin/scratch/repeats/input_chromosomes/Chicken

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/557/775/GCA_027557775.1_bGalGal4.pri/GCA_027557775.1_bGalGal4.pri_genomic.fna.gz

# Unzip
gunzip GCA_027557775.1_bGalGal4.pri_genomic.fna.gz 

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_027557775.1_bGalGal4.pri_genomic.fna  > Chicken_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Chicken_2.fasta
find -type f -size -10000000c -delete
rename CM0503 Chicken_H0_Chr *
rename .1.fasta .fasta *


# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


#############################################
# Wheat 
# wheat cultivar Aikang58 (AK58)
# https://pubmed.ncbi.nlm.nih.gov/36739481/

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018294505.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Wheat

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/294/505/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna.gz

# Unzip
gunzip GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna  > Wheat_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Wheat_2.fasta
find -type f -size -10000000c -delete
rename NC_057 Wheat_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


#########################
# try another Tomato species - with less gaps
# https://www.ncbi.nlm.nih.gov/genome/?term=txid4081[Organism:exp]

# main reference - Solanum lycopersicum 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000188115.5/

#     Solanum pimpinellifolium (currant tomato)
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_014964335.1/

#     Solanum pennellii 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001406875.1/

#    Solanum arcanum 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027789145.1/

#     Solanum lycopersicoides 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022817965.1/

#----------------------------------------

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Tomatov2

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.5_SL3.1/GCF_000188115.5_SL3.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/964/335/GCA_014964335.1_ASM1496433v1/GCA_014964335.1_ASM1496433v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/406/875/GCF_001406875.1_SPENNV200/GCF_001406875.1_SPENNV200_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/789/145/GCA_027789145.1_ASM2778914v1/GCA_027789145.1_ASM2778914v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/817/965/GCA_022817965.1_SlydLA2951_v2.0/GCA_022817965.1_SlydLA2951_v2.0_genomic.fna.gz

# Unzip
gunzip GCF_000188115.5_SL3.1_genomic.fna.gz
gunzip GCA_014964335.1_ASM1496433v1_genomic.fna.gz 
gunzip GCF_001406875.1_SPENNV200_genomic.fna.gz
gunzip GCA_027789145.1_ASM2778914v1_genomic.fna.gz
gunzip GCA_022817965.1_SlydLA2951_v2.0_genomic.fna.gz


# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCF_000188115.5_SL3.1_genomic.fna  > Sly_2.fasta
seqkit seq -w 60 GCA_014964335.1_ASM1496433v1_genomic.fna > Spim_2.fasta
seqkit seq -w 60 GCF_001406875.1_SPENNV200_genomic.fna > Spen_2.fasta
seqkit seq -w 60 GCA_027789145.1_ASM2778914v1_genomic.fna > Sarc_2.fasta
seqkit seq -w 60 GCA_022817965.1_SlydLA2951_v2.0_genomic.fna > Slys_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Sly_2.fasta
find -type f -size -10000000c -delete
rename 	NC_0154 Tomato_H0_ChrSly *
rename .3.fasta .fasta *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Spim_2.fasta
find -type f -size -10000000c -delete
rename 	CM0265 Tomato_H0_ChrSpim *
rename .1.fasta .fasta *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Spen_2.fasta
find -type f -size -10000000c -delete
rename 	NC_0286 Tomato_H0_ChrSpen *
rename .1.fasta .fasta *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Sarc_2.fasta
find -type f -size -10000000c -delete
rename CM0509 Tomato_H0_ChrSarc *
rename .1.fasta .fasta *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Slys_2.fasta
find -type f -size -10000000c -delete
rename 	CM0406 Tomato_H0_ChrSlys *
rename .1.fasta .fasta *


# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


#############################
Sly
Chr ID
38-49
1:12

38 1
39 2
40 3
41 4

#############################
Spim
Chr ID


#############################
Spen
Chr ID


#############################
Sarc
Chr ID


#############################
Slys
Chr ID


#######################################
# Zebrafish
# https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0195-y
# https://en.wikipedia.org/wiki/Compositional_domain
# TEs are massively represented in large and GC-poor genomes such zebrafish 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_944039275.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Zebrafish

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/944/039/275/GCA_944039275.1_fDanRer4.1/GCA_944039275.1_fDanRer4.1_genomic.fna.gz

# Unzip
gunzip GCA_944039275.1_fDanRer4.1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  GCA_944039275.1_fDanRer4.1_genomic.fna > Zebrafish_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Zebrafish_2.fasta
find -type f -size -10000000c -delete
rename OX0633 Zebrafish_H0_Chr *
rename .1.fasta .fasta *



#####################################
# Pufferfish
# TEs are clearly depleted in compact and GC-rich genomes

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Pufferfish

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/901/000/725/GCA_901000725.3_fTakRub1.3/GCA_901000725.3_fTakRub1.3_genomic.fna.gz

# Unzip
gunzip GCA_901000725.3_fTakRub1.3_genomic.fna.gz


# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_901000725.3_fTakRub1.3_genomic.fna  > Pufferfish_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Pufferfish_2.fasta
find -type f -size -10000000c -delete
rename LR5842 Pufferfish_H0_Chr *
rename .1.fasta .fasta *
rename .2.fasta .fasta *


##########################################
# Anolis sagrei
# https://www.nature.com/articles/s42003-022-04074-5
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/TTKBFU

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Anolis

# copied from website link above

# Unzip
gunzip AnoSag2.1.fa.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  AnoSag2.1.fa > Anolis_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Anolis_2.fasta
find -type f -size -10000000c -delete
rename scaffold Anolis_H0_Chr *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


#############################################
# Green Algae
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_023205875.1/
#     Chloropicon primus 

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Algae

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/205/875/GCA_023205875.1_ASM2320587v1/GCA_023205875.1_ASM2320587v1_genomic.fna.gz

# Unzip
gunzip GCA_023205875.1_ASM2320587v1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60 GCA_023205875.1_ASM2320587v1_genomic.fna > Algae_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Algae_2.fasta


rename CP0607 Algae_H0_Chr *
rename .1.fasta .fasta *


##########################################
# Ant
# Monomorium pharaonis
# https://academic.oup.com/gigascience/article/9/12/giaa143/6034789
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013373865.1/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Ant

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/373/865/GCF_013373865.1_ASM1337386v2/GCF_013373865.1_ASM1337386v2_genomic.fna.gz

# Unzip
gunzip GCF_013373865.1_ASM1337386v2_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  GCF_013373865.1_ASM1337386v2_genomic.fna  > Ant_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Ant_2.fasta
find -type f -size -10000000c -delete
rename  NC_0504 Ant_H0_Chr *
rename .1.fasta .fasta *

############################################
# Zebra finch
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003957565.2/


# download
cd /home/celphin/scratch/repeats/input_chromosomes/Finch

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz

# Unzip
gunzip GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  GCF_003957565.2_bTaeGut1.4.pri_genomic.fna  > Finch_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Finch_2.fasta
find -type f -size -5000000c -delete
rename   NC_0442 Finch_H0_Chr *
rename .2.fasta .fasta *

###############################
# Duck
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_008746955.2/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Duck

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/746/955/GCA_008746955.2_CAU-Wild1.1/GCA_008746955.2_CAU-Wild1.1_genomic.fna.gz

# Unzip
gunzip GCA_008746955.2_CAU-Wild1.1_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  GCA_008746955.2_CAU-Wild1.1_genomic.fna  > Duck_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Duck_2.fasta
find -type f -size -5000000c -delete
rename   CM0183 Duck_H0_Chr *
rename .1.fasta .fasta *

# needed if over 100Mbp
# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

####################################
# Ornithorhynchus anatinus (Platypus) 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004115215.2/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Ha412

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  Ha412HOv2.0-20181130.fasta  > Ha412_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Ha412_2.fasta
find -type f -size -5000000c -delete
rename Ha412HOChr Ha412_H0_Chr *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

####################################
# Sunflower

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Sunflower
ls
HA89_2.fasta                    
HanXRQv2.fasta                      
Hneg2702r1.0-20210702.genome.fasta
Hanom1.0-20201016.genome.fasta  
Harg2202r1.0-20210824.genome.fasta
HanPSC8.fasta                   
Hdeb2414r1.0-20210622.genome.fasta

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  HanHA89r1.0-20210811.genome.fasta    > HA89_2.fasta
seqkit seq -w 60  HanXRQv2.fasta                       > HanXRQ_2.fasta
seqkit seq -w 60  Hneg2702r1.0-20210702.genome.fasta   > Hneg_2.fasta
seqkit seq -w 60  Hanom1.0-20201016.genome.fasta       > Hanom_2.fasta
seqkit seq -w 60  Harg2202r1.0-20210824.genome.fasta   > Harg_2.fasta
seqkit seq -w 60  HanPSC8.fasta                        > PSC8_2.fasta
seqkit seq -w 60  Hdeb2414r1.0-20210622.genome.fasta   > Hdeb_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../HA89_2.fasta
find -type f -size -5000000c -delete
rename  HanHA89Chr Sunflower_H0_ChrHa89 *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../HanXRQ_2.fasta
find -type f -size -5000000c -delete
rename  HanXRQChr Sunflower_H0_ChrXRQ *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../PSC8_2.fasta
find -type f -size -5000000c -delete
rename  HanPSC8Chr Sunflower_H0_ChrPSC8 *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


####################################
# Ornithorhynchus anatinus (Platypus) 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004115215.2/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Platypus

# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/115/215/GCF_004115215.2_mOrnAna1.pri.v4/GCF_004115215.2_mOrnAna1.pri.v4_genomic.fna.gz

# Unzip
gunzip GCF_004115215.2_mOrnAna1.pri.v4_genomic.fna.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  GCF_004115215.2_mOrnAna1.pri.v4_genomic.fna  > Platypus_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Platypus_2.fasta
find -type f -size -5000000c -delete
rename   NC_0417 Platypus_H0_Chr *
rename .1.fasta .fasta *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

####################################
# Strawberry
# http://eplant.njau.edu.cn/strawberry/
# http://eplantftp.njau.edu.cn/Fragaria/F._vesca/F._vesca_v6.0/

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Strawberry

# go to FTP directory for GenBank assembly
wget http://eplantftp.njau.edu.cn/Fragaria/F._vesca/F._vesca_v6.0/Fragaria_vesca_v6_genome.fasta

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  Fragaria_vesca_v6_genome.fasta  > Strawberry_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Strawberry_2.fasta
find -type f -size -5000000c -delete
rename   chr Strawberry_H0_Chr *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


####################################
# Einkorn Wheat
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.v41ns1rxj
# https://www.nature.com/articles/s41586-023-06389-7

# download
cd /home/celphin/scratch/repeats/input_chromosomes/Einkorn

# go to FTP directory for GenBank assembly
wget https://datadryad.org/stash/downloads/file_stream/2350888  # Triticum monococcum monococcum
wget https://datadryad.org/stash/downloads/file_stream/2350887  # Triticum monococcum aegilpoides

# rename 
mv 2350888 monococcum.fasta.gz
mv 2350887 aegilpoides.fasta.gz

# Unzip
gunzip monococcum.fasta.gz
gunzip aegilpoides.fasta.gz

# if need to fold 
module load seqkit/2.3.1
seqkit seq -w 60  monococcum.fasta  > Einkorn_monococcum_2.fasta
seqkit seq -w 60  aegilpoides.fasta  > Einkorn_aegilpoides_2.fasta

mkdir chromosome_files/
cd ./chromosome_files/

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Einkorn_monococcum_2.fasta
find -type f -size -5000000c -delete
rename   chr Einkorn_H0_Chrmono *_TA10622.fasta
rename _TA10622. . *

awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Einkorn_aegilpoides_2.fasta
find -type f -size -5000000c -delete
rename   chr Einkorn_H0_Chraeg *_TA299.fasta
rename _TA299. . *

# loop to split long chromosomes
ls > ../list_chromosomes.txt
sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

while IFS= read -r Chr; 
do
split -l 1666666 ${Chr}.fasta ${Chr}part --numeric-suffixes=1
rm ${Chr}.fasta
done < ../list_chrnum.txt

wc -c *

cd ..
find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;


#############################################
# Reindeer
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-023-09189-5 


#############################################
# kiwifruit Actinidia chinensis 
# https://academic.oup.com/hr/article/10/2/uhac264/6865344 


####################################
# Banana: 
# https://www.nature.com/articles/s42003-021-02559-3 


###################################
# Watermelon 
# https://www.sciencedirect.com/science/article/pii/S1674205222001927 


##########################
# Fungi 
# https://academic.oup.com/dnaresearch/article/30/3/dsad006/7142856


###################################################
# code to run on chromosome files by making the R and bash  scripts separately 
# https://www.taniarascia.com/how-to-create-and-use-bash-scripts/

nano /home/celphin/scratch/repeats/pre-repeats.sh

#!/bin/bash

# input directory set to just above chromosome_files directory
# inlude / at the end of pathname
pathname=$1
SPP=$2

cd ${pathname}
echo ${pathname}
echo $(pwd)

# get chromosome lengths
wc -c ./chromosome_files/* > chromosome_lengths.txt

# make plots of CG isochores
# list of chromosomes
ls ./chromosome_files/* > chr_list.txt

while IFS= read -r chr; 
do
isochore -sequence ${chr} -outfile ${chr}.isochore  -window 10000 -shift 5000 -graph png -goutfile ${chr}_png
done < chr_list.txt

mkdir isochore
mv ./chromosome_files/*.isochore isochore
mv ./chromosome_files/*.png isochore

# identify telomere locations - plants
grep -n "TTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTT" ./chromosome_files/* > plant_telomere_positions.txt

# identify telomere locations - animals
grep -n "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT" ./chromosome_files/* > animal_telomere_positions.txt

qq='"'

# make R script with correct input
cat << EOF1 > repeats_fourier.R

library(RepeatOBserverV1)

inpath=${qq}${pathname}/chromosome_files/${qq}
fname=${qq}${SPP}${qq}

#----------------------------------------

outpath="~/scratch/repeats/output_chromosomes"
x_cpu=19
pflag=FALSE
writeflag=FALSE
plotflag=FALSE

nam_list0 <- list.files(inpath)
nam_list1 <- tools::file_path_sans_ext(nam_list0)
nam_list1
nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
chr_list <- nam_list3[,1]
chr_list

#--------------------
# run intial spectra if directories do not exist
for (nam in nam_list1){
 run_plot_chromosome_parallel(nam=nam, fname=fname, inpath=inpath, outpath=outpath, pflag=pflag, plotflag=plotflag,  writeflag=writeflag, x_cpu=x_cpu)
}

#-----------------------
# try to run writing of summary files on different cpu

nam_write_all_spec <- function(x, nam_list1=nam_list1, chr_list=chr_list, fname=fname, inpath=inpath, outpath=outpath, pflag=pflag, plotflag=plotflag,  writeflag=writeflag, x_cpu=x_cpu){
  library(RepeatOBserverV1)
  nam <<- nam_list1[x]
  print(nam)
  chromosome <- chr_list[x]
  print(chromosome)
  write_All_spec_DNAwalk(nam=nam, fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
}

x_cpu=2

cl <- parallel::makeCluster(x_cpu)
results1 <- parallel::parSapply(cl, base::seq_along(c(1:length(nam_list1))), nam_write_all_spec, nam_list1=nam_list1, chr_list=chr_list, fname=fname, inpath=inpath, outpath=outpath, pflag=pflag, plotflag=plotflag,  writeflag=writeflag, x_cpu=x_cpu)

#-----------------------

# run chromosomes on different cpu
uni_chr_list <- unique(chr_list)

chromosome_summary <-  function(x, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath){
  library(RepeatOBserverV1)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  merge_spectra(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
  join_chromosome_parts(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
  run_summary_hist(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
  run_diversity_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
  run_summary_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
}
x_cpu=5
cl <- parallel::makeCluster(x_cpu)
results2 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), chromosome_summary,uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)

EOF1


# make R script with correct input
cat << EOF2 > centromere_histogram.R

library(RepeatOBserverV1)

inpath=${qq}${pathname}/chromosome_files/${qq}
fname=${qq}${SPP}${qq}

#----------------------------------------

outpath="~/scratch/repeats/output_chromosomes"
x_cpu=19
pflag=FALSE
writeflag=FALSE
plotflag=FALSE

nam_list0 <- list.files(inpath)
nam_list1 <- tools::file_path_sans_ext(nam_list0)
nam_list1
nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
chr_list <- nam_list3[,1]
chr_list


# run chromosomes on different cpu
uni_chr_list <- unique(chr_list)

chromosome_summary <-  function(x, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath){
  library(RepeatOBserverV1)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  #run_summary_hist(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
  run_summary_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
#  run_diversity_plots_no_telomere(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
#  run_diversity_plots_35_2000_no_telomere(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
#  run_diversity_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
#  run_diversity_plots_35_2000(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
}

cl <- parallel::makeCluster(x_cpu)
results2 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), chromosome_summary,uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)


EOF2


EOF

chmod +x /home/celphin/scratch/repeats/pre-repeats.sh


#------------------------------

nano /home/celphin/scratch/repeats/post_repeats.sh

#!/bin/bash

SPP_Hap=$1

cd /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/
mkdir Summary_output; cd Summary_output

cd /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Summary_output

#-----------

mkdir spectra_parts_35-2000
cp -v -u  /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp35_2000seq2501_*TRUE.png ./spectra_parts_35-2000

mkdir spectra_parts_15-35
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp15_35seq2501_*TRUE.png ./spectra_parts_15-35

mkdir spectra_parts_2-8
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp2_8seq2501_*TRUE.png ./spectra_parts_2-8

#------------------------------
mkdir spectra_total_merged
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/*All_spec_bp35_2000*.png ./spectra_total_merged

mkdir histograms
cp -v -u  /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/histograms/*POWER_SUM*s_0.5std_1*.pdf ./histograms/

mkdir DNAwalks
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/*DNAwalk*.png ./DNAwalks

#----------

mkdir Shannon_div
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/${SPP_Hap}_Chr*_Shannon_plot_norm.png ./Shannon_div


mkdir Shannon_div_100
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/*roll_mean_Shannon_100.png ./Shannon_div_100

mkdir Shannon_div_1000
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/*roll_mean_Shannon_1000.png ./Shannon_div_1000

mkdir Shannon_div_250
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/*roll_mean_Shannon_250.png ./Shannon_div_250

mkdir Shannon_div_500
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/*roll_mean_Shannon_500.png ./Shannon_div_500

#-------------

rm Centromere_summary.txt
cat /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/histograms/Centromere_*.txt > Centromere_summary.txt
#more Centromere_summary.txt

rm Centromere_summary_Shannon.txt
cat /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/*Centromere_MIN_Shannon.txt > Centromere_summary_Shannon.txt

grep "cent25 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_25.txt
grep "cent100 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_100.txt
grep "cent250 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_250.txt
grep "cent500 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_500.txt
grep "cent1000 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_1000.txt
grep "centwind " Centromere_summary_Shannon_35_no_telo.txt > Centromere_summary_Shannon_wind_35_no_telo.txt

mkdir Shannon_div_window
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/*Shannon_div_window*.png ./Shannon_div_window

EOF

chmod +x /home/celphin/scratch/repeats/post_repeats.sh


#------------------------
# slurm scripts to launch runs using the R scripts made in the bash scripts above

cd /home/celphin/scratch/repeats/scripts

cat << EOF > Arabidopsis_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Arabidopsis" "Arab_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Arabidopsis 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Arabidopsis/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Arabidopsis/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Arab_H0"

EOF

sbatch Arabidopsis_cedar_repeats.sh

#---------------------------------
# Arab_COLCEN
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Arab_COLCEN_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Arab_COLCEN" "COLCEN_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Arab_COLCEN 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Arab_COLCEN/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Arab_COLCEN/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "COLCEN_H0"

EOF

sbatch Arab_COLCEN_cedar_repeats.sh


#---------------------------------
# Orchid
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Orchid_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Orchid" "Orchid_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Orchid 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Orchid/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Orchid/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Orchid_H0"

EOF

sbatch Orchid_cedar_repeats.sh


#---------------------------------
# Maize
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Maize_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Maize" "Maize_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Maize 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Maize/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Maize/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Maize_H0"

EOF

sbatch Maize_cedar_repeats.sh


#---------------------------------
# Pea
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Pea_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Pea" "Pea_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Pea 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Pea/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Pea/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Pea_H0"

EOF

sbatch Pea_cedar_repeats.sh

#---------------------------------
# Mouse
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Mouse_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Mouse" "Mouse_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Mouse 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Mouse/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Mouse/centromere_histogram.R


srun /home/celphin/scratch/repeats/post_repeats.sh "Mouse_H0"

EOF

sbatch Mouse_cedar_repeats.sh


#---------------------------------
# Fern
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Fern_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Fern" "Fern_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Fern 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Fern/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Fern/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Fern_H0"

EOF

sbatch Fern_cedar_repeats.sh

#--------------------------------
# Human
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Human_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Human" "Human_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Human 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Human/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Human/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Human_H0"

EOF

sbatch Human_cedar_repeats.sh

#---------------------------------
# Rice
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Rice_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Rice" "Rice_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Rice 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Rice/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Rice/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Rice_H0"

EOF

sbatch Rice_cedar_repeats.sh

#---------------------------------
# Chicken
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Chicken_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Chicken" "Chicken_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Chicken 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Chicken/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Chicken/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Chicken_H0"

EOF

sbatch Chicken_cedar_repeats.sh

#---------------------------------
# Wheat
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Wheat_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Wheat" "Wheat_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Wheat 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Wheat/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Wheat/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Wheat_H0"

EOF

sbatch Wheat_cedar_repeats.sh


#######################################
# Zebrafish
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Zebrafish_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Zebrafish" "Zebrafish_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Zebrafish 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Zebrafish/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Zebrafish/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Zebrafish_H0"

EOF

sbatch Zebrafish_cedar_repeats.sh

#---------------------------------
# Pufferfish
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Pufferfish_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Pufferfish" "Pufferfish_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Pufferfish 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Pufferfish/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Pufferfish/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Pufferfish_H0"

EOF

sbatch Pufferfish_cedar_repeats.sh

#---------------------------------
# Anolis
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Anolis_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Anolis" "Anolis_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Anolis 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Anolis/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Anolis/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Anolis_H0"

EOF

sbatch Anolis_cedar_repeats.sh


#---------------------------------
# Algae
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Algae_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Algae" "Algae_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Algae 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Algae/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Algae/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Algae_H0"

EOF

sbatch Algae_cedar_repeats.sh


#---------------------------------
# Ant
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Ant_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Ant" "Ant_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Ant 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Ant/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Ant/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Ant_H0"

EOF

sbatch Ant_cedar_repeats.sh

#---------------------------------
# Finch
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Finch_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Finch" "Finch_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Finch
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Finch/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Finch/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Finch_H0"

EOF

sbatch Finch_cedar_repeats.sh

#---------------------------------
# Duck
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Duck_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Duck" "Duck_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Duck
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Duck/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Duck/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Duck_H0"

EOF

sbatch Duck_cedar_repeats.sh

#---------------------------------
# Platypus
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Platypus_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Platypus" "Platypus_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Platypus
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Platypus/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Platypus/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Platypus_H0"

EOF

sbatch Platypus_cedar_repeats.sh

#---------------------------------
# Strawberry
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Strawberry_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Strawberry" "Strawberry_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Strawberry
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Strawberry/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Strawberry/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Strawberry_H0"

EOF

sbatch Strawberry_cedar_repeats.sh


#---------------------------------
# Einkorn
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Einkorn_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Einkorn" "Einkorn_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Einkorn
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Einkorn/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Einkorn/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Einkorn_H0"

EOF

sbatch Einkorn_cedar_repeats.sh


#---------------------------------
# Tomato
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Tomato_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Tomato" "Tomato_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Tomato 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Tomato/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Tomato/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Tomato_H0"

EOF

sbatch Tomato_cedar_repeats.sh

#####################################
# Corn
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Corn_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Corn" "Corn_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Corn 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Corn/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Corn/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Corn_H0"

EOF

sbatch Corn_cedar_repeats.sh

#---------------------------------
# Fly
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Fly_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Fly" "Fly_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Fly 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Fly/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Fly/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Fly_H0"

EOF

sbatch Fly_cedar_repeats.sh

#---------------------------------
# Cotton
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Cotton_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Cotton" "Cotton_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Cotton 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Cotton/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Cotton/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Cotton_H0"

EOF

sbatch Fly_cedar_repeats.sh
#---------------------------------
# Brassica
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Brassica_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Brassica" "Brassica_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Brassica 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Brassica/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Brassica/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Brassica_H0"

EOF

sbatch Brassica_cedar_repeats.sh


#---------------------------------
# Eel
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Eel_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Eel" "Eel_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Eel 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Eel/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Eel/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Eel_H0"

EOF

sbatch Eel_cedar_repeats.sh

#---------------------------------
# Ginkgo
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Ginkgo_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Ginkgo" "Ginkgo_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Ginkgo 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Ginkgo/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Ginkgo/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Ginkgo_H0"

EOF

sbatch Ginkgo_cedar_repeats.sh

#---------------------------------
# Polygonum
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Polygonum_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Polygonum" "Polygonum_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Polygonum 
module load StdEnv/2020 r/4.1.2
Rscript /home/celphin/scratch/repeats/input_chromosomes/Polygonum/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Polygonum/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Polygonum_H0"

EOF

sbatch Polygonum_cedar_repeats.sh

#---------------------------------
# Fagopyrum
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Fagopyrum_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Fagopyrum" "Fagopyrum_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Fagopyrum 
module load StdEnv/2020 r/4.1.2
Rscript /home/celphin/scratch/repeats/input_chromosomes/Fagopyrum/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Fagopyrum/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Fagopyrum_H0"

EOF

sbatch Fagopyrum_cedar_repeats.sh

#---------------------------------
# Alligator
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Alligator_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Alligator" "Alligator_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Alligator 
module load StdEnv/2020 r/4.1.2
Rscript /home/celphin/scratch/repeats/input_chromosomes/Alligator/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Alligator/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Alligator_H0"

EOF

sbatch Alligator_cedar_repeats.sh

#---------------------------------
# Yew
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Yew_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Yew" "Yew_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Yew 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Yew/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Yew/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Yew_H0"

EOF

sbatch Yew_cedar_repeats.sh
#---------------------------------
# Pineapple
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Pineapple_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Pineapple" "Pineapple_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Pineapple 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Pineapple/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Pineapple/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Pineapple_H0"

EOF

sbatch Pineapple_cedar_repeats.sh

#---------------------------------
# Luzula
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Luzula_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Luzula" "Luzula_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Luzula 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Luzula/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Luzula/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Luzula_H0"

EOF

sbatch Luzula_cedar_repeats.sh

#---------------------------------
# Chionographis
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Chionographis_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Chionographis" "Chionographis_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Chionographis 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Chionographis/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Chionographis/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Chionographis_H0"

EOF

sbatch Chionographis_cedar_repeats.sh

#---------------------------------
# Ha412
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Ha412_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Ha412" "Ha412_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Ha412 
module load StdEnv/2020 r/4.1.2
#Rscript /home/celphin/scratch/repeats/input_chromosomes/Ha412/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Ha412/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Ha412_H0"

EOF

sbatch Ha412_cedar_repeats.sh

#----------------------------------------
# Sunflower
cd /home/celphin/scratch/repeats/scripts

cat << EOF > Sunflower_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128000M

module load StdEnv/2020
module load emboss/6.6.0
srun /home/celphin/scratch/repeats/pre-repeats.sh "/home/celphin/scratch/repeats/input_chromosomes/Sunflower" "Sunflower_H0"

cd /home/celphin/scratch/repeats/input_chromosomes/Sunflower 
module load StdEnv/2020 r/4.1.2
Rscript /home/celphin/scratch/repeats/input_chromosomes/Sunflower/repeats_fourier.R
Rscript /home/celphin/scratch/repeats/input_chromosomes/Sunflower/centromere_histogram.R

srun /home/celphin/scratch/repeats/post_repeats.sh "Sunflower_H0"

EOF

sbatch Sunflower_cedar_repeats.sh

#---------------------------------

#####################
# to edit all to suite Graham nodes
# https://docs.alliancecan.ca/wiki/Graham#Node_characteristics 

cd /home/celphin/scratch/repeats/scripts

sed -i 's/#SBATCH --cpus-per-task=20/#SBATCH --cpus-per-task=15/g' *
# change the x_cpu to 10 and 14 for the various processes

sed -i 's/#SBATCH --mem=192000M/#SBATCH --mem=128000M/g' *.sh

##################################
# submit all

cd /home/celphin/scratch/repeats/scripts

sbatch Arabidopsis_cedar_repeats.sh
sbatch Arab_COLCEN_cedar_repeats.sh
sbatch Orchid_cedar_repeats.sh
sbatch Maize_cedar_repeats.sh
sbatch Pea_cedar_repeats.sh
sbatch Mouse_cedar_repeats.sh
sbatch Fern_cedar_repeats.sh
sbatch Human_cedar_repeats.sh
sbatch Rice_cedar_repeats.sh
sbatch Chicken_cedar_repeats.sh
sbatch Zebrafish_cedar_repeats.sh
sbatch Pufferfish_cedar_repeats.sh
sbatch Anolis_cedar_repeats.sh
sbatch Ant_cedar_repeats.sh
sbatch Finch_cedar_repeats.sh
sbatch Duck_cedar_repeats.sh
sbatch Platypus_cedar_repeats.sh
sbatch Strawberry_cedar_repeats.sh
sbatch Tomato_cedar_repeats.sh
sbatch Einkorn_cedar_repeats.sh
sbatch Corn_cedar_repeats.sh
sbatch Fly_cedar_repeats.sh
sbatch Brassica_cedar_repeats.sh
sbatch Eel_cedar_repeats.sh
sbatch Ginkgo_cedar_repeats.sh
sbatch Yew_cedar_repeats.sh
sbatch Pineapple_cedar_repeats.sh
sbatch Luzula_cedar_repeats.sh
sbatch Chionographis_cedar_repeats.sh
sbatch Polygonum_cedar_repeats.sh
sbatch Fagopyrum_cedar_repeats.sh
sbatch Alligator_cedar_repeats.sh
sbatch Ha412_cedar_repeats.sh
sbatch Sunflower_cedar_repeats.sh


######################################
# Make Summary folder for all chromosomes

cd ~/scratch/repeats
mkdir Summary

cd ~/scratch/repeats/Summary

mkdir spectra_total_merged
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*_All_spec_bp35_2000*.png ./spectra_total_merged
cd ~/scratch/repeats/Summary/spectra_total_merged
mkdir norm_diff
mv *diff_norm.png ./norm_diff


mkdir histograms
cp -v -u  /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/histograms/*POWER_SUM*s_0.5std_1*.pdf ./histograms/

mkdir DNAwalks
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*DNAwalk*.png ./DNAwalks

#-----------
cd ~/scratch/repeats/Summary

mkdir spectra_parts_35-2000
cp -v -u  /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/largeimages.png/All_spec1_*_Chr*_bp35_2000seq2501_*TRUE.png ./spectra_parts_35-2000

mkdir spectra_parts_15-35
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/largeimages.png/All_spec1_*_Chr*_bp15_35seq2501_*TRUE.png ./spectra_parts_15-35

mkdir spectra_parts_2-8
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/largeimages.png/All_spec1_*_Chr*_bp2_8seq2501_*TRUE.png ./spectra_parts_2-8

#----------
cd ~/scratch/repeats/Summary

mkdir Pielou_div
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*Pielou_div_norm.png ./Pielou_div

mkdir  Simpson_div
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*Simpson_div_norm.png ./Simpson_div

mkdir Shannon_div
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*Shannon_plot_norm.png ./Shannon_div

mkdir Shannon_div_25
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*roll_mean_Shannon_25.png ./Shannon_div_25

mkdir Shannon_div_100
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*roll_mean_Shannon_100.png ./Shannon_div_100

mkdir Shannon_div_1000
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*roll_mean_Shannon_1000.png ./Shannon_div_1000

mkdir Shannon_div_250
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*roll_mean_Shannon_250.png ./Shannon_div_250

mkdir Shannon_div_500
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*roll_mean_Shannon_500.png ./Shannon_div_500


#----------
cd ~/scratch/repeats/Summary

mkdir isochores
cp -v -u /home/celphin/scratch/repeats/input_chromosomes/*/isochore/*.fasta_png.1.png ./isochores

#----------
cd ~/scratch/repeats/Summary

mkdir DNAwalks_coloured_genes
cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/Genes_coloured_*DNAwalk*.png ./DNAwalks_coloured_genes

#-------------
cd ~/scratch/repeats/Summary

rm Centromere_summary.txt
cat /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/histograms/Centromere_*.txt > Centromere_summary.txt
more Centromere_summary.txt

rm Centromere_summary_Shannon.txt
cat /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*Centromere_MIN_Shannon.txt > Centromere_summary_Shannon.txt
more Centromere_summary_Shannon.txt

grep "cent25 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_25.txt
grep "cent100 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_100.txt
grep "cent250 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_250.txt
grep "cent500 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_500.txt
grep "cent1000 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_1000.txt
grep "centwind " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_wind.txt


#######################################################
# some example code

library(RepeatOBserverV1)
inpath="/home/celphin/scratch/repeats/input_chromosomes/Arab_COLCEN/chromosome_files/"
fname= "COLCEN_H0"
outpath="/home/celphin/scratch/repeats/output_chromosomes"
chromosome="Chr4"

x_cpu=1
pflag=FALSE
writeflag=FALSE
plotflag=FALSE

nam_list0 <- list.files(inpath)
nam_list1 <- tools::file_path_sans_ext(nam_list0)
nam_list1
nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
chr_list <- nam_list3[,1]
chr_list

# run chromosomes on different cpu
uni_chr_list <- unique(chr_list)
uni_chr_list 

# to run from initial fasta files
for (nam in chr_list){
print(nam)
run_plot_chromosome_parallel(nam=nam, fname=fname, inpath=inpath, outpath=outpath, pflag=pflag, plotflag=plotflag,  writeflag=writeflag, x_cpu=x_cpu)
}

for (nam in chr_list){
print(nam)
  write_All_spec_DNAwalk(nam=nam, fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
}

for (chromosome in uni_chr_list){
  merge_spectra(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
  join_chromosome_parts(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
}

#----------------------------------------
# To just plot the histograms, shannon diversity or heatmaps
# the following can be only run successfully after the above code or  the
# default functions (run_plot_chromosome_parallel, write_All_spec_DNAwalk, merge_spectra, join_chromosome_parts) 
# have run properly

# to re-plot the DNAwalks and Fourier heatmaps for a specific chromosome (can also wrap in loop as above)
chromosome = "ChrXX"

run_summary_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)

# to remake specific shannon diversity plots
run_diversity_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)

# to remake specific chromosome histogram plots
run_summary_hist(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath, atflag=TRUE)

################################
# to make plots of longer repeats > 2000bp up to 10 000bp
# or plots for a more specific range 

library(RepeatOBserverV1)
inpath="~/scratch/repeats/input_chromosomes/Einkorn/chromosome_files/"
fname= "Einkorn_H0"

outpath="~/scratch/repeats/output_chromosomes"
x_cpu=1
pflag=FALSE
writeflag=FALSE
plotflag=FALSE

nam_list0 <- list.files(inpath)
nam_list1 <- tools::file_path_sans_ext(nam_list0)
nam_list1
nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
chr_list <- nam_list3[,1]
chr_list

# run chromosomes on different cpu
uni_chr_list <- unique(chr_list)
uni_chr_list 

for (nam in chr_list){
print(nam)
run_long_repeats_NEW(nam=nam, fname=fname, inpath=inpath, outpath=outpath)
}

# or for a specific chromosome
chromosome="Chraeg1A"
run_long_repeats_NEW(nam=nam, fname=fname, inpath=inpath,outpath=outpath)

#################################
# to make 20Mbp plots of the Fourier transform data in a pdf

for (nam in chr_list){
print(nam)
run_20Mbpimag(nam=nam, fname=fname, inpath=inpath, outpath=outpath,  pflag=FALSE, plotflag=FALSE,writeflag=FALSE)
}

#############################################
# to make Zoomed in plots of Fourier spectra and DNAwalks

#-----------------------
# HA412 Chr 5 Zoom 35-40Mbp

# HA 89 Chr 5 35-40Mbp

# DNA walks for HA412 and HA89
# HA 412 Chr 5 37 130 000 - 37 150 000
# HA 89 Chr 5 35 520 000 - 35 550 000
# HA 89 Chr 5 38 255 000 - 38 275 000

inpath="~/scratch/repeats/input_chromosomes/Ha412/chromosome_files/"
fname= "Ha412_H0"
chromosome="Chr05"

# spectra
startbp=35005001
endbp=39000000

#----------------------

inpath="~/scratch/repeats/input_chromosomes/Sunflower/chromosome_files/"
fname= "Sunflower_H0"
chromosome="ChrHa8905"

# spectra
startbp=35005001
endbp=39000000

#----------------------
library(RepeatOBserverV1)

outpath="~/scratch/repeats/output_chromosomes"
x_cpu=1
pflag=FALSE
writeflag=FALSE
plotflag=FALSE

nam_list0 <- list.files(inpath)
nam_list1 <- tools::file_path_sans_ext(nam_list0)
nam_list1
nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
chr_list <- nam_list3[,1]
chr_list

# to plot the DNAwalks and Fourier transform heatmaps in the startbp to endbp range above
run_summary_plots_range(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath, startbp=startbp, endbp=endbp)

################################################
# to make Fourier transform heatmap plots with your own range of repeat lengths

# Example Fruit Fly

inpath="~/scratch/repeats/input_chromosomes/Fly/chromosome_files/"
fname= "Fly_H0"
chromosome="ChrDlow77"
part="part01"
nam="Fly_H0_ChrDlow77part01"
#----------------------
library(RepeatOBserverV1)

outpath="~/scratch/repeats/output_chromosomes"
x_cpu=1
pflag=FALSE
writeflag=FALSE
plotflag=FALSE

nam_list0 <- list.files(inpath)
nam_list1 <- tools::file_path_sans_ext(nam_list0)
nam_list1
nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
chr_list <- nam_list3[,1]
chr_list

# writes out file for DNA walk and All_spec (the Fourier transform)
write_All_spec_DNAwalk(nam=nam, fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)

# read in the All_spec part of interest

All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome,part,"_All_spec.txt"), check.names = FALSE))

# remove blank/zero columns in All_spec
colnames(All_spec0)[ncol(All_spec0)]

# https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
All_spec <- All_spec0[,which(!is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0)]

# https://r-graph-gallery.com/heatmap
#stats::heatmap(All_spec[,c(1:100)], Rowv=FALSE, Colv=FALSE)
pngflag=TRUE
ofi <- paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,part,"_All_spec_")

# set repeat lengths of interest : 17-30 bp
r1=base::c(17,30)

base::cat("\nofi:",ofi,"\n")
chromosome <- paste0(chromosome, part)
if(!base::is.null(r1))largeimagesub_NEW(All_spec,fname=fname, chromosome=chromosome,ofi,rangebp=r1,pngflag = pngflag, repround=1, part=1)


