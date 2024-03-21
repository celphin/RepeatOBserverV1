#################################
# Notes for running R code on the server
# June 2023 
# Repeat runs for centromere paper
###################################

# Code to install the package in R

cd ~/scratch/repeats/

module load StdEnv/2020 r/4.1.2
R


#############################################
# Code to download and prepare the reference genome assembly from NCBI

# Maize Example 
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02029-9
# https://www.ncbi.nlm.nih.gov/assembly/GCA_029775835.1

# download the data to a folder called input_chromosomes

cd /~/scratch/

mkdir repeats/; cd repeats
mkdir output_chromosomes
mkdir input_chromosomes; cd input_chromosomes
mkdir Maize; cd Maize

# go to FTP directory for GenBank assembly and download the genomic.fna.gz - any fasta file for a reference genome assembly can work
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/775/835/GCA_029775835.1_ASM2977583v1/GCA_029775835.1_ASM2977583v1_genomic.fna.gz

# Unzip if the file ends in .gz
gunzip GCA_029775835.1_ASM2977583v1_genomic.fna.gz


#--------------------
# maybe script from here?

# To help R read in long sequences fold the data to be 60bp per line using seqkit
module load seqkit/2.3.1
seqkit seq -w 60 GCA_029775835.1_ASM2977583v1_genomic.fna > Maize_2.fasta

# make a directory to store the individual chromosome files in
mkdir chromosome_files/
cd ./chromosome_files/

# split the assembly fasta file into one fasta file per chromosome
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../Maize_2.fasta

# remove any chromosomes or scaffolds that are less than 5Mbp _ this can be set as low as 1Mbp
find -type f -size -5000000c -delete

# rename the chromosomes in the form Refspp_Haplotype_Chr#
# the 3 part naming of the files here is important for the reading into the program
# Refspp can be any name you use to identify this species but should not include any underscores
# Haplotype can be any way to record the haplotype, H0 if the haplotype is unknown, again this should not include any underscores
# Chromosome number should be your ID for the chromosome and should not include any underscores

# in this case "CP1167" is the start of the NCBI chromosome names and we renmae this to "Maize_H0_Chr"
# note if you do this each chromosome will have the NCBI ID number and not the chromosome number starting from 1 in the output
rename CP1167 Maize_H0_Chr *

# the files should also end in .fasta not .fna or .fa
rename .1.fasta .fasta *

# loop to split long chromosomes
# any chromosomes over 100Mbp long will be run in 100Mbp sections and rejoined
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

# all chromosomes should be labelled in the form:
# Maize_H0_Chr39part01.fasta

#########################
# Chr ID
# 1	39
# 2	40
# 3	41
# 4	42
# 5	43
# 6	44
# 7	45
# 8	46
# 9	47
# 10 48
	

###################################################
# The following code creates a script that runs the code needed pre-repeats analysis

nano /~/scratch/repeats/pre-repeats.sh

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
isochore -sequence ${chr} -outfile ${chr}.isochore -graph png -goutfile ${chr}_png
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

inpath=${qq}${pathname}/chromosome_files/${qq}
fname=${qq}${SPP}${qq}

# set the pathname to your computer/server
outpath="~/scratch/repeats/output_chromosomes"

# set the number of cpu you have access to -1
x_cpu=19

########################
# below here nothing should need to be changed unless you want to run something specific

library(RepeatObserver)

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
  library(RepeatObserver)
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
  library(RepeatObserver)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  merge_spectra(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
  join_chromosome_parts(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
  run_summary_hist(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
  run_summary_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
}

x_cpu=5

cl <- parallel::makeCluster(x_cpu)
results2 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), chromosome_summary,uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)


EOF1

# make R script with correct input
cat << EOF2 > centromere_histogram.R

inpath=${qq}${pathname}/chromosome_files/${qq}
fname=${qq}${SPP}${qq}

# set the pathname to your computer/server
outpath="~/scratch/repeats/output_chromosomes"

# set the number of cpu you have access to -1
# if you run out of memory increase the memory or decrease the number of cpu here
x_cpu=19

########################
# below here nothing should need to be changed unless you want to run something specific

library(RepeatObserver)

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
  library(RepeatObserver)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  run_summary_hist(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
}

cl <- parallel::makeCluster(x_cpu)
results2 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), chromosome_summary,uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)


EOF2


EOF

chmod +x /~/scratch/repeats/pre-repeats.sh

#------------------------------
# this code creates the post repeats script

nano /~/scratch/repeats/post_repeats.sh

#!/bin/bash

SPP_Hap=$1

cd /~/scratch/repeats/output_chromosomes/${SPP_Hap}
mkdir Summary_output
cd Summary_output

mkdir histograms
cp -v -u /~/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/histograms/*POWER_SUM*.pdf ./histograms/

mkdir spectra_total_merged
cp -v -u /~/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/${SPP_Hap}*All_spec*.png ./spectra_total_merged

mkdir boxplots
cp -v -u /~/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/5000bp_barplots/*_Chr*_POWER_SUM_seqval_*_s_1std_*_35_2000.pdf ./boxplots/

mkdir DNAwalks
cp -v -u /~/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/${SPP_Hap}*DNAwalk*.png ./DNAwalks

mkdir spectra_parts
cp -v -u  /~/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp35_2000seq2501_*TRUE.png ./spectra_parts
cp -v -u /~/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp15_35seq2501_*TRUE.png ./spectra_parts
cp -v -u /~/scratch/repeats/output_chromosomes/${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp2_8seq2501_*TRUE.png ./spectra_parts

EOF

chmod +x /~/scratch/repeats/post_repeats.sh



############################################
# slurm scripts

# make a script to submit the pre, post and R scripts to your server
# code below is setup for slurm

# Maize Example
cd /~/scratch/repeats/
mkdir scripts
cd /~/scratch/repeats/scripts

cat << EOF > Maize_cedar_repeats.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=192000M

module load StdEnv/2020
module load emboss/6.6.0
srun /~/scratch/repeats/pre-repeats.sh "/~/scratch/repeats/input_chromosomes/Maize" "Maize_H0"

cd /~/scratch/repeats/input_chromosomes/Maize 
module load StdEnv/2020 r/4.1.2
#Rscript /~/scratch/repeats/input_chromosomes/Maize/repeats_fourier.R
Rscript /~/scratch/repeats/input_chromosomes/Maize/centromere_histogram.R

srun /~/scratch/repeats/post_repeats.sh "Maize_H0"

EOF

sbatch Maize_cedar_repeats.sh


######################################
# Make Summary folder for all chromosomes

cd ~/scratch/repeats
mkdir Summary

cd ~/scratch/repeats/Summary

mkdir histograms
cp -v -u  /~/scratch/repeats/output_chromosomes/*/Chr*/histograms/*POWER_SUM*s_0.5std_1*.pdf ./histograms/

cat /~/scratch/repeats/output_chromosomes/*/Chr*/histograms/Centromere_*.txt > Centromere_summary.txt
more Centromere_summary.txt

mkdir spectra_total_merged
cp -v -u /~/scratch/repeats/output_chromosomes/*/Chr*/*All_spec_bp35_2000*.png ./spectra_total_merged

mkdir DNAwalks
cp -v -u /~/scratch/repeats/output_chromosomes/*/Chr*/*DNAwalk*.png ./DNAwalks

mkdir spectra_parts_35-2000
cp -v -u  /~/scratch/repeats/output_chromosomes/*/Chr*/largeimages.png/All_spec1_*_Chr*_bp35_2000seq2501_*TRUE.png ./spectra_parts_35-2000

mkdir spectra_parts_15-35
cp -v -u /~/scratch/repeats/output_chromosomes/*/Chr*/largeimages.png/All_spec1_*_Chr*_bp15_35seq2501_*TRUE.png ./spectra_parts_15-35

mkdir spectra_parts_2-8
cp -v -u /~/scratch/repeats/output_chromosomes/*/Chr*/largeimages.png/All_spec1_*_Chr*_bp2_8seq2501_*TRUE.png ./spectra_parts_2-8


########################################
