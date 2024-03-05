#!/bin/bash

#Modifications: 
    #Need to hand modify the renaming section, search for: DO: RENAME
run() {
  #Get input:
  get_input "$@"
  echo "Running RepeatOBserverV1 on: $species"
  echo "with fasta file: $fasta_file"
  echo "With $cpu cpus over $memory MB"
  set_up_directories
  split_fasta
  build_pre_repeats
  build_post_repeats
  run_all 
}


#get_input:
  #Reads input from original send, 
  #TO DO: modify to account for CPU
get_input() {
  while getopts ":i:f:h:c:m:g:" opt; do
    case $opt in
     i)
        species="$OPTARG"
        ;;
     f)
        fasta_file="$OPTARG"
        ;;
     h) 
        Haplotype="$OPTARG"
        ;;
     c)
        cpu="$OPTARG"
        ;;
     m)
        memory="$OPTARG"
        ;;
     g)
        CGflag="$OPTARG"
        ;;
      \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
      :)
        echo "Option -$OPTARG requires an argument." >&2
        exit 1
      ;;
    esac
  done
}

#Set up directories:
set_up_directories() {
  # change haplotype name if CGflag=TRUE
  if [ ${CGflag} = TRUE ]
  then
    haplotype=${Haplotype}-CG
  else
    haplotype=${Haplotype}-AT
  fi

  path_name=$(pwd)
  mkdir output_chromosomes
  mkdir input_chromosomes; cd input_chromosomes
  mkdir ${species}_${haplotype}; cd ${species}_${haplotype}
  cp $path_name/$fasta_file .
  echo "Running in directory ${path_name}"
}

split_fasta() {
  # To help R read in long sequences fold the data to be 60bp per line using seqkit

  #seqkit seq -w 60 GCA_029775835.1_ASM2977583v1_genomic.fna > Maize_2.fasta
  seqkit seq -w 60 $fasta_file > ${species}_2.fasta

  if [ -d chromosome_files ]
  then
    echo "chromosome_files directory exists"
    cd ./chromosome_files/
  else
  # make a directory to store the individual chromosome files in
  mkdir chromosome_files/
  cd ./chromosome_files/
  
  # split the assembly fasta file into one fasta file per chromosome
  awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }'  ../${species}_2.fasta
  # remove any chromosomes or scaffolds that are less than 5Mbp 
  find -type f -size -5000000c -delete
  echo "The resulting chromosome files are: "
  ls    
  
  # rename chromosomes
  touch ../chromosome_renaming.txt
  echo "Old_Name	New_Name" >> ../chromosome_renaming.txt
  count=1
  for filename in *; do
    new_name="${species}_${haplotype}_Chr$count"
    # Rename the file
    mv "$filename" "$new_name"
    echo "$filename	$new_name" >> ../chromosome_renaming.txt
    ((count++))
  done
  echo "Files renamed"
  
  # check chromosome lengths, find all greater than 400Mbp and split into 400Mbp parts
  find -type f -size +400000000c  > ../list_long_chrnum.txt

  while IFS= read -r Chr; 
  do
  split -l 6666660 ${Chr} ${Chr}- --numeric-suffixes=1
  rm ${Chr}
  done < ../list_long_chrnum.txt
  wc -c *
  
  # rename parts
  # loop to split long (>100Mbp) chromosomes into parts
  # any chromosomes over 100Mbp long will be run in 100Mbp sections and rejoined
  ls > ../list_chromosomes.txt
  sed 's/.fasta//g' ../list_chromosomes.txt > ../list_chrnum.txt

  while IFS= read -r Chr; 
  do
  split -l 1666666 ${Chr} ${Chr}part --numeric-suffixes=1
  rm ${Chr}
  done < ../list_chrnum.txt

  wc -c *

  cd ..

  find ./chromosome_files/ -type f -exec mv {} {}".fasta" \;

  fi
}


#Writes and builds pre-repeats script in three parts:
  #1: General set up
  #2: Make Repeats_Fourier.R
  #3: Make Summary_plots.R
build_pre_repeats() {
    # The following code creates a script that runs the code needed pre-repeats analysis
  #cat << EOF done in several parts due to max length
cd ${path_name}

  #PART 1: Set up everything until Fourier transform script
cat << EOF > pre-repeats.sh
#!/bin/bash
## input directory set to just above chromosome_files directory
# inlude / at the end of pathname
pathname=\$1
SPP=\$2
cpu=\$3
CGflag=\$4

cd \${pathname}
echo "need to be in: \${pathname}"
echo "Currently in: \$(pwd)"

# get chromosome lengths
wc -c \${pathname}/chromosome_files/* > chromosome_lengths.txt

# make plots of CG isochores
# list of chromosomes
ls \${pathname}/chromosome_files/* > chr_list.txt

while IFS= read -r chr; 
do
isochore -sequence \${chr} -outfile \${chr}.isochore -window 100000 -shift 50000  -graph png -goutfile \${chr}_png
done < chr_list.txt

mkdir isochore
mv ./chromosome_files/*.isochore isochore
mv ./chromosome_files/*.png isochore

#TO DO: add if statements to file existence
# identify telomere locations - plants
grep -n "TTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTT" ./chromosome_files/* > plant_telomere_positions.txt

# identify telomere locations - animals
grep -n "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT" ./chromosome_files/* > animal_telomere_positions.txt

qq='"'

EOF

#Part 2: make fourier script
cat << EOF >> pre-repeats.sh

# make R script with correct input
cat << EOF1 > repeats_fourier.R

print("repeats_fourier.R starting")

library(RepeatOBserverV1)
inpath=\${qq}\${pathname}/chromosome_files/\${qq}
fname=\${qq}\${SPP}\${qq}
#----------------------------------------

if (\${CGflag}) {
AT_flag=FALSE
} else {
AT_flag=TRUE
}

outpath="${path_name}/output_chromosomes"
x_cpu=\${cpu} 
pflag=FALSE
writeflag=FALSE
plotflag=FALSE
atflag=TRUE

nam_list0 <- list.files(inpath)
nam_list1 <- tools::file_path_sans_ext(nam_list0)
nam_list1
nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
nam_list2
nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
nam_list3
chr_list <- nam_list3[,1]
chr_list

#--------------------
# run intial spectra if directories do not exist
for (nam in nam_list1){
 run_plot_chromosome_parallel(nam=nam, fname=fname, inpath=inpath, outpath=outpath, atflag=atflag, AT_flag=AT_flag, pflag=pflag, plotflag=plotflag,  writeflag=writeflag, x_cpu=x_cpu)
}

#-----------------------
# try to run writing of summary files on different cpu

print("write_all_spec starting")
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

print("write_all_spec complete")

#-----------------------
# merge chromosome parts on different cpu
uni_chr_list <- unique(chr_list)

print("join_parts starting")
parts_join <-  function(x, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath){
  library(RepeatOBserverV1)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  merge_spectra(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
  join_chromosome_parts(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath)
}

x_cpu=2
cl <- parallel::makeCluster(x_cpu)
results2 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), parts_join, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)

print("join_parts complete")

#----------------------------------
# plot centromeres 
print("centromere summary starting")

centromere_summary <-  function(x, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath){
  library(RepeatOBserverV1)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  run_diversity_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
  run_summary_hist(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
}

x_cpu=\${cpu}
cl <- parallel::makeCluster(x_cpu)
results3 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), centromere_summary, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)

print("centromere summary complete")

#----------------------------------
# plot chromosomes - summary plots

print("summary plots starting")

chromosome_summary <-  function(x, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath){
  library(RepeatOBserverV1)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  run_summary_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
}

x_cpu=\${cpu}
cl <- parallel::makeCluster(x_cpu)
results4 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), chromosome_summary, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)

print("summary plots complete")

print("repeats_fourier.R finished running")
EOF1
echo "repeats_fourier.R made"
EOF

#Part 3: make Summary plots script
cat << EOF >> pre-repeats.sh

# make R script for summary plots
cat << EOF2 > Summary_plots.R

print("Summary_plots.R starting")

library(RepeatOBserverV1)
inpath=\${qq}\${pathname}/chromosome_files/\${qq}
fname=\${qq}\${SPP}\${qq}
#----------------------------------------

outpath="${path_name}/output_chromosomes"
x_cpu=\${cpu} 
pflag=FALSE
writeflag=FALSE
plotflag=FALSE

nam_list0 <- list.files(inpath)
nam_list1 <- tools::file_path_sans_ext(nam_list0)
nam_list1
nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
nam_list2
nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
nam_list3
chr_list <- nam_list3[,1]
chr_list

uni_chr_list <- unique(chr_list)

#----------------------------------
# Plot each chromosome DNA walk and fourier transforms again

print("summary plots 2 starting")
chromosome_summary <-  function(x, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath){
  library(RepeatOBserverV1)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  run_summary_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
}

x_cpu=\${cpu}
cl <- parallel::makeCluster(x_cpu)
results4 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), chromosome_summary, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)

print("summary plots 2 complete")
#----------------------------------
# try to plot centromeres again

print("centromere summary 2 starting")

centromere_summary <-  function(x, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath){
  library(RepeatOBserverV1)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  run_summary_hist(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
  run_diversity_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
}

x_cpu=\${cpu}
cl <- parallel::makeCluster(x_cpu)
results3 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), centromere_summary, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)

print("centromere summary 2 complete")

#---------------------------------
print("plot all chromosomes starting")

# plot all chromosomes shannon diversity and histograms in one plot
plot_all_chromosomes(fname=fname, inpath=inpath, outpath=outpath)

print("plot all chromosomes complete")

print("Summary_plots.R finished running")
EOF2

echo "Summary_plots.R made"
EOF

chmod +x pre-repeats.sh
echo "pre-repeats.sh file made"
}

build_post_repeats() {
cd ${path_name}
touch post-repeats.sh
cat << EOF > post-repeats.sh
#!/bin/bash

SPP_Hap=\$1

cd ${path_name}/output_chromosomes/\${SPP_Hap}/
mkdir Summary_output; cd Summary_output
cd ${path_name}/output_chromosomes/\${SPP_Hap}/Summary_output
#-----------
mkdir spectra; cd spectra
mkdir spectra_parts_35-2000
mkdir spectra_parts_15-35
mkdir spectra_parts_2-8
mkdir spectra_total_merged; cd ..

cp -v -u  ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp35_2000seq2501_*TRUE.png ./spectra/spectra_parts_35-2000
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp15_35seq2501_*TRUE.png ./spectra/spectra_parts_15-35
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp2_8seq2501_*TRUE.png ./spectra/spectra_parts_2-8
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*All_spec_bp35_2000*.png ./spectra/spectra_total_merged

#--------------------------
mkdir histograms
cp -v -u  ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/histograms/*_histogram_POWER_SUM_seqval_*s_0.5std_1*.png ./histograms/

rm Centromere_summary.txt
cat ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/histograms/Centromere_*.txt > ./histograms/Centromere_histograms_summary.txt

#-----------------
mkdir DNAwalks; cd DNAwalks
mkdir 2D
mkdir 1D; cd ..
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*DNAwalk2D*.png ./DNAwalks/2D
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*DNAwalk1D*.png ./DNAwalks/1D

#----------
#Make shannon_div directories: 100,1000, 250, 500

mkdir Shannon_div; cd Shannon_div

mkdir Shannon_div_5kbp
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/\${SPP_Hap}_Chr*_Shannon_plot_norm.png ./Shannon_div_5kbp
mkdir Shannon_div_500kbp
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*roll_mean_Shannon_100.png ./Shannon_div_500kbp
mkdir Shannon_div_5Mbp
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*roll_mean_Shannon_1000.png ./Shannon_div_5Mbp
mkdir Shannon_div_1.25Mbp
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*roll_mean_Shannon_250.png ./Shannon_div_1.25Mbp
mkdir Shannon_div_2.5Mbp
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*roll_mean_Shannon_500.png ./Shannon_div_2.5Mbp

rm Centromere_summary_Shannon.txt
cat ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*Centromere_MIN_Shannon.txt > Centromere_summary_Shannon.txt

grep "cent25 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_25.txt
grep "cent100 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_100.txt
grep "cent250 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_250.txt
grep "cent500 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_500.txt
grep "cent1000 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_1000.txt

grep "centwind " Centromere_summary_Shannon_35_no_telo.txt > Centromere_summary_Shannon_wind_35_no_telo.txt
mkdir Shannon_div_window
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*Shannon_div_window*.png ./Shannon_div_window

cd ..

#-------------
# copy over the raw Shannon and Histogram data

mkdir output_data
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/Total_dnawalk_every50_\${SPP_Hap}_Chr*.txt ./output_data
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/Total_\${SPP_Hap}_Chr*_All_spec_merged.txt ./output_data
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/\${SPP_Hap}_Chr*_Shannon_div.txt ./output_data
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/histograms/\${SPP_Hap}_Chr*_Histogram_input*.txt ./output_data

#------------------------
# copy over chromosome renaming table
cp -v -u ${path_name}/input_chromosomes/\${SPP_Hap}/chromosome_renaming.txt .

#-------------------------
mkdir isochores
cp -v -u ${path_name}/input_chromosomes/${SPP_Hap}/isochore/*.png ./isochores

EOF

chmod +x post-repeats.sh
echo "post-repeats.sh file made"

}

############################################

run_all(){

  cd ${path_name}

  ${path_name}/pre-repeats.sh "${path_name}/input_chromosomes/${species}_${haplotype}" "${species}_${haplotype}" "${cpu}" "${CGflag}"

  cd ${path_name}/input_chromosomes/${species}_${haplotype}

  Rscript ${path_name}/input_chromosomes/${species}_${haplotype}/repeats_fourier.R 

  ${path_name}/post-repeats.sh "${species}_${haplotype}"

  Rscript ${path_name}/input_chromosomes/${species}_${haplotype}/Summary_plots.R 

  ${path_name}/post-repeats.sh "${species}_${haplotype}"

  # list and remove empty folders
  cd ${path_name}/output_chromosomes/${species}_${haplotype}
  find . -type d -empty -print
  find . -type d -empty -delete

  echo "RepeatOBserverV1 complete"
}


run "$@"


