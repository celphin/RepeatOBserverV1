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
  while getopts ":i:f:h:c:m:" opt; do
    case $opt in
     i)
        species="$OPTARG"
        ;;
     f)
        fasta_file="$OPTARG"
        ;;
     h) 
        haplotype="$OPTARG"
        ;;
     c)
        cpu="$OPTARG"
        ;;
     m)
        memory="$OPTARG"
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
  path_name=$(pwd)
  mkdir output_chromosomes
  mkdir input_chromosomes; cd input_chromosomes
  mkdir $species; cd $species
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
  # remove any chromosomes or scaffolds that are less than 1Mbp 
  find -type f -size -1000000c -delete
  echo "The resulting chromosome files are: "
  ls    
  
  # rename parts
  count=1
  for filename in *; do
    new_name="${species}_${haplotype}_Chr$count.fasta"
    # Rename the file
    mv "$filename" "$new_name"
    ((count++))
  done
  echo "Files renamed"

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

  fi
}


#Writes and builds pre-repeats script in three parts:
  #1: General set up
  #2: Make Repeats_Fourier.R
  #3: Make Centromere_Histogram.R
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
x_cpu=1
cl <- parallel::makeCluster(x_cpu)
results2 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), chromosome_summary,uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)
print("repeats_fourier.R finished running")
EOF1
echo "repeats_fourier.R made"
EOF

#Part3: make centromere histogram script
cat << EOF >> pre-repeats.sh
cat << EOF2 > centromere_histogram.R
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
nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
chr_list <- nam_list3[,1]
chr_list

# run chromosomes on different cpu
uni_chr_list <- unique(chr_list)

chromosome_summary <-  function(x, uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath){
  library(RepeatOBserverV1)
  chromosome <<- uni_chr_list[x]
  print(chromosome)
  run_summary_hist(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
  run_summary_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
  run_diversity_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
}

cl <- parallel::makeCluster(x_cpu)
results2 <- parallel::parSapply(cl, base::seq_along(uni_chr_list), chromosome_summary,uni_chr_list=uni_chr_list, fname=fname, inpath=inpath, outpath=outpath)
print("centromere_histogram.R finished running")
EOF2
echo "centromere_histogram.R made"
echo "pre-repeats.sh finished running"

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
mkdir spectra_parts_35-2000
cp -v -u  ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp35_2000seq2501_*TRUE.png ./spectra_parts_35-2000

mkdir spectra_parts_15-35
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp15_35seq2501_*TRUE.png ./spectra_parts_15-35

mkdir spectra_parts_2-8
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/largeimages.png/All_spec1_*_Chr*_bp2_8seq2501_*TRUE.png ./spectra_parts_2-8

#------------------------------
mkdir spectra_total_merged
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*All_spec_bp35_2000*.png ./spectra_total_merged
mkdir histograms
cp -v -u  ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/histograms/*POWER_SUM*s_0.5std_1*.pdf ./histograms/
mkdir DNAwalks
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*DNAwalk*.png ./DNAwalks
#----------
#Make shannon_div directories: 100,1000, 250, 500
mkdir Shannon_div
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/\${SPP_Hap}_Chr*_Shannon_plot_norm.png ./Shannon_div
mkdir Shannon_div_100
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*roll_mean_Shannon_100.png ./Shannon_div_100
mkdir Shannon_div_1000
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*roll_mean_Shannon_1000.png ./Shannon_div_1000
mkdir Shannon_div_250
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*roll_mean_Shannon_250.png ./Shannon_div_250
mkdir Shannon_div_500
cp -v -u ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/*roll_mean_Shannon_500.png ./Shannon_div_500

#-------------

rm Centromere_summary.txt
cat ${path_name}/output_chromosomes/\${SPP_Hap}/Chr*/histograms/Centromere_*.txt > Centromere_summary.txt
#more Centromere_summary.txt

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
EOF

chmod +x post-repeats.sh
echo "post-repeats.sh file made"
}

run_all(){
  
############################################

  cd ${path_name}

  ${path_name}/pre-repeats.sh "${path_name}/input_chromosomes/${species}" "${species}_H0" "${cpu}"

  cd ${path_name}/input_chromosomes/${species}

  Rscript ${path_name}/input_chromosomes/${species}/repeats_fourier.R 
  Rscript ${path_name}/input_chromosomes/${species}/centromere_histogram.R

  ${path_name}/post-repeats.sh "${species}_H0"

  echo "RepeatOBserverV1 complete"
}


run "$@"

