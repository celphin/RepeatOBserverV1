# Test code

##############################
# run code

#-----------------------
# Custom settings
inpath="/home/celphin/scratch/repeats/input_chromosomes/Arabidopsis/chromosome_files/"
fname= "Arab_H0"
outpath="/home/celphin/scratch/repeats/output_chromosomes"
pflag=FALSE
writeflag=FALSE
plotflag=FALSE
x_cpu=19

#--------------------------
library(RepeatObserver)

chr_list0 <- list.files(inpath)
chr_list <- tools::file_path_sans_ext(chr_list0)
chr_list <- chr_list[c(4)]
chr_list

for (nam in chr_list){
  print(nam)
  run_plot_chromosome_parallel(nam=nam, fname=fname, inpath=inpath, outpath=outpath, pflag=pflag, plotflag=plotflag,  writeflag=writeflag, x_cpu=x_cpu)
}
# done

#----------------------------

for (nam in chr_list){
  print(nam)
  runs_run_barplots(nam=nam, fname=fname, inpath=inpath, outpath=outpath, pflag=pflag, plotflag=plotflag,  writeflag=writeflag)
}

# works - maybe comment out bedGraph making in run_barplots? still prints?
#--------------------------

for (nam in chr_list){
  print(nam)
  run_20Mbpimag(nam=nam, fname=fname, inpath=inpath, outpath=outpath,  pflag=pflag, plotflag=plotflag,  writeflag=writeflag)
}

# works
#--------------------------

for (nam in chr_list){
  print(nam)
  run_long_repeats_NEW(nam=nam, fname=fname, inpath=inpath, outpath=outpath, pflag=pflag, plotflag=plotflag,  writeflag=writeflag)
}

# works
# 5000bp_spectra.pdf/ is blank?
# 5000bp_spectra.txt/ is blank?

###################################
