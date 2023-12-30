# #---------------------------------------
# # update the code for summary files and run new code
#
# tmux new-session -s Repeats1
# tmux attach-session -t Repeats1
#
# cd /home/celphin/scratch/repeats/
#
#   #salloc -c20 --time 2:50:00 --mem 192000M --account def-cronk
#
#   module load StdEnv/2020 r/4.1.2
# R
#
# #------------------------------
# # test barplots, test shannon
#
# library(RepeatObserver)
# inpath="~/scratch/repeats/input_chromosomes/Ha412/chromosome_files/"
# fname= "Ha412_H0"
#
# outpath="~/scratch/repeats/output_chromosomes"
# x_cpu=1
# pflag=FALSE
# writeflag=FALSE
# plotflag=FALSE
#
# nam_list0 <- list.files(inpath)
# nam_list1 <- tools::file_path_sans_ext(nam_list0)
# nam_list1
# nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
# nam_list3 <- stringr::str_split(nam_list2[,3], "part", simplify =TRUE)
# chr_list <- nam_list3[,1]
# chr_list
#
# # run chromosomes on different cpu
# uni_chr_list <- unique(chr_list)
# uni_chr_list
# chromosome="Chr15"
#
# run_diversity_plots_no_telomere(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
# run_diversity_plots_35_2000_no_telomere(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
# run_diversity_plots(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
# run_diversity_plots_35_2000(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath)
#
# #--------------
# # look at output
#
# cd ~/scratch/repeats/output_chromosomes/Ha412_H0/Chr01/histograms/
#   more Centromere_Ha412_H0_Chr01_MIN_POWER_SUM_403_s_0.5std_1_35_2000.txt
# Ha412_H0 Chr01 2.5e+07 155215000
#
# cd ~/scratch/repeats/output_chromosomes/Ha412_H0/Chr01/
#   more Ha412_H0_Chr01_Centromere_MIN_Shannon.txt
# Ha412_H0_Chr01_cent25 57282501 159500000 Ha412_H0 Chr01
#
# ############################################
# # Summary join output files
#
# cd ~/scratch/repeats/Summary/Centromere_prediction
#
# rm Centromere_summary.txt
# cat /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/histograms/Centromere_*.txt > Centromere_summary.txt
# more Centromere_summary.txt
#
# awk 'NF>3' Centromere_summary.txt > Centromere_summary_1.txt
# awk 'NF<5' Centromere_summary_1.txt > Centromere_summary_2.txt
#
# #----------------------------------------
#
# rm Centromere_summary_Shannon.txt
# cat /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*Centromere_MIN_Shannon.txt > Centromere_summary_Shannon.txt
#
# rm Centromere_summary_Shannon_*
#   grep "cent25 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_25.txt
# grep "cent100 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_100.txt
# grep "cent250 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_250.txt
# grep "cent500 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_500.txt
# grep "cent1000 " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_1000.txt
# grep "centwind " Centromere_summary_Shannon.txt > Centromere_summary_Shannon_wind.txt
#
# awk 'NF>4' Centromere_summary_Shannon_500.txt > Centromere_summary_Shannon_500_1.txt
# awk 'NF>4' Centromere_summary_Shannon_wind.txt > Centromere_summary_Shannon_wind_1.txt
#
# #-----------------------
# cd ~/scratch/repeats/Summary/Centromere_prediction
#
# rm Centromere_summary_Shannon_35.txt
# cat /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*Centromere_MIN_Shannon_35.txt > Centromere_summary_Shannon_35.txt
#
# rm Centromere_summary_Shannon_*
#   grep "cent25 " Centromere_summary_Shannon_35.txt > Centromere_summary_Shannon_25_35.txt
# grep "cent100 " Centromere_summary_Shannon_35.txt > Centromere_summary_Shannon_100_35.txt
# grep "cent250 " Centromere_summary_Shannon_35.txt > Centromere_summary_Shannon_250_35.txt
# grep "cent500 " Centromere_summary_Shannon_35.txt > Centromere_summary_Shannon_500_35.txt
# grep "cent1000 " Centromere_summary_Shannon_35.txt > Centromere_summary_Shannon_1000_35.txt
# grep "centwind " Centromere_summary_Shannon_35.txt > Centromere_summary_Shannon_wind_35.txt
#
# awk 'NF>4' Centromere_summary_Shannon_500_35.txt > Centromere_summary_Shannon_500_35_1.txt
# awk 'NF>4' Centromere_summary_Shannon_wind_35.txt > Centromere_summary_Shannon_wind_35_1.txt
#
# #--------------------------------------
# rm Centromere_summary_Shannon_no_telo.txt
# cat /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*Centromere_MIN_Shannon_no_telo.txt > Centromere_summary_Shannon_no_telo.txt
# more Centromere_summary_Shannon.txt
#
# rm Centromere_summary_Shannon_*
#   grep "cent25 " Centromere_summary_Shannon_no_telo.txt > Centromere_summary_Shannon_25_no_telo.txt
# grep "cent100 " Centromere_summary_Shannon_no_telo.txt > Centromere_summary_Shannon_100_no_telo.txt
# grep "cent250 " Centromere_summary_Shannon_no_telo.txt > Centromere_summary_Shannon_250_no_telo.txt
# grep "cent500 " Centromere_summary_Shannon_no_telo.txt > Centromere_summary_Shannon_500_no_telo.txt
# grep "cent1000 " Centromere_summary_Shannon_no_telo.txt > Centromere_summary_Shannon_1000_no_telo.txt
# grep "centwind " Centromere_summary_Shannon_no_telo.txt > Centromere_summary_Shannon_wind_no_telo.txt
#
#
# awk 'NF>4' Centromere_summary_Shannon_wind_no_telo.txt> Centromere_summary_Shannon_500_no_telo_1.txt
# awk 'NF>4' Centromere_summary_Shannon_wind_no_telo.txt> Centromere_summary_Shannon_wind_no_telo_1.txt
#
#
# #-----------------------
# cd ~/scratch/repeats/Summary/Centromere_prediction
#
# rm Centromere_summary_Shannon_35_no_telo.txt
# cat /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*Centromere_MIN_Shannon_35_no_telo.txt > Centromere_summary_Shannon_35_no_telo.txt
# #more Centromere_summary_Shannon_35.txt
#
# rm Centromere_summary_Shannon_*
#   grep "cent25 " Centromere_summary_Shannon_35_no_telo.txt > Centromere_summary_Shannon_25_35_no_telo.txt
# grep "cent100 " Centromere_summary_Shannon_35_no_telo.txt > Centromere_summary_Shannon_100_35_no_telo.txt
# grep "cent250 " Centromere_summary_Shannon_35_no_telo.txt > Centromere_summary_Shannon_250_35_no_telo.txt
# grep "cent500 " Centromere_summary_Shannon_35_no_telo.txt > Centromere_summary_Shannon_500_35_no_telo.txt
# grep "cent1000 " Centromere_summary_Shannon_35_no_telo.txt > Centromere_summary_Shannon_1000_35_no_telo.txt
# grep "centwind " Centromere_summary_Shannon_35_no_telo.txt > Centromere_summary_Shannon_wind_35_no_telo.txt
#
# awk 'NF>4' Centromere_summary_Shannon_wind_35_no_telo.txt > Centromere_summary_Shannon_500_35_no_telo_1.txt
# awk 'NF>4' Centromere_summary_Shannon_wind_35_no_telo.txt > Centromere_summary_Shannon_wind_35_no_telo_1.txt
#
#
# #----------------------------
# cd ~/scratch/repeats/Summary/
#   mkdir Shannon_div_window
# cp -v -u /home/celphin/scratch/repeats/output_chromosomes/*/Chr*/*Shannon_div_window*.png ./Shannon_div_window
#
# #############################################
# # test join
# tmux new-session -s Repeats1
# tmux attach-session -t Repeats1
#
# module load StdEnv/2020 r/4.1.2
# R
#
# #####################################

library(RepeatObserver)

# setwd("~/scratch/repeats/Summary/Centromere_prediction/")

setwd("~/Cedar_transfers/Centromere_prediction")

#-----------------------------
# load barplot and shannon data
barplot_cent <- as.data.frame(utils::read.table("./Centromere_summary_2.txt"))

shannon_cent <- as.data.frame(utils::read.table("./Centromere_summary_Shannon_500_1.txt"))
shannon_35_cent <- as.data.frame(utils::read.table("./Centromere_summary_Shannon_500_35_no_telo_1.txt"))
shannon_no_telo_cent <- as.data.frame(utils::read.table("./Centromere_summary_Shannon_500_35_no_telo_1.txt"))
shannon_35_no_telo_cent <- as.data.frame(utils::read.table("./Centromere_summary_Shannon_500_35_no_telo_1.txt"))

shannon_cent_wind <- as.data.frame(utils::read.table("./Centromere_summary_Shannon_wind_1.txt"))
shannon_35_cent_wind <- as.data.frame(utils::read.table("./Centromere_summary_Shannon_wind_35_no_telo_1.txt"))
shannon_no_telo_cent_wind <- as.data.frame(utils::read.table("./Centromere_summary_Shannon_wind_35_no_telo_1.txt"))
shannon_35_no_telo_cent_wind <- as.data.frame(utils::read.table("./Centromere_summary_Shannon_wind_35_no_telo_1.txt"))

#---------------------------
# edit column names
# barplot
colnames(barplot_cent) <- c("fname", "chromosome", "centromere_barplot", "full_length_barplot")
barplot_cent$Mix <- paste0(barplot_cent$fname, "_", barplot_cent$chromosome)

#----------
shannon_cent <- shannon_cent[,-1]
colnames(shannon_cent) <- c("centromere_shannon", "full_length_shannon", "fname", "chromosome")
shannon_cent$Mix <- paste0(shannon_cent$fname, "_", shannon_cent$chromosome)

shannon_35_cent <- shannon_35_cent[,-1]
colnames(shannon_35_cent) <- c("centromere_shannon_35", "full_length_shannon", "fname", "chromosome")

shannon_no_telo_cent <- shannon_no_telo_cent[,-1]
colnames(shannon_no_telo_cent) <- c("centromere_shannon_notelo", "full_length_shannon", "fname", "chromosome")

shannon_35_no_telo_cent <- shannon_35_no_telo_cent[,-1]
colnames(shannon_35_no_telo_cent) <- c("centromere_shannon_35_notelo", "full_length_shannon", "fname", "chromosome")

#---------
shannon_cent_wind <- shannon_cent_wind[,-1]
colnames(shannon_cent_wind) <- c("centromere_shannon_wind", "full_length_shannon", "fname", "chromosome")

shannon_35_cent_wind <- shannon_35_cent_wind[,-1]
colnames(shannon_35_cent_wind ) <- c("centromere_shannon_wind_35", "full_length_shannon", "fname", "chromosome")

shannon_no_telo_cent_wind <- shannon_no_telo_cent_wind[,-1]
colnames(shannon_no_telo_cent_wind) <- c("centromere_shannon_wind_notelo", "full_length_shannon", "fname", "chromosome")

shannon_35_no_telo_cent_wind <- shannon_35_no_telo_cent_wind[,-1]
colnames(shannon_35_no_telo_cent_wind) <- c("centromere_shannon_wind_35_notelo", "full_length_shannon", "fname", "chromosome")

#---------

# merge
centromere_data <- as.data.frame(merge(barplot_cent, shannon_cent,  by=c("fname", "chromosome")))

# write out centromere data
utils::write.table(centromere_data, "./Centromeres_total.txt", append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = TRUE)


centromere_data1 <- merge(centromere_data, shannon_35_cent, by=c("fname", "full_length_shannon", "chromosome"))
centromere_data2 <- merge(centromere_data1, shannon_no_telo_cent, by=c("fname", "full_length_shannon","chromosome"))
centromere_data3 <- merge(centromere_data2, shannon_35_no_telo_cent ,by=c("fname","full_length_shannon", "chromosome"))
centromere_data4 <- merge(centromere_data3, shannon_cent_wind, by=c("fname", "full_length_shannon","chromosome"))
centromere_data5 <- merge(centromere_data4, shannon_35_cent_wind ,by=c("fname", "full_length_shannon","chromosome"))
centromere_data6 <- merge(centromere_data5, shannon_no_telo_cent_wind,  by=c("fname", "full_length_shannon","chromosome"))
centromere_data7 <- merge(centromere_data6, shannon_35_no_telo_cent_wind, by=c("fname", "full_length_shannon","chromosome"))

#-----------
# add in the known centromere data
known_cent <- as.data.frame(utils::read.table("./Centromeres_known.txt", sep = "\t", header=TRUE))

centromere_data7$Mix0 <- paste0(centromere_data2$fname, "_", centromere_data7$chromosome)
colnames(known_cent) <- c("Mix0", "centromere_known")

centromere_datax <- merge(centromere_data7, known_cent)

centromere_datax[c(1:3),]

centromere_data2 <- dplyr::distinct(centromere_datax)
#--------------
# make predictions in Mbp
centromere_data2$centromere_shannon <- centromere_data2$centromere_shannon/1e6
centromere_data2$centromere_shannon_35 <- centromere_data2$centromere_shannon_35/1e6
centromere_data2$centromere_shannon_notelo <- centromere_data2$centromere_shannon_notelo/1e6
centromere_data2$centromere_shannon_35_notelo <- centromere_data2$centromere_shannon_35_notelo/1e6
centromere_data2$centromere_shannon_wind <- centromere_data2$centromere_shannon_wind/1e6
centromere_data2$centromere_shannon_wind_35 <- centromere_data2$centromere_shannon_wind_35/1e6
centromere_data2$centromere_shannon_wind_notelo <- centromere_data2$centromere_shannon_wind_notelo/1e6
centromere_data2$centromere_shannon_wind_35_notelo <- centromere_data2$centromere_shannon_wind_35_notelo/1e6

centromere_data2$centromere_barplot <- centromere_data2$centromere_barplot/1e6
centromere_data2$full_length_shannon <- centromere_data2$full_length_shannon/1e6

centromere_data2[c(1:3),]

colnames(centromere_data2)[c(6,7,9)]
centromere_data2 <- centromere_data2[,-c(6,7,9)]

#--------------------
# write out centromere data
utils::write.table(centromere_data2 , "./Centromeres_Table_final.txt", append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = TRUE)

#---------------------------
# gather shannnon and barplot columns
 colnames(centromere_data2)
#install.packages("tidyr")
library(tidyr)
gathered_cent<- tidyr::gather(centromere_data2, "Technique", "Cent_Pos", 5:13)

centromere_data2 <- gathered_cent

#----------------------
# subtract known from predicted

centromere_data2$Cent_Pos <- as.numeric(centromere_data2$Cent_Pos)
centromere_data2$centromere_known <- as.numeric(centromere_data2$centromere_known)

centromere_data2$cent_accuracy <- abs(centromere_data2$Cent_Pos - centromere_data2$centromere_known)
centromere_data2$accuracy_norm <- (centromere_data2$cent_accuracy/centromere_data2$full_length_shannon)*100

#---------------------------
# calc how many get within certain percent of total size
centromere_data2$accuracy_count <- centromere_data2$accuracy_norm

centromere_data2$accuracy_count[which(centromere_data2$accuracy_norm<=15)] <-  1
centromere_data2$accuracy_count[which(centromere_data2$accuracy_norm>15)] <-  0

#---------------------------------
# calc how many total are right

centromere_data2$Mix <- paste0( centromere_data2$fname, "_", centromere_data2$Technique)

centromere_by_spp <- dplyr::group_by(centromere_data2, by=Mix)
max_error_by_spp <- dplyr::summarise(centromere_by_spp, max = max(cent_accuracy, na.rm = T),
                                     min = min(cent_accuracy, na.rm = T),
                                     mean = mean(cent_accuracy, na.rm = T),
                                     count = sum(accuracy_count, na.rm = T))

centromere_by_Tech <- dplyr::group_by(centromere_data2, by=Technique)
max_error_by_Tech <- dplyr::summarise(centromere_by_Tech, max = max(cent_accuracy, na.rm = T),
                                     min = min(cent_accuracy, na.rm = T),
                                     mean = mean(cent_accuracy, na.rm = T),
                                     count = sum(accuracy_count, na.rm = T))


# change count to percent of total number of chromosomes


#--------------------------------------
# keep only the two shannon and barplot techniques
unique(centromere_data3$Technique)
centromere_data2 <- centromere_data2[which(centromere_data2$Technique=="centromere_shannon" | centromere_data2$Technique=="centromere_barplot"),]

#-------------------------------
# remove Maize and Arab since doubled

centromere_data3 <- centromere_data2
centromere_data2 <- centromere_data3[-which(centromere_data3$fname=="Arab_H0"),]
centromere_data2 <- centromere_data2[-which(centromere_data2$fname=="Maize_H0"),]
centromere_data2 <- centromere_data2[-which(centromere_data2$fname=="Finch_H0"),]

centromere_data2 <- centromere_data2[-grep("ChrDath", centromere_data2$chromosome),]
centromere_data2 <- centromere_data2[-grep("ChrDpse", centromere_data2$chromosome),]
centromere_data2 <- centromere_data2[-grep("ChrDsub", centromere_data2$chromosome),]

#---------------------------------
# plot
col <- sort(unique(centromere_data2$Mix))
col[grep("centromere_barplot",col)] <- "red"
col[grep("centromere_shannon", col)] <- "lightblue"


png("./Centromere_by_spp.png", width = 1000, height = 600)
# sets the bottom, left, top and right margins
par(mar=c(15,4.1,4.1,2.1))
boxplot(cent_accuracy~Mix, data=centromere_data2, las=2, col=col, xlab="Spp_Technique", ylab="Centromere prediction accuracy")
stripchart(cent_accuracy~Mix, vertical = TRUE, method = "jitter", pch = 16,las = 2, col = "black", cex = 1.5, data = centromere_data2, add=TRUE)
dev.off()

png("./Centromere_by_spp_norm.png", width = 1000, height = 600)
# sets the bottom, left, top and right margins
par(mar=c(15,4.1,4.1,2.1))
boxplot(accuracy_norm~Mix, data=centromere_data2, las=2, col=col, xlab="Spp_Technique", ylab="Centromere prediction accuracy")
stripchart(accuracy_norm~Mix, vertical = TRUE, method = "jitter", pch = 16,las = 2, col = "black", cex = 1.5, data = centromere_data2, add=TRUE)
dev.off()


png("./Centromere_by_spp_norm_colour_TEch.png", width = 1000, height = 600)
library(ggplot2)
ggplot(data=centromere_data2, aes(x=Mix, y=accuracy_norm))+
  geom_boxplot(aes(x=fname, y=accuracy_norm))+
  facet_wrap(~Technique, scales = "free")+ theme_classic()+
  labs(fill="Technique")
dev.off()

png("./Centromere_by_spp_norm_colour.png", width = 1000, height = 600)
library(ggplot2)
ggplot(data=centromere_data2, aes(x=Mix, y=accuracy_norm))+
  geom_boxplot(aes(x=Technique, y=accuracy_norm))+
  facet_wrap(~fname, scales = "free")+ theme_classic()+
  labs(fill="Technique")
dev.off()

centromere_by_spp <- dplyr::group_by(centromere_data2, by=Mix)
max_error_by_spp2 <- dplyr::summarise(centromere_by_spp, max = max(cent_accuracy, na.rm = T),
                                      min = min(cent_accuracy, na.rm = T),
                                      mean = mean(cent_accuracy, na.rm = T),
                                      count = sum(accuracy_count, na.rm = T),
                                      total = length(accuracy_count))

centromere_by_Tech <- dplyr::group_by(centromere_data2, by=Technique)
max_error_by_Tech2 <- dplyr::summarise(centromere_by_Tech, max = max(cent_accuracy, na.rm = T),
                                       min = min(cent_accuracy, na.rm = T),
                                       mean = mean(cent_accuracy, na.rm = T),
                                       count = sum(accuracy_count, na.rm = T),
                                       total = nrow(accuracy_count))

#------------------------
# order by success in shannon or barplot

unique(centromere_data2$Mix)
centromere_data2$Mix2 <- factor(centromere_data2$Mix , levels=c("COLCEN_H0_centromere_barplot",    "COLCEN_H0_centromere_shannon",
                                                               "Human_H0_centromere_barplot",     "Human_H0_centromere_shannon",
                                                               "Fly_H0_centromere_barplot" ,      "Fly_H0_centromere_shannon",
                                                               "Rice_H0_centromere_barplot",      "Rice_H0_centromere_shannon",
                                                               "Brassica_H0_centromere_barplot",  "Brassica_H0_centromere_shannon",
                                                               "Strawberry_H0_centromere_barplot", "Strawberry_H0_centromere_shannon",
                                                               "Ha412_H0_centromere_barplot" ,    "Ha412_H0_centromere_shannon",
                                                               "Tomato_H0_centromere_barplot",    "Tomato_H0_centromere_shannon",
                                                               "Chicken_H0_centromere_barplot",   "Chicken_H0_centromere_shannon",
                                                               "Mouse_H0_centromere_barplot",     "Mouse_H0_centromere_shannon",
                                                               "Corn_H0_centromere_barplot",      "Corn_H0_centromere_shannon"))

col <- sort(unique(centromere_data2$Mix))
col[grep("centromere_barplot",col)] <- "red"
col[grep("centromere_shannon", col)] <- "lightblue"


png("./Centromere_by_spp_norm_subset_ordered.png", width = 1000, height = 600)
# sets the bottom, left, top and right margins
par(mar=c(15,4.1,4.1,2.1))
boxplot(accuracy_norm~Mix2, data=centromere_data2, las=2, col=col, xlab="Spp_Technique", ylab="Centromere prediction accuracy")
stripchart(accuracy_norm~Mix2, vertical = TRUE, method = "jitter", pch = 16,las = 2, col = "black", cex = 1.5, data = centromere_data2, add=TRUE)
dev.off()


#####################################
# Plot a subset for the Spp grouped by the techniques that worked

centromere_data2 <- centromere_data2[-which((centromere_data2$fname=="Corn_H0") & (centromere_data2$Technique=="centromere_barplot")),]
centromere_data2 <- centromere_data2[-which((centromere_data2$fname=="Mouse_H0") & (centromere_data2$Technique=="centromere_barplot")),]
centromere_data2 <- centromere_data2[-which((centromere_data2$fname=="Strawberry_H0") & (centromere_data2$Technique=="centromere_barplot")),]
centromere_data2 <- centromere_data2[-which((centromere_data2$fname=="Chicken_H0") & (centromere_data2$Technique=="centromere_barplot")),]


centromere_data2 <- centromere_data2[-which((centromere_data2$fname=="Tomato_H0") & (centromere_data2$Technique=="centromere_shannon")),]
centromere_data2 <- centromere_data2[-which((centromere_data2$fname=="Ha412_H0") & (centromere_data2$Technique=="centromere_shannon")),]

col <- sort(unique(centromere_data2$Mix))
col[grep("centromere_barplot",col)] <- "red"
col[grep("centromere_shannon", col)] <- "lightblue"

png("./Centromere_by_spp_norm_subset.png", width = 1000, height = 600)
# sets the bottom, left, top and right margins
par(mar=c(15,4.1,4.1,2.1))
boxplot(accuracy_norm~Mix, data=centromere_data2, las=2, col=col, xlab="Spp_Technique", ylab="Centromere prediction accuracy")
stripchart(accuracy_norm~Mix, vertical = TRUE, method = "jitter", pch = 16,las = 2, col = "black", cex = 1.5, data = centromere_data2, add=TRUE)
dev.off()



#-----------------------------

centromere_by_spp <- dplyr::group_by(centromere_data2, by=Mix)
max_error_by_spp3 <- dplyr::summarise(centromere_by_spp, max = max(cent_accuracy, na.rm = T),
                                     min = min(cent_accuracy, na.rm = T),
                                     mean = mean(cent_accuracy, na.rm = T),
                                     count = sum(accuracy_count, na.rm = T),
                                     total = length(accuracy_count))

centromere_by_Tech <- dplyr::group_by(centromere_data2, by=Technique)
max_error_by_Tech3 <- dplyr::summarise(centromere_by_Tech, max = max(cent_accuracy, na.rm = T),
                                      min = min(cent_accuracy, na.rm = T),
                                      mean = mean(cent_accuracy, na.rm = T),
                                      count = sum(accuracy_count, na.rm = T),
                                      total = nrow(accuracy_count))




