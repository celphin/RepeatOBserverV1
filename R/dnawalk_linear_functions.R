
#' build_DNAwalk
#'
#' Builds a DNA walk for any fasta file of a DNA sequnece
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
build_DNAwalk <- function(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath) {

  # get total length of chromosome
  dna1 <- utils::read.table(paste0(inpath, fname,"/chromosome_files/",fname,"_",chromosome, ".fasta"), sep = "\t", header = FALSE, check.names = FALSE)

  # Load DNA sequence
  dna1 <- dna1[-1, ]
  waxy <- paste0(dna1, collapse = "\n")
  waxy <- stringr::str_replace_all(waxy, "([\n])", "")
  waxy <- unlist(strsplit(waxy, split = ""))
  x <- base::casefold(waxy, upper = FALSE)

  # Crop region
  #x <- x[start:end]

  # Build DNA walk
  ATval <- rep(0, length(x))
  CGval <- rep(0, length(x))

  Aindices <- grep("a", x)
  Tindices <- grep("t", x)
  Cindices <- grep("c", x)
  Gindices <- grep("g", x)

  ATval[Aindices] <- 1
  ATval[Tindices] <- -1

  CGval[Cindices] <- 1
  CGval[Gindices] <- -1

  atwalk <- cumsum(ATval)
  cgwalk <- cumsum(CGval)

  DNAwalk_long <- cbind(atwalk, cgwalk)
  dnawalk <- as.data.frame(DNAwalk_long)
  DNAwalk_long_rows <- cbind(dnawalk, c(1:nrow(dnawalk)))
  colnames(DNAwalk_long_rows) <- c("AT", "CG", "Position")

  # Optional: Save result for each chunk
  out_file <- paste0(outpath, "/", fname,"_", chromosome, "_total_dnawalk.txt")
  write.table(DNAwalk_long_rows, file = out_file, sep = "\t", row.names = FALSE, col.names = TRUE)

  grDevices::png(paste0(outpath, "/png/", fname,"_", chromosome, "_DNAwalk.png"), height=1000, width=1000)
  plot(DNAwalk_long_rows$AT, DNAwalk_long_rows$CG, main=fname )
  grDevices::dev.off()

  print(max(DNAwalk_long_rows$Position))

}

#' calc_R2_DNAwalk
#'
#' Calculates the slope and R^2 for any DNAwalk
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
calc_R2_DNAwalk <- function(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath) {

  # read in total dnawalk
  DNAwalk_long_rows <- utils::read.table(paste0(outpath, "/",  fname, "_", chromosome, "_total_dnawalk.txt"), sep = "\t", header = TRUE)

  head(DNAwalk_long_rows)
  tail(DNAwalk_long_rows)

  #----------------
  # run R^2 for total chromosome

  # sort by x
  orderedAT_DNAwalk <- DNAwalk_long_rows[order(DNAwalk_long_rows[,1],decreasing=FALSE),]
  head(orderedAT_DNAwalk)
  gc()

  # calculate R^2
  orderedAT_DNAwalk <- as.data.frame(orderedAT_DNAwalk)
  ml = lm(AT~CG, data=orderedAT_DNAwalk)
  summary(ml)$r.squared
  # SLOPE
  slope <- coef(ml)["CG"]
  #print(slope)

  # Finch2 Chr7 = 0.9017472
  # Finch3 Chr7 = 0.9343643
  # Slope 1.810583

  # arabidopsis 0.4400775
  # slope 0.8166362

  gc()

  #----------------------
  # run R^2 for each 1Mbp piece

  # Define the window size (1Mb) and the total length
  window_size_bps <- 1000000  # 1Mbp
  total_length_bps <- DNAwalk_long_rows[nrow(DNAwalk_long_rows), 3]  # 100Mb in base pairs

  # Number of data points (since data is spaced every 50bp)
  points_per_window <- window_size_bps / 1
  num_windows <- total_length_bps / window_size_bps  # Number of 5Mb windows

  # Store the R^2 values for each window
  r_squared_values <- numeric(num_windows)
  slope_values <- numeric(num_windows)

  # Loop through each 5Mb window
  for (i in 1:num_windows) {

    # Define the start and end positions for the current window in terms of data points
    start_index <- (i - 1) * points_per_window + 1
    end_index <- i * points_per_window

    # Subset the data for this window
    window_data <- DNAwalk_long_rows[start_index:end_index, ]

    # Sort the data by the first column (x-values, which is presumably the position or coordinate)
    ordered_window_data <- window_data[order(window_data[,1], decreasing = FALSE), ]

    # Convert the data into a data frame for modeling
    df_window <- as.data.frame(ordered_window_data)

    # Run the linear regression for AT ~ CG
    ml <- lm(AT ~ CG, data = df_window)

    # Store the R^2 value for this window
    r_squared_values[i] <- summary(ml)$r.squared
    slope_values[i] <- coef(ml)["CG"]
    gc()
  }

  # Print out the R^2 values for each window
  print("R^2")
  print(r_squared_values)

  print("Slope")
  print(slope_values)

  R2_slopes <- cbind(r_squared_values, slope_values)

  write.table(R2_slopes, file = paste0(outpath, "/",  fname, "_", chromosome,"_r2_slope.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)

  non_zero_r2 <- r_squared_values[r_squared_values != 0]
  average_r2 <- mean(non_zero_r2)
  #print(average_r2)

  non_zero_slope <- slope_values[slope_values != 0]
  average_slope <- mean(non_zero_slope)
  #print(average_slope)

}


#' calc_R2_genes
#'
#' Calculates an R^2 value for each gene in a gff3 file
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
calc_R2_genes <- function(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath,
                          gff3_path, gff3_chromosome) {

  # import the gff3 annotation file
  genes_gff30 <- utils::read.table(gff3_path, sep = "\t", header=FALSE, check.names = FALSE)
  colnames(genes_gff30) <- c("chrnum", "Program", "genetype", "Start", "End", "dot1", "Strand", "dot2", "ID")
  genes_gff3 <- genes_gff30[,-c(2,6,8,9)]
  genes_gff3 <- as.data.frame(genes_gff3)

  unique(genes_gff3$genetype)
  unique(genes_gff3$chrnum)
  unique(genes_gff3$Strand)

  # extract only the chromosome's genes
  chr_genes_gff3 <- genes_gff3[which(genes_gff3$chrnum==gff3_chromosome),]
  chr_genes_gff3$ID <- paste0(chr_genes_gff3$genetype, "_", chr_genes_gff3$Start, "_", chr_genes_gff3$End)
  head(chr_genes_gff3)
  ordered_genes_gff3 <- chr_genes_gff3[order(chr_genes_gff3[,3],decreasing=FALSE),]
  head(ordered_genes_gff3)
  # remove duplicates
  ordered_genes_gff3 <- ordered_genes_gff3[!duplicated(ordered_genes_gff3$ID), ]


  #-------------------------------
  # First, convert DNAwalk_long to a data frame and name the columns
  DNAwalk_long_rows <- utils::read.table(paste0(outpath, "/",  fname, "_", chromosome, "_total_dnawalk.txt"), sep = "\t", header = TRUE)
  DNAwalk_long_rows <- as.data.frame(DNAwalk_long_rows)
  colnames(DNAwalk_long_rows) <- c("AT", "CG", "Position")

  dnawalk_df <- DNAwalk_long_rows[seq(0, nrow(DNAwalk_long_rows), by = 50), ]
  head(dnawalk_df)

  # Prepare an empty data frame to store R² results
  gene_r2_results <- data.frame(
    GeneID = character(),
    R_squared = numeric(),
    Slope = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop over each gene in the gene list
  for (i in 1:nrow(ordered_genes_gff3)) {
    gene_id <- ordered_genes_gff3$ID[i]
    start_pos <- ordered_genes_gff3$Start[i]
    end_pos <- ordered_genes_gff3$End[i]

    # Subset DNA walk data for this gene region
    gene_data <- dnawalk_df[dnawalk_df$Position >= start_pos & dnawalk_df$Position <= end_pos, ]

    # If there's enough data (at least 5 rows), perform regression
    if (nrow(gene_data) > 5) {
      model <- lm(AT ~ CG, data = gene_data)
      r2_value <- summary(model)$r.squared
      slope_value <- coef(model)["CG"]
    } else {
      r2_value <- NA
      slope_value <- NA
    }

    # Add result to the results data frame
    gene_r2_results <- rbind(gene_r2_results, data.frame(GeneID = gene_id, R_squared = r2_value, Slope = slope_value))
    # print(gene_id)
    # print(r2_value)
    # print(slope_value)
  }

  # View the results
  print(gene_r2_results)

  # Working!
  gene_r2_results_unique <- gene_r2_results[!duplicated(gene_r2_results$GeneID), ]

  # Write results to a file
  utils::write.table(x=gene_r2_results_unique, file=paste0(outpath, "/",  fname, "_", chromosome, "_gene_r2_results.txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)

  #-----------------
  # Figure out averages
  library(dplyr)
  library(stringr)

  # Step 1: Extract gene type from GeneID
  gene_r2_results_unique <- gene_r2_results_unique %>%
    mutate(gene_type = str_extract(GeneID, "^[^_]+"))

  # Step 2: Group by gene type and calculate mean R²
  gene_avg_r2 <- gene_r2_results_unique %>%
    group_by(gene_type) %>%
    summarise(mean_R_squared = mean(R_squared, na.rm = TRUE)) %>%
    arrange(desc(mean_R_squared))

  # View the result
  print(gene_avg_r2)

  # Step 3: Group by gene type and calculate mean slope
  gene_avg_slope <- gene_r2_results_unique %>%
    group_by(gene_type) %>%
    summarise(mean_slope = mean(Slope, na.rm = TRUE)) %>%
    arrange(desc(mean_slope))

  # get a count for each category

  # View the result
  print(gene_avg_slope)

}


#' calc_R2_TEs
#'
#' Calculates an R^2 value for each TE in a gff3 file from EDTA
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
calc_R2_TEs <- function(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath,
                        EDTA_path, EDTA_chromosome) {

  #Import the EDTA file
  EDTA_gff30 <- utils::read.table(EDTA_path, sep = "\t", header=FALSE, check.names = FALSE)
  colnames(EDTA_gff30) <- c("chrnum", "Program", "TE_type", "Start", "End", "bp", "Strand", "dot", "ID")
  EDTA_gff3 <- EDTA_gff30[,-c(2,6,8,9)]
  EDTA_gff3 <- as.data.frame(EDTA_gff3)

  unique(EDTA_gff3$genetype)
  unique(EDTA_gff3$chrnum)
  unique(EDTA_gff3$Strand)

  # extract only the chromosome's genes and TEs
  chr_TE_gff3 <- EDTA_gff3[which(EDTA_gff3$chrnum==gff3_chromosome),]
  chr_TE_gff3$ID <- paste0(chr_TE_gff3$TE_type, "_", chr_TE_gff3$Start, "_", chr_TE_gff3$End)
  head(chr_TE_gff3)
  ordered_TE_gff3 <- chr_TE_gff3[order(chr_TE_gff3[,3],decreasing=FALSE),]
  head(ordered_TE_gff3)
  # remove duplicates
  ordered_TE_gff3 <- ordered_TE_gff3[!duplicated(ordered_TE_gff3$ID), ]

  #-------------------------------

  # Prepare an empty data frame to store R² results
  TE_r2_results <- data.frame(
    TEID = character(),
    R_squared = numeric(),
    Slope = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop over each gene in the gene list
  for (i in 1:nrow(ordered_TE_gff3)) {
    TE_id <- ordered_TE_gff3$ID[i]
    start_pos <- ordered_TE_gff3$Start[i]
    end_pos <- ordered_TE_gff3$End[i]

    # Subset DNA walk data for this gene region
    TE_data <- dnawalk_df[dnawalk_df$Position >= start_pos & dnawalk_df$Position <= end_pos, ]

    # If there's enough data (at least 5 rows), perform regression
    if (nrow(TE_data) > 5) {
      model <- lm(AT ~ CG, data = TE_data)
      r2_value <- summary(model)$r.squared
      slope_value <- coef(model)["CG"]
    } else {
      r2_value <- NA
      slope_value <- NA
    }

    # Add result to the results data frame
    TE_r2_results <- rbind(TE_r2_results, data.frame(TEID = TE_id, R_squared = r2_value, Slope = slope_value))
    # print(gene_id)
    # print(r2_value)
    # print(slope_value)
  }

  # View the results
  print(TE_r2_results)

  TE_r2_results_unique <- TE_r2_results[!duplicated(TE_r2_results$TEID), ]

  # Write results to a file
  utils::write.table(x=TE_r2_results_unique, file=paste0(outpath, "/",  fname, "_", chromosome, "_TE_r2_results_unique.txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)

  #-------------------------------
  # Figure out averages

  library(dplyr)
  library(stringr)

  # Step 1: Extract TE type from GeneID
  TE_r2_results_unique <- TE_r2_results_unique %>%
    mutate(TE_type = str_extract(TEID, "^[^_]+"))

  # Step 2: Group by TE type and calculate mean R²
  TE_avg_r2 <- TE_r2_results_unique %>%
    group_by(TE_type) %>%
    summarise(mean_R_squared = mean(R_squared, na.rm = TRUE)) %>%
    arrange(desc(mean_R_squared))

  # View the result
  print(TE_avg_r2)

  # Step 3: Group by TE type and calculate mean slope
  TE_avg_slope <- TE_r2_results_unique %>%
    group_by(TE_type) %>%
    summarise(mean_slope = mean(Slope, na.rm = TRUE)) %>%
    arrange(desc(mean_slope))

  # View the result
  print(TE_avg_slope)

}





#' DNAwalks_with_genes
#'
#' Creates 2D and 1D  DNA walks of segements of chromosomes showing the genes plotted in blue (forward), red (reverse), black (intergenic) and yellow (overlapping forward and reverse).
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
DNAwalks_with_genes <- function(chromosome=chromosome, fname=fname, inpath=inpath,
                                outpath=outpath, gff3_path=gff3_path, cytchr=cytchr,
                                start=start, end=end){

  # import the gff3 annotation file
  cytobands0 <-utils::read.table(gff3_path, sep = "\t", header=FALSE, check.names = FALSE)

  colnames(cytobands0) <- c("chrnum", "Program", "genetype", "Start", "End", "dot1", "Strand", "dot2", "ID")

  cytobands <- cytobands0[,-c(2,6,8)]

  cytobands<- as.data.frame(cytobands)
  unique(cytobands$genetype)
  unique(cytobands$chrnum)
  unique(cytobands$Strand)

  # extract only the CDS
  CDS_cytobands <- cytobands[which(cytobands$chrnum==cytchr),]
  chr_cytobands <- CDS_cytobands[which(CDS_cytobands$genetype=="mRNA"),]

  #nrow(cytobands)
  #nrow(CDS_cytobands)
  #nrow(chr_cytobands)
  # 4133

  # select only genes in the start to end range
  chr_cytobands <- chr_cytobands[which(chr_cytobands$Start>=start),]
  chr_cytobands <- chr_cytobands[which(chr_cytobands$Start<=end),]

  # how many genes?
  nrow(chr_cytobands)
  # 4000
  #--------------------------------
  # set these regions as colours
  #walk_colours <- c("#00FFFF", "#FF00FF", "black", "#CCFF00")
  walk_colours <- c("red", "blue", "black", "#CCFF00")
  cytoband_types <- c("+", "-", "not_gene", "both")
  cyto_type_colours <- as.data.frame(cbind(cytoband_types, walk_colours))

  chr_cytobands
  #chr_cytobands <- Pos_chr_cytobands
  # order by starting position
  chr_cytobands$Start <- as.numeric(chr_cytobands$Start)
  chr_cytobands$End <- as.numeric(chr_cytobands$End)
  chr_cytobands <- chr_cytobands[order(chr_cytobands$Start),]

  walk_bands <- NULL
  skip=FALSE

  # add in the starting values
  chr_cytobands$Start0 <- chr_cytobands$Start
  chr_cytobands$Start <- chr_cytobands$Start-start
  chr_cytobands$End0 <- chr_cytobands$End
  chr_cytobands$End <- chr_cytobands$End-start
  START0 = chr_cytobands$Start[1]
  reptimes <- round(START0)
  colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == "not_gene")]
  coltmp <- rep(colour, reptimes)
  walk_bands <- c(coltmp, walk_bands)

  for (i in 1:nrow(chr_cytobands)) {
    START = chr_cytobands$Start[i]
    END = chr_cytobands$End[i]
    START2 = chr_cytobands$Start[i+1]
    END2 = chr_cytobands$End[i+1]
    print(skip)
    print(i)
    if (!is.na(START2)){
      if (skip==FALSE){
        reptimes <- END-START
        print(paste0("Gene = ",reptimes))
        colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == chr_cytobands$Strand[i])]
        coltmp <- rep(colour, reptimes)
        walk_bands <- c(walk_bands, coltmp)
        print(START)
        print(END)
        print(length(walk_bands))
        skip=FALSE
        if (END < START2){  # genes have space between them
          reptimes2 <- START2-END
          print(paste0("Space = ",reptimes2))
          colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == "not_gene")]
          coltmp <- rep(colour, reptimes2)
          walk_bands <- c(walk_bands, coltmp)
          print(START)
          print(END)
          print(length(walk_bands))
          skip=FALSE
        }
        if (END == START2){ # no overlap but no space, nothing to add but don't skip
          skip=FALSE
        }
        if (END > START2){
          if (END >= END2){	  # genes overlap completely
            skip=TRUE
          }
          if (END < END2){ # genes overlap partially
            reptimes3 <- START2-END
            print(paste0("Overlap = ",reptimes3))
            colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == chr_cytobands$Strand[i])]
            coltmp <- rep(colour, abs(reptimes3))
            walk_bands <- walk_bands[-c((length(walk_bands)- abs(reptimes3)+1):length(walk_bands))]
            walk_bands <- c(walk_bands, coltmp)
            print(START)
            print(END)
            print(length(walk_bands))
            skip=FALSE
          }
        }
      } else{
        START = chr_cytobands$Start[i]
        END = chr_cytobands$End[i]
        START2 = chr_cytobands$Start[i+1]
        END2 = chr_cytobands$End[i+1]
        if (END2 > length(walk_bands)) {
          skip=FALSE
          print("Next Partial overlap")
          if (START2 >= length(walk_bands)) {
            reptimes4 <- START2-length(walk_bands)
            print(paste0("Space = ",reptimes4))
            colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == "not_gene")]
            coltmp <- rep(colour, reptimes4)
            walk_bands <- c(walk_bands, coltmp)
            print(START)
            print(END)
            print(length(walk_bands))
            skip=FALSE
          }
          if (START2 < length(walk_bands)) {
            reptimes5 <- START2-length(walk_bands)
            print(paste0("Overlap = ",reptimes5))
            colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == chr_cytobands$Strand[i])]
            coltmp <- rep(colour, abs(reptimes5))
            walk_bands <- walk_bands[-c((length(walk_bands)- abs(reptimes5)+1):length(walk_bands))]
            walk_bands <- c(walk_bands, coltmp)
            print(START)
            print(END)
            print(length(walk_bands))
            skip=FALSE
          }
        }
        else {
          skip=TRUE
          print("Another complete overlap")
          print(START)
          print(END)
          print(length(walk_bands))
        }
      }
    }else{ # START2 is NA and loop is done
      reptimes6 <- END-length(walk_bands)
      if (reptimes6 >0){
        print(paste0("Space = ",reptimes6))
        colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == "not_gene")]
        coltmp <- rep(colour, reptimes6)
        walk_bands <- c(walk_bands, coltmp)
        print(START)
        print(END)
        print(length(walk_bands))
        print("Loop complete")
      }else{
        print(START)
        print(END)
        print(length(walk_bands))
        print("Loop complete")
      }

    }
  }

  #utils::write.table(walk_bands, file=paste0(outpath, "/", fname, "/", chromosome,"/", fname, "_", chromosome,"_DNAwalkcolours.txt"), append = FALSE, sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)

  walk_bands <- as.data.frame(walk_bands)

  #----------------------
  # load in DNA sequence
  dna1 <- utils::read.table(paste0(inpath, fname, "_", chromosome,".fasta"), sep = "\t", header=FALSE, check.names = FALSE)

  # remove header
  dna1 <- dna1[-1,]

  # concatenate
  waxy <- paste0(dna1, collapse="\n")
  waxy<-stringr::str_replace_all(waxy, "([\n])", "")
  waxy<-(base::unlist(base::strsplit(waxy,split="")))

  # fold upper case
  x <- base::casefold(waxy, upper=F)

  # crop to only keep start-end values
  x <- x[c(start:end)]
  # test with random generated sequence
  #dna1 <- utils::read.table(paste0(inpath,"Random_sequence.csv"), sep = ".", header=FALSE, check.names = FALSE)
  #x <- c(dna1$V1)
  #fname="TEST"

  # build the DNA walk

  ATval<- base::rep(0,base::length(base::c(x)))
  CGval<- base::rep(0,base::length(base::c(x)))
  ACval<- base::rep(0,base::length(base::c(x)))
  TGval<- base::rep(0,base::length(base::c(x)))
  gc()
  MAGval<- base::rep(0,base::length(base::c(x)))
  MATval<- base::rep(0,base::length(base::c(x)))

  Cindices<-base::grep("c",x=x)
  Gindices<-base::grep("g",x=x)
  Tindices<-base::grep("t",x=x)
  Aindices<-base::grep("a",x=x)

  gc()

  ATval[Cindices]<-  0
  ATval[Gindices]<-  0
  ATval[Aindices]<-  1;  MAGval[Aindices]<-  1
  ATval[Tindices]<- (-1);  MAGval[Tindices]<- -1

  CGval[Cindices]<-  1;  MAGval[Cindices]<-  1
  CGval[Gindices]<- (-1);  MAGval[Gindices]<- -1
  CGval[Aindices]<-  0
  CGval[Tindices]<-  0

  gc()
  ACval[Cindices]<-  (-1);  MATval[Cindices]<- -1
  ACval[Gindices]<-  0
  ACval[Aindices]<-  1;  MATval[Aindices]<-  1
  ACval[Tindices]<-  0

  TGval[Cindices]<-  0
  TGval[Gindices]<- (-1);  MATval[Gindices]<- -1
  TGval[Aindices]<-  0
  TGval[Tindices]<-  1;  MATval[Tindices]<-  1

  #Walklist<-base::list(atwalk=ATval,cgwalk=CGval,dnawalk=MAGval)
  Walklist<-base::list(atwalk=ATval,cgwalk=CGval, acwalk=ACval,tgwalk=TGval)

  # convert to cumulative walk

  atwalk<-base::cumsum(Walklist$atwalk)   #cumulative sum of AT single walk values
  cgwalk<-base::cumsum(Walklist$cgwalk)    #cumulative sum of CG single walk values
  acwalk<-base::cumsum(Walklist$acwalk)   #cumulative sum of AT single walk values
  tgwalk<-base::cumsum(Walklist$tgwalk)    #cumulative sum of CG single walk values

  #DNAwalk_long <- cbind(atwalk, cgwalk) #, acwalk, tgwalk)

  #---------------------------------------
  #dnawalk <- as.data.frame(DNAwalk_long)

  gc()
  # set region to plot
  length(atwalk)
  #chr_cytobands$Start
  walk_bands$x <- walk_bands$walk_bands
  length(walk_bands$x)

  # write 2D walk to png
  ofile <-  paste0(outpath, "/", fname, "/", chromosome,"/", "Genes_coloured_", fname, "_", chromosome,"_", start, "_", end)
  # Dryas specific
  # ofile <-  paste0(outpath, "/",  "Genes_coloured_", fname, "_", chromosome,"_", start, "_", end)

  ofile_AT_CG <- paste0(ofile, "_2D_DNAwalk_AT_CG.png")
  base::cat("\n ouput to", ofile_AT_CG)
  grDevices::png(file=ofile_AT_CG, width = 700, height = 500)
  plot(atwalk, cgwalk, cex=0.05, xlab="AT DNA Walk",  ylab="CG DNA Walk", col = walk_bands$x[c(1:length(atwalk))])
  grDevices::dev.off()
  gc()

  genomepos <- c(1:length(atwalk))+start

  # plot 1D DNAwalks
  ofile_AT <-  paste0(ofile, "_1D_DNAwalk_AT.png")
  base::cat("\n ouput to", ofile_AT)
  grDevices::png(file=ofile_AT, width = 700, height = 500)
  plot(genomepos, atwalk, cex=0.05, xlab="Genome Position",  ylab="AT DNA Walk", col =  walk_bands$x[c(1:length(atwalk))])
  grDevices::dev.off()

  ofile_CG <-  paste0(ofile, "_1D_DNAwalk_CG.png")
  base::cat("\n ouput to", ofile_CG)
  grDevices::png(file=ofile_CG, width = 700, height = 500)
  plot(genomepos, cgwalk, cex=0.05, xlab="Genome Position",  ylab="CG DNA Walk", col =  walk_bands$x[c(1:length(atwalk))])
  grDevices::dev.off()

  #------------------
  # plot the smoothed spline for AT
  spar=0.6
  lowpass1.spline <- stats::smooth.spline(genomepos, atwalk, spar = spar) ## Control spar for amount of smoothing
  # plot 1D DNAwalks
  ofile_AT <-  paste0(ofile, "_1D_DNAwalk_AT_with_spline.png")
  base::cat("\n ouput to", ofile_AT)
  grDevices::png(file=ofile_AT, width = 700, height = 500)
  plot(genomepos, atwalk,  cex=0.05, xlab="Genome Position",  ylab="AT DNA Walk", col =  walk_bands$x[c(1:length(atwalk))])
  graphics::lines(stats::predict(lowpass1.spline, genomepos), col = "red", lwd = 1.5)
  grDevices::dev.off()

  ofile_AT <-  paste0(ofile, "_1D_DNAwalk_AT_smoothed.png")
  base::cat("\n ouput to", ofile_AT)
  grDevices::png(file=ofile_AT, width = 700, height = 500)
  highpass1 <- atwalk - stats::predict(lowpass1.spline, genomepos)$y
  #graphics::lines(genomepos, highpass1, lwd =  2)
  base::plot(genomepos, highpass1, type="l", pch=15, lwd =  3, xlab="Genome Position", ylab="Smoothed AT DNA Walk")
  grDevices::dev.off()

  #------------------
  # other

  # ofile_AC_TG <- paste0(ofile, "_2D_DNAwalk_AC_TG.png")
  # base::cat("\n ouput to", ofile_AC_TG)
  # grDevices::png(file=ofile_AC_TG, width = 3000, height = 3000)
  # plot(acwalk, tgwalk, col = walk_bands$x[c(1:length(atwalk))]) #"black") # walk_bands$x[c(start:end)])
  # grDevices::dev.off()
  # gc()
  #
  # ofile_AC <-  paste0(ofile, "_1D_DNAwalk_AC.png")
  # base::cat("\n ouput to", ofile_AC)
  # grDevices::png(file=ofile_AC, width = 3000, height = 1000)
  # plot(genomepos, acwalk, col = walk_bands$x[c(1:length(atwalk))]) #"black") # walk_bands$x[c(start:end)])
  # grDevices::dev.off()
  #
  # ofile_TG <-  paste0(ofile, "_1D_DNAwalk_TG.png")
  # base::cat("\n ouput to", ofile_TG)
  # grDevices::png(file=ofile_TG, width = 3000, height = 1000)
  # plot(genomepos, tgwalk, col = walk_bands$x[c(1:length(atwalk))]) #"black") # walk_bands$x[c(start:end)])
  # grDevices::dev.off()

  # should also try
  # AG and TC

}

#' gene_slopes
#'
#' Plots graphs of the slopes of each DNA walk for each gene
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
gene_slopes <- function(chromosome=chromosome, fname=fname, inpath=inpath,
                        outpath=outpath, gff3_path=gff3_path, cytchr=cytchr,
                        start=start, end=end){


  # import the gff3 annotation file
  cytobands0 <-utils::read.table(gff3_path, sep = "\t", header=FALSE, check.names = FALSE)

  colnames(cytobands0) <- c("chrnum", "Program", "genetype", "Start", "End", "dot1", "Strand", "dot2", "ID")

  cytobands <- cytobands0[,-c(2,6,8)]

  cytobands<- as.data.frame(cytobands)
  unique(cytobands$genetype)
  unique(cytobands$chrnum)
  unique(cytobands$Strand)

  # extract only the CDS
  CDS_cytobands <- cytobands[which(cytobands$chrnum==cytchr),]
  chr_cytobands <- CDS_cytobands[which(CDS_cytobands$genetype=="mRNA"),]

  #nrow(cytobands)
  #nrow(CDS_cytobands)
  #nrow(chr_cytobands)
  # 4133

  # select only genes in the start to end range
  chr_cytobands <- chr_cytobands[which(chr_cytobands$Start>=start),]
  chr_cytobands <- chr_cytobands[which(chr_cytobands$Start<=end),]

  # how many genes?
  nrow(chr_cytobands)
  # 4000
  #--------------------------------
  # set these regions as colours
  #walk_colours <- c("#00FFFF", "#FF00FF", "black", "#CCFF00")
  walk_colours <- c("red", "blue", "black", "#CCFF00")
  cytoband_types <- c("+", "-", "not_gene", "both")
  cyto_type_colours <- as.data.frame(cbind(cytoband_types, walk_colours))

  chr_cytobands
  #chr_cytobands <- Pos_chr_cytobands
  # order by starting position
  chr_cytobands$Start <- as.numeric(chr_cytobands$Start)
  chr_cytobands$End <- as.numeric(chr_cytobands$End)
  chr_cytobands <- chr_cytobands[order(chr_cytobands$Start),]

  walk_bands <- NULL
  skip=FALSE

  #--------------------------------------------------
  # load in DNA sequence
  dna1 <- utils::read.table(paste0(inpath, fname, "_", chromosome,".fasta"), sep = "\t", header=FALSE, check.names = FALSE)

  # remove header
  dna1 <- dna1[-1,]

  # concatenate
  waxy <- paste0(dna1, collapse="\n")
  waxy<-stringr::str_replace_all(waxy, "([\n])", "")
  waxy<-(base::unlist(base::strsplit(waxy,split="")))

  # fold upper case
  x <- base::casefold(waxy, upper=F)

  # crop to only keep start-end values
  x <- x[c(start:end)]
  # test with random generated sequence
  #dna1 <- utils::read.table(paste0(inpath,"Random_sequence.csv"), sep = ".", header=FALSE, check.names = FALSE)
  #x <- c(dna1$V1)
  #fname="TEST"

  # build the DNA walk

  ATval<- base::rep(0,base::length(base::c(x)))
  CGval<- base::rep(0,base::length(base::c(x)))
  ACval<- base::rep(0,base::length(base::c(x)))
  TGval<- base::rep(0,base::length(base::c(x)))
  AGval<- base::rep(0,base::length(base::c(x)))
  TCval<- base::rep(0,base::length(base::c(x)))

  gc()
  MAGval<- base::rep(0,base::length(base::c(x)))
  MATval<- base::rep(0,base::length(base::c(x)))
  MACval<- base::rep(0,base::length(base::c(x)))

  Cindices<-base::grep("c",x=x)
  Gindices<-base::grep("g",x=x)
  Tindices<-base::grep("t",x=x)
  Aindices<-base::grep("a",x=x)

  gc()

  ATval[Cindices]<-  0
  ATval[Gindices]<-  0
  ATval[Aindices]<-  1;  MAGval[Aindices]<-  1
  ATval[Tindices]<- (-1);  MAGval[Tindices]<- -1

  CGval[Cindices]<-  1;  MAGval[Cindices]<-  1
  CGval[Gindices]<- (-1);  MAGval[Gindices]<- -1
  CGval[Aindices]<-  0
  CGval[Tindices]<-  0

  gc()
  ACval[Cindices]<-  (-1);  MATval[Cindices]<- -1
  ACval[Gindices]<-  0
  ACval[Aindices]<-  1;  MATval[Aindices]<-  1
  ACval[Tindices]<-  0

  TGval[Cindices]<-  0
  TGval[Gindices]<- (-1);  MATval[Gindices]<- -1
  TGval[Aindices]<-  0
  TGval[Tindices]<-  1;  MATval[Tindices]<-  1

  gc()
  AGval[Cindices]<-  0
  AGval[Gindices]<-  (-1); MACval[Gindices]<- -1
  AGval[Aindices]<-  1;  MACval[Aindices]<-  1
  AGval[Tindices]<-  0

  TCval[Cindices]<-  (-1);  MACval[Cindices]<- -1
  TCval[Gindices]<-  0
  TCval[Aindices]<-  0
  TCval[Tindices]<-  1;  MACval[Tindices]<-  1

  #Walklist<-base::list(atwalk=ATval,cgwalk=CGval,dnawalk=MAGval)
  Walklist<-base::list(atwalk=ATval,cgwalk=CGval, acwalk=ACval,tgwalk=TGval, agwalk=AGval,tcwalk=TCval)

  # convert to cumulative walk

  atwalk<-base::cumsum(Walklist$atwalk)   #cumulative sum of AT single walk values
  cgwalk<-base::cumsum(Walklist$cgwalk)    #cumulative sum of CG single walk values

  acwalk<-base::cumsum(Walklist$acwalk)   #cumulative sum of AC single walk values
  tgwalk<-base::cumsum(Walklist$tgwalk)    #cumulative sum of TG single walk values

  agwalk<-base::cumsum(Walklist$agwalk)   #cumulative sum of AG single walk values
  tcwalk<-base::cumsum(Walklist$tcwalk)    #cumulative sum of TC single walk values

  #DNAwalk_long <- cbind(atwalk, cgwalk) #, acwalk, tgwalk)

  #------------------------------------------
  # calculate slopes at and cg

  genomepos <- c(1:length(atwalk))+start
  DNAwalk_long<- as.data.frame(cbind(genomepos, atwalk, cgwalk))

  Genes_long <- as.data.frame(chr_cytobands[,c(1,3,4,5)])

  colnames(DNAwalk_long) <- c("Start", "atwalk", "cgwalk")
  Gene_starts <- dplyr::left_join(Genes_long, DNAwalk_long, by="Start")

  colnames(DNAwalk_long) <- c("End", "atwalk", "cgwalk")
  Gene_ends <- dplyr::left_join(Genes_long, DNAwalk_long, by="End")

  colnames(Gene_starts) <- c("chrnum", "Start", "End", "Strand", "Start_atpos", "Start_cgpos")
  colnames(Gene_ends) <- c("chrnum", "Start", "End", "Strand", "End_atpos", "End_cgpos")

  Gene_walks <- as.data.frame(merge(Gene_starts, Gene_ends))

  Gene_walks <- dplyr::distinct(Gene_walks)

  Gene_walks$Length <- (Gene_walks$End-Gene_walks$Start)


  Gene_walks$AT_slope <- (Gene_walks$End_atpos - Gene_walks$Start_atpos)/(Gene_walks$End-Gene_walks$Start)
  Gene_walks$CG_slope <- (Gene_walks$End_cgpos - Gene_walks$Start_cgpos)/(Gene_walks$End-Gene_walks$Start)
  Gene_walks$AT_CG_slope <- (Gene_walks$End_atpos - Gene_walks$Start_atpos)/(Gene_walks$End_cgpos-Gene_walks$Start_cgpos)

  Slope_by_strand <- dplyr::group_by(Gene_walks, by=Strand)
  Slopes_per_strand_AT <- dplyr::summarise(Slope_by_strand, ATmax = max(AT_slope, na.rm = T),
                                           ATmin = min(AT_slope, na.rm = T),
                                           ATmean = mean(AT_slope , na.rm = T))

  Slopes_per_strand_CG <- dplyr::summarise(Slope_by_strand, CGmax = max(CG_slope, na.rm = T),
                                           CGmin = min(CG_slope, na.rm = T),
                                           CGmean = mean(CG_slope , na.rm = T))
  #------------------------------------
  # think about slope of best fit line through walk - might fix the short genes that get affected by a few values
  colnames(DNAwalk_long) <- c("genomepos", "atwalk", "cgwalk")
  Gene_walks$best_fit_AT_slope <- Gene_walks$AT_slope
  Gene_walks$best_fit_CG_slope <- Gene_walks$AT_slope
  for (i in 1:nrow(Gene_walks)){
    if (Gene_walks$Start[i]<Gene_walks$End[i]){
      DNAwalk_subset <- DNAwalk_long[which(DNAwalk_long$genomepos >= Gene_walks$Start[i] & DNAwalk_long$genomepos <=Gene_walks$End[i]),]
      Gene_walks$best_fit_AT_slope[i] <- coef(lm(DNAwalk_subset$atwalk~DNAwalk_subset$genomepos))[2]
      Gene_walks$best_fit_CG_slope[i] <- coef(lm(DNAwalk_subset$cgwalk~DNAwalk_subset$genomepos))[2]
    } else {
      DNAwalk_subset <- DNAwalk_long[which(DNAwalk_long$genomepos <= Gene_walks$Start[i] & DNAwalk_long$genomepos >=Gene_walks$End[i]),]
      Gene_walks$best_fit_AT_slope[i] <- coef(lm(DNAwalk_subset$atwalk~DNAwalk_subset$genomepos))[2]
      Gene_walks$best_fit_CG_slope[i] <- coef(lm(DNAwalk_subset$cgwalk~DNAwalk_subset$genomepos))[2]
    }
  }

  #-----------------------------------
  # plots
  png(paste0(outpath, "/",fname, "_gene_fitslopes_at.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_AT_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="AT Slope")
  stripchart(best_fit_AT_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_fitslopes_cg.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_CG_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="CG Slope")
  stripchart(best_fit_CG_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_fitslopes_at_cg.png"), width = 1000, height = 1000)
  plot(best_fit_CG_slope ~ best_fit_AT_slope, data=Gene_walks, pch=16)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_fitslopes_AT_length.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_AT_slope ~ Length, data=Gene_walks, las=2, xlab="Gene Length", ylab="AT Slope")
  stripchart(best_fit_AT_slope ~ Length, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_fitslopes_AT_length_strand.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_AT_slope ~ Length+Strand, data=Gene_walks, las=2, xlab="Length/Strand", ylab="AT Slope")
  stripchart(best_fit_AT_slope ~ Length+Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()


  png(paste0(outpath, "/",fname, "test.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_AT_slope ~ AT_slope, data=Gene_walks, las=2, xlab="Strand", ylab="AT Slope")
  stripchart(best_fit_AT_slope ~ AT_slope, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()
  #---------
  png(paste0(outpath, "/",fname, "_Gene_start_pos.png"), width = 300, height = 15000)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  #boxplot(Start ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="Genome Position")
  stripchart(Start ~ Strand, data=Gene_walks, vertical = TRUE, pch = 16,las = 1, col = "black", cex = 1.5)#, add=TRUE)
  dev.off()

  #-------------------------------------------
  # plots
  png(paste0(outpath, "/",fname, "_gene_slopes_at.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(AT_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="AT Slope")
  stripchart(AT_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  # this was wrong
  png(paste0(outpath, "/",fname, "_gene_slopes_cg.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(CG_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="CG Slope")
  stripchart(CG_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_slopes_at_cg.png"), width = 1000, height = 1000)
  plot(CG_slope ~ AT_slope, data=Gene_walks, pch=16)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_slopes_AT_length.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(AT_slope ~ Length, data=Gene_walks, las=2, xlab="Gene Length", ylab="AT Slope")
  stripchart(AT_slope ~ Length, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_slopes_AT_length_strand.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(AT_slope ~ Length+Strand, data=Gene_walks, las=2, xlab="Gene Length/Strand", ylab="AT Slope")
  stripchart(AT_slope ~ Length+Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()


  # #------------------------------------------
  # # calculate slopes ag and tc
  #
  # genomepos <- c(1:length(agwalk))+start
  # DNAwalk_long<- as.data.frame(cbind(genomepos, agwalk, tcwalk))
  #
  # Genes_long <- as.data.frame(chr_cytobands[,c(1,3,4,5)])
  #
  # colnames(DNAwalk_long) <- c("Start", "agwalk", "tcwalk")
  # Gene_starts <- dplyr::left_join(Genes_long, DNAwalk_long, by="Start")
  #
  # colnames(DNAwalk_long) <- c("End", "agwalk", "tcwalk")
  # Gene_ends <- dplyr::left_join(Genes_long, DNAwalk_long, by="End")
  #
  # colnames(Gene_starts) <- c("chrnum", "Start", "End", "Strand", "Start_agpos", "Start_tcpos")
  # colnames(Gene_ends) <- c("chrnum", "Start", "End", "Strand", "End_agpos", "End_tcpos")
  #
  # Gene_walks <- as.data.frame(merge(Gene_starts, Gene_ends))
  #
  # Gene_walks <- dplyr::distinct(Gene_walks)
  #
  #
  # Gene_walks$ag_slope <- (Gene_walks$End_agpos - Gene_walks$Start_agpos)/(Gene_walks$End-Gene_walks$Start)
  # Gene_walks$tc_slope <- (Gene_walks$End_tcpos - Gene_walks$Start_tcpos)/(Gene_walks$End-Gene_walks$Start)
  # Gene_walks$ag_tc_slope <- (Gene_walks$End_agpos - Gene_walks$Start_agpos)/(Gene_walks$End_tcpos-Gene_walks$Start_tcpos)
  #
  # Gene_walks$Length <- (Gene_walks$End-Gene_walks$Start)
  #
  # Slope_by_strand <- dplyr::group_by(Gene_walks, by=Strand)
  # Slopes_per_strand_ag <- dplyr::summarise(Slope_by_strand, agmax = max(ag_slope, na.rm = T),
  #                                          agmin = min(ag_slope, na.rm = T),
  #                                          agmean = mean(ag_slope , na.rm = T))
  #
  # Slopes_per_strand_tc <- dplyr::summarise(Slope_by_strand, tcmax = max(tc_slope, na.rm = T),
  #                                          tcmin = min(tc_slope, na.rm = T),
  #                                          tcmean = mean(tc_slope , na.rm = T))
  #
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ag.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(ag_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="ag Slope")
  # stripchart(ag_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_tc.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(tc_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="tc Slope")
  # stripchart(ag_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ag_tc.png"), width = 1000, height = 1000)
  # plot(tc_slope ~ ag_slope, data=Gene_walks, pch=16)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ag_length.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(ag_slope ~ Length, data=Gene_walks, las=2, xlab="Strand", ylab="ag Slope")
  # stripchart(ag_slope ~ Length, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # #----------------------
  # LinearAG <- aov(ag_slope ~ Strand, data = Gene_walks)
  # summary(LinearAG)
  #
  # LinearTC <- aov(tc_slope ~ Strand, data = Gene_walks)
  # summary(LinearTC)
  #
  # #------------------------------------------
  # # calculate slopes ac and tg
  #
  # genomepos <- c(1:length(acwalk))+start
  # DNAwalk_long<- as.data.frame(cbind(genomepos, acwalk, tgwalk))
  #
  # Genes_long <- as.data.frame(chr_cytobands[,c(1,3,4,5)])
  #
  # colnames(DNAwalk_long) <- c("Start", "acwalk", "tgwalk")
  # Gene_starts <- dplyr::left_join(Genes_long, DNAwalk_long, by="Start")
  #
  # colnames(DNAwalk_long) <- c("End", "acwalk", "tgwalk")
  # Gene_ends <- dplyr::left_join(Genes_long, DNAwalk_long, by="End")
  #
  # colnames(Gene_starts) <- c("chrnum", "Start", "End", "Strand", "Start_acpos", "Start_tgpos")
  # colnames(Gene_ends) <- c("chrnum", "Start", "End", "Strand", "End_acpos", "End_tgpos")
  #
  # Gene_walks <- as.data.frame(merge(Gene_starts, Gene_ends))
  #
  # Gene_walks <- dplyr::distinct(Gene_walks)
  #
  #
  # Gene_walks$ac_slope <- (Gene_walks$End_acpos - Gene_walks$Start_acpos)/(Gene_walks$End-Gene_walks$Start)
  # Gene_walks$tg_slope <- (Gene_walks$End_tgpos - Gene_walks$Start_tgpos)/(Gene_walks$End-Gene_walks$Start)
  # Gene_walks$ac_tg_slope <- (Gene_walks$End_acpos - Gene_walks$Start_acpos)/(Gene_walks$End_tgpos-Gene_walks$Start_tgpos)
  #
  # Slope_by_strand <- dplyr::group_by(Gene_walks, by=Strand)
  # Slopes_per_strand_ac <- dplyr::summarise(Slope_by_strand, acmax = max(ac_slope, na.rm = T),
  #                                          acmin = min(ac_slope, na.rm = T),
  #                                          acmean = mean(ac_slope , na.rm = T))
  #
  # Slopes_per_strand_tg <- dplyr::summarise(Slope_by_strand, tgmax = max(tg_slope, na.rm = T),
  #                                          tgmin = min(tg_slope, na.rm = T),
  #                                          tgmean = mean(tg_slope , na.rm = T))
  #
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ac.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(ac_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="ac Slope")
  # stripchart(ac_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_tg.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(tg_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="tg Slope")
  # stripchart(ac_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ac_tg.png"), width = 1000, height = 1000)
  # plot(tg_slope ~ ac_slope, data=Gene_walks, pch=16)
  # dev.off()
  #
}

#----------------------------------------
# Create documentations for functions above
# devtools::document()

