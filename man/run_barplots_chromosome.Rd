% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnawalk_fft_functions.R
\name{run_barplots_chromosome}
\alias{run_barplots_chromosome}
\title{run_barplots_chromosome}
\usage{
run_barplots_chromosome(
  All_spec = All_spec,
  chromosome = chromosome,
  fname = fname,
  inpath = inpath,
  outpath = outpath,
  full_length = full_length,
  atflag = TRUE,
  chromnum = 0,
  numstd = 3,
  numstd2 = 1,
  numrange = 200,
  repeat_range = NULL,
  samplesize = 1,
  fftlength = 5000,
  binnum = 20
)
}
\arguments{
\item{nam}{input dataset}
}
\value{
output dataset
}
\description{
Makes barplots of minimums and maximums for repeats  and bedGraphs of repeat locations
}
\examples{
function()
}
