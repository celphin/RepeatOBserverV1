% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnawalk_fft_functions.R
\name{run_chromosomelistNEW}
\alias{run_chromosomelistNEW}
\title{run_chromosomelistNEW}
\usage{
run_chromosomelistNEW(
  nam = nam,
  fname = fname,
  inpath = inpath,
  outpath = outpath,
  startgroup = 1,
  endgroup = NULL,
  majorlen = 5e+06,
  length_majorgroup = 1e+06,
  name_majorgroup = "1500Kbp",
  length_submajor = 5e+05,
  name_submajorgroup = "500Kbp",
  length_fftgroup = 5000,
  name_fftgroup = "5Kbp",
  length_minor = 25000,
  submajor_nam = "25Kb",
  full_length = NULL
)
}
\arguments{
\item{nam}{input dataset}
}
\value{
output dataset
}
\description{
rebuilds twolists if already run
}
\examples{
function()
}
