% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{tmm_matrix}
\alias{tmm_matrix}
\title{Use edgeR to do TMM normalization with or without full library size}
\usage{
tmm_matrix(counts, lib.size = NULL, group = c(rep(1, 3), rep(2, 3)))
}
\arguments{
\item{lib.size}{Numeric vector representing the full library size for all samples.
If it is null, use the sum of read count in each sample as the library size.}

\item{group}{vector or factor giving the experimental group/condition for each sample/library.
If it is null,no control/condition are used.}

\item{rawcount}{Matrix representing the raw read counts for samples. Row are the peaks, column are the samples}
}
\value{
normalized read counts \code{(log2(read counts / sum of read counts + 1))} .
}
\description{
Use edgeR to do TMM normalization with or without full library size
}
\examples{
cpmD = cpmNorm(countAll)

# The robustness.
cpmD = cpmNorm(countAll, libsize=fd$lsIP)
cpmD = cpmNorm(countAll)
}
