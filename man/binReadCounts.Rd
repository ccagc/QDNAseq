%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  binReadCounts.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{binReadCounts}
\alias{binReadCounts}


\title{Calculate binned read counts from a set of BAM files}

\usage{
binReadCounts(bins, bamfiles=NULL, path=".", ext="bam", bamnames=NULL, phenofile=NULL,
  genome="hg19", cache=TRUE, samtools="samtools", f="", F="0x0404", q=37, maxChunk=1e+08,
  allosomeBins="flag", incompleteBins="flag", blacklistedBins="flag")
}

\description{
 Calculate binned read counts from a set of BAM files.
}

\arguments{
  \item{bins}{...}
  \item{bamfiles}{...}
  \item{path}{...}
  \item{ext}{...}
  \item{bamnames}{...}
  \item{phenofile}{...}
  \item{genome}{...}
  \item{cache}{...}
  \item{samtools}{...}
  \item{f}{...}
  \item{F}{...}
  \item{q}{...}
  \item{maxChunk}{...}
  \item{allosomeBins}{...}
  \item{incompleteBins}{...}
  \item{blacklistedBins}{...}
}

\value{
  Returns a named \code{\link[base]{list}} with elements ...
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link{.binReadCountsPerSample}} is used.
}


\keyword{IO}
