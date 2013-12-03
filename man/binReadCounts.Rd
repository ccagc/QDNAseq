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
binReadCounts(bins, bamfiles=NULL, path=NULL, ext="bam", bamnames=NULL, phenofile=NULL,
  cache=getOption("QDNAseq::cache", FALSE), force=!cache, isPaired=NA, isProperPair=NA,
  isUnmappedQuery=FALSE, hasUnmappedMate=NA, isMinusStrand=NA, isMateMinusStrand=NA,
  isFirstMateRead=NA, isSecondMateRead=NA, isNotPrimaryRead=NA,
  isNotPassingQualityControls=FALSE, isDuplicate=FALSE, minMapq=37)
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
  \item{cache=getOption("QDNAseq::cache", FALSE)}{...}
  \item{force=!cache}{...}
  \item{isPaired=NA}{...}
  \item{isProperPair=NA}{...}
  \item{isUnmappedQuery=FALSE}{...}
  \item{hasUnmappedMate=NA}{...}
  \item{isMinusStrand=NA}{...}
  \item{isMateMinusStrand=NA}{...}
  \item{isFirstMateRead=NA}{...}
  \item{isSecondMateRead=NA}{...}
  \item{isNotPrimaryRead=NA}{...}
  \item{isNotPassingQualityControls=NA}{...}
  \item{isDuplicate=FALSE}{...}
  \item{minMapq}{If quality scores exists, the minimum quality score required
    in order to keep a read, otherwise all reads are kept.}
}

\value{
  Returns a named \code{\link[base]{list}} with elements ...
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link{.binReadCountsPerSample}} is used.
}


\keyword{IO}
