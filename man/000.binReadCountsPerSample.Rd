%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  binReadCounts.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{.binReadCountsPerSample}
\alias{.binReadCountsPerSample}


\title{Calculate binned read counts from a BAM file}

\usage{
.binReadCountsPerSample(bins, bamfile, cache=TRUE, maxChunk=100000000L, isPaired=NA,
  isProperPair=NA, isUnmappedQuery=FALSE, hasUnmappedMate=NA, isMinusStrand=NA,
  isMateMinusStrand=NA, isFirstMateRead=NA, isSecondMateRead=NA, isNotPrimaryRead=NA,
  isNotPassingQualityControls=FALSE, isDuplicate=FALSE, minMapq=37)
}

\description{
 Calculate binned read counts from a BAM file.
}

\arguments{
  \item{bins}{...}
  \item{bamfile}{...}
  \item{cache}{...}
  \item{maxChunk}{...}
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
  \item{minMapq=37}{...}
}

\value{
  Returns ...
}

\author{Ilari Scheinin}

\seealso{
  To retrieve bin read counts for a set of BAM files,
  see \code{\link{binReadCounts}}().
}



\keyword{IO}
\keyword{internal}
