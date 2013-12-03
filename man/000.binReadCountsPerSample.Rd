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
.binReadCountsPerSample(bins, bamfile, cache, force, isPaired, isProperPair,
  isUnmappedQuery, hasUnmappedMate, isMinusStrand, isMateMinusStrand, isFirstMateRead,
  isSecondMateRead, isNotPrimaryRead, isNotPassingQualityControls, isDuplicate, minMapq)
}

\description{
 Calculate binned read counts from a BAM file.
}

\arguments{
  \item{bins}{...}
  \item{bamfile}{...}
  \item{cache}{...}
  \item{force}{...}
  \item{isPaired}{...}
  \item{isProperPair}{...}
  \item{isUnmappedQuery}{...}
  \item{hasUnmappedMate}{...}
  \item{isMinusStrand}{...}
  \item{isMateMinusStrand}{...}
  \item{isFirstMateRead}{...}
  \item{isSecondMateRead}{...}
  \item{isNotPrimaryRead}{...}
  \item{isNotPassingQualityControls}{...}
  \item{isDuplicate}{...}
  \item{minMapq}{...}
}

\value{
  Returns an \code{\link[base]{integer}} \code{\link[base]{vector}}.
}

\author{Ilari Scheinin}

\seealso{
  To retrieve bin read counts for a set of BAM files,
  see \code{\link{binReadCounts}}().
}



\keyword{IO}
\keyword{internal}
