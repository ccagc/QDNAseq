%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{segmentBins}
\alias{segmentBins}


\alias{segmentBins,QDNAseqReadCounts-method}

\title{Segments and calls total copy numbers}

\usage{
segmentBins(object, smoothBy=FALSE, alpha=0.0000000001, undo.splits="sdundo", undo.SD=1,
  normalize=TRUE, inter=c(-0.1, 0.1), force=FALSE, ...)
}

\description{
 Segments and calls total copy numbers.
}

\arguments{
  \item{object}{...}
  \item{smoothBy}{...}
  \item{alpha}{...}
  \item{undo.splits}{...}
  \item{undo.SD}{...}
  \item{normalize}{...}
  \item{inter}{...}
  \item{force}{...}
  \item{...}{Additional arguments passed to \code{\link[DNAcopy]{segment}}.}
}

\value{
  Returns ...
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[DNAcopy]{segment}} of the \pkg{DNAcopy} package,
  which implements the CBS method, is used to segment the data.
}

