%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentDataWithWeights.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{segmentDataWithWeights}
\alias{segmentDataWithWeights}


\title{Segments and calls total copy numbers}

\usage{
segmentDataWithWeights(input, weights=TRUE, ...)
}

\description{
 Segments and calls total copy numbers.
}

\arguments{
  \item{input}{...}
  \item{weights}{Either \code{\link[base:logical]{TRUE}} or a vector of weights. If \code{\link[base:logical]{TRUE}},
    the 1000 Genomes residuals are used as weights.}
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

