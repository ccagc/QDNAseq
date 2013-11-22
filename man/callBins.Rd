%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  callBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{callBins}
\alias{callBins}


\alias{callBins,QDNAseqReadCounts-method}

\title{Call aberrations from segmented copy number data}

\usage{
callBins(object, ...)
}

\description{
 Call aberrations from segmented copy number data.
}

\arguments{
  \item{object}{An object of class QDNAseqReadCounts}
  \item{...}{Additional arguments passed to \code{\link[CGHcall]{CGHcall}}.}
}

\value{
  Returns an object of class QDNAseqReadCounts with segmentation results
    added.
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[CGHcall]{CGHcall}} and \code{\link[CGHcall]{ExpandCGHcall}} of
    the \pkg{DNAcopy} package are used.
}

