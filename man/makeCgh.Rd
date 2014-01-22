%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  makeCgh.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{makeCgh}
\alias{makeCgh}


\alias{makeCgh,QDNAseqReadCounts-method}

\title{Constructs a 'cghRaw', 'cghSeg', or 'cghCall' object}

\usage{
makeCgh(object)
}

\description{
 Constructs a 'cghRaw', 'cghSeg', or 'cghCall' object.
}

\arguments{
  \item{object}{A \code{\link{QDNAseqReadCounts}} object.}
}

\value{
  Returns a \code{\link[CGHbase]{cghRaw}} if the object has not been segmented,
  a \code{\link[CGHbase]{cghSeg}} if it has been segmented but not called,
  or \code{\link[CGHbase]{cghCall}} if it has been called as well.
}

\author{Ilari Scheinin}

