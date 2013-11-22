%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  correctBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{correctBins}
\alias{correctBins}


\alias{correctBins,QDNAseqReadCounts-method}

\title{Correct binned read counts for GC content and mappability}

\usage{
correctBins(object, span=0.65, family="symmetric", adjustIncompletes=TRUE,
  type="relative", keepCounts=TRUE, storeResiduals=TRUE, storeLoess=FALSE, force=FALSE,
  ...)
}

\description{
 Correct binned read counts for GC content and mappability.
}

\arguments{
  \item{object}{...}
  \item{span}{...}
  \item{family}{...}
  \item{adjustIncompletes}{...}
  \item{type}{...}
  \item{keepCounts}{...}
  \item{storeResiduals}{...}
  \item{storeLoess}{...}
  \item{...}{Further aguments to loess.}
  \item{force}{...}
}

\value{
  Returns a named \code{\link[base]{list}} with elements ...
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[stats]{loess}} is used to fit the regression model.
}

