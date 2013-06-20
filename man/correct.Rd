%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  correct.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{correct}
\alias{correct}


\alias{correct,QDNAseqReadCounts-method}

\title{Correct binned read counts for GC content and mappability}

\usage{
correct(object, span=0.65, family="symmetric", adjustIncompletes=TRUE, keepCounts=TRUE,
  storeResiduals=TRUE, ...)
}

\description{
 Correct binned read counts for GC content and mappability.
}

\arguments{
  \item{object}{...}
  \item{span}{...}
  \item{family}{...}
  \item{adjustIncompletes}{...}
  \item{keepCounts}{...}
  \item{storeResiduals}{...}
  \item{...}{Further aguments to loess.}
}

\value{
  Returns a named \code{\link[base]{list}} with elements ...
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[stats]{loess}} is used to fit the regression model.
}

