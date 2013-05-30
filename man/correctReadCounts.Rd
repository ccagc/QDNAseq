%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  correctReadCounts.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{correctReadCounts}
\alias{correctReadCounts}


\title{Correct binned read counts for GC content and mappability}

\usage{
correctReadCounts(obj, span=0.65, family="symmetric", plotting=FALSE, ...)
}

\description{
 Correct binned read counts for GC content and mappability.
}

\arguments{
  \item{obj}{...}
  \item{span}{...}
  \item{family}{...}
  \item{plotting}{...}
  \item{...}{Not used.}
}

\value{
  Returns a named \code{\link[base]{list}} with elements ...
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[stats]{loess}} is used to fit the regression model.
}

