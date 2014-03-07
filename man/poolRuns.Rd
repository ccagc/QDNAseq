%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  poolRuns.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{poolRuns}
\alias{poolRuns}


\alias{poolRuns,QDNAseqSignals,character-method}

\title{Pools binned read counts across samples}

\usage{
poolRuns(object, samples, force=FALSE)
}

\description{
 Pools binned read counts across samples.
}

\arguments{
  \item{object}{A \code{\link{QDNAseqReadCounts}} or \code{\link{QDNAseqCopyNumbers}}
    object.}
  \item{samples}{A character vector of new sample names. Samples with
    identical names will be pooled together. Must be the same length as there
    are samples in \code{object}.}
  \item{force}{Whether to force the operation even when downstream data will
    be lost.}
}

\value{
  Returns a \code{\link{QDNAseqReadCounts}} or \code{\link{QDNAseqCopyNumbers}} object.
}

\author{Ilari Scheinin}

