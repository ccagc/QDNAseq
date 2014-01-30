%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeSegmentedBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeSegmentedBins}
\alias{normalizeSegmentedBins}


\alias{normalizeSegmentedBins,QDNAseqCopyNumbers-method}

\title{Normalize segmented bins}

\usage{
normalizeSegmentedBins(object, inter=c(-0.1, 0.1), force=FALSE)
}

\description{
 Normalize segmented bins.
}

\arguments{
  \item{object}{An object of class QDNAseqCopyNumbers.}
  \item{inter}{The interval in which the function should search for the
    normal level.}
  \item{force}{Whether to force execution when it causes removal of
    downstream calling results.}
}
\details{
  This function recursively searches for the interval containing the
  most segmented data, decreasing the interval length in each
  recursion.  The recursive search makes the post-segmentation
  normalization robust against local maxima. This function is
  particularly useful for profiles for which, after segmentation,
  the 0-level does not coincide with many segments. It is more or
  less harmless to other profiles. We advise to keep the search
  interval (inter) small, in particular at the positive (gain) side
  to avoid that the 0-level is set to a common gain level.
}
\value{
  Returns an object of class QDNAseqCopyNumbers with re-normalized data.
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[CGHcall]{postsegnormalize}} of the \pkg{CGHcall} package
  is used.
}

