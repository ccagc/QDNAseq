%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{segmentBins}
\alias{segmentBins}


\alias{segmentBins,QDNAseqCopyNumbers-method}

\title{Segments normalized copy number data}

\usage{
segmentBins(object, smoothBy=FALSE, alpha=0.0000000001, undo.splits="sdundo", undo.SD=1,
  force=FALSE, ...)
}

\description{
 Segments normalized copy number data.
}

\arguments{
  \item{object}{An object of class QDNAseqCopyNumbers.}
  \item{smoothBy}{An optional integer value to perform smoothing before
    segmentation by taking the mean of every smoothBy bins, and then segment
    those means. Default is to perform no smoothing.}
  \item{alpha}{Significance levels for the test to accept change-points.}
  \item{undo.splits}{A character string specifying how change-points are to be
    undone, if at all.  Default is "none".  Other choices are
    "prune", which uses a sum of squares criterion, and "sdundo",
    which undoes splits that are not at least this many SDs
    apart.}
  \item{undo.SD}{The number of SDs between means to keep a split if
    undo.splits="sdundo".}
  \item{force}{Whether to force execution when it causes removal of
    downstream calling results.}
  \item{...}{Additional arguments passed to \code{\link[DNAcopy]{segment}}.}
}

\value{
  Returns an object of class QDNAseqCopyNumbers with segmentation results
    added.
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[DNAcopy]{segment}} of the \pkg{DNAcopy} package,
  which implements the CBS method, is used to segment the data.
}

