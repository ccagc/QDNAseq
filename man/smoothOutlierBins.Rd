%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  smoothOutlierBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{smoothOutlierBins}
\alias{smoothOutlierBins}


 \alias{smoothOutlierBins,QDNAseqCopyNumbers-method}

 \title{Smooth outlier bins after normalization}

 \usage{
smoothOutlierBins(object, logTransform=TRUE, force=FALSE, ...)
}

 \description{
     Smooth outlier bins after normalization.
 }

 \arguments{
     \item{object}{A \code{\link{QDNAseqCopyNumbers}} object with \code{copynumber}
         data.}
     \item{logTransform}{If \code{\link[base:logical]{TRUE}} (default), data will be log2-transformed.}
     \item{force}{Running this function will remove possible segmentation and
         calling results. When they are present, running requires specifying
         \code{force} is \code{\link[base:logical]{TRUE}}.}
     \item{...}{Additional arguments passed to \code{\link[DNAcopy]{smooth.CNA}}.}
%     \item{verbose}{If \code{\link[base:logical]{TRUE}}, verbose messages are produced.}
 }

 \value{
     Returns a \code{\link{QDNAseqCopyNumbers}} object with the values for outliers
     smoothed. See \code{\link[DNAcopy]{smooth.CNA}} for more details. If
     \code{logTransform} is \code{\link[base:logical]{TRUE}}, these signals are log2-transformed prior
     to smoothing, but afterwards back-transformed..
 }

 \examples{
 data(LGG150)
 readCounts <- LGG150
 readCountsFiltered <- applyFilters(readCounts)
 readCountsFiltered <- estimateCorrection(readCountsFiltered)
 copyNumbers <- correctBins(readCountsFiltered)
 copyNumbersNormalized <- normalizeBins(copyNumbers)
 copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
 }

 \author{Ilari Scheinin}


\keyword{manip}
