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
segmentBins(object, smoothBy=FALSE, alpha=1e-10, undo.splits="sdundo", undo.SD=1,
  force=FALSE, transformFun="log2", ...)
}

 \description{
     Segments normalized copy number data.
 }

 \arguments{
     \item{object}{An object of class QDNAseqCopyNumbers.}
     \item{smoothBy}{An optional integer value to perform smoothing before
         segmentation by taking the mean of every smoothBy bins, and then
         segment those means. Default (\code{\link[base:logical]{FALSE}}) is to perform no smoothing.
         \code{smoothBy=1L} is a special case that will not perform smoothing,
         but will split the segmentation process by chromosome instead of by
         sample.}
     \item{alpha}{Significance levels for the test to accept change-points.
         Default is 1e-10.}
     \item{undo.splits}{A character string specifying how change-points are to
         be undone, if at all. Default is "sdundo", which undoes splits that
         are not at least this many SDs apart. Other choices are
         "prune", which uses a sum of squares criterion, and "none".}
     \item{undo.SD}{The number of SDs between means to keep a split if
         undo.splits="sdundo". Default is 1.0.}
     \item{force}{Whether to force execution when it causes removal of
         downstream calling results.}
     \item{transformFun}{A function to transform the data with. This can be
         the default "log2" for log2(x + .Machine$double.xmin),
         "sqrt" for the Anscombe transform of sqrt(x * 3/8) which
         stabilizes the variance, "none" for no transformation, or any
         R function that performs the desired transformation and also its
         inverse when called with parameter \code{inv=TRUE}.}
%     \item{segmentStatistic}{A character vector specifying which segment
%         statistic to use.}
%     \item{storeSegmentObjects}{A boolean to indicate whether to store the raw
%         DNAcopy objects within the QDNAseq objects. Segment objects can be
%         retrieved as segmentObject in assayData
%         eg. "assayDataElement(object, 'segmentObject')"}
     \item{...}{Additional arguments passed to \code{\link[DNAcopy]{segment}}.}
%     \item{verbose}{If \code{\link[base:logical]{TRUE}}, verbose messages are produced.}
 }

 \value{
     Returns an object of class QDNAseqCopyNumbers with segmentation results
         added.
 }

 \section{Numerical reproducibility}{
  This method make use of random number generation (RNG) via the
  \code{\link[DNAcopy]{segment}} used internally.  Because of this, calling the
  method with the same input data multiple times will each time give slightly
  different results.  To get numerically reproducible results, the random
  seed must be fixed, e.g. by using `set.seed()` at the top of the script.
 }

 \section{Parallel processing}{
   This function uses \pkg{future} to segment samples in parallel.
 }

 \examples{
 data(LGG150)
 readCounts <- LGG150
 readCountsFiltered <- applyFilters(readCounts)
 readCountsFiltered <- estimateCorrection(readCountsFiltered)
 copyNumbers <- correctBins(readCountsFiltered)
 copyNumbersNormalized <- normalizeBins(copyNumbers)
 copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
 copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
 }

 \author{Ilari Scheinin}

 \seealso{
     Internally, \code{\link[DNAcopy]{segment}} of the \pkg{DNAcopy} package,
     which implements the CBS method [1,2], is used to segment the data.
 }

 \references{
  [1] A.B. Olshen, E.S. Venkatraman (aka Venkatraman E. Seshan), R. Lucito
      and M. Wigler, \emph{Circular binary segmentation for the analysis of
      array-based DNA copy number data}, Biostatistics, 2004 \cr
  [2] E.S. Venkatraman and A.B. Olshen, \emph{A faster circular binary
      segmentation algorithm for the analysis of array CGH data},
      Bioinformatics, 2007 \cr
 }



\keyword{manip}
\keyword{smooth}
