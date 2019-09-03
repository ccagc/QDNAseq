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
  force=FALSE, transformFun="log2", seed=FALSE, ...)
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
         sample. This has an effect when using parallel computing or specifying
         seeds for random number generation.}
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
     \item{seed}{If \code{\link[base:logical]{TRUE}}, random number generation (RNG) is fully
         reproducible. Normally, the segmentation process is
         split by sample, and provided seeds also used per sample. But when
         smoothing is performed (or in the the special case of
         \code{smoothBy=1L}), the process is split by chromosome, seeds used
         per chromosome, and results not necessarily reproducible across
         samples.}
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
     which implements the CBS method, is used to segment the data.
 }



\keyword{manip}
\keyword{smooth}
