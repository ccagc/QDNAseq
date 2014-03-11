%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  compareToReference.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{compareToReference}
\alias{compareToReference}


\alias{compareToReference,QDNAseqCopyNumbers,numeric-method}

\title{Divide binned read counts with those of reference samples}

\usage{
compareToReference(object, references, force=FALSE)
}

\description{
 Divide binned read counts with those of reference samples.
}

\arguments{
  \item{object}{An object of class \code{\link{QDNAseqCopyNumbers}}.}
  \item{references}{A numeric vector of indexes of the reference sample. Must
    be the same length as there are samples in object. When \code{\link[base]{NA}}, the sample
    will be kept as is. When \code{\link[base:logical]{FALSE}}, the sample will be removed from the
    output. As an example, object contains three samples: tumor1, tumor2, and
    normal2. There is no reference for tumor1, but normal2 is a matched normal
    sample from the same patient as tumor2. The keep tumor1 as is, but to
    divide tumor2 with normal2, argument references should be \code{c(NA, 3,
    FALSE)}.}
  \item{force}{Whether to force the operation even when downstream data will
    be lost.}
}

\value{
  Returns a \code{\link{QDNAseqCopyNumbers}} object with the desired samples divided
  by the signal of their reference samples.
}

\examples{
  data(LGG150)
  readCounts <- LGG150
  readCountsFiltered <- applyFilters(readCounts)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  # Note: the following command will compare the sample to itself, which
  # does not really make sense:
  tumorVsNormal <- compareToReference(copyNumbersSmooth, c(1))
}

\author{Ilari Scheinin}

