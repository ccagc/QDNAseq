%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  makeCgh.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{makeCgh}
\alias{makeCgh}


\alias{makeCgh,QDNAseqCopyNumbers-method}

\title{Constructs a 'cghRaw', 'cghSeg', or 'cghCall' object}

\usage{
makeCgh(object, filter=TRUE, chromosomeReplacements=c(X = 23, Y = 24, MT = 25), ...)
}

\description{
    Constructs a 'cghRaw', 'cghSeg', or 'cghCall' object.
}

\arguments{
    \item{object}{A \code{\link{QDNAseqCopyNumbers}} object.}
    \item{filter}{If \code{\link[base:logical]{TRUE}}, bins are filtered, otherwise not.}
    \item{chromosomeReplacements}{A named integer vector of chromosome name
        replacements to be done. QDNAseq stores chromosome names as
        characters, but CGHcall expects them to be integers. Defaults to
        \code{c(X=23, Y=24, MT=25)} for human. Value of "auto" will use
        running numbers in order of appearance in the bin annotations.}
    \item{...}{Not used.}
}

\value{
    Returns a \code{\link[CGHbase]{cghRaw}} if the object has not been segmented,
    a \code{\link[CGHbase]{cghSeg}} if it has been segmented but not called,
    or \code{\link[CGHbase]{cghCall}} if it has been called as well.
}

\examples{
data(LGG150)
readCounts <- LGG150
readCountsFiltered <- applyFilters(readCounts)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
cgh <- makeCgh(copyNumbersSmooth)
}

\author{Ilari Scheinin}


\keyword{manip}
