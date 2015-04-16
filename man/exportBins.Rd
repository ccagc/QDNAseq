%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  exportBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{exportBins}
\alias{exportBins}


\alias{exportBins,QDNAseqSignals-method}

\title{Exports to a file}

\usage{
exportBins(object, file, format=c("tsv", "igv", "bed"), type=c("copynumber", "segments",
  "calls"), filter=TRUE, logTransform=TRUE, digits=3, ...)
}

\description{
    Exports to a file.
}

\arguments{
    \item{object}{A \code{\link{QDNAseqReadCounts}} or \code{\link{QDNAseqCopyNumbers}}
        object.}
    \item{file}{Filename. For formats that support only one sample per file,
        such as BED, '\%s' can be used as a placeholder for sample name or
        '\%d' for sample number.}
    \item{format}{Format to export in. Currently supported ones are "tsv" (tab
        separated values), "igv" (Integrative Genomics Viewer), and "bed" (BED
        file format).}
    \item{type}{Type of data to export, options are "copynumber" (corrected or
        uncorrected read counts), "segments", or "calls".}
    \item{filter}{If \code{\link[base:logical]{TRUE}}, bins are filtered, otherwise not.}
    \item{logTransform}{If \code{\link[base:logical]{TRUE}} (default), data will be log2-transformed.}
    \item{digits}{The number of digits to round to. If not \code{\link[base]{numeric}}, no
        no rounding is performed.}
    \item{...}{Additional arguments passed to \code{\link[utils]{write.table}}.}
}

\details{
    Exports \code{object} to a file.
}

\examples{
\dontrun{
data(LGG150)
readCounts <- LGG150
readCountsFiltered <- applyFilters(readCounts)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
exportBins(copyNumbersSmooth, file="LGG150.igv", format="igv")
}
}

\author{Ilari Scheinin}



\keyword{IO}
\keyword{file}
