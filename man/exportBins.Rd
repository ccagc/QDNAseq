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
exportBins(object, file, format=c("tsv", "igv", "bed", "vcf", "seg"),
  type=c("copynumber", "segments", "calls"), filter=TRUE, logTransform=TRUE, digits=3,
  chromosomeReplacements=c(`23` = "X", `24` = "Y", `25` = "MT"), ...)
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
    \item{logTransform}{If \code{\link[base:logical]{TRUE}} (default), exported data will be log2
        transformed for \code{format} in \code{"tsv"}, \code{"igv"}, and
        \code{"bed"}.  This argument is ignored if code{type = "calls"}.}
    \item{digits}{The number of digits to round to. If not \code{\link[base]{numeric}}, no
        no rounding is performed.}
    \item{chromosomeReplacements}{A named character vector of chromosome name
        replacements to be done. Only used when \code{object} is of class
        \code{\link{cghRaw}}, \code{\link{cghSeg}}, \code{\link{cghCall}}, or \code{\link{cghRegions}},
       since these classes store chromosome names as integers, whereas all
       QDNAseq object types use character vectors. Defaults to
        \code{c("23"="X", "24"="Y", "25"="MT")} for human.}
    \item{...}{Additional arguments passed to \code{\link[utils]{write.table}}.}
}

\value{
Returns the pathnames of the files written.
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
