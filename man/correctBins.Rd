%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  correctBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{correctBins}
\alias{correctBins}


\alias{correctBins,QDNAseqReadCounts-method}

\title{Correct binned read counts for GC content and mappability}

\usage{
correctBins(object, fit=NULL, method="ratio", adjustIncompletes=TRUE, ...)
}

\description{
    Correct binned read counts for GC content and mappability.
}

\arguments{
    \item{object}{An \code{\link{QDNAseqReadCounts}} object with \code{counts} data.}
    \item{fit}{An optional matrix of values to use for the correction. If
        NULL (default), assay data \code{fit} from object is used. If it is
        missing, it is generated with a call to \code{\link{estimateCorrection}}().}
    \item{method}{A character(1) string speficying the correction method.
        \code{ratio} (default) divides \code{counts} with \code{fit}.
        \code{median} calculates the median \code{fit}, and defines the
        correction for bins with GC content \code{gc} and mappability
        \code{map} as \code{median(fit) - fit(gc,map)}, which is added to
        \code{counts}. Method \code{none} leaves \code{counts} untouched.}
    \item{adjustIncompletes}{A boolean(1) specifying whether \code{counts} for
        bins with uncharacterized nucleotides (N's) in their reference genome
        sequence should be adjusted by dividing them with the percentage of
        characterized (A, C, G, T) nucleotides. Defaults to \code{\link[base:logical]{TRUE}}.}
    \item{...}{Additional arguments passed to \code{\link{estimateCorrection}}().}
}

\value{
    Returns a \code{\link{QDNAseqCopyNumbers}} object with assay data element
    \code{copynumber}.
}

\examples{
data(LGG150)
readCounts <- LGG150
readCountsFiltered <- applyFilters(readCounts)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
}

\author{Ilari Scheinin}

\seealso{
    Internally, \code{\link[stats]{loess}} is used to fit the regression model.
}


\keyword{manip}
\keyword{loess}
