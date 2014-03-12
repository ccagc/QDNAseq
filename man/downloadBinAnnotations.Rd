%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  downloadBinAnnotations.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{downloadBinAnnotations}
\alias{downloadBinAnnotations}


\title{Gets bin annotation data for a particular bin size}

\usage{
downloadBinAnnotations(binSize, genome="hg19", type="SR50", force=FALSE)
}

\description{
    Gets bin annotation data for a particular bin size.
}

\arguments{
    \item{binSize}{A \code{\link[base]{numeric}} scalar specifying the width of the bins
       in units of kbp (1000 base pairs), e.g. \code{binSize=15} corresponds
       to 15 kbp bins.}
    \item{genome}{A \code{\link[base]{character}} string specify the genome and genome version
       to be used.}
    \item{type}{A \code{\link[base]{character}} string specify the experiment type, e.g. "SR50"
       or "PE1000".}
    \item{force}{If \code{\link[base:logical]{TRUE}}, the bin anonnation data is retrieved/calculated
       regardless of it already exists in the cache or not.}
}

\value{
    Returns a \code{\link[Biobase]{AnnotatedDataFrame}} object.
}

\examples{
\dontrun{
bins <- downloadBinAnnoattions(15)
}
}

\author{Ilari Scheinin}

\seealso{
    \code{\link{createBins}}().
}


\keyword{IO}
