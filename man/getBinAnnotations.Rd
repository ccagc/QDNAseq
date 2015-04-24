%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  getBinAnnotations.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getBinAnnotations}
\alias{getBinAnnotations}

\title{Gets bin annotation data for a particular bin size}

\usage{
getBinAnnotations(binSize, genome="hg19", type="SR50", force=FALSE,
  path=getOption("QDNAseq::binAnnotationPath", "http://qdnaseq.s3.amazonaws.com"))
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
       or "PE100".}
    \item{force}{If \code{\link[base:logical]{TRUE}}, the bin anonnation data is retrieved/calculated
       regardless of it already exists in the cache or not.}
    \item{path}{A \code{\link[base]{character}} string specifying the path for the bin
        annotation files. Defaults to downloading from the Internet, but can
        also be a local path. Can also be defined by setting the
        \code{QDNAseq::binAnnotationPath} option.}
}

\details{
    Gets bin annotation data for a particular bin size
}

\value{
    Returns a \code{\link[Biobase]{AnnotatedDataFrame}} object.
}

\examples{
\dontrun{
bins <- getBinAnnotations(15)
}
}

\author{Ilari Scheinin}

\seealso{
    \code{\link{createBins}}().
}


\keyword{IO}
