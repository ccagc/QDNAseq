%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  addPhenodata.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{addPhenodata}
\alias{addPhenodata}


\title{Adds phenotype data from a file to a QDNAseqReadCounts or a
    QDNAseqCopyNumbers object}

\usage{
addPhenodata(object, phenofile)
}

\description{
    Adds phenotype data from a file to a QDNAseqReadCounts or a
    QDNAseqCopyNumbers object.
}

\arguments{
    \item{object}{A \code{\link{QDNAseqReadCounts}} or \code{\link{QDNAseqCopyNumbers}}
        object.}
    \item{phenofile}{A file name with phenotypic data for samples in
        \code{object}.}
}

\value{
    Returns a \code{\link{QDNAseqReadCounts}} or \code{\link{QDNAseqCopyNumbers}} object
    with phenotype data added.
}

\examples{
data(LGG150)
readCounts <- LGG150
\dontrun{
readCounts <- addPhenodata(readCounts, "phenodata.txt")
}
}
\author{Ilari Scheinin}



\keyword{IO}
\keyword{file}
