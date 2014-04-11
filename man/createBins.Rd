%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  createBinAnnotations.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{createBins}
\alias{createBins}


\alias{calculateMappability}
\alias{calculateBlacklist}
\alias{iterateResiduals}

\title{Builds bin annotation data for a particular bin size}

\usage{
createBins(bsgenome, binSize, ignoreMitochondria=TRUE)
}

\description{
    Builds bin annotation data for a particular bin size.
}

\arguments{
    \item{bsgenome}{A BSgenome package.}
    \item{binSize}{A \code{\link[base]{numeric}} scalar specifying the width of the bins
        in units of kbp (1000 base pairs), e.g. \code{binSize=15} corresponds
        to 15 kbp bins.}
    \item{ignoreMitochondria}{Wheter to ignore the mitochondria.}
}

\value{
    Returns a \code{\link[base]{data.frame}} with columns \code{chromosome, start, end, bases,
    and gc}, which correspond to the chromosome name, positions of the first
    and last base pair in the bin, the percentage of characterized nucleotides
    (A, C, G, or T, i.e. non-N), and GC content (percentage of C and G
    nucleotides of non-N nucleotides).
}

\examples{
\dontrun{# NOTE: These take a very long time to run.
library(BSgenome.Hsapiens.UCSC.hg19)
bins <- createBins(BSgenome.Hsapiens.UCSC.hg19, 15)
bins$mappability <- calculateMappability(bins,
    bigWigFile='/path/to/wgEncodeCrgMapabilityAlign50mer.bigWig',
    bigWigAverageOverBed='/path/to/bigWigAverageOverBed')
bins$blacklist <- calculateBlacklist(bins,
    bedFiles=c('/path/to/wgEncodeDacMapabilityConsensusExcludable.bed',
    '/path/to/wgEncodeDukeMapabilityRegionsExcludable.bed'))
bins$residual <- iterateResiduals(readCountsG1K)
}
}

\author{Ilari Scheinin}

\seealso{
    \code{\link{downloadBinAnnotations}}().
}

