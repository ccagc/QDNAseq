%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plot-methods.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{isobarPlot}
\alias{isobarPlot}


 \alias{isobarPlot,QDNAseqReadCounts,missing-method}

 \title{Plot median read counts as a function of GC content and mappability}

 \usage{
isobarPlot(x, y, ...)
}

 \description{
     Plot median read counts as a function of GC content and mappability.
 }

 \arguments{
     \item{x}{A \code{\link{QDNAseqReadCounts}} object.}
     \item{y}{missing}
     \item{...}{...}
%     \item{verbose}{If \code{\link[base:logical]{TRUE}}, verbose messages are produced.}
 }

 \examples{
 data(LGG150)
 readCounts <- LGG150
 isobarPlot(readCounts)
 }

 \author{Ilari Scheinin}


\keyword{hplot}
