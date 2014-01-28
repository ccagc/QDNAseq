%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QDNAseqReadCounts-accessors.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{QDNAseqReadCounts}
\docType{class}
\alias{QDNAseqReadCounts}
\alias{QDNAseqReadCounts-class}


\title{Container for QDNAseq read count data}

\description{
  Container for QDNAseq read count data
}

\section{Assay data elements}{
  An object of this class contains (a subset) the following elements:
  \describe{
    \item{\code{counts}}{(\code{\link[base]{numeric}}) Binned uncorrected read counts
      as non-negative integers in \eqn{\{0,1,2,...\}}.
      An object with only this field is created by \code{\link{binReadCounts}}().}
    \item{\code{fit}}{(\code{\link[base]{numeric}}; optional) Loess fit of "count" signals
      as non-negative doubles in \eqn{[0,+\infty)}.
      This element is added after calling \code{\link{estimateCorrection}}().}
  }
}

\section{Missing values}{
  The bin data (assay data) may contain missing values.
}


\author{Ilari Scheinin}
\keyword{classes}
