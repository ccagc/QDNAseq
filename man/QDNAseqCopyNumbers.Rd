%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QDNAseqCopyNumbers-accessors.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{QDNAseqCopyNumbers}
\docType{class}
\alias{QDNAseqCopyNumbers}
\alias{QDNAseqCopyNumbers-class}


\title{Container for QDNAseq read count data}

\description{
  Container for QDNAseq read count data
}

\section{Assay data elements}{
  An object of this class contains the following elements:
  \describe{
    \item{\code{copynumber}}{(\code{\link[base]{numeric}}) Corrected
      "count" signals in \eqn{[0,+\infty)} An object with only this field is
      created by \code{\link{correctBins}}().}
    \item{\code{segmented}}{(\code{\link[base]{numeric}}; optional) Segmented data in
      \eqn{[0,+\infty)}, added by calling \code{\link{segmentBins}}().}
    \item{\code{calls}}{(\code{\link[base]{integer}}; optional) Calls as -2=deletion, -1=loss,
      0=normal, 1=gain, 2=amplification, added by calling \code{\link{callBins}}().}
    \item{\code{probdloss}}{(\code{\link[base]{numeric}}; optional) Probabilities of deletions in
      \eqn{[0,1]}, added by calling \code{\link{callBins}}().}
    \item{\code{probloss}}{(\code{\link[base]{numeric}}; optional) Probabilities of losses in
      \eqn{[0,1]}, added by calling \code{\link{callBins}}().}
    \item{\code{probnorm}}{(\code{\link[base]{numeric}}; optional) Probabilities of normal copy
      number in \eqn{[0,1]}, added by calling \code{\link{callBins}}().}
    \item{\code{probgain}}{(\code{\link[base]{numeric}}; optional) Probabilities of gains in
      \eqn{[0,1]}, added by calling \code{\link{callBins}}().}
    \item{\code{probamp}}{(\code{\link[base]{numeric}}; optional) Probabilities of amplifications
      in \eqn{[0,1]}, added by calling \code{\link{callBins}}().}
  }
}

\section{Missing values}{
  The bin data (assay data) may contain missing values.
}

\author{Ilari Scheinin}
\keyword{classes}
