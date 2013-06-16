%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalize.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{normalize}
\alias{normalize}


 \alias{normalize,QDNAseqReadCounts-method}
 \alias{normalize,cghRaw-method}

 \title{Normalizes binned read counts}

 \usage{
normalize(object, method="median", smoothOutliers=TRUE, logTransform=TRUE, ...)
}

 \description{
  Normalizes binned read counts.
 }

 \arguments{
   \item{object}{...}
   \item{method}{A \code{\link[base]{character}} string specifying ...}
   \item{logTransform}{If \code{\link[base:logical]{TRUE}}, ..., otherwise, ...}
   \item{smoothOutliers}{If \code{\link[base:logical]{TRUE}}, ..., otherwise, ...}
   \item{...}{Further arguments to DNAcopy::smooth.CNA}
 }

 \value{
   Returns a named \code{\link[base]{list}} containing elements ...
 }

 \author{Ilari Scheinin}

 \seealso{
   Internally, ...
 }

