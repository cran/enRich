\name{enrich.mrf}
\alias{enrich.mrf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Detection of enriched and differentially bound regions for fitting results of \code{mrf} and \code{mrf.joint}.}
\description{
\code{enrich.mrf} returns the enriched regions or differentially bound regions using the \code{mrf} or the \code{mrf.joint} model, by controlling a given FDR level. 
\code{enrich.mrf} also calculates the IP efficiencies for each experiment.  
}
\usage{
enrich.mrf(object, analysis = "joint", differential = FALSE, 
    diff.vec = NULL, cr = 0.05, crdiff = 0.05)
}
\arguments{
  \item{object}{
The output of \code{mrf} if analysis="separate" or of \code{mrf.joint} if analysis="joint". 
}
  \item{analysis}{
A character variable. Default value is "joint" and the object should be the output of \code{mrf.joint}. If analysis="separate", then the object should be the output of \code{mrf}. 
}
  \item{differential}{
A logical variable. If TRUE, the function will compute the posterior probability of differential binding of any two experiments or two conditions, as specified by \code{diff.vec}. Default value is FALSE.
}
  \item{diff.vec}{
A numeric vector. If differential = TRUE, diff.vec must be given to show which experiments are to be used in the comparison. At the moment, this is restricted to two conditions
(e.g. two proteins at the same time point), so the value for diff.vec should be only 0, 1, 2, where 0 indicates which experiments are not to be used in the analysis, 
1 and 2 stand for conditions 1 and 2, respectively. diff.vec should be of the same length as the number of experiments in object.
}
  \item{cr}{
A numeric variable. The level of FDR for identifying the enriched regions.
}
  \item{crdiff}{
A numeric variable. The level of FDR for identifying the differentially bound regions.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{enrich}{The list of enriched regions for each condition at the chosen FDR. Note that there is only one list of enriched regions for replicates, if a joint model is used.}
\item{diffenrich1}{The list of regions bound only by condition 1.}
\item{diffenrich2}{The list of regions bound only by condition 2.}
\item{X}{A n x p matrix of enrichment for each region and each condition, at the given FDR cutoff (1: enriched, 0: not-enriched).}
\item{diffprob1}{A n-dimensional vector of posterior probabilities of differential binding for the two conditions under study; diffprob0=1-diffprob1.}
\item{diffX1}{A n-dimensional index of regions bound only by condition 1 (0: not bound, 1: bound).}
\item{diffX2}{A n-dimensional index of regions bound only by condition 2.}
\item{IPE}{A p-dimensional vector of estimated IP efficiency values for each experiment.}
}
\references{
Bao et al. Joint modelling of ChIP-seq data via a Markov random field model, Biostatistics 2014, 15(2):296-310 DOI:10.1093/biostatistics/kxt047.
}
\author{
Yanchun Bao and Veronica Vinciotti
}
\seealso{
See also \code{\link{mrf}, \link{mrf.joint}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistics}

