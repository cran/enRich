\name{enrich.mix}
\alias{enrich.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Detection of enriched and differentially bound regions for fitting results of \code{mix} and \code{mix.joint}.}
\description{
\code{enrich.mix} returns the enriched regions or differentially bound regions using the \code{mix} or the \code{mix.joint} model, by controlling a given FDR level. 
\code{enrich.mix} also calculates the IP efficiencies for each experiment.  
}
\usage{
enrich.mix(object, analysis = "joint", differential = FALSE, 
    diff.vec = NULL, cr = 0.05, crdiff = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
The output of \code{mix} if analysis="separate" or of \code{mix.joint} if analysis="joint". 
}
  \item{analysis}{
A character variable. Default value is "joint" and the object should be the output of \code{mix.joint}. If analysis="separate", then the object should be the output of \code{mix}. 
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
\item{enrich}{The list of enriched regions for each condition at the chosen FDR. Note that there is only one list of enriched regions for replicates, if a joint model is used.}
\item{diffenrich1}{The list of regions bound only by condition 1.}
\item{diffenrich2}{The list of regions bound only by condition 2.}
\item{ppx1}{A n x p matrix of posterior probabilities of enrichment for each region and each condition. ppx0=1-ppx1.}
\item{X}{A n x p matrix of enrichment for each region and each condition, at the given FDR cutoff (1: enriched, 0: not-enriched).}
\item{diffprob1}{A n-dimensional vector of posterior probabilities of differential binding for the two conditions under study; diffprob0=1-diffprob1.}
\item{diffX1}{A n-dimensional index of regions bound only by condition 1 (0: not bound, 1: bound).}
\item{diffX2}{A n-dimensional index of regions bound only by condition 2.}
\item{IPE}{A p-dimensional vector of estimated IP efficiency values for each experiment.}
}
\references{
Bao et al. Accounting for immunoprecipitation efficiencies in the statistical analysis of ChIP-seq data. BMC Bioinformatics 2013, 14:169 DOI:10.1186/1471-2105-14-169.
}
\author{
Yanchun Bao and Veronica Vinciotti
}
\seealso{
See also \code{\link{mix}, \link{mix.joint}}
}
\examples{
tempdir()
data(p300cbp.1000bp)
exp.label=c("CBPT0", "CBPT301", "CBPT302", "p300T0", 
    "p300T301", "p300T302", "WangCBP", "Wangp300")

## Simple examples -- only two experiments and first 5000 observations
CBPT30=list()
CBPT30$region=p300cbp.1000bp$region[1:5000,]
CBPT30$count=p300cbp.1000bp$count[1:5000,2:3]
Poisfit.simple<-mix(CBPT30, method="Poisson", exp.label=exp.label[c(2,3)])
enrich.mix.simple<-enrich.mix(Poisfit.simple, analysis="separate")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistics}

