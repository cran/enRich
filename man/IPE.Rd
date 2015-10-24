\name{IPE}
\alias{IPE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimating the ImmunoPrecipitation (IP) efficiency of a ChIP-seq experiment.
}
\description{
Calculate the IP efficiency of an experiment by using the mixture model parameters. 
}
\usage{
IPE(para, method = NULL)
}
\arguments{
  \item{para}{
A numeric vector. The parameters estimated by the mixture model. 
}
  \item{method}{
A charater variable. Can be "poisson" or "NB" and it refers to the densities of the mixture distribution.
}
}
\value{
IPE estimated value.
}
\references{
Bao et al. Accounting for immunoprecipitation efficiencies in the statistical analysis of ChIP-seq data. BMC Bioinformatics 2013, 14:169 DOI:10.1186/1471-2105-14-169.
}
\author{
Yanchun Bao and Veronica Vinciotti 
}
\seealso{
See also \code{\link{mix}, \link{mix.joint},\link{enrich.mix}}
}
\keyword{Statistics}

