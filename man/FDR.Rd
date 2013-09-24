\name{FDR}
\alias{FDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identification of enriched regions by controlling a given FDR level 
}
\description{
Identify enriched regions by controlling false discovery rate at a specified level. 
}
\usage{
FDR(prob0, cr = 0.05)
}
\arguments{
  \item{prob0}{
A numeric vector. The probability of the null hypothesis being true (i.e. a region not-enriched). 
}
  \item{cr}{
A numeric variable. The level of FDR for identifying the enriched regions.
}
}
\value{
\item{X}{The index of enriched regions (1: enriched, 0: not-enriched).}
}
\references{
Bao et al. Accounting for immunoprecipitation efficiencies in the statistical analysis of ChIP-seq data. BMC Bioinformatics 2013, 14:169 doi:10.1186/1471-2105-14-169.
}
\author{
Yanchun Bao and Veronica Vinciotti 
}
\seealso{
See also \code{\link{mix}, \link{mix.joint}, \link{enrich.mix}}
}
\keyword{ Statistics}
