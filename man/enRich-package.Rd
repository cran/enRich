\name{enRich-package}
\alias{enRich-package}
\alias{enRich}
\docType{package}
\title{
An R package for the analysis of multiple ChIP-seq data.
}
\description{
enRich is an R package that performs a joint statistical modelling of ChIP-seq data, accounting for technical/biological replicates, multiple conditions and
the different IP efficiencies of individual experiments. 
}
\details{
\tabular{ll}{
Package: \tab enRich\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2013-09-24\cr
Depends: \tab R(>= 2.15.3), parallel\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
mainfunctions: mix, mix.joint, mrf, mrf.joint, enrich.mix, enrich.mrf
}
\author{
Yanchun Bao \email{yanchun.bao@brunel.ac.uk} and 
\cr Veronica Vinciotti \email{veronica.vinciotti@brunel.ac.uk}

Maintainer: Yanchun Bao \email{yanchun.bao@brunel.ac.uk} 
}
\references{
Bao et al. Accounting for immunoprecipitation efficiencies in the statistical analysis of ChIP-seq data. BMC Bioinformatics 2013, 14:169 doi:10.1186/1471-2105-14-169.\cr
Bao et al. Joint modelling of ChIP-seq data via a Markov random field model, 2013.
}
\keyword{ package}


