\name{enRich-package}
\alias{enRich-package}
\alias{enRich}
\docType{package}
\title{
Analysis of multiple ChIP-seq data.
}
\description{
enRich is an R package that performs a joint statistical modelling of ChIP-seq data, accounting for technical/biological replicates, multiple conditions and
the different IP efficiencies of individual experiments. 
}
\details{
\tabular{ll}{
Package: \tab enRich\cr
Type: \tab Package\cr
Version: \tab 3.1\cr
Date: \tab 2020-02-03\cr
Depends: \tab R(>= 3.2.2), parallel\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
mainfunctions: mix, mix.joint, mrf, mrf.joint, enrich.mix, enrich.mrf
}
\author{
Yanchun Bao \email{ybaoa@essex.ac.uk} and 
\cr Veronica Vinciotti \email{veronica.vinciotti@brunel.ac.uk}

Maintainer: Yanchun Bao \email{ybaoa@essex.ac.uk} 
}
\references{
Bao et al. Accounting for immunoprecipitation efficiencies in the statistical analysis of ChIP-seq data. BMC Bioinformatics 2013, 14:169 DOI:10.1186/1471-2105-14-169.\cr
Bao et al. Joint modelling of ChIP-seq data via a Markov random field model, Biostatistics 2014, 15(2):296-310 DOI:10.1093/biostatistics/kxt047.
}
\keyword{ package}
