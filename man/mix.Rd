\name{mix}
\alias{mix}
%\link{mix.joint}, \link{enrich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting mixture of two densities, either Poisson or Negative Binomial, to ChIP-seq data.
}
\description{
\code{mix} uses an EM algorithm to fit ChIP-seq count data by a latent mixture model with two components. One component is the signal density and the other is the background density.
\code{mix} can deal with more than one experiment at the same time. In this case, it fits individual models to each experiment. The output of this function can be used for further analysis by \code{mix.joint} or \code{enrich.mix}.  
}
\usage{
mix(data, method = NULL, initialpara=NULL, fixoffset=FALSE, fixk=3,krange=c(0:10), 
     exp.label=NULL, stopdiff=1e-04, parallel=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A list, whose first argument is a n x 3 matrix with information on the bins. The three columns should contain  "Chromosome", "Start" and "Stop" information. 
The second list contains the counts of ChIP-seq experiments. This is a n x p matrix, where n is the number of bins and p is the number of experiments.
Count data for at least one experiment should be given.
}
  \item{method}{
A character variable. Can be "Poisson" or "NB" and it refers to the densities of the mixture distribution.
}
\item{initialpara}{
A numeric matrix or vector. The initial parameters given for EM algorithm. In form of c("p", "lambda_S", "lambda_B") if method="Poisson" or c("p", "mu_S", "phi_S", "mu_B", "phi_B") if method="NB". 
Could be a matrix if initial values are the different for multiple experiments or a vector if initial values are the same. 
If not given, then a default value of (0.1, 10, 1) or (0.1, 10, 1, 1, 1) for method="Poisson" or "NB" respectively. 
}
  \item{fixoffset}{
A logical variable. If TRUE, the offset of the signal distribution is fixed by the user and is the same for all experiments. If FALSE, the offset is estimated empirically for each experiment. Default value is FALSE.
}
  \item{fixk}{
A numeric variable. The value of the offset, when fixoffset = TRUE.
}
  \item{krange}{
A numeric vector. The range of the offset, when fixoffset = FALSE. Default range is from 0 to 10. 
}  
\item{exp.label}{
A charater vector, giving a label for each experiments.
%If exp.label=NULL, then each experiment will be represented by Experiment and the order of that Experiment in data. 
}
  \item{stopdiff}{
A numeric variable. A prescribed small quantity for determining the convergence of the EM algorithm. Default value is 1e-04.  
}
\item{parallel}{
A logical variable. If TRUE, then the individual experiments will be processed in parallel, using the \code{clusterApplyLB} function in package \code{parallel}. Default value is TRUE.
}
}
\value{
\item{data}{ The data provided as input.}
\item{parameters}{The parameters estimated by the mixture model. The parameters are (p, lambda_S, lambda_B, k) when method="Poisson" or (p, mu_S, phi_S, mu_B, phi_B, k) when method="NB". 
p is the proportion of signal in the mixture model. For a Poisson mixture model, lambda_S and lambda_B represent the mean of the signal and mean of the background, respectively. 
For a NB mixture model, mu_S and phi_S are the mean and overdispersion of the signal density, respectively, whereas mu_B and phi_B are the mean and overdispersion of the background density, respectively.}
\item{method}{The method used for the analysis}
}
\references{
Bao et al. Accounting for immunoprecipitation efficiencies in the statistical analysis of ChIP-seq data. BMC Bioinformatics 2013, 14:169 DOI:10.1186/1471-2105-14-169.
}
\author{
Yanchun Bao and Veronica Vinciotti 
}
\seealso{
See also \code{\link{mix.joint}, \link{enrich.mix}}
}
\examples{
data(p300cbp.1000bp)
exp.label=c("CBPT0", "CBPT301", "CBPT302", "p300T0", 
"p300T301", "p300T302", "WangCBP", "Wangp300")
## Simple examples -- only two experiments and first 5000 observations
CBPT30=list()
CBPT30$region=p300cbp.1000bp$region[1:5000,]
CBPT30$count=p300cbp.1000bp$count[1:5000,2:3]
Poissonfit.simple<-mix(CBPT30, method="Poisson", exp.label=exp.label[c(2,3)])

\dontrun{
##separate analysis of single experiments using a Poisson latent mixture model for all 8 data sets
Poisfit.separate<-mix(p300cbp.1000bp, method="Poisson", exp.label=exp.label)
##separate analysis of single experiments using a NB latent mixture model for all 8 data sets
NB.separate<-mix(p300cbp.1000bp, method="NB", exp.label=exp.label)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Model }% __ONLY ONE__ keyword per line
