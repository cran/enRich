\name{mix.joint}
\alias{mix.joint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint fitting of mixture of Poisson or NB densities to ChIP-seq experiments.
}
\description{
\code{mix.joint} uses an EM algorithm to jointly fit ChIP-seq data for two or more experiments. 
Technical replicates are accounted for in the model as well as individual ChIP efficiencies for each experiment. 
Prior biological knowledge, such as the expectation of a similar number of binding profiles for the same protein
under two similar conditions, can also be included in the model to aid robustness in the detection of enriched and differentially bound regions. 
The output of \code{mix.joint} can be further analysed by \code{enrich.mix}. 
}
\usage{
mix.joint(data, method = NULL, para.sep = NULL, rep.vec = NULL, 
    p.vec = NULL, exp.label = NULL, stopdiff = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A list, whose first argument is a n x 3 matrix with information on the regions. The three columns should contain  "Chromosome", "Start" and "Stop" information. 
The second list contains the counts of ChIP-seq experiments. This is a n x p matrix, where n is the number of regions and p is the number of experiments.
Count data for at least one experiment should be given.
}
 \item{method}{
A character variable. Can be "Poisson" or "NB" and it refers to the densities of the mixture distribution.
}
\item{para.sep}{
A p x q matrix, where p is the number of experiments and q is the number of parameters in the mixture model. 
This is used as initial parameters for the joint modelling function.  We recommend using the parameters of \code{mix}, as these are optimized for each experiment. 
If there are no technical replicates, then the parameters of the mix function are automatically used for the mix.joint output. 
}
  \item{rep.vec}{
A non-zero integer vector. The vector of replicate indices, of length equal to the number of experiments. Technical replicates share the same index, 
e.g c(1,2,2,3,4,4,5,6) for 8 experiments where the 2nd and 3rd are two technical replicates and similarly the 5th and 6th.
}
  \item{p.vec}{
A non-zero integer vector. The vector of p indices, with p=P(X_s=1) for any region s. This vector is of length equal to the number of experiments. 
Experiments with the same probability of enrichment share the same p index, such as technical replicates and/or proteins with
a similar number of binding sites, e.g. c(1,1,1,2,3,3,3,4) if the first three experiments have the same p and similarly the 5th, 6th and 7th experiments.
This allows to propertly account for the different IP efficiencies in the joint analysis. 
At least one of rep.vec or p.vec should be given. For those experiments which do not share the same index (p.vec or rep.vec) with any other experiment, a single mixture model will be fitted. 
}
\item{exp.label}{
A charater vector, giving the labels for each experiment. 
}
  \item{stopdiff}{
A numeric variable. A prescribed small quantity for determining the convergence of the EM algorithm. Default value is 1e-04. 
}
}

\value{
\item{data}{ The data provided as input.}
\item{parameters}{The parameters estimated by the mixture model. The parameters are (p, lambda_S, lambda_B, k) when method="poisson" or (p, mu_S, phi_S, mu_B, phi_B, k) when method="NB". 
p is the proportion of signal in the mixture model. For a Poisson mixture model, lambda_S and lambda_B represent the mean of the signal and mean of the background, respectively. 
For a NB mixture model, mu_S and phi_S are the mean and overdispersion of the signal density, respectively, whereas mu_B and phi_B are the mean and overdispersion of the background density, respectively.}
\item{rep.vec}{The rep.vec used for the analysis}
\item{p.vec}{The p.vec used for the analysis}
\item{method}{The method used for the analysis}
}
\references{
Bao et al. Accounting for immunoprecipitation efficiencies in the statistical analysis of ChIP-seq data. BMC Bioinformatics 2013, 14:169 DOI:10.1186/1471-2105-14-169.
}
\author{
Yanchun Bao and Veronica Vinciotti 
}

\seealso{
See also \code{\link{mix}, \link{enrich.mix}}
}
\examples{
data(p300cbp.1000bp)
exp.label=c("CBPT0", "CBPT301", "CBPT302", "p300T0", 
    "p300T301", "p300T302", "WangCBP", "Wangp300")
## Simple examples -- only two experiments and first 5000 observations
CBPT30=list()
CBPT30$region=p300cbp.1000bp$region[1:5000,]
CBPT30$count=p300cbp.1000bp$count[1:5000,2:3]
Poisfit.simple<-mix(CBPT30, method="Poisson", exp.label=exp.label[c(2,3)])
## Joint analysis combining technical replicates 
##   (CBPT301,CBPT302)
Poisfit.joint<-mix.joint(CBPT30, Poisfit.simple$parameters, method="Poisson", 
    rep.vec=c(1,1), p.vec=c(1,1), exp.label=exp.label[c(2,3)]) 

\dontrun{
##separate analysis of single experiments using a Poisson latent mixture model for all 8 data sets
Poisfit.separate<-mix(p300cbp.1000bp, method="Poisson", exp.label=exp.label)

## Joint analysis combining technical replicates 
##   (CBPT301,CBPT302 and p300T301,p300T302)
Poisfit.joint1<-mix.joint(p300cbp.1000bp, Poisfit.separate$parameters, method="Poisson", 
    rep.vec=c(1,2,2,3,4,4,5,6), p.vec=c(1,2,2,3,4,4,5,6), 
    exp.label=exp.label) 

## Joint analysis, assuming the same p for two different proteins 
##   at the same time point
Poisfit.joint2<-mix.joint(p300cbp.1000bp, Poisfit.separate$parameters, method="Poisson", 
    rep.vec=c(1,2,2,3,4,4,5,6), p.vec=c(1,2,2,1,2,2,5,5), exp.label=exp.label) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Model }% __ONLY ONE__ keyword per line
