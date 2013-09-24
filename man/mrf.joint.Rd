\name{mrf.joint}
\alias{mrf.joint}
\title{
Joint fitting of a one-dimensional Markov random field model to multiple ChIP-seq datasets.
}
\description{
\code{mrf.joint} uses an MCMC algorithm to fit one-dimensional Markov random field models to multiple ChIP-seq datasets. 
These datasets could contain technical and biological replicates. If a single experiment is given, then the function \code{mrf} is used.
The emission distribution of the  enriched state (signal) could be either Poisson or Negative Binomial (NB), while the emission distribution of the non-enriched state (background) could be either a Zero-inflated Poisson (ZIP) or a Zero-inflated Negative Binomial (ZINB). 
}
\usage{
mrf.joint(data, method = NULL, rep.vec = NULL, p.vec = NULL, exp.label = NULL, 
      Niterations = 10000, Nburnin = 5000, Poisprior = NULL, NBprior = NULL, 
      PoisNBprior = NULL, var.NB = NULL, var.q=NULL, parallel=TRUE)
}
\arguments{
 \item{data}{
 A list, whose first argument is a n x 3 matrix with information on the regions. The three columns should contain  "Chromosome", "Start" and "Stop" information. 
The second list contains the counts of ChIP-seq experiments. This is a n x p matrix, where n is the number of regions and p is the number of experiments.
Count data for at least one experiment should be given.
}
  \item{method}{
  A character variable. Can be "Poisson", "PoisNB" or "NB" and it refers to the densities of the mixture distribution. "Poisson" means that a ZIP distribution is used for the background (with parameters pi and mean lambda_B) and a Poisson distribution for the signal (with parameter lambda_S); 
"PoisNB" means that a ZIP distribution is used for the bacground (with parameter pi and lambda_B) and a NB distribution for the signal (with mean mu_S and overdispersion phi_S) ; 
"NB" means that a ZINB distribution is used for the background (with parameters pi, mu_B and phi_B) and a  NB distribution for the signal (with mean mu_S and overdispersion phi_S).
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
At least one of rep.vec or p.vec should be given. For those experiments which do not share the same index (p.vec or rep.vec) with any other experiment, function \code{mrf} will be used. 
}
  \item{exp.label}{
A charater vector, giving a label for each experiments.
}
  \item{Niterations}{
An integer value, giving the number of MCMC iteration step. Default value is 10000.
}
  \item{Nburnin}{An integer value, giving the number of burn-in step. Default value is 5000.}
  \item{Poisprior}{The gamma priors for the parameter lambda in the Poisson-Poisson mixture: the first two elements are the priors for signal and the second two are priors for background. If p experiments are given, 
  then the prior should be a matrix 4 x p, where each column represents the priors for each experiment. Default values are (5,1, 0.5, 1) for each single experiment.}
  \item{NBprior}{The gamma priors for the mean mu and overdispersion parameter phi in the NB-NB mixture: the first two elements are the priors for mu_S for the signal; the third and fourth elements are priors for phi_S; the fifth and sixth elements are priors for mu_B for the background and the seventh and eighth are priors for phi_B. If p experiments are given then the prior should be a matrix 8 x p, where each column represents the priors for each experiment. Default values are (5, 1, 1, 1, 0.5, 1, 1, 1) for each single experiment. }
  \item{PoisNBprior}{The gamma priors for lambda_B and mu_S, phi_S in Poisson-NB mixture, the first two are priors for mu_S, the third and the fourth are priors for phi_S, the fifth and the sixth are priors for lambda_B. If more than one experiments are given then the prior should be a matrix 6 x p, each column represents the priors for each experiment. Default values are (5, 1,1,1, 0.5, 1) for each single experiment. }
  \item{var.NB}{The variances used in Metropolis-Hastling algorithm for estimates of (mu_S, phi_S, mu_B, phi_B) for NB mixture or for estimates of (mu_S, phi_S) for PoisNB mixture. If p experiments are given then var.NB should be 4 x p or 2 x p matrix for NB and PoisNB respectively, each column represents the variances used for each experiment. Default values for each single experiment are (0.1, 0.1, 0.1, 0.1) or (0.1, 0.1) for NB and PoisNB respectively. }
  \item{var.q}{the variances used in Metropolis-Hastling algorithm for estimates of q_0 and common ratio parameter when assume same p condition for multiple experiments. 
The number of components of var.q equals to number of experiment +1. Default values are 0.001 for each experiment and 0.005 for common ratio parameter. 
For example,  var.q=(0.001, 0.001, 0.005) for two experiments.} 
  \item{parallel}{A logical variable. If TRUE and the experiment has more than one chromosome, then the individual chromosomes will be processed in parallel, using the \code{clusterApplyLB} function in package \code{parallel}. Default value is TRUE.}
}
\value{
\item{data}{ The data provided as input.}
\item{parameters}{The list of parameters for each experiment, where each list contains the samples matrix drawing from the posterior distributions of the parameters. 
The samples are collected one from every ten steps right after burn-in step. The column names for the matrix are (q_1, q_0, lambda_S, pi, lambda_B) if method="Poisson" or (q_1, q_0, mu_S, phi_S, pi, mu_B, phi_B) 
if method ="NB" or (q_1, q_0, mu_S, phi_S, pi, lambda_B) if method="PoisNB", where q_1 and q_0 are the transition probabilities that the current region is enriched given the previous region is enriched or not enriched, respectively. 
}
\item{PP}{The list of posterior probabilities for each experiment, where each list contains a vector of posterior probabilities that regions are enriched.}
\item{rep.vec}{The rep.vec used for the analysis.}
\item{p.vec}{The p.vec used for the analysis.}
\item{method}{The method used for the analysis.}
\item{acrate.NB}{The acceptance rate of Metropolis-Hastling method.}
}
\references{
Bao et al. Joint modelling of ChIP-seq data via a Markov random field model, 2013.
}
\author{
Yanchun Bao and Veronica Vinciotti 
}
\seealso{
#See also \code{\link{mrf.joint}, \link{enrich.mrf}}
}
\examples{
data(p300cbp.200bp)
exp.label=c("CBPT0", "CBPT301", "CBPT302", "p300T0", 
    "p300T301", "p300T302", "WangCBP", "Wangp300")
CBPT30=list()
CBPT30$region=p300cbp.200bp$region[200001:210000,]
CBPT30$count=p300cbp.200bp$count[200001:210000,2:3]
\dontrun{
NBfit.simple<-mrf.joint(CBPT30, method="NB", rep.vec=c(1,1),
    p.vec=c(1,1), exp.label=exp.label[c(2,3)])

## Joint analysis combining technical replicates 
##   (CBPT301,CBPT302) and (p300T301, p300T302)
p300cbp.mrf<-mrf.joint(p300cbp.200bp, method="NB", 
    rep.vec=c(1,2,2,3,4,4,5,6), p.vec=c(1,2,2,3,4,4,5,6), exp.label=exp.label)
    
## Joint analysis, assuming the same p for two different proteins 
p300cbp.mrf.p<-mrf.joint(p300cbp.200bp, method="NB", 
    rep.vec=c(1,2,2,3,4,4,5,6), p.vec=c(1,2,2,1,2,2,3,3), exp.label=exp.label)
}
}

