\name{mrf}
\alias{mrf}
\title{
Fitting a one-dimensional Markov random field mixture model to ChIP-seq data.
}
\description{
\code{mrf} uses an MCMC algorithm to fit a one-dimensional Markov random field model for the latent binding profile from ChIP-seq data. The emission distribution of the  enriched state (signal) can be either Poisson or Negative Binomial (NB), while the emission distribution of the non-enriched state (background) can be either a Zero-inflated Poisson (ZIP) or a Zero-inflated Negative Binomial (ZINB). 
}
\usage{
mrf(data, method=NULL, exp.label = NULL, Niterations=10000, Nburnin=5000,  
     Poisprior=c(5, 1, 0.5, 1), NBprior=c(5, 1, 1, 1, 0.5, 1, 1, 1), 
     PoisNBprior=c(5,1,1,1, 0.5,1), var.NB=c(0.1, 0.1, 0.1, 0.1), parallel=TRUE)
}
\arguments{
 \item{data}{A list, whose first argument is a n x 3 matrix with information on the bins. The three columns should contain  "Chromosome", "Start" and "Stop" information. 
The second argument contains the counts of a single ChIP-seq experiment. This is a n x 1 matrix, where n is the number of bins. 
}
  \item{method}{
A character variable. Can be "Poisson", "PoisNB" or "NB" and it refers to the densities of the mixture distribution. "Poisson" means that a ZIP distribution is used for the background (with parameters pi and mean lambda_B) and a Poisson distribution for the signal (with parameter lambda_S); 
"PoisNB" means that a ZIP distribution is used for the bacground (with parameter pi and lambda_B) and a NB distribution for the signal (with mean mu_S and overdispersion phi_S); 
"NB" means that a ZINB distribution is used for the background (with parameters pi, mu_B and phi_B) and a  NB distribution for the signal (with mean mu_S and overdispersion phi_S).
}
  \item{exp.label}{
A charater vector, giving a label for experiment.
}
  \item{Niterations}{
An integer value, giving the number of MCMC iteration steps. Default value is 10000.
}
  \item{Nburnin}{An integer value, giving the number of burn-in steps. Default value is 5000.}
  \item{Poisprior}{The gamma priors for the parameter lambda in the Poisson-Poisson mixture: the first two elements are the priors for signal and the second two are priors for background. Default values are (5,1, 0.5, 1). }
  \item{NBprior}{The gamma priors for the mean mu and overdispersion parameter phi in the NB-NB mixture: the first two elements are the priors for mu_S for the signal; the third and fourth elements are priors for phi_S; the fifth and sixth elements are priors for mu_B for the background and the seventh and eighth are priors for phi_B. Default values are (5, 1, 1, 1, 0.5, 1, 1, 1).}
  \item{PoisNBprior}{The gamma priors for lambda_B and mu_S, phi_S in Poisson-NB mixture, the first two are priors for mu_S, the third and the fourth are priors for phi_S, the fifth and the sixth are priors for lambda_B. Default values are (5, 1,1,1, 0.5, 1). }
  \item{var.NB}{The variances used in the Metropolis-Hastling algorithm for estimating (mu_S, phi_S, mu_B, phi_B) for NB mixture or for estimating (mu_S, phi_S) for PoisNB mixture. 
  Default values are (0.1, 0.1, 0.1, 0.1) or (0.1, 0.1) for NB and PoisNB respectively. }
  \item{parallel}{A logical variable. If TRUE and the experiment has more than one chromosome, then the individual chromosomes will be processed in parallel, using the \code{clusterApplyLB} function in package \code{parallel}. Default value is TRUE.}
}
\value{
\item{data}{ The data provided as input.}
\item{parameters}{The estimates of parameters which are the mean of samples of parameters.}
\item{parameters.sample}{The samples matrix drawing from the posterior distributions of the parameters. 
The samples are collected one from every ten steps right after burn-in step. The column names for the matrix are (q_1, q_0, lambda_S, pi, lambda_B) if method="Poisson" or (q_1, q_0, mu_S, phi_S, pi, mu_B, phi_B) if method ="NB" or (q_1, q_0, mu_S, phi_S, pi, lambda_B) if method="PoisNB", where q_1 and q_0 are the transition probabilities that the current bin is enriched given the previous bin is enriched or not enriched, respectively. 
}
\item{PP}{The posterior probabilities that bins are enriched. }
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
CBPT0=list()
CBPT0$region=p300cbp.200bp$region
CBPT0$count=p300cbp.200bp$count[,1]

\dontrun{
##separate analysis of a single experiment using a ZINB-NB Markov random field model
CBPT0.mrf<-mrf(CBPT0, method="NB", exp.label="CBPT0")
}
}

