\name{p300cbp.200bp}
\alias{p300cbp.200bp}
\docType{data}
\title{
Example Data
}
\description{
\code{p300cbp.200bp} contains ChIP-seq counts in 200bp length bins on chromosome 21 for 8 experiments. The names of the 8 experiments are CBPT0, CBPT301, CBPT302, p300T0, p300T301, p300T302, WangCBP, Wangp300. The data are consistu
}
\usage{data(p300cbp.200bp)}
\format{
List of 2:\cr
 region:'data.frame', 234721 obs. of  3 variables, Chromosome, Start, Stop\cr  
 count: 'numeric', 234721 obs. of 8 variables,  CBPT0, CBPT301, CBPT302, p300T0, p300T301, p300T302, WangCBP, Wangp300
}
\source{
\describe{
\item{}{The first 6 datasets are from the GEO database, accession number GSE21026.}
\item{}{The last 2 data sets are from the GEO database, accession number GSE15735.}}
}
\references{
\describe{
\item{}{Romas, et.al, 2010. Genome-wide assessment of differential roles for p300 and CBP in transcription
regulation. Nucleic Acids Research, 38(16):5396-5408.}
\item{}{Wang, et.al, 2009. Genome-wide Mapping of HATs and
HDACs Reveals Distinct Functions in Active and Inactive Genes. Cell, 138:1019-1031.}
}
}
\examples{
data(p300cbp.200bp)
}
\keyword{datasets}
