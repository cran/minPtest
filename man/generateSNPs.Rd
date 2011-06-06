\name{generateSNPs}
\alias{generateSNPs}
\alias{print.generateSNPs}
\title{
Simulation of SNP data
}
\description{
Simulates SNP data with genotypes coded by 0, 1 and 2 as well as a binary and a continuous covariate, together with case-control status specified by logistic regression.
}
\usage{
generateSNPs(n, gene.no, block.no, block.size, p.same,
             p.different = NULL, p.minor, n.sample, SNPtoBETA)
}
\arguments{
  \item{n}{
an integer specifying the number of observations (cases and controls with 1:1 match) that should be generated. \code{n} should be an even number.
}
  \item{gene.no}{
an integer specifying the number of genes that should be generated.
}
  \item{block.no}{
an integer specifying the number of blocks per gene.
}
  \item{block.size}{
an integer specifying the number of SNPs per block.
}
  \item{p.same}{
either a numeric value specifying the probability for neighborhood SNPs within a block or a numeric vector of length \code{block.size}. In the latter case the argument \code{p.different} is ignored and has to be specified in the first item of \code{p.same}. The remaining items in the \code{p.same} vector specify the probabilities for neighborhood SNPs within the blocks, i.e. the probability that two neighboring alleles are equal within a block. If a numeric value, all SNPs, except the first item of each a block, will have the same neighborhood probability. If a vector of length \code{block.size}, each SNP of each block will have the neighborhood probability specified in the corresponding entry in \code{p.same}.
}
  \item{p.different}{
a numeric value specifying the probability for neighborhood blocks within a gene which is used if \code{p.same} is a scalar. The argument is ignored if \code{p.same} is a numeric vector and has to be specified in the first entry in \code{p.same}.
}
  \item{p.minor}{
a vector of length \code{block.no} containing the allele frequencies of the SNPs within a block. All SNPs in a block will have the same allele frequency.
}
  \item{n.sample}{
an integer specifying the number of simulated subjects from which the observations (case-control status) \code{n} are drawn.
}
  \item{SNPtoBETA}{
a matrix of non-negative numeric values of dimension \code{m * 2} consisting of the SNP index (first column) with \code{m <= snp.no} and the parameters (size of effect) of these SNPs (second column) for generating of case-control status.
}
}
\details{
\code{generateSNPs} generates a matrix consisting of \code{n} observations,
\code{snp.no=gene.no*block.no*block.size} SNPs with genotypes coded by 0, 1 and 2, two automatically generated covariates for adjustment or matching and the matchset numbers. The neighborhood probabilities for SNPs is given by \code{p.different} and/or \code{p.same} and the allele frequencies for SNPs is given by \code{p.minor}. The allele frequencies (\code{p.minor}) and the probabilities for neighborhood blocks (\code{p.different}) and/or \code{p.same}, respectively, can be different between the blocks on one gene but is repeated similar over all genes \code{gene.no}. The simulated SNP data structure is similar as in Schwender et al. (2011).

The response is determined by a logistic regression model given the SNPs, the binary covariate and the continuous covariate in the sim.cov matrix:

\code{P(Y=1|sim.cov)=exp(sim.cov*beta)/(1+exp(sim.cov*beta))}

Using the above model \code{P(Y=1|sim.cov)} is computed for each subject in \code{n.sample}, then the case and one control status for each of the \code{n.sample} subjects are determined by drawn randomly from a Bernoulli distribution using the probability \code{P(Y=1|sim.cov)}. From these \code{n.sample} subjects one case and one control observation is randomly drawn.  This algorithm is repeated \code{n/2} times for each randomly sampled value from the continuous covariate, i.e. one case and one control is randomly drawn from each of \code{n/2} times to generate the complete response vector of \code{length(n)}.

As output \code{generateSNPs} provides a response vector \code{y}, a SNP matrix \code{x}, a covariate matrix \code{cov} and a matchset vector \code{matchset} which can directly be used as input for the \code{\link{minPtest}}, see the example of the \code{\link{minPtest}} function.
}
\value{
An object of class 'generateSNPs', which is a list containing the following components:
\item{sim.data}{a matrix with \code{n} rows and \code{(snp.no+4)} columns containing response (case-control status) values, simulated SNP values, continuous matching covariate, binary matching covariate and matchset numbers.}
\item{y}{a numeric response vector coded with 0 (coding for controls) and 1 (coding for cases) of length \code{n}.}
\item{x}{a numeric \code{n * snp.no} matrix containing the simulated SNP data with genotypes coded by 0, 1 and 2.}
\item{cov}{a \code{n * 2} matrix containing the continuous matching covariate (likewise to age) and the binary matching covariate (likewise to gender).}
\item{matchset}{a numeric vector of length \code{n} containing the matching numbers (1:1 match).}
\item{snp.no}{number of SNPs in the simulated data set.}
\item{SNPtoGene}{the mapping matrix of dimension \code{p x 2} comprising of SNP names (first column) and the name of the genes (second column) on which the SNPs are located.}
\item{call}{call.}
}
\references{
Schwender, H. et al. (2011). Testing SNPs and sets of SNPs for importance in association studies. Biostatistics, 12, 18-32.
}
\author{
Stefanie Hieke \email{hieke@fdm.uni-freiburg.de}
}
\seealso{
\code{\link{minPtest}}
}
\examples{
# Generate a data set consisting of 100 subjects and 200 SNPs on 5 genes,
# with 4 blocks per gene with block size of 10, i.e. 10 SNPs per block
# yielding 40 SNPs per gene:

# specifying the matrix for 6 SNPs and corresponding parameters (effect size)
# for the generation of case-control status

SNP <- c(6,26,54,135,156,186)
BETA <- c(0.9,0.7,1.5,0.5,0.6,0.8)
SNPtoBETA <- matrix(c(SNP,BETA),ncol=2,nrow=6)
colnames(SNPtoBETA) <- c("SNP.item","SNP.beta")

set.seed(191)
sim1 <- generateSNPs(n=100,gene.no=5,block.no=4,block.size=10,p.same=0.9,
p.different=0.75,p.minor=c(0.1,0.4,0.1,0.4),n.sample=80,SNPtoBETA=SNPtoBETA)

# to reconstruct how to adopt the output from generateSNPs,
# see the example of the minPtest function.
}

\keyword{simulated SNP data}
