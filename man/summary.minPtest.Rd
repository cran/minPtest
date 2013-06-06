\name{summary.minPtest}
\alias{summary.minPtest}
\alias{print.summary.minPtest}
\title{
Summary method for a "minPtest" object
}
\description{
Summary method for objects of class "minPtest"
}
\usage{
\S3method{summary}{minPtest}(object, level = 0.05, sign.SNP = FALSE, ...)
}
\arguments{
  \item{object}{
an object of class \code{minPtest}, i.e. the output of a \code{\link{minPtest}} call.
}
  \item{level}{
a numeric threshold that specifies which genes are shown in the summary, i.e. the genes with adjusted permutation-based p-values smaller than or equal to that threshold are printed. Default is 0.05.
}
\item{sign.SNP}{
  a logical value; if TRUE, print, in addition to the genes selected by a \code{level}, the SNPs with adjusted marginal p-values smaller than or equal to the \code{level} (same as for genes) located on these genes. Default is FALSE, all SNP located on these genes, selected according to the \code{level}, are shown. Default is FALSE.
}
\item{...}{Further arguments for the \code{summary} method. Not used.}
}
\details{
Prints the genes with adjusted permutation-based p-value smaller than or equal to a \code{level}, the corresponding permutation-based p-values, the adjusted permutation-based p-values and the SNPs located on these genes, either all SNPs or SNPs selected by the \code{level}, sorted by the adjusted marginal p-values, with marginal p-values and adjusted marginal p-values.
}
\value{
\code{summary.minPtest} returns a list. Each item characterizes a gene, selected according to a \code{level}, list items are named by means of these genes. Each gene item contains a list of data frames, a data frame for the permutation-based p-values and adjusted permutation-based p-values for this gene and a data frame for the marginal p-values and adjusted marginal p-values for the SNPs located on that gene, either all SNPs or SNPs selected by the \code{level}.
}
\author{
Stefanie Hieke \email{hieke@imbi.uni-freiburg.de}
}
\seealso{
\code{\link{generateSNPs}}, \code{\link{minPtest}}
}
\examples{
## Continuing the example from minPtest and generateSNPs:
# generate a data set consisting of 100 subjects and 200 SNPs on 5 genes.

SNP <- c(6,26,54,135,156,186)
BETA <- c(0.9,0.7,1.5,0.5,0.6,0.8)
SNPtoBETA <- matrix(c(SNP,BETA),ncol=2,nrow=6)
colnames(SNPtoBETA) <- c("SNP.item","SNP.beta")

set.seed(191)
sim1 <- generateSNPs(n=100,gene.no=5,block.no=4,block.size=10,p.same=0.9,
	             p.different=0.75,p.minor=c(0.1,0.4,0.1,0.4),
                     n.sample=80,SNPtoBETA=SNPtoBETA)

# Cochran Armitage Trend Test without covariates and default permutations.
# Example: Run R sequential

### Seed
set.seed(10)
seed1 <- sample(1:1e7,size=1000)
###
minPtest.object <- minPtest(y=sim1$y, x=sim1$x, SNPtoGene=sim1$SNPtoGene,
	                    seed=seed1)
###
summary(minPtest.object)
}
