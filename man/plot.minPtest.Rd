\name{plot.minPtest}
\alias{plot.minPtest}
\title{
Plot method for a "minPtest" object
}
\description{
plot method for an object of class 'minPtest'. This function plots \code{(-log_{10})} transformed marginal p-values of each SNP in a basic scatterplot. In addition, horizontal lines of \code{(-log_{10})} transformed permutation-based p-values of each gene, covering all SNPs located on that gene, are plotted. The plot is indicated by two separated y-axes (\code{(-log_{10}(psnp))} at left hand side and \code{(-lambda*log_{10}(minp))} at the right hand side). After the Bonferroni correction depending on the \code{level}, significant genes and SNPs are highlighted in red, i.e. each gene line and each SNP point with permutation-based p-value and marginal p-value smaller than or equal to the \code{level}, respectively, is highlighted in red.
}
\usage{
\S3method{plot}{minPtest}(x, level = 0.05, lambda = 1, gene.name = FALSE, ...)
}
\arguments{
  \item{x}{
an object of class \code{minPtest}.
}
  \item{level}{
a numeric threshold that specifies which genes and SNPs are highlighted in the plot, i.e. the genes and SNPs with adjusted permutation-based p-values and marginal p-value, respectively, smaller than or equal to that threshold are highlighted in red. Default is 0.05.
}
  \item{lambda}{
a numeric value to scale the y-axis for the permutation-based p-values of the genes (indicated at the right hand side). Default is 1.
}
\item{gene.name}{
  a logical value, if TRUE, the gene names are shown at the
  x-axis. Default is FALSE.
}
\item{...}{Further arguments for the plot function.}
}
\details{
The plot method allows to compare the permutation-based p-values of each gene with marginal p-values of each SNP located on these genes. The y-axes are \code{(-log_{10})} transformed to obtain a disposition as a Manhattan plot for the points of the marginal p-values of the SNPs. The horizontal lines of the permutation-based p-values of each gene covering all SNPs located on this gene. Permutation-based p-values and marginal p-values, respectively, are highlighted in red, if they are smaller than or equal to \code{level} after Bonferroni correction.
}

\note{
The default for \code{gene.name=FALSE} should retained for performance reasons, if a large number of genes are included in the fit.
}

\value{
No value returned
}

\author{
Stefanie Hieke \email{hieke@fdm.uni-freiburg.de}
}
\seealso{
  \code{\link{minPtest}}, \code{\link{generateSNPs}}
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
	p.different=0.75,p.minor=c(0.1,0.4,0.1,0.4),n.sample=80,SNPtoBETA=SNPtoBETA)


# Cochran Armitage Trend Test without covariates and default permutations.
# Example: Run R sequential

### Seed
set.seed(10)
seed1 <- sample(1:1e7,size=1000)
###
minPtest.object <- minPtest(y=sim1$y, x=sim1$x, SNPtoGene=sim1$SNPtoGene,
	seed=seed1)
###
plot(minPtest.object,lambda=0.5,gene.name=TRUE)
}
