\name{plot.minPtest}
\alias{plot.minPtest}
\title{
Plot method for a "minPtest" object
}
\description{
plot method for an object of class 'minPtest'. This function plots \code{(-log_{10})} transformed marginal p-values of each SNP in a basic scatterplot. In addition, horizontal lines of \code{(-log_{10})} transformed permutation-based p-values of each gene, covering all SNPs located on that gene, are plotted. The plot is indicated by two separated y-axes (\code{(-log_{10}(psnp))} at left hand side and \code{(-lambda*log_{10}(minp))} at the right hand side). After the Bonferroni correction depending on the \code{level}, significant genes and SNPs are by default highlighted in red, i.e. each gene line and each SNP point with permutation-based p-value and marginal p-value smaller than or equal to the \code{level}, respectively, is highlighted in red.
}
\usage{
\S3method{plot}{minPtest}(x, level=0.05, lambda=1, gene.name=FALSE,
     sigPch=pch, nonsigPch=pch, pch=20,
     sigLty=lty, nonsigLty=lty, lty=1,
     sigCol=col, nonsigCol=col, col=NULL, xlab="SNP", ...)
}
\arguments{
  \item{x}{
an object of class \code{minPtest}.
}
  \item{level}{
a numeric threshold that specifies which genes and SNPs are highlighted in the plot, i.e. the genes and SNPs with adjusted permutation-based p-values and marginal p-value, respectively, smaller than or equal to that threshold are by default highlighted in red. Default is 0.05.
}
  \item{lambda}{
a numeric value to scale the y-axis for the permutation-based p-values of the genes (indicated at the right hand side). Default is 1.
}
\item{gene.name}{
  a logical value, if TRUE, the gene names are shown at the
  x-axis. Default is FALSE.
}
\item{sigPch}{
Type of plotting for significant marginal p-values of the SNPs (if neither \code{sigPch} nor \code{pch} set: points)
}
\item{nonsigPch}{
Type of plotting for not significant marginal p-values of the SNPs (if neither \code{nonsigPch} nor \code{pch} set: points)
}
\item{pch}{
Set type of plotting for both \code{sigPch} and \code{nonsigPch} (but can be overwritten by \code{sigPch} and \code{nonsigPch} if set)
}
\item{sigLty}{
Type of plotting for significant permutation-based p-values of the genes (if neither \code{sigLty} nor \code{lty} set: solid lines)
}
\item{nonsigLty}{
Type of plotting for not significant permutation-based p-values of the genes (if neither \code{nonsigLty} nor \code{lty} set: solid lines)
}
\item{lty}{
Set type of plotting for both \code{sigLty} and \code{nonsigLty} (but can be overwritten by \code{sigLty} and \code{nonsigLty} if set)
}
\item{sigCol}{
Color for significant genes and SNPs (if neither \code{sigCol} nor \code{col} set: red)
}
\item{nonsigCol}{
Color for not significant genes and not significant SNPs (if neither \code{nonsigCol} nor \code{col} set: black)
}
\item{col}{
Set color for both \code{sigCol} and \code{nonsigCol} (but can be overwritten by \code{sigCol} and \code{nonsigCol} if set)
}
\item{xlab}{
xlab (Default: SNP)
}
\item{...}{
Further arguments for the plot function.
}
}
\details{
The plot method allows to compare the permutation-based p-values of each gene with marginal p-values of each SNP located on these genes. The y-axes are \code{(-log_{10})} transformed to obtain a disposition as a Manhattan plot for the points of the marginal p-values of the SNPs. The horizontal lines of the permutation-based p-values of each gene covering all SNPs located on this gene. Permutation-based p-values and marginal p-values, respectively, are by default highlighted in red, if they are smaller than or equal to \code{level} after Bonferroni correction.
}

\note{
The default for \code{gene.name=FALSE} should retained for performance reasons, if a large number of genes are included in the fit. No \code{ylim} should be specified as the plot is indicated by two separate y-axes.
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

### Plot significant permutation-based p-values of the genes and significant marginal p-values of the SNPs as blue dotted lines and blue points
## plot(minPtest.object,lambda=0.5,gene.name=TRUE,sigCol="blue",sigLty=2)
}
