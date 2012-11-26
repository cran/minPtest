\name{plot.minPtest}
\alias{plot.minPtest}
\title{
Plot method for "minPtest" object
}
\description{
plot method for an object of class 'minPtest'. Plots allowing to get an impression of important genes or/and SNPs.
}
\usage{
\S3method{plot}{minPtest}(x, type=c("gene" ,"SNP" ,"both"), level=0.05, lambda=1, gene.name=FALSE,
     sigPch=pch, nonsigPch=pch, pch=20,
     sigLty=lty, nonsigLty=lty, lty=1,
     sigCol=col, nonsigCol=col, col=NULL, xlab, ...)
}
\arguments{
  \item{x}{
an object of class \code{minPtest}.
}
\item{type}{
by default, permutation-based p-values for each gene are plotted (\code{"gene"}). \code{"SNP"}: marginal p-values for each SNP are plotted. \code{"both"}: marginal p-values for each SNP and the transformed permutation-based p-values for each gene are displayed in a combined plot, see Details.
}
  \item{level}{
a numeric threshold that specifies which genes or/and SNPs are highlighted in the plot. I.e. not depending on the used \code{type} argument, the genes or/and SNPs with adjusted permutation-based p-values or/and marginal p-value, respectively, which are smaller than or equal to that threshold are by default highlighted in red. Default is 0.05.
}
  \item{lambda}{
only useful for \code{type="both"}. A numeric value to scale the y-axis for the permutation-based p-values of the genes (indicated at the right hand side). Default is 1.
}
\item{gene.name}{
  only useful for \code{type="SNP"} and \code{type="both"}. A logical value, if TRUE, the gene names are shown at the
  x-axis. Default is FALSE.
}
\item{sigPch}{
Type of plotting for significant permutation-based p-values (\code{type="gene"}) or for significant marginal p-values (\code{type="SNP"} and \code{type="both"}) (if neither \code{sigPch} nor \code{pch} set: points)
}
\item{nonsigPch}{
Type of plotting for non significant permutation-based p-values (\code{type="gene"}) or non significant marginal p-values (\code{type="SNP"} and \code{type="both"}) (if neither \code{nonsigPch} nor \code{pch} set: points)
}
\item{pch}{
Set type of plotting for both \code{sigPch} and \code{nonsigPch} (but can be overwritten by \code{sigPch} and \code{nonsigPch} if set)
}
\item{sigLty}{
only used for \code{type="both"}. Type of plotting for significant permutation-based p-values (if neither \code{sigLty} nor \code{lty} set: solid lines)
}
\item{nonsigLty}{
only used for \code{type="both"}. Type of plotting for non significant permutation-based p-values (if neither \code{nonsigLty} nor \code{lty} set: solid lines)
}
\item{lty}{
only used for \code{type="both"}. Set type of plotting for both \code{sigLty} and \code{nonsigLty} (but can be overwritten by \code{sigLty} and \code{nonsigLty} if set)
}
\item{sigCol}{
Color for significant genes or/and significant SNPs (if neither \code{sigCol} nor \code{col} set: red)
}
\item{nonsigCol}{
Color for non significant genes or/and non significant SNPs (if neither \code{nonsigCol} nor \code{col} set: black)
}
\item{col}{
Set color for both \code{sigCol} and \code{nonsigCol} (but can be overwritten by \code{sigCol} and \code{nonsigCol} if set)
}
\item{xlab}{
xlab (Default: Gene if \code{type="both"}) and SNP if \code{type=SNP} or \code{type="both"}, respectively
}
\item{...}{
Further arguments for the plot function.
}
}
\details{
The function plots either \code{(-log_{10})} transformed permutation-based p-values for each gene or \code{(-log_{10})} transformed marginal p-values for each SNP in a basic scatterplot. The y-axis is \code{(-log_{10})} transformed to obtain a disposition as a Manhattan plot for the points of the marginal p-values of the SNPs. Furthermore, an alternative given by the function is to display the marginal p-values for each SNP and the transformed permutation-based p-values for each gene in a combined plot. The \code{(-log_{10})} transformed marginal p-values for each SNP are plotted as points. In addition, horizontal lines of \code{(-lambda*log_{10})} transformed permutation-based p-values of each gene, covering all SNPs located on that gene, are plotted. The composed plot is indicated by two separated y-axes (\code{(-log_{10}(psnp))} at left hand side and \code{(-lambda*log_{10}(minp))} at the right hand side). After correction for multiple hypothesis testing depending on the \code{level} and the argument \code{adj.method} in the \code{\link{minPtest}} function, but not depending on the used type of plot, significant genes and SNPs are by default highlighted in red, i.e. each permutation-based p-value or/and marginal p-value smaller than or equal to the \code{level}, respectively, is highlighted in red.
}

\note{
The default for \code{gene.name=FALSE}, used for \code{type="SNP"} and \code{type="both"}, should kept for performance reasons, if a large number of genes are included in the fit. For \code{type="both"} no \code{ylim} should be specified as the plot is indicated by two separate y-axes.
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
### Combined plot for permutation-based p-values and marginal p-values.
plot(minPtest.object,type="both",lambda=0.5,gene.name=TRUE)

### Combined plot for permutation-based p-values and marginal p-values. Plot significant permutation-based p-values and significant marginal p-values as blue dotted lines and blue points
## plot(minPtest.object,type="both",lambda=0.5,gene.name=TRUE,sigCol="blue",sigLty=2)
}
