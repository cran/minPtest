\name{minPtest}
\alias{minPtest}
\alias{print.minPtest}
\title{
A gene region-level testing procedure for each candidate gene based on resampling using the min P test
}
\description{
Permutation-based p-values estimation via min P test, a gene region-level summary for each candidate gene. The gene region-level summary assesses the smallest p-trend within each gene region comparing cases and controls. The min P test is permutation-based method that can be based on different univariate tests per SNP. Inference is based on the permutation distribution of the ordered p-values from the marginal tests of each SNP. Potentially accelerated by parallelization, if a compute cluster or a multicore computer is available.
}
\usage{
minPtest(y, x, SNPtoGene, formula = NULL, cov = NULL, matchset = NULL,
         permutation = 1000, seed = NULL, subset = NULL,
         multicore = FALSE, parallel = FALSE, trace = FALSE,
         aggregation.fun = min, ...)
}
\arguments{
  \item{y}{
a numeric response vector coded with 0 (coding for controls) and 1 (coding for cases) of length \code{n}.
}
  \item{x}{
a numeric \code{n * p} matrix of covariates (i.e. SNPs) containing the genotypes coded by 0, 1 and 2. Thus, each column is assumed to represent one of the SNPs with corresponding column names. Detail for SNP coding are given below.
}
  \item{SNPtoGene}{
a mapping matrix of dimension \code{p x 2} comprising of SNP names (first column) same as the column names of x, and the names of the genes (second column) on which the SNPs are located.
}
  \item{formula}{
(optional) for unconditional or conditional logistic regression, respectively, with or without covariates other than SNPs. A symbolic description of the model to be fitted including covariates only, see \code{\link{glm}} or \code{\link[Epi]{clogistic}}, the latter requires library \pkg{Epi}, else the default method Cochran Armitage Trend Test, which requires library \pkg{scrime}, is fitted. Details of model specification are given below.
}
  \item{cov}{
(optional) a \code{n * q} matrix containing the covariates for adjustment with corresponding column names.
}
  \item{matchset}{
(optional) a numeric vector of length \code{n} containing matching numbers, needed for conditional logistic regression.
}
  \item{permutation}{
number of permutations employed for obtaining a null distribution.
}
  \item{seed}{
(optional) vector of length \code{permutation}. Allows reproducibility even when running in parallel and for different numbers of parallel processes.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{multicore}{
indicates whether computation in the permuted data sets should be performed in parallel, using package \pkg{multicore}. If TRUE, the default numbers of cores are used. A value larger than 1 is taken to be the number of cores. Note that, \pkg{multicore} is currently not available with Microsoft Windows.
}
  \item{parallel}{
logical value indicates whether computation should be performed in parallel on a compute cluster, using package \pkg{snowfall}. \code{parallel=FALSE} sequential execution is performed. If TRUE the initialization function of this package, \code{sfInit()}, should be called before calling minPtest. See Details.
}
  \item{trace}{
logical value indicating whether progress in estimation should be indicated by printing the number of permutation that is currently used. (ignored if running in parallel via \pkg{snowfall}).
}
\item{aggregation.fun}{
  function that is used to combine the trend p-values over multiple loci within a gene region. By default the minimum  (\code{"min"}) is applied to obtain candidate gene region-level summaries. Any other function to integrate the p-values into one single test statistic can be used, e.g., \code{median} or different functions designed by the user.
}
\item{...}{
Further arguments for \code{aggregation.fun}. In case of \code{NA/NaN} within the evaluated marginal trend p-values for the SNPs from the original data (\code{psnp}) or within the permuted trend p-values in the permutation samples (\code{psnpperm}), \code{na.rm=TRUE} has to be specified.
}
}
\details{
The idea of the gene region-level summary, using the min P test procedure (Westfall and Young, 1993; Westfall et al., 2002; Chen et al., 2006), is to identify candidate genes by assessing the statistical significance of the smallest p-trend from a set of SNPs (single nucleotide polymorphisms) within each gene region comparing cases and controls by permutation-based resampling methods. A SNP, i.e. genetic variation occurring at a single base pair position in the genome, can take three possible values (genotypes); either there is no SNP variant in comparison to some reference coding (homozygous reference (0)) or the SNP variant occurs on one of the two base pair positions (heterozygous (1)), or both base pairs have a variant comparing to the reference coding. \pkg{minPtest} permits to include, instead of the genotypes 0, 1 and 2, also combined carrier SNPs, e.g. coding 0 and 1 (1 + 2).

Computation of min P test is based on marginal trend p-values for a set of univariate SNP disease association and trend p-values for permutation samples for each SNP. The \pkg{minPtest} package brings together three different kinds of tests to compute such p-values that are scattered over several R packages, and automatically selects the one most appropriate for the design at hand. In any case a response vector \code{y}, a SNP matrix \code{x} and a mapping matrix \code{SNPtoGene} are required. Then the default, a Cochran Armitage Trend Test (Cochran, 1954; Armitage, 1955), is automatically fitted to compute p-values. The Cochran Armitage Trend Test does not depend on covariates and matching scenario. Additionally adding a formula, see also \code{\link{glm}} from package \pkg{base}, and a covariate matrix \code{cov} an unconditional logistic regression is fitted. Unconditional logistic regression can be used without or with covariates for adjustment; either \code{formula=y~1} or \code{formula=y~cov1+cov2+...}. The former does not need any information relative to covariates and matching scenario. However, the latter is general for frequency matching with the inclusion of matching variables for adjustment specified in the covariate matrix \code{cov}. Providing a matchset, as in the case of 1:1; 1:2 etc. matching, and a formula, see also \code{\link[Epi]{clogistic}} from package \pkg{Epi}, a conditional logistic regression is fitted. Conditional logistic regression can be used without or with covariates for adjustment; either \code{formula=y~1} or \code{formula=y~cov1+cov2...}. In the latter case covariates other than matching variables can be used and have to be specified in the covariate matrix \code{cov}. In general, there are two possibilities to specify the formula, first if no covariates are used for adjustment, the formula has to be written as \code{y~1} without specifying the covariate matrix \code{cov}. Second if covariates other than SNPs are used for adjustment, the formula has to be written as response vector \code{y} on the left of a \code{~} operator, and the clinical covariates on the right, as well as a covariate matrix has to be specified.

If SNPs genotypes are coded by 0, 1 and 2, they are included as continuous variables in the logistic regression models. If SNPs are coded as carrier SNPs 0 and 1, they are included as binary variables in the logistic regression models. If covariates are used for adjustment, the column names of the covariate matrix \code{cov} have to be specified same as in the formula specification, to link the formula with the covariate matrix \code{cov}.

Missing SNP genotypes in \code{x} or, if used, missing values in \code{cov} are accounted for, as each marginal test makes use of only the available data for that SNP in \code{x} and for that covariate in \code{cov}. \pkg{minPtest} uses all subjects with available data for each SNP (and covariates) when fitting Cochran Armitage Trend Test or unconditional logistic regression. Note that in conditional logistic regression, the matched subjects are removed together in case of 1:1 matching. In the 1:2 matching scenario, matched subjects are removed when the missing occurs in a case, otherwise when a missing occurs in one control, only that control is removed.

Concerning parallelization on a compute cluster, i.e. with argument \code{parallel=TRUE}, there are two possibilities to run \pkg{minPtest}:
\enumerate{
\item Start R on a commandline with sfCluster (Knaus et al., 2009) and preferred options, e.g. number of cpus. The initialization function of package \pkg{snowfall}, \code{sfInit()}, should be called before calling minPtest.

\item Use any other solutions supported by \pkg{snowfall}. Argument \code{parallel} has to be set to TRUE and number of cpus can be chosen in the \code{sfInit()} function.
}

sfCluster is a Unix tool for convenient management of R parallel processes. It is available at
\url{www.imbi.uni-freiburg.de/parallel}, with detailed information.

A print function returns a short overviews of the results. The print function describes the number of subjects included in the analysis, which method is used by the package, briefing of the number of genes, the number of SNPs, the number of missings in the SNP matrix \code{x} and the number of permutations used for the fit. A \code{\link{summary.minPtest}} and a \code{\link{plot.minPtest}} function are available.
}
\value{
An object of class 'minPtest', which is a list containing the following components:
\item{minp}{\code{nrgene * 1} matrix of permutation-based p-values of the min P test for each candidate gene.}
\item{p.adj.minp}{\code{nrgene * 1} matrix of corrected permutation-based p-values via Bonferroni correction method for each candidate gene.}
\item{psnp}{\code{nrsnp * 1} matrix of marginal trend p-values for each SNP from the original data set.}
\item{p.adj.psnp}{\code{nrsnp * 1} matrix of corrected marginal trend p-values via Bonferroni correction method for each SNP from the original data set.}
\item{psnpperm}{\code{nrsnp * n.permute} matrix of permuted trend p-values for each SNP in each permutation step.}
\item{zgen}{\code{nrgene * 1} matrix of min P test statistics for each candidate gene from the original data set.}
\item{zgenperm}{\code{nrgene * n.permute} matrix of permuted min P test statistics for each candidate gene in each permutation step.}
\item{n}{number of subjects in the original data set.}
\item{nrsnp}{number of SNPs in the original data set.}
\item{nrgene}{number of genes in the original data set.}
\item{snp.miss}{number of missings in the SNP matrix \code{x}.}
\item{n.permute}{number of permutations.}
\item{method}{used method.}
\item{call}{call.}
\item{SNPtoGene}{the mapping matrix of dimension \code{p x 2} comprising of SNP names (first column) and names of the genes (second column) on which the SNPs are located.}
}
\references{
Armitage,P. (1955). Tests for linear trends in proportions and frequencies. Biometrics, 11(3), 375-386.

Chen,B.E. et al. (2006). Resampling-based multiple hypothesis testing procedures for genetic case-control association studies. Genetic Epidemiology, 30, 495-507.

Cochran,W.G. (1954). Some methods for strengthening the common chi-squared tests. Biometrics, 10(4), 417-451.

Knaus,J. et al. (2009). Easier parallel computing in R with snowfall and sfCluster. The R Journal, 1, 54-59.

Westfall,P.H. et al.(2002). Multiple tests for genetic effects in association studies. Methods Mol Biol, 184, 143-168.

Westfall,P.H. and Young,S.S. (1993). Resampling-Based Multiple Testing: Example and Methods for p-Value Adjustment. Wiley, New York.
}
\author{
Stefanie Hieke \email{hieke@fdm.uni-freiburg.de}
}
\seealso{
\code{\link{summary.minPtest}}, \code{\link{plot.minPtest}}
}
\examples{
# generate a simulated data set as in the example of the function generateSNPs
# consisting of 100 subjects and 200 SNPs on 5 genes.

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

}
\keyword{single nucleotide polymorphisms}
\keyword{min P test}
\keyword{permutation-based resampling}
