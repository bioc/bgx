\name{bgx}
\alias{bgx}
\alias{standalone.bgx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fully Bayesian integrated approach to the analysis of Affymetrix GeneChip data}
\description{
  'bgx' estimates Bayesian Gene eXpression (BGX) measures from an AffyBatch object.

  'standalone.bgx' creates various files needed by the bgx standalone binary and places them in a directory. One of these files is 'infile.txt'. In order to run standalone BGX, compile it and run 'bgx <path\_to\_infile.txt>' from the command line.
}
\usage{
bgx(aData, samplesets = NULL, genes = NULL, genesToWatch = NULL,
  burnin = 16384, iter = 65536, output = c("minimal","trace","all"), serial = FALSE, 
  probeAff = TRUE, probecat_threshold = 100, adaptive = TRUE, basepath = file.path(tempdir(), "bgx"))

standalone.bgx(aData, samplesets = NULL, genes = NULL, genesToWatch = NULL,
  burnin = 16384, iter = 65536, output = c("minimal", "trace", "all"),
  probeAff = TRUE, probecat_threshold = 100,
  adaptive = TRUE, batch_size = 50, optimalAR = 0.44, dirname = "input")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aData}{An \code{AffyBatch} object.}
  \item{samplesets}{A numeric vector specifying which condition each array belongs to. E.g. if samplesets=c(2,2), then the first two replicates belong to one condition and the last two replicates belong to another condition. If NULL, each array is assumed to belong to a different condition. If the aData object contains information about the experiment design in its \code{phenoData} slot, this argument is not required.}
  \item{genes}{A numeric vector specifying which genes to analyse. If NULL, all genes are analysed.}
  \item{genesToWatch}{A numeric vector specifying which genes to monitor closely amongst those chosen to be analysed (see below for details).}
  \item{burnin}{Number of burn-in iterations.}
  \item{iter}{Number of post burn-in iterations.}
  \item{output}{One of "minimal", "trace" or "all". See below for details.}
  \item{serial}{Run a separate simulation for each individual array.}
  \item{probeAff}{Stratify the mean (lambda) of the cross-hybridisation parameter (H) by categories according to probe-level sequence information.}
  \item{probecat_threshold}{Minimum amount of probes per probe affinity category.}
  \item{adaptive}{Adapt the variance of the proposals for Metropolis Hastings objects (that is: S, H, Lambda, Eta, Sigma and Mu).}
  \item{batch_size}{Size of batches for calculating acceptance ratios and updating jumps.}
  \item{optimalAR}{Optimal acceptance ratio.}
  \item{basepath}{The directory in which to place the output run directories. By default this is the "bgx" directory in the R temporary directory, which you can obtain by using \code{tempdir()}.}
  \item{dirname}{The name of the directory in which to place the files for standalone binary.}
}
\details{
  \itemize{
    \item{genesToWatch}{Specify the subset of genes for which thinned samples from the full posterior distributions of log(S+1) (x) and log(H+1) (y) are collected.}
    \item{output}{Output the following to disk:}
    \itemize{
      \item{"minimal"}{The gene expression measure (muave), thinned samples from the full posterior distributions of mu (mu.[1..c]), where 'c' is the number of conditions, the integrated autocorrelation time (IACT) and the Markov chain Monte Carlo Standard Error (MCSE) for each gene under each condition. Note that the IACT and MCSE are calculated from the thinned samples of mu.}
      \item{"trace"}{The same as "minimal" plus thinned samples from the full posterior distributions of sigma2 (sigma2.[1..c]), lambda (lambda.[1..s]), eta2 (eta2), phi (phi) and tau2 (tau2), where 's' is the number of samples. If there are probes with unknown sequences, output a thinned trace of their categorisation.}
      \item{"all"}{The same as "trace" plus acceptance ratios for S (sacc), H (hacc), mu (muacc), sigma (sigmaacc), eta (etaacc) and lambda (lambdasacc).}
    }
  }
}
\value{
  'bgx' returns an \code{ExpressionSet} object containing gene expression information for each gene under each condition (not each replicate).

  'standalone.bgx' returns the path to the BGX input files.
}
\references{
Hein, A., Richardson, S., Causton, H., Ambler, G., Green., P. (2005) BGX: a fully Bayesian integrated approach to the analysis of Affymetrix GeneChip data. Biostatistics (2005), 6, 3, pp. 349-373.

Hekstra, D., Taussig, A. R., Magnasco, M., and Naef, F. (2003) Absolute mRNA concentrations from sequence-specific calibration of oligonucleotide array. Nucleic Acids Research, 31. 1962-1968.

G.O. Roberts and J.S. Rosenthal (September, 2006) Examples of Adaptive MCMC.}
\author{Ernest Turro}
\note{The bgx() method and the bgx standalone binary create a directory in the working directory called 'run.x' (x:1,2,3,...), wherein files are placed for further detailed analysis.
}
\examples{
  if(require(affydata) && require(hgu95av2cdf)) {
    data(Dilution)
    eset <- bgx(Dilution, samplesets=c(2,2), probeAff=F, burnin=2048, iter=8192,
      genes=c(12500:12599), output="all", basepath=".")
  }
}
\keyword{manip}% at least one, from doc/KEYWORDS
