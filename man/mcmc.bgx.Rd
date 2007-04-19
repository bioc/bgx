\name{mcmc.bgx}
\alias{mcmc.bgx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal wrapper function for calling the bgx C++ function.}
\description{
  This internal function calls the \code{bgx} method in a loaded bgx shared object (bgx.so/bgx.dll)
}
\usage{
mcmc.bgx(pm, mm, samplesets, probesets, numberCategories, categories, unknownProbeSeqs, numberOfUnknownProbeSeqs, 
  numberGenesToWatch, whichGenesToWatch, whichProbesToWatch, iter, burnin,
  adaptive, batch_size=50, optimalAR=0.44, output, samplenames = "unknown",
  subsample = ifelse(iter > 1024, iter/1024, 1), seed = 192492, basepath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pm}{Perfect Match probes}
  \item{mm}{MisMatch probes}
  \item{samplesets}{A numeric vector specifying which condition each array belongs to. E.g. if samplesets=c(2,2), then the first two replicates belong to one condition and the last two replicates belong to another condition. If NULL, each array is assumed to belong to a different condition.}
  \item{probesets}{A numeric vector specifying how probes are grouped into probesets.}
  \item{numberCategories}{Number of probe affinity categories.}
  \item{categories}{A numeric vector specifying which category each probe belongs to.}
  \item{unknownProbeSeqs}{A numeric vector specifying which probes lack sequence information.}
  \item{numberOfUnknownProbeSeqs}{Number of probes lacking sequence information.}
  \item{numberGenesToWatch}{How many genes to monitor closely.}
  \item{whichGenesToWatch}{A numeric vector specifying which genes to monitor closely.}
  \item{whichProbesToWatch}{The starting position for each probe in each gene to monitor closely.}
  \item{iter}{Number of post burn-in iterations.}
  \item{burnin}{Number of burn-in iterations.}
  \item{adaptive}{Use adaptive MCMC for better mixing.}
  \item{batch_size}{Batch size for adaptive MCMC.}
  \item{optimalAR}{Optimal acceptance ratio.}
  \item{output}{One of "minimal", "trace", "diagnostic" or "mcse".}
  \item{samplenames}{Vector of names for each array.}
  \item{subsample}{Subsampling interval.}
  \item{seed}{Seed for PRNG.}
  \item{basepath}{The directory in which to place the output run directories. By default this is the "bgx" directory in the R temporary directory, which you can obtain by using \code{tempdir()}.}
}
\details{
See \code{\link{bgx}} for more details.
}
\value{
  The name of the output directory.
}
\author{Ernest Turro}
\note{You shouldn't call this function directly, but if you do, make sure the appropriate shared object is loaded.
}
\seealso{\code{\link{bgx}}, \code{\link{standalone.bgx}}}
\examples{
}
\keyword{internal}% at least one, from doc/KEYWORDS
