\name{setupVars.bgx}
\alias{setupVars.bgx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initialise variables needed to run BGX simulation.}
\description{
  This internal function initialises several variables, which it returns in a list.
}
\usage{
setupVars.bgx(data, samplesets, genes, genesToWatch, probeAff, probecat_threshold, rounding_dec_places = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An \code{AffyBatch} object.}
  \item{samplesets}{A numeric vector specifying which condition each array belongs to. E.g. if samplesets=c(2,2), then the first two replicates belong to one condition and the last two replicates belong to another condition. If NULL, each array is assumed to belong to a different condition.}
  \item{genes}{A numeric vector specifying which genes to analyse. If NULL, all genes are analysed.}
  \item{genesToWatch}{A numeric vector specifying which genes to monitor closely amongst those chosen to be analysed (see below for details).}
  \item{probeAff}{Stratify the mean (lambda) for the cross-hybridisation parameter (H) by categories according to probe-level sequence information.}
  \item{probecat_threshold}{Minimum amount of probes per probe affinity category.}
  \item{rounding_dec_places}{The initial probe categorisation is done by rounding affinities to the nearest \code{rounding_dec_places} decimal places. 1 is a good value.}
}
\details{
}
\value{
  A list:
  \item{pm}{Perfect Match probes.}
  \item{mm}{MisMatch probes.}
  \item{samplesets}{A numeric vector specifying which condition each array belongs to. E.g. if samplesets=c(2,2), then the first two replicates belong to one condition and the last two replicates belong to another condition. If NULL, each array is assumed to belong to a different condition.}
  \item{probesets}{A numeric vector specifying how probes are grouped into probesets.}
  \item{numberOfCategories}{Number of probe affinity categories.}
  \item{categories}{A numeric vector specifying which category each probe belongs to.}
  \item{unknownProbeSeqs}{A numeric vector specifying which probes lack sequence information.}
  \item{numberOfUnknownProbeSeqs}{Number of probes lacking sequence information.}
  \item{genesToWatch}{A numeric vector specifying which genes to monitor closely.}
  \item{firstProbeInEachGeneToWatch}{The starting position for each probe in each gene to monitor closely.}
  \item{numArrays}{Number of arrays.}
}
\references{See \code{\link{bgx}}}
\author{Ernest Turro}
\note{
  This function shouldn't be called directly.
}
\seealso{\code{\link{bgx}}}
\examples{
}
\keyword{internal}% at least one, from doc/KEYWORDS
