\name{analysis.bgx}
\alias{analysis.bgx}
\alias{plotExpressionDensity}
\alias{plotDEDensity}
\alias{plotDEHistogram}
\alias{plotDERank}
\alias{plotDiffRank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analyse BGX output.}
\description{
  Functions for plotting expression densities, differential expression densities, histogram of proportion of differentially expressed genes, etc.
}
\usage{
  plotExpressionDensity(bgxOutput, gene=NULL, normalize=c("none","mean","loess"),...)
  plotDEDensity(bgxOutput, gene=NULL, conditions=c(1,2), normalize=c("none","mean","loess"), ...)
  plotDEHistogram(bgxOutput, conditions=c(1,2), normalize=c("none","mean","loess"), df=floor(1.8 * log10(length(bgxOutput[["geneNames"]]))))
  plotDERank(bgxOutput, conditions=c(1,2),normalize=c("none", "mean", "loess"))
  plotDiffRank(bgxOutput, conditions=c(1,2),normalize=c("none", "mean", "loess"),ymax=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgxOutput}{A list obtained from running \code{\link{readOutput.bgx}} on a BGX output directory.}
  \item{gene}{Which gene to analyse. This can either be an index or a name.}
  \item{conditions}{Indices of conditions to compare.}
  \item{normalize}{"none": do not normalise posterior distributions of mu. "mean": normalise by scaling posterior distributions of mu for conditions > 1 to have the same mean as the posterior distribution of mu for condition 1. "loess": same as "mean" but use loess normalisation.}
  \item{df}{Residual degrees of freedom. Decrease to 6 if the histogram fit goes haywire.}
  \item{ymax}{Specify upper limit of y axis.}
  \item{...}{Parameters to pass to density function (where applicable).}
}
\details{
  See \code{\link{bgx}} for more details.
}
\value{
  None.
}
\author{Ernest Turro}
\seealso{\code{\link{bgx}}, \code{\link{standalone.bgx}}, \code{\link{readOutput.bgx}}, \code{\link{plotExpressionDensity}}, \code{\link{plotDEDensity}}, \code{\link{plotDEHistogram}}}
\examples{
}
\keyword{IO}% at least one, from doc/KEYWORDS

