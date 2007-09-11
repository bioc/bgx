\name{readOutput.bgx}
\alias{readOutput.bgx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read in the output from a BGX run.}
\description{
  readOutput.bgx reads in output from BGX which can then be fed into BGX analysis functions.
}
\usage{
  readOutput.bgx(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{Paths of BGX output directories. If you specify more than one path, then the runs will be combined such that each condition from each run is treated as a different different from all the others.}
}
\details{
  See \code{\link{bgx}} for more details.
}
\value{
  A list containing data from the BGX output.
}
\author{Ernest Turro}
\seealso{\code{\link{bgx}}, \code{\link{standalone.bgx}}, \code{\link{plotExpressionDensity}}, \code{\link{plotDEDensity}}, \code{\link{plotDEHistogram}}}
\examples{
}
\keyword{IO}% at least one, from doc/KEYWORDS

