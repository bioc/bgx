#  This file is part of BGX, the Bayesian Gene eXpression program.
#  Copyright 2007 Ernest Turro <ernest.turro@ic.ac.uk>
#
#  BGX is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License, version 2, as
#  published by the Free Software Foundation.
#
#  BGX is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


saveAffinityPlot.bgx <- function(originalAffinities, categories, dir, probecat_threshold) {
  path <- file.path(dir,"affinityPlot.pdf")
  
  pdf(path)
  par(mfrow=c(2,2))

  plot(as.numeric(levels(factor(originalAffinities))), table(factor(originalAffinities)),
    xlab="Rounded affinity", ylab="Number of probes", main="Before collapse")
  plot(as.numeric(levels(factor(originalAffinities))), table(factor(originalAffinities)), ylim=c(0,probecat_threshold*2),
    xlab="Rounded affinity", ylab="Number of probes", main="Before collapse\n(blow-up)")


  plot(as.numeric(levels(factor(categories))), table(factor(categories)),
    xlab="Affinity category", ylab="Number of probes", main=paste("After collapse (threshold=", probecat_threshold,")",sep=""))
  plot(as.numeric(levels(factor(categories))), table(factor(categories)), ylim=c(0,probecat_threshold*2),
    xlab="Affinity category", ylab="Number of probes", main=paste("After collapse (threshold=", probecat_threshold,")\n(blow-up)", sep=""))

  dev.off()
  return(path)

}
