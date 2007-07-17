########################################################################
#
#  This file is part of BGX, the Bayesian Gene eXpression program.
#  Copyright (c) 2003-2004  Graeme Ambler <graeme@ambler.me.uk>
#                2004       Anne-Mette Hein <a.hein@imperial.ac.uk>
#                2007       Ernest Turro <ernest.turro@ic.ac.uk>
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
#
########################################################################

"mcmc.bgx" <-
function(pm,mm,samplesets,probesets,numberCategories,categories,unknownProbeSeqs,numberOfUnknownProbeSeqs,
  numberGenesToWatch, whichGenesToWatch,whichProbesToWatch,iter,burnin,adaptive, batch_size=50, optimalAR=0.44,output,
  samplenames="unknown", subsample=ifelse(iter>1024,iter/1024,1),seed=192492, rundir) {

  #make an indicator variable for what type of output we want
  if(output=="minimal") out.ind<-0
  else if(output=="trace") out.ind<-1
  else if(output=="all") out.ind<-2
  else stop("Invalid value for \"output\" parameter.")

  cat("Starting MCMC simulation...\n")
 
  # free allocated memory on user interrupt/end of simulation
  on.exit(.C("freeBGXMemory", as.integer(out.ind), as.integer(numberGenesToWatch), PACKAGE = "bgx"))

  a<-.C("bgx",
        as.double(pm),#pm
        as.double(mm),#mm
        as.integer(dim(pm)[2]),#samples
        as.integer(length(samplesets)),#conditions
        as.integer(dim(pm)[1]),#probes
        as.integer(length(probesets)),#genes
        as.integer(numberCategories),#numberCategories
        as.integer(numberGenesToWatch),#number genes to monitor
        as.integer(samplesets),#structure of samples
        as.integer(probesets),#structure of probes
        as.integer(categories),#array of probe affinity categories (one for each probe pair) 
        as.integer(unknownProbeSeqs), #array of indices for probes with no sequence information
        as.integer(numberOfUnknownProbeSeqs), #number of probes with no sequenece information
        as.integer(whichGenesToWatch),#gene numbers among analysed to monitor fully
        as.integer(whichProbesToWatch),#first probe numbers for genes to monitor fully
        as.integer(iter),#sweeps
        as.integer(burnin),#burnin
        as.double(30),#s_jump 
        as.double(350),#h_jump
        as.double(1.1),#mu_jump
        as.double(1.5),#tau_jump
        as.double(0.04),#lambda_jump
        as.double(0.1),#eta_jump
        as.logical(adaptive),#adaptive mcmc
        as.integer(batch_size),#batch size for adapting
        as.double(optimalAR),#optimal acceptance ratio
        as.integer(subsample),#subsampling interval
        as.integer(out.ind),#indicator of what output we want
        as.character(paste(character(length=length(rundir)+16), collapse="x")), # output path returned in here
        as.character(rundir),# base path in which to put run directories 
        as.integer(seed),#seed
        as.character(samplenames), # chip names 
        as.character(NULL),# no need to copy affinity plot to output directory 
        as.character(NULL), # no need to copy geneNames file to output directory
        PACKAGE = "bgx")
  #return the name of the output directory
  return(a[[29]])
}

