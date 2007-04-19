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


"bgx" <-
function(aData,samplesets=NULL,genes=NULL,genesToWatch=NULL,burnin=16384,iter=65536,output=c("minimal","trace","all"), serial=FALSE, probeAff=TRUE, probecat_threshold = 100, adaptive=TRUE, basepath=file.path(tempdir(),"bgx")) {
#  if(burnin %% 1024 != 0 || iter %% 1024 != 0)
#    stop("\"iter\" and \"burnin\" must be a multiple of 1024")
  # create directory where runs will be saved if necessary
  if(!file.exists(basepath)) dir.create(basepath)

  # do not analyse the same gene more than once, even if specified in the arguments
  genes <- unique(genes); genesToWatch <- unique(genesToWatch)

  output <- match.arg(output)

  # not being able to pass by reference sucks.  
  vars <- setupVars.bgx(data=aData,samplesets=samplesets,genes=genes,genesToWatch=genesToWatch,probeAff=probeAff, probecat_threshold=probecat_threshold)
  pm<-vars$pm
  mm<-vars$mm
  samplesets<-vars$samplesets
  probesets<-vars$probesets
  categories<-vars$categories
  numberOfCategories<-vars$numberOfCategories
  unknownProbeSeqs<-vars$unknownProbeSeqs
  numberOfUnknownProbeSeqs<-vars$numberOfUnknownProbeSeqs
  originalAffinities<-vars$originalAffinities
  genes<-vars$genes
  genesToWatch<-vars$genesToWatch
  numberOfGenesToWatch=vars$numberOfGenesToWatch
  firstProbeInEachGeneToWatch = vars$firstProbeInEachGeneToWatch
  numArrays<-vars$numArrays

  if(serial) {
    samplesets=1
    for(c in 1:numArrays) {
      outdir=mcmc.bgx(as.matrix(pm[,c]),as.matrix(mm[,c]),samplesets,probesets,numberOfCategories, categories,
        unknownProbeSeqs, numberOfUnknownProbeSeqs, numberOfGenesToWatch, genesToWatch, 
        firstProbeInEachGeneToWatch,iter,burnin,adaptive, output=output,samplenames=sampleNames(aData)[c], basepath=basepath) 
      cat("CEL file ",sampleNames(aData)[c]," analysed.\n")
    }
  } else {
    outdir=mcmc.bgx(pm,mm,samplesets,probesets,numberOfCategories, categories, unknownProbeSeqs, 
    numberOfUnknownProbeSeqs,
    numberOfGenesToWatch,genesToWatch,firstProbeInEachGeneToWatch,iter,burnin,adaptive, output=output,samplenames=sampleNames(aData), basepath=basepath)
    cat("CEL files analysed: ")
    for(c in 1:numArrays) cat(sampleNames(aData)[c]," ")
  }
  
  if(probeAff) saveAffinityPlot.bgx(originalAffinities, categories, outdir, probecat_threshold)
  write(geneNames(aData)[genes], file=file.path(outdir,"geneNames.txt"))
  
  cat("\nFiles in '",outdir,"'\n")
  cat("Run(s) finished.\n")

  # this should be an argument
  conditionnames <- paste("condition",c(1:length(samplesets)))

  return(new('ExpressionSet',
    exprs = matrix(data=scan(file.path(outdir,"muave"),quiet=TRUE),nrow=length(genes),ncol=length(samplesets),
      dimnames=list(geneNames(aData)[genes],conditionnames)),
    phenoData = new("AnnotatedDataFrame", data=as.data.frame(matrix(c(1:length(samplesets)),dimnames=list(c(conditionnames),c("condition")))), varMetadata=as.data.frame(list(labelDescription="Arbitrary numbering"), row.names="condition")),
    annotation = annotation(aData),
    experimentData = description(aData)))
}

