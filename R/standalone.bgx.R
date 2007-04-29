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


"standalone.bgx" <-
function(aData,samplesets=NULL,genes=NULL,genesToWatch=NULL,burnin=16384,iter=65536, output=c("minimal","trace", "all"), probeAff=TRUE, probecat_threshold=100, adaptive = TRUE, batch_size=50, optimalAR=0.44, dirname="input") {
#  if(burnin %% 1024 != 0 || iter %% 1024 != 0) stop("\"iter\" and \"burnin\" must be a multiple of 1024")
  output <- match.arg(output)
  if(length(grep("[[:space:]]", dirname)) > 0) {
    cat("Warning: removing space characters from dirname.")
    dirname <- gsub("[[:space:]]","", dirname)
  }
  # not being able to pass by reference sucks.  
  vars <- setupVars.bgx(data=aData,samplesets=samplesets,genes=genes,genesToWatch=genesToWatch,probeAff=probeAff,probecat_threshold=probecat_threshold)
  pm<-vars$pm
  mm<-vars$mm
  samplesets<-vars$samplesets
  probesets<-vars$probesets
  categories<-vars$categories
  unknownProbeSeqs<-vars$unknownProbeSeqs
  numberOfUnknownProbeSeqs<-vars$numberOfUnknownProbeSeqs
  numberOfCategories<-vars$numberOfCategories
  originalAffinities<-vars$originalAffinities
  genes<-vars$genes
  genesToWatch<-vars$genesToWatch
  numberOfGenesToWatch<-vars$numberOfGenesToWatch
  firstProbeInEachGeneToWatch<-vars$firstProbeInEachGeneToWatch
  numArrays<-vars$numArrays

  while(!is.na(file.info(dirname)$size)) {
    cat("I tried to create a new directory named \"",dirname,"\" in which to place the standalone BGX files,
    but a file or directory with that name already exists.\n
    Would you like to:
    \t* overwrite it (o)
    \t* specify a different directory name (d)
    \t* abort (A)?\n",sep="")
    answer = scan(what=character(1), n=1)
    if(answer=="o" || answer=="O") {
      cat("Overwriting \"", dirname, "\"\n",sep="")
      unlink(dirname,recursive=TRUE)
    } else if(answer=="d" || answer=="D") {
      cat("Please choose a different directory name:\n")
      dirname = scan(what=character(), n=1)
    } else stop("Aborting.")
  }
  
  dir.create(dirname) 

  write(pm,file.path(dirname, "PM.txt"))
  write(mm,file.path(dirname,"MM.txt"))
  write(probesets,file.path(dirname,"PS.txt"))

  write(samplesets,file.path(dirname,"SS.txt"))

  write(categories,file.path(dirname,"categories.txt"))
  write(unknownProbeSeqs,file.path(dirname,"unknownProbeSeqs.txt"))
  
  write(genesToWatch,file.path(dirname,"genesToWatch.txt"))
  write(firstProbeInEachGeneToWatch,file.path(dirname,"firstProbeInEachGeneToWatch.txt"))

  # create category summary table
  cats=sort(unique(categories))
  numInEachCat=vector(mode="integer",length=length(cats))
  for( c in 1: length(cats)) numInEachCat[c]=length(categories[categories==cats[c]])
  categorySummaryTable=matrix(nrow=length(cats),ncol=1)
  categorySummaryTable[,1]=numInEachCat
  rownames(categorySummaryTable)=cats
  colnames(categorySummaryTable)="occurances"
  write.table(categorySummaryTable,file.path(dirname,"categorySummaryTable.csv"))
  if(probeAff) affinityPlotFile <- saveAffinityPlot.bgx(originalAffinities, categories, dirname, probecat_threshold)
  write(geneNames(aData)[genes], file=file.path(dirname,"geneNames.txt"))
    
  cat(paste("samples ",numArrays,"\n",
    "conditions ",length(samplesets),"\n",
    "probes ",dim(pm)[1],"\n",
    "genes ",length(probesets),"\n",
    "geneNamesFile ", file.path(dirname, "geneNames.txt"), "\n",
    "numberOfCategories ",length(unique(categories)),"\n",
    "numberOfGenesToWatch ",numberOfGenesToWatch,"\n",
    "SampleSets ",file.path(dirname,"SS.txt"),"\n",
    "ProbeSets ",file.path(dirname,"PS.txt"),"\n",
    "Categories ",file.path(dirname,"categories.txt"),"\n",
    "UnknownProbeSeqs ", file.path(dirname,"unknownProbeSeqs.txt"),"\n",
    "numberOfUnknownProbeSeqs ", numberOfUnknownProbeSeqs,"\n",
    "genesToWatch ",file.path(dirname,"genesToWatch.txt"),"\n",
    "firstProbeInEachGeneToWatch ",file.path(dirname,"firstProbeInEachGeneToWatch.txt"),"\n",
    "PM ",file.path(dirname,"PM.txt"),"\n",
    "MM ",file.path(dirname,"MM.txt"),"\n",
    "seed ",192492,"\n",
    "sweeps ",iter,"\n",
    "burn-in ",burnin,"\n",
    "Output ",output,"\n",
    "Adaptive ", as.integer(adaptive), "\n",
    "BatchSize ", as.integer(batch_size), "\n",
    "OptimalAR ", as.double(optimalAR), "\n",
    "S_jmp 30\n",
    "H_jmp 350\n",
    "Mu_jmp 1.1\n",
    "Tau_jmp 1.5\n",
    "Lambda_jmp 0.04\n",
    "Eta_jmp 0.1\n",
    "CELfiles ", paste(sampleNames(aData), collapse=" "), "\n",
    if(probeAff)"affinityPlotFile ", if(probeAff)affinityPlotFile, sep=""),
    file=file.path(dirname,"infile.txt"),append=TRUE)

    cat("Standalone BGX files are in \"",dirname,"\".\n",sep="")
    return(dirname)
}

