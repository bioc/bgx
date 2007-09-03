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


"setupVars.bgx" <-
function(data,samplesets,genes,genesToWatch,probeAff,probecat_threshold, rounding_dec_places=1) {
    calcProbeAffCategories <- function(data, theProbes) {
      affinitiesOfMMprobes<- as.vector(round(mm(gcrma::compute.affinities(cdfName(data))),rounding_dec_places))[theProbes]
      affinitiesOfMMprobesOrig <- affinitiesOfMMprobes
      affinitiesOfMMprobesOrig[is.na(affinitiesOfMMprobesOrig)] <- median(affinitiesOfMMprobesOrig, na.rm=TRUE)

      sumUpTo <- function(t,c) {
        sum <- 0
        for(i in 1:c) sum <- sum + t[[i]]
        return(sum)
      }
    
      sumDownTo <- function(t,c) {
        sum <- 0
        for(i in length(t):c) sum <- sum + t[[i]]
        return(sum)
      }
      
      ## collapse categories with fewer than probecat_threshold probes.
      startAffinity <- names(which.max(table(factor(affinitiesOfMMprobes))))
      
      # bottom half 
      affinityTable <- table(factor(affinitiesOfMMprobes))
      numberOfProbesInEachCategory <- as.vector(affinityTable)
      distinctRoundedAffinities <- as.numeric(levels(factor(affinitiesOfMMprobes)))
      start <- (1:length(affinityTable))[names(affinityTable)==startAffinity]

      temp <- 0
      offset <- 0.0
      concat <- FALSE
      toofewleft <- FALSE
      for(c in start:1){
        temp <- temp + affinityTable[[c]]
        if(toofewleft || sumUpTo(affinityTable,c-1) < probecat_threshold) toofewleft <- TRUE
        if((concat || toofewleft) && c != start) {
          offset <- offset + as.numeric(names(affinityTable))[c+1] - as.numeric(names(affinityTable))[c]
          affinitiesOfMMprobes[affinitiesOfMMprobes == as.numeric(names(affinityTable)[c]) & !is.na(affinitiesOfMMprobes)] <- as.numeric(names(affinityTable[c])) + offset
        }
        if(!toofewleft && temp>=probecat_threshold) {
          temp <- 0
          offset <- 0.0
          concat <- FALSE
        } else concat <- TRUE
      }
      affinitiesOfMMprobes <- round(affinitiesOfMMprobes, rounding_dec_places)
     
      # top half 
      affinityTable <- table(factor(affinitiesOfMMprobes))
      numberOfProbesInEachCategory <- as.vector(affinityTable)
      distinctRoundedAffinities <- as.numeric(levels(factor(affinitiesOfMMprobes)))
      start <- (1:length(affinityTable))[names(affinityTable)==startAffinity]
      
      temp <- 0
      offset <- 0.0
      concat <- FALSE
      toofewleft <- FALSE
      for(c in start:length(affinityTable)){
        temp <- temp + affinityTable[[c]]
        if(toofewleft || sumDownTo(affinityTable,c+1) < probecat_threshold) toofewleft <- TRUE
        if((concat || toofewleft) && c != start) {
          offset <- offset + (as.numeric(names(affinityTable))[c] - as.numeric(names(affinityTable))[c-1])
          affinitiesOfMMprobes[affinitiesOfMMprobes == as.numeric(names(affinityTable)[c]) & !is.na(affinitiesOfMMprobes)] <- as.numeric(names(affinityTable[c])) - offset
        }
        if(!toofewleft && temp>=probecat_threshold) {
          temp <- 0
          offset <- 0.0
          concat <- FALSE
        } else concat <- TRUE
      }
      affinitiesOfMMprobes <- round(affinitiesOfMMprobes, rounding_dec_places)

      distinctRoundedAffinities <- as.numeric(levels(factor(affinitiesOfMMprobes)))
      categories <- vector(mode="numeric", length=length(affinitiesOfMMprobes))
      for(i in 1:length(distinctRoundedAffinities)){
        categories[affinitiesOfMMprobes == distinctRoundedAffinities[i] & !is.na(affinitiesOfMMprobes)] <- i-1
      }

      unknownProbeSeqs <- (1:length(affinitiesOfMMprobes))[is.na(affinitiesOfMMprobes)] - 1# no sequence information; indexed from 0, not 1
      categories[is.na(affinitiesOfMMprobes)] <- median(categories, na.rm=TRUE) 
       # floor(runif(length(unknownProbeSeqs))*length(unique(categories)))

      return(list(categories=categories,unknownProbeSeqs=unknownProbeSeqs, originalAffinities=affinitiesOfMMprobesOrig))
    }

    integersToNiceString <- function(vec) {
      if(length(vec)< 2) return(as.character(vec))
    
      c <- as.character(vec[1])
      prev <- vec[1]
      len <- length(vec)
      
      for(i in 2:len){
        if(prev!=vec[i]-1){
          if(i==2 || vec[i-2] != prev-1) c <- paste(c,",",sep="") 
          else c <- paste(c,":",prev,",",sep="")
        }
        if(i==len){
            if(vec[i-1] == vec[i]-1) c <- paste(c,":",vec[i],sep="")
            else c <- paste(c,vec[i],sep="")
        } else if(substr(c,nchar(c),nchar(c))==",")
            c <- paste(c,vec[i],sep="")
        prev <- vec[i]
      }
      return(c)
    }
    
    # START HERE
    numArrays<-ncol(pm(data))
    if(is.null(samplesets)){
      if(!is.null(pData(data)[,1])) {
        samplesets <- c()
        cat("Note: Samplesets not specified. Using phenoData.\n")
        for(con in levels(factor(pData(data)[,1]))) {
          samplesets <- append(samplesets, sum(pData(data)[,1] == con))
        }
      } else {
        samplesets=rep(1,numArrays)
      }
    }
    else if(sum(samplesets)!=numArrays)
      stop("Samplesets not sensible.  You have ", numArrays, " arrays, but sum(samplesets) = ",sum(samplesets))
    cat("Analyse",numArrays,"array(s) in", length(samplesets), "condition(s):\n")
    for(a in 1:length(samplesets)){
      cat("\t- condition ",a,": ",samplesets[a], " array(s)\n",sep="")
    }

    if(is.null(genes)) {
      cat("Analyse all genes\n")
      geneNames<-geneNames(data)
      genes<-seq(1,length(geneNames))
    } else {
      cat("Analyse genes ",integersToNiceString(genes),"\n")
      geneNames<-geneNames(data)[genes]
    }

    if(any(!(genesToWatch %in% genes)))
      stop("Sorry, 'genesToWatch' contains genes not in c(1:genes).")
    
    pm<-pm(data, geneNames)
    mm<-mm(data, geneNames)
    
    ip<-indexProbes(data,"pm")[genes]
    probesets <- vector(mode="integer", length=length(ip))
    for(i in 1:length(ip)) probesets[i]<-length(ip[[i]])
    
    # to calculate the probe affinities of the analysed probes we need 
    # to know which ones, among the full set of probes, they are:
    ipAll=indexProbes(data,"pm")
    probesetsAll <- vector(mode="integer", length=length(ipAll))
    for(i in 1:length(ipAll)) probesetsAll[i]=length(ipAll[[i]])
    
    if(probeAff){
      cat("Take into account probe affinities (threshold is ", probecat_threshold, " probes per category)\n")
      countProbes=0
      theProbes=vector(mode="integer", length=sum(probesetsAll[genes]))
      for(g in 1:length(genes)){
        theProbes[(countProbes+1):(countProbes+probesetsAll[genes[g]])]=c((sum(probesetsAll[1:(genes[g])])-probesetsAll[genes[g]]+1):sum(probesetsAll[1:genes[g]]))
        countProbes = countProbes + probesetsAll[genes[g]]
      }
      cats <- calcProbeAffCategories(data, theProbes)
      categories <- cats$categories
      numberOfCategories=length(unique(categories)) 
      unknownProbeSeqs <- cats$unknownProbeSeqs
      originalAffinities <- cats$originalAffinities
      numberOfUnknownProbeSeqs <- length(unknownProbeSeqs)
      cat("There are ",numberOfCategories, " probe affinity categories\n")
      cat("There are ",numberOfUnknownProbeSeqs, " probes with no sequence information (initial category set to median)\n")
    } else {
      categories=vector(mode="integer",length=length(pm[,1]))
      numberOfCategories=1
      unknownProbeSeqs=c()
      numberOfUnknownProbeSeqs=0
    }
    
    numberOfGenesToWatch <- length(genesToWatch)
    if(is.null(genesToWatch)) {
      genesToWatch <- numeric(0) 
      firstProbeInEachGeneToWatch <- numeric(0)
    } else {
      cat("Parameter 'genesToWatch' set. Monitor genes ", integersToNiceString(genesToWatch)," in detail.\n")
      firstProbeInEachGeneToWatch <- vector(mode="integer", length=length(genesToWatch))
      for(g in 1:length(genesToWatch)){
        firstProbeInEachGeneToWatch[g] <- sum(probesets[1:genesToWatch[g]-1]) + 1
      }
      genesToWatch <- genesToWatch - 1 # In C we start from 0
      firstProbeInEachGeneToWatch <- firstProbeInEachGeneToWatch - 1 # ibid
    }
    
    return(list(pm=pm,mm=mm,samplesets=samplesets,probesets=probesets,
      categories=categories,numberOfCategories=numberOfCategories,unknownProbeSeqs=unknownProbeSeqs,
      numberOfUnknownProbeSeqs=numberOfUnknownProbeSeqs,originalAffinities=if(probeAff)originalAffinities, genes=genes,
      genesToWatch=genesToWatch, numberOfGenesToWatch=numberOfGenesToWatch,
      firstProbeInEachGeneToWatch = firstProbeInEachGeneToWatch, numArrays=numArrays))
}

