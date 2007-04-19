#  This file is part of BGX, the Bayesian Gene eXpression program.
#  Copyright 2006, 2007 Anne-Mette K Hein, Ernest Turro <ernest.turro@ic.ac.uk>
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



### Read in BGX output
readOutput.bgx <- function(path) {
  cat("Reading '",path, "'\n", sep="")
  summary <- read.delim(file.path(path,"summary.txt"))
  print(summary)
  cat('\n')
  
  summary <- as.vector(summary[[1]])
  noOfConditions <- as.numeric(summary[[2]])
  noOfGenes <- as.numeric(summary[[4]])
  
  mu <- list() 
  for(m in 1:noOfConditions) {
    cat("Reading mu under condition ", m, "...", sep="")
    mu[[m]] <- matrix(scan(file.path(path, paste("mu.",m,sep=""))), ncol=1024)
  }
  cat("Reading mu average...")
  muave <- matrix(scan(file.path(path, "muave")), ncol=noOfConditions)
  cat("Reading gene names...")
  geneNames <- scan(file.path(path, "geneNames.txt"), what="character")
  
  cat("Done.\n")

  return(list(mu=mu, muave=muave, geneNames=geneNames))
}

### Plot densities of mu's under all conditions
plotExpressionDensity <- function(bgxOutput, gene=NULL, normalize=c("none","mean","loess"), ...) {
  if(is.null(gene)) stop("Please specify a genes index/name.")
  if(is.character(gene)) gene <- (1:length(bgxOutput$geneNames))[gene == bgxOutput$geneNames]
  normalize <- match.arg(normalize)
  if(normalize=="mean") bgxOutput$mu <- meanNorm(bgxOutput$mu)
  else if(normalize=="loess") bgxOutput$mu <- loessNorm(bgxOutput$mu)

  noOfConditions <- length(bgxOutput$mu)
  ymax <- 0
  xmin <- 9999
  xmax <- -9999
  densities <- list(length=noOfConditions)
  
  for(i in 1:noOfConditions) {
    densities[[i]] <- density(bgxOutput$mu[[i]][gene,], ...)
    yma <- max(densities[[i]]$y)
    xmi <- min(densities[[i]]$x)
    xma <- max(densities[[i]]$x)
    if(ymax < yma) ymax <- yma
    if(xmax < xma) xmax <- xma
    if(xmin > xmi) xmin <- xmi
  }

  plot(NA, type="n", xlab="Expression", ylab="Density", main=paste("Densities of mu for gene", bgxOutput$geneNames[gene]), xlim=c(xmin,xmax), ylim=c(0,ymax+0.2))# leave some space for legend

  for(i in 1:noOfConditions) lines(densities[[i]], lty=i)

  legend("topright", paste("Cond", c(1:noOfConditions)), lty=c(1:noOfConditions), cex=0.75)

}

### Plot densities of differential expression between two conditions
plotDEDensity <- function(bgxOutput, gene=NULL, conditions=c(1,2),normalize=c("none","mean","loess"), ...) {
  if(is.null(gene)) stop("Please specify a genes index/name.")
  if(is.character(gene)) gene <- (1:length(bgxOutput$geneNames))[gene == bgxOutput$geneNames]
  if(length(conditions) !=2 ) stop("Please specify exactly two conditions (e.g. conditions=c(1,2)).")
  normalize <- match.arg(normalize)
  if(normalize=="mean") bgxOutput$mu <- meanNorm(bgxOutput$mu, target=conditions[1])
  else if(normalize=="loess") bgxOutput$mu <- loessNorm(bgxOutput$mu, target=conditions[1])

  diff <- bgxOutput$mu[[conditions[2]]][gene,] - bgxOutput$mu[[conditions[1]]][gene,] 
  
  plot(density(diff, ...), main=paste("Density of mu for gene", bgxOutput$geneNames[gene], "\ncondition ", conditions[2] , " - condition", conditions[1]), xlab="Differential Expression")
  abline(v=0)

}

### Plot sorted DE of genes 
plotDERank <- function(bgxOutput, conditions=c(1,2),normalize=c("none", "mean", "loess")) {
  normalize <- match.arg(normalize)
  if(normalize=="mean") bgxOutput$mu <- meanNorm(bgxOutput$mu, target=conditions[1])
  else if(normalize=="loess") bgxOutput$mu <- loessNorm(bgxOutput$mu, target=conditions[1])

  mu_diff <- bgxOutput$mu[[conditions[2]]] - bgxOutput$mu[[conditions[1]]]
  de <- c()
  for(i in 1:nrow(bgxOutput$mu[[1]])) de[i] <- abs(mean(mu_diff[i,]))/sd(mu_diff[i,])
  order <- sort(de,decreasing=TRUE, index.return=TRUE)$ix

  plot(c(1:nrow(bgxOutput$mu[[1]])), de[order], cex=0.5, pch=19, xlab="rank", ylab="Differential expression", 
    main=paste("Sorted differential expression between cond", conditions[1], " and cond",conditions[2],"\nnormalisation: ",normalize))
}

plotDiffRank <- function(bgxOutput, conditions=c(1,2),normalize=c("none", "mean", "loess"), ymax=NULL) {
  normalize <- match.arg(normalize)
  if(normalize=="mean") bgxOutput$mu <- meanNorm(bgxOutput$mu, target=conditions[1])
  else if(normalize=="loess") bgxOutput$mu <- loessNorm(bgxOutput$mu, target=conditions[1])

  mu_diff <- bgxOutput$mu[[conditions[2]]] - bgxOutput$mu[[conditions[1]]]
  
  if(is.null(ymax)) ymax <- nrow(mu_diff)
  
  quantiles <- matrix(ncol=2, nrow=nrow(mu_diff))
  for(r in 1:nrow(mu_diff)){
    q <- quantile(mu_diff[r,], probs=c(0.025,0.975))
    quantiles[r,1] <- q[[1]]
    quantiles[r,2] <- q[[2]]
   }
  
  plot(0,type="n", xlim=c(min(quantiles),max(quantiles)),ylim=c(0,ymax), xlab="CI 2.5%, 97.5%", ylab="Ranked gene index")
  
  diff <- c()
  for(i in 1:nrow(bgxOutput$mu[[1]])) diff[i] <- mean(mu_diff[i,])/sd(mu_diff[i,])
  order <- sort(diff,decreasing=TRUE, index.return=TRUE)$ix
  
  for(r in 1:ymax) segments(quantiles[order[r],1],r, quantiles[order[r],2],r)

  abline(v=0, col=2)
  
}

### Estimate proportion of DE genes. Method similar to Efron
plotDEHistogram <- function(bgxOutput, conditions=c(1,2), normalize=c("none", "mean", "loess"), df=7) {
  normalize <- match.arg(normalize)
  if(normalize=="mean") bgxOutput$mu <- meanNorm(bgxOutput$mu, target=conditions[1])
  else if(normalize=="loess") bgxOutput$mu <- loessNorm(bgxOutput$mu, target=conditions[1])

  diff <- bgxOutput$mu[[conditions[2]]] - bgxOutput$mu[[conditions[1]]]
  pp <- c()
  for(i in 1:length(bgxOutput$geneNames)) pp[i] <- mean(diff[i,] < 0)

  bre = seq(0,1,by=0.025)
#  df = 7
  pct = 0
  pct0 = 2/3
  nulltype = 1
  type = 1

  require("splines")

  if(length(bre) > 1) {
    lo <- min(bre)
    up <- max(bre)
    bre <- length(bre)
  }
  else {
    if(length(pct) > 1) {
      lo <- pct[1]
      up <- pct[2]
    }
    else {
      if(pct == 0) {
        lo <- min(pp)
        up <- max(pp)
      }
      else {
        v <- quantile(pp, c(pct, 1 - pct))
        lo <- v[1]
        up <- v[2]
      }
    }
  }
  ppp <- pmax(pmin(pp, up), lo)
  breaks <- seq(lo, up, length = bre)
  ph <- hist(ppp, breaks = breaks, plot = FALSE)

  # fitting curve 'f' to histogram:

  x <- (breaks[-1] + breaks[ - length(breaks)])/2
  y <- ph$counts
  N <- length(y)
  if(pct > 0) {
    y[1.] <- min(y[1.], 1.)
    y[N] <- min(y[N], 1.)
  }
  if(type == 0) {
    f <- glm(y ~ ns(x, df = df), poisson)$fit
  } 
  else {
    f <- glm(y ~ poly(x, df = df), poisson)$fit
  }

#  cat("The histogram:\n")
#  cat("x: ")
#  cat(x,"\n")
#  cat("y: ")
#  cat(y,"\n")
#  cat("The curve fitted to the full histogram using poisson regression, f:\n")
#  cat(round(f),"\n")
  
  # identifying local maxima and minima on curve 'f': 

  help=sign(f[1:(length(f)-1)]-f[2:(length(f))])
 
  i = 1
  signChanges=0
  numberSignChanges = 0
  localMaxima = 0
  numberLocalMaxima = 0
  localMinima = 0
  numberLocalMinima = 0

  if( help[1] > 0.0 ){
    numberLocalMaxima = numberLocalMaxima+1
    localMaxima[numberLocalMaxima] = 1
  }
  else{
    numberLocalMinima = numberLocalMinima+1
    localMinima[numberLocalMinima] = 1
  }
  while(i < (length(f)-1)){
    if( sign(help[i]) != sign(help[i+1]) ){
       if( sign(help[i]) < 0.0 ){
         numberLocalMaxima = numberLocalMaxima+1
         localMaxima[numberLocalMaxima] = i+1
       }
       else{
         numberLocalMinima = numberLocalMinima+1
         localMinima[numberLocalMinima] = i+1
       }
       numberSignChanges = numberSignChanges+1
       signChanges[numberSignChanges] = i+1
    } 
    i = i+1
  }
  if( sign(help[(length(f)-1)]) < 0.0 ){
    numberLocalMaxima = numberLocalMaxima+1
    localMaxima[numberLocalMaxima] = length(f)
  }
  else{
    numberLocalMinima = numberLocalMinima+1
    localMinima[numberLocalMinima] = length(f)
  }

  specialPoints=union(union(1,signChanges),length(f))

#  cat("Number localMaxima: ",numberLocalMaxima,"\n")
#  cat("localMaxima (red): x's:",x[localMaxima],"  y's:",f[localMaxima],"\n")
#  cat("Number localMinima: ",numberLocalMinima,"\n")
#  cat("localMinima (green): x's:",x[localMinima],"  y's:",f[localMinima],"\n")
#  cat("When fitting f0 we give zero weight to the categories below/above (both incl) these\n")
#  cat("two minimas: x:",x[localMinima[c(1,length(localMinima))]]," f:",f[localMinima[c(1,length(localMinima))]],"\n")
#  cat("apart from the first and last categories which are given count 0 and weight 1\n")
  innermaxval=max(f[setdiff(localMaxima,c(1,40))])
  innermaxentry=localMaxima[f[localMaxima]==innermaxval]
 
  # fitting two curves to histogram: one to the left (f0left) and one to the right
  # (f0right) of the inner max. These curves make up the null.

  df0=df

  # fitting the curve left of innermax:
  #   - first fixing count y[1] to zero and giving it weight 1:

  y <- ph$counts    
  y[1]=0
  weights0=rep(1,times=length(y))

  #   - giving weight 0 to categories between the first and the minimum right of first:

  if( localMinima[1] > 1 ){
    weights0[2:(localMinima[1])]=0
  }

 # cat("initial weights of categories for f0 fit:\n")
 # cat(weights0[1:innermaxentry],"\n")

  # - testing whether histogram counts are decreasing for categories right of
  #   'localMinima[1]' category. If they are the fit of f0 goes haywire. To avoid
  #   this we further fix the weights of the categories with decreasing counts immediately
  #   right of localMinima[1] to zero.

#  k=localMinima[1]+1

 # while( k<(innermaxentry-1) & ph$counts[k]>ph$counts[k+1]){
 #   weights0[k]=0
 #   k=k+1
 # }

#  cat("final weights of categories for f0 fit (after possibly setting weights of\n")
#  cat("additional categories to zero, in order to obtain required histogram shape for f0):\n")
#  cat(weights0[1:innermaxentry],"\n")

  if(type == 0) {
    f0left <- glm(y[1:innermaxentry] ~ ns(x[1:innermaxentry], df = df0),weights=weights0[1:innermaxentry], poisson)$fit
  } 
  else {
    f0left <- glm(y[1:innermaxentry] ~ poly(x[1:innermaxentry], df = df0),weights=weights0[1:innermaxentry], poisson)$fit
  }   

  NumberNullGenesf0left=round(sum(f0left),0)

  # fitting the curve right of innermax:
  #   - first fixing count y[length(y)] to zero and giving it weight 1: 

  y <- ph$counts    
  y[length(y)]=0
  weights0=rep(1,times=length(y))

  #   - giving weight 0 to categories between last and minimum left of last category

  if( localMinima[length(localMinima)] < length(y) ){
    weights0[localMinima[length(localMinima)]:(length(y)-1)]=0 
  }

  # testing whether histogram counts are decreasing for categories left of
  # 'localMinima[length(localMinima)]' category. If they are, the fit of f0 goes haywire. 
  # To avoid this we further fix the weights of the categories with decreasing counts immediately
  # left of localMinima[length(localMinima)] to zero.

#  k=localMinima[length(localMinima)]-1

#  while( k>(innermaxentry+1) & ph$counts[k-1]<ph$counts[k]){
#    weights0[k]=0
#    k=k+1
#  }

  if(type == 0) {
    f0right <- glm(y[innermaxentry:length(y)] ~ ns(x[innermaxentry:length(y)], df = df0),weights=weights0[innermaxentry:length(y)], poisson)$fit
  } 
  else {
    f0right <- glm(y[innermaxentry:length(y)] ~ poly(x[innermaxentry:length(y)], df = df0),weights=weights0[innermaxentry:length(y)], poisson)$fit
  }   
  
  NumberNullGenesf0right=round(sum(f0right[2:length(f0right)]),0)
 
  NumberNullGenesf0=min(NumberNullGenesf0left+NumberNullGenesf0right,sum(ph$counts))

  # taking care of roundings:

#  cat("the f0left curve:\n")
#  cat(round(f0left),"\n")
#  cat("the f0right curve:\n")
#  cat(round(f0right),"\n")

  NumberDiffGenesLeft=max(sum(ph$counts[1:innermaxentry])-NumberNullGenesf0left,0)
  NumberDiffGenesRight=max(sum(ph$counts[(innermaxentry+1):length(ph$counts)])-NumberNullGenesf0right,0)
  NumberDiffGenes=NumberDiffGenesLeft+NumberDiffGenesRight

  cat("number diff genes left: ",NumberDiffGenesLeft,"\n")
  cat("number diff genes right: ",NumberDiffGenesRight,"\n")
  cat("numberDiff genes: ",NumberDiffGenes,"\n")
  
  # plotting the histogram, fitted curve and intersting points:

  plot(ph,main="Choe: between",xlab=paste("#DEG:",NumberDiffGenes," up:",NumberDiffGenesLeft," down:",NumberDiffGenesRight,sep=""),cex.main=1.6,cex.lab=1.5,cex.axis=1.4)
  lines(x,f,col=8,lwd=2)
  points(x[localMaxima],f[localMaxima],col=2,pch=1,cex=2)  # maxima: red circle
  points(x[localMinima],f[localMinima],col=4,pch=2,cex=2)  # minima: blue
  points(x[innermaxentry],f[innermaxentry],col=1,pch=20,cex=2)
  # points(x[localMinima[c(1,length(localMinima))]],f[localMinima[c(1,length(localMinima))]],col=1,pch=20,cex=2)
  # points(x[innermaxentry],f[innermaxentry],col=1,pch=5,cex=2)
  #points(x[localMinima[c(1,length(localMinima))]],f[localMinima[c(1,length(localMinima))]],col=4,pch=5,cex=2)
  #points(x[postonegreflpoints],f[postonegreflpoints],col=1,pch=20,cex=2)
  #points(x[negtoposreflpoints],f[negtoposreflpoints],col=1,pch=20,cex=2)
  lines(x[1:innermaxentry],f0left,col=1,lwd=2)
  lines(x[innermaxentry:length(y)],f0right,col=1,lwd=2)

#  return(c(NumberDiffGenesLeft,NumberDiffGenesRight))

}

### helper functions
meanNorm <- function(mu, target=1) {
  if(length(mu)>1) for(i in c(1:length(mu))[!c(1:length(mu))==target]) mu[[i]]<- mu[[i]] * mean(mu[[target]])/mean(mu[[i]])
  return(mu)
}

loessNorm <- function(mu,target=1) {
  if(length(mu)>1) {
    for(i in c(1:length(mu))[!c(1:length(mu))==target]){
      muave <- matrix(nrow=nrow(mu[[1]]), ncol=2)
      muave[,1] <- apply(mu[[target]],1,mean)
      muave[,2] <- apply(mu[[i]],1,mean)
      muaveSubset <- data.frame(M = muave[,2] - muave[,1], A = 0.5*(muave[,1] + muave[,2]))
      muaveSubsetLoess <- loess(M~A, muaveSubset, span=0.1)
      predMuaveSubsetLoess <- predict(muaveSubsetLoess, data.frame(A=0.5*(muave[,1]+muave[,2])),se=FALSE)
      mu[[i]] <- mu[[i]] - predMuaveSubsetLoess 
    }
  }
  return(mu)
}

