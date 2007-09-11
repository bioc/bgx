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

 for(i in 1:nrow(summary)) { cat("***", row.names(summary)[i], "\n   ", paste(summary[i,1]),"\n")}
  
  summary <- as.vector(summary[[1]])
  noOfConditions <- as.numeric(summary[[2]])
  noOfGenes <- as.numeric(summary[[4]])
  subsample <- as.numeric(summary[[8]])
  iter <- as.numeric(summary[[10]])
  traceLength <- iter/subsample
  
  mu <- list() 
  for(m in 1:noOfConditions) {
    cat("Reading mu under condition ", m, "...", sep="")
    mu[[m]] <- matrix(scan(file.path(path, paste("mu.",m,sep=""))), ncol=traceLength)
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
plotDEDensity <- function(bgxOutput, gene=NULL, conditions=c(1,2),normalize=c("none","mean","loess"), normgenes=c(1:length(bgxOutput[["geneNames"]])), ...) {
  if(is.null(gene)) stop("Please specify a genes index/name.")
  if(is.character(gene)) gene <- (1:length(bgxOutput$geneNames))[gene == bgxOutput$geneNames]
  if(length(conditions) !=2 ) stop("Please specify exactly two conditions (e.g. conditions=c(1,2)).")
  normalize <- match.arg(normalize)
  if(normalize=="mean") bgxOutput$mu <- meanNorm(bgxOutput$mu, target=conditions[1])
  else if(normalize=="loess") bgxOutput$mu <- loessNorm(bgxOutput$mu, target=conditions[1], normgenes=normgenes)

  diff <- bgxOutput$mu[[conditions[2]]][gene,] - bgxOutput$mu[[conditions[1]]][gene,] 
  
  plot(density(diff, ...), main=paste("Density of mu for gene", bgxOutput$geneNames[gene], "\ncondition ", conditions[2] , " - condition", conditions[1]), xlab="Differential Expression")
  abline(v=0)

}

### Rank genes by deifferential expression
rankByDE <- function(bgxOutput, conditions=c(1,2),normalize=c("none", "mean", "loess"), normgenes=c(1:length(bgxOutput[["geneNames"]])), absolute=TRUE) {
  normalize <- match.arg(normalize)
  if(normalize=="mean") bgxOutput$mu <- meanNorm(bgxOutput$mu, target=conditions[1])
  else if(normalize=="loess") bgxOutput$mu <- loessNorm(bgxOutput$mu, target=conditions[1],normgenes=normgenes)

  mu_diff <- bgxOutput$mu[[conditions[2]]] - bgxOutput$mu[[conditions[1]]]
  de <- c()
  tau <- var <- m <- 0
  for(i in 1:nrow(bgxOutput$mu[[1]])) {
    temp <- mu_diff[i,]

    # using sokal function
    var <- tau <- 0.0
    m <- 0
    sok <- .C("sokal", as.integer(1024), as.double(temp), as.double(var), as.double(tau), as.integer(m), PACKAGE="bgx")
    mcse <- sqrt(1023*sok[[3]]*sok[[4]]/(1024*(1024-sok[[4]])))

    # slow
    #tau <- 1
    #a <- acf(temp, plot=FALSE)
    #maxlag <- max(a$lag)
    #for(j in 2:maxlag) tau <- tryCatch(tau + 2 * (1-a$lag[j]/maxlag) * a$acf[j], error= function(e) NaN)
    #mcse <- sqrt((a$n.used-1)*tau*var(temp)/(a$n.used*(a$n.used-tau)))
    
    if(absolute) de[i] <- abs(mean(mu_diff[i,]))/mcse
    else de[i] <- mean(mu_diff[i,])/mcse
  }

  de[which(is.nan(de))] <- max(de, na.rm=T) # Get rid of NaNs returned in sokal()
  order <- sort(de,decreasing=TRUE, index.return=TRUE)$ix

  return(matrix(c(order,de[order]), dimnames=list( c(bgxOutput$geneNames[order]), c("Position", "DiffExpression")), ncol=2))
}

### Plot sorted 2.5-97.5% CI of DE (not absolute value)
plotDiffRank <- function(bgxOutput, conditions=c(1,2),normalize=c("none", "mean", "loess"), normgenes=c(1:length(bgxOutput[["geneNames"]])), ymax=NULL) {
  normalize <- match.arg(normalize)
  if(normalize=="mean") bgxOutput$mu <- meanNorm(bgxOutput$mu, target=conditions[1])
  else if(normalize=="loess") bgxOutput$mu <- loessNorm(bgxOutput$mu, target=conditions[1], normgenes=normgenes)

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
plotDEHistogram <- function(bgxOutput, conditions=c(1,2), normalize=c("none", "mean", "loess"), 
  normgenes=c(1:length(bgxOutput[["geneNames"]])), df=floor(1.8*log10(length(bgxOutput[["geneNames"]])))) {
  normalize <- match.arg(normalize)
  if(normalize=="mean") bgxOutput$mu <- meanNorm(bgxOutput$mu, target=conditions[1])
  else if(normalize=="loess") bgxOutput$mu <- loessNorm(bgxOutput$mu, target=conditions[1], normgenes)

  diff <- bgxOutput$mu[[conditions[2]]] - bgxOutput$mu[[conditions[1]]]
  pp <- c()
  for(i in 1:length(bgxOutput$geneNames)) pp[i] <- mean(diff[i,] < 0)

  if(length(pp) < 100) stop("Sorry, your sample of genes is too small. This method only works well with 100 or more genes.\n")

  cat("Degrees of freedom: ",df, "(decrease for smoother curve).\n")

#  require("splines")
  ppp <- pmax(pmin(pp, 1), 0)
  breaks <- seq(0, 1, by=1/ifelse(length(pp)<1000,20,40))
  ph <- hist(ppp, breaks = breaks, plot = FALSE)

  # fitting curve 'f' to histogram:

  x <- (breaks[-1] + breaks[ - length(breaks)])/2
  y <- ph$counts
  N <- length(y)
  f <- glm(y ~ poly(x, df = df), poisson)$fit

  # identifying local maxima and minima on curve 'f': 

  help=sign(f[1:(length(f)-1)]-f[2:(length(f))])
 
  i = 1
  signChanges=0
  numberSignChanges = 0
  localMaxima = 0
  numberLocalMaxima = 0
  localMinima = 0
  numberLocalMinima = 0

  if( help[1] > 0.0 ) {
    numberLocalMaxima = numberLocalMaxima+1
    localMaxima[numberLocalMaxima] = 1
  } else {
    numberLocalMinima = numberLocalMinima+1
    localMinima[numberLocalMinima] = 1
  }
  while(i < (length(f)-1)){
    if( sign(help[i]) != sign(help[i+1]) ){
       if( sign(help[i]) < 0.0 ){
         numberLocalMaxima = numberLocalMaxima+1
         localMaxima[numberLocalMaxima] = i+1
       } else{
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
  } else{
    numberLocalMinima = numberLocalMinima+1
    localMinima[numberLocalMinima] = length(f)
  }

  specialPoints=union(union(1,signChanges),length(f))

  innermaxval=max(f[setdiff(localMaxima,c(1,40))])
  innermaxentry=localMaxima[f[localMaxima]==innermaxval]
 
  # fitting two curves to histogram: one to the left (f0left) and one to the right
  # (f0right) of the inner max. These curves make up the null.

  # fitting the curve left of innermax:
  #   - first fixing count y[1] to zero and giving it weight 1:

  y <- ph$counts    
  y[1]=0
  weights0=rep(1,times=length(y))

  #   - giving weight 0 to categories between the first and the minimum right of first:

  if( localMinima[1] > 1 ){
    weights0[2:(localMinima[1])]=0
  }

  # Abort if the innermaxentry is too near the edges
  if(length(1:innermaxentry) < df || length(innermaxentry:length(y)) < df) stop("Sorry, differential expression appears to be systematically biased. Maybe you need to normalize (see help(\"plotDEHistogram\")) or maybe you are not using a representative subset of genes.\n")

  f0left <- glm(y[1:innermaxentry] ~ poly(x[1:innermaxentry], df = df),weights=weights0[1:innermaxentry], poisson)$fit

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

  f0right <- glm(y[innermaxentry:length(y)] ~ poly(x[innermaxentry:length(y)], df = df),weights=weights0[innermaxentry:length(y)], poisson)$fit
  
  NumberNullGenesf0right=round(sum(f0right[2:length(f0right)]),0)
 
  NumberNullGenesf0=min(NumberNullGenesf0left+NumberNullGenesf0right,sum(ph$counts))

  # taking care of roundings:

  NumberDiffGenesLeft=max(sum(ph$counts[1:innermaxentry])-NumberNullGenesf0left,0)
  NumberDiffGenesRight=max(sum(ph$counts[(innermaxentry+1):length(ph$counts)])-NumberNullGenesf0right,0)
  NumberDiffGenes=NumberDiffGenesLeft+NumberDiffGenesRight

  cat("Number of differentially expressed  genes (left): ",NumberDiffGenesLeft,"\n")
  cat("Number of differentially expressed  genes (right): ",NumberDiffGenesRight,"\n")
  cat("Total number of differentially expressed genes: ",NumberDiffGenes,"\n")
  
  # plotting the histogram, fitted curve and intersting points:

  plot(ph,main=substitute(paste("Histogram of P(",mu[two] - mu[one]," < 0); df=",df),list(two=paste("g",conditions[2],sep=""),one=paste("g",conditions[1],sep=""),df=df)) ,xlab=paste("#DEG:",NumberDiffGenes," up:",NumberDiffGenesLeft," down:",NumberDiffGenesRight,sep=""))
  lines(x,f,col=8,lwd=2)
  points(x[localMaxima],f[localMaxima],col=2,pch=1,cex=2)  # maxima: red circle
  points(x[localMinima],f[localMinima],col=4,pch=2,cex=2)  # minima: blue
  points(x[innermaxentry],f[innermaxentry],col=1,pch=20,cex=2)
  lines(x[1:innermaxentry],f0left,col=1,lwd=2)
  lines(x[innermaxentry:length(y)],f0right,col=1,lwd=2)

#  return(c(NumberDiffGenesLeft,NumberDiffGenesRight))

}

### helper functions
meanNorm <- function(mu, target=1) {
  if(length(mu)>1) for(i in c(1:length(mu))[!c(1:length(mu))==target]) mu[[i]]<- mu[[i]] * mean(mu[[target]])/mean(mu[[i]])
  return(mu)
}

loessNorm <- function(mu,target=1, normgenes=c(1:nrow(mu[[target]]))) {
  if(length(mu)>1) {
    for(i in c(1:length(mu))[!c(1:length(mu))==target]){
      muave <- matrix(nrow=nrow(mu[[1]]), ncol=2)
      muave[,1] <- apply(mu[[target]],1,mean)
      muave[,2] <- apply(mu[[i]],1,mean)
      muaveSubset <- data.frame(M = muave[normgenes,2] - muave[normgenes,1], A = 0.5*(muave[normgenes,1] + muave[normgenes,2]))
      muaveSubsetLoess <- loess(M~A, muaveSubset, span=0.1)
      predMuaveSubsetLoess <- predict(muaveSubsetLoess, data.frame(A=0.5*(muave[,1]+muave[,2])),se=FALSE)
      mu[[i]] <- mu[[i]] - predMuaveSubsetLoess 
    }
  }
  return(mu)
}

