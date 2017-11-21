## Does the modeling on abundances, Shannon diversity and MDS axes
## parametric and nonparametric
## And the Shannon diversity
## And modeling on MDS axes
## bestAt doesn't work for RDP because of duplicate names across taxa
rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")
library("Kendall")
library("coin")

## TODO: Make this more generic
## TODO: Add in Kendall for the "Age" variable (and others)
baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging"
setwd(baseDir)

## dataType <- "closedQIIMER1"
dataType <- "RDPR1"
metadataDir <- paste(baseDir, "metadata", sep="/")
dataDir <- paste(baseDir, dataType, sep="/")
processedDir <- paste(dataDir, "processed", sep="/")
analysisDir <- paste(baseDir, dataType, "analysis", sep="/")

## baseDataFileName <- "closed_reference_otu_table_L"
baseDataFileName <- "Samples.classified_"

## It would be great to make the name change the variables in place
## This should replace the each modeling code
#filePrefix <- "Age"
#filePrefix <- "Mucosal.Gene.Ct"
filePrefix <- "Group"
#filePrefix <- "Fecal.Gene.Ct"

## taxaLevels <- c(2:7)
taxaLevels <- c("phylum", "class", "order", "family", "genus")
bestAt <- data.frame()
for (taxa in taxaLevels){
  setwd(processedDir)
  inFileName <- paste(baseDataFileName, taxa, "_LogNormwithMetadata.txt", sep="")
  myT <- read.table(inFileName, header=TRUE, sep="\t")
  numCols <- ncol(myT)
  
  ## endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
  ## The 9 could be done better...
  myColClasses <- metaColClasses<- c("character", "character", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", rep("numeric", numCols - 9))
  myT <- read.table(inFileName, header=TRUE, sep="\t", colClasses = myColClasses)
  
  endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
  setwd(analysisDir)
  ## Need to find a way to do all tissues
  ## No, that is another scripts job as it involves a second variable.
  if (filePrefix == "Fecal.Gene.Ct") {
    myT <- myT[complete.cases(myT),]
  }
  tissues <- unique(myT$Sample.Type)
  for(tissue in tissues){
    subT <- myT[myT$Sample.Type == tissue,]
    
    names <- vector()
    pValuesfilePrefix<- vector()
    pValuesfilePrefixWilcox<- vector()
    pValuesCage<- vector()
    iccCage <- vector()
    pValuesSampleType <- vector()
    pValuesAge<- vector()
    pValuesAnimal<- vector()
    pValuesInteraction <- vector()
    allpvals <- list()
    ## Additional checks
    pValuesDepth <- vector()
    ## pValuesDepthWilcox <- vector()
    pValuesDepthwithfilePrefix <- vector()
    pValuesfilePrefixwithDepth <- vector()
    
    index <- 1
    
    ## Over the non-metadata columns
    ## Will replace with a function call without the loop
    ## Add in some timing information.
    for ( j in c((endMetadataIndex + 1):ncol(subT))) {
      bug <- subT[,j]
      ## Removing rare taxa
      ## This should be a bit more flexible
      if( sum(bug != 0 ) > nrow(subT) / 4 ) {
        names[index] <- names(subT)[j]

        ifelse(class(subT[filePrefix][,1]) == "character",
               ##class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character",
               varCare <- factor(as.numeric(as.factor(unlist(subT[filePrefix])))),
               varCare <- as.numeric(as.factor(unlist(subT[filePrefix])))
        )
        sampleType <- factor(subT$Sample.Type)
        
        cage <-  factor(subT$Pen.Location)
        depth <- subT$sampleDepth
        ## depth <- subT$depthAtLevel
        
        myFrame <- data.frame(bug, varCare, sampleType, depth, cage)
        
        varCareModel <- lm(bug~varCare, data= myFrame)
        #print("Didn't get through the first one")
        depthModel <- lm(bug~depth, data = myFrame)
        depthWithCareModel <- lm(bug~varCare + depth, data=myFrame)
        mixedAnova <- anova(varCareModel)
        ## This method of getting p-values only works for the univariate case
        pValuesfilePrefix[index] <- mixedAnova$"Pr(>F)"[1]
        ## ifelse( class(subT[filePrefix][,1]) == "character",
        ##class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character",
        ## pValuesfilePrefixWilcox[index] <- wilcox.test( bug[subT[filePrefix] == "Young"],    bug[subT[filePrefix] == "Old"])$p.value## ,
        ## This might need to become more complicated
        if(class(varCare) == "numeric") {
          pValuesfilePrefixWilcox[index] <- cor.test( bug, varCare, method="kendall")$p.value[1]
        } else {
          pValuesfilePrefixWilcox[index] <- pvalue(wilcox_test(bug~varCare, data=myFrame))
        }
        
        ##pValuesfilePrefixWilcox[index] <- wilcox_test( bug[subT[filePrefix] == "Young"]~bug[subT[filePrefix] == "Old"], )$p.value,
        ## So this shouldn't actually appear
        ## pValuesfilePrefixWilcox[index] <- cor.test(bug, as.numeric(as.factor(unlist(subT[filePrefix]))), method="kendall")$p.value[1]
        ## )
        pValuesDepth[index] <- anova(depthModel)$"Pr(>F)"[1]
        ## pValuesDepthwithfilePrefix[index] <-
        pValuesfilePrefixwithDepth[index] <- anova(depthWithCareModel, depthModel)$"Pr(>F)"[2]
        
        ## Model for depth
        ## Model for varCare and depth
        
        index <- index + 1
      }
    }
    
    dFrame <- data.frame(taxa, tissue, names, pValuesfilePrefix, pValuesfilePrefixWilcox, pValuesDepth, pValuesfilePrefixwithDepth)
    dFrame <- dFrame [order(dFrame$pValuesfilePrefix),]
    rownames(dFrame) <- dFrame$names
    dFrame$adjustedpValuesfilePrefix <- p.adjust( dFrame$pValuesfilePrefix, method = "BH" )
    dFrame$adjustedpValuesfilePrefixWilcox <- p.adjust( dFrame$pValuesfilePrefixWilcox, method = "BH" )
    dFrame$adjustedpValuesDepth <- p.adjust(dFrame$pValuesDepth, method="BH")
    dFrame$adjustedpValuesfilePrefixwithDepth <- p.adjust(dFrame$pValuesfilePrefixwithDepth, method="BH")
    oldColnames <- colnames(dFrame)
    ## This can probably be better written
    if(class(varCare) == "numeric") {
      colnames(dFrame) <- c("TaxonomicLevel", "tissue", "names", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Kendall"), paste0("p-valuesDepth"), paste0("p-valuesDepthwith", filePrefix), paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Kendall"), paste0("adjusted_p-valuesDepth"), paste0("adjusted_p-valuesDepthwith", filePrefix))
      
    } else {
      colnames(dFrame) <- c("TaxonomicLevel", "tissue", "names", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Wilcox"), paste0("p-valuesDepth"), paste0("p-valuesDepthwith", filePrefix), paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Wilcox"), paste0("adjusted_p-valuesDepth"), paste0("adjusted_p-valuesDepthwith", filePrefix))
    }
    write.table(dFrame, file=paste(dataType, "_L_", taxa, "_", tissue, "_", filePrefix, ".txt",sep=""), sep="\t", row.names = FALSE)
    ## Save results to temporary before merging across taxonomic levels if they are below the 0.05
    bestAt <- rbind(bestAt, dFrame[which(dFrame[grep("adjust", colnames(dFrame))] < 0.10, arr.ind=TRUE)[,1],])
    colnames(dFrame) <- oldColnames
    ## This plots out everything. I'll need another loop to do just the significant ones.
    pdf( paste(dataType, "_L_", taxa, "_", tissue, "_", filePrefix, "_boxplots.pdf", sep=""))
    index <- 1
    for (nameIter in rownames(dFrame)) {
      ## This should be stored once near the top and altered there
      par(mfrow=c(3,1),
          oma = c(1, 3, 2, 0) + 0.1,
          mar = c(3, 2, 3, 1) + 0.1,
          mgp = c(2, 0.75, 0))
      
      bug <- subT[,nameIter]
      ## ifelse(class(subT[filePrefix][,1]) == "character",
      ##        ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character",
      ##        varCare <- factor(as.numeric(as.factor(unlist(subT[filePrefix])))),
      ##        varCare <- as.numeric(as.factor(unlist(subT[filePrefix])))
      ##        )
      
      ## depth <- subT$depthAtLevel
      varCare <- subT[filePrefix][[1]]
      depth <- subT$sampleDepth
      sampleType <- factor(subT$Sample.Type)
      cage <-  factor(subT$Pen.Location)
      myFrame <- data.frame(bug, varCare, sampleType, cage, depth)
      ## Rerunning the model to test for the visuals
      ## varCareModel <- lm(bug~varCare, data= myFrame)
      ## depthModel <- lm(bug~depth, data = myFrame)
      ## depthWithCareModel <- lm(bug~varCare + depth, data=myFrame)
      
      if(##class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character"
        class(subT[filePrefix][,1]) == "character"){
        
        boxplot( bug ~ varCare , ylab="Log normalized abundance", pch=19, xlab=filePrefix,
                 main = paste(filePrefix," linear model p-value", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), " Wilcoxon test p-value", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ), outline=FALSE)
        
        stripchart(bug ~ varCare ,
                   data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "Log normalized abundance", method="jitter", col=c("red", "blue"))
        ## axis(1, at=c(1,2), labels=c("Old", "Young")##, col=c("red", "blue"))
      } else {
        plot( bug ~ varCare , ylab="Log normalized abundance", pch=19, xlab=filePrefix,
              main = paste(filePrefix," linear model p-value ", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), " Kendall correlation test p-value ", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) )
      }
      ## Turning off the coloration so you know it only cares about depth
      plot( bug ~ depth , ylab="Log normalized abundance", pch=19,
            main = paste("Depth linear model p-value ", format(dFrame$adjustedpValuesDepth[index],digits=3)), xlab="Sequences in Sample")##, col=ifelse(varCare == "Young", "blue", "red"))
      
      stripchart(bug ~ depth ,
                 data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "Log normalized abundance", method="jitter", xlab="Sequences in Sample", col=c("red", "blue"))
      ## title(xlab = "Sequences in Sample")
      ## axis(1, at=c(1,2), labels=c("Old", "Young")##, col=c("red", "blue"))
      if (## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character"
        class(subT[filePrefix][,1]) == "character"){
        plot(depth[varCare == "Young"], bug[varCare == "Young"], ylab="Log normalized abundance", xlab="Sequencing depth at this taxonomic Level", pch=19,
             main = paste(filePrefix, " and Depth linear model p-value ", format(dFrame$adjustedpValuesfilePrefixwithDepth[index],digits=3)), col="blue",
             xlim = c(min(depth), max(depth)), ylim = c(min(bug), max(bug)))
        points(depth[varCare == "Old"], bug[varCare == "Old"], col="red", pch=19)
        
        deptho <- depth[order(depth)]
        varo <- varCare[order(depth)]
        bugo <- bug[order(depth)]
        lines(deptho[varo == "Young"], bugo[varo == "Young"], col = "blue")
        lines(deptho[varo == "Old"], bugo[varo == "Old"], col = "red")
      } else{
        ## TODO: Need a new plot here!
        plot.new()
      }
      
      ## plot( bug ~ varCare + depth , ylab="Log normalized abundance", pch=19,
      ##      main = paste(filePrefix, "and Depth p-value linear model", format(dFrame$adjustedpValuesfilePrefixwithDepth[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))
      
      ## stripchart(bug ~ varCare + depth ,
      ##            data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "Log normalized abundance", method="jitter")
      
      mtext(paste(taxa, nameIter), outer=TRUE, cex = 1.0)
      mtext("Log normalized abundance", outer=TRUE, side=2, line=1)
      
      index = index + 1
    }
    
    hist(pValuesfilePrefix, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", filePrefix, " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    
    if(class(varCare) == "numeric") {
      hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Kendall-tau correlation test) for\n", filePrefix, " at Taxonomic Level:", taxa),
           ylab = "Count", xlab="Unadjusted p-values")
    } else {
      hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Wilcoxon test) for\n", filePrefix, " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    }
    hist(pValuesDepth, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", "Depth", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesfilePrefixwithDepth, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", filePrefix, " after depth", " at Taxonomic Level:", taxa),       ylab = "Count", xlab="Unadjusted p-values")
    
    dev.off()
  }
}

bestAt <- unique(bestAt)
write.table(bestAt, file=paste0("SigHits_",filePrefix,".txt"), sep="\t", row.names=FALSE)

pdf( paste(dataType, "_", filePrefix, "_SIGHITboxplots.pdf", sep=""))
par(mfrow=c(3,1),
    oma = c(1, 3, 2, 0) + 0.1,
    mar = c(3, 2, 3, 1) + 0.1,
    mgp = c(2, 0.75, 0))

colnames(bestAt) <- oldColnames
for (taxa in taxaLevels){
  setwd(processedDir)
  inFileName <- paste(baseDataFileName, taxa, "_LogNormwithMetadata.txt", sep="")
  myT <- read.table(inFileName, header=TRUE, sep="\t")
  numCols <- ncol(myT)
  
  ## endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
  ## The 9 could be done better...
  myColClasses <- metaColClasses<- c("character", "character", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", rep("numeric", numCols - 9))
  myT <- read.table(inFileName, header=TRUE, sep="\t", colClasses = myColClasses)
  
  endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
  setwd(analysisDir)
  
  if (filePrefix == "Fecal.Gene.Ct") {
    myT <- myT[complete.cases(myT),]
  }
  
  ## pdf( paste(dataType, "_L_", taxa, "_", tissue, "_", filePrefix, "_SIGHITboxplots.pdf", sep=""))
  index <- 1
  ## for (nameIter in rownames(bestAt)) {
  for (nameIter in bestAt$names) {
    ## This should be stored once near the top and altered there
    ## par(mfrow=c(3,1),
    ##     oma = c(0, 3, 2, 0) + 0.1,
    ##     mar = c(3, 2, 3, 1) + 0.1,
    ##     mgp = c(3, 1.25, 0))
    subT <- myT[myT$Sample.Type == bestAt[nameIter,]$tissue,]
    if (bestAt[nameIter,]$taxa == taxa) {
      bug <- subT[,nameIter]
      ## ifelse(class(subT[filePrefix][,1]) == "character",
      ##        ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character",
      ##        varCare <- factor(as.numeric(as.factor(unlist(subT[filePrefix])))),
      ##        varCare <- as.numeric(as.factor(unlist(subT[filePrefix])))
      ##        )
      
      ## depth <- subT$depthAtLevel
      varCare <- subT[filePrefix][[1]]
      depth <- subT$sampleDepth
      sampleType <- factor(subT$Sample.Type)
      cage <-  factor(subT$Pen.Location)
      myFrame <- data.frame(bug, varCare, sampleType, cage, depth)
      ## Looking over the models for problems
      varCareModel <- lm(bug~varCare, data= myFrame)
      #print("Didn't get to not sure where we are")
      depthModel <- lm(bug~depth, data = myFrame)
      depthWithCareModel <- lm(bug~varCare + depth, data=myFrame)
      
      if(class(subT[filePrefix][,1]) == "character"
         ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character"
      ){
        boxplot( bug ~ varCare , ylab="Log normalized abundance", pch=19,
                 main = paste(filePrefix," linear model p-value ", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), " Wilcoxon test p-value ", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) , outline=FALSE)
        
        stripchart(bug ~ varCare ,
                   data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "Log normalized abundance", method="jitter", col=c("red", "blue"))
      } else {
        plot( bug ~ varCare , ylab="Log normalized abundance", pch=19,
              main = paste(filePrefix," linear model p-value ", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), " Kendall-tau correlation test p-value ", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) )
      }
      
      hist(residuals(varCareModel), col="darkgrey")
      plot(fitted(varCareModel),
           residuals(varCareModel))
      
      plot( bug ~ depth , ylab="Log normalized abundance", pch=19,
            main = paste("Depth linear model p-value ", format(bestAt[nameIter,]$adjustedpValuesDepth,digits=3)), col=ifelse(varCare == "Young", "blue", "red"))
      
      stripchart(bug ~ depth ,
                 data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "Log normalized abundance", method="jitter", col=c("red", "blue"))
      hist(residuals(depthModel), col="darkgrey")
      plot(fitted(depthModel),
           residuals(depthModel))
      
      ## plot( bug ~ varCare + depth , ylab="Log normalized abundance", pch=19,
      ##     main = paste(filePrefix, "and Depth p-value linear model", format(bestAt[nameIter,]$adjustedpValuesfilePrefixwithDepth,digits=3)), col=ifelse(varCare == "Young", "blue", "red"))
      
      if( class(subT[filePrefix][,1]) == "character"
          ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character"
      ){
        
        plot(depth[varCare == "Young"], bug[varCare == "Young"], ylab="Log normalized abundance", xlab="Sequencing depth at this taxonomic Level", pch=19,
             main = paste(filePrefix, " and Depth linear model p-value", format(bestAt[nameIter,]$adjustedpValuesfilePrefixwithDepth,digits=3)), col="blue",
             xlim = c(min(depth), max(depth)), ylim = c(min(bug), max(bug)))
        points(depth[varCare == "Old"], bug[varCare == "Old"], col="red", pch=19)
        
        deptho <- depth[order(depth)]
        varo <- varCare[order(depth)]
        bugo <- bug[order(depth)]
        lines(deptho[varo == "Young"], bugo[varo == "Young"], col = "blue")
        lines(deptho[varo == "Old"], bugo[varo == "Old"], col = "red")
      } else{
        ## TODO: Need a new plot here!
        plot.new()
      }
      
      ## stripchart(bug ~ varCare + depth ,
      ##            data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "Log normalized abundance", method="jitter")
      hist(residuals(depthWithCareModel), col="darkgrey")
      plot(fitted(depthWithCareModel),
           residuals(depthWithCareModel))
      
      ## Try to do this iteratively over the models
      ## Plotting to check for normality of the residuals
      
      ## Plotting of residuals versus predicted values
      
      mtext(paste(taxa, nameIter), outer=TRUE, cex = 1.0)
      mtext("Log normalized abundance", outer=TRUE, side=2, line=1)
    }
    
    ## index = index + 1
  }
}
dev.off()
## This is for the MDS analysis
for (taxa in taxaLevels){
  setwd(processedDir)
  inFileName <- paste(baseDataFileName, taxa, "_LogNormwithMetadata.txt", sep="")
  myT <- read.table(inFileName, header=TRUE, sep="\t")
  numCols <- ncol(myT)
  
  ## endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
  ## The 9 could be done better...

  myColClasses <- metaColClasses<- c("character", "character", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", rep("numeric", numCols - 9))
  myT <- read.table(inFileName, header=TRUE, sep="\t", colClasses = myColClasses)
  
  endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
  setwd(analysisDir)
  if (filePrefix == "Fecal.Gene.Ct") {
    myT <- myT[complete.cases(myT),]
    tissues <- c("Feces")
  } else {
    tissues <- unique(myT$Sample.Type)
  }
  ## Need to find a way to do all tissues
  ## No, that is another scripts job as it involves a second variable.
  #tissues <- unique(myT$Sample.Type)
  ##tissues <- c("Feces")
  for(tissue in tissues){
    subT <- myT[myT$Sample.Type == tissue,]
    inFileName <- paste(baseDataFileName, taxa, "_", gsub('([[:punct:]])|\\s+','_',tissue), "_pcoa.txt",sep="")
    ## This will not work the first time through
    myPCoA <- read.table(inFileName, sep="\t", header=TRUE)
    inFileName2 <- paste(baseDataFileName, taxa, "_", gsub('([[:punct:]])|\\s+','_',tissue), "_eigenValues.txt",sep="")
    myEigs <- read.table(inFileName2, sep="\t", header=TRUE)
    myPercent <- myEigs/sum(myEigs)
    names <- vector()
    pValuesfilePrefix<- vector()
    pValuesfilePrefixWilcox<- vector()
    pValuesCage<- vector()
    iccCage <- vector()
    pValuesSampleType <- vector()
    pValuesAge<- vector()
    pValuesAnimal<- vector()
    pValuesInteraction <- vector()
    allpvals <- list()
    ## Additional checks
    pValuesDepth <- vector()
    ## pValuesDepthWilcox <- vector()
    pValuesDepthwithfilePrefix <- vector()
    pValuesfilePrefixwithDepth <- vector()
    
    index <- 1
    
    ## Over the non-metadata columns
    ## Will replace with a function call without the loop
    ## Add in some timing information.
    for ( j in 1:10) {
      bug <- myPCoA[,j]
      ## Removing rare taxa
      ## This should be a bit more flexible
      names[index] <- names(myPCoA)[j]
      
      ## ifelse(class(subT[filePrefix][,1]) == "character",
      ##        ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character",
      ##        varCare <- factor(as.numeric(as.factor(unlist(subT[filePrefix])))),
      ##        varCare <- as.numeric(as.factor(unlist(subT[filePrefix])))
      ##        )
      varCare <- subT[filePrefix][[1]]
      sampleType <- factor(subT$Sample.Type)
      
      cage <-  factor(subT$Pen.Location)
      ## depth <- subT$depthAtLevel
      depth <- subT$sampleDepth
      
      myFrame <- data.frame(bug, varCare, sampleType, depth, cage)
      
      varCareModel <- lm(bug~varCare, data= myFrame)
      #print("Doesn't get here MDS")
      depthModel <- lm(bug~depth, data = myFrame)
      depthWithCareModel <- lm(bug~varCare + depth, data=myFrame)
      mixedAnova <- anova(varCareModel)
      ## This method of getting p-values only works for the univariate case
      pValuesfilePrefix[index] <- mixedAnova$"Pr(>F)"[1]
      ## ifelse(class(subT[filePrefix][,1]) == "character",
      ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character",
      ##pValuesfilePrefixWilcox[index] <- wilcox.test( bug[subT[filePrefix] == "Young"],    bug[subT[filePrefix] == "Old"])$p.value##,
      
      if(class(varCare) == "numeric") {
        pValuesfilePrefixWilcox[index] <- cor.test( bug, varCare, method="kendall")$p.value[1]
      } else {
        pValuesfilePrefixWilcox[index] <- pvalue(wilcox_test(bug~varCare, data=myFrame))
      }
      ## pValuesfilePrefixWilcox[index] <- cor.test(bug, as.numeric(as.factor(unlist(subT[filePrefix]))), method="kendall")$p.value[1]
      ##   )
      
      pValuesDepth[index] <- anova(depthModel)$"Pr(>F)"[1]
      ## pValuesDepthwithfilePrefix[index] <-
      pValuesfilePrefixwithDepth[index] <- anova(depthWithCareModel, depthModel)$"Pr(>F)"[2]
      
      ## Model for depth
      ## Model for varCare and depth
      
      index <- index + 1
    }
    
    dFrame <- data.frame(names, myPercent[1:length(names),], pValuesfilePrefix, pValuesfilePrefixWilcox, pValuesDepth, pValuesfilePrefixwithDepth)
    
    ##dFrame <- dFrame [order(dFrame$pValuesfilePrefix),]
    rownames(dFrame) <- dFrame$names
    dFrame$adjustedpValuesfilePrefix <- p.adjust( dFrame$pValuesfilePrefix, method = "BH" )
    dFrame$adjustedpValuesfilePrefixWilcox <- p.adjust( dFrame$pValuesfilePrefixWilcox, method = "BH" )
    dFrame$adjustedpValuesDepth <- p.adjust(dFrame$pValuesDepth, method="BH")
    dFrame$adjustedpValuesfilePrefixwithDepth <- p.adjust(dFrame$pValuesfilePrefixwithDepth, method="BH")
    oldColnames <- colnames(dFrame)
    ## This can probably be better written
    if (class(varCare) =="numeric"){
      colnames(dFrame) <- c("names", "percent_variance", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Kendall"), paste0("p-valuesDepth"), paste0("p-valuesDepthwith", filePrefix), paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Kendall"), paste0("adjusted_p-valuesDepth"), paste0("adjusted_p-valuesDepthwith", filePrefix))
    } else {
    colnames(dFrame) <- c("names", "percent_variance", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Wilcox"), paste0("p-valuesDepth"), paste0("p-valuesDepthwith", filePrefix), paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Wilcox"), paste0("adjusted_p-valuesDepth"), paste0("adjusted_p-valuesDepthwith", filePrefix))
    }
    dFrame <- dFrame [order(dFrame$"percent_variance", decreasing = TRUE),]
    
    write.table(dFrame, file=paste(dataType, "_L_", taxa, "_", tissue, "_", filePrefix, "MDS.txt",sep=""), sep="\t", row.names = FALSE)
    colnames(dFrame) <- oldColnames
    pdf( paste(dataType, "_L_", taxa, "_", tissue, "_", filePrefix, "_MDSboxplots.pdf", sep=""))
    index <- 1
    for (j in 1:10) {
      ## This should be stored once near the top and altered there
      par(mfrow=c(3,1),
          oma = c(1, 3, 2, 0) + 0.1,
          mar = c(3, 2, 3, 1) + 0.1,
          mgp = c(2, 0.75, 0))
      
      
      bug <- myPCoA[,j]
      ## ifelse(class(subT[filePrefix][,1]) == "character",
      ##        ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character",
      ##        varCare <- factor(as.numeric(as.factor(unlist(subT[filePrefix])))),
      ##        varCare <- as.numeric(as.factor(unlist(subT[filePrefix])))
      ##        )
      
      ## depth <- subT$depthAtLevel
      varCare <- subT[filePrefix][[1]]
      depth <- subT$sampleDepth
      sampleType <- factor(subT$Sample.Type)
      cage <-  factor(subT$Pen.Location)
      myFrame <- data.frame(bug, varCare, sampleType, cage, depth)
      
      if(class(subT[filePrefix][,1]) == "character"
         ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character"
      ){
        boxplot( bug ~ varCare , ylab="Log normalized abundance", pch=19,
                 main = paste(filePrefix, " linear model p-value ", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), " Wilcoxon test p-value ", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) , outline=FALSE)
        
        stripchart(bug ~ varCare ,
                   data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "Log normalized abundance", method="jitter", col=c("red", "blue"))
      } else {
        plot( bug ~ varCare , ylab="Log normalized abundance", pch=19,
              main = paste(filePrefix, " linear model p-value", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), " Kendall-tau correlation test p-value ", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) )
      }
      
      plot( bug ~ depth , ylab="MDS Axis", pch=19,
            main = paste("Depth linear model p-value ", format(dFrame$adjustedpValuesDepth[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))
      
      if(class(subT[filePrefix][,1]) == "character"
         ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character"
      ){
        
        plot(depth[varCare == "Young"], bug[varCare == "Young"], ylab="Log normalized abundance", xlab="Sequencing depth at this taxonomic Level", pch=19,
             main = paste(filePrefix, " and Depth  linear model p-value ", format(dFrame$adjustedpValuesfilePrefixwithDepth[index],digits=3)), col="blue",
             xlim = c(min(depth), max(depth)), ylim = c(min(bug), max(bug)))
        points(depth[varCare == "Old"], bug[varCare == "Old"], col="red", pch=19)
        
        deptho <- depth[order(depth)]
        varo <- varCare[order(depth)]
        bugo <- bug[order(depth)]
        lines(deptho[varo == "Young"], bugo[varo == "Young"], col = "blue")
        lines(deptho[varo == "Old"], bugo[varo == "Old"], col = "red")
      } else{
        ## TODO: Need a new plot here!
        plot.new()
      }
      
      ## stripchart(bug ~ varCare + depth ,
      ##           data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "MDS Axis", method="jitter")
      
      ## mtext(nameIter, outer=TRUE, cex = 0.7)
      mtext(paste0(taxa, " MDS Axis ", j), outer=TRUE, side=2, line=1)
      
      index = index + 1
    }
    
    hist(pValuesfilePrefix, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", filePrefix, " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    if (class(varCare) == "numeric"){
      hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Kendall-tau correlation test) for\n", filePrefix, " at Taxonomic Level:", taxa),
           ylab = "Count", xlab="Unadjusted p-values")
    } else {
      hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Wilcoxon test) for\n", filePrefix, " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    }
    hist(pValuesDepth, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", "Depth", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesfilePrefixwithDepth, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", filePrefix, " after depth", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    
    dev.off()
  }
}

##for (taxa in taxaLevels){
## tissues <- c("LI Lumen", "Feces", "LI Mucosa")

for(tissue in tissues){
  pdf( paste(dataType, "_", tissue, "_", filePrefix, "_ShannonDiversity.pdf", sep=""))
  par(mfrow=c(3,1),
      oma = c(1, 3, 2, 0) + 0.1,
      mar = c(3, 2, 3, 1) + 0.1,
      mgp = c(2, 0.75, 0))
  
  names <- vector()
  pValuesfilePrefix<- vector()
  pValuesfilePrefixWilcox <- vector()
  pValuesCage<- vector()
  iccCage <- vector()
  pValuesSampleType <- vector()
  pValuesAge<- vector()
  pValuesAnimal<- vector()
  pValuesInteraction <- vector()
  allpvals <- list()
  ShannonDiversityList <- list()
  ## Additional checks
  pValuesDepth <- vector()
  ## pValuesDepthWilcox <- vector()
  pValuesDepthwithfilePrefix <- vector()
  pValuesfilePrefixwithDepth <- vector()
  
  index <- 1
  for(taxa in taxaLevels){
    setwd(processedDir)
    inFileName <- paste(baseDataFileName, taxa, "_LogNormwithMetadata.txt", sep="")
    myT <- read.table(inFileName, header=TRUE, sep="\t")
    numCols <- ncol(myT)
    
    ## endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    ## The 9 could be done better...
    myColClasses <- metaColClasses<- c("character", "character", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", rep("numeric", numCols - 9))
    myT <- read.table(inFileName, header=TRUE, sep="\t", colClasses = myColClasses)
    
    endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    setwd(analysisDir)
    if (filePrefix == "Fecal.Gene.Ct") {
      myT <- myT[complete.cases(myT),]
    }
    ## names <- vector()
    ## pValuesfilePrefix<- vector()
    ## pValuesfilePrefixWilcox <- vector()
    ## pValuesCage<- vector()
    ## iccCage <- vector()
    ## pValuesSampleType <- vector()
    ## pValuesAge<- vector()
    ## pValuesAnimal<- vector()
    ## pValuesInteraction <- vector()
    ## allpvals <- list()
    ## ShannonDiversityList <- list()
    ## ## Additional checks
    ## pValuesDepth <- vector()
    ## ## pValuesDepthWilcox <- vector()
    ## pValuesDepthwithfilePrefix <- vector()
    ## pValuesfilePrefixwithDepth <- vector()
    
    
    subT <- myT[myT$Sample.Type == tissue,]
    
    subT$ShannonDiversity <- apply(subT[,(endMetadataIndex+1):ncol(subT)],1,diversity)
    bug <- subT$ShannonDiversity
    ShannonDiversityList[[index]] <- bug
    
    ## ifelse(class(subT[filePrefix][,1]) == "character",
    ##        ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character",
    ##        varCare <- factor(as.numeric(as.factor(unlist(subT[filePrefix])))),
    ##        varCare <- as.numeric(as.factor(unlist(subT[filePrefix])))
    ##        )
    
    varCare <- subT[filePrefix][[1]]
    sampleType <- factor(subT$Sample.Type)
    cage <-  factor(subT$Pen.Location)
    ## depth <- subT$depthAtLevel
    depth <- subT$sampleDepth
    myFrame <- data.frame(bug, varCare, sampleType, cage, depth)
    
    varCareModel <- lm(bug~varCare, data= myFrame)
    #print("Doesn't get here Shannon")
    depthModel <- lm(bug~depth, data = myFrame)
    depthWithCareModel <- lm(bug~varCare + depth, data=myFrame)
    mixedAnova <- anova(varCareModel)
    ## This method of getting p-values only works for the univariate case
    pValuesfilePrefix[index] <- mixedAnova$"Pr(>F)"[1]
    ## How would I change the Wilcoxon test to suit this purpose?
    ## ifelse(class(subT[filePrefix][,1]) == "character",
    ## class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character",
    ## pValuesfilePrefixWilcox[index] <- wilcox.test( bug[subT[filePrefix] == "Young"],    bug[subT[filePrefix] == "Old"])$p.value##,
    if(class(varCare) == "numeric") {
      pValuesfilePrefixWilcox[index] <- cor.test( bug, varCare, method="kendall")$p.value[1]
    } else {
      pValuesfilePrefixWilcox[index] <- pvalue(wilcox_test(bug~varCare, data=myFrame))
    }

    
    ##pValuesfilePrefixWilcox[index] <- cor.test(bug, as.numeric(as.factor(unlist(subT[filePrefix]))), method="kendall")$p.value[1]
    ##)
    pValuesDepth[index] <- anova(depthModel)$"Pr(>F)"[1]
    ## pValuesDepthwithfilePrefix[index] <-
    pValuesfilePrefixwithDepth[index] <- anova(depthWithCareModel, depthModel)$"Pr(>F)"[2]
    
    if( ##class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character"
      class(subT[filePrefix][,1]) == "character"){
      boxplot( bug ~ varCare , ylab="Shannon Diversity", pch=19,
               main = paste(filePrefix, " linear model p-value ", format(pValuesfilePrefix[index],digits=3), " Wilcoxon test p-value", format(pValuesfilePrefixWilcox[index],digits=3) ) , outline=FALSE)
      
      stripchart(bug ~ varCare ,
                 data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "Log normalized abundance", method="jitter")
    } else {
      plot( bug ~ varCare , ylab="Log normalized abundance", pch=19,
            main = paste(filePrefix," linear model p-value ", format(pValuesfilePrefix[index],digits=3), " Kendall-tau correlation test p-value ", format(pValuesfilePrefixWilcox[index],digits=3) ) )
    }
    
    plot( bug ~ depth , ylab="Shannon Diversity", pch=19,
          main = paste("Depth linear model p-value ", format(pValuesDepth[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))
    
    stripchart(bug ~ depth ,
               data = myFrame, vertical = TRUE, pch=19, add=TRUE, ylab = "Log normalized abundance", method="jitter")
    
    if (##class(as.numeric(as.factor(unlist(subT[filePrefix])))) == "character"
      class(subT[filePrefix][,1]) == "character"){
      
      plot(depth[varCare == "Young"], bug[varCare == "Young"], ylab="Log normalized abundance", xlab="Sequencing depth at this taxonomic Level", pch=19,
           main = paste(filePrefix, " and Depth linear model p-value ", format(pValuesfilePrefixwithDepth[index],digits=3)), col="blue",
           xlim = c(min(depth), max(depth)), ylim = c(min(bug), max(bug)))
      points(depth[varCare == "Old"], bug[varCare == "Old"], col="red", pch=19)
      
      deptho <- depth[order(depth)]
      varo <- varCare[order(depth)]
      bugo <- bug[order(depth)]
      lines(deptho[varo == "Young"], bugo[varo == "Young"], col = "blue")
      lines(deptho[varo == "Old"], bugo[varo == "Old"], col = "red")
    } else{
      ## TODO: Need a new plot here!
      plot.new()
    }
    
    mtext(paste0("Shannon diversity at level ",taxa), outer=TRUE, cex = 1.0)
    ## mtext("Shannon Diversity", outer=TRUE, side=2, line=1, cex = 1.1)
    
    index <- index + 1
  }
  
  ## Will have to do something about labeling the taxonomic levels
  ## hist(pValuesfilePrefix, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for ", filePrefix, "\nat Taxonomic Level:", taxa))
  ## hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Wilcoxon test) for ", filePrefix, "\nat Taxonomic Level:", taxa))
  
  dev.off()
  
  dSFrame <- data.frame(taxaLevels, pValuesfilePrefix, pValuesfilePrefixWilcox, pValuesDepth, pValuesfilePrefixwithDepth)
  
  dSFrame <- dSFrame[order(dSFrame$pValuesfilePrefix),]
  dSFrame$adjustedpValuesfilePrefix <- p.adjust( dSFrame$pValuesfilePrefix, method = "BH" )
  dSFrame$adjustedpValuesfilePrefixWilcox <- p.adjust(dSFrame$pValuesfilePrefixWilcox, method = "BH")
  rownames(dSFrame) <- dSFrame$names
  dSFrame$adjustedpValuesDepth <- p.adjust(dSFrame$pValuesDepth, method="BH")
  dSFrame$adjustedpValuesfilePrefixwithDepth <- p.adjust(dSFrame$pValuesfilePrefixwithDepth, method="BH")
  oldColnames <- colnames(dSFrame)
  ## This can probably be better written
  if(class(varCare) == "numeric"){
    colnames(dSFrame) <- c("names", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Kendall"), paste0("p-valuesDepth"), paste0("p-valuesDepthwith", filePrefix), paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Kendall"), paste0("adjusted_p-valuesDepth"), paste0("adjusted_p-valuesDepthwith", filePrefix))
  } else {
    colnames(dSFrame) <- c("names", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Wilcox"), paste0("p-valuesDepth"), paste0("p-valuesDepthwith", filePrefix), paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Wilcox"), paste0("adjusted_p-valuesDepth"), paste0("adjusted_p-valuesDepthwith", filePrefix))
  }
  
  write.table(dSFrame, file=paste(dataType, tissue, "_", filePrefix, "_ShannonDiversity.txt", sep=""), sep="\t", row.names=FALSE)
}
