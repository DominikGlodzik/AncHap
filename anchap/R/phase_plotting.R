

plotPhasing <- function(probandId,plotFolder, genomeMap, cumEnds, iMap, noInds, plotFileExt='') {
  # requires
  #
  # noSurrs1 - number of surrogates at a locus: genome-leng vector
  # genomeLength
  # selectedSequencesH1
  # selectedSequencesH2
  # cumEnds
  # stitches
  # betHapMismatches
  # i1
  # i2
  
# finds position of a snip on the genome
# requires a map which localises snps on the genome
scaleToGenome <- function (snpNo) {
  snp <- map[snpNo,]
  snpChrom <- as.integer(snp$chromosome);
  if (snpChrom==1) {pos <- snp$position}
  else {pos <- cumEnds[snpChrom-1] + as.integer(snp$position)}
  pos
}

# draws the chromosome borders, based on the map file
drawChromosomes2 <- function (max, cumEnds) {
  for (c in 1:22) { endP <- cumEnds[c]; lines(c(endP,endP),c(0,max),col= 'black', lwd=0.2)}  
}

drawSeqs <- function (seqs, c, l, width) {

  
  if ((dim(seqs)[1])>0) {
    for (s in 1:dim(seqs)[1]) {
        snp1 <- seqs[s,'start'];
        snp2 <- seqs[s,'end'];   
        
        id <- seqs[s,'id'];

        xs <- c(iMap[snp1], iMap[snp2-1])
        ys <- c(id,id)
 
        lines(xs, ys, col= c, lty=l, lwd=width)
    }
  }
}

drawStitches <- function (seqs, height, c) {

  noSeqs <- dim(seqs)[1]

  if (noSeqs>0) {

    for (s in 1:noSeqs) {

      snp1 <- seqs[s,'start'];
      snp2 <- seqs[s,'end'];
     
      xs <- c(iMap[snp1], iMap[snp2-1])
      ys <- c(height,height)
      lines(xs, ys, col= c, lwd=0.2)
      points(xs,  ys, lwd=0.05, pch="x", cex=0.2, col=c)
    }
  }
}

plotPoints <- function (pnts, height, c) {
    x <- iMap[pnts]
    height <- integer(length(x)) + height
    points(x,  height, lwd=0.05, pch="o", cex=0.2, col=c)
}

plotPointsData <- function (xs, ys, c) {
    for (i in 1:length(xs)) {
      xs[i] <- iMap[xs[i]]
    }  
    lines(xs,  ys, lwd=0.1, pch=".", cex=0.1, col=c)
  
}
  

  print("generating a plot ...")
  fileName <- paste(plotFolder, "phasing","_", probandId,'-',plotFileExt,'.pdf', sep="")
  
  pdf(fileName)
  print(paste("Writing the plot to", fileName))

  indx <- sample(1:length(noSurrs1), min(3000, length(iMap)), replace=FALSE)            
  indx<-indx[order(indx)]
  allSurrs <- (noSurrs1[indx] + noSurrs2[indx])
  maxSurrCount <- max(allSurrs)
  
  par(fig=c(0,1,0.22,1))
   
  plot(c(min(genomeMap$position),max(cumEnds)*1.15), c(-1 ,noInds+10),'n', ylab='Surrogate parent', xlab='', bty='o',  xaxt="n")
  mtext(paste("Orkney , proband ",  probandId, sep="", " - sequences shared with surrogate parents"), 3, cex=0.7)
  drawChromosomes2(800, cumEnds)

  drawSeqs(selectedSequencesH1, 'red', 'solid', 0.2);
  drawSeqs(selectedSequencesH2, 'green', 'solid', 0.2);


  drawSeqs(rejectedSequences, 'dodgerblue2', 'solid', 0.2);
  drawSeqs(matchingBothSequences, 'turquoise3', 'solid', 0.2);
  drawSeqs( noOverlapSequences, 'goldenrod2', 'solid', 0.2);
  drawSeqs(partialSequencesH1, 'red4', 'solid', 0.15)
  drawSeqs(partialSequencesH2, 'green4', 'solid', 0.1)

    legend("topright", c("sequences sharing \n with haplotype 1", "sequences sharing  \n with  haplotype 2", 'sequences not matching \n any haplotype', 'sequences matching \n both haplotypes', 'sequences with \n not enough overlap', 'seqs. partially \n matching haplotype 1', 'seqs. partially \n matching haplotype 2'), col = c('red','green', 'dodgerblue2', 'turquoise3', 'goldenrod2', 'red4', 'green4'),
           text.col = "black", lty = c('solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid'), lwd=c(0.2, 0.2, 0.2, 0.2, 0.2, 0.15, 0.1) , pch = c(-1, -1), cex = 0.3,  merge = TRUE, y.intersp=2)
    
    par(fig=c(0,1,0,0.45), new=TRUE) 
    plot(c(min(genomeMap$position),max(cumEnds)*1.15), c(- 18 ,maxSurrCount),'n', ylab="Count of surr seqs, inconsistencies", xlab='Position on genome', bty='o')
    drawChromosomes2(maxSurrCount, cumEnds)
  
    ##plotPointsData(indx, allSurrs-maxSurrCount, 'blue')
    plotPointsData(indx, noSurrs1[indx], 'red')
    plotPointsData(indx, noSurrs2[indx], 'green')
    lines(c(1,max(genomeMap$position)), c(0,0), col='black', lwd=0.2)

    drawStitches(stitches,   - 3, 'black')

  
  plotRest <- TRUE
    if (plotRest) {
                                        # plotting inconsistencies between the haplotypes
    plotPoints(which(betHapMismatches), -10, 'black')
    
    if (length(i1>0)) {  plotPoints(which(i1>0),  - 13, 'red')}
    if (length(i2>0)) { plotPoints(which(i2>0),  - 16, 'green')}

    
    legend("topright", c("number of sequences \n sharing with hap 1", "number of sequences \n sharing with hap 2", "stitches", "between hap \n inconsistencies", "within hap 1 \n inconsistencies", "within hap 2 \n inconsistencies"), col = c('red','green', 'black', 'black', 'red', 'green'),
           text.col = "black", lty = c(1, 1, 1, -1, -1, -1), lwd=0.2 , pch = c(-1, -1, 4, 1, 1, 1), cex = 0.3,  pt.cex=0.5, pt.lwd=1, merge = TRUE,  y.intersp=2)
  }
  
  dev.off()
  print(paste("Generated plot ", fileName))
  

}
