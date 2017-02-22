data("toyData")

individualsPhased <- 1:5


#####
# run the program
result <- phase( individualsPhased = individualsPhased,
                data = toyData,
                )

avgSurrs <- 2*sum(result$surrMatDf$end - result$surrMatDf$start)/(max(individualsPhased)*length(toyData$iGene[[1]]))
cat(paste('On average there were ',avgSurrs,'surrogate parents at a locus. \n'))
