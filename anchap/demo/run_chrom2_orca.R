data("chr2Data")

individualsPhased <- 1:30
                                        # outputs

outputPaths <- list()
                                        # path, where the surrogacy matrix is stored
outputPaths$surrmatFile <- "../Data/results/orca749chr2/Surrs1000.Rdata";
                                        # path, to where the plots illustrating stitches and inconsistences go
outputPaths$plotFolder <-  "../Data/results/orca749chr2/plots/"
                                        # path where the phased haplotypes will be saved
outputPaths$phasingOutputFolder <- "../Data/results/orca749chr2/phasings/"
                                        # path where the results of the trio analysis will go
outputPaths$trioResultPath <- "../Data/results/orca749chr2/"
                                        # path to the plot that summarises ancestral haplotypes

#####
# run the program

result <- phase( individualsPhased = individualsPhased,
                data = chr2Data,
                outputPaths = outputPaths
                )
cat(names(result)); cat("\n")
cat(names(result$stitchingResults[[1]])); cat("\n")
cat(names(result$phasingResults[[1]])); cat("\n")
cat(names(result$trioResults)); cat("\n")
