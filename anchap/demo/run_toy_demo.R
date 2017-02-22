data("toyData")

individualsPhased <- 1:5

                                        # outputs
outputPaths <- list()
                                        # path, where the surrogacy matrix is stored
outputPaths$surrmatFile <- "toy_dataset/results/Surrs1000.Rdata"
                                        # path, to where the plots illustrating stitches and inconsistences go
outputPaths$plotFolder <-  "toy_dataset/results/plots/"
                                        # path where the phased haplotypes will be saved
outputPaths$phasingOutputFolder <- "toy_dataset/results/phasings/"
                                        # path where the results of the trio analysis will go
outputPaths$trioResultPath <- "toy_dataset/results/"
                                        # path to the plot that summarises ancestral haplotypes



#####
# run the program
result <- phase( individualsPhased = individualsPhased,
                data = toyData,
                outputPaths = outputPaths
                )
cat(names(result)); cat("\n")
cat(names(result$stitchingResults[[1]])); cat("\n")
cat(names(result$phasingResults[[1]])); cat("\n")
cat(names(result$trioResults)); cat("\n")
