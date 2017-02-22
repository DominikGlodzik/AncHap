# inputs
#### uncomment to read in the full Orkney dataset

genotypesFile <- "../Data/data/Orca/orca749chr2/orca749chr2";
mapFile <- "../Data/data/Orca/orca749chr2/orca749chr2.map";
pedigreeFilepath <- "../Data/data/Orca/orca749chr2/Yurii.xls"
iGeneFilepath <- "../Data/data/Orca/orca749chr2/iGene.Rdata"

chr2Data <- loadData(mapFile = mapFile,
                     genotypesFile = genotypesFile,
                     pedigreeFilepath = pedigreeFilepath,
                     iGeneFile = iGeneFilepath,
                     iGeneCreate = TRUE)
