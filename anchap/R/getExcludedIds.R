getExcludedIds <- function (individualPhased, excludeParents, excludeChildren, pedigree, orcaIds ) {
 # id of father of the proband, if he was genotyped

  faID <- which(pedigree[orcaIds[individualPhased],]$fa == orcaIds);
  # id of mother of the proband, if she was genotyped;
  moID <- which(pedigree[orcaIds[individualPhased],]$mo == orcaIds);

  # ids of children of the proband, where the proband is a father
  childrenMo <- which(orcaIds %in% which(pedigree$mo==orcaIds[individualPhased]))
  # ids of children of the proband, where the proband is a mother
  childrenFa <- which(orcaIds %in% which(pedigree$fa==orcaIds[individualPhased]))

  excluded<-c();
  if (excludeParents) {
    excluded <- c(excluded, faID, moID);
  }
  if (excludeChildren) {
    excluded <- c(excluded, childrenMo,childrenFa);
  }

  excluded


}
