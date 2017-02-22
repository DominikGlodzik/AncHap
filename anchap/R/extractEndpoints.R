extractEndpoints <- function(map) {

  endpoints <- integer(23)
  snpChromosome <- map$chromosome

  for (c in 1:23) {
    if (c %in% snpChromosome) {
      lastChromSnp <- max(which(snpChromosome==c))
      endpoints[c] <- map[lastChromSnp,'position']
    }
  }
  endpoints
}

extractEndpointsSnp<- function(map) {

  endpoints <- integer(22)
  snpChromosome <- map$chromosome

  for (c in 1:23) {
    if (c %in% snpChromosome) {
      lastChromSnp <- max(which(snpChromosome==c))
      endpoints[c] <- lastChromSnp
    }
  }
  endpoints
}
