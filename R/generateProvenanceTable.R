#' Generate a provenance table for reforestation. This details what species are planted in what ecoregionGroups
#' and optionally allows for a climate-sensitive genetic modifier
#' @param transferTable a data.table containing height reductions by seed source and ecoregion. At the moment this is
#' BECvarfut_plantation, BECvar_seed, and HTp_pred. If not provided, this will return a table where
#' every species is planted wherever possible, according to ecoregionDT
#' @param ecoregionCol the column name in ecoregion to match with transferTable - this matching is hardcoded
#' @param ecoregionKey a table for mapping ecoregion codes to the BEC zone/subzone/variant
#' column used by transferTable
#' @param sppEquiv the LandR sppEquivalency table for matchign species names
#' @param sppEquivCol the column in sppEquiv by which provenanceTable species should be named
#' @param method experimental - this will change how particular species are prioritized for planting
#' @return a provenance table object to be used with LandR_reforestation
#' @export
generateBCProvenanceTable <- function(transferTable = NULL, ecoregionCol, ecoregionKey,
                                      method = 'default', sppEquiv, sppEquivCol) {

  transferTable <- copy(transferTable)
  if (!is.null(transferTable)) {
    #due to the factorization of the zsv column in BEC zone data,
    splitUpZones <- strsplit(x = as.character(ecoregionKey[, get(ecoregionCol)]), split = "-")
    newCol <- lapply(splitUpZones, FUN = function(x){
      newX <- gsub(x, pattern = "^NA", replacement = "") %>%
        paste(., collapse = "")
    }) %>%
      unlist(.)

    ecoregionKey$newCol <- newCol
    spp <- c('BC_Forestry', eval(sppEquivCol))
    sppEquiv <- sppEquiv[, .SD, .SDcols = spp]
    transferTable <- transferTable[sppEquiv, on = c("species" = "BC_Forestry")]
    transferTable[, species := NULL]
    setnames(transferTable, eval(sppEquivCol), 'speciesCode')

    setkey(ecoregionKey, newCol)
    setkey(transferTable, BECvarfut_plantation)

    if (method == "default") {

      #subset each ecoregionGroup/species by the minimum height reduction (ie max when expressed as proportion) among provenances
      optimalProvenance <- transferTable[, score := rank(HTp_pred, ties.method = 'random'), by = .(BECvarfut_plantation, speciesCode)]# %>%
        .[, best := max(score), .(BECvarfut_plantation, speciesCode)] %>%
        .[score == best,]

      optimalProvenance <- optimalProvenance[ecoregionKey, on = c('BECvarfut_plantation' = 'newCol'), allow.cartesian = TRUE]
      setnames(optimalProvenance, old = "BECvar_seed", "Provenance")

      provenanceTable <- optimalProvenance[, .(ecoregionGroup, Provenance, speciesCode)]

    } else if (method == "Elizabeth's other ideas") {
      #implement Elizabeth's other ideas here
    } else {
      stop("unrecognized method")
    }
    #When time is available, scenarios that aren't relevant to CS - e.g. no hardwood reforestation, no migration, etc
  } else {
    stop("for now, I need the transfer tables!")
    #decide some way to make a provenance table.
  }

  provenanceTable[, ecoregionGroup := as.factor(ecoregionGroup)]
  provenanceTable[, Provenance := as.factor(Provenance)]
  #these must be factors.

  return(provenanceTable)
}

