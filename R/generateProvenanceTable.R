#' Generate a provenance table for reforestation. This details what species are planted in what ecoregionGroups
#' and optionally allows for a climate-sensitive genetic modifier
#' @param transferTable a data.table containing height reductions by seed source and ecoregion. At the moment this is
#' BECvarfut_plantation, BECvar_seed, and HTp_pred. If not provided, this will return a table where
#' every species is planted wherever possible, according to ecoregionDT
#' @param ecoregionDT the table saved to the attributes of sim$ecoregionRst
#' @param ecoregionKey a table for mapping ecoregion codes to the BEC zone/subzone/variant
#' column used by transferTable
#' @param sppEquiv the LandR sppEquivalency table for matchign species names
#' @param sppEquivCol the column in sppEquiv by which provenanceTable species should be named
#' @param method experimental - this will change how particular species are prioritized for planting
#' @return a provenance table object to be used with LandR_reforestation
#' @export
generateBCProvenanceTable <- function(transferTable = NULL, ecoregionDT, ecoregionKey,
                                      method = 'default', sppEquiv, sppEquivCol) {
  browser()
  if (is.factor(ecoregionDT$ecoregion)) {
    #ecoregions below 10 have 0 in front of them. I prefer to change this than change every other use of ecoregion
    ecoregionDT[, ecoregion := as.numeric(as.character(ecoregion))]
  }

  if (!is.null(transferTable)) {
    #due to the factorization of the zsv column in BEC zone data,
    splitUpZones <- strsplit(x = as.character(ecoregionKey$zone_subzone_variant), split = "-")
    newCol <- lapply(splitUpZones, FUN = function(x){
      newX <- gsub(x, pattern = "^NA", replacement = "") %>%
        paste(., collapse = "")
    }) %>%
      unlist(.)

    ecoregionKey$zone_subzone_variant <- newCol

    setkey(ecoregionKey, zone_subzone_variant)
    setkey(transferTable, BECvarfut_plantation)

    provenanceTable <- transferTable[ecoregionKey, on = c('BECvarfut_plantation' = 'zone_subzone_variant')]
    setnames(provenanceTable, 'speciesEcoregionCode', new = 'ecoregion')
    provenanceTable <- provenanceTable[ecoregionKey, on = c("BECvar_seed" = 'zone_subzone_variant')]
    setnames(provenanceTable, 'speciesEcoregionCode', new = 'Provenance')

    #rename spp col to match sppEquivalencies
    spp <- c('BC_Forestry', eval(sppEquivCol))
    sppEquiv <- sppEquiv[, .SD, .SDcols = spp]
    provenanceTable <- provenanceTable[sppEquiv, on = c("species" = "BC_Forestry")]
    provenanceTable[, species := NULL]
    setnames(provenanceTable, eval(sppEquivCol), 'speciesCode')

    if (method == "default") {

      #subset each ecoregionGroup/species by the minimum height reduction (ie max when expressed as proportion) among provenances
      optimalProvenance <- provenanceTable[ecoregion %in% ecoregionDT$ecoregion] %>%
        .[, score := rank(HTp_pred), by = .(ecoregionGroup, speciesCode)] %>%
        .[score == 1, .(speciesCode, ecoregionGroup, Provenance)]

      #These are ecoregions - need to add back in 'speciesEcoregion' aka the ecoregion/landcover grouping
      provenanceTable <- optimalProvenance[ecoregionDT, on = c("ecoregionGroup" = "ecoregion"), allow.cartesian = TRUE]
      #ecoregionGroup represents the ecoregionGroup here - because there is no 'ecoregionGroup' in the AMAT tables,
      provenanceTable[, ecoregionGroup := ecoregionGroup]
      provenanceTable[, c('ecoregionGroup', 'ID') := c(NULL, NULL)]

      #You don't need to know about land cover with provenance - growth modifications will only be done by ecoregion
      #We could optionally pass the ecoregionMap table to harvestPixelCohortData - join ecoregions to ecoregionGroup,
      #and join that resulting table to provenance Table via one-to-many
      #It is irrelevant whether species-ecoregionGroup combinations exist that aren't in speciesEcoregion - they will never be planted
      #ID is the value inside sim$ecoregionMap - It isn't necessary might be better to leave it in for analyses?

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

