
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "LandR_reforestation",
  description = "LandR-ecosystem module for simulating regeneration following harvest",
  keywords = c("harvest", "reforestation", "assisted migration"),
  authors = person("Ian", "Eddy", email = "ieddy@canada.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.5.9005", LandR_reforestation = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "LandR_reforestation.Rmd"),
  reqdPkgs = list("raster", "PredictiveEcology/LandR@development", "magrittr"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = "Should this entire module be run with caching activated? This is generally
                    intended for data-type modules, where stochasticity and time are not relevant"),
    defineParameter("cohortDefinitionCols", 'character', c("pixelGroup", 'age', 'speciesCode'), NA, NA,
                    desc = 'columns in cohortData that determine unique cohorts'),
    defineParameter("reforestInitialTime", "numeric", start(sim), NA, NA, "Time of first reforest. Set to NA if no
                    reforestation is desired. Harvest will still occur and map reclassified, with natural regen"),
    defineParameter("reforestInterval", "numeric", 1, NA, NA, "Time between reforest events"),
    defineParameter("simulateHarvest", 'logical', FALSE, NA, NA,
                    desc = 'generate a random 50 pixel harvest layer from pixelGroupMap for testing purposes only'),
    defineParameter("successionTimestep", "numeric", 10, NA, NA,
                    desc = "succession time step used by biomass succession module"),
    defineParameter("trackPlanting", 'logical', FALSE, NA, NA, 'add "harvest" column to cohortData that tracks planted cohorts')
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = 'rstCurrentHarvest', objectClass = "RasterLayer",
                 desc = "Binary raster layer with locations of harvest (represented as 1)", sourceURL = NA),
    expectsInput(objectName = 'pixelGroupMap', objectClass = "RasterLayer",
                 desc = "Location of pixel groups"),
    expectsInput(objectName = "cohortData", objectClass = "data.table",
                 desc = "table with attributes of cohorts that are harvested"),
    expectsInput(objectName = "speciesEcoregion", objectClass = "data.table",
                 desc = "data table with maxB and maxANPP estimates..."),
    expectsInput(objectName = "provenanceTable", objectClass = "data.table",
                 desc = "Data table with 3 columns. Location is reforestation location (as ecoregion),
                 Provenance is the provenance to be planted (as ecoregion),
                 speciesCode is the species to be planted")
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = 'harvestedBiomass', objectClass = "RasterLayer",
                  desc = "raster with harvested biomass from cohortData")
  )
))

## event types
#   - type `init` is required for initialization

doEvent.LandR_reforestation = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$reforestInitialTime, "LandR_reforestation", "reforest", eventPriority = 7)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandR_reforestation", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "LandR_reforestation", "save")

      if (P(sim)$simulateHarvest) {
        sim <- scheduleEvent(sim, P(sim)$reforestInitialTime, "LandR_reforestation", 'simulateHarvest', eventPriority = 6)
      }
    },
    plot = {
      # ! ----- EDIT BELOW ----- ! #
      Plot(sim$harvestedBiomass, new = TRUE, title = "Harvested Biomass")
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "LandR_reforestation", "plot")
    },

    save = {
      Save(sim)
      scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "LandR_reforestation", "save")
    },

    simulateHarvest = {
      sim$rstCurrentHarvest <- makeHarvestRaster(pixelGroupMap = sim$pixelGroupMap,
                                                 time = time(sim))
      scheduleEvent(sim, time(sim) + P(sim)$reforestInterval, "LandR_reforestation", 'simulateHarvest', eventPriority = 6)
    },

    reforest = {

      #planting new cohorts only scheduled if there is a non-null disturbance layer with current year
      if (!LandR::scheduleDisturbance(sim$rstCurrentHarvest, time(sim))) {
        sim <- plantNewCohorts(sim)
      }

      sim <- scheduleEvent(sim, time(sim) + P(sim)$reforestInterval, "LandR_reforestation", "reforest", eventPriority = 7)

    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {

  #initiate harvested biomass
  sim$harvestedBiomass <- raster(sim$pixelGroupMap)
  sim$harvesteBiomass[!is.na(sim$harvestedBiomass[])] <- 0

  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}


### template for your event2
plantNewCohorts <- function(sim) {

  cohortData <- copy(sim$cohortData)
  cols <- c("pixelGroup", 'speciesCode', 'ecoregionGroup', 'age', 'B') %>%
    .[. %in% colnames(cohortData)] #originally had Provenance
  cohortData <- cohortData[, ..cols]
  pixelGroupMap <- sim$pixelGroupMap

  newcols <- c(cols, "pixelIndex")
  harvestPixelCohortData <- sim$harvestedCohorts[, ..newcols]

  #this prevents harvest from occuring on pixels that were burned.
  #Need better coordination of scheduling of Biomass_regeneration, harvest, and fire
  harvestPixelCohortData <- harvestPixelCohortData[pixelGroup != 0 & !is.na(speciesCode)]

  thpt <- unique(harvestPixelCohortData[, .(pixelGroup, pixelIndex)])

  #Remove biomass from cohortData
  #treeHarvestPixelTable used to 0 pixelGroupMap

  outs <- updateCohortDataPostHarvest (newPixelCohortData = harvestPixelCohortData,
                                       cohortData = sim$cohortData,
                                       pixelGroupMap = sim$pixelGroupMap,
                                       currentTime = round(time(sim)),
                                       speciesEcoregion = sim$speciesEcoregion,
                                       cohortDefinitionCols = P(sim)$cohortDefinitionCols,
                                       treedHarvestPixelTable = thpt,
                                       provenanceTable = sim$provenanceTable,
                                       successionTimestep = P(sim)$successionTimestep,
                                       trackPlanting = P(sim)$trackPlanting)

  if (is.null(outs$cohortData$Provenance)) {
    warning("LandR_reforestation Provenance is NULL during plantNewCohorts")
  }

  sim$cohortData <- outs$cohortData
  sim$pixelGroupMap <- outs$pixelGroupMap
  sim$pixelGroupMap[] <- as.integer(sim$pixelGroupMap[])

  return(invisible(sim))
}

makeHarvestRaster <- function(pixelGroupMap, time){

  index <- 1:ncell(pixelGroupMap)
  index <- index[!is.na(pixelGroupMap[])]
  rstCurrentHarvest <- pixelGroupMap
  harvestLoc <- sample(x = index, size = 50, replace = FALSE)
  rstCurrentHarvest[!is.na(rstCurrentHarvest)] <- 0
  rstCurrentHarvest[harvestLoc] <- 1
  rstCurrentHarvest@data@attributes$Year <- time

  return(rstCurrentHarvest)
}



.inputObjects <- function(sim) {
  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("pixelGroupMap", sim)) {
    message("You should probably run LandR_BorealLBMRDataPrep")
    sim$pixelGroupMap <- raster(extent(-5,5,-5,5), res = c(1,1))
    sim$pixelGroupMap[] <- replicate(n = 10, expr = sample(1:10, 10, replace = FALSE), simplify = TRUE)
  }

  if (!suppliedElsewhere("cohortData", sim)) {
    #Figure out what cohort data looks like
    B <- runif(10, 1, 100)
    sim$cohortData <- data.table("speciesCode" = rep(c("Pice_gla", "Pice_mar", "Popu_tre", "Pinu_con", "Abie_las"), 2),
                                 "pixelGroup" = 1:10,
                                 "ecoregionGroup" = rep(c(1,2), each = 5),
                                 "age" = round(runif(10, 10,100), digits = 0),
                                 "B" = B,
                                 "totalB" = B,
                                 "mortality" = 0,
                                 "aNPPAct" = 1,
                                 'speciesProportion' = 1
                                 )
    #Add second cohort to a pixelGroup so it is a little more realistic for now
    twoCohort <- data.table("Abie_las", 1, 1, 20, sim$cohortData$B[1], sim$cohortData$B[1]*2, 0, 1, 0.5)
    names(twoCohort) <- names(sim$cohortData)
    sim$cohortData[1,]$totalB <- sim$cohortData[1,]$totalB * 2
    sim$cohortData[1,]$speciesProportion <- 0.5
    sim$cohortData <- rbind(sim$cohortData, twoCohort)
  }

  if (!suppliedElsewhere("rstCurrentHarvest", sim)) {
    message("No harvest layer supplied. Simulating one year of harvest")
    sim$rstCurrentHarvest <- sim$pixelGroupMap %>%
      setValues(., round(runif(n = ncell(.), min = 0, max = 0.55), digits = 0))
    sim$rstCurrentHarvest[is.na(sim$pixelGroupMap)] <- NA
    sim$rstCurrentHarvest@data@attributes <- list("Year" = P(sim)$reforestInitialTime)
  }

  if (!suppliedElsewhere('treedHarvestPixelTableSinceLastDisp', sim)) {

    sim$treedHarvestPixelTableSinceLastDisp <- data.table(pixelIndex = integer(0), pixelGroup = integer(0),
                                                       harvestTime = numeric(0))
  }

  if (!suppliedElsewhere('speciesEcoregion', sim)) {
    sim$speciesEcoregion <- data.table(factor(ecoregionGroup = rep(x = 1:2, each = 5)),
                                       speciesCode = rep(x = c("Pice_gla", "Pice_mar", "Popu_tre", "Pinu_con", "Abie_las"),
                                                         times = 2, ),
                                       establishprob = 0.5,
                                       maxB = runif(10, 5000, 10000),
                                       year = 0)
    sim$speciesEcoregion[, maxANPP := maxB/30]
  }

  if (!suppliedElsewhere("provenanceTable", sim)) {
    sim$provenanceTable <- data.table("ecoregionGroup" = sim$speciesEcoregion$ecoregionGroup,
                                      "Provenance" = sim$speciesEcoregion$ecoregionGroup,
                                      "speciesCode" = sim$speciesEcoregion$speciesCode)
  }
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
