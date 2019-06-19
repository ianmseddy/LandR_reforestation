
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
  reqdPkgs = list("raster", "PredictiveEcology/LandR@development"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant"),
    defineParameter("harvestInitialTime", "numeric", start(sim), NA, NA, "Time of first harvest"),
    defineParameter("reforestInitialTime", "numeric", start(sim), NA, NA, "Time of first reforest"),
    defineParameter("harvestInterval", "numeric", 1, NA, NA, "Time between harvest events"),
    defineParameter("reforestInterval", "numeric", 1, NA, NA, "Time between reforest events")
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = 'rstCurrentHarvest', objectClass = "RasterLayer", desc = "Binary raster layer with locations of harvest (represented as 1)", sourceURL = NA),
    expectsInput(objectName = 'pixelGroupMap', objectClass = "RasterLayer", desc = "Location of pixel groups"),
    expectsInput(objectName = "cohortData", objectClass = "data.table", desc = "table with attributes of cohorts that are harvested")
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = 'harvestedBiomass', objectClass = "RasterLayer", desc = "raster with harvested biomass from cohortData")
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
      sim <- scheduleEvent(sim, P(sim)$harvestInitialTime, "LandR_reforestation", "harvest")
      sim <- scheduleEvent(sim, P(sim)$reforestInitialTime, "LandR_reforestation", "reforest")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandR_reforestation", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "LandR_reforestation", "save")
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

    harvest = {
      #I want to prevent harvest from occuring if rstCurrentHarvest is unchanged
      #This will schedule it only if rstCurrentHarvest is from this year... I think
      if (!LandR::scheduleDisturbance(sim$rstCurrentHarvest, time(sim)))
      sim <- harvestBiomass(sim)

      # REVIEW HARVEST PRIORITY
      sim <- scheduleEvent(sim, time(sim) + P(sim)$harvestInterval, "LandR_reforestation", "harvest", priority = 4)

      # ! ----- STOP EDITING ----- ! #
    },
    reforest = {

      sim <- plantNewCohorts(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$reforestInterval, "LandR_reforestation", "reforest")

      # ! ----- STOP EDITING ----- ! #
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


  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}


### template for your event1
harvestBiomass <- function(sim) {
  #Calculate sum of B at harvest locations
  cohortDataLong <- LandR::addPixels2CohortData(cohortData = sim$cohortData,
                                                pixelGroupMap = sim$pixelGroupMap)



  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event2
plantNewCohorts <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {

  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("pixelGroupMap", sim)) {
    message("You should probably run LandR_BorealLBMRDataPrep")
    sim$pixelGroupMap <- raster(extent(-5,5,-5,5), res = c(1,1))
    sim$pixelGroupMap[] <- 1:100
  }

  if (!suppliedElsewhere("cohortData", sim)) {
    #Figure out what cohort data looks like
    sim$cohortData <- data.table("pixelGroup" = 1:100,
                                 "biomass" = runif(100, 1,100),
                                 "maxB" = 100,
                                 "aNPP" = 3)
  }
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
