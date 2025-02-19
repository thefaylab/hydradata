# Read in all of the initial csv files and create RData files to lazily load in package

create_RData <- function() {
  path <- "data-raw"
  # read in data files associated with pin and dat file
  d <- get_DatData(path)
  p <- get_PinData(path)

  d$recName <- rep("segmented",d$Nspecies) # this is a hack. NEED to sort out data inputs

  # add all of fields to hydraData
  hydraDataList <- d
  hydraDataList <- modifyList(hydraDataList,p)

  # save each one as an RData file to lazily load with package
  # for (eachName in names(d)){
  #   #assigns the list variable to its own variable
  #   assign(eachName,d[[eachName]])
  #   save(list=eachName,file=paste0("data/",eachName,".rda"))
  # }
  # for (eachName in names(p)){
  #   #assigns the list variable to its own variable
  #   assign(eachName,p[[eachName]])
  #   save(list=eachName,file=paste0("data/",eachName,".rda"))
  # }

  usethis::use_data(hydraDataList,overwrite = TRUE)

}

## subfunctions get_pinData, get_DatData

get_DatData <- function(path){
  # We stipulate all of the input data needed to write to the .dat file
  # Eventually it will be better to read all of these in from text files or a GUI. For now this will suffice
  d <- list() # set up list for data storage

  # path to data
  #path <- paste0(getwd(),"/",options$pathToDataFiles)

  # list of species and guilds (Functional Groups)hydr
  speciesList <- read.csv(paste0(path,"/speciesList_NOBA.csv"),header=TRUE)
  # add species number to be used in writing hydra_est file
  speciesList["speciesNum"] <- as.numeric(as.factor(speciesList$species))
  d$speciesList <- as.character(speciesList$species) 
  d$guildNames <- as.character(speciesList$guild)
  d$numGuilds <- length(unique(d$guildNames))
  d$guildMembership <- speciesList$guildmember
  d$predOrPrey <- speciesList$predOrPrey
  d$speciesNum <- speciesList$speciesNum

  singletons <- read.csv(paste0(path,"/singletons_NOBA.csv"),header=FALSE,row.names = 1)
  d$debugState <- singletons["debug",]
  d$Nyrs <- singletons["Nyrs",]
  d$Nspecies <- singletons["Nspecies",]
  d$Nsizebins <- singletons["Nsizebins",]
  d$Nareas <- singletons["Nareas",]
  d$Nfleets <- singletons["Nfleets",]
  d$Nsurveys <- singletons["Nsurveys",]
  d$wtconv <- singletons["wtconv",]
  d$yr1Nphase <- singletons["yr1Nphase",]
  d$recphase <- singletons["recphase",]
  d$avgRecPhase <- singletons["avg_rec_phase",]
  d$recsigmaphase <- singletons["recsigmaphase",]
  d$avgFPhase <- singletons["avg_F_phase",]
  d$devRecPhase <- singletons["dev_rec_phase",]
  d$devFPhase <- singletons["dev_F_phase",]
  d$fqphase <- singletons["fqphase",]
  d$sqphase <- singletons["sqphase",]
  d$ssigPhase <- singletons["ssig_phase",]
  d$csigPhase <- singletons["csig_phase",]
  d$phimax <- singletons["phimax",]
  d$scaleInitialN <- singletons["scaleInitialN",]
  d$otherFood <- singletons["otherFood",]
  d$LFISize <- singletons["LFI_size",]
  d$bandwidthMetric <- singletons["bandwidth_metric",]
  d$assessmentPeriod <- unlist(singletons["assessmentPeriod",])
  d$flagMSE<- unlist(singletons["flagMSE",])
  d$flagLinearRamp <- unlist(singletons["flagLinearRamp",])
  d$baselineThreshold <- unlist(singletons["baseline_threshold",])

  # sizebin lengths
  binwidth <- read.csv(paste0(path,"/length_sizebins_NOBA.csv"),header=TRUE)
  d$binwidth <- binwidth[1:d$Nspecies,1:d$Nsizebins]
  row.names(d$binwidth) <- binwidth[,ncol(binwidth)]

  # length to weight coefficients/parameters
  lenwt <- read.csv(paste0(path,"/lengthweight_species_NOBA.csv"),header=TRUE)
  d$lenwta <- lenwt$a
  d$lenwtb <- lenwt$b

  # covariate information relating to recruitment, growth and maturity
  recruitmentCovs <- read.csv(paste0(path,"/recruitment_covariates_NOBA.csv"),header=TRUE)
  maturityCovs <- read.csv(paste0(path,"/maturity_covariates_NOBA.csv"),header=TRUE)
  growthCovs <- read.csv(paste0(path,"/growth_covariates_NOBA.csv"),header=TRUE)

  d$recruitmentCov <- t(recruitmentCovs)
  d$maturityCov <- t(maturityCovs)
  d$growthCov <- t(growthCovs)
  # number of covariates
  d$NrecruitmentCov <- dim(d$recruitmentCov)[1]
  d$NmaturityCov <- dim(d$maturityCov)[1]
  d$NgrowthCov <- dim(d$growthCov)[1]

  # observed survey biomass
  #obsBio <- read.csv(paste0(path,"/observation_biomass_NOBA_BTS_fall_allbox_effic1.csv"),header=TRUE)
  #d$observedBiomass <- t(obsBio)
  # new long format
  obsBio <- read.csv(paste0(path,"/observation_biomass_NOBA_allsurvs.csv"),header=TRUE)
  d$Nsurvey_obs <- dim(obsBio)[1]
  obsBio$survey <- as.numeric(as.factor(obsBio$survey))
  obsBio$species <- speciesList$speciesNum[match(unlist(obsBio$species), speciesList$species)]
  d$observedBiomass <- obsBio 
  
  # observed survey size composition
  obsSurvSize <- read.csv(paste0(path,"/observation_lengths_NOBA_allsurvs.csv"),header=TRUE)
  d$Nsurvey_size_obs <- dim(obsSurvSize)[1]
  obsSurvSize$survey <- as.numeric(as.factor(obsSurvSize$survey))
  obsSurvSize$species <- speciesList$speciesNum[match(unlist(obsSurvSize$species), speciesList$species)]
  d$observedSurvSize <- obsSurvSize 
  
  # observed catch biomass
  #obsCatch <- read.csv(paste0(path,"/observation_catch_NOBA_census.csv"),header=TRUE)
  #d$observedCatch <- t(obsCatch)
  # new long format
  obsCatch <- read.csv(paste0(path,"/observation_catch_NOBA_allfisheries.csv"),header=TRUE)
  d$Ncatch_obs <- dim(obsCatch)[1]
  obsCatch$fishery <- as.numeric(as.factor(obsCatch$fishery))
  obsCatch$species <- speciesList$speciesNum[match(unlist(obsCatch$species), speciesList$species)]
  d$observedCatch <- obsCatch
  
  # observed catch size composition
  obsCatchSize <- read.csv(paste0(path,"/observation_lengths_NOBA_allfisheries.csv"),header=TRUE)
  d$Ncatch_size_obs <- dim(obsCatchSize)[1]
  obsCatchSize$fishery <- as.numeric(as.factor(obsCatchSize$fishery))
  obsCatchSize$species <- speciesList$speciesNum[match(unlist(obsCatchSize$species), speciesList$species)]
  d$observedCatchSize <- obsCatchSize 
  
  # observed survey diet proportion by weight
  obsSurvDiet <- read.csv(paste0(path,"/observation_diets_NOBA_allsurvs.csv"),header=TRUE)
  d$Ndietprop_obs <- dim(obsSurvDiet)[1]
  obsSurvDiet$survey <- as.numeric(as.factor(obsSurvDiet$survey))
  obsSurvDiet$species <- speciesList$speciesNum[match(unlist(obsSurvDiet$species), speciesList$species)]
  obsSurvDiet$sizebin <- as.numeric(regmatches(obsSurvDiet$sizebin, regexpr("\\d+", obsSurvDiet$sizebin)))#take the number portion of sizebinN
  d$observedSurvDiet <- obsSurvDiet
  
  # observed effort by fleet (dummy, not used in estimation model)
  obsEffort <- read.csv(paste0(path,"/observation_effort_NOBA.csv"),header=TRUE)
  d$observedEffort <- t(obsEffort)
  d$fleetNames <- (names(obsEffort)[2:(d$Nfleets+1)])

  # observed temperature
  obsTemp <- read.csv(paste0(path,"/observation_temperature_NOBA_BTS_fall_allbox_effic1.csv"),header=TRUE)
  d$observedTemperature <- t(obsTemp)

  # stomach weight
  stomachContent <- read.csv(paste0(path,"/intake_stomachContent_NOBA.csv"),header=TRUE)
  d$intakeStomach <- as.matrix(stomachContent[,2:(d$Nsizebins+1)])

  # recruitment parameters
  stockRecruit <- read.csv(paste0(path,"/recruitment_species_NOBA.csv"),header=TRUE,row.names=1)

  d$alphaEggRicker <- unlist(stockRecruit["eggRicker_alpha",])
  d$shapeEggRicker <- unlist(stockRecruit["eggRicker_shape",])
  d$betaEggRicker <- unlist(stockRecruit["eggRicker_beta",])
  d$alphaDS <- unlist(stockRecruit["DS_alpha",])
  d$shapeDS <- unlist(stockRecruit["DS_shape",])
  d$betaDS <- unlist(stockRecruit["DS_beta",])
  d$alphaGamma <- unlist(stockRecruit["gamma_alpha",])
  d$shapeGamma <- unlist(stockRecruit["gamma_shape",])
  d$betaGamma <- unlist(stockRecruit["gamma_beta",])
  d$alphaRicker <- unlist(stockRecruit["ricker_alpha",])
  d$shapeRicker <- unlist(stockRecruit["ricker_shape",])
  d$betaRicker <- unlist(stockRecruit["ricker_beta",])
  d$alphaBH <- unlist(stockRecruit["BH_alpha",])
  d$shapeBH <- unlist(stockRecruit["BH_shape",])
  d$betaBH <- unlist(stockRecruit["BH_beta",])
  d$alphaShepherd <- unlist(stockRecruit["shepherd_alpha",])
  d$shapeShepherd <- unlist(stockRecruit["shepherd_shape",])
  d$betaShepherd<- unlist(stockRecruit["shepherd_beta",])
  d$alphaHockey <- unlist(stockRecruit["hockey_alpha",])
  d$shapeHockey <- unlist(stockRecruit["hockey_shape",])
  d$betaHockey <- unlist(stockRecruit["hockey_beta",])
  d$alphaSegmented <- unlist(stockRecruit["segmented_alpha",])
  d$shapeSegmented <- unlist(stockRecruit["segmented_shape",])
  d$betaSegmented <- unlist(stockRecruit["segmented_beta",])

  d$recSigma <- unlist(stockRecruit["sigma",])
  d$recType <- unlist(stockRecruit["type",])
  d$recStochastic <- unlist(stockRecruit["stochastic",])


  # sex ratio
  sexRatio <- read.csv(paste0(path,"/sexratio_NOBA.csv"),header=TRUE,row.names=1)
  d$sexRatio <- unlist(sexRatio)

  # recruitment covariate effects. # columns = d$Nrecruitment_cov
  rec_covEffects <- read.csv(paste0(path,"/recruitment_covariateEffects_NOBA.csv"),header=TRUE,row.names=1)
  d$recruitCovEffects <- as.matrix(rec_covEffects)

  # fecundity
  fecundity_d_h <- read.csv(paste0(path,"/fecundity_species_NOBA.csv"),header=TRUE,row.names=1)
  d$fecundityd <- unlist(fecundity_d_h["d",])
  d$fecundityh <- unlist(fecundity_d_h["h",])

  fecundity_Theta <- read.csv(paste0(path,"/fecundity_Theta_NOBA.csv"),header=TRUE,row.names=1)
  d$fecundityTheta <- format(as.matrix(fecundity_Theta),digits=5)

  # maturity
  maturity <- read.csv(paste0(path,"/maturity_species_NOBA.csv"),header=TRUE,row.names=1)
  d$maturityNu <- unlist(maturity["nu",])
  d$maturityOmega <- unlist(maturity["omega",])

  maturity_covEffects <- read.csv(paste0(path,"/maturity_covariateEffects_NOBA.csv"),header=TRUE,row.names=1)
  d$maturityCovEffects <- as.matrix(maturity_covEffects)

  # growth
  growth <- read.csv(paste0(path,"/growth_species_NOBA_era1_BTS_fall_allbox_effic1.csv"),header=TRUE,row.names=1)
  d$growthPsi <- unlist(growth["psi",])
  d$growthKappa <- unlist(growth["kappa",])
  d$growthLinf <- unlist(growth["Linf",])
  d$growthK <- unlist(growth["k",])
  d$growthType <- unlist(growth["growthType",])

  growth_covEffects <- read.csv(paste0(path,"/growth_covariateEffects_NOBA.csv"),header=TRUE,row.names=1)
  d$growthCovEffects <- as.matrix(growth_covEffects)

  # intake
  intake <- read.csv(paste0(path,"/intake_species_NOBA.csv"),header=TRUE,row.names=1)
  d$intakeAlpha <- unlist(intake["alpha",])
  d$intakeBeta <- unlist(intake["beta",])

  # M1
  M1 <- read.csv(paste0(path,"/mortality_M1_NOBA.csv"),header=TRUE,row.names = 1)
  d$M1 <- as.matrix(M1)

  #foodweb
  foodweb <- read.csv(paste0(path,"/foodweb_NOBA.csv"),header=TRUE,row.names = 1)
  d$foodweb <- as.matrix(foodweb)

  #M2 size preference
  M2sizePref <- read.csv(paste0(path,"/mortality_M2sizePreference_NOBA.csv"),header=TRUE,row.names = 1)
  d$M2sizePrefMu <- as.matrix(M2sizePref["mu",])
  d$M2sizePrefSigma <- as.matrix(M2sizePref["sigma",])

  #fishery/fleet selectivity
  fisherySelectivityc<- read.csv(paste0(path,"/fishing_selectivityc_NOBA.csv"),header=TRUE,row.names = 1)
  d$fisherySelectivityc <- unlist(as.matrix(t(fisherySelectivityc)))
  fisherySelectivityd<- read.csv(paste0(path,"/fishing_selectivityd_NOBA.csv"),header=TRUE,row.names = 1)
  d$fisherySelectivityd <- unlist(as.matrix(t(fisherySelectivityd)))

  # B0 - equilibrium biomass
  B0 <- read.csv(paste0(path,"/B0_NOBA.csv"),header=TRUE,row.names = 1)
  d$B0 <- unlist(B0)

  # assessment thresholds + exploitations
  assessmentThresholds <-  read.csv(paste0(path,"/assessmentThresholds.csv"),header=TRUE)
  d$thresholds <- assessmentThresholds$thresholds
  thresholds <- d$thresholds[d$thresholds< 1]
  d$Nthresholds <- length(d$thresholds)
  d$exploitationOptions <- assessmentThresholds[,2:dim(assessmentThresholds)[2]]
  d$minMaxThresholds <- c(min(thresholds),max(thresholds))

  # additionl level added to species specific threshold
  assessmentThresholdsSpecies <-  read.csv(paste0(path,"/assessmentThresholdsSpecies_NOBA.csv"),header=TRUE,row.names = 1)
  d$thresholdSpecies <- unlist(assessmentThresholdsSpecies)

  #scaled Efort - not used
  scaledEffort <-  read.csv(paste0(path,"/observation_effortScaling_NOBA.csv"),header=TRUE)
  d$scaledEffort <- unlist(scaledEffort)

  # discard coefficient - prob of discard
  discardCoef <-  read.csv(paste0(path,"/fishing_discards_NOBA.csv"),header=TRUE,#skip=3,
                           row.names = 1)
  d$discardCoef <- (unlist(as.matrix(discardCoef)))

  # discard survival probability | discard
  discardSurvival <-  read.csv(paste0(path,"/fishing_discardsSurvival_NOBA.csv"),header=TRUE,#skip=3,
                               row.names = 1)
  d$discardSurvival <- (unlist(as.matrix(discardSurvival)))


  # fishery catchability indicator (q's)
  indicatorFisheryqs<- read.csv(paste0(path,"/fishing_q_indicator_NOBA.csv"),header=TRUE,row.names = 1)
  d$indicatorFisheryq<- unlist(as.matrix(t(indicatorFisheryqs)))

  # Autoregressive parameters for alternative error structure (AR) for survey, recruitment, catch
  ARparameters<- read.csv(paste0(path,"/observation_error.csv"),header=TRUE)
  d$ARParameters <- unlist(ARparameters)
  
  
  # residence time in the modeled area
  residenceTime <- read.csv(paste0(path,"/restime_NOBA.csv"),header=TRUE)
  d$residenceTime <- unlist(residenceTime)
  
  # mortality rate outside the modeled area
  areaMortality <- read.csv(paste0(path,"/outsidemort_NOBA.csv"),header=TRUE)
  d$areaMortality <- unlist(areaMortality)
  
  # add missing vectors from sim not used in est
  if(is.null(d$fleetMembership)) d$fleetMembership <- rep(1, d$numGuilds) #WARNING: maps guilds to fleet 1
  if(is.null(d$minExploitation)) d$minExploitation <- rep(1e-05, d$Nfleets)
  if(is.null(d$maxExploitation)) d$maxExploitation <- rep(1e-05, d$Nfleets)
    

  # replace NA
  d <- lapply(d, function(x) {x[is.na(x)] <- -999; x})
  
  
  return(d)

}



get_PinData <- function(path){
  # Stipulate all information required for the pin data file
  p <- list()
  # path to data
  # list of species and guilds (Functional Groups)
  Y1N <- read.csv(paste0(path,"/observation_Y1N_NOBA.csv"),header=TRUE,row.names=1)
  p$Y1N <- Y1N


  # Avg recruitemnt and deviations
  redundantAvgRec <- read.csv(paste0(path,"/AvgRecPinData_NOBA.csv"),header=TRUE,row.names=1)
  p$redundantAvgRec <- unlist(redundantAvgRec)
  redundantRecDevs <- read.csv(paste0(path,"/RecDevsPinData_NOBA.csv"),header=TRUE,row.names=1)
  p$redundantRecDevs <- unlist(t(as.matrix(redundantRecDevs)))

  

  # fishery catchability (q's)
  fisheryqs<- read.csv(paste0(path,"/fishing_q_NOBA.csv"),header=TRUE,row.names = 1)
  p$fisheryq<- fisheryqs


  # survey q and obs error
  survey<- read.csv(paste0(path,"/survey_info_NOBA.csv"),header=TRUE,row.names = 1)
  p$surveyq<- unlist(survey["q",])
  p$surveySigma<- unlist(survey["obs_error",])

  # fishing error
  fishery_sigma <- read.csv(paste0(path,"/fishing_error_NOBA.csv"),header=TRUE,row.names = 1)
  p$fisherySigma <- fishery_sigma
  
  # replace NA
  p <- lapply(p, function(d) {d[is.na(d)] <- -999; d})
  
  return(p)
}
