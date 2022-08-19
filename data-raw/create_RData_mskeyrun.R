# Reads data inputs from mskeyrun data package
# does *not* make interim csv files

# switch for simulated or real data
# user specifies n bins, n fleets?

library(magrittr)
library(FSA)

# Reuse Andy's functions
create_RData_mskeyrun <- function(dattype = c("sim", "real"), 
                                  nlenbin = 5) {
  
  if(dattype == "sim"){
    focalspp <- mskeyrun::simFocalSpecies
    survindex <- mskeyrun::simSurveyIndex
    survlen <- mskeyrun::simSurveyLencomp
    survdiet <- mskeyrun::simSurveyDietcomp
    survtemp <- mskeyrun::simSurveyBottemp
    fishindex <- mskeyrun::simCatchIndex
    fishlen <- mskeyrun::simFisheryLencomp
  }

  if(dattype == "real"){
    focalspp <- mskeyrun::focalSpecies %>%
      dplyr::mutate(Name = modelName)
    
    survindex <- mskeyrun::surveyIndex %>%
      dplyr::left_join(focalspp) %>%
      dplyr::mutate(Name = modelName,
                    year = YEAR,
                    survey = SEASON)
    
    survlen <- mskeyrun::realSurveyLencomp
    
    survdiet <- get_survDiet()
    
    survtemp <- get_bottemp() #use ecodata
    
    fishindex <- mskeyrun::catchIndex%>%
      dplyr::left_join(focalspp) %>%
      dplyr::mutate(Name = modelName,
                    year = YEAR,
                    survey = SEASON)
    
    fishlen <- mskeyrun::realFisheryLencomp
  }
  
  
  # read in data files associated with pin and dat file
  d <- get_DatData_msk(nlenbin,
                       focalspp,
                       survindex,
                       survlen,
                       survdiet,
                       survtemp,
                       fishindex,
                       fishlen)
  p <- get_PinData_msk(nlenbin,
                       focalspp,
                       survindex,
                       survlen,
                       survdiet,
                       survtemp,
                       fishindex,
                       fishlen)
  
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

get_DatData_msk <- function(nlenbin,
                            focalspp,
                            survindex,
                            survlen,
                            survdiet,
                            survtemp,
                            fishindex,
                            fishlen){
  # We stipulate all of the input data needed to write to the .dat file
  # Eventually it will be better to read all of these in from text files or a GUI. For now this will suffice
  d <- list() # set up list for data storage
  
  # path to data
  #path <- paste0(getwd(),"/",options$pathToDataFiles)
  
  predPrey <- survdiet %>%
    dplyr::mutate(species = Name) %>%
    dplyr::filter(prey %in% unique(species)) %>%
    dplyr::select(species, agecl, year, prey, value) %>%
    dplyr::group_by(species, prey) %>%
    dplyr::summarise(avgpreyprop = mean(value)) 
  
  foodweb <- dplyr::left_join(focalspp, predPrey, by=c("Name" = "species")) %>%
    dplyr::select(species=Name, prey, avgpreyprop) %>%
    dplyr::mutate(avgpreyprop = ceiling(avgpreyprop)) %>%
    tidyr::spread(species, avgpreyprop) %>%
    dplyr::filter(!is.na(prey)) %>%
    replace(is.na(.), 0)
  
  predOrPrey <- ifelse(colSums(foodweb[-1])>0, 1, 0)
  
  # list of species and guilds (Functional Groups)hydr
  speciesList <- focalspp[order(focalspp$Name),] %>%
    dplyr::mutate(species = Name,
                  guild = seq(1:length(Name)),
                  guildmember = seq(1:length(Name)),)
  
  # add species number to be used in writing hydra_est file
  speciesList["speciesNum"] <- as.numeric(as.factor(speciesList$species))
  d$speciesList <- as.character(speciesList$species) 
  d$guildNames <- as.character(speciesList$guild)
  d$numGuilds <- length(unique(d$guildNames))
  d$guildMembership <- speciesList$guildmember
  d$predOrPrey <- unname(predOrPrey)
  d$speciesNum <- speciesList$speciesNum
  
  # these are defaults, may be changed by user for different runs
  #singletons <- read.csv(paste0(path,"/singletons_NOBA.csv"),header=FALSE,row.names = 1)
  d$debugState <- 4
  d$Nyrs <- max(length(unique(survindex$year)), length(unique(fishindex$year)))
  d$Nspecies <- length(focalspp$Name)
  d$Nsizebins <- nlenbin
  d$Nareas <- 1
  d$Nfleets <- 2
  d$Nsurveys <- length(unique(survindex$survey))
  d$wtconv <- 1
  d$yr1Nphase <- 1
  d$recphase <- -1
  d$avgRecPhase <- 1
  d$recsigmaphase <- -1
  d$avgFPhase <- 1
  d$devRecPhase <- 4
  d$devFPhase <- 2
  d$fqphase <- 1
  d$fsphase <- 2
  d$sqphase <- 1
  d$ssphase <- 1
  d$ssigPhase <- -1
  d$csigPhase <- -1
  d$phimax <- 1
  d$scaleInitialN <- 1
  d$otherFood <- 1e+09
  d$LFISize <- 60
  d$bandwidthMetric <- 5
  d$assessmentPeriod <- 3
  d$flagMSE<- 0
  d$flagLinearRamp <- 1
  d$baselineThreshold <- 0.2
  
  # sizebin lengths--at present equal nlenbin increments across length range
  # pull observed min and max size from data to optimize start and end bin
  minmaxsurv <- survlen %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(minlen = min(lenbin, na.rm = TRUE),
                     maxlen = max(lenbin, na.rm = TRUE),
                     range = maxlen-minlen) 
  minmaxfish <- fishlen %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(minlen = min(lenbin, na.rm = TRUE),
                     maxlen = max(lenbin, na.rm = TRUE),
                     range = maxlen-minlen) 
  
  maxLrange <- dplyr::bind_rows(minmaxsurv, minmaxfish) %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(minminL = min(minlen),
                     maxmaxL = max(maxlen))
  
  equalbinwidths <- function(Lmax, Nsizebins){
    bindef <- max(1,ceiling(Lmax/Nsizebins))
    binwidths <- rep(bindef, Nsizebins)
    return(binwidths)
  }
  
  newbins <- maxLrange %>%
    dplyr::group_by(Name) %>%
    dplyr::group_modify(~broom::tidy(equalbinwidths(.$maxmaxL, d$Nsizebins))) %>%
    dplyr::mutate(lb = paste0("sizebin",seq(1:d$Nsizebins))) %>%
    tidyr::pivot_wider(names_from = lb, values_from = x)
  
  binwidth <- as.data.frame(newbins) %>%
    dplyr::select(-Name, everything()) %>%
    dplyr::rename('commentField-notused' = Name)
  
  d$binwidth <- binwidth[1:d$Nspecies,1:d$Nsizebins]
  row.names(d$binwidth) <- binwidth[,ncol(binwidth)]
  
  # length to weight coefficients/parameters
  lenwt <- read.csv(paste0(path,"/lengthweight_species_NOBA.csv"),header=TRUE)
  d$lenwta <- lenwt$a
  d$lenwtb <- lenwt$b
  
  # covariate information relating to recruitment, growth and maturity
  # recruitmentCovs <- read.csv(paste0(path,"/recruitment_covariates_NOBA.csv"),header=TRUE)
  # maturityCovs <- read.csv(paste0(path,"/maturity_covariates_NOBA.csv"),header=TRUE)
  # growthCovs <- read.csv(paste0(path,"/growth_covariates_NOBA.csv"),header=TRUE)
  
  # not used in mskeyrun, dummy variables  
  d$recruitmentCov <- rep(1, d$Nyrs)
  d$maturityCov <- rep(1, d$Nyrs)
  d$growthCov <- rep(1, d$Nyrs)
  # number of covariates
  d$NrecruitmentCov <- dim(d$recruitmentCov)[1]
  d$NmaturityCov <- dim(d$maturityCov)[1]
  d$NgrowthCov <- dim(d$growthCov)[1]
  
  # observed survey biomass
  # new long format
  obsBio <- survindex %>%
    dplyr::mutate(species = Name) %>%
    tidyr::pivot_wider(-units, names_from = "variable", values_from = "value") %>%
    dplyr::select(survey, year, species, biomass, cv)
  
  d$Nsurvey_obs <- dim(obsBio)[1]
  obsBio$survey <- as.numeric(as.factor(obsBio$survey))
  obsBio$species <- speciesList$speciesNum[match(unlist(obsBio$species), speciesList$species)]
  d$observedBiomass <- obsBio 
  
  # observed survey size composition
  modbins <- d$binwidth %>%
    tibble::rownames_to_column("species") %>%
    tidyr::pivot_longer(!species, names_to = "sizebin", values_to = "binwidth") %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(modbin.min = pcumsum(binwidth),
                  modbin.max = cumsum(binwidth)) %>%
    dplyr::mutate(sizebin = factor(sizebin, levels = unique(sizebin))) 
  
  obsSurvSize <- survlen %>%
    dplyr::mutate(species = Name) %>%
    dplyr::select(species,  year, survey, lenbin, value) %>%
    dplyr::left_join(modbins) %>%
    dplyr::filter(modbin.min <= lenbin & lenbin < modbin.max) %>% #lenbin defined as lower
    dplyr::group_by(species, year, survey, sizebin) %>%
    dplyr::summarise(sumlen = sum(value)) %>%
    tidyr::spread(sizebin, sumlen) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(inpN = rowSums(.[,-(1:3)], na.rm = TRUE)) %>%
    dplyr::mutate(type = 0) %>%
    dplyr::mutate(dplyr::across(c(dplyr::contains("sizebin")), ~./inpN)) %>%
    dplyr::select(survey, year, species, type, inpN, everything()) %>%
    dplyr::arrange(survey)
  
  # use -999 for missing value
  obsSurvSize <- obsSurvSize %>%
    replace(is.na(.),-999)
  
  # WARNING currently hardcoded cap inpN at 1000
  obsSurvSize$inpN[obsSurvSize$inpN > 1000] <- 1000
  
  d$Nsurvey_size_obs <- dim(obsSurvSize)[1]
  obsSurvSize$survey <- as.numeric(as.factor(obsSurvSize$survey))
  obsSurvSize$species <- speciesList$speciesNum[match(unlist(obsSurvSize$species), speciesList$species)]
  d$observedSurvSize <- obsSurvSize 
  
  # observed catch biomass
  # need fleet defs first. the user may need to set these up
  # WARNING currently hardcoded for 2 NOBA sim fleets
  
  fleetdef <- focalspp %>%
    dplyr::mutate(species = Name) %>%
    dplyr::arrange(species) %>%
    dplyr::mutate(allcombined = dplyr::case_when(species %in% c("Long_rough_dab", "Polar_cod") ~ 0,
                                          TRUE ~ 1)) %>%
    dplyr::mutate(allbutcod = dplyr::case_when(species=="North_atl_cod" ~ 0,
                                        TRUE ~ allcombined),
                  codfleet = dplyr::case_when(species=="North_atl_cod" ~ 2,
                                       TRUE ~ 0),
                  qind = allbutcod + codfleet) %>%
    dplyr::select(species, qind)
  
  # new long format
  obsCatch <- fishindex %>%
    dplyr::mutate(species = Name) %>%
    tidyr::pivot_wider(-units, names_from = "variable", values_from = "value") %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(totcatch = sum(catch)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(totcatch > 0) %>%
    dplyr::select(-totcatch) %>%
    dplyr::mutate(area = 1) %>%
    dplyr::left_join(fleetdef) %>%
    dplyr::mutate(fishery = qind) %>%
    dplyr::select(fishery, area, year, species, catch, cv) %>%
    dplyr::arrange(fishery, species, year)
  
  d$Ncatch_obs <- dim(obsCatch)[1]
  obsCatch$fishery <- as.numeric(as.factor(obsCatch$fishery))
  obsCatch$species <- speciesList$speciesNum[match(unlist(obsCatch$species), speciesList$species)]
  d$observedCatch <- obsCatch
  
  # observed catch size composition
  obsCatchSize <- fishlen %>%
    dplyr::mutate(species = Name) %>%
    dplyr::select(species,  year, lenbin, value) %>%
    dplyr::left_join(modbins) %>%
    dplyr::filter(modbin.min <= lenbin & lenbin < modbin.max) %>% #lenbin defined as lower
    dplyr::group_by(species, year, sizebin) %>%
    dplyr::summarise(sumlen = sum(value)) %>%
    tidyr::spread(sizebin, sumlen) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(inpN = rowSums(.[,-(1:3)], na.rm = TRUE)) %>%
    dplyr::mutate(type = 0) %>%
    dplyr::left_join(fleetdef) %>%
    dplyr::rename(fishery = qind) %>%
    dplyr::mutate(area = 1) %>%
    dplyr::mutate(dplyr::across(c(dplyr::contains("sizebin")), ~./inpN)) %>%
    dplyr::select(fishery, area, year, species, type, inpN, everything()) %>%
    dplyr::arrange(fishery, area, species, year)
  
  #add back columns with no lengths as NAs
  missing <- setdiff(modbins$sizebin, names(obsCatchSize))
  obsCatchSize[missing] <- NA_real_ 
  colorder <- c("fishery", "area", "year", "species", "type", "inpN", levels(modbins$sizebin))
  obsCatchSize <- obsCatchSize[colorder]
  
  # use -999 for missing value
  obsCatchSize <- obsCatchSize %>%
    replace(is.na(.),-999)
  
  # cap inpN at 1000
  obsCatchSize$inpN[obsCatchSize$inpN > 1000] <- 1000
  
  d$Ncatch_size_obs <- dim(obsCatchSize)[1]
  obsCatchSize$fishery <- as.numeric(as.factor(obsCatchSize$fishery))
  obsCatchSize$species <- speciesList$speciesNum[match(unlist(obsCatchSize$species), speciesList$species)]
  d$observedCatchSize <- obsCatchSize 
  
  # observed survey diet proportion by weight
  # need length at age to get diet at length
  obsSurvDiet <- survdiet %>%
    dplyr::mutate(species = Name) %>%
    
    
    
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
