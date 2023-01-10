# Reads data inputs from mskeyrun data package
# does *not* make interim csv files

# switch for simulated or real data
# user specifies n bins, n fleets?

library(magrittr)
library(FSA)

# Reuse Andy's functions
create_RData_mskeyrun <- function(dattype = c("sim", "real"), 
                                  nlenbin = 5) {
  
  path <- "data-raw" #will still read in some csvs
  
  if(dattype == "sim"){
    focalspp <- mskeyrun::simFocalSpecies
    survindex <- mskeyrun::simSurveyIndex
    survlen <- mskeyrun::simSurveyLencomp
    survagelen <- mskeyrun::simSurveyAgeLencomp
    survlensamp <- NULL
    focalprey <- NULL
    survdiet <- mskeyrun::simSurveyDietcomp
    survdietlen <- NULL
    survtemp <- mskeyrun::simSurveyBottemp
    survbiopar <- mskeyrun::simBiolPar
    fishindex <- mskeyrun::simCatchIndex
    fishlen <- mskeyrun::simFisheryLencomp
    fishlensamp <- NULL
    percapcons <- mskeyrun::simPerCapCons
    startpars <- mskeyrun::simStartPars
    
    fleetdef <- focalspp %>%
      dplyr::mutate(species = Name) %>%
      dplyr::arrange(species) %>%
      dplyr::mutate(allcombined = dplyr::case_when(species %in% c("Long_rough_dab", "Polar_cod") ~ 0,
                                                   TRUE ~ 1)) %>%
      dplyr::mutate(allbutcod = dplyr::case_when(species=="North_atl_cod" ~ 0,
                                                 TRUE ~ allcombined),
                    codfleet = dplyr::case_when(species=="North_atl_cod" ~ 2,
                                                TRUE ~ 0),
                    qind = allbutcod + codfleet,
                    codfleet = codfleet/2) %>%
      dplyr::select(species, allbutcod, codfleet, qind)
    
  }

  if(dattype == "real"){
    modyears <- 1978:2019  # 1968 agreed for project, testing with shorter dataset
  
    focalspp <- mskeyrun::focalSpecies %>%
      dplyr::filter(modelName != "Pollock") %>% # not using in these models
      dplyr::mutate(Name = modelName) 
    
    # single index, use if fitting to two not working
    survindex.all <- mskeyrun::surveyIndex %>%
      dplyr::filter(YEAR %in% modyears) %>%
      dplyr::left_join(focalspp %>% dplyr::select(-NESPP3) %>% dplyr::distinct()) %>%
      dplyr::filter(!is.na(modelName)) %>%
      dplyr::mutate(ModSim = "Actual",
                    year = YEAR,
                    Code = SPECIES_ITIS,
                    Name = modelName,
                    survey = paste0("Combined ", SEASON)) %>%
      dplyr::select(ModSim, year, Code, Name, survey, variable, value, units)

    survindex.AL <- mskeyrun::surveyIndexA4 %>%
      dplyr::filter(YEAR %in% modyears) %>%
      dplyr::left_join(focalspp %>% dplyr::select(-NESPP3) %>% dplyr::distinct()) %>%
      dplyr::filter(!is.na(modelName)) %>%
      dplyr::mutate(ModSim = "Actual",
                    year = YEAR,
                    Code = SPECIES_ITIS,
                    Name = modelName,
                    survey = paste0("Albatross ", SEASON)) %>%
      dplyr::select(ModSim, year, Code, Name, survey, variable, value, units)

    survindex.HB <- mskeyrun::surveyIndexHB %>%
      dplyr::filter(YEAR %in% modyears) %>%
      dplyr::left_join(focalspp %>% dplyr::select(-NESPP3) %>% dplyr::distinct()) %>%
      dplyr::filter(!is.na(modelName)) %>%
      dplyr::mutate(ModSim = "Actual",
                    year = YEAR,
                    Code = SPECIES_ITIS,
                    Name = modelName,
                    survey = paste0("Bigelow ", SEASON)) %>%
      dplyr::select(ModSim, year, Code, Name, survey, variable, value, units)
    
    # try to fit separately
    #survindex <- dplyr::bind_rows(survindex.AL, survindex.HB) 
    
    # try to simplify
    survindex <- survindex.all
    
    # rename for functions below
    survindex$variable <- stringr::str_replace_all(survindex$variable, "strat.biomass", "biomass")
    
    # observed survey biomass
    towarea <- 0.0384
    
    GBsurvstrata  <- c(1090, 1130:1210, 1230, 1250, 3460, 3480, 3490, 3520:3550)
    
    GBarea <- FishStatsUtils::northwest_atlantic_grid %>%
      dplyr::filter(stratum_number %in% GBsurvstrata) %>%
      dplyr::summarise(area = sum(Area_in_survey_km2)) %>%
      as.double()
    
    # get cv = stdev/bio = sqrt(variance)/bio
    survindex <- survindex %>%
      dplyr::filter(variable %in% c("biomass", "biomass.var")) %>%
      dplyr::select(-units) %>%
      tidyr::pivot_wider(names_from = "variable", values_from = "value") %>%
      dplyr::mutate(cv = sqrt(biomass.var)/biomass) %>%
      dplyr::select(-biomass.var) %>%
      # expand: area of tow =  0.0384 and GB area = above and kg to tons
      dplyr::mutate(biomass = biomass * (GBarea/towarea) /1000) %>% 
      tidyr::pivot_longer(c(biomass, cv), names_to = "variable", values_to = "value") %>%
      dplyr::mutate(units = ifelse(variable=="biomass", "tons", "unitless"))
    
    # need to add sample size for lengths from mskeyrun::surveyLenSampN
    survlen <- mskeyrun::realSurveyLennumcomp %>%
      dplyr::filter(year %in% modyears,
                    variable == "numbers") %>%
      dplyr::mutate(Code = as.character(Code)) %>%
      dplyr::left_join(focalspp %>% dplyr::select(-NESPP3) %>% dplyr::distinct(), by=c("Name" = "LongName", "Code" = "SPECIES_ITIS")) %>%
      dplyr::mutate(Name = modelName) %>%
      dplyr::mutate(vessel = dplyr::case_when(year<2009 ~ "Albatross",
                                              year>=2009 ~ "Bigelow", 
                                              TRUE ~ as.character(NA)),
                    #survey = paste(vessel, season)) %>%
                    survey = paste0("Combined ", season)) %>%
      dplyr::select(ModSim, year, Code, Name, survey, lenbin, variable, value, units)
    
    survagelen <- NULL #diet already has length in it
    
    survlensamp <- mskeyrun::surveyLenSampN %>%
      dplyr::ungroup() %>%
      dplyr::mutate(year = as.integer(YEAR),
                    Code = as.character(SPECIES_ITIS),
                    vessel = dplyr::case_when(year<2009 ~ "Albatross",
                                              year>=2009 ~ "Bigelow", 
                                              TRUE ~ as.character(NA)),
                    #survey = paste(vessel, SEASON)) %>%
                    survey = paste0("Combined ", SEASON)) %>%
      dplyr::select(year, Code, survey, lensampsize, ntows)
    
    focalprey <- focalspp %>%
      dplyr::select(-NESPP3) %>%
      dplyr::distinct() %>%
      dplyr::select(SCIENTIFIC_NAME, SPECIES_ITIS, Name)
    
    survdiet <- mskeyrun::surveyDietcomp %>%
      dplyr::filter(year %in% modyears,
                    season %in% c("FALL", "SPRING")) %>%
      dplyr::mutate(vessel = dplyr::case_when(year<2009 ~ "Albatross",
                                              year>=2009 ~ "Bigelow", 
                                              TRUE ~ as.character(NA)),
                    #survey = paste(vessel, season),
                    survey = paste0("Combined ", season),
                    agecl = NA) %>%
      dplyr::left_join(focalprey, by=c("prey" = "SCIENTIFIC_NAME")) %>%
      dplyr::rename(preySci = prey,
                    prey = Name.y,
                    Name = Name.x) %>%
      dplyr::filter(variable == "relmsw") %>%
      dplyr::mutate(value = value/100)
    
    survdietlen <- mskeyrun::surveyLenDietcomp %>%
      dplyr::filter(year %in% modyears,
                    season %in% c("FALL", "SPRING")) %>%
      dplyr::mutate(vessel = dplyr::case_when(year<2009 ~ "Albatross",
                                              year>=2009 ~ "Bigelow", 
                                              TRUE ~ as.character(NA)),
                    #survey = paste(vessel, season))
                    survey = paste0("Combined ", season))
    
    survtempdat <- ecodata::bottom_temp %>% 
      dplyr::filter(EPU == "GB",
                    Time %in% modyears,
                    Var %in% c("bottom temp anomaly in situ",
                               "reference bt in situ (1981-2010)")) %>%  
      tidyr::pivot_wider(names_from = "Var", values_from = "Value") %>%  
      dplyr::mutate(ModSim = "Actual",
                    year = Time,
                    survey = paste0("2022ecodata::bottom_temp ",Time),
                    variable = "bottomtemp",
                    value = `bottom temp anomaly in situ` +`reference bt in situ (1981-2010)`,
                    units = Units) %>%
      dplyr::select(ModSim, year, survey, variable, value, units)
    
    if(modyears[1]<1978){
    survtempfill <- data.frame(ModSim = "Actual",
                           year = setdiff(survindex$year, survtempdat$year),
                           survey = "2022ecodata::bottom_temp mean 1977-1981",
                           variable = "bottomtemp",
                           value = mean(survtempdat$value[survtempdat$year %in% 1977:1981]),
                           units = "degreesC")
    
    survtemp <- dplyr::bind_rows(survtempfill, survtempdat)
    }else{survtemp <- survtempdat}
    
    survbiopar <- mskeyrun::realBiolPar

    #fix for real species simplest allocation: mackerel herring fleet 2 all else fleet 1
    fleetdef <- focalspp %>%
      dplyr::select(-NESPP3) %>% 
      dplyr::distinct() %>%
      dplyr::mutate(species = Name) %>%
      dplyr::arrange(species) %>%
      # #2 fleet config used for WGSAM model fitting fall 2022
      # dplyr::mutate(pelagic = dplyr::case_when(species %in% c("Atlantic_herring",
      #                                                         "Atlantic_mackerel") ~ 2,
      #                                          TRUE ~ 0),
      #               demersal = dplyr::case_when(!species %in% c("Atlantic_herring",
      #                                                           "Atlantic_mackerel") ~ 1,
      #                                           TRUE ~ 0),
      #               qind = pelagic + demersal,
      #               pelagic = pelagic/2) %>%
      # dplyr::select(species, demersal, pelagic, qind)
      #
      # #species = fleet config suggested at WGSAM Dec 2022
      dplyr::mutate(fleetnames = Name,
                    value = 1,
                    qind = which(Name==Name)) %>%
      dplyr::select(species, fleetnames, value, qind) %>%
      tidyr::pivot_wider(names_from = fleetnames, values_from = value, values_fill = 0)
      
    
    # need to have landings + discards = catch
    fishindex <- mskeyrun::catchIndex %>%
      dplyr::filter(YEAR %in% modyears) %>%
      dplyr::left_join(focalspp %>% dplyr::mutate(NESPP3 = as.integer(NESPP3))) %>%
      dplyr::filter(!is.na(Name)) %>%
      dplyr::select(-units) %>%
      tidyr::pivot_wider(names_from = "variable", values_from = "value") %>%
      dplyr::group_by(YEAR, Name) %>%
      dplyr::mutate(catch = sum(`commercial landings`,`commercial discards`, na.rm = TRUE),
                    cv = 0.05) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_longer(c(`commercial landings`, `commercial discards`, catch, cv), names_to = "variable", values_to = "value") %>%
      dplyr::mutate(units = ifelse(variable=="cv", "unitless", "metric tons")) %>%
      dplyr::left_join(fleetdef, by=c("Name" = "species")) %>%
      dplyr::mutate(ModSim = "Actual",
                    year = YEAR,
                    Code = SPECIES_ITIS,
                    Name = modelName,
                    fishery = qind) %>%
      dplyr::select(ModSim, year, Code, Name, fishery, variable, value, units)
    
    #fishlen <- mskeyrun::realFisheryLencomp %>% #identical column names to sim, single fishery in real data
    fishlen <- mskeyrun::realFisheryLencompRaw %>% #this dataset does not borrow lengths  
      dplyr::filter(year %in% modyears,
                    variable == "abundance") %>%
      dplyr::mutate(Code = as.character(Code)) %>%
      dplyr::left_join(focalspp %>% dplyr::select(-NESPP3) %>% dplyr::distinct(), by=c("Name" = "LongName", "Code" = "SPECIES_ITIS")) %>%
      dplyr::mutate(Name = modelName) %>%
      dplyr::left_join(fleetdef, by=c("Name" = "species")) %>%
      dplyr::mutate(fishery = qind) %>%
      dplyr::select(ModSim, year, Code, Name, fishery, lenbin, variable, value, units) #%>%
      #dplyr::filter(lenbin < 250) #temporary fix to remove goosefish lenbin 347
    
    fishlensamp <- mskeyrun::fisheryLenSampN %>%
      dplyr::mutate(year = as.integer(YEAR),
                    Code = as.character(SPECIES_ITIS)) %>%
      dplyr::select(year, Code, lensampsize, ntrips)
    
    percapcons <- NULL # read existing csv of GB stomwt 
    
    startpars <- NULL # read existing csvs of Y1N and avgrec
  }
  
  # calculate model length bins for both pin and dat
  bindef <- get_modbins(nlenbin,
                        focalspp,
                        survlen,
                        fishlen)
  
  # read in data files associated with pin and dat file
  d <- get_DatData_msk(dattype,
                       nlenbin,
                       focalspp,
                       bindef,
                       survindex,
                       survlen,
                       survlensamp,
                       survagelen,
                       focalprey,
                       survdiet,
                       survdietlen,
                       survtemp,
                       survbiopar,
                       fleetdef,
                       fishindex,
                       fishlen,
                       fishlensamp,
                       percapcons,
                       path)
  p <- get_PinData_msk(dattype,
                       nlenbin,
                       speciesList = d$speciesList,
                       bindef,
                       startpars,
                       Nyrs = d$Nyrs,
                       Nfleets = d$Nfleets,
                       Nsurveys = d$Nsurveys,
                       Nareas = d$Nareas,
                       Npreds = sum(d$predOrPrey),
                       fqind = d$indicatorFisheryq,
                       observedCatch = d$observedCatch,
                       observedBiomass = d$observedBiomass,
                       path)
  
  d$recName <- rep("avgplusdevs",d$Nspecies) # this is a hack. NEED to sort out data inputs
  
  # add all of fields to hydraData
  hydraDataList_msk <- d
  hydraDataList_msk <- modifyList(hydraDataList_msk,p)
  
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
  
  usethis::use_data(hydraDataList_msk,overwrite = TRUE)
  
}

## subfunctions get_modbins,get_pinData, get_DatData

get_modbins <- function(nlenbin,
                        focalspp,
                        survlen,
                        fishlen){
  # sizebin lengths--at present equal nlenbin increments across length range
  # pull observed min and max size from data to optimize start and end bin
  bindef <- list()
  
  Nsizebins <- nlenbin
  Nspecies <- length(unique(focalspp$Name))
  
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
    dplyr::group_modify(~broom::tidy(equalbinwidths(.$maxmaxL, Nsizebins))) %>%
    dplyr::mutate(lb = paste0("sizebin",seq(1:Nsizebins))) %>%
    tidyr::pivot_wider(names_from = lb, values_from = x)
  
  binwidth <- as.data.frame(newbins) %>%
    dplyr::select(-Name, everything()) %>%
    dplyr::rename('commentField-notused' = Name)
  
  bindef$binwidth <- binwidth[1:Nspecies,1:Nsizebins]
  row.names(bindef$binwidth) <- binwidth[,ncol(binwidth)]
  
  modbins <- bindef$binwidth %>%
    tibble::rownames_to_column("species") %>%
    tidyr::pivot_longer(!species, names_to = "sizebin", values_to = "binwidth") %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(modbin.min = pcumsum(binwidth),
                  modbin.max = cumsum(binwidth)) %>%
    dplyr::mutate(sizebin = factor(sizebin, levels = unique(sizebin))) 
  
  bindef$modbins <- modbins
  
  return(bindef)
  
}

get_DatData_msk <- function(dattype,
                            nlenbin,
                            focalspp,
                            bindef,
                            survindex,
                            survlen,
                            survlensamp,
                            survagelen,
                            focalprey,
                            survdiet,
                            survdietlen,
                            survtemp,
                            survbiopar,
                            fleetdef,
                            fishindex,
                            fishlen,
                            fishlensamp,
                            percapcons,
                            path){
  # We stipulate all of the input data needed to write to the .dat file
  # Eventually it will be better to read all of these in from text files or a GUI. For now this will suffice
  d <- list() # set up list for data storage
  
  # path to data
  #path <- paste0(getwd(),"/",options$pathToDataFiles)
  
  dietyrs <- survdiet %>% dplyr::summarise(count = dplyr::n_distinct(year))
  
  predPrey <- survdiet %>%
    dplyr::mutate(species = Name) %>%
    dplyr::filter(prey %in% unique(species)) %>%
    dplyr::select(species, agecl, year, prey, value) %>%
    dplyr::group_by(species, prey) %>%
    dplyr::summarise(yrct = dplyr::n_distinct(year),
                     avgpreyprop = mean(value, na.rm=T)) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(maxyrct =  max(yrct)) %>%
    dplyr::filter(maxyrct > 0.3 * dietyrs$count) # drops preds with isolated obs 
                           #0.15 keeps haddock 0.2 keeps mackerel
  
  neverprey <- setdiff(predPrey$species, predPrey$prey)
  
  foodweb <- dplyr::left_join(focalspp %>% dplyr::select(Name) %>% dplyr::distinct(), 
                              predPrey, by=c("Name" = "species")) %>%
    dplyr::select(species=Name, prey, avgpreyprop) %>%
    dplyr::mutate(avgpreyprop = ceiling(avgpreyprop)) %>%
    tidyr::spread(species, avgpreyprop) %>%
    dplyr::mutate(prey = dplyr::case_when(is.na(prey) ~ neverprey,
                                   TRUE ~ prey)) %>%
    dplyr::arrange(prey) %>%
    replace(is.na(.), 0)
  
  predOrPrey <- ifelse(colSums(foodweb[-1])>0, 1, 0)
  
  # list of species and guilds (Functional Groups)hydr
  #speciesList <- focalspp[order(focalspp$Name),] %>%
  speciesList <- focalspp %>% 
    dplyr::select(Name) %>% 
    dplyr::distinct() %>% 
    dplyr::arrange(Name) %>%
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
  d$debugState <- 0
  d$Nyrs <- max(length(unique(survindex$year)), length(unique(fishindex$year)))
  d$Nspecies <- length(unique(focalspp$Name))
  d$Nsizebins <- nlenbin
  d$Nareas <- 1
  d$Nfleets <- length(unique(fishindex$fishery)) 
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
  d$m1phase <-  -1
  d$oF1phase <- -1
  d$oFdevphase <- -1
  d$phimax <- 1
  d$scaleInitialN <- 1
  d$otherFood <- 1e+09
  d$LFISize <- 60
  d$bandwidthMetric <- 5
  d$assessmentPeriod <- 3
  d$flagMSE<- 0
  d$flagLinearRamp <- 1
  d$baselineThreshold <- 0.2
  
  # # sizebin lengths--at present equal nlenbin increments across length range
  # # pull observed min and max size from data to optimize start and end bin
  # minmaxsurv <- survlen %>%
  #   dplyr::group_by(Name) %>%
  #   dplyr::summarise(minlen = min(lenbin, na.rm = TRUE),
  #                    maxlen = max(lenbin, na.rm = TRUE),
  #                    range = maxlen-minlen) 
  # minmaxfish <- fishlen %>%
  #   dplyr::group_by(Name) %>%
  #   dplyr::summarise(minlen = min(lenbin, na.rm = TRUE),
  #                    maxlen = max(lenbin, na.rm = TRUE),
  #                    range = maxlen-minlen) 
  # 
  # maxLrange <- dplyr::bind_rows(minmaxsurv, minmaxfish) %>%
  #   dplyr::group_by(Name) %>%
  #   dplyr::summarise(minminL = min(minlen),
  #                    maxmaxL = max(maxlen))
  # 
  # equalbinwidths <- function(Lmax, Nsizebins){
  #   bindef <- max(1,ceiling(Lmax/Nsizebins))
  #   binwidths <- rep(bindef, Nsizebins)
  #   return(binwidths)
  # }
  # 
  # newbins <- maxLrange %>%
  #   dplyr::group_by(Name) %>%
  #   dplyr::group_modify(~broom::tidy(equalbinwidths(.$maxmaxL, d$Nsizebins))) %>%
  #   dplyr::mutate(lb = paste0("sizebin",seq(1:d$Nsizebins))) %>%
  #   tidyr::pivot_wider(names_from = lb, values_from = x)
  # 
  # binwidth <- as.data.frame(newbins) %>%
  #   dplyr::select(-Name, everything()) %>%
  #   dplyr::rename('commentField-notused' = Name)
  # 
  # d$binwidth <- binwidth[1:d$Nspecies,1:d$Nsizebins]
  # row.names(d$binwidth) <- binwidth[,ncol(binwidth)]
  
  d$binwidth <- bindef$binwidth
  
  # length to weight coefficients/parameters
  lenwt <- survbiopar %>%
    dplyr::mutate(species = Name) %>%
    dplyr::arrange(species)
    
  d$lenwta <- lenwt$WLa
  d$lenwtb <- lenwt$WLb
  
  # covariate information relating to recruitment, growth and maturity
  # recruitmentCovs <- read.csv(paste0(path,"/recruitment_covariates_NOBA.csv"),header=TRUE)
  # maturityCovs <- read.csv(paste0(path,"/maturity_covariates_NOBA.csv"),header=TRUE)
  # growthCovs <- read.csv(paste0(path,"/growth_covariates_NOBA.csv"),header=TRUE)
  
  # not used in mskeyrun, dummy variables  
  d$recruitmentCov <- matrix(rep(1, d$Nyrs), nrow = 1, byrow = TRUE)
  d$maturityCov <- matrix(rep(1, d$Nyrs), nrow = 1, byrow = TRUE)
  d$growthCov <- matrix(rep(1, d$Nyrs), nrow = 1, byrow = TRUE)
  # number of covariates
  d$NrecruitmentCov <- dim(d$recruitmentCov)[1]
  d$NmaturityCov <- dim(d$maturityCov)[1]
  d$NgrowthCov <- dim(d$growthCov)[1]
  
  fitstartyr <- survindex$year[1]-1
  
  # new long format
  obsBio <- survindex %>%
    dplyr::mutate(species = Name) %>%
    dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
    tidyr::pivot_wider(-units, names_from = "variable", values_from = "value") %>%
    dplyr::select(survey, year, species, biomass, cv) %>%
    dplyr::filter(!is.na(cv)) #no 0 catch or missing cv
  
  d$Nsurvey_obs <- dim(obsBio)[1]
  obsBio$survey <- as.numeric(as.factor(obsBio$survey))
  obsBio$species <- speciesList$speciesNum[match(unlist(obsBio$species), speciesList$species)]
  d$observedBiomass <- obsBio 
  
  # observed survey size composition
  # modbins <- d$binwidth %>%
  #   tibble::rownames_to_column("species") %>%
  #   tidyr::pivot_longer(!species, names_to = "sizebin", values_to = "binwidth") %>%
  #   dplyr::group_by(species) %>%
  #   dplyr::mutate(modbin.min = pcumsum(binwidth),
  #                 modbin.max = cumsum(binwidth)) %>%
  #   dplyr::mutate(sizebin = factor(sizebin, levels = unique(sizebin))) 
  
  modbins <- bindef$modbins
  
  if(dattype == "sim") {
    
    obsSurvSize <- survlen %>%
      dplyr::mutate(species = Name) %>%
      dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
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
    
    # use 0 for missing value
    obsSurvSize <- obsSurvSize %>%
      replace(is.na(.),0)
    
    # WARNING currently hardcoded cap inpN at 1000
    obsSurvSize$inpN[obsSurvSize$inpN > 1000] <- 1000
    
  }
  
  if(dattype == "real") {
    
    obsSurvSize <- survlen %>%
      dplyr::mutate(species = Name) %>%
      dplyr::select(species,  Code, year, survey, lenbin, value) %>%
      dplyr::left_join(modbins) %>%
      dplyr::filter(modbin.min <= lenbin & lenbin < modbin.max) %>% #lenbin defined as lower
      dplyr::group_by(species, Code, year, survey, sizebin) %>%
      dplyr::summarise(sumlen = sum(value)) %>%
      tidyr::spread(sizebin, sumlen) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(totN = rowSums(.[,-(1:4)], na.rm = TRUE)) %>%
      dplyr::left_join(survlensamp) %>%
      dplyr::rename(inpN = ntows) %>%
      dplyr::mutate(type = 0) %>%
      dplyr::mutate(dplyr::across(c(dplyr::contains("sizebin")), ~./totN)) %>%
      dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
      dplyr::select(-Code, -totN, -lensampsize) %>%
      dplyr::select(survey, year, species, type, inpN, everything()) %>%
      dplyr::arrange(survey)
    
    # use 0 for missing value
    obsSurvSize <- obsSurvSize %>%
      replace(is.na(.),0)
    
  }
  
  d$Nsurvey_size_obs <- dim(obsSurvSize)[1]
  obsSurvSize$survey <- as.numeric(as.factor(obsSurvSize$survey))
  obsSurvSize$species <- speciesList$speciesNum[match(unlist(obsSurvSize$species), speciesList$species)]
  d$observedSurvSize <- obsSurvSize 
  
  # observed catch biomass
  # need fleet defs first. the user will need to set these up above
  # WARNING currently hardcoded for 2 NOBA sim fleets
  
  
  # new long format
  obsCatch <- fishindex %>%
    dplyr::mutate(species = Name) %>%
    dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
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
  
  
  if(dattype == "sim"){
    
  # observed catch size composition
  obsCatchSize <- fishlen %>%
    dplyr::mutate(species = Name) %>%
    dplyr::select(species, year, lenbin, value) %>%
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
    dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
    dplyr::select(fishery, area, year, species, type, inpN, everything()) %>%
    dplyr::arrange(fishery, area, species, year)
  
  # cap inpN at 1000
  obsCatchSize$inpN[obsCatchSize$inpN > 1000] <- 1000
  }
  
  
  if(dattype == "real"){
    
    # observed catch size composition
    obsCatchSize <- fishlen %>%
      dplyr::mutate(species = Name) %>%
      dplyr::select(species, Code, year, lenbin, value) %>%
      dplyr::left_join(modbins) %>%
      dplyr::filter(modbin.min <= lenbin & lenbin < modbin.max) %>% #lenbin defined as lower
      dplyr::group_by(species, Code, year, sizebin) %>%
      dplyr::summarise(sumlen = sum(value)) %>%
      tidyr::spread(sizebin, sumlen) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(totN = rowSums(.[,-(1:3)], na.rm = TRUE)) %>%
      dplyr::left_join(fishlensamp) %>%
      dplyr::rename(inpN = ntrips) %>%
      dplyr::mutate(type = 0) %>%
      dplyr::left_join(fleetdef %>% dplyr::select(species, qind)) %>%
      dplyr::rename(fishery = qind) %>%
      dplyr::mutate(area = 1) %>%
      dplyr::mutate(dplyr::across(c(dplyr::contains("sizebin")), ~./totN)) %>%
      dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
      dplyr::select(-totN, -Code) %>%
      dplyr::select(fishery, area, year, species, type, inpN, everything()) %>%
      dplyr::arrange(fishery, area, species, year)
    
    # some catch at length has no matching sample size (multiple species)
    # temporary fix
    obsCatchSize$inpN[is.na(obsCatchSize$inpN)] <- 5
    
    # some sample size has no matching catch at length (winter skate)
    # temporary fix
    obsCatchSize <- obsCatchSize %>%
      dplyr::filter(!dplyr::if_all(c(dplyr::contains("sizebin")), is.na))
    
    # there are sizebins with 0 proportions
    # temporary fix, change these to NAs
    obsCatchSize <- obsCatchSize %>%
      dplyr::mutate(dplyr::across(c(dplyr::contains("sizebin")), ~replace(., . ==  0 , NA)))
    
  }
  
  
  #add back columns with no lengths as NAs
  missing <- setdiff(modbins$sizebin, names(obsCatchSize))
  obsCatchSize[missing] <- NA_real_ 
  colorder <- c("fishery", "area", "year", "species", "type", "inpN", levels(modbins$sizebin))
  obsCatchSize <- obsCatchSize[colorder]
  
  # use 0 for missing value
  obsCatchSize <- obsCatchSize %>%
    replace(is.na(.),0)
  
  
  d$Ncatch_size_obs <- dim(obsCatchSize)[1]
  obsCatchSize$fishery <- as.numeric(as.factor(obsCatchSize$fishery))
  obsCatchSize$species <- speciesList$speciesNum[match(unlist(obsCatchSize$species), speciesList$species)]
  d$observedCatchSize <- obsCatchSize 
  
  # observed survey diet proportion by weight
  # for simulated diet! need different function for real diet that is by length
  # need length comp by age to get diet at length
  
  if(dattype == "sim") {
    
  svagelenbin <- survagelen %>%
    dplyr::mutate(species = Name) %>%
    dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
    dplyr::select(species, year, survey, agecl, lenbin, value) %>%
    dplyr::left_join(modbins) %>%
    dplyr::filter(modbin.min <= lenbin & lenbin < modbin.max) %>% #lenbin defined as lower
    dplyr::group_by(species, survey, year, agecl, sizebin) %>%
    dplyr::summarise(sumlen = sum(value)) %>%
    dplyr::group_by(species, year, sizebin) %>%
    dplyr::mutate(propage = sumlen/sum(sumlen)) #proportion of each agecl contributing to lengthbin
  
  obsSurvDiet <- survdiet %>%
    dplyr::mutate(species = Name) %>%
    dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
    dplyr::left_join(svagelenbin) %>%
    dplyr::mutate(dietpropage = value*propage) %>% #reweight diets for lengthbins
    dplyr::group_by(species, survey, year, sizebin, prey) %>%
    dplyr::summarise(dietsize = sum(dietpropage)) %>%
    dplyr::filter(prey %in% unique(modbins$species)) %>% #drops prey that aren't our modeled species
    tidyr::spread(prey, dietsize) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(sizebin)) 
  
  }
  
  # need code for real diet data here
  if(dattype == "real"){
    
    lenbintots <- survdietlen %>%
      dplyr::mutate(species = Name) %>%
      dplyr::filter(variable %in% c("totwt")) %>%
      dplyr::mutate(prey = "ALL") %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = variable, values_from = value) %>%
      dplyr::left_join(modbins) %>%
      dplyr::filter(modbin.min <= pdlen & pdlen < modbin.max) %>%
      dplyr::group_by(species, survey, year, sizebin) %>%
      dplyr::summarise(totwtsize = sum(totwt, na.rm = TRUE))
    
    obsSurvDiet <- survdietlen %>%
      dplyr::mutate(species = Name) %>%
      dplyr::filter(variable %in% c("meansw")) %>%
      tidyr::pivot_wider(names_from = variable, values_from = value) %>%
      dplyr::left_join(modbins) %>%
      dplyr::filter(modbin.min <= pdlen & pdlen < modbin.max) %>%
      dplyr::group_by(species, survey, year, sizebin, prey) %>%
      dplyr::summarise(dietwtsize = sum(meansw, na.rm = TRUE)) %>%
      dplyr::left_join(lenbintots) %>%
      dplyr::mutate(dietsize = dietwtsize/totwtsize) %>%
      dplyr::select(-dietwtsize, -totwtsize) %>%
      dplyr::left_join(focalprey, by=c("prey" = "SCIENTIFIC_NAME")) %>%
      dplyr::mutate(prey = Name) %>%
      dplyr::filter(prey %in% unique(modbins$species)) %>% #drops prey that aren't our modeled species
      dplyr::select(-SPECIES_ITIS, -Name) %>%
      tidyr::spread(prey, dietsize) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
      dplyr::filter(!is.na(sizebin)) 
  }
  
  # add back any modeled species that weren't prey so they show up as columns
  missedpreds <- setdiff(unique(modbins$species), names(obsSurvDiet)[-(1:4)])
  obsSurvDiet[missedpreds] <- NA_real_
  dcolorder <- c("species", "survey", "year", "sizebin", (sort.default(d$speciesList)))
  obsSurvDiet <- obsSurvDiet[dcolorder]
  
  obsSurvDiet <- obsSurvDiet %>% 
    dplyr::mutate(allotherprey = 1-rowSums(.[,-(1:4)], na.rm = TRUE)) %>%
    #dplyr::left_join(svalphamultlook) %>% #what is effective sample size for Dirichlet?
    #dplyr::mutate(inpN = max(surv_alphamult_n/100000, 5)) %>% #rescale this input for now
    dplyr::mutate(inpN = 100) %>% #hardcoded simulated sample size, revisit
    dplyr::select(survey, year, species, sizebin, inpN, everything()) %>% #, -surv_alphamult_n)
    dplyr::filter(species %in% names(predOrPrey)[predOrPrey==1]) #include only kept predators
  
  # use 0 for missing value
  obsSurvDiet <- obsSurvDiet %>%
    replace(is.na(.),0)
    
  d$Ndietprop_obs <- dim(obsSurvDiet)[1]
  obsSurvDiet$survey <- as.numeric(as.factor(obsSurvDiet$survey))
  obsSurvDiet$species <- speciesList$speciesNum[match(unlist(obsSurvDiet$species), speciesList$species)]
  obsSurvDiet$sizebin <- as.numeric(regmatches(obsSurvDiet$sizebin, regexpr("\\d+", obsSurvDiet$sizebin)))#take the number portion of sizebinN
  d$observedSurvDiet <- obsSurvDiet
  
  # observed effort by fleet (dummy, not used in estimation model) but fleenames used
  # also hardcoded for current 2 fleets--fixed dec 2022
  #fleetnames <- c("allbutcod", "codfleet")
  
  fleetnames <- names(fleetdef)[!names(fleetdef) %in% c("species", "qind")]
  
  obsEffort <- data.frame(1:d$Nyrs,
                          matrix(1, d$Nyrs, d$Nfleets))
  names(obsEffort) <- c("year", fleetnames)
    
  d$observedEffort <- t(obsEffort)
  d$fleetNames <- fleetnames #(names(obsEffort)[2:(d$Nfleets+1)])
  
  # observed temperature
  
  if(dattype == "sim") {
  # simulation using only fall: NOBA_BTS_fall_allbox_effic1
  obsTemp <- survtemp %>%
    dplyr::filter(survey=="BTS_fall_allbox_effic1") %>%
    dplyr::mutate(year = as.integer(year-fitstartyr)) %>% #year starts at 1
    dplyr::select(year, meantemp=value)
  }
  
  if(dattype == "real") {
    obsTemp <- survtemp %>%
      dplyr::mutate(year = as.integer(year-fitstartyr)) %>% #year starts at 1
      dplyr::select(year, meantemp=value)
  }
  
  d$observedTemperature <- t(obsTemp)
  
  # stomach weight
  # uses same length-age lookup as above to convert to length specific percapcons
  # takes mean across all years to get single vector
  
  if(dattype == "sim"){
  stomachContent <- percapcons %>%
    dplyr::mutate(species = Name) %>%
    dplyr::mutate(year = year-fitstartyr) %>% #year starts at 1
    tidyr::pivot_wider(-units, names_from = "variable", values_from = "value") %>%
    dplyr::left_join(svagelenbin) %>%
    dplyr::filter(!is.na(sizebin)) %>% # see note below
    dplyr::mutate(conspropage = totconsagecl*propage) %>% #reweight cons for lengthbins
    dplyr::mutate(totNpropage = totNagecl*propage) %>% #reweight nums for lengthbins
    dplyr::group_by(species, year, sizebin) %>%
    dplyr::summarise(intakesize = (sum(conspropage)/sum(totNpropage))*1000000) %>%
    tidyr::complete(sizebin) %>% #adds back any missing sizebin factor levels
    tidyr::spread(sizebin, intakesize) %>%
    dplyr::ungroup() %>% # stop here for individual year input
    dplyr::group_by(species) %>%
    dplyr::summarise(across(everything(), ~mean(.x, na.rm=T))) %>%
    dplyr::select(-year) %>%
    dplyr::rename_with(~gsub("bin", "class", .x)) %>% #stop here for unfilled
    dplyr::mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    tidyr::pivot_longer(-species) %>%
    dplyr::group_by(species) %>%
    tidyr::fill(value, .direction = "updown") %>%
    tidyr::pivot_wider()
  
  }
  
  if(dattype == "real") {
    # original hydra_sim parameterization, reorder with mskeyrun names
    stomachContent <- read.csv(paste0(path,"/intake_stomachContent.csv"),header=TRUE)
    stomachContent$species <- c("Spiny_dogfish",
                                "Winter_skate",
                                "Atlantic_herring",
                                "Atlantic_cod",
                                "Haddock",
                                "Yellowtail_flounder",
                                "Winter_flounder",
                                "Atlantic_mackerel",
                                "Silver_hake",
                                "Goosefish")
    
    stomachContent <-stomachContent[order(stomachContent$species),]
  }
  
  d$intakeStomach <- as.matrix(stomachContent[,2:(d$Nsizebins+1)])
  
  # recruitment parameters, all dummy; not used in estimation model
  stockRecruit <- matrix(1, 27, d$Nspecies,
                         dimnames = list(c("eggRicker_alpha",
                                           "eggRicker_shape",
                                           "eggRicker_beta",
                                           "DS_alpha",
                                           "DS_shape",
                                           "DS_beta",
                                           "gamma_alpha",
                                           "gamma_shape",
                                           "gamma_beta",
                                           "ricker_alpha",
                                           "ricker_shape",
                                           "ricker_beta",
                                           "BH_alpha",
                                           "BH_shape",
                                           "BH_beta",
                                           "shepherd_alpha",
                                           "shepherd_shape",
                                           "shepherd_beta",
                                           "hockey_alpha",
                                           "hockey_shape",
                                           "hockey_beta",
                                           "segmented_alpha",
                                           "segmented_shape",
                                           "segmented_beta",
                                           "sigma",
                                           "type",
                                           "stochastic"
                         ), c(sort.default(d$speciesList))
                         ))
  # change rectype input to 9=avg+devs
  stockRecruit[c("type"),] <- 9
  
  # change stochastic recruitment switch to off
  stockRecruit[c("stochastic"),] <- 0
  
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
  
  
  # sex ratio assumed 50:50 for all
  sexRatio <- matrix(0.5, 1, d$Nspecies,
                     dimnames = list("ratio", c(sort.default(d$speciesList))))
  d$sexRatio <- unlist(sexRatio)
  
  # recruitment covariate effects. # columns = d$Nrecruitment_cov
  # not used in estimation model
  rec_covEffects <- data.frame(names =c (sort.default(d$speciesList)), speciesEffect1 = 0)
  rownames(rec_covEffects) <- rec_covEffects$names
  rec_covEffects <- subset(rec_covEffects, select=-names)
  # ensure 0 numeric not character

  d$recruitCovEffects <- as.matrix(rec_covEffects)
  
  # fecundity all dummy variables not used in estimation model
  fecundity_d_h <- matrix(1, 2, d$Nspecies,
                          dimnames = list(c("d", "h"), c(sort.default(d$speciesList))))
  d$fecundityd <- unlist(fecundity_d_h["d",])
  d$fecundityh <- unlist(fecundity_d_h["h",])
  
  fecundity_Theta <- as.data.frame(matrix(0, d$Nspecies, d$Nsizebins, 
                                          dimnames = list(c(sort.default(d$speciesList)),  
                                                          c(paste0("sizeclass",seq(1:d$Nsizebins))))))
  d$fecundityTheta <- format(as.matrix(fecundity_Theta),digits=5)
  
  # maturity all dummy variables not used in estimation model
  maturity <- matrix(1, 2, d$Nspecies,
                     dimnames = list(c("nu", "omega"), c(sort.default(d$speciesList))))
  d$maturityNu <- unlist(maturity["nu",])
  d$maturityOmega <- unlist(maturity["omega",])
  
  # dummy variables, no covariate effects in estimation model
  maturity_covEffects <- data.frame(name = c(sort.default(d$speciesList)), speciesEffect1 = 0)
  rownames(maturity_covEffects) <- maturity_covEffects$name
  maturity_covEffects <- subset(maturity_covEffects, select=-name)
  
  d$maturityCovEffects <- as.matrix(maturity_covEffects)
  
  # growth FIX THIS TO BE ESTIMATED IN A FUNCTION HERE
  if(dattype == "sim"){
    growth <- read.csv(paste0(path,"/growth_species_NOBA_era1_BTS_fall_allbox_effic1.csv"),header=TRUE,row.names=1)
  }
  
  if(dattype == "real"){
    #original hydra_sim parameters, change to mskeyrun names and order
    growth <- read.csv(paste0(path,"/growth_species.csv"),header=FALSE,row.names=1)
    names(growth) <- c("Spiny_dogfish",
                       "Winter_skate",
                       "Atlantic_herring",
                       "Atlantic_cod",
                       "Haddock",
                       "Yellowtail_flounder",
                       "Winter_flounder",
                       "Atlantic_mackerel",
                       "Silver_hake",
                       "Goosefish")
    growth <- growth[-1,]
    growth <- growth[,order(names(growth))]
    
  }
  d$growthPsi <- unlist(growth["psi",])
  d$growthKappa <- unlist(growth["kappa",])
  d$growthLinf <- unlist(growth["Linf",])
  d$growthK <- unlist(growth["k",])
  d$growthType <- unlist(growth["growthType",])
  
  d$growthType[d$growthType==4] <- 3 #change to no covariates
  d$growthType[d$growthType==2] <- 1 #change to no covariates
  
  # dummy variables, no covariate effects in estimation model
  growth_covEffects <- data.frame(name = c(sort.default(d$speciesList)), speciesEffect1 = 0)
  rownames(growth_covEffects) <- growth_covEffects$name
  growth_covEffects <- subset(growth_covEffects, select=-name)
  
  d$growthCovEffects <- as.matrix(growth_covEffects)
  
  # intake: fixed parameters for all species
  intake <- rbind(alpha = rep(0.004, d$Nspecies),
                  beta = rep(0.11, d$Nspecies))
  colnames(intake) <- sort.default(d$speciesList)
  intake <- as.data.frame(intake)
  
  d$intakeAlpha <- unlist(intake["alpha",])
  d$intakeBeta <- unlist(intake["beta",])
  
  # M1, annual value to be rescaled in hydra after Nstepsyr is calculated
  # Dec 2022 moved to pin and reformatted for estimation
  # M1 <- as.data.frame(matrix(0.1, d$Nspecies, d$Nsizebins, 
  #                            dimnames = list(c(sort.default(d$speciesList)),  
  #                                            c(paste0("sizeclass",seq(1:d$Nsizebins)))))) 
  # d$M1 <- as.matrix(M1)
  
  #foodweb, object created above
  foodweb <- foodweb %>% tibble::column_to_rownames("prey")
  #rownames(foodweb) <- foodweb[,1]
  #foodweb[,1] <- NULL
  d$foodweb <- as.matrix(foodweb)
  
  #M2 size preference
  M2sizePref <- rbind(mu = as.numeric(rep(0.5, d$Nspecies)),
                      sigma = as.numeric(rep(2), d$Nspecies))
  colnames(M2sizePref) <- sort.default(d$speciesList)
  M2sizePref <- as.data.frame(M2sizePref)
  
  d$M2sizePrefMu <- as.matrix(M2sizePref["mu",])
  d$M2sizePrefSigma <- as.matrix(M2sizePref["sigma",])
  
  # fishery catchability indicator (q's)
  indicatorFisheryqs<- fleetdef %>%
    dplyr::select(dplyr::all_of(fleetnames))
  d$indicatorFisheryq<- as.matrix(t(indicatorFisheryqs))
  
  
  #fishery/fleet selectivity MOVED TO PIN FILE AND REFORMATTED
  # fisherySelectivityc<- read.csv(paste0(path,"/fishing_selectivityc_NOBA.csv"),header=TRUE,row.names = 1)
  # d$fisherySelectivityc <- unlist(as.matrix(t(fisherySelectivityc)))
  # fisherySelectivityd<- read.csv(paste0(path,"/fishing_selectivityd_NOBA.csv"),header=TRUE,row.names = 1)
  # d$fisherySelectivityd <- unlist(as.matrix(t(fisherySelectivityd)))
  
  # NOTHING BELOW USED IN ESTIMATION, READ IN PLACEHOLDERS FROM CSV
  # But these will break with different data sources so update with calculations
  # B0 - equilibrium biomass faked in with initial year survey biomass
  B0 <- survindex %>%
    dplyr::filter(year==min(year),
                  variable=="biomass") %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(B0 = mean(value)) %>%
    dplyr::select(Name, B0) 
  
  d$B0 <- B0$B0
  
  # assessment thresholds + exploitations (fixed for sims, not used in est)
  assessmentThresholds <-  read.csv(paste0(path,"/assessmentThresholds.csv"),header=TRUE)
  d$thresholds <- assessmentThresholds$thresholds
  thresholds <- d$thresholds[d$thresholds< 1]
  d$Nthresholds <- length(d$thresholds)
  d$exploitationOptions <- assessmentThresholds[,2:dim(assessmentThresholds)[2]]
  d$minMaxThresholds <- c(min(thresholds),max(thresholds))
  
  # additional level added to species specific threshold
  assessmentThresholdsSpecies <-  matrix(0, 1, d$Nspecies, 
                                         dimnames = c("additionalThreshold", list(sort.default(d$speciesList))))
  d$thresholdSpecies <- unlist(assessmentThresholdsSpecies)
  
  #scaled Efort - not used
  scaledEffort <-  as.data.frame(matrix(1, 1, d$Nspecies, 
                                        dimnames = c("", list(c(sort.default(d$speciesList))))))
  d$scaledEffort <- unlist(scaledEffort)
  
  # discard coefficient - prob of discard
  discardCoef <-  matrix(0, d$Nspecies*d$Nfleets, d$Nsizebins)
  
  d$discardCoef <- (unlist(as.matrix(discardCoef)))
  
  # discard survival probability | discard
  discardSurvival <-  matrix(1, d$Nspecies*d$Nfleets, d$Nsizebins)
  
  d$discardSurvival <- (unlist(as.matrix(discardSurvival)))
  
  # Autoregressive parameters for alternative error structure (AR) for survey, recruitment, catch
  ARparameters<- read.csv(paste0(path,"/observation_error.csv"),header=TRUE)
  d$ARParameters <- unlist(ARparameters)
  
  
  # residence time in the modeled area
  residenceTime <- as.data.frame(matrix(1, 1, d$Nspecies, 
                                        dimnames = c("", list(c(sort.default(d$speciesList))))))
  d$residenceTime <- unlist(residenceTime)
  
  # mortality rate outside the modeled area
  areaMortality <- as.data.frame(matrix(0, 1, d$Nspecies, 
                                        dimnames = c("", list(c(sort.default(d$speciesList))))))
  d$areaMortality <- unlist(areaMortality)
  
  # add missing vectors from sim not used in est
  if(is.null(d$fleetMembership)) d$fleetMembership <- fleetdef$qind #WARNING:assumes guilds=species
  if(is.null(d$minExploitation)) d$minExploitation <- rep(1e-05, d$Nfleets)
  if(is.null(d$maxExploitation)) d$maxExploitation <- rep(0.999, d$Nfleets)
  
  
  # replace NA
  d <- lapply(d, function(x) {x[is.na(x)] <- -999; x})
  
  
  return(d)
  
}



get_PinData_msk <- function(dattype,
                            nlenbin,
                            speciesList = d$speciesList,
                            bindef,
                            startpars,
                            Nyrs = d$Nyrs,
                            Nfleets = d$Nfleets,
                            Nsurveys = d$Nsurveys,
                            Nareas = d$Nareas,
                            Npreds = sum(d$predOrPrey),
                            fqind = d$indicatorFisheryq,
                            observedCatch = d$observedCatch,
                            observedBiomass = d$observedBiomass,
                            path){
  # Stipulate all information required for the pin data file
  p <- list()
  # path to data
  # list of species and guilds (Functional Groups)'
  #speciesList <- focalspp[order(unique(focalspp$Name)),]
  Nspecies <- length(speciesList)
  
  modbins <- bindef$modbins
  
  # starting length structure
  if(dattype == "sim"){
  Y1N <- startpars %>%
    dplyr::filter(variable=="Natlen") %>%
    dplyr::left_join(modbins, by=c("Name" = "species")) %>%
    dplyr::filter(modbin.min <= lenbin & lenbin < modbin.max) %>% #lenbin defined as lower
    dplyr::group_by(Name, sizebin) %>%
    dplyr::summarise(sumlen = log(sum(value)/1000000)) %>%
    tidyr::complete(sizebin) %>% #adds back any missing sizebin factor levels
    tidyr::spread(sizebin, sumlen) %>%
    dplyr::ungroup() %>%
    replace(is.na(.),0) %>% # use 0 for missing value in pin file
    dplyr::select(-Name) %>%
    as.matrix()
  }
  
  if(dattype == "real"){
    # try original parameterization
    #Y1N <- read.csv(paste0(path,"/observation_Y1N.csv"),header=TRUE,row.names=1)
    # this is in the paper supplement, the above must have come afterwards?
    # initial N millions
    Y1Npaper <- list("Spiny dogfish" = c(
                      44.0358,
                      39.76466,
                      29.46116,
                      18.65818,
                      8.19665),
                      "Winter skate" = c(
                      116.86965,
                      131.41996,
                      118.27867,
                      102.31659,
                      112.13381),
                      "Atlantic herring" = c(
                      1229.2325,
                      1229.2325,
                      1307.77136,
                      1307.77136,
                      2724.31378),
                      "Atlantic cod" = c(
                      75.24274,
                      64.75629,
                      42.76626,
                      50.49759,
                      16.06118),
                      "Haddock" = c(
                      129.57077,
                      320.35115,
                      596.82624,
                      588.06217,
                      725.45692),
                      "Yellowtail flounder" = c(
                      119.02818,
                      202.48411,
                      227.44696,
                      272.49665,
                      675.00921),
                      "Winter flounder" = c(
                      57.89208,
                      86.18024,
                      94.94266,
                      112.72898,
                      161.48656),
                      "Atlantic mackerel" = c(
                      172.99207,
                      220.51843,
                      198.87867,
                      177.85968,
                      161.15555),
                      "Silver hake" = c(
                      135.37161,
                      153.75246,
                      111.32931,
                      76.79369,
                      101.00855),
                      "Monkfish" = c(
                      10.54952,
                      10.73453,
                      8.87141,
                      9.40031,
                      4.09328)
    )
    
    Y1N <- as.data.frame(do.call(rbind, Y1Npaper))
    names(Y1N) <- c("sizeclass1","sizeclass2","sizeclass3","sizeclass4","sizeclass5")
    
    rownames(Y1N) <- c("Spiny_dogfish",
                     "Winter_skate",
                     "Atlantic_herring",
                     "Atlantic_cod",
                     "Haddock",
                     "Yellowtail_flounder",
                     "Winter_flounder",
                     "Atlantic_mackerel",
                     "Silver_hake",
                     "Goosefish")
    
    Y1N <-Y1N[order(rownames(Y1N)),]
    
    Y1N <- log(Y1N[,] + 1e-8)
  }
  
  rownames(Y1N) = sort.default(speciesList)
  
  p$Y1N <- Y1N
  
  p$recalpha <- rep(0.0, Nspecies)
  
  p$recshape <- rep(0.0, Nspecies)
  
  p$recbeta <- rep(0.0, Nspecies)
  
  # Avg recruitment and deviations
  if(dattype == "sim"){
  redundantAvgRec <- startpars %>%
    dplyr::filter(variable=="AvgRec") %>%
    dplyr::select(Name, ln_avgrec=value) %>%
    tibble::column_to_rownames("Name") 
  }
  
  if(dattype == "real"){
    redundantAvgRec <- read.csv(paste0(path,"/redundantAvgRecPinData.csv"),header=FALSE,row.names=1)
    names(redundantAvgRec) <- c("Spiny_dogfish",
                       "Winter_skate",
                       "Atlantic_herring",
                       "Atlantic_cod",
                       "Haddock",
                       "Yellowtail_flounder",
                       "Winter_flounder",
                       "Atlantic_mackerel",
                       "Silver_hake",
                       "Goosefish")
    redundantAvgRec <- redundantAvgRec[-1,]
    redundantAvgRec <- redundantAvgRec[,order(names(redundantAvgRec))]
    
  }
  
  
  p$redundantAvgRec <- t(redundantAvgRec)
  
  # recdevs start at 0 and are estimated for year 2 on, so Nrys-1
  redundantRecDevs <-  matrix(0, Nspecies, Nyrs-1,
                              dimnames = list(c(sort.default(speciesList), NULL)))
  
  p$redundantRecDevs <- unlist(redundantRecDevs)
  
  p$ln_recsigma <- rep(1.0, Nspecies)
  
  # starting values for F are log exploitation rates based on observed fleets
  sumcatchfleet <- observedCatch %>%
    dplyr::group_by(fishery) %>%
    dplyr::summarise(totcatch = sum(catch))
  
  fleetspp <- observedCatch %>%
    dplyr::distinct(species, fishery, .keep_all = TRUE) %>%
    dplyr::select(species, fishery)
  
  sumbiofleet <- observedBiomass %>%
    dplyr::left_join(fleetspp) %>%
    dplyr::group_by(fishery) %>%
    dplyr::summarise(sppbio = sum(biomass))
  
  ln_avgF <- merge(sumcatchfleet, sumbiofleet) %>%
    dplyr::mutate(avgF = totcatch/sppbio,
                  lnavgF = log(avgF))
  
  p$ln_avgF <- t(ln_avgF$lnavgF)
   
  # Fdevs start at 0 
  Fdevs <- matrix(0, Nfleets, Nyrs)
  
  p$Fdevs <- Fdevs
  
  # fishsel_pars
  fishsel_c <- rep(1, Nfleets)
  fishsel_d <- rep(1, Nfleets)
  
  p$fishsel_pars <- matrix(c(fishsel_c, fishsel_d), ncol = Nfleets, byrow = TRUE)
  
  
  # ln_fishery_q: ln_fishery_q(1,Nqpars) (= sum of fishery_indicator_q - Nfleet*Narea):
  # formula  needs adjustment now that we have 2 fisheries in qind, can't sum
  # count of non-zero entries in fquind would do it
  Nqpars <-length(fqind[fqind>0]) - Nfleets*Nareas
  
  # fishery catchability (q's)
  fisheryqs<- rep(1, Nqpars)
  
  p$fisheryq <- read.csv(paste0(path,"/fishing_q_NOBA.csv"),header=TRUE,row.names = 1) #old for -sim.pin
  
  p$ln_fishery_q <- log(fisheryqs) # for estimation pin file
  
  
  # ln_survey_q:
  surveyqs<- matrix(1, Nsurveys, Nspecies)
  
  p$surveyq <-  surveyqs
  
  p$ln_survey_q <- log(surveyqs)
  

  # survey_selpars
  survsel_c <- rep(1, Nsurveys)
  survsel_d <- rep(1, Nsurveys)
  
  p$survsel_pars <- matrix(c(survsel_c,survsel_d), ncol = Nsurveys, byrow = TRUE)
  
  # M1 init_matrix ln_M1ann(1,Nareas,1,Nspecies,m1_phase)
  M1ann <- rep(0.1, Nspecies)
  p$ln_M1ann <- log(M1ann)
  
  # ln_otherFood base // amount of other food included in the M2 term for the base (predator 1)
  p$ln_otherFood <- 21
  
  # other food devs //deviation from base other food for predators 2+  (same other food for all size classes of each predator)
  p$otherFood_dev <- rep(0, Npreds-1)
  
  # fishery and survey sigmas not used in estimation
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
