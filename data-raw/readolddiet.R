# look at original GB Hydra diet data GB_diets_10pred_allYears.csv
library(magrittr)

dietdat <- readr::read_csv(here::here("data-raw/GB_diets_10pred_allYears.csv"))

lookup <- mskeyrun::focalSpecies %>% 
  dplyr::select(SVSPP, SPECIES_ITIS, modelName, SCIENTIFIC_NAME) %>% 
  dplyr::distinct()

preycat <- unique(dietdat$preycat)
# [1] "CLUHAR" "GADMOR" "LIMFER" "LOPAME" "MELAEG" 
#      "MERBIL" "OTHER"  "OTHFIS" "PLEAME" "RAJERI" 
#      "SCOSCO" "SQUACA" "UNIFIS"

preyname <- c("Atlantic_herring", "Atlantic_cod", "Yellowtail_flounder", "Goosefish", "Haddock", 
              "Silver_hake", "*Otherfood*", "*Otherfood*", "Winter_flounder", "*Little skate*",
              "Atlantic_mackerel", "Spiny_dogfish", "*Otherfood*")

preylook <- data.frame(preycat, preyname)

predpreyorig <- dietdat %>%
  dplyr::left_join(lookup, by=c("svspp"="SVSPP")) %>%
  dplyr::left_join(preylook) %>%
  dplyr::select(SEASON, preyname, relmsw, modelName) %>%
  dplyr::arrange(modelName) %>%
  dplyr::group_by(SEASON, modelName, preyname) %>%
  dplyr::summarise(relmsw = sum(relmsw)) %>%
  tidyr::pivot_wider(names_from = modelName, values_from = relmsw) %>%
  dplyr::arrange(SEASON, factor(preyname, levels = c(sort(lookup$modelName), 
                                                     "*Little skate*", 
                                                     "*Otherfood*")))
  
  
flextable::flextable(predpreyorig) %>% flextable::autofit()


