# test mskeyrun:catchAtLengthProportions vs aggregated catchIndex

focalspp <- mskeyrun::focalSpecies %>%
  dplyr::filter(modelName != "Pollock") %>% # not using in these models
  dplyr::mutate(Name = modelName) 

totcatch <- mskeyrun::catchIndex %>%
  dplyr::left_join(focalspp %>% dplyr::mutate(NESPP3 = as.integer(NESPP3)))

gearcatch <- mskeyrun::catchAtLengthProportions %>%
  dplyr::left_join(focalspp, by=c("species_itis" = "SPECIES_ITIS"))
