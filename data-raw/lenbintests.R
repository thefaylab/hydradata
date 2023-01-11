library(dplyr)
library(ggplot2)
library(mskeyrun)



quantsurv <- realSurveyLennumcomp %>%
  filter(variable == "numbers") %>%
  group_by(Name) %>%
  summarise(minlen = min(lenbin, na.rm = TRUE),
            maxlen = max(lenbin, na.rm = TRUE),
            lenqant = quantile(lenbin, weights=value, 
                          c(0.05, 0.10, 0.5, 0.90, 0.95)), 
            quant = c(0.05, 0.10, 0.5, 0.90, 0.95)) %>%
  tidyr::pivot_wider(names_from = quant, names_prefix = "q",values_from = lenqant)

quantfish <- realFisheryLencompRaw %>%
  filter(variable == "abundance",
         !is.na(value)) %>%
  group_by(Name) %>%
  summarise(minlen = min(lenbin, na.rm = TRUE),
            maxlen = max(lenbin, na.rm = TRUE),
            lenqant = quantile(lenbin, weights=value, 
                               c(0.05, 0.10, 0.5, 0.90, 0.95)), 
            quant = c(0.05, 0.10, 0.5, 0.90, 0.95)) %>%
  tidyr::pivot_wider(names_from = quant, names_prefix = "q",values_from = lenqant)

# use q.0.10 survey for lower limit, q.0.90 fishery for upper limit
# divide that range into 5 equal bins, allow bin 1 to start at 0, bin 5 to extend to max
# so binwidths object still starts at 0 and goes to max, bin 1 and last different widths

Nspecies <- length(unique(quantsurv$Name))
Nsizebins <- 5

# testing
# mixquant <- quantsurv %>%
#   select(Name, q0.1, maxlens) %>%
#   left_join(quantfish %>% select(Name, q0.9, maxlenf)) %>%
#   mutate(range = q0.9 - q0.1,
#          binwidth = ceiling(range/Nsizebins),
#          binwidth1 = ceiling(q0.1 + binwidth),
#          binwidthN = ceiling(binwidth + max(maxlens, maxlenf) - q0.9),
#          check = binwidth1 + binwidth * (Nsizebins - 2) + binwidthN)

# summary object for function
# takes min lower quant and max upper quant
maxLrange <- dplyr::bind_rows(quantsurv, quantfish) %>%
  dplyr::group_by(Name) %>%
  dplyr::summarise(minqlo = min(q0.1),
                   maxqhi = max(q0.9),
                   maxmaxL = max(maxlen))

#function to sub in
quantbinwidths <- function(qlo, qhi, Lmax, Nsizebins){
  range  <- qhi-qlo
  bindef <- max(1,ceiling(range/Nsizebins))
  binwidths <- c(ceiling(qlo + bindef),
                 rep(bindef, Nsizebins-2),
                 ceiling(bindef + Lmax-qhi))
  return(binwidths)
}

newbins <- maxLrange %>%
  dplyr::group_by(Name) %>%
  dplyr::group_modify(~broom::tidy(quantbinwidths(.$minqlo,
                                                  .$maxqhi,
                                                  .$maxmaxL, 
                                                  Nsizebins))) %>%
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

plotbins <- modbins %>%
  rename(Name = species)

plotsurvlen <- realSurveyLennumcomp %>% 
  filter(variable == "numbers") %>%
  left_join(quantsurv) %>%
  ggplot(aes(x=lenbin, y=value)) +
  geom_bar(stat = "identity", width = 1, fill = "darkgreen") +
  theme_bw() +
  facet_wrap(~Name, scales = "free") +
  ggtitle("Survey lengths with proposed bins") +
  #sapply(quantsurv$q0.1, function(xint) geom_vline(aes(xintercept = xint))) 
  
  geom_vline(data=filter(plotbins, Name=="Atlantic cod"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Atlantic herring"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Atlantic mackerel"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Goosefish"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Haddock"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Atlantic herring"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Silver hake"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Spiny dogfish"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Winter flounder"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Winter skate"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Yellowtail flounder"),
             aes(xintercept = modbin.max)) 
  


plotfishlen <- realFisheryLencompRaw %>% 
  filter(variable == "abundance") %>%
  left_join(quantfish) %>%
  ggplot(aes(x=lenbin, y=value)) +
  geom_bar(stat = "identity", width = 1, fill = "blue") +
  theme_bw() +
  facet_wrap(~Name, scales = "free") +
  ggtitle("Fishery lengths with proposed bins") +
  
  geom_vline(data=filter(plotbins, Name=="Atlantic cod"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Atlantic herring"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Atlantic mackerel"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Goosefish"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Haddock"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Atlantic herring"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Silver hake"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Spiny dogfish"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Winter flounder"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Winter skate"),
             aes(xintercept = modbin.max)) +
  geom_vline(data=filter(plotbins, Name=="Yellowtail flounder"),
             aes(xintercept = modbin.max)) 

