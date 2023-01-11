library(dplyr)
library(ggplot2)
library(mskeyrun)



quantsurv <- realSurveyLennumcomp %>%
  filter(variable == "numbers") %>%
  group_by(Name) %>%
  summarise(lenqant = quantile(lenbin, weights=value, 
                          c(0.05, 0.25, 0.5, 0.75, 0.95)), 
            q = c(0.05, 0.25, 0.5, 0.75, 0.95))

quantfish <- realFisheryLencompRaw %>%
  filter(variable == "abundance",
         !is.na(value)) %>%
  group_by(Name) %>%
  summarise(lenqant = quantile(lenbin, weights=value, 
                               c(0.05, 0.25, 0.5, 0.75, 0.95)), 
            q = c(0.05, 0.25, 0.5, 0.75, 0.95))

fishlen <- ggplot(realFisheryLencompRaw %>% filter(variable == "abundance"),
                  aes(x=lenbin, y=value)) +
  geom_bar(stat = "identity", width = 1) +
  #sapply(quantfish, function(lenquant) geom_vline(aes(xintercept = lenquant))) +
  facet_wrap(~Name, scales = "free")


survlen <- ggplot(realSurveyLennumcomp %>% filter(variable == "numbers"),
                  aes(x=lenbin, y=value)) +
  geom_bar(stat = "identity", width = 1) +
  #sapply(quantsurv, function(lenquant) geom_vline(aes(xintercept = lenquant))) +
  facet_wrap(~Name, scales = "free")
