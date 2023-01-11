library(dplyr)
library(ggplot2)
library(mskeyrun)



quantsurv <- realSurveyLennumcomp %>%
  filter(variable == "numbers") %>%
  group_by(Name) %>%
  summarise(lenqant = quantile(lenbin, weights=value, 
                          c(0.05, 0.10, 0.5, 0.90, 0.95)), 
            quant = c(0.05, 0.10, 0.5, 0.90, 0.95)) %>%
  tidyr::pivot_wider(names_from = quant, names_prefix = "q",values_from = lenqant)

plotsurvlen <- realSurveyLennumcomp %>% 
  filter(variable == "numbers") %>%
  left_join(quantsurv) %>%
  ggplot(aes(x=lenbin, y=value)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~Name, scales = "free") #+ 
  # geom_vline(data=filter(quantsurv, Name=="Atlantic cod"),
  #            xintercept = c("q0.05", "q0.95")) +
  # geom_vline(data=filter(quantsurv, Name=="Atlantic herring"),
  #            xintercept = c("q0.05", "q0.95")) +
  # geom_vline(data=filter(quantsurv, Name=="Atlantic mackerel"),
  #            xintercept = c("q0.05", "q0.95"))


quantfish <- realFisheryLencompRaw %>%
  filter(variable == "abundance",
         !is.na(value)) %>%
  group_by(Name) %>%
  summarise(lenqant = quantile(lenbin, weights=value, 
                               c(0.05, 0.10, 0.5, 0.90, 0.95)), 
            quant = c(0.05, 0.10, 0.5, 0.90, 0.95)) %>%
  tidyr::pivot_wider(names_from = quant, names_prefix = "q",values_from = lenqant)

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

# use q.0.10 survey for lower limit, q.0.90 fishery for upper limit
# divide that range into 5 equal bins, allow bin 1 to start at 0, bin 5 to extend to max
# so binwidths object still starts at 0 and goes to max, bin 1 and last different widths

Nsizebins <- 5

mixquant <- quantsurv %>%
  select(Name, q0.1) %>%
  left_join(quantfish %>% select(Name, q0.9)) %>%
  mutate(range = q0.9 - q0.1,
         binwidth = ceiling(range/Nsizebins),
         binwidth1 = ceiling(q0.1 + binwidth))
