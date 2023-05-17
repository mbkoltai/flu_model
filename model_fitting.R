# source functions, call libraries
source("fcns/fcns.R")
# load data
source("fcns/load_flunet.R")

# FluEvidenceSynthesis
library(fluEvidenceSynthesis)

# Engl ILI dataset
as.data.frame(ili$ili) %>% mutate(week=1:52) %>% pivot_longer(!week,names_to = "age") %>% 
  ggplot() + geom_line(aes(x=week,y=value)) + facet_wrap(~age,scales="free_y") + 
  theme_bw() + standard_theme

### ####

cromer_samples_by_age <- qread("/home/lshmk17/Downloads/cromer_samples_by_age.qs")
cromer_summ_stat <- cromer_samples_by_age %>% group_by(risk_group,subtype,outcome,age) %>% 
  summarise(median=median(value),mean=mean(value),
            ci95_l=quantile(value,probs=0.025),ci95_u=quantile(value,probs=0.975),
            ci50_l=quantile(value,probs=0.25),ci50_u=quantile(value,probs=0.75))

cromer_summ_stat %>% filter(age<99) %>% # mutate(across(c(mean,ci95_l,ci95_u,ci50_l,ci50_u), ~ ifelse(. == 0, NA, .))) %>%
  ggplot(aes(x=age)) + facet_grid(outcome~risk_group,scales="free_y") +
  geom_line(aes(y=mean*100,color=subtype),size=1.2) + 
  geom_ribbon(aes(ymin=ci50_l*100,ymax=ci50_u*100,fill=subtype),alpha=1/5) +
  ylab("probability (%)") + # scale_y_log10() + 
  theme_bw() + standard_theme # 
