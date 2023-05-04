# FluEvidenceSynthesis

# source functions, call libraries
source("fcns/fcns.R")

as.data.frame(ili$ili) %>% mutate(week=1:52) %>% pivot_longer(!week,names_to = "age") %>% 
  ggplot() + geom_line(aes(x=week,y=value)) + facet_wrap(~age,scales = "free_y") + 
  theme_bw() + standard_theme
