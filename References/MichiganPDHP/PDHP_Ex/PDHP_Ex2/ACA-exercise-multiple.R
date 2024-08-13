#---------------------------------------------------
#      Empirical Exercise
#      Effect of ACA on Health insurance
#---------------------------------------------------
#-----------------------------------------------------------------------------

# Load packages
#-----------------------------------------------------------------------------
# Libraries
library(readstata13)
library(statar)
library(ggplot2)
library(ggpubr)

#Download latest version of did package and load it
#remotes::install_github("bcallaway11/did")
library(did)
# Use here package to facilitate relative paths
library(here)
# Use these for data manipulation, and plots
library(tidyverse)
library(ggplot2)
library(scales)
library(ggrepel)
library(dplyr)
library(bacondecomp) 
library(TwoWayFEWeights)
library(fixest)
library(glue)
#---------------------------------------------------------------------------------------
# Set ggplot theme
theme_set(
  #theme_clean() + 
  theme_classic() +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      legend.background = element_rect(color = "white"),
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.spacing = unit(10, "lines"))
)
#---------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Load data
acafs <- data.frame(read.dta13(here("data",'ehec_data.dta')))
#-----------------------------------------------------------------------------
# Do some data manipulations
# Ensure that stfips are numeric
acafs$stfips2 <- as.numeric(acafs$stfips)
# Make sure year is numeric (and correct)
acafs$year2 <- 1985 + as.numeric(acafs$year)
# If year of adoption is missing, set it to Infinity ( zero or any large number  would also work with the did package)
acafs$yexp2[is.na(acafs$yexp2)] <- Inf
# Create treatment dummy - 1 if treated by that year, 0 otherwise
acafs$treated <- as.numeric(acafs$year2 >= acafs$yexp2)
#-----------------------------------------------------------------------------
# Create subset of data withour never-treted
acafs_no_never <- subset(acafs, acafs$yexp2!=Inf)

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Start the analysis
#---------------------------------------------------------------------------------------
# Get TWFE coefficient
twfe <- fixest::feols(dins ~ treated| stfips2 + year2, 
                      data = acafs,
                      weights = ~W, 
                      cluster = ~stfips2)

summary(twfe)
#---------------------------------------------------------------------------------------
# Get Bacon decomposition (without weights bc not supported in the R package)
df_bacon <- bacon(dins ~ treated,
                  data = acafs,
                  id_var = "stfips2",
                  time_var = "year2")
#---------------------------------------------------------------------------------------
# Get de Chaisemartin and D'Haultfoeuille Decomposition
dCDH_decomp <- twowayfeweights(
  df = acafs, 
  Y = "dins", 
  G = "stfips2",
  T = "year2", 
  D ="treated",
  cmd_type =  "feTR"
)

# Weakly Positive weights
dCDH_positive <- sum(dCDH_decomp$weight[dCDH_decomp$weight>=0])

# Negative weights
dCDH_negative <- sum(dCDH_decomp$weight[dCDH_decomp$weight<0])
#---------------------------------------------------------------------------------------
# Get TWFE event study coefficients

# create event times (rel_year)
acafs <- acafs %>% 
  mutate(
    rel_year = year2 - yexp2
  )

# Create min_year to omit it in the TWFE regression
min_year <- min(acafs$rel_year * (acafs$rel_year != -Inf), na.rm = T)
# Formula we will use
formula_twfe_es <- as.formula(glue::glue("dins ~ i(rel_year, ref=c(-1, {min_year}, -Inf)) | stfips2 + year2"))

# estimate the TWFE coefficients
twfe_es <- fixest::feols(formula_twfe_es, data = acafs, cluster = ~stfips2)
summary(twfe_es)

# Put the TWFE coefficients in a tibble that is easy to plot
twfe_es <- broom::tidy(twfe_es) %>%
  filter(str_detect(term, "rel_year::")) %>% 
  mutate(
    rel_year = as.numeric(str_remove(term, "rel_year::")),
  ) %>%
  filter(rel_year <= 5 & rel_year >= -5) %>%
  select(event.time = rel_year, 
         estimate,
         std.error) %>%
  add_row(event.time = -1, 
          estimate = 0,
          std.error =0)  %>%
  mutate( point.conf.low  = estimate - stats::qnorm(1 - .05/2) * std.error,
          point.conf.high = estimate + stats::qnorm(1 - .05/2) * std.error)


twfe_es


#---------------------------------------------------------------------------------------
# Plots for TWFE
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
p_es_twfe<- ggplot(data = twfe_es,
                   mapping = aes(x = event.time, y = estimate)) +
  geom_line(size = 0.5, alpha = 2, colour = "black") +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = 0, colour="black",  linetype = "dotted")+
  #geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = FALSE, linetype= 1, size = 1.1,
  #                color = "red")+
  geom_pointrange(aes(ymin = point.conf.low, ymax = point.conf.high), show.legend = FALSE, size = 1.1)+
  xlab("Event time") +
  ylab("Event Study Estimate") +
  ylim(range(-0.03, 0.11))+
  scale_x_continuous(breaks = c(-5:5)) +
  #scale_y_continuous(breaks = seq(-600, 200, 200), limits = c(-700,200))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  ) + 
  annotate(geom="text", x=3, y=-0.01, label="Static TWFE Coefficient:
           0.0740 (0.013)",
           color="black")
p_es_twfe




#---------------------------------------------------------------------------------------
# Callaway and Sant'Anna (2021) procedure
#---------------------------------------------------------------------------------------
# Use never-treated as comparison group
atts_never <- did::att_gt(yname = "dins", # name of the LHS variable
                          tname = "year2", # name of the time variable
                          idname = "stfips2", # name of the id variable
                          gname = "yexp2", # name of the first treatment period variable
                          data = acafs, # name of the data
                          xformla = NULL,
                          weightsname = "W",
                          est_method = "reg", # estimation method. "dr" means doubly robust
                          control_group = "nevertreated", # set the control group which is either "nevertreated" or "notyettreated" 
                          bstrap = TRUE, # if TRUE compute boostrapped SE
                          biters = 1000, # number of boostrap interations
                          print_details = FALSE, # if TRUE, print detailed results
                          panel = TRUE) # whether the data is panel or repeated cross-sectional
summary(atts_never)
agg_effects_es_never <- aggte(atts_never, type = "dynamic", min_e = -5, max_e = 5)
summary(agg_effects_es_never)


### Control groups: Not yet treated

# Use not-yet-treated as comparison group
atts_ny <- did::att_gt(yname = "dins", # name of the LHS variable
                       tname = "year2", # name of the time variable
                       idname = "stfips2", # name of the id variable
                       gname = "yexp2", # name of the first treatment period variable
                       data = acafs, # name of the data
                       xformla = NULL,
                       weightsname = "W",
                       est_method = "reg", # estimation method. "dr" means doubly robust
                       control_group = "notyettreated", # set the control group which is either "nevertreated" or "notyettreated" 
                       bstrap = TRUE, # if TRUE compute boostrapped SE
                       biters = 1000, # number of boostrap interations
                       print_details = FALSE, # if TRUE, print detailed results
                       panel = TRUE) # whether the data is panel or repeated cross-sectional
summary(atts_ny)
agg_effects_es_ny <- aggte(atts_ny, type = "dynamic", min_e = -5, max_e = 5)
summary(agg_effects_es_ny)

#---------------------------------------------------------------------------------------

# Use not-yet-treated as comparison group (drop all never-treated)
atts_ny2 <- did::att_gt(yname = "dins", # name of the LHS variable
                        tname = "year2", # name of the time variable
                        idname = "stfips2", # name of the id variable
                        gname = "yexp2", # name of the first treatment period variable
                        data = acafs_no_never, # name of the data
                        xformla = NULL,
                        weightsname = "W",
                        est_method = "reg", # estimation method. "dr" means doubly robust
                        control_group = "notyettreated", # set the control group which is either "nevertreated" or "notyettreated" 
                        bstrap = TRUE, # if TRUE compute boostrapped SE
                        biters = 1000, # number of boostrap interations
                        print_details = FALSE, # if TRUE, print detailed results
                        panel = TRUE) # whether the data is panel or repeated cross-sectional
summary(atts_ny2)
agg_effects_es_ny2 <- aggte(atts_ny2, type = "dynamic", min_e = -5, max_e = 5)
summary(agg_effects_es_ny2)
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Now we are ready to go! Let me put all the outputs into a table

event_study_diff_ny <-   data.frame(
  type          = "dynamic",
  term = paste0('ATT(', agg_effects_es_ny$egt, ")"),
  event.time= agg_effects_es_ny$egt,
  estimate  = agg_effects_es_ny$att.egt,
  std.error = agg_effects_es_ny$se.egt,
  conf.low  = agg_effects_es_ny$att.egt - agg_effects_es_ny$crit.val.egt * agg_effects_es_ny$se.egt,
  conf.high = agg_effects_es_ny$att.egt + agg_effects_es_ny$crit.val.egt  * agg_effects_es_ny$se.egt,
  point.conf.low  = agg_effects_es_ny$att.egt - stats::qnorm(1 - agg_effects_es_ny$DIDparams$alp/2) * agg_effects_es_ny$se.egt,
  point.conf.high = agg_effects_es_ny$att.egt + stats::qnorm(1 - agg_effects_es_ny$DIDparams$alp/2) * agg_effects_es_ny$se.egt
)

# Now we are ready to go! Let me put all this into a table
event_study_diff_never <-   data.frame(
  type          = "dynamic",
  term = paste0('ATT(', agg_effects_es_never$egt, ")"),
  event.time= agg_effects_es_never$egt,
  estimate  = agg_effects_es_never$att.egt,
  std.error = agg_effects_es_never$se.egt,
  conf.low  = agg_effects_es_never$att.egt - agg_effects_es_never$crit.val.egt * agg_effects_es_never$se.egt,
  conf.high = agg_effects_es_never$att.egt + agg_effects_es_never$crit.val.egt  * agg_effects_es_never$se.egt,
  point.conf.low  = agg_effects_es_never$att.egt - stats::qnorm(1 - agg_effects_es_never$DIDparams$alp/2) * agg_effects_es_never$se.egt,
  point.conf.high = agg_effects_es_never$att.egt + stats::qnorm(1 - agg_effects_es_never$DIDparams$alp/2) * agg_effects_es_never$se.egt
)

event_study_diff_ny2 <-   data.frame(
  type          = "dynamic",
  term = paste0('ATT(', agg_effects_es_ny2$egt, ")"),
  event.time= agg_effects_es_ny2$egt,
  estimate  = agg_effects_es_ny2$att.egt,
  std.error = agg_effects_es_ny2$se.egt,
  conf.low  = agg_effects_es_ny2$att.egt - agg_effects_es_ny2$crit.val.egt * agg_effects_es_ny2$se.egt,
  conf.high = agg_effects_es_ny2$att.egt + agg_effects_es_ny2$crit.val.egt  * agg_effects_es_ny2$se.egt,
  point.conf.low  = agg_effects_es_ny2$att.egt - stats::qnorm(1 - agg_effects_es_ny2$DIDparams$alp/2) * agg_effects_es_ny2$se.egt,
  point.conf.high = agg_effects_es_ny2$att.egt + stats::qnorm(1 - agg_effects_es_ny2$DIDparams$alp/2) * agg_effects_es_ny2$se.egt
)


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# I will use two different styles of plot.
# You pick the one you like best!
# I usually use the first one when the number of event times is large
#  and the second one when then number event times is smallish.
#---------------------------------------------------------------------------------------
# Plots for never-treated
#---------------------------------------------------------------------------------------
# First option
p_es_never1 <- ggplot(data = event_study_diff_never,
                      mapping = aes(x = event.time, y = estimate)) +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_ribbon(aes(ymin= point.conf.low, ymax=  point.conf.high), alpha = 0.5, size = 1, fill = "steelblue")+
  geom_ribbon(aes(ymin=  conf.low, ymax =  conf.high), alpha =  0.3, size = 1, fill = "steelblue")+
  geom_line(mapping = aes(x = event.time, y=estimate), colour = "black", size = 0.6, linetype = "dashed") +
  geom_line(size = 1.2, alpha = 2, colour = "darkblue") +
  
  geom_hline(yintercept = 0, colour="black", size = 0.25, linetype = "dotted")+
  xlab('Event time') +
  ylab("Event-Study Estimate") +
  scale_x_continuous(breaks = c(-5:5)) +
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  )+
  annotate(geom="text", x=3, y=-0.01, label="Average of post-treatment ES coef's:
           0.0751 (0.013)",
           color="black")

p_es_never1


# This is another option

p_es_never2<- ggplot(data = event_study_diff_never,
                     mapping = aes(x = event.time, y = estimate)) +
  geom_line(size = 0.5, alpha = 2, colour = "black") +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = 0, colour="black",  linetype = "dotted")+
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = FALSE, linetype= 1, size = 1.1,
                  color = "red")+
  geom_pointrange(aes(ymin = point.conf.low, ymax = point.conf.high), show.legend = FALSE, size = 1.1)+
  xlab("Event time") +
  ylab("Event Study Estimate") +
  ylim(range(event_study_diff_ny$conflow,
             event_study_diff_ny$confhigh,
             event_study_diff_ny2$conflow,
             event_study_diff_ny2$confhigh,
             twfe_es$point.conf.low,
             twfe_es$point.conf.high,
             event_study_diff_never$conf.low,
             event_study_diff_never$conf.high
  ))+
  scale_x_continuous(breaks = c(-5:5)) +
  #scale_y_continuous(breaks = seq(-600, 200, 200), limits = c(-700,200))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  )+
  annotate(geom="text", x=3, y=-0.01, label="Average of post-treatment ES coef's:
           0.0751 (0.013)",
           color="black")

p_es_never2

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Plots with not-yet-treated
#---------------------------------------------------------------------------------------
# First option
p_es_ny1 <- ggplot(data = event_study_diff_ny,
                   mapping = aes(x = event.time, y = estimate)) +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_ribbon(aes(ymin= point.conf.low, ymax=  point.conf.high), alpha = 0.5, size = 1, fill = "steelblue")+
  geom_ribbon(aes(ymin=  conf.low, ymax =  conf.high), alpha =  0.3, size = 1, fill = "steelblue")+
  geom_line(mapping = aes(x = event.time, y=estimate), colour = "black", size = 0.6, linetype = "dashed") +
  geom_line(size = 1.2, alpha = 2, colour = "darkblue") +
  geom_hline(yintercept = 0, colour="black", size = 0.25, linetype = "dotted")+
  xlab('Event time') +
  ylab("Event-Study Estimate") +
  scale_x_continuous(breaks = c(-5:5)) +
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  )+
  annotate(geom="text", x=3, y=-0.01, label="Average of post-treatment ES coef's:
           0.0758 (0.0124)",
           color="black")

p_es_ny1

# This is another option

p_es_ny2<- ggplot(data = event_study_diff_ny,
                  mapping = aes(x = event.time, y = estimate)) +
  geom_line(size = 0.5, alpha = 2, colour = "black") +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = 0, colour="black",  linetype = "dotted")+
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = FALSE, linetype= 1, size = 1.1,
                  color = "red")+
  geom_pointrange(aes(ymin = point.conf.low, ymax = point.conf.high), show.legend = FALSE, size = 1.1)+
  xlab("Event time") +
  ylab("Event Study Estimate") +
  scale_x_continuous(breaks = c(-5:5)) +
  ylim(range(event_study_diff_ny$conflow,
             event_study_diff_ny$confhigh,
             event_study_diff_ny2$conflow,
             event_study_diff_ny2$confhigh,
             twfe_es$point.conf.low,
             twfe_es$point.conf.high,
             event_study_diff_never$conf.low,
             event_study_diff_never$conf.high
  ))+
  #scale_y_continuous(breaks = seq(-600, 200, 200), limits = c(-700,200))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  )+
  annotate(geom="text", x=3, y=-0.01, label="Average of post-treatment ES coef's:
          0.0758 (0.0124)",
           color="black")

p_es_ny2

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Plots for not-yet-treated, dropping never treated
#---------------------------------------------------------------------------------------
# This is one options
p_es_ny12 <- ggplot(data = event_study_diff_ny2,
                    mapping = aes(x = event.time, y = estimate)) +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_ribbon(aes(ymin= point.conf.low, ymax=  point.conf.high), alpha = 0.5, size = 1, fill = "steelblue")+
  geom_ribbon(aes(ymin=  conf.low, ymax =  conf.high), alpha =  0.3, size = 1, fill = "steelblue")+
  geom_line(mapping = aes(x = event.time, y=estimate), colour = "black", size = 0.6, linetype = "dashed") +
  geom_line(size = 1.2, alpha = 2, colour = "darkblue") +
  geom_hline(yintercept = 0, colour="black", size = 0.25, linetype = "dotted")+
  xlab('Event time') +
  ylab("Event-Study Estimate") +
  scale_x_continuous(breaks = c(-5:5)) +
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  )+
  annotate(geom="text", x=3, y=-0.01, label="Average of post-treatment ES coef's:
          0.0803 (0.0123)",
           color="black")

p_es_ny12

# This is another one that I use when event times is not that large

p_es_ny22<- ggplot(data = event_study_diff_ny2,
                   mapping = aes(x = event.time, y = estimate)) +
  geom_line(size = 0.5, alpha = 2, colour = "black") +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = 0, colour="black",  linetype = "dotted")+
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = FALSE, linetype= 1, size = 1.1,
                  color = "red")+
  geom_pointrange(aes(ymin = point.conf.low, ymax = point.conf.high), show.legend = FALSE, size = 1.1)+
  xlab("Event time") +
  ylab("Event Study Estimate") +
  scale_x_continuous(breaks = c(-5:5)) +
  ylim(range(event_study_diff_ny$conflow,
             event_study_diff_ny$confhigh,
             event_study_diff_ny2$conflow,
             event_study_diff_ny2$confhigh,
             twfe_es$point.conf.low,
             twfe_es$point.conf.high,
             event_study_diff_never$conf.low,
             event_study_diff_never$conf.high
  ))+
  #scale_y_continuous(breaks = seq(-600, 200, 200), limits = c(-700,200))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  )+
  annotate(geom="text", x=3, y=-0.01, label="Average of post-treatment ES coef's:
          0.0803 (0.0123)",
           color="black")

p_es_ny22


#---------------------------------------------------------------------------------------
# Plots for TWFE
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
p_es_twfe<- ggplot(data = twfe_es,
                   mapping = aes(x = event.time, y = estimate)) +
  geom_line(size = 0.5, alpha = 2, colour = "black") +
  geom_vline(xintercept = 0-0.05, color = 'grey', size = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = 0, colour="black",  linetype = "dotted")+
  #geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = FALSE, linetype= 1, size = 1.1,
  #                color = "red")+
  geom_pointrange(aes(ymin = point.conf.low, ymax = point.conf.high), show.legend = FALSE, size = 1.1)+
  xlab("Event time") +
  ylab("Event Study Estimate") +
  ylim(range(event_study_diff_ny$conflow,
             event_study_diff_ny$confhigh,
             event_study_diff_ny2$conflow,
             event_study_diff_ny2$confhigh,
             twfe_es$point.conf.low,
             twfe_es$point.conf.high,
             event_study_diff_never$conf.low,
             event_study_diff_never$conf.high
  ))+
  scale_x_continuous(breaks = c(-5:5)) +
  #scale_y_continuous(breaks = seq(-600, 200, 200), limits = c(-700,200))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",  size = 12))+
  theme(plot.title=ggtext::element_markdown(size = 12,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  ) + 
  annotate(geom="text", x=3, y=-0.01, label="Static TWFE Coefficient:
           0.0740 (0.013)",
           color="black")
p_es_twfe

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Plots for the ATT(g,t) using not-yet-treated
#---------------------------------------------------------------------------------------
p_attgt <- ggdid(atts_ny, 
                 title  = "ATT(g,t)'s with not-yet-treated comparison groups", 
                 ncol = 3, ylim = range(0.15, -0.15) ) + 
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12))
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Save all plots

ggsave(here("plots", "es-not-yet-treated.pdf"),
       plot = p_es_ny2, width = 12, height = 5,  bg = "transparent")


ggsave(here("plots", "es-never-treated.pdf"),
       plot = p_es_never2, width = 12, height = 5,  bg = "transparent")

ggsave(here("plots", "es-no-untreated.pdf"),
       plot = p_es_ny22, width = 12, height = 5,  bg = "transparent")

ggsave(here("plots", "es-twfe.pdf"),
       plot = p_es_twfe, width = 12, height = 5,  bg = "transparent")

ggsave(here("plots", "attgt.pdf"),
       plot = p_attgt, width = 20, height = 9,  bg = "transparent")

#---------------------------------------------------------------------------------------