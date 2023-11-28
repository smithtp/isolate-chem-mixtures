#########################################
# basic-visualisations.R
#
# Create various visualisations of our
# results for figs:
# 1b, 2a,b,c, 4b,c,e,f
# supplementary figs:
# S1, S2, S5a,b, S7
#
# Tom Smith 2023
#########################################

library(ggplot2)
library(grafify)
library(tidyr)
library(dplyr)
library(forcats)
library(ggbeeswarm) # to avoid overplotting

# set main plotting theme
plotting_theme <- theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 16))

# read in the data
emergent_interactions <- read.csv("results/bootstrapped-emergent-interactions.csv")
net_interactions <- read.csv("results/bootstrapped-net-interactions.csv")

# add a new set of strain names
strain_matching <- data.frame(Strain = c("100", "302", "306", "331", "371", "419", "448", "487", "527",
                                         "74", "Ecoli", "mixture", "vfischeri"),
                              Strain_nice = c("P. baetica 1", "F. glaciei", "A. humicola 1", "P. baetica 2", "R. herbae",
                                              "S. faeni", "C. gallinarum", "A. popoffii", "A. humicola 2", "N. soli", 
                                              "E. coli", "mixture", "A. fischeri"))

emergent_interactions <- left_join(emergent_interactions, strain_matching)
net_interactions <- left_join(net_interactions, strain_matching)

# also consider whether the responses are positive or negative
additive_effects <- read.csv("results/bootstrapped-no-response-test.csv")


########################
# ---    Fig 1B    --- #
########################
# Side plot for the heatmap
# Interactions vs additive effects
additive_summary <- additive_effects %>%
  group_by(Complexity, response) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# choose what levels we want for plotting
additive_summary$response <- factor(additive_summary$response, levels = c("None", "Negative", "Positive"))

additive_summary$Complexity <- factor(additive_summary$Complexity, levels = c(8, 7, 6, 5, 4, 3, 2, 1, 0))

sideplot <- ggplot(additive_summary %>% filter(Complexity != 0), aes(x = response, y = freq)) +
  geom_col(aes(fill = response), col = "black") +
  scale_fill_manual(values = c("white", "firebrick", "slateblue")) +
  labs(x = "Response", y = "Frequency", title = "Complexity") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  facet_wrap(~Complexity, ncol = 1) +
  plotting_theme +
  theme(legend.position = "none",
        #axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        strip.text = element_text(size = 16))
sideplot  

ggsave("results/response-side-plot.svg", sideplot, width = 3, height = 12)

# quantify it
#
negative.lm <- lm(freq ~ as.numeric(as.character(Complexity)), data = additive_summary %>% filter(response == "Negative"))
summary(negative.lm)

########################
# ---    Fig 2A    --- #
########################

# do a little model selection
linmod <- lm(mean_response ~ Complexity, data = net_interactions %>% filter(Strain != "mixture"  & Complexity > 0))
logmod <- lm(mean_response ~ log(Complexity), data = net_interactions %>% filter(Strain != "mixture"  & Complexity > 0))
polymod <- lm(mean_response ~ poly(Complexity, 2), data = net_interactions %>% filter(Strain != "mixture"  & Complexity > 0))

summary(linmod)
summary(logmod)
summary(polymod)

AIC(linmod)
AIC(logmod)
AIC(polymod)
# the straight-line model has lower AIC

# models with oxytet as a factor
linmod_tet <- lm(mean_response ~ Complexity + Oxytetracycline, data = net_interactions %>% filter(Strain != "mixture"  & Complexity > 0))
summary(linmod_tet)

# create the plot
chems_nooxy <- seq(1, 7, 1)
dauc_nooxy <- coefficients(linmod_tet)[["(Intercept)"]] + coefficients(linmod_tet)[["Complexity"]] * chems_nooxy

chems_oxy <- seq(1, 8, 1)
dauc_oxy <- coefficients(linmod_tet)[["(Intercept)"]] + coefficients(linmod_tet)[["Complexity"]] * chems_oxy +
  coefficients(linmod_tet)[["Oxytetracycline"]]

line_nooxy <- data.frame(chems_nooxy, dauc_nooxy)
line_oxy <- data.frame(chems_oxy, dauc_oxy)

p_beeswarm_tet <- ggplot(net_interactions %>% filter(!(Strain %in% c("mixture", "Ecoli", "vfischeri"))
                                                     & Complexity > 0), aes(x = Complexity, y = mean_response)) +
  geom_beeswarm(aes(fill = as.character(Oxytetracycline)), shape = 21, cex = 0.5, alpha = 0.6) +
  labs(x = "Number of chemicals", y = "Relative Growth (G)", fill = "Oxytetracycline") +
  scale_fill_manual(values = c("black", "orange")) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  lims(y = c(0, 1.5)) +
  geom_line(data = line_nooxy, aes(x = chems_nooxy, y = dauc_nooxy), lwd = 1) +
  geom_line(data = line_oxy, aes(x = chems_oxy, y = dauc_oxy), col = "orange", lwd = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  plotting_theme +
  theme(legend.position = "bottom")
p_beeswarm_tet

ggsave("results/fig-2a.svg", p_beeswarm_tet, width = 7, height = 7)


########################
# ---    Fig 2B    --- #
########################

# How similar is mixture to the mean of other IL strains?
mixture_summary <- net_interactions %>%
  filter(Strain == "mixture") %>%
  select(id, r_mean, mean_response, Complexity, Amoxicillin, Chlorothalonil, Diflufenican, Glyphosate,
         Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole, response) %>%
  rename(Mixture.AUC = r_mean, Mixture.proportion = mean_response)

iceland_average <- net_interactions %>%
  filter(Strain %in% c("74", "100", "302", "306", "331", "371", "419", "448", "487", "527")) %>%
  select(id, r_mean, mean_response) %>%
  group_by(id) %>%
  summarise(Mean.AUC = mean(r_mean), Mean.proportion = mean(mean_response))

comparison_df <- left_join(mixture_summary, iceland_average, by = "id")

# repeat a beeswarm plot for the mixture only?

mixture_beeswarm <- ggplot(mixture_summary %>% filter(Complexity > 0), aes(x = Complexity, y = Mixture.proportion)) +
  geom_beeswarm(aes(fill = as.character(Oxytetracycline)), shape = 21, cex = 0.5, alpha = 0.6, size = 2) +
  labs(x = "Number of chemicals", y = "Relative Growth (G)", fill = "Oxytetracycline") +
  scale_fill_manual(values = c("black", "orange")) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  lims(y = c(0, 1.5)) +
  #geom_line(data = line_nooxy, aes(x = chems_nooxy, y = dauc_nooxy), lwd = 1) +
  #geom_line(data = line_oxy, aes(x = chems_oxy, y = dauc_oxy), col = "orange", lwd = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  plotting_theme +
  theme(legend.position = "bottom")
mixture_beeswarm

ggsave("results/fig-2b.svg", mixture_beeswarm, width = 7, height= 7)

# get the regression lines
linmod_tet_mix <- lm(Mixture.proportion ~ Complexity + Oxytetracycline, data = mixture_summary)
summary(linmod_tet_mix)


########################
# ---    Fig 2C    --- #
########################

mixture_vs_isolates <- ggplot(comparison_df, aes(x = Mean.proportion, y = Mixture.proportion)) +
  geom_point(size = 5, alpha = 0.8, aes(shape = as.factor(Oxytetracycline), fill = Complexity)) +
  scale_fill_viridis_c(name = "# Chemicals") +
  scale_shape_manual(values = c(21, 22), name = "Oxytetracycline") +
  geom_abline(slope = 1, linetype = "dashed") +
  labs(x = "Relative Growth (Mean of Isolates)", y = "Relative Growth (Isolate Mixture)") +
  #guides(shape = guide_legend(title = "Oxytetracycline"),) +
  plotting_theme +
  theme(legend.position = c(0.25, 0.7))
mixture_vs_isolates

ggsave("results/fig-2c.svg", mixture_vs_isolates, width = 7, height = 7)


#####################################
# ---   Supplementary Fig S2    --- #
#####################################

p_linear_strains <- ggplot(net_interactions %>% filter(Strain != "mixture"), aes(x = Complexity, y = mean_response)) +
  geom_point(aes(col = Strain_nice)) +
  scale_colour_grafify(palette = "muted") +
  labs(x = "Number of chemicals", y = "Relative Growth (G)", title = "All Monoculture Strains", col = "Strain") +
  geom_smooth(method = lm, col = "black") +
  theme_bw() +
  facet_wrap(~Strain_nice)
p_linear_strains

svg("results/supp-fig-s2.svg", width = 9, height = 6)
p_linear_strains
dev.off()

### --------------------------------------------------------------------- ###

##################################
## --- INTERACTIONS FIGURES --- ##
##################################

# need to incorporate no response into the interactions
# add the test of responses being additive, or just no response
noresponse_test <- additive_effects %>%
  select(Strain, id, significant) %>%
  rename(additive.significant = significant)

net_interactions <- left_join(net_interactions, noresponse_test)

net_interactions$response.tested <- net_interactions$response
net_interactions[net_interactions$response == "Additive" & 
                   net_interactions$additive.significant == FALSE,]$response.tested <- "None"

# same for emergent
emergent_interactions <- left_join(emergent_interactions, noresponse_test)

emergent_interactions$response.tested <- emergent_interactions$response
emergent_interactions[emergent_interactions$response == "Additive" & 
                   emergent_interactions$additive.significant == FALSE,]$response.tested <- "None"

net_summary <- net_interactions %>%
  filter(Strain != "mixture") %>%
  group_by(Complexity, response.tested) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  rename(response = response.tested)

emergent_summary <- emergent_interactions %>%
  filter(Strain != "mixture") %>%
  group_by(Complexity, response.tested) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  rename(response = response.tested)


#####################################
# ---   Supplementary Fig S5    --- #
#####################################

# like a STRUCTURE plot, showing antagonisms/synergism % horizontally across 255 regions
# get the factor ordered for plotting
net_interactions_multi <- net_interactions %>% filter(Complexity > 1)
net_interactions_multi$response.tested <- factor(net_interactions_multi$response.tested, 
                                                 levels = c("None", "Additive", "Antagonism", "Synergism"))
emergent_interactions_multi <- emergent_interactions %>% filter(Complexity > 1)
emergent_interactions_multi$response.tested <- factor(emergent_interactions_multi$response.tested, 
                                                 levels = c("None", "Additive", "Antagonism", "Synergism"))

# make an empty theme removing a load of the standard stuff
facet_plot_theme <- theme(axis.text = element_blank(),
                          axis.title = element_blank(),
                          axis.ticks = element_blank(),
                          panel.background = element_blank(),
                          legend.title = element_blank(),
                          strip.background = element_blank(),
                          strip.text = element_blank(),
                          legend.text = element_blank())

net_interactions_multi$id <- factor(net_interactions_multi$id, levels = net_interactions_multi$id[1:247])

all_net_interactions <- ggplot(net_interactions_multi %>% filter(Strain != "mixture"), aes(x = id, fill = response.tested)) + 
  geom_bar() +
  scale_fill_viridis_d() +
  facet_grid(~Complexity, scales = "free_x", space = "free") +
  facet_plot_theme +
  theme(legend.position = "bottom",
        legend.spacing.x = unit(2.0, 'cm'))
all_net_interactions

emergent_interactions_multi$id <- factor(emergent_interactions_multi$id, levels = emergent_interactions_multi$id[1:247])

all_emergent_interactions <- ggplot(emergent_interactions_multi %>% filter(Strain != "mixture"), aes(x = id, fill = response.tested)) + 
  geom_bar() +
  scale_fill_viridis_d(name = "Interaction") +
  facet_grid(~Complexity, scales = "free_x", space = "free") +
  facet_plot_theme +
  theme(legend.position = "none")
all_emergent_interactions

ggsave("results/supp-fig-S5a.svg", all_net_interactions, width = 12, height = 3)
ggsave("results/supp-fig-S5b.svg", all_emergent_interactions, width = 12, height = 3)


#########################
# ---   Figure 4    --- #
#########################

# summarise the above into a more digestible figure
summary_net_interactions <- ggplot(net_interactions_multi  %>% filter(Strain != "mixture"), aes(x = Complexity, fill = response.tested)) + 
  geom_bar(position = "fill") +
  scale_fill_viridis_d() +
  #labs(x = "Number of Chemicals", y = "Response Proportion") +
  scale_x_continuous(breaks=seq(2,8,1)) +
  # facet_plot_theme +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank(),
        legend.position = "none")
summary_net_interactions

summary_emergent_interactions <- ggplot(emergent_interactions_multi  %>% filter(Strain != "mixture"), aes(x = Complexity, fill = response.tested)) + 
  geom_bar(position = "fill") +
  scale_fill_viridis_d() +
  labs(x = "Number of Chemicals", y = "Response Proportion") +
  scale_x_continuous(breaks=seq(2,8,1)) +
  # facet_plot_theme +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank(),
        legend.position = "none")
summary_emergent_interactions

ggsave("results/fig-4b.svg", summary_net_interactions, width = 6, height = 4)
ggsave("results/fig-4e.svg", summary_emergent_interactions, width = 6, height = 4)

# contrast these with the mixture responses
mixture_net_interactions <- ggplot(net_interactions_multi  %>% filter(Strain == "mixture"), aes(x = Complexity, fill = response.tested)) + 
  geom_bar(position = "fill") +
  scale_fill_viridis_d() +
  #labs(x = "Number of Chemicals", y = "Response Proportion") +
  scale_x_continuous(breaks=seq(2,8,1)) +
  # facet_plot_theme +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank(),
        legend.position = "none")
mixture_net_interactions

mixture_emergent_interactions <- ggplot(emergent_interactions_multi  %>% filter(Strain == "mixture"), aes(x = Complexity, fill = response.tested)) + 
  geom_bar(position = "fill") +
  scale_fill_viridis_d() +
  labs(x = "Number of Chemicals", y = "Response Proportion") +
  scale_x_continuous(breaks=seq(2,8,1)) +
  # facet_plot_theme +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank(),
        legend.position = "none")
mixture_emergent_interactions

ggsave("results/fig-4c.svg", mixture_net_interactions, width = 6, height = 4)
ggsave("results/fig-4f.svg", mixture_emergent_interactions, width = 6, height = 4)




####################################
# --- Supplementarry Figure S7 --- #
####################################
# try to produce a supplementary plot
# that shows something useful about which
# interactions are prevalent - which
# mixture combinations commonly produce
# an interaction?

net_interactions_multi %>% 
  filter(Strain != "mixture") %>%
  group_by(id, response.tested) %>%
  summarise(n = n()) %>%
  #filter(n < 11) %>%
  pivot_wider(names_from = response.tested, values_from = n) %>% 
  replace(is.na(.), 0) %>% # replace NAs with 0
  mutate(total.interactions = Antagonism + Synergism) %>%
  filter(total.interactions >= 6)

emergent_interactions_multi %>% 
  filter(Strain != "mixture") %>%
  group_by(id, response.tested) %>%
  summarise(n = n()) %>%
  #filter(n < 11) %>%
  pivot_wider(names_from = response.tested, values_from = n) %>% 
  replace(is.na(.), 0) %>%
  mutate(total.interactions = Antagonism + Synergism) %>%
  filter(total.interactions >= 4)

# create a heatmap plot from this.
interesting.interactions <- emergent_interactions_multi %>% 
  filter(Strain != "mixture") %>%
  group_by(id, response.tested) %>%
  summarise(n = n()) %>%
  #filter(n < 11) %>%
  pivot_wider(names_from = response.tested, values_from = n) %>% 
  replace(is.na(.), 0) %>%
  mutate(total.interactions = Antagonism + Synergism) %>%
  filter(total.interactions >= 6) %>%
  rename(Multiplicative = Additive) %>%
  select(-total.interactions) %>%
  tibble::column_to_rownames("id")

# turn it into proportions
interesting.interactions <- interesting.interactions/12

hmp <- Heatmap(as.matrix(interesting.interactions),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = c("white", "black"),
        name = "Prevalence",
        #show_row_names = FALSE,
        show_column_names = FALSE)

svg("results/supplementary-fig-s7.svg", 5, 4)
hmp
dev.off()


####################################
# --- Supplementarry Figure S1 --- #
####################################

# response of each strain to the single stressors only
library(readr)

# load data with all replicates
all_data <- read_csv("results/spline_fits.csv", show_col_types = FALSE) %>%
  arrange(Strain, Complexity)

# single stressors
single_stressor_data <- all_data %>% 
  filter(Complexity == 1)

control_data <- all_data %>% 
  filter(Complexity == 0)

control_means <- control_data %>%
  select(Strain, AUC) %>%
  group_by(Strain) %>%
  summarise(control_AUC = mean(AUC))

single_stressor_data <- left_join(single_stressor_data, control_means) %>%
  mutate(dAUC = AUC/control_AUC) %>%
  filter(Strain != "mixture")

strain_matching <- data.frame(Strain = c("100", "302", "306", "331", "371", "419", "448", "487", "527",
                                         "74", "Ecoli", "mixture", "vfischeri"),
                              Strain_nicer = c("P. baetica 1", "F. glaciei", "A. humicola 1", "P. baetica 2", "R. herbae",
                                               "S. faeni", "C. gallinarum", "A. popoffii", "A. humicola 2", "N. soli", 
                                               "E. coli", "mixture", "A. fischeri"))

single_stressor_data <- left_join(single_stressor_data, strain_matching) %>%
  mutate(Stressor = NA)

single_stressor_data[single_stressor_data$Amoxicillin == 1,]$Stressor <- "Amoxicillin"
single_stressor_data[single_stressor_data$Chlorothalonil == 1,]$Stressor <- "Chlorothalonil"
single_stressor_data[single_stressor_data$Diflufenican == 1,]$Stressor <- "Diflufenican"
single_stressor_data[single_stressor_data$Glyphosate == 1,]$Stressor <- "Glyphosate"
single_stressor_data[single_stressor_data$Imidacloprid == 1,]$Stressor <- "Imidacloprid"
single_stressor_data[single_stressor_data$Metaldehyde == 1,]$Stressor <- "Metaldehyde"
single_stressor_data[single_stressor_data$Oxytetracycline == 1,]$Stressor <- "Oxytetracycline"
single_stressor_data[single_stressor_data$Tebuconazole == 1,]$Stressor <- "Tebuconazole"


stressor_plot <- ggplot(single_stressor_data, aes(x = Strain_nicer, y = dAUC, fill = Strain_nicer)) +
  geom_boxplot() +
  scale_fill_grafify(palette = "muted") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(x = "Strain", y = "Relative Growth (G)") +
  facet_wrap(~Stressor, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
stressor_plot

ggsave("results/supplementary-fig-s1.svg", stressor_plot)