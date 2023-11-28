#######################################
# interaction-networks.R
#
# Visualising the interaction results
# as a network
# 
# Tom Smith 2023
#######################################

library(ggnetwork)
library(dplyr)

# read the data
emergent_interactions <- read.csv("results/bootstrapped-emergent-interactions.csv")
net_interactions <- read.csv("results/bootstrapped-net-interactions.csv")

strain_matching <- data.frame(Strain = c("100", "302", "306", "331", "371", "419", "448", "487", "527",
                                         "74", "Ecoli", "mixture", "vfischeri"),
                              Strain_nice = c("P. baetica 1", "F. glaciei", "A. humicola 1", "P. baetica 2", "R. herbae",
                                              "S. faeni", "C. gallinarum", "A. popoffii", "A. humicola 2", "N. soli", 
                                              "E. coli", "mixture", "A. fischeri"))

emergent_interactions <- left_join(emergent_interactions, strain_matching)
net_interactions <- left_join(net_interactions, strain_matching)

# create a field describing the chemicals
# in mixture using their first letters
emergent_interactions$ChemCode <- NA
for(i in 1:nrow(emergent_interactions)){
  subs_dat <- emergent_interactions[i,]
    
  chem_codes <- names(subs_dat[,2:9])[apply(subs_dat[,2:9], 1, function(i) which(i == 1))]
  chem_result <- paste(gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1', chem_codes, perl = TRUE), collapse = "")
    
  emergent_interactions[i,]$ChemCode <- chem_result
}

net_interactions$ChemCode <- NA
for(i in 1:nrow(net_interactions)){
  subs_dat <- net_interactions[i,]
  
  chem_codes <- names(subs_dat[,2:9])[apply(subs_dat[,2:9], 1, function(i) which(i == 1))]
  chem_result <- paste(gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1', chem_codes, perl = TRUE), collapse = "")
  
  net_interactions[i,]$ChemCode <- chem_result
}


###########################################
# create the necessary parts for a network

chem_matrix <- read.csv("data/chemical-matrix-plus.csv")

chem_matrix <- chem_matrix %>%
  arrange(Complexity, ChemCode)

# create edges
# edges will be from the complexity level below, which contain all but one of the chems in the upper node

from_vector <- c()
to_vector <- c()

for(i in 2:8){
  from_options <- chem_matrix[chem_matrix$Complexity == i-1,]$ChemCode
  to_options <- chem_matrix[chem_matrix$Complexity == i,]$ChemCode
  for(j in 1:length(from_options)){
    for(k in 1:length(to_options)){
      condition <- unlist(strsplit(from_options[j], split = "")) %in% unlist(strsplit(to_options[k], split = ""))
      if(!(FALSE %in% condition)){
        from_vector <- c(from_vector, from_options[j])
        to_vector <- c(to_vector, to_options[k])
      }
    }
  }
  
}


all_edges <- data.frame(source = from_vector, destination = to_vector)

# create nodes which can be merged with the interaction dataset
all_nodes <- data.frame(label = unique(c(all_edges$source, all_edges$destination)))
all_nodes$id <- rownames(all_nodes)

nodes <- data.frame(id = all_nodes$id, label = all_nodes$label)

# now need to link the edges to the ID labels in nodes
edges <- all_edges %>% 
  left_join(nodes, by = c("source" = "label")) %>% 
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("destination" = "label")) %>% 
  rename(to = id)

edges <- select(edges, from, to)

# The layout is defined in a matrix with 2 columns and a row for each node. 
# The first column indicates its x position and the second its y position, and scale is not relevant 
# (it is always rescaled to fit a -1 to 1 plotting area.)

# nodes$x <- c(seq(28,42, by = 2), seq(21,48), seq(7,62), seq(1,70), seq(7,62), seq(21,48), seq(32,39), 35) 
nodes$x <- c(seq(21,49, by = 4), seq(8,62, by = 2), seq(2,68, by = 1.2), seq(1,70), seq(2,68, by = 1.2), seq(8,62, by = 2), seq(21,49, by = 4), 35) 
nodes$y <- c(rep(1, 8), rep(2, 28), rep(3, 56), rep(4, 70), rep(5, 56), rep(6, 28), rep(7, 8), 8)

# create the network from the nodes and edges
chem_network_gg <- ggnetwork(edges, layout = as.matrix(nodes[,c("x", "y")]), arrow.gap = 0)

# now going to create a network for each strain by adding
# metadata to the above network
full_network <- data.frame()

net_interactions <- net_interactions %>%
  mutate(Net_interaction = ifelse(Complexity > 1, response, NA))
emergent_interactions <- emergent_interactions %>%
  mutate(Emergent_interaction = ifelse(Complexity > 1, response, NA))

# this bit to visualise everything
interaction_data <- net_interactions

interaction_data$Emergent_interaction <- emergent_interactions$Emergent_interaction

interaction_data <- interaction_data %>% filter(Complexity > 0) %>%
  rename(chem_id = id)

strains <- unique(interaction_data$Strain)

for(i in seq_along(strains)){
  strain_nodes <- left_join(nodes, interaction_data %>% filter(Strain == strains[i]), by = c("label" = "ChemCode"))
  # merge in the metadata from before
  strain_network_gg <- left_join(chem_network_gg, strain_nodes[,c("Strain", "Strain_nice", "id", "label", "mean_response", "Net_interaction", "Emergent_interaction")],
                               by = c("vertex.names" = "id"))
  
  full_network <- bind_rows(full_network, strain_network_gg)
}




###############################################
## --- Create the network visualisations --- ##
###############################################

# Networks with only the interaction nodes
# and edges plotted

all_networks <- data.frame()

strain_list <- c("306", "331", "371", "419", "448", "487", "527", "74", "Ecoli", "vfischeri", "mixture")
#nice_strain_list <- c("306", "331", "371", "419", "448", "487", "527", "74", "E. coli", "A. fischeri")

for(strain in seq_along(strain_list)){
  strain_network <- full_network %>% filter(Strain == strain_list[strain])
  # take only edges that go between nodes with interactions, and the single-stressors
  interaction_nodes <- strain_network %>% filter(Net_interaction != "Additive" |
                                                label %in% c("A", "C", "D", "G", "I", "M", "O", "T")) # %>%
    # mutate(Strain_nice = nice_strain_list[strain])
  # now we just filter out all nodes whose xend/yend don't match an x/y in the data
  interaction_edges <- interaction_nodes[paste0(interaction_nodes$xend, interaction_nodes$yend) %in% paste0(interaction_nodes$x, interaction_nodes$y),]
  # bind into the data
  all_networks <- bind_rows(all_networks, interaction_edges)
}

# make Strain nice a factor to order them how we like
# can we roughly order it by phylogeny??
all_networks$Strain_nice <- factor(all_networks$Strain_nice, 
                                   levels = c("P. baetica 2", "E. coli", "A. fischeri", "A. popoffii", "S. faeni", 
                                              "R. herbae", "N. soli", "C. gallinarum", "A. humicola 1", "A. humicola 2",
                                              "mixture"))

all_interactions_network <- ggplot(all_networks, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodes(data = full_network %>% filter(!Strain %in% c("100", "302")), shape = 21) +
  geom_nodes(data = all_networks %>% filter(!label %in% c("A", "C", "D", "G", "I", "M", "O", "T")),
             aes(fill = Net_interaction), shape = 21, size = 3) +
  geom_nodelabel(data = full_network[full_network$vertex.names %in% c(1:8),] %>% filter(!Strain %in% c("100", "302")),
                 aes(label = label), nudge_y = -0.05, fontface = "bold") +
  # geom_nodetext(data = full_network[full_network$vertex.names %in% c(1:8),] %>% filter(!Strain %in% c("100", "302", "mixture")), 
  #               aes(label = label), fontface = "bold") +
  # geom_nodelabel_repel(aes(label = label), fontface = "bold", box.padding = unit(1, "lines")) +
  theme_blank() +
  scale_fill_manual(values = c("#21918c", "#fde725"), name = "Net Interaction") +
  guides(fill = guide_legend(override.aes = list(size=6))) +
  facet_wrap(~as.factor(Strain_nice)) +
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.position = c(0.9, 0.2),
        legend.direction = "vertical", legend.box = "horizontal")
all_interactions_network

ggsave("results/all_interaction_networks.svg", all_interactions_network, width = 16, height = 12)


# repeat this to generate two plots which we can layout next to each other in inkscape

all_interactions_network_top <- ggplot(all_networks %>% filter(Strain_nice %in% c("P. baetica 2", "E. coli", "A. fischeri", "A. popoffii",
                                                                                  "S. faeni", "R. herbae")), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodes(data = full_network %>% filter(Strain_nice %in% c("P. baetica 2", "E. coli", "A. fischeri", "A. popoffii",
                                                               "S. faeni", "R. herbae")), shape = 21) +
  geom_nodes(data = all_networks %>% filter(!label %in% c("A", "C", "D", "G", "I", "M", "O", "T"),
                                            Strain_nice %in% c("P. baetica 2", "E. coli", "A. fischeri", "A. popoffii",
                                                               "S. faeni", "R. herbae")),
             aes(fill = Net_interaction), shape = 21, size = 3) +
  geom_nodelabel(data = full_network[full_network$vertex.names %in% c(1:8),] %>% filter(Strain_nice %in% c("P. baetica 2", "E. coli", "A. fischeri", "A. popoffii",
                                                                                                           "S. faeni", "R. herbae")),
                 aes(label = label), nudge_y = -0.05, fontface = "bold") +
  # geom_nodetext(data = full_network[full_network$vertex.names %in% c(1:8),] %>% filter(!Strain %in% c("100", "302", "mixture")), 
  #               aes(label = label), fontface = "bold") +
  # geom_nodelabel_repel(aes(label = label), fontface = "bold", box.padding = unit(1, "lines")) +
  theme_blank() +
  scale_fill_manual(values = c("#35b779", "#fde725"), name = "Net Interaction") +
  guides(fill = guide_legend(override.aes = list(size=6))) +
  facet_wrap(~as.factor(Strain_nice), ncol = 4) +
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.position = c(0.7, 0.25),
        legend.direction = "vertical", legend.box = "horizontal")
all_interactions_network_top

ggsave("results/all_interaction_networks_top.svg", all_interactions_network_top, width = 16, height = 8)

all_interactions_network_bot <- ggplot(all_networks %>% filter(Strain_nice %in% c("N. soli", "C. gallinarum", "A. humicola 1", "A. humicola 2")), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodes(data = full_network %>% filter(Strain_nice %in% c("N. soli", "C. gallinarum", "A. humicola 1", "A. humicola 2")), shape = 21) +
  geom_nodes(data = all_networks %>% filter(!label %in% c("A", "C", "D", "G", "I", "M", "O", "T"),
                                            Strain_nice %in% c("N. soli", "C. gallinarum", "A. humicola 1", "A. humicola 2")),
             aes(fill = Net_interaction), shape = 21, size = 3) +
  geom_nodelabel(data = full_network[full_network$vertex.names %in% c(1:8),] %>% filter(Strain_nice %in% c("N. soli", "C. gallinarum", "A. humicola 1", "A. humicola 2")),
                 aes(label = label), nudge_y = -0.05, fontface = "bold") +
  # geom_nodetext(data = full_network[full_network$vertex.names %in% c(1:8),] %>% filter(!Strain %in% c("100", "302", "mixture")), 
  #               aes(label = label), fontface = "bold") +
  # geom_nodelabel_repel(aes(label = label), fontface = "bold", box.padding = unit(1, "lines")) +
  theme_blank() +
  scale_fill_manual(values = c("#35b779", "#fde725"), name = "Net Interaction") +
  guides(fill = guide_legend(override.aes = list(size=6))) +
  facet_wrap(~as.factor(Strain_nice), ncol = 4) +
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.position = "none")
all_interactions_network_bot

ggsave("results/all_interaction_networks_bot.svg", all_interactions_network_bot, width = 16, height = 4)

#################################################
## Oxytetracycline + Tebuconazole interactions ##
## i.e. supplementary figure S6                ##
#################################################
#
# that network plot is cool and all... 
# but what if we use it to show something more specifically
# that this abundant interaction works its way up through
# the levels of complexity. 

# how to do it...
# very light shade of grey for whole network
# subset to just those strains with the O/T interaction
# then colour the O/T edges/nodes by
# what interaction they represent

# nb strains = A. fischeri; A. popoffii; S. faeni; R. herbae; C. gallinarum; A. humicola 2.

OT_networks <- data.frame()

strain_list <- c("vfischeri", "487", "371", "419", "448", "527")

# first subset it too all the O/T interactions
OT_labels <- unique(net_interactions[net_interactions$Oxytetracycline == 1 & net_interactions$Tebuconazole == 1,]$ChemCode)

for(strain in seq_along(strain_list)){
  strain_network <- full_network %>% filter(Strain == strain_list[strain])
  # take only edges that go between nodes with interactions, and the single-stressors
  interaction_nodes <- strain_network %>% filter(label %in% OT_labels |
                                                   label %in% c("A", "C", "D", "G", "I", "M", "O", "T")) # %>%
  # mutate(Strain_nice = nice_strain_list[strain])
  # now we just filter out all nodes whose xend/yend don't match an x/y in the data
  interaction_edges <- interaction_nodes[paste0(interaction_nodes$xend, interaction_nodes$yend) %in% paste0(interaction_nodes$x, interaction_nodes$y),]
  # bind into the data
  OT_networks <- bind_rows(OT_networks, interaction_edges)
}

# now subset it specifically to the antagonisms or synergisms
antagonism_strains <- c("vfischeri", "487", "448", "527")
antagonism_networks <- data.frame()
for(strain in seq_along(antagonism_strains)){
  strain_network <- OT_networks %>% filter(Strain == antagonism_strains[strain])
  # take only edges that go between nodes with interactions, and the single-stressors
  interaction_nodes <- strain_network %>% filter(Net_interaction == "Antagonism" |
                                                   label %in% c("A", "C", "D", "G", "I", "M", "O", "T")) # %>%
    # mutate(Strain_nice = nice_strain_list[strain])
  # now we just filter out all nodes whose xend/yend don't match an x/y in the data
  interaction_edges <- interaction_nodes[paste0(interaction_nodes$xend, interaction_nodes$yend) %in% paste0(interaction_nodes$x, interaction_nodes$y),]
  # bind into the data
  antagonism_networks <- bind_rows(antagonism_networks, interaction_edges)
}

synergism_strains <- c("371", "419")
synergism_networks <- data.frame()
for(strain in seq_along(synergism_strains)){
  strain_network <- OT_networks %>% filter(Strain == synergism_strains[strain])
  # take only edges that go between nodes with interactions, and the single-stressors
  interaction_nodes <- strain_network %>% filter(Net_interaction == "Synergism" |
                                                   label %in% c("A", "C", "D", "G", "I", "M", "O", "T")) # %>%
  # mutate(Strain_nice = nice_strain_list[strain])
  # now we just filter out all nodes whose xend/yend don't match an x/y in the data
  interaction_edges <- interaction_nodes[paste0(interaction_nodes$xend, interaction_nodes$yend) %in% paste0(interaction_nodes$x, interaction_nodes$y),]
  # bind into the data
  synergism_networks <- bind_rows(synergism_networks, interaction_edges)
}

# make the plot
OT_interaction_network <-  ggplot(full_network %>% filter(Strain %in% c(antagonism_strains, synergism_strains)), aes(x = x, y = y, xend = xend, yend = yend)) +
  #geom_edges(color = "grey95") +
  geom_edges(data = OT_networks, color = "grey80") +
  geom_edges(data = antagonism_networks, color = "grey20") +
  geom_edges(data = synergism_networks, color = "grey20") +
  geom_nodes(data = full_network %>% filter(Strain %in% c(antagonism_strains, synergism_strains)), shape = 21) +
  geom_nodes(data = antagonism_networks %>% filter(!label %in% c("A", "C", "D", "G", "I", "M", "O", "T")),
             fill = "#35b779", shape = 21, size = 3) +
  geom_nodes(data = synergism_networks %>% filter(!label %in% c("A", "C", "D", "G", "I", "M", "O", "T")),
             fill = "#fde725", shape = 21, size = 3) +
  geom_nodelabel(data = OT_networks[OT_networks$vertex.names %in% c(1:8),],
                 aes(label = label), nudge_y = -0.05, fontface = "bold") +
  theme_blank() +
  facet_wrap(~as.factor(Strain_nice)) +
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.position = c(0.9, 0.2),
        legend.direction = "vertical", legend.box = "horizontal")
OT_interaction_network


svg("results/supplementary-fig-s6.svg", width = 16, height = 8)
OT_interaction_network
dev.off()
