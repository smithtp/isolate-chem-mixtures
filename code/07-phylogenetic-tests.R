#######################################
# phylogenetic-tests.R
#
# Script to perform phylogenetic analyses
# and produce phylogenetic visualisations
# 
# Tom Smith 2023
#######################################

library(dplyr)
library(tidyr)
library(phytools)
library(vegan)
library(ggplot2)

rm(list = ls())

## read in the tree
tree <- read.tree("data/phylogenetics/RAxML_result.constrained_output")
## plot tree
plotTree(tree)
nodelabels()

# reroot
rerooted_tree <- reroot(tree, 16)
# make it essentially a cladogram:
rerooted_tree$edge.length <- NULL

plot(rerooted_tree)

clustered_tree <- read.tree("results/species_clustering.nwk")

plotTree(clustered_tree)

# neaten up the tip labels
rerooted_tree$tip.label <- c("C. gallinarum", "A. humicola 2", "A. humicola 1", "N. soli", "P. baetica 2",
                             "P. baetica 1", "A. popoffii", "A. fischeri", "E. coli", "R. herbae", "S. faeni",
                             "F. glaciei")
clustered_tree$tip.label <- c("P. baetica 1", "F. glaciei", "P. baetica 2", "E. coli", "N. soli", "A. fischeri", 
                              "S. faeni", "A. popoffii", "C. gallinarum", "R. herbae", "A. humicola 1",
                              "A. humicola 2")

# plot them side-by-side
plot(rerooted_tree)
nodelabels()
rerooted_tree <- rotate(rerooted_tree, node = 21)
plot(rerooted_tree)

# make the clustered tree also a cladogram
clustered_tree$edge.length <- NULL

plot(clustered_tree)
nodelabels()
clustered_tree <- rotate(clustered_tree, node = 19)
clustered_tree <- rotate(clustered_tree, node = 22)

svg("results/tree-topology-phylo.svg")
plotTree(rerooted_tree)
dev.off()

svg("results/tree-topology-pheno.svg")
plotTree(clustered_tree, direction = "leftwards")
dev.off()

# now without labels
svg("results/tree-topology-phylo-nolabels.svg")
plot(rerooted_tree, show.tip.label = FALSE, edge.width = 3)
dev.off()

svg("results/tree-topology-pheno-nolabels.svg")
plot(clustered_tree, direction = "leftwards", show.tip.label = FALSE, edge.width = 3)
dev.off()


#########################################################
# Phylogenetic test
# Is there phylogenetic signal in the responses
# of the bacteria to chems and chem mixes?
#
# can we just directly compare the distance matrices
# from the phylogenetic tree and the phenotypic clustering?
#
# ----- MANTEL TEST ----- #

# give the tree some sensible tip labels for this:
phylo_tree <- reroot(tree, 16)
phylo_tree$tip.label <- c("448 C. gallinarum", "527 A. humicola", "306 A. humicola", "74 N. soli", "331 P. baetica",
                          "100 P. baetica", "487 A. popoffii", "A. fischeri", "E. coli", "371 R. herbae", "419 S. faeni",
                          "302 F. glaciei")

phylo_dist <- cophenetic.phylo(phylo_tree)
# clustered_dist <- cophenetic.phylo(clustered_tree)

# read in the data to compute phenotypic distances
net_interactions <- read.csv("results/bootstrapped-net-interactions.csv")
emergent_interactions <- read.csv("results/bootstrapped-emergent-interactions.csv")

net_interactions[net_interactions$Strain == "Ecoli",]$Strain <- "E. coli"
net_interactions[net_interactions$Strain == "vfischeri",]$Strain <- "A. fischeri"

interaction.df <- net_interactions %>%
  filter(Complexity > 0) %>%
  filter(Strain != "mixture") %>%
  select(id, mean_response, Strain) %>%
  pivot_wider(names_from = Strain, values_from = mean_response)

df.1strain <- net_interactions %>% filter(Complexity > 0, Strain == 100)
interaction.df <- cbind(interaction.df, df.1strain[,c(2:9,20)])

interaction.matrix <- as.matrix(interaction.df[,2:13])
rownames(interaction.matrix) <- interaction.df$id

interaction.matrix <- t(interaction.matrix)
rownames(interaction.matrix) <- c("100 P. baetica", "302 F. glaciei", "306 A. humicola", "331 P. baetica", "371 R. herbae",
                                  "419 S. faeni", "448 C. gallinarum",  "487 A. popoffii",  "527 A. humicola", "74 N. soli", 
                                  "E. coli", "A. fischeri"
                                 )

phenotypic.distance <- dist(interaction.matrix, diag = TRUE, upper = TRUE)
phenotypic.distance <- as.matrix(phenotypic.distance, labels=TRUE)

# do we need to get the columns in the same order?
col.order <- c("448 C. gallinarum", "527 A. humicola", "306 A. humicola", "74 N. soli", "331 P. baetica",
               "100 P. baetica", "487 A. popoffii", "A. fischeri", "E. coli", "371 R. herbae", "419 S. faeni",
               "302 F. glaciei")
# clustered_dist <- clustered_dist[col.order,col.order]
phenotypic.distance <- phenotypic.distance[col.order,col.order]

# vegan lets us choose a method for the mantel test
mantel(phylo_dist, phenotypic.distance, method = "kendall", permutations = 9999)

# no correlation between matrices

########################################
## -- Supplementary Figs S3 and S4 -- ##
########################################

# display the comparisons of distance matrices

comparison_df <- data.frame(phylo.distance = phylo_dist[upper.tri(phylo_dist, diag = FALSE)],
                            stressor.distance = phenotypic.distance[upper.tri(phenotypic.distance, diag = FALSE)])


distance_comparison <- ggplot(comparison_df, aes(x = phylo.distance, y = stressor.distance)) +
  geom_point(size = 3) +
  labs(x = "Phylogenetic Distance", y = "Phenotypic Distance") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))
distance_comparison
ggsave("results/supplementary-fig-S3.svg", distance_comparison, width = 7, height = 7)

# --- Repeat without Oxytetracycline --- #

interaction.nooxy <- net_interactions %>%
  filter(Complexity > 0) %>%
  filter(Strain != "mixture" & Oxytetracycline == 0) %>%
  select(id, mean_response, Strain) %>%
  pivot_wider(names_from = Strain, values_from = mean_response)

df.1strain.nooxy <- net_interactions %>% filter(Complexity > 0, Strain == 100,
                                          Oxytetracycline == 0)
interaction.nooxy <- cbind(interaction.nooxy, df.1strain.nooxy[,c(2:9,20)])

interaction.matrix.nooxy <- as.matrix(interaction.nooxy[,2:13])
rownames(interaction.matrix.nooxy) <- interaction.nooxy$id

interaction.matrix.nooxy <- t(interaction.matrix.nooxy)
rownames(interaction.matrix.nooxy) <- c("100 P. baetica", "302 F. glaciei", "306 A. humicola", "331 P. baetica", "371 R. herbae",
                                  "419 S. faeni", "448 C. gallinarum",  "487 A. popoffii",  "527 A. humicola", "74 N. soli", 
                                  "E. coli", "A. fischeri"
)

phenotypic.distance.nooxy <- dist(interaction.matrix.nooxy, diag = TRUE, upper = TRUE)
phenotypic.distance.nooxy <- as.matrix(phenotypic.distance.nooxy, labels=TRUE)

phenotypic.distance.nooxy <- phenotypic.distance.nooxy[col.order,col.order]

# vegan lets us choose a method for the mantel test
mantel(phylo_dist, phenotypic.distance.nooxy, method = "kendall", permutations = 9999)

comparison_df_nooxy <- data.frame(phylo.distance = phylo_dist[upper.tri(phylo_dist, diag = FALSE)],
                            stressor.distance = phenotypic.distance.nooxy[upper.tri(phenotypic.distance.nooxy, diag = FALSE)])


distance_comparison_nooxy <- ggplot(comparison_df_nooxy, aes(x = phylo.distance, y = stressor.distance)) +
  geom_point(size = 3) +
  labs(x = "Phylogenetic Distance", y = "Phenotypic Distance") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))
distance_comparison_nooxy
ggsave("results/supplementary-fig-S4.svg", distance_comparison_nooxy, width = 7, height = 7)


###########################################################
## --- Testing phylogenetic signal with lambda and K --- ##
###########################################################

# try the phylogenetic signal of each single chem response
# maybe the antibiotics have phylo signal but the others don't?
# or something like that.

library(geiger)
library(phytools)

# need the data
single_stressor_data <- net_interactions %>% 
  filter(Complexity == 1) %>%
  mutate(Stressor = "")

# add "stressor" data
# Loop across the single stressor info to make sure they all get the right stressor assigned to the right effect 
for (p in 1:nrow(single_stressor_data)){
  for (q in 2:9){ # these are the columns where we've got the stressor names
    if (single_stressor_data[p,q] == 1){
      single_stressor_data$Stressor[p] <- colnames(single_stressor_data[q])
    }
  }
}


trait_data <- data.frame(single_stressor_data %>%
                           select(Strain, mean_response, Stressor) %>%
                           pivot_wider(names_from = Stressor, values_from = mean_response) %>%
                           filter(Strain != "mixture"))

row.names(trait_data) <- trait_data$Strain

# fix tip labels again...
phylo_tree$tip.label <- c(448, 527, 306, 74, 331, 100, 487, "A. fischeri", "E. coli", 371, 419, 302)

name.check(phylo_tree, trait_data)

Amoxicillin <- trait_data[,"Amoxicillin"]
Chlorothalonil <- trait_data[,"Chlorothalonil"]
Diflufenican <- trait_data[,"Diflufenican"]
Glyphosate <- trait_data[,"Glyphosate"]
Imidacloprid <- trait_data[,"Imidacloprid"]
Metaldehyde <- trait_data[,"Metaldehyde"]
Oxytetracycline <- trait_data[,"Oxytetracycline"]
Tebuconazole <- trait_data[,"Tebuconazole"]

names(Amoxicillin) <- names(Chlorothalonil) <- names(Diflufenican) <- names(Glyphosate) <-
  names(Imidacloprid) <- names(Metaldehyde) <- names(Oxytetracycline) <- names(Tebuconazole) <- rownames(trait_data)

# check for normality
ggplot(single_stressor_data %>% filter(Strain != "Mixture"), aes(x = mean_response)) +
  geom_density() +
  facet_wrap(~Stressor, scales = "free")
# A, C, D and I are a bit left-skewed
histogram(Amoxicillin^4)
histogram(Chlorothalonil^4)
histogram(Diflufenican^3)
histogram(Imidacloprid^10)

# signal using lambda
Amox_lambda <- phylosig(phylo_tree, Amoxicillin^4, method = "lambda", test = T)
Chlor_lambda <- phylosig(phylo_tree, Chlorothalonil^4, method = "lambda", test = T)
Dif_lambda <- phylosig(phylo_tree, Diflufenican^3, method = "lambda", test = T)
Gly_lambda <- phylosig(phylo_tree, Glyphosate, method = "lambda", test = T)
Imi_lambda <- phylosig(phylo_tree, Imidacloprid^10, method = "lambda", test = T)
Met_lambda <- phylosig(phylo_tree, Metaldehyde, method = "lambda", test = T)
Oxy_lambda <- phylosig(phylo_tree, Oxytetracycline, method = "lambda", test = T)
Teb_lambda <- phylosig(phylo_tree, Tebuconazole, method = "lambda", test = T)

# signal using K
Amox_K <- phylosig(phylo_tree, Amoxicillin^4, method = "K", test = T)
Chlor_K <- phylosig(phylo_tree, Chlorothalonil^4, method = "K", test = T)
Dif_K <- phylosig(phylo_tree, Diflufenican^3, method = "K", test = T)
Gly_K <- phylosig(phylo_tree, Glyphosate, method = "K", test = T)
Imi_K <- phylosig(phylo_tree, Imidacloprid^10, method = "K", test = T)
Met_K <- phylosig(phylo_tree, Metaldehyde, method = "K", test = T)
Oxy_K <- phylosig(phylo_tree, Oxytetracycline, method = "K", test = T)
Teb_K <- phylosig(phylo_tree, Tebuconazole, method = "K", test = T)

# stick them together into a dataframe for plotting
signal_results <- data.frame(test.statistic = c(Amox_lambda$lambda,
                                                Chlor_lambda$lambda,
                                                Dif_lambda$lambda,
                                                Gly_lambda$lambda,
                                                Imi_lambda$lambda,
                                                Met_lambda$lambda,
                                                Oxy_lambda$lambda,
                                                Teb_lambda$lambda,
                                                Amox_K$K,
                                                Chlor_K$K,
                                                Dif_K$K,
                                                Gly_K$K,
                                                Imi_K$K,
                                                Met_K$K,
                                                Oxy_K$K,
                                                Teb_K$K),
                             p.val = c(Amox_lambda$P,
                                                Chlor_lambda$P,
                                                Dif_lambda$P,
                                                Gly_lambda$P,
                                                Imi_lambda$P,
                                                Met_lambda$P,
                                                Oxy_lambda$P,
                                                Teb_lambda$P,
                                                Amox_K$P,
                                                Chlor_K$P,
                                                Dif_K$P,
                                                Gly_K$P,
                                                Imi_K$P,
                                                Met_K$P,
                                                Oxy_K$P,
                                                Teb_K$P))

signal_results$test <- rep(c("lambda", "K"), each = 8)
signal_results$stressor <- rep(c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate",
                                 "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole"), 2)

# make it wider for a radar plot
# (this isn't very ggplot like!)
signal_results_wide <- signal_results %>%
  select(-p.val) %>%
  pivot_wider(values_from = test.statistic, names_from = stressor)

library(ggradar)
svg("results/fig-3a.svg", width = 8, height = 8)
ggradar(signal_results_wide,
        values.radar = element_blank(),
        grid.min = 0,
        grid.mid = 0.5,
        grid.max = 1,
        legend.position = "bottom", 
        legend.text.size = 0,
        axis.label.size = 0)
dev.off()


##########################################
# Plot the phylogenetic tree onto a map! #
##########################################

# turns out its too hard without a proper shapefile
# so just going to plot the tree nicely and do the lines manually!
iceland_tree <- drop.tip(rerooted_tree, c("E. coli", "A. fischeri"))

svg("results/phylo-iceland.svg", width = 4, height= 10)
plot(iceland_tree)
dev.off()