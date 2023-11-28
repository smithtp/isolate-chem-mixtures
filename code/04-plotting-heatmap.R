#########################################
# plotting-heatmap.R
#
# Create a large heatmap detailing the
# responses of every bacteria to every
# mixture of chemicals. Shown in Fig 1A
#
# Tom Smith 2023
#########################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap) # package for fancy heatmap plots
library(circlize) # do this for making colour palettes

# read the data
net_interactions <- read.csv("results/bootstrapped-net-interactions.csv")

# add the test of responses being additive, or just no response
noresponse_test <- read.csv("results/bootstrapped-no-response-test.csv") %>%
  select(Strain, id, significant) %>%
  rename(additive.significant = significant)

net_interactions <- left_join(net_interactions, noresponse_test)

net_interactions$response.tested <- net_interactions$response
net_interactions[net_interactions$response == "Additive" & 
                   net_interactions$additive.significant == FALSE,]$response.tested <- "None"

# where there is no response, make the dAUC value 1
net_interactions$dAUC <- net_interactions$mean_response
net_interactions[net_interactions$response.tested == "None",]$dAUC <- 1

strain_matching <- data.frame(Strain = c("100", "302", "306", "331", "371", "419", "448", "487", "527",
                                         "74", "Ecoli", "mixture", "vfischeri"),
                              Strain_nice = c("P. baetica 1", "F. glaciei", "A. humicola 1", "P. baetica 2", "R. herbae",
                                              "S. faeni", "C. gallinarum", "A. popoffii", "A. humicola 2", "N. soli", 
                                              "E. coli", "Mixture", "A. fischeri"))

net_interactions <- left_join(net_interactions, strain_matching)

####################################
# ----- MAKE HEATMAP FIGURES ----- #
####################################

# make a colour function
col_fun <- colorRamp2(c(0, 1, 2), c("firebrick", "white", "slateblue"))

# data needs to be in a matrix format
interaction.df <- net_interactions %>%
  filter(Complexity > 0) %>%
  #filter(Strain != "mixture") %>% # nah, lets keep the mixture now
  select(id, dAUC, Strain_nice) %>%
  pivot_wider(names_from = Strain_nice, values_from = dAUC)

df.1strain <- net_interactions %>% filter(Complexity > 0, Strain == 100)
interaction.df <- cbind(interaction.df, df.1strain[,c(2:9,20)])

# figure out how to reorder the rows the way we want them
interaction.matrix <- as.matrix(interaction.df[,2:14])
rownames(interaction.matrix) <- interaction.df$id

# second heatmap of simply the stressor matrix
stressor.matrix <- as.matrix(interaction.df[,15:22])
pal <- c("white", "black")

# do something cunning to highlight the E. Coli and A. fischeri columns
# Rows to highlight
labstrainCols <- c("E. coli", "A. fischeri")
mixCol <- "Mixture"

# Set stylings for column names and make our selected columns unique
lab_col_idx <- which(colnames(interaction.matrix) %in% labstrainCols)
mix_col_idx <- which(colnames(interaction.matrix) %in% mixCol)
fontcolors <- rep("black", ncol(interaction.matrix))
fontcolors[lab_col_idx] <- "purple"
fontfaces <- rep("plain", ncol(interaction.matrix))
fontfaces[mix_col_idx] <- "bold"

# Create text annotation object for displaying column names
colAnno <- columnAnnotation(cols = anno_text(colnames(interaction.matrix), gp = gpar(col = fontcolors, fontface = fontfaces)))

## Create the main heatmap
h1 <- Heatmap(interaction.matrix, col = col_fun, 
              name = "Growth Reduction",
              #cluster_rows = FALSE, # turn off row clustering
              heatmap_legend_param = list(title = "G"),
              bottom_annotation = colAnno, # use our column annotation
              show_column_names = FALSE, # and don't use the normal column names
              # cluster_row_slices = FALSE,
              column_title = "Bacterial Strain",  column_title_side = "bottom",
              row_title = NULL,
              show_row_names = FALSE)

# create the stressor matrix heatmap
h2 <- Heatmap(stressor.matrix, col = c("white", "black"),
              name = "Present/Absent",
              cluster_rows = FALSE, cluster_columns = FALSE, 
              column_order = c("Amoxicillin", "Oxytetracycline", "Chlorothalonil", "Tebuconazole",
                               "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde"),
              column_title = "Chemical Stressor",  column_title_side = "bottom",
              row_title = NULL)

# create the row annotation
new_cols <- colorRamp2(c(0, 8), c("white", "orange"))
h3 <- rowAnnotation(Complexity = interaction.df[,c("Complexity")],  col = list(Complexity = new_cols))

ht_list <- h1 + h3 + h2

# Export the heatmap list figure
svg("results/heatmap-figure.svg", width = 10, height = 8)
draw(ht_list)
dev.off()