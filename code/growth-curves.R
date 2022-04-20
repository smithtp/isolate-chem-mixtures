#############################
# growth curves from multiple plates

library(tidyr)
library(dplyr)
library(ggplot2)


# read the data

tidy_data <- read.csv("data/tidy-data.csv")


###############################
# Model fitting functions
###############################

# note that gompertz is expressed here with log10 (i.e. log10 transform the OD data before use)
gompertz_model <- function(N_max, N_0, r_max, t_lag, t){ # Modified gompertz growth model (Zwietering 1990)
  return(N_0 + (N_max - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t)/((N_max - N_0) * log(10)) + 1)))
}

logistic_model <- function(N_max, N_0, r_max, t) {
  return( N_max / (1 + ((N_max - N_0) / N_0) * exp(-r_max * t)))
}


AreaUnderGompertz <- function(N_max, N_0, r_max, t_lag, t_min = 0, t_max) {
  auc_l <- integrate(function(x) gompertz_model(N_max, N_0, r_max, t_lag, x), t_min, t_max)
  return(auc_l)
}

AreaUnderLogistic <- function(N_max, N_0, r_max, t_min = 0, t_max) {
  auc_l <- integrate(function(x) logistic_model(N_max, N_0, r_max, x), t_min, t_max)
  return(auc_l)
}

# function to integrate spline curve, with minT and maxT
# the time range to integrate across
AreaUnderSpline <- function(Time, OD600, minT, maxT){
  spline.area <- stats::integrate(stats::splinefun(Time, OD600, method = "natural"), 
                                  lower = minT, upper = maxT)$value
  return(spline.area)
}


# -------------------- Model Fitting ----------------------- #

# create a place to store the data
tidy_growth_data <- tidy_data[0,] %>%
  select(-Time, -Well, -Plate, -OD600)

fitting_wells <- unique(tidy_data$PlateWell)
plate_reps <- unique(tidy_data$Rep)
strains <- unique(tidy_data$Strain)

# crazy loop to fit model for every well
for(i in seq_along(strains)){
  for(j in seq_along(fitting_wells)){
    # print something to say where we are in the loop
    print(paste("Fitting #", j, "of", length(fitting_wells), "for", strains[i], sep = " "))
    for(k in seq_along(plate_reps)){
      # pull out a single growth curve
      fit_data <- tidy_data[tidy_data$Strain == strains[i] &
                              tidy_data$PlateWell == fitting_wells[j] &
                              tidy_data$Rep == plate_reps[k],]
      if(nrow(fit_data) > 0){ # check we actually have this combo in the dataframe
        ## Forget fully modelling growth curves - just fit spline function and calculate area under curve
        # need to make sure to integrate across the same time frame for each
        # species, now we've altered the times of later plates
        
        # take away the lowest OD measurement from all other measurements
        # to account for any differences in volume or differences in OD
        # caused by the chemicals directly
        min.OD <- min(fit_data$OD600)
        fit_data$OD.norm <- fit_data$OD600-min.OD
        
        AUC <- AreaUnderSpline(fit_data$Time, fit_data$OD.norm, 1, 72)
      
        temp_data <- bind_cols(PlateWell = fit_data$PlateWell[[1]],
                               PlateTrue = fit_data$PlateTrue[[1]],
                               TruePlateWell = fit_data$TruePlateWell[[1]],
                               Rep = fit_data$Rep[[1]],
                               Amoxicillin = fit_data$Amoxicillin[[1]],
                               Chlorothalonil = fit_data$Chlorothalonil[[1]],
                               Diflufenican = fit_data$Diflufenican[[1]],
                               Glyphosate = fit_data$Glyphosate[[1]],
                               Imidacloprid = fit_data$Imidacloprid[[1]],
                               Metaldehyde = fit_data$Metaldehyde[[1]],
                               Oxytetracycline = fit_data$Oxytetracycline[[1]],
                               Tebuconazole = fit_data$Tebuconazole[[1]],
                               Complexity = fit_data$Complexity[[1]],
                               Strain = fit_data$Strain[[1]],
                               AUC = AUC)
      
        # tidy the results into one place
        tidy_growth_data <- bind_rows(tidy_growth_data, temp_data)
      }
    }
  }
}

# write to csv
write.csv(tidy_growth_data, "results/spline_fits.csv", row.names = FALSE)

# plot stuff to see whats happened
ggplot(tidy_growth_data, aes(x = Complexity, y = AUC)) + 
  geom_point() +
  facet_wrap(~Strain)

########################################
# create summaries of the growth curves
########################################

# Pick a growth metric. I'm using AUC, but it could be some other parameter
# if we had fit a growth model to calculate that
growth_metric <- "AUC"

# Pull out our controls so we can average them by Isolate rather than location
control_means <- tidy_growth_data %>%
  dplyr::filter(Complexity == 0) %>%
  group_by(Strain) %>%
  summarise(Mean = mean(UQ(rlang::sym(growth_metric))),
            SD = sd(UQ(rlang::sym(growth_metric))),
            n = n(),
            Amoxicillin = mean(Amoxicillin),
            Chlorothalonil = mean(Chlorothalonil),
            Diflufenican = mean(Diflufenican),
            Glyphosate = mean(Glyphosate),
            Imidacloprid = mean(Imidacloprid),
            Metaldehyde = mean(Metaldehyde),
            Oxytetracycline = mean(Oxytetracycline),
            Tebuconazole = mean(Tebuconazole),
            Complexity = mean(Complexity))

# Average the remaining treatment wells by location
growth_means <- tidy_growth_data %>%
  dplyr::filter(Complexity != 0) %>%
  group_by(Strain, PlateWell) %>%
  summarise(Mean = mean(UQ(rlang::sym(growth_metric))),
            SD = sd(UQ(rlang::sym(growth_metric))),
            n = n(),
            Amoxicillin = mean(Amoxicillin),
            Chlorothalonil = mean(Chlorothalonil),
            Diflufenican = mean(Diflufenican),
            Glyphosate = mean(Glyphosate),
            Imidacloprid = mean(Imidacloprid),
            Metaldehyde = mean(Metaldehyde),
            Oxytetracycline = mean(Oxytetracycline),
            Tebuconazole = mean(Tebuconazole),
            Complexity = mean(Complexity))

# And merge the summaries back together
summarised_growth_data <- bind_rows(growth_means, control_means)

write.csv(summarised_growth_data, "results/model_summaries.csv", row.names = FALSE)
