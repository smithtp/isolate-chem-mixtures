#######################################
# Combine the raw data from each plate
# and do a little cleaning up
#
#

library(dplyr)
library(tidyr)

# functions to apply to plates
remove.temp <- function(x){
  # remove the temperature column
  y <- x[,c(1,3:98)]
  return(y)
}

# clean up the time column
clean.time <- function(x){
  x$Time <- sapply(strsplit(x$Time, ":"),
                   function(z) {
                     z <- as.numeric(z)
                     z[1]+(z[2]+z[3]/60)/60
                   })
  return(x) 
}

# some plates got put in backwards by mistake, lets fix that here
re.number.cols <- function(x){
  the.names <- names(x)
  new.names <- c(the.names[1], rev(the.names[2:97])) 
  names(x) <- new.names
  return(x)
}

# read in list of files from different plates


clean.raw.data <- function(isolate.name, plate.name){
  filepath <- paste("data/", "isolate-", isolate.name, "-plate-", plate.name, ".txt", sep = "")
  if(file.exists(filepath)){ # check the file exists
    plate_data <- clean.time(remove.temp(read.table(filepath,  sep = "\t", header = TRUE, check.names = FALSE)))
    # edit the time column of the plates
    # because it takes the stacker time to run each plate...
    # we're adding 100 seconds to each sequential plate (1min 40)
    plate_data$Time <- plate_data$Time+((plate.name*100)-100)/3600
    
    # renumber the plates that we fucked up putting into the robot
    if(isolate.name == 100 & plate.name %in% c(12:14)){
      plate_data <- re.number.cols(plate_data)
    }
    # if(isolate.name == 527 & plate.name %in% c(4,9,14,19)){
    #   plate_data <- re.number.cols(plate_data)
    # } # I think this was for Jenny's data, but not Tom M-Ls!
    
    # number the plates for the chemical matrix
    if(plate.name %in% seq(1, 16, by = 5)){
      plate.num <- 1
    } else if(plate.name %in% seq(2, 17, by = 5)){
      plate.num <- 2
    } else if(plate.name %in% seq(3, 18, by = 5)){
      plate.num <- 3
    } else if(plate.name %in% seq(4, 19, by = 5)){
      plate.num <- 4
    } else if(plate.name %in% seq(5, 20, by = 5)){
      plate.num <- 5
    } else{
      print("Error: plate number can't be resolved")
    }
    
    # number the replicate plates
    rep.num <- ceiling(plate.name/5)
    
    clean_data <- pivot_longer(plate_data, cols = c(2:97), names_to = "Well", values_to = "OD600") %>%
      arrange(Well) %>%
      mutate(Plate = plate.num) %>%
      mutate(PlateTrue = plate.name) %>%
      mutate(Rep = rep.num) %>%
      mutate(Strain = isolate.name)
    
    return(clean_data)
  } else{
    return(NA)
  }
}



# run through a big ol' loop to collect all the data
clean_data <- data.frame()

isolate.names <- c(74, 100, 302, 306, 331, 371, 419, 448, 487, 527, "Ecoli", "vfischeri", "mixture")
plate.names <- seq(1, 20, 1)


for(i in seq_along(isolate.names)){
  for(j in seq_along(plate.names)){
    # need some error catching for if some plates don't exist
    cleaned_data <- clean.raw.data(isolate.names[i], plate.names[j])
    if(is.data.frame(cleaned_data)){
      clean_data <- bind_rows(clean_data, cleaned_data)
    }
  }
}

# read the (previously created) chemical matrix
chem.plates <- read.csv("data/chemical-matrix.csv")

# combine OD readings with chemical matrix
combined_data <- left_join(clean_data, chem.plates, by = c("Plate", "Well"))

# create a column combining plate and well
combined_data$PlateWell <- paste(combined_data$Plate, combined_data$Well, sep = "")

# create a column for the True plate well
combined_data$TruePlateWell <- paste(combined_data$PlateTrue, combined_data$Well, sep = "")

# remove the wells without bacteria
tidy_data <- combined_data[!is.na(combined_data$Complexity),]

# and do a quick fix for Strain 419 where the Tebuconazole and Oxytet got swapped
# around when the plates were made
data_419 <- tidy_data %>% 
  filter(Strain == 419) %>%
  rename(Tebuconazole = Oxytetracycline, Oxytetracycline = Tebuconazole)

data_remaining <- tidy_data %>% 
  filter(Strain != 419)

tidy_data <- bind_rows(data_remaining, data_419)

rm(data_419, data_remaining)

######################################
###   Data cleaning
############################

# remove contaminated plates
tidy_data <- tidy_data[!(tidy_data$Strain == 302 &
                           tidy_data$PlateTrue %in% c(4, 9, 18)),]

# something weird happened in a row of one of the #302 plates
tidy_data <- tidy_data[!(tidy_data$Strain == 302 &
                           tidy_data$TruePlateWell %in% c("6D1", "6D2", "6D3", "6D4", "6D5", "6D6", 
                                                          "6D7", "6D8", "6D9", "6D10", "6D11", "6D12")),]

# remove dodgy well
tidy_data <- tidy_data[!(tidy_data$Strain == 331 &
                           tidy_data$TruePlateWell == "20B2"),]

# Filter out some dodgy data from isolate #100 which had weird starting values
tidy_data <- tidy_data[!(tidy_data$Strain == 100 &
                           tidy_data$TruePlateWell %in% tidy_data[tidy_data$Strain == 100 &
                                                                    tidy_data$Time < 1 &
                                                                    tidy_data$OD600 > 0.2,]$TruePlateWell),]

# Filter out some dodgy data from isolate #448 with weird values (evaporation?)
tidy_data <- tidy_data[!(tidy_data$Strain == 448 &
                           tidy_data$TruePlateWell %in% c("1G2", "7C4", "7C5", "7E10", "8B2", "8B7",
                                                          "9D7", "12E10", "13B7")),]

# # Filter a some dodgy wells from isolate # 527. 6B5 got no bacteria
# tidy_data <- tidy_data[!(tidy_data$Strain == 527 &
#                            tidy_data$TruePlateWell == "6B5"),]

# Many problem wells from Tom M-Ls first #527 run -
# looks like multichannel pipetting error mostly
# 1B6-G6 ; 16C2-G2
# 2B10; 7B2-G2; 17B3-F3; 12C11-G11; 17C9-E9
# 8G7, 9B4
# 20B11-G11; 10B3-G3; 20B3-G3
tidy_data <- tidy_data[!(tidy_data$Strain == 527 &
                           tidy_data$TruePlateWell %in% c("1B6", "1C6", "1D6", "1E6", "1F6", "1G6",
                                                          "16C2", "16D2", "16E2", "16F2", "16G2",
                                                          "2B10", "7B2", "7C2", "7D2", "7E2", "7F2", "7G2",
                                                          "17B3", "17C3", "17D3", "17E3", "17F3", "12C11", "12D11",
                                                          "12E11", "12F11", "12G11", "17C9", "17D9", "17E9",
                                                          "8G7", "9B4", "20B11", "20C11", "20D11", "20E11",
                                                          "20F11", "20G11", "10B3", "10C3", "10D3", "10E3", "10F3",
                                                          "10G3", "20B3", "20C3", "20D3", "20E3", "20F3", "20G3")),]

# Filter out a couple of wells from the Ecoli run (one OD starts too high, one looks like it failed to get any bacteria in)
# 13G8 got no bacteria # 8F8 got a double-helping, or had a bubble
tidy_data <- tidy_data[!(tidy_data$Strain == "Ecoli" &
                           tidy_data$TruePlateWell %in% c("13G8", "8F8")),]

# make the strain names nicer
tidy_data$Strain_nice <- NA

isolate.names <- c(74, 100, 302, 306, 331, 371, 419, 448, 487, 527, "Ecoli", "vfischeri", "mixture")

tidy_data[tidy_data$Strain == 74,]$Strain_nice <- "IL 74"
tidy_data[tidy_data$Strain == 100,]$Strain_nice <- "IL 100"
tidy_data[tidy_data$Strain == 302,]$Strain_nice <- "IL 302"
tidy_data[tidy_data$Strain == 306,]$Strain_nice <- "IL 306"
tidy_data[tidy_data$Strain == 331,]$Strain_nice <- "IL 331"
tidy_data[tidy_data$Strain == 371,]$Strain_nice <- "IL 371"
tidy_data[tidy_data$Strain == 419,]$Strain_nice <- "IL 419"
tidy_data[tidy_data$Strain == 448,]$Strain_nice <- "IL 448"
tidy_data[tidy_data$Strain == 487,]$Strain_nice <- "IL 487"
tidy_data[tidy_data$Strain == 527,]$Strain_nice <- "IL 527"
tidy_data[tidy_data$Strain == "Ecoli",]$Strain_nice <- "E. coli"
tidy_data[tidy_data$Strain == "vfischeri",]$Strain_nice <- "V. fischeri"
tidy_data[tidy_data$Strain == "mixture",]$Strain_nice <- "IL Mixture"

tidy_data$Strain_nice <- factor(tidy_data$Strain_nice, levels = c("IL 74", "IL 100", "IL 302", "IL 306", "IL 331", "IL 371",
                                                                     "IL 419", "IL 448", "IL 487", "IL 527", "E. coli",
                                                                     "V. fischeri", "IL Mixture"))

write.csv(tidy_data, "data/tidy-data.csv", row.names = FALSE)
