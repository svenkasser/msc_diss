# Primary ABM code

# Set up packages

## For data handling

library(tidyverse)
library(dplyr)

## For parallel processing

library(parallel)
library(doParallel)

## Set up parallel processing

list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "kableExtra"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

### loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

n.cores <- parallel::detectCores() - 1

### create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

### check cluster definition (optional)
print(my.cluster)

### register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

### check if it is registered (optional)
foreach::getDoParRegistered()

foreach::getDoParWorkers()


## Set up ABM


ABMmodel6 <- function(N, t_max, r_max, pc_0 = 0.5, pb_0 = 0.5, b_c = 0.5, b_b = 0.5, cb_bh = 0.5, pf_0 = 0.5, b_f = 0.5, mono_max = 3, poly_max = 2) {
  
  # Set parameter values for the parallel loop to access
  
  N = N
  t_max = t_max
  r_max = r_max
  pc_0 = pc_0
  pb_0 = pb_0
  b_c = b_c
  b_b = b_b
  cb_bh = cb_bh
  pf_0 = pf_0
  b_f = b_f
  mono_max = mono_max
  poly_max = poly_max
  
  # Loop for each run (essentially each simulated population) - uses parallel processing capabilities
  
  final_output <- foreach(r = 1:r_max, .packages = 'tidyverse', .export = c("matingmarket_sim2", "trait_determination"), .combine = 'rbind') %dopar% {
    
    r = r
    
    # Create first generation
    population <- tibble(sex = sample(c("Female", "Male"), N, replace = TRUE, prob = c(pf_0, 1 - pf_0)),
                         trait_c = sample(c("Monogamy", "Polygamy"), N, replace = TRUE, prob = c(pc_0, 1 - pc_0)), 
                         trait_b = sample(c("Monogamy", "Polygamy"), N, replace = TRUE, prob = c(pb_0, 1 - pb_0)),
                         bh = as.character(rep(NA, N)))
    
    # Determining behavioural tendencies based on cultural and biological trait combinations
    
    # Both traits favour monogamy, thus behavioural preference is monogamous
    full_mono <- population$trait_c == "Monogamy" & population$trait_b == "Monogamy"
    if (sum(full_mono) > 0) {
      population[full_mono, ]$bh <- "Monogamy"  
    }
    
    # Both traits favour polygamy, thus behavioural preference is polygamous
    full_poly <- population$trait_c == "Polygamy" & population$trait_b == "Polygamy"
    if (sum(full_poly) > 0) {
      population[full_poly, ]$bh <- "Polygamy" 
    }
    
    # If any empty NA slots (i.e. mixed cultural/biological traits) are present
    if (anyNA(population$bh)) {  
      # They will adopt the behavioural variant corresponding to their cultural trait with probability cb_bh 
      population[is.na(population$bh)& population$trait_c == "Monogamy",]$bh <- 
        sample(c("Monogamy", "Polygamy"), sum(is.na(population$bh)&population$trait_c == "Monogamy"), prob = c(cb_bh, 1 - cb_bh), replace = TRUE)
      population[is.na(population$bh)& population$trait_c == "Polygamy",]$bh <-
        sample(c("Polygamy", "Monogamy"), sum(is.na(population$bh)&population$trait_c == "Polygamy"), prob = c(cb_bh, 1 - cb_bh), replace = TRUE)
    }
    
    # Create output file - this is he data what we want to save from each run and generation
    
    output <- tibble(generation = 1, 
                     p_c = as.numeric(rep(NA, 1)),
                     p_b = as.numeric(rep(NA, 1)),
                     p_bh = as.numeric(rep(NA, 1)),
                     p_f = as.numeric(rep(NA, 1)),
                     p_mono = as.numeric(rep(NA, 1)),
                     p_polyg = as.numeric(rep(NA,1)),
                     p_polya = as.numeric(rep(NA,1)),
                     poly_partners = as.numeric(rep(NA, 1)),
                     poly_partners_m = as.numeric(rep(NA, 1)),
                     poly_partners_f = as.numeric(rep(NA, 1)),
                     RS_mono = as.numeric(rep(NA,1)),
                     RS_poly = as.numeric(rep(NA,1)),
                     RS_mono_m = as.numeric(rep(NA,1)),
                     RS_poly_m = as.numeric(rep(NA,1)),
                     RS_mono_f = as.numeric(rep(NA,1)),
                     RS_poly_f = as.numeric(rep(NA,1)),
                     run = as.factor(rep(r, 1)))
    
    
    # Add first generation's p_c (proportion of cultural Monogamy trait in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p_c <- 
      sum(population$trait_c == "Monogamy") / N 
    
    # Add first generation's p_b (proportion of biological Monogamy trait in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p_b <- 
      sum(population$trait_b == "Monogamy") / N 
    
    # Add first generation's p_bh (proportion of Monogamous behaviour in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p_bh <- 
      sum(population$bh == "Monogamy") / N 
    
    # Add first generation's p_f (proportion of female individuals, analogous to sex-ratio) for run r
    output[output$generation == 1 & output$run == r, ]$p_f<- 
      sum(population$sex == "Female") / N 
    
    # While the first generation is always generated in the same way, each subsequent generation needs to go through the following loop
    
    for (t in 2:t_max) {
      
      # Copy current population to previous_population tibble
      
      previous_population <- population
      
      # Create matng pool from previous population
      # This is the input fpr te mating market sim
      
      matingpool <- previous_population %>%
        add_column(matesearched = "No", available = "Yes", bond_ID1 = as.character("NA"), bond_ID2 = as.character("NA"), partners = 0, offspring = 0, ID = 1:N)
      
      # Create empty offspring tibble, which includes data on the parents traits
      # This is the empty table of offspring that will constitute the next generation once filled by the mating loop
      
      matepairings <- tibble(sex = sample(c("Female", "Male"), N, replace = TRUE, prob = c(b_f, 1 - b_f)),
                             ID_p1 = as.character(rep(NA, N)),
                             trait_c_p1 = as.character(rep(NA, N)),
                             trait_b_p1 = as.character(rep(NA, N)),
                             ID_p2= as.character(rep(NA, N)),
                             trait_c_p2 = as.character(rep(NA, N)),
                             trait_b_p2 = as.character(rep(NA, N)))
      
      # This is the output for the mating market submodel
      
      matingmodel_output <- matingmarket_sim2(matingpool, matepairings, N = N, mono_maxRS = mono_max, poly_maxRS = poly_max)
      
      # The mating market simulation is self-contained and does not automatically update the matingpool sitting outside of it
      # Here, we copy the reproductive success data recorded within the mating market simulation back into the main loop, to be recorded in the output table
      
      matingpool <- as_tibble(cbind(matingmodel_output[2]$input$sex,
                                    matingmodel_output[2]$input$trait_c,
                                    matingmodel_output[2]$input$trait_b,
                                    matingmodel_output[2]$input$bh,
                                    matingmodel_output[2]$input$matesearched,
                                    matingmodel_output[2]$input$available,
                                    matingmodel_output[2]$input$bond_ID1,
                                    matingmodel_output[2]$input$bond_ID2,
                                    matingmodel_output[2]$input$partners,
                                    matingmodel_output[2]$input$offspring,
                                    matingmodel_output[2]$input$ID))
      
      matingpool <- rename(matingpool, sex = V1, trait_c = V2, trait_b = V3, bh = V4, matesearched = V5, availabe = V6, bond_ID1 = V7, bond_ID2 = V8, partners = V9, offspring = V10, ID = V11)
      
      # The same applies to the offspring table, labelled "nextgen"
      
      nextgen <- as_tibble(cbind(matingmodel_output[1]$output$sex, 
                                 matingmodel_output[1]$output$trait_c_p1, 
                                 matingmodel_output[1]$output$trait_b_p1, 
                                 matingmodel_output[1]$output$trait_c_p2, 
                                 matingmodel_output[1]$output$trait_b_p2))
      
      nextgen <- rename(nextgen, sex = V1, trait_c_p1 = V2, trait_b_p1 = V3, trait_c_p2 = V4 , trait_b_p2 = V5)
      
      # Create an empty population tibble for the nextgen tibble to feed into
      
      population <- select(nextgen, -c(trait_c_p1, trait_c_p2, trait_b_p1, trait_b_p2)) %>%
        add_column(trait_c = as.character(rep(NA, N)), 
                   trait_b = as.character(rep(NA, N)),
                   bh = as.character(rep(NA, N)))
      
      # Run trait determination submodel on nextgen tibble, then record otcome in population tibble
      
      population <- trait_determination(population, nextgen, b_c = b_c, b_b = b_b, cb_bh = cb_bh)
      
      # Add new row of output
      
      new_output <- tibble(generation = t, 
                           p_c = as.numeric(rep(NA, 1)),
                           p_b = as.numeric(rep(NA, 1)),
                           p_bh = as.numeric(rep(NA, 1)),
                           p_f = as.numeric(rep(NA, 1)),
                           p_mono = as.numeric(rep(NA, 1)),
                           p_polyg = as.numeric(rep(NA,1)),
                           p_polya = as.numeric(rep(NA,1)),
                           poly_partners = as.numeric(rep(NA, 1)),
                           poly_partners_m = as.numeric(rep(NA, 1)),
                           poly_partners_f = as.numeric(rep(NA, 1)),
                           RS_mono = as.numeric(rep(NA,1)),
                           RS_poly = as.numeric(rep(NA,1)),
                           RS_mono_m = as.numeric(rep(NA,1)),
                           RS_poly_m = as.numeric(rep(NA,1)),
                           RS_mono_f = as.numeric(rep(NA,1)),
                           RS_poly_f = as.numeric(rep(NA,1)),
                           run = as.factor(rep(r, 1)))
      
      # Get p_c (cultural trait frequency for "Monogamy") and put it into output slot for this generation t and run r
      new_output[new_output$generation == t & new_output$run == r, ]$p_c <- 
        sum(population$trait_c == "Monogamy") / N 
      
      # Get p_b (Biological trait frequency for "Monogamy") and put it into output slot for this generation t and run r
      new_output[new_output$generation == t & new_output$run == r, ]$p_b <- 
        sum(population$trait_b == "Monogamy") / N 
      
      # Get p_bh (Behavioural trait frequency for "Monogamy") and put it into output slot for this generation t and run r
      new_output[new_output$generation == t & new_output$run == r, ]$p_bh <- 
        sum(population$bh == "Monogamy") / N 
      
      # Get the sex ratio p_f (as proportion of female individuals) for this generation t in run r
      new_output[new_output$generation == t & new_output$run == r, ]$p_f <- 
        sum(population$sex == "Female") / N 
      
      # Get observable monogamy outcome from mating pool
      output[output$generation == t-1 & output$run == r, ]$p_mono <-
        sum(matingpool$partners == 1) / N
      
      # Get average reproductive success (number of offspring) for polygamous agents
      output[output$generation == t-1 & output$run == r, ]$RS_poly <-
        mean(as.integer(matingpool[matingpool$bh=="Polygamy",]$offspring))
      
      # Get average reproductive success (number of offspring) for monogamous agents
      output[output$generation == t-1 & output$run == r, ]$RS_mono <-
        mean(as.integer(matingpool[matingpool$bh=="Monogamy",]$offspring))
      
      # Get freqency of observable polygyny (i.e. polygynous outcomes regardless of traits) among male agents
      output[output$generation == t-1 & output$run == r, ]$p_polyg <-
        sum(matingpool$partners == 2 & matingpool$sex == "Male") / sum(matingpool$sex=="Male")
      
      # Get freqency of observable polyandry (i.e. polyandrous outcomes regardless of traits) among female agents
      output[output$generation == t-1 & output$run == r, ]$p_polya <-
        sum(matingpool$partners == 2 & matingpool$sex == "Female") / sum(matingpool$sex=="Female")
      
      # Get average reproductive success (number of offspring) for male polygamous agents
      output[output$generation == t-1 & output$run == r, ]$RS_poly_m <-
        mean(as.integer(matingpool[matingpool$bh=="Polygamy" & matingpool$sex=="Male",]$offspring))
      
      # Get average reproductive success (number of offspring) for female monogamous agents
      output[output$generation == t-1 & output$run == r, ]$RS_mono_f <-
        mean(as.integer(matingpool[matingpool$bh=="Monogamy" & matingpool$sex=="Female",]$offspring))
      
      # Get average reproductive success (number of offspring) for female polygamous agents
      output[output$generation == t-1 & output$run == r, ]$RS_poly_f <-
        mean(as.integer(matingpool[matingpool$bh=="Polygamy" & matingpool$sex=="Female",]$offspring))
      
      # Get average reproductive success (number of offspring) for male monogamous agents
      output[output$generation == t-1 & output$run == r, ]$RS_mono_m <-
        mean(as.integer(matingpool[matingpool$bh=="Monogamy" & matingpool$sex=="Male",]$offspring))
      
      # Get average number of mating partners (mating success) for polygamous agents
      output[output$generation == t-1 & output$run == r, ]$poly_partners <-
        mean(as.integer(matingpool[matingpool$bh=="Polygamy",]$partners))
      
      # Get average number of mating partners (mating success) for male polygamous agents
      output[output$generation == t-1 & output$run == r, ]$poly_partners_m <-
        mean(as.integer(matingpool[matingpool$bh=="Polygamy" & matingpool$sex=="Male",]$partners))
      
      # Get average number of mating partners (mating success) for female polygamous agents
      output[output$generation == t-1 & output$run == r, ]$poly_partners_f <-
        mean(as.integer(matingpool[matingpool$bh=="Polygamy" & matingpool$sex=="Female",]$partners))
      
      # Add new row of output data to existing output tibble
      
      output <- rbind(output,new_output)
      
    }
    return(output)
  }
  
  # Export data from function
  final_output
  
}

