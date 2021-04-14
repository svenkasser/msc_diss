library(tidyverse)
library(dplyr)
library(ggplot2)

# First, I'm gonna follow the Acerbi book closely to build a simple vertical transmission model of Marriage system preference as a cultural belief

# This initial model is not going to have biological sex, or mating, and only has a features a singular cultural trait

# Only allows for a handful of parameters to vary - pop size, proportion at gen 1, likelihood of Monogamy trait in the case of mixed inheritance (I could fix this at 0.5), max number of generations simulated, and max number of populations simulated

## ABM version 1 - Vertical transmission, no mutation

vertical_transmission <- function(N, p_0, b, t_max, r_max) {
  output <- tibble(generation = rep(1:t_max, r_max), 
                   p = as.numeric(rep(NA, t_max * r_max)), 
                   run = as.factor(rep(1:r_max, each = t_max)))
  
  for (r in 1:r_max) {
    # Create first generation
    population <- tibble(trait = sample(c("Monogamy", "Polygamy"), N, 
                                        replace = TRUE, prob = c(p_0, 1 - p_0)))
    
    # Add first generation's p (proportion of Monogamy trait in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p <- 
      sum(population$trait == "Monogamy") / N 
    
    for (t in 2:t_max) {
      # Copy individuals to previous_population tibble
      previous_population <- population
      
      # Randomly pick mothers and fathers (this will have to be rewritten when sex is added)
      mother <- tibble(trait = sample(previous_population$trait, N, replace = TRUE))  
      father <- tibble(trait = sample(previous_population$trait, N, replace = TRUE)) 
      
      # Prepare next generation
      population <- tibble(trait = as.character(rep(NA, N))) 
      
      # Both parents are monogamous, thus child adopts Monogamy
      both_A <- mother$trait == "Monogamy" & father$trait == "Monogamy"
      if (sum(both_A) > 0) {
        population[both_A, ]$trait <- "Monogamy"  
      }
      
      # Both parents are polygamous, thus child adopts Polygamy
      both_B <- mother$trait == "Polygamy" & father$trait == "Polygamy"
      if (sum(both_B) > 0) {
        population[both_B, ]$trait <- "Polygamy" 
      }
      # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
      if (anyNA(population)) {  
        # They adopt Monogamy with probability b
        population[is.na(population)[,1],]$trait <- 
          sample(c("Monogamy", "Polygamy"), sum(is.na(population)), prob = c(b, 1 - b), replace = TRUE)
      }
      
      # Get p and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p <- 
        sum(population$trait == "Monogamy") / N 
    }
  }
  # Export data from function
  output 
}

run1 <- vertical_transmission(1000, 0.5, 0.5, 100, 20)

# Again, the plotting function is adapted from Acerbi et al., mostly for convenience at this point

plot_multiple_runs <- function(data_model) {
  ggplot(data = data_model, aes(y = p, x = generation)) +
    geom_line(aes(colour = run)) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    ylim(c(0, 1)) +
    theme_bw() +
    labs(y = "p (proportion of individuals with Monogamy trait)", x = "Generation")
}

plot_multiple_runs(run1)

# As a simple mechanism, I'll try to add random, unbiased mutation to this - so after individuals inherit the trait from their parent, there is a small chance it will flip (this chance is the mutation rate "mu")

## ABM Version 2 - Vertical transmission, variable mutation rate

vertical_transmission2 <- function(N, mu, p_0, b, t_max, r_max) {
  output <- tibble(generation = rep(1:t_max, r_max), 
                   p = as.numeric(rep(NA, t_max * r_max)), 
                   run = as.factor(rep(1:r_max, each = t_max)))
  
  for (r in 1:r_max) {
    # Create first generation
    population <- tibble(trait = sample(c("Monogamy", "Polygamy"), N, 
                                        replace = TRUE, prob = c(p_0, 1 - p_0)))
    
    # Add first generation's p (proportion of Monogamy trait in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p <- 
      sum(population$trait == "Monogamy") / N 
    
    for (t in 2:t_max) {
      # Copy individuals to previous_population tibble
      previous_population <- population
      
      # Randomly pick mothers and fathers (this will have to be rewritten when sex is added)
      mother <- tibble(trait = sample(previous_population$trait, N, replace = TRUE))  
      father <- tibble(trait = sample(previous_population$trait, N, replace = TRUE)) 
      
      # Prepare next generation
      population <- tibble(trait = as.character(rep(NA, N))) 
      
      # Both parents are monogamous, thus child adopts Monogamy
      both_A <- mother$trait == "Monogamy" & father$trait == "Monogamy"
      if (sum(both_A) > 0) {
        population[both_A, ]$trait <- "Monogamy"  
      }
      
      # Both parents are polygamous, thus child adopts Polygamy
      both_B <- mother$trait == "Polygamy" & father$trait == "Polygamy"
      if (sum(both_B) > 0) {
        population[both_B, ]$trait <- "Polygamy" 
      }
      # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
      if (anyNA(population)) {  
        # They adopt Monogamy with probability b
        population[is.na(population)[,1],]$trait <- 
          sample(c("Monogamy", "Polygamy"), sum(is.na(population)), prob = c(b, 1 - b), replace = TRUE)
      }
      
      # Determine 'mutant' individuals
      mutate <- sample(c(TRUE, FALSE), N, prob = c(mu, 1 - mu), replace = TRUE)
      
      # If there are 'mutants' from Monogamy to Polygamy
      if (nrow(population[mutate & population$trait == "Monogamy", ]) > 0) { 
        # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
        population[mutate & population$trait == "Monogamy", ]$trait <- "Polygamy_m" 
      }
      
      # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
      if (nrow(population[mutate & population$trait == "Polygamy", ]) > 0) { 
        # Then flip them to Monogamy (no need to rename the Monogamy trait here)
        population[mutate & population$trait == "Polygamy", ]$trait <- "Monogamy" 
      }
      
      # To make sure there is no confusion for the next generation, we will now correct the name of the Polygamy trait
      
      population[population$trait == "Polygamy_m", ]$trait <- "Polygamy"
      
      # Get p and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p <- 
        sum(population$trait == "Monogamy") / N 
    }
  }
  # Export data from function
  output 
}

vertical_transmission2(1000,0.05,0.5,0.5,100,10) -> testrun

plot_multiple_runs(testrun)

# The unbiased mutation rate functions as a kind of "stabilising factor", reinforcing an even split of cultural traits

vertical_transmission2(1000,0.05,0.2,0.5,100,10) -> testrun2

plot_multiple_runs(testrun2)

# I Wrote out the function to test it step by step, to check everything is working as intended

population <- tibble(trait = sample(c("Monogamy", "Polygamy"), 1000, 
                                    replace = TRUE, prob = c(0.5, 1 - 0.5)))
previous_population <- population
mother <- tibble(trait = sample(previous_population$trait, 1000, replace = TRUE))  
father <- tibble(trait = sample(previous_population$trait, 1000, replace = TRUE)) 

population <- tibble(trait = as.character(rep(NA, 1000))) 

both_A <- mother$trait == "Monogamy" & father$trait == "Monogamy"
if (sum(both_A) > 0) {
  population[both_A, ]$trait <- "Monogamy"  
}

both_B <- mother$trait == "Polygamy" & father$trait == "Polygamy"
if (sum(both_B) > 0) {
  population[both_B, ]$trait <- "Polygamy" 
}

if (anyNA(population)) {  
  population[is.na(population)[,1],]$trait <- 
    sample(c("Monogamy", "Polygamy"), sum(is.na(population)), prob = c(0.5, 0.5), replace = TRUE)
}

unmutated <- population

mutate <- sample(c(TRUE, FALSE), 1000, prob = c(0.1, 0.9), replace = TRUE)

nrow(unmutated[mutate & unmutated$trait == "Monogamy", ])
nrow(unmutated[mutate & unmutated$trait == "Polygamy", ])

if (nrow(population[mutate & population$trait == "Monogamy", ]) > 0) { 
  population[mutate & population$trait == "Monogamy", ]$trait <- "Polygamy_m" 
}

if (nrow(population[mutate & population$trait == "Polygamy", ]) > 0) { 
  population[mutate & population$trait == "Polygamy", ]$trait <- "Monogamy" 
}

population[population$trait == "Polygamy_m", ]$trait <- "Polygamy"

view(unmutated)
view(population)

# Looks to me as if the population mutated at the level we'd expect. The code seems to be doing what it is supposed to do.

## ABM version 3

# Eventually, we want to simulate the co-evolution of two traits. This will require a slightly modified model

twotrait_vertical_transmission <- function(N, mu_c, mu_b, pc_0, pb_0, b_c, b_b, t_max, r_max) {
  output <- tibble(generation = rep(1:t_max, r_max), 
                   p_c = as.numeric(rep(NA, t_max * r_max)),
                   p_b = as.numeric(rep(NA, t_max * r_max)),
                   run = as.factor(rep(1:r_max, each = t_max)))
  
  for (r in 1:r_max) {
    # Create first generation
    population <- tibble(trait_c = sample(c("Monogamy", "Polygamy"), N, 
                                        replace = TRUE, prob = c(pc_0, 1 - pc_0)), 
                                        trait_b = sample(c("Monogamy", "Polygamy"), N, replace = TRUE, prob = c(pb_0, 1 - pb_0)))
    
    # Add first generation's p_c (proportion of cultural Monogamy trait in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p_c <- 
      sum(population$trait_c == "Monogamy") / N 
    
    # Add first generation's p_b (proportion of biological Monogamy trait in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p_b <- 
      sum(population$trait_b == "Monogamy") / N 
    
    for (t in 2:t_max) {
      # Copy individuals to previous_population tibble
      previous_population <- population
      
      # Randomly pick mothers and fathers (this will have to be rewritten when sex is added)
      mother <- slice_sample(previous_population, n = N, replace = TRUE) 
      father <- slice_sample(previous_population, n = N, replace = TRUE)
      
      # Prepare next generation
      population <- tibble(trait_c = as.character(rep(NA, N)), trait_b = as.character(rep(NA, N))) 
      
      # Both parents are culturally monogamous, thus child adopts cultural Monogamy
      both_mono_c <- mother$trait_c == "Monogamy" & father$trait_c == "Monogamy"
      if (sum(both_mono_c) > 0) {
        population[both_mono_c, ]$trait_c <- "Monogamy"  
      }
      
      # Both parents are culturally polygamous, thus child adopts cultural Polygamy
      both_poly_c <- mother$trait_c == "Polygamy" & father$trait_c == "Polygamy"
      if (sum(both_poly_c) > 0) {
        population[both_poly_c, ]$trait_c <- "Polygamy" 
      }
      # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
      if (anyNA(population$trait_c)) {  
        # They adopt cultural Monogamy with probability b_c
        population[is.na(population$trait_c),]$trait_c <- 
          sample(c("Monogamy", "Polygamy"), sum(is.na(population$trait_c)), prob = c(b_c, 1 - b_c), replace = TRUE)
      }
      
      # Determine 'mutant' individuals
      mutate_c <- sample(c(TRUE, FALSE), N, prob = c(mu_c, 1 - mu_c), replace = TRUE)
      
      # If there are 'mutants' from cultural Monogamy to cultural Polygamy
      if (nrow(population[mutate_c & population$trait_c == "Monogamy", ]) > 0) { 
        # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
        population[mutate_c & population$trait_c == "Monogamy", ]$trait_c <- "Polygamy_m" 
      }
      
      # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
      if (nrow(population[mutate_c & population$trait_c == "Polygamy", ]) > 0) { 
        # Then flip them to Monogamy (no need to rename the Monogamy trait here)
        population[mutate_c & population$trait_c == "Polygamy", ]$trait_c <- "Monogamy" 
      }
      
      # To make sure there is no confusion for the next generation, we will now correct the name of the Polygamy trait
      
      population[population$trait_c == "Polygamy_m", ]$trait_c <- "Polygamy"
      
      # Get p and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p_c <- 
        sum(population$trait_c == "Monogamy") / N 
      
      # We then repeat this process for the biological trait
      
      # Both parents are culturally monogamous, thus child adopts cultural Monogamy
      both_mono_b <- mother$trait_b == "Monogamy" & father$trait_b == "Monogamy"
      if (sum(both_mono_b) > 0) {
        population[both_mono_b, ]$trait_b <- "Monogamy"  
      }
      
      # Both parents are culturally polygamous, thus child adopts cultural Polygamy
      both_poly_b <- mother$trait_b == "Polygamy" & father$trait_b == "Polygamy"
      if (sum(both_poly_b) > 0) {
        population[both_poly_b, ]$trait_b <- "Polygamy" 
      }
      # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
      if (anyNA(population$trait_b)) {  
        # They adopt cultural Monogamy with probability b_b
        population[is.na(population$trait_b),]$trait_b <- 
          sample(c("Monogamy", "Polygamy"), sum(is.na(population$trait_b)), prob = c(b_b, 1 - b_b), replace = TRUE)
      }
      
      # Determine 'mutant' individuals
      mutate_b <- sample(c(TRUE, FALSE), N, prob = c(mu_b, 1 - mu_b), replace = TRUE)
      
      # If there are 'mutants' from cultural Monogamy to cultural Polygamy
      if (nrow(population[mutate_b & population$trait_b == "Monogamy", ]) > 0) { 
        # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
        population[mutate_b & population$trait_b == "Monogamy", ]$trait_b <- "Polygamy_m" 
      }
      
      # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
      if (nrow(population[mutate_b & population$trait_b == "Polygamy", ]) > 0) { 
        # Then flip them to Monogamy (no need to rename the Monogamy trait here)
        population[mutate_b & population$trait_b == "Polygamy", ]$trait_b <- "Monogamy" 
      }
      
      # To make sure there is no confusion for the next generation, we will now correct the name of the Polygamy trait
      
      population[population$trait_b == "Polygamy_m", ]$trait_b <- "Polygamy"
      
      # Get p and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p_b <- 
        sum(population$trait_b == "Monogamy") / N 
      
      
    }
  }
  # Export data from function
  output 
}

plot_2t_multiple_runs <- function(data_model) {
  ggplot(data = data_model, aes(x = generation)) +
    geom_line(aes(y = p_c, colour = run)) +
    geom_line(aes(y = p_b, colour = run), linetype = "dashed") +
    stat_summary(aes(y = p_c), fun = mean, geom = "line", size = 1) +
    stat_summary(aes(y = p_b), fun = mean, geom = "line", size = 1, linetype = "dashed") +
    ylim(c(0, 1)) +
    theme_bw() +
    labs(y = "p (proportion of individuals with Monogamy trait)", x = "Generation")
}

testrun3 <- twotrait_vertical_transmission(1000,0.05,0.05,0.5,0.5,0.5,0.5,10,10)

plot_2t_multiple_runs(testrun3)


## ABM version 4

# The next step is to use the combination of biological and cultural trait to generate a behavioural expression of the underlying phenogenotype

twotrait_vertical_transmission2 <- function(N, mu_c, mu_b, pc_0, pb_0, b_c, b_b, b_bh, t_max, r_max) {
  output <- tibble(generation = rep(1:t_max, r_max), 
                   p_c = as.numeric(rep(NA, t_max * r_max)),
                   p_b = as.numeric(rep(NA, t_max * r_max)),
                   p_bh = as.numeric(rep(NA, t_max * r_max)),
                   run = as.factor(rep(1:r_max, each = t_max)))
  
  for (r in 1:r_max) {
    # Create first generation
    population <- tibble(trait_c = sample(c("Monogamy", "Polygamy"), N, 
                                          replace = TRUE, prob = c(pc_0, 1 - pc_0)), 
                         trait_b = sample(c("Monogamy", "Polygamy"), N, replace = TRUE, prob = c(pb_0, 1 - pb_0)),
                         bh = as.character(rep(NA, N)))
   
    # Determining behavioural tendencies based on cultural and biological trait combination
    
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
      # They will adopt the behavioural variant corresponding to their cultural trait with probability b_bh
      population[is.na(population$bh)& population$trait_c == "Monogamy",]$bh <- 
        sample(c("Monogamy", "Polygamy"), sum(is.na(population$bh)&population$trait_c == "Monogamy"), prob = c(b_bh, 1 - b_bh), replace = TRUE)
      population[is.na(population$bh)& population$trait_c == "Polygamy",]$bh <- 
        sample(c("Polygamy", "Monogamy"), sum(is.na(population$bh)&population$trait_c == "Polygamy"), prob = c(b_bh, 1 - b_bh), replace = TRUE)
    }
    
    # Add first generation's p_c (proportion of cultural Monogamy trait in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p_c <- 
      sum(population$trait_c == "Monogamy") / N 
    
    # Add first generation's p_b (proportion of biological Monogamy trait in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p_b <- 
      sum(population$trait_b == "Monogamy") / N 
    
    # Add first generation's p_bh (proportion of Monogamous behaviour in entire population) for run r
    output[output$generation == 1 & output$run == r, ]$p_bh <- 
      sum(population$bh == "Monogamy") / N 
    
    for (t in 2:t_max) {
      # Copy individuals to previous_population tibble
      previous_population <- population
      
      # Randomly pick mothers and fathers (this will have to be rewritten when sex is added)
      mother <- slice_sample(previous_population, n = N, replace = TRUE) 
      father <- slice_sample(previous_population, n = N, replace = TRUE)
      
      # Prepare next generation
      population <- tibble(trait_c = as.character(rep(NA, N)), trait_b = as.character(rep(NA, N)), bh = as.character(rep(NA, N))) 
      
      # Both parents are culturally monogamous, thus child adopts cultural Monogamy
      both_mono_c <- mother$trait_c == "Monogamy" & father$trait_c == "Monogamy"
      if (sum(both_mono_c) > 0) {
        population[both_mono_c, ]$trait_c <- "Monogamy"  
      }
      
      # Both parents are culturally polygamous, thus child adopts cultural Polygamy
      both_poly_c <- mother$trait_c == "Polygamy" & father$trait_c == "Polygamy"
      if (sum(both_poly_c) > 0) {
        population[both_poly_c, ]$trait_c <- "Polygamy" 
      }
      # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
      if (anyNA(population$trait_c)) {  
        # They adopt cultural Monogamy with probability b_c
        population[is.na(population$trait_c),]$trait_c <- 
          sample(c("Monogamy", "Polygamy"), sum(is.na(population$trait_c)), prob = c(b_c, 1 - b_c), replace = TRUE)
      }
      
      # Determine 'mutant' individuals (for cultural traits, this is roughly analogous to "individual learners")
      mutate_c <- sample(c(TRUE, FALSE), N, prob = c(mu_c, 1 - mu_c), replace = TRUE)
      
      # If there are 'mutants' from cultural Monogamy to cultural Polygamy
      if (nrow(population[mutate_c & population$trait_c == "Monogamy", ]) > 0) { 
        # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
        population[mutate_c & population$trait_c == "Monogamy", ]$trait_c <- "Polygamy_m" 
      }
      
      # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
      if (nrow(population[mutate_c & population$trait_c == "Polygamy", ]) > 0) { 
        # Then flip them to Monogamy (no need to rename the Monogamy trait here)
        population[mutate_c & population$trait_c == "Polygamy", ]$trait_c <- "Monogamy" 
      }
      
      # To make sure there is no confusion for the next generation, we will now correct the name of the Polygamy trait
      
      population[population$trait_c == "Polygamy_m", ]$trait_c <- "Polygamy"
      
      # Get p and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p_c <- 
        sum(population$trait_c == "Monogamy") / N 
      
      # We then repeat this process for the biological trait
      
      # Both parents are biologically monogamous, thus child adopts biological Monogamy
      both_mono_b <- mother$trait_b == "Monogamy" & father$trait_b == "Monogamy"
      if (sum(both_mono_b) > 0) {
        population[both_mono_b, ]$trait_b <- "Monogamy"  
      }
      
      # Both parents are biologically polygamous, thus child adopts biological Polygamy
      both_poly_b <- mother$trait_b == "Polygamy" & father$trait_b == "Polygamy"
      if (sum(both_poly_b) > 0) {
        population[both_poly_b, ]$trait_b <- "Polygamy" 
      }
      # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
      if (anyNA(population$trait_b)) {  
        # They inherit biological Monogamy with probability b_b
        population[is.na(population$trait_b),]$trait_b <- 
          sample(c("Monogamy", "Polygamy"), sum(is.na(population$trait_b)), prob = c(b_b, 1 - b_b), replace = TRUE)
      }
      
      # Determine 'mutant' individuals (for the biological trait, this is closer to genetic mutation). NB this mutation does not generate novel traits outside our poly/mono dichotimy
      mutate_b <- sample(c(TRUE, FALSE), N, prob = c(mu_b, 1 - mu_b), replace = TRUE)
      
      # If there are 'mutants' from biological Monogamy to biological Polygamy
      if (nrow(population[mutate_b & population$trait_b == "Monogamy", ]) > 0) { 
        # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
        population[mutate_b & population$trait_b == "Monogamy", ]$trait_b <- "Polygamy_m" 
      }
      
      # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
      if (nrow(population[mutate_b & population$trait_b == "Polygamy", ]) > 0) { 
        # Then flip them to Monogamy (no need to rename the Monogamy trait here)
        population[mutate_b & population$trait_b == "Polygamy", ]$trait_b <- "Monogamy" 
      }
      
      # To make sure there is no confusion for the next generation, we will now correct the name of the Polygamy trait
      
      population[population$trait_b == "Polygamy_m", ]$trait_b <- "Polygamy"
      
      # Get p and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p_b <- 
        sum(population$trait_b == "Monogamy") / N 
      
      # Finally we determine the behavioural variant resulting from the biological/cultural trait combination as we did in Gen 1
      
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
        # They will adopt the behavioural variant corresponding to their cultural trait with probability b_bh
        population[is.na(population$bh)& population$trait_c == "Monogamy",]$bh <- 
          sample(c("Monogamy", "Polygamy"), sum(is.na(population$bh)&population$trait_c == "Monogamy"), prob = c(b_bh, 1 - b_bh), replace = TRUE)
        population[is.na(population$bh)& population$trait_c == "Polygamy",]$bh <- 
          sample(c("Polygamy", "Monogamy"), sum(is.na(population$bh)&population$trait_c == "Polygamy"), prob = c(b_bh, 1 - b_bh), replace = TRUE)
      }
      
      # Get p and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p_bh <- 
        sum(population$bh == "Monogamy") / N
      
      
    }
  }
  # Export data from function
  output 
}


testrun4 <- twotrait_vertical_transmission2(1000, 0.05, 0.05, 0.5, 0.5, 0.5, 0.5, 0.5, 100, 10)

plot_multiple_runs_bh <- function(data_model) {
  ggplot(data = data_model, aes(y = p_bh, x = generation)) +
    geom_line(aes(colour = run)) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    ylim(c(0, 1)) +
    theme_bw() +
    labs(y = "p (proportion of individuals with Monogamy trait)", x = "Generation")
}

plot_multiple_runs_bh(testrun4)

testrun5 <- twotrait_vertical_transmission2(1000, 0.05, 0.05, 0.1, 0.1, 0.51, 0.5, 0.7, 500, 10)

plot_multiple_runs_bh(testrun5)


