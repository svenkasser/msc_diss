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

## Looks to me as if the population mutated at the level we'd expect. The code seems to be doing what it is supposed to do.