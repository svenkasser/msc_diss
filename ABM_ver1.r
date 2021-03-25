library(tidyverse)

# First, I'm gonna follow the Acerbi book closely to build a simple vertical transmission model of Marriage system preference as a cultural belief

# This initial model is not going to have biological sex, or mating, and only has a features a singular cultural trait

# Only allows for a handful of parameters to vary - pop size, proportion at gen 1, likelihood of Monogamy trait in the case of mixed inheritance (I could fix this at 0.5), max number of generations simulated, and max number of populations simulated

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

## Again, the plotting function is adapted from Acerbi et al., mostly for convenience at this point

plot_multiple_runs <- function(data_model) {
  ggplot(data = data_model, aes(y = p, x = generation)) +
    geom_line(aes(colour = run)) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    ylim(c(0, 1)) +
    theme_bw() +
    labs(y = "p (proportion of individuals with Monogamy trait)", x = "Generation")
}

plot_multiple_runs(run1)
