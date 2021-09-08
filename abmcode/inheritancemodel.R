## Load required packages for data handling

library(tidyverse)
library(dplyr)

## Custom function to determine traits for the next gen

trait_determination <- function(input1, input2, b_c, b_b, cb_bh) {
  
  population <- input1
  nextgen <- input2
  
  
  # Determining cultural traits
  
  # Both traits favour monogamy, thus cultural preference is monogamous
  full_mono <- nextgen$trait_c_p1 == "Monogamy" & nextgen$trait_c_p2 == "Monogamy"
  if (sum(full_mono) > 0) {
    population[full_mono, ]$trait_c <- "Monogamy"  
  }
  
  # Both traits favour polygamy, thus cultural preference is polygamous
  full_poly <- nextgen$trait_c_p1 == "Polygamy" & nextgen$trait_c_p2 == "Polygamy"
  if (sum(full_poly) > 0) {
    population[full_poly, ]$trait_c <- "Polygamy" 
  }
  
  # If they have mixed parental cultural traits (i.e. one Monogamy and one Polygamy parent)
  if (anyNA(population$trait_c)) {  
    # They adopt cultural Monogamy with probability b_c (parameter allows for a kind of dominance/recessiveness)
    population[is.na(population$trait_c),]$trait_c <- 
      sample(c("Monogamy", "Polygamy"), sum(is.na(population$trait_c)), prob = c(b_c, 1 - b_c), replace = TRUE)
  }
  
  
  # Determining biological traits
  
  # Both traits favour monogamy, thus cultural preference is monogamous
  full_mono <- nextgen$trait_b_p1 == "Monogamy" & nextgen$trait_b_p2 == "Monogamy"
  if (sum(full_mono) > 0) {
    population[full_mono, ]$trait_b <- "Monogamy"  
  }
  
  # Both traits favour polygamy, thus cultural preference is polygamous
  full_poly <- nextgen$trait_b_p1 == "Polygamy" & nextgen$trait_b_p2 == "Polygamy"
  if (sum(full_poly) > 0) {
    population[full_poly, ]$trait_b <- "Polygamy" 
  }
  
  # If they have mixed parental biological traits (i.e. one Monogamy and one Polygamy parent)
  if (anyNA(population$trait_b)) {  
    # They adopt cultural Monogamy with probability b_c (parameter allows for a kind of dominance/recessiveness)
    population[is.na(population$trait_b),]$trait_b <- 
      sample(c("Monogamy", "Polygamy"), sum(is.na(population$trait_b)), prob = c(b_b, 1 - b_b), replace = TRUE)
  }
  
  # Determining behavioural tendencies based on cultural and biological trait combinations
  
  # Both traits favour monogamy, thus behaviour is monogamous
  full_mono <- population$trait_c == "Monogamy" & population$trait_b == "Monogamy"
  if (sum(full_mono) > 0) {
    population[full_mono, ]$bh <- "Monogamy"  
  }
  
  # Both traits favour polygamy, thus cultural preference is polygamous
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
  
  return(population)
  
}
