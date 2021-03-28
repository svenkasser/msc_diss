### WIP - Mating model

library(tidyverse)
library(dplyr)
library(ggplot2)

# Create first generation
population <- tibble(sex = sample(c("Female", "Male"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)),
                     trait_c = sample(c("Monogamy", "Polygamy"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)), 
                     trait_b = sample(c("Monogamy", "Polygamy"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)),
                     bh = as.character(rep(NA, 100)))

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
    sample(c("Monogamy", "Polygamy"), sum(is.na(population$bh)&population$trait_c == "Monogamy"), prob = c(0.5, 1 - 0.5), replace = TRUE)
  population[is.na(population$bh)& population$trait_c == "Polygamy",]$bh <- 
    sample(c("Polygamy", "Monogamy"), sum(is.na(population$bh)&population$trait_c == "Polygamy"), prob = c(0.5, 1 - 0.5), replace = TRUE)
}

# Create mating pool

matingpool <- population %>%
  add_column(mated = "No", offspring = 0, ID = 1:100)

# Create offspring tibble

nextgen <- tibble(sex = sample(c("Female", "Male"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)),
                     trait_c = as.character(rep(NA, 100)), 
                     trait_b = as.character(rep(NA, 100)),
                     bh = as.character(rep(NA, 100)))

for (i in 1:100) {
  
  mate1 <- slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$mated=="No")|(matingpool$bh=="Polygamy"),], n=1)
  
  if (mate1$sex=="Female") {
    mate2 <- 
      slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$mated=="No" & matingpool$sex=="Male")|(matingpool$bh=="Polygamy" & matingpool$sex=="Male"),], n=1)
  }
  if (mate1$sex=="Male") {
    mate2 <- 
      slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$mated=="No" & matingpool$sex=="Female")|(matingpool$bh=="Polygamy" & matingpool$sex=="Female"),], n=1)
  }
  
  matingpool[mate1$ID,]$mated <- "Yes"
  matingpool[mate2$ID,]$mated <- "Yes"
  
  matingpool[mate1$ID,]$offspring <- matingpool[mate1$ID,]$offspring + 1
  matingpool[mate2$ID,]$offspring <- matingpool[mate2$ID,]$offspring + 1
  
  # Both parents are culturally monogamous, thus child adopts cultural Monogamy
  if (mate1$trait_c == "Monogamy" & mate2$trait_c == "Monogamy") {
    nextgen$trait_c[i] <- "Monogamy"  
  }
  
  # Both parents are culturally polygamous, thus child adopts cultural Polygamy
  if (mate1$trait_c == "Polygamy" & mate2$trait_c == "Polygamy") {
    nextgen$trait_c[i] <- "Polygamy" 
  }
  
  # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
  if (is.na(nextgen$trait_c[i])) {  
    # They adopt cultural Monogamy with probability b_c
    nextgen$trait_c[i] <- 
      sample(c("Monogamy", "Polygamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  }
  
  # Determine 'mutant' individuals (for cultural traits, this is roughly analogous to "individual learners")
  mutate_c <- sample(c(TRUE, FALSE), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  
  # If there are 'mutants' from cultural Monogamy to cultural Polygamy
  if (mutate_c & nextgen$trait_c[i] == "Monogamy") { 
    # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
    nextgen$trait_c[i] <- "Polygamy_m" 
  }
  
  # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
  if (mutate_c & nextgen$trait_c[i] == "Polygamy") { 
    # Then flip them to Monogamy (no need to rename the Monogamy trait here)
    nextgen$trait_c[i] <- "Monogamy" 
  }
  
  if (nextgen$trait_c[i]=="Polygamy_m"){
    nextgen$trait_c[i] <- "Polygamy"
  }
  
  ## Rinse and repeat with biological traits
  
  # Both parents are biologically monogamous, thus child inherits biological Monogamy
  if (mate1$trait_b == "Monogamy" & mate2$trait_b == "Monogamy") {
    nextgen$trait_b[i] <- "Monogamy"  
  }
  
  # Both parents are culturally polygamous, thus child inherits biological Polygamy
  if (mate1$trait_b == "Polygamy" & mate2$trait_b == "Polygamy") {
    nextgen$trait_b[i] <- "Polygamy" 
  }
  
  # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
  if (is.na(nextgen$trait_b[i])) {  
    # They adopt biological Monogamy with probability b_b
    nextgen$trait_b[i] <- 
      sample(c("Monogamy", "Polygamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  }
  
  # Determine 'mutant' individuals (for cultural traits, this is roughly analogous to "individual learners")
  mutate_b <- sample(c(TRUE, FALSE), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  
  # If there are 'mutants' from cultural Monogamy to cultural Polygamy
  if (mutate_b & nextgen$trait_b[i] == "Monogamy") { 
    # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
    nextgen$trait_b[i] <- "Polygamy_m" 
  }
  
  # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
  if (mutate_b & nextgen$trait_b[i] == "Polygamy") { 
    # Then flip them to Monogamy (no need to rename the Monogamy trait here)
    nextgen$trait_b[i] <- "Monogamy" 
  }
  
  if (nextgen$trait_b[i]=="Polygamy_m"){
    nextgen$trait_b[i] <- "Polygamy"
  }
  
  # Determining behavioural tendencies based on cultural and biological trait combination
  
  # Both traits favour monogamy, thus behavioural preference is monogamous
  if (nextgen$trait_c[i] == "Monogamy" & nextgen$trait_b[i] == "Monogamy") {
    nextgen$bh[i] <- "Monogamy"  
  }
  
  # Both traits favour polygamy, thus behavioural preference is polygamous
  if (nextgen$trait_c[i] == "Polygamy" & nextgen$trait_b[i] == "Polygamy") {
    population$bh[i] <- "Polygamy" 
  }
  
  # If any empty NA slots (i.e. mixed cultural/biological traits) are present
  # They will adopt the behavioural variant corresponding to their cultural trait with probability b_bh
  if (is.na(nextgen$bh[i]) & nextgen$trait_c[i] == "Monogamy") {  
    nextgen$bh[i] <- 
      sample(c("Monogamy", "Polygamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  } else {
    nextgen$bh[i] <- 
      sample(c("Polygamy", "Monogamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  }
  

  
}

######################################## Test Space ####################################################

i <- 1

nextgen <- tibble(sex = sample(c("Female", "Male"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)),
                  trait_c = as.character(rep(NA, 100)), 
                  trait_b = as.character(rep(NA, 100)),
                  bh = as.character(rep(NA, 100)))
  
  mate1 <- slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$mated=="No")|(matingpool$bh=="Polygamy"),], n=1)
  
  if (mate1$sex=="Female") {
    mate2 <- 
      slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$mated=="No" & matingpool$sex=="Male")|(matingpool$bh=="Polygamy" & matingpool$sex=="Male"),], n=1)
  }
  if (mate1$sex=="Male") {
    mate2 <- 
      slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$mated=="No" & matingpool$sex=="Female")|(matingpool$bh=="Polygamy" & matingpool$sex=="Female"),], n=1)
  }
  
  matingpool[mate1$ID,]$mated <- "Yes"
  matingpool[mate2$ID,]$mated <- "Yes"
  
  matingpool[mate1$ID,]$offspring <- matingpool[mate1$ID,]$offspring + 1
  matingpool[mate2$ID,]$offspring <- matingpool[mate2$ID,]$offspring + 1
  
  # Both parents are culturally monogamous, thus child adopts cultural Monogamy
  if (mate1$trait_c == "Monogamy" & mate2$trait_c == "Monogamy") {
    nextgen$trait_c[i] <- "Monogamy"  
  }
  
  # Both parents are culturally polygamous, thus child adopts cultural Polygamy
  if (mate1$trait_c == "Polygamy" & mate2$trait_c == "Polygamy") {
    nextgen$trait_c[i] <- "Polygamy" 
  }
  
  # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
  if (is.na(nextgen$trait_c[i])) {  
    # They adopt cultural Monogamy with probability b_c
    nextgen$trait_c[i] <- 
      sample(c("Monogamy", "Polygamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  }
  
  # Determine 'mutant' individuals (for cultural traits, this is roughly analogous to "individual learners")
  mutate_c <- sample(c(TRUE, FALSE), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  
  # If there are 'mutants' from cultural Monogamy to cultural Polygamy
  if (mutate_c & nextgen$trait_c[i] == "Monogamy") { 
    # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
    nextgen$trait_c[i] <- "Polygamy_m" 
  }
  
  # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
  if (mutate_c & nextgen$trait_c[i] == "Polygamy") { 
    # Then flip them to Monogamy (no need to rename the Monogamy trait here)
    nextgen$trait_c[i] <- "Monogamy" 
  }
  
  if (nextgen$trait_c[i]=="Polygamy_m"){
    nextgen$trait_c[i] <- "Polygamy"
  }
  
  ## Rinse and repeat with biological traits
  
  # Both parents are biologically monogamous, thus child inherits biological Monogamy
  if (mate1$trait_b == "Monogamy" & mate2$trait_b == "Monogamy") {
    nextgen$trait_b[i] <- "Monogamy"  
  }
  
  # Both parents are culturally polygamous, thus child inherits biological Polygamy
  if (mate1$trait_b == "Polygamy" & mate2$trait_b == "Polygamy") {
    nextgen$trait_b[i] <- "Polygamy" 
  }
  
  # If any empty NA slots (i.e. one Monogamy and one Polygamy parent) are present
  if (is.na(nextgen$trait_b[i])) {  
    # They adopt biological Monogamy with probability b_b
    nextgen$trait_b[i] <- 
      sample(c("Monogamy", "Polygamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  }
  
  # Determine 'mutant' individuals (for cultural traits, this is roughly analogous to "individual learners")
  mutate_b <- sample(c(TRUE, FALSE), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  
  # If there are 'mutants' from cultural Monogamy to cultural Polygamy
  if (mutate_b & nextgen$trait_b[i] == "Monogamy") { 
    # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
    nextgen$trait_b[i] <- "Polygamy_m" 
  }
  
  # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
  if (mutate_b & nextgen$trait_b[i] == "Polygamy") { 
    # Then flip them to Monogamy (no need to rename the Monogamy trait here)
    nextgen$trait_b[i] <- "Monogamy" 
  }
  
  if (nextgen$trait_b[i]=="Polygamy_m"){
    nextgen$trait_b[i] <- "Polygamy"
  }
  
  # Determining behavioural tendencies based on cultural and biological trait combination
  
  # Both traits favour monogamy, thus behavioural preference is monogamous
  if (nextgen$trait_c[i] == "Monogamy" & nextgen$trait_b[i] == "Monogamy") {
    nextgen$bh[i] <- "Monogamy"  
  }
  
  # Both traits favour polygamy, thus behavioural preference is polygamous
  if (nextgen$trait_c[i] == "Polygamy" & nextgen$trait_b[i] == "Polygamy") {
    population$bh[i] <- "Polygamy" 
  }
  
  # If any empty NA slots (i.e. mixed cultural/biological traits) are present
  # They will adopt the behavioural variant corresponding to their cultural trait with probability b_bh
  if (is.na(nextgen$bh[i]) & nextgen$trait_c[i] == "Monogamy") {  
    nextgen$bh[i] <- 
      sample(c("Monogamy", "Polygamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  } else {
    nextgen$bh[i] <- 
      sample(c("Polygamy", "Monogamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  }
