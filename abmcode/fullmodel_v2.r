library(tidyverse)
library(dplyr)
library(ggplot2)

#### First complete model

### Quick definitions for all parameters

## Simulation parameters

# r_max = number of separate model populations simulated
# Possible values: 0 < r_max
# Default value: N/A
# t_max = number of distinct, discrete generations simulated for each population
# Possible values: 0 < t_max
# Default value: N/A

## Demographic parameters

# N = population size. Population size is static for the time being.
# Possible values: 0 < N
# Default value: N/A
# pf_0 = avg. percentage of female individuals >in the starting population<. This is our starting sex-ratio stand-in.
# Possible values: 0 =< pf_0 =< 1
# Default value: 0.5 - equal starting sex-ratio
# b_f = female demographic bias in subsequent generations. This can be used to model stable sex ratio biases in the adult population (e.g. sex differences in infant mortality)
# Possible values: 0 =< b_f =< 1
# Default value: 0.5 - no non-randomly biased sex-ratios

## Transmission parameters

# pc_0 = avg. starting trait frequency for cultural "Monogamy"
# Possible values: 0 =< pc_0 =< 1
# Default value: 0.5 - equal starting cultural Poly/Mono split
# pb_0 = avg. starting trait frequency for biological "Monogamy"
# Possible values: 0 =< pb_0 =< 1
# Default value: 0.5 - equal starting biological Poly/Mono split
# b_c = transmission bias for cultural "Monogamy". This only plays a role for individuals with mixed parental traits.
# Possible values: 0 =< b_c =< 1
# Default value: 0.5 - no biased transmission. Higher/lower values increase/decrease the likelihood of adopting "Monogamy".
# b_b = transmission bias for biological "Monogamy". This only plays a role for individuals that inherit mixed parental traits.
# Possible values: 0 =< b_b =< 1
# Default value: 0.5 - no biased transmission. Higher/lower values increase/decrease the likelihood of inheriting "Monogamy".
# cb_bh = culture bias in behaviour. This only plays a role for individuals with mixed cultural / biological traits.
# Possible values: 0 =< b_bh =< 1
# Default value: 0.5 - no bias toward cultural trait in behaviour. Higher/lower values increase/decrease the likelihood of adopting the cultural trait as the behavioural phenotype.

## Mutation parameters
# NB. Mutation in this model cannot generate "new" trait variants, but only flip the existing ones within the Mono/Poly dichotomy

# mu_c = mutation rate of cultural traits. This is roughly analogous to the frequency of individual learning "overriding" a vertically inherited trait.
# Possible values: 0 =< mu_c =< 1
# Default value: 0 - no mutation
# mu_b = mutation rate of biological traits. This is roughly analogous to genetic mutation.
# Possible values: 0 =< mu_c =< 1
# Default value: 0 - no mutation

## Reproduction parameters

# mono_max = maximum number of offspring for one monogamous pair
# Possible values: 0 < mono_max
# Default value: 100000 - effectively no reproductive limits on monogamous pairs
# poly_max = maximum number of offspring for a polygamous individual
# Possible values: 0 < poly_max
# Default value: 100000 - effectively no reproductive limits on monogamous pairs

## The explanations about the inner workings of this model are laid out with comments throughout the code below
# N.B. This is the first try combining the transmission models from Acerbi et al. with a mating system.
# As such, this model is currently painfully slow and in dire need of optimisation.
# The next step would probably be to replace the current mating system loop which generates offspring one by one with a noon-looping system


ABMmodel2 <- function(N, t_max, r_max, mu_c = 0, mu_b = 0, pc_0 = 0.5, pb_0 = 0.5, b_c = 0.5, b_b = 0.5, cb_bh = 0.5, pf_0 = 0.5, b_f = 0.5, mono_max = 100000, poly_max = 100000) {
  
  # Create output file - this is he data what we want to save from each run and generation
  
  output <- tibble(generation = rep(1:t_max, r_max), 
                   p_c = as.numeric(rep(NA, t_max * r_max)),
                   p_b = as.numeric(rep(NA, t_max * r_max)),
                   p_bh = as.numeric(rep(NA, t_max * r_max)),
                   p_f = as.numeric(rep(NA, t_max * r_max)),
                   p_mono = as.numeric(rep(NA, t_max * r_max)),
                   run = as.factor(rep(1:r_max, each = t_max)))
  
  # Loop for each run (essentially each simulated population)
  
  for (r in 1:r_max) {
    
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
      
      # Create mating pool tibble
      # The mating pool tibble stores a bit more information that the normal population table, incl. mating status, number of offspring
      # Each individual is also given a unique ID so that each successful mating from the loop below can be recorded in the mating pool tibble
      
      matingpool <- previous_population %>%
        add_column(bonded = "No", bond_ID = as.character("NA"), offspring = 0, ID = 1:N)
      
      # Create empty offspring tibble
      # This is the empty table of offspring that will constitute the next generation once filled by the mating loop. 
      # Without population growth, the next generation will consist of N individuals
      
      nextgen <- tibble(sex = sample(c("Female", "Male"), N, replace = TRUE, prob = c(b_f, 1 - b_f)),
                        trait_c = as.character(rep(NA, N)), 
                        trait_b = as.character(rep(NA, N)),
                        bh = as.character(rep(NA, N)))
      
      # The next loop is the mating process script that I wrote separately (this is not from the Acerbi book, though it used concepts an structures introduced there)
      
      # Each iteration of the loop is a single mating encounter with a 100% chance of reproduction. This is an assumption that could eventually be relaxed.
      
      # Mating is random - the only rule is that monogamous individuals mate only with one partner. Polygamous individuals are not limited in their mating opportunities, except that they can't mate with "bonded" monogamists
      
      # Because I wrote this myself, this is likelihood the least optimised bit of this entire script for now.
      
      for (i in 1:N) {
        
        # Randomly draw the first participant of this mating encounter from the mating pool
        # NB. Has to be Monogamous with less than X (mono_max) offspring or Polygamous with less than X (poly_max) offspring.
        # Setting mono_max and poly_max allows us to set hard constraints on the number of offspring a monogamous pair or polygamous parent can have
        # This can approximate time or resource constraints in the absence of modelling those directly
        
        mate1 <- slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$offspring < mono_max)|(matingpool$bh=="Polygamy" & matingpool$offspring < poly_max),], n=1)
        
        # Based on the sex of the first participant, an opposite sex partner is drawn from the same population (same rules)
        
        if (mate1$sex=="Female" & mate1$bonded=="No") {
          mate2 <- 
            slice_sample(matingpool[matingpool$sex=="Male" & matingpool$bonded == "No" & matingpool$offspring < poly_max & matingpool$ID != mate1$ID,], n=1)
        } 
        
        if (mate1$sex=="Male" & mate1$bonded=="No") {
          mate2 <- 
            slice_sample(matingpool[matingpool$sex=="Female" & matingpool$bonded == "No" & matingpool$offspring < poly_max & matingpool$ID != mate1$ID,], n=1)
        } 
        
        # If the first participant of this mating encounter is already monogamously bonded (from a previous encounter), then Mate 2 will simply be his bondmate from previous mating encounters
        
        if (mate1$bonded=="Yes") {
          mate2 <- 
            matingpool[mate1$bond_ID,]
        }
        
        # Record mating encounter in matingpool tibble
        
        if (mate1$bh=="Monogamy"){
          matingpool[mate1$ID,]$bonded <- "Yes"
          matingpool[mate1$ID,]$bond_ID <- as.character(mate2$ID)
        }
        
        if (mate2$bh=="Monogamy") {
          matingpool[mate2$ID,]$bonded <- "Yes"
          matingpool[mate2$ID,]$bond_ID <- as.character(mate1$ID)
        }
        
        matingpool[mate1$ID,]$offspring <- matingpool[mate1$ID,]$offspring + 1
        matingpool[mate2$ID,]$offspring <- matingpool[mate2$ID,]$offspring + 1
        
        ## Next, we start to populate the data row for the offspring that this encounter has produced, starting with the cultural trait
        
        # Both parents are culturally monogamous, thus child adopts cultural Monogamy
        if (mate1$trait_c == "Monogamy" & mate2$trait_c == "Monogamy") {
          nextgen$trait_c[i] <- "Monogamy"  
        }
        
        # Both parents are culturally polygamous, thus child adopts cultural Polygamy
        if (mate1$trait_c == "Polygamy" & mate2$trait_c == "Polygamy") {
          nextgen$trait_c[i] <- "Polygamy" 
        }
        
        # If they have mixed parental cultural traits (i.e. one Monogamy and one Polygamy parent)
        if (is.na(nextgen$trait_c[i])) {  
          # They adopt cultural Monogamy with probability b_c (parameter allows for a kind of dominance/recessiveness)
          nextgen$trait_c[i] <- 
            sample(c("Monogamy", "Polygamy"), 1, prob = c(b_c, 1 - b_c), replace = TRUE)
        }
        
        # Determine if this is a 'mutant' individual (for cultural traits, this is roughly analogous to "individual learners")
        mutate_c <- sample(c(TRUE, FALSE), 1, prob = c(mu_c, 1 - mu_c), replace = TRUE)
        
        # If they are a 'mutant' from cultural Monogamy to cultural Polygamy
        if (mutate_c & nextgen$trait_c[i] == "Monogamy") { 
          # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
          nextgen$trait_c[i] <- "Polygamy_m" 
        }
        
        # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
        if (mutate_c & nextgen$trait_c[i] == "Polygamy") { 
          # Then flip them to Monogamy (no need to rename the Monogamy trait here)
          nextgen$trait_c[i] <- "Monogamy" 
        }
        
        # Rename to standardize (note: I'll have to check the code, but I think these problems with Polygamy_m, etc., are an artifact of previous problems in non-mating transmission models - I'll have to test whether they are still strictly necessar here)
        
        if (nextgen$trait_c[i]=="Polygamy_m"){
          nextgen$trait_c[i] <- "Polygamy"
        }
        
        ## We repeat this with the biological traits
        
        # Both parents are biologically monogamous, thus child inherits biological Monogamy
        if (mate1$trait_b == "Monogamy" & mate2$trait_b == "Monogamy") {
          nextgen$trait_b[i] <- "Monogamy"  
        }
        
        # Both parents are culturally polygamous, thus child inherits biological Polygamy
        if (mate1$trait_b == "Polygamy" & mate2$trait_b == "Polygamy") {
          nextgen$trait_b[i] <- "Polygamy" 
        }
        
        # If they have mixed parental biological traits (i.e. one Monogamy and one Polygamy parent)
        if (is.na(nextgen$trait_b[i])) {  
          # They adopt biological Monogamy with probability b_b
          nextgen$trait_b[i] <- 
            sample(c("Monogamy", "Polygamy"), 1, prob = c(b_b, 1 - b_b), replace = TRUE)
        }
        
        # Determine 'mutant' individuals (for biological traits, this is roughly analogous to random genetic mutation, but is constrained to flipping a trait rather than generating new variants)
        mutate_b <- sample(c(TRUE, FALSE), 1, prob = c(mu_b, 1 - mu_b), replace = TRUE)
        
        # If there are 'mutants' from biological Monogamy to biological Polygamy
        if (mutate_b & nextgen$trait_b[i] == "Monogamy") { 
          # Then flip them to Polygamy (for the time being, we call this mutation-induced Polygamy "Polygamy_m", because otherwise the next operation is going to reverse this entire step)
          nextgen$trait_b[i] <- "Polygamy_m" 
        }
        
        # If there are 'mutants' from Polygamy to Monogamy (note: if the mutation induced Polygamy trait was simply called "Polygamy" here, this would reverse the previous step as this operation selects all "Polygamy" labelled traits, including mutated ones)
        if (mutate_b & nextgen$trait_b[i] == "Polygamy") { 
          # Then flip them to Monogamy (no need to rename the Monogamy trait here)
          nextgen$trait_b[i] <- "Monogamy" 
        }
        
        # Rename to standardize (note: as above)
        
        if (nextgen$trait_b[i]=="Polygamy_m"){
          nextgen$trait_b[i] <- "Polygamy"
        }
        
        ## Determining behavioural tendencies based on cultural and biological trait combination
        
        # Both traits favour monogamy, thus behavioural preference is monogamous
        if (nextgen$trait_c[i] == "Monogamy" & nextgen$trait_b[i] == "Monogamy") {
          nextgen$bh[i] <- "Monogamy"  
        }
        
        # Both traits favour polygamy, thus behavioural preference is polygamous
        if (nextgen$trait_c[i] == "Polygamy" & nextgen$trait_b[i] == "Polygamy") {
          population$bh[i] <- "Polygamy" 
        }
        
        # If the individual has mixed cultural/biological traits
        if (is.na(nextgen$bh[i]) & nextgen$trait_c[i] == "Monogamy") {  
          # They will adopt the behavioural variant corresponding to their cultural trait with probability cb_bh
          nextgen$bh[i] <- 
            sample(c("Monogamy", "Polygamy"), 1, prob = c(cb_bh, 1 - cb_bh), replace = TRUE)
        } else {
          nextgen$bh[i] <- 
            sample(c("Polygamy", "Monogamy"), 1, prob = c(cb_bh, 1 - cb_bh), replace = TRUE)
        }
        
        
        
      }
      
      # Now that mating is complete we have a new population of individuals, from which we can extract the trait frequencies from time t in r
      
      population <- nextgen
      
      # Get p_c (cultural trait frequency for "Monogamy") and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p_c <- 
        sum(population$trait_c == "Monogamy") / N 
      
      # Get p_b (Biological trait frequency for "Monogamy") and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p_b <- 
        sum(population$trait_b == "Monogamy") / N 
      
      # Get p_bh (Behavioural trait frequency for "Monogamy") and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p_bh <- 
        sum(population$bh == "Monogamy") / N 
      
      # Get the sex ratio p_f (as proportion of female individuals) for this generation t in run r
      output[output$generation == t & output$run == r, ]$p_f <- 
        sum(population$sex == "Female") / N 
      
      # Get observable monogamy outcome from mating pool
      output[output$generation == t-1 & output$run == r, ]$p_mono <-
        sum(matingpool$bonded == "Yes") / N
      
      
    }
    
  }
  
  # Export data from function
  output
  
}




tstrun <- ABMmodel2(100, 50, 10)


## Some plots and exploratory analysis

ggplot(data = tstrun, aes(y = p_bh, x = generation)) +
  geom_line(aes(colour = run)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with Monogamy behaviour)", x = "Generation")

ggplot(data = tstrun, aes(y = p_c, x = generation)) +
  geom_line(aes(colour = run)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with Monogamy cultural trait)", x = "Generation")

ggplot(data = tstrun, aes(y = p_b, x = generation)) +
  geom_line(aes(colour = run)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with Monogamy biological trait)", x = "Generation")

ggplot(data = tstrun, aes(y = p_mono, x = generation)) +
  geom_line(aes(colour = run)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with observable Monogamy)", x = "Generation")



tstrun2 <- ABMmodel2(100, 100, 10, b_f = 0.3)

ggplot(data = tstrun2, aes(y = p_bh, x = generation)) +
  geom_line(aes(colour = run)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with Monogamy behaviour)", x = "Generation")

ggplot(data = tstrun2, aes(y = p_c, x = generation)) +
  geom_line(aes(colour = run)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with Monogamy cultural trait)", x = "Generation")

ggplot(data = tstrun2, aes(y = p_b, x = generation)) +
  geom_line(aes(colour = run)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with Monogamy biological trait)", x = "Generation")

ggplot(data = tstrun2, aes(y = p_mono, x = generation)) +
  geom_line(aes(colour = run)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with observable Monogamy)", x = "Generation")









########################################### Step-by-step testing ###############################################################

# Create first generation
population <- tibble(sex = sample(c("Female", "Male"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)),
                     trait_c = sample(c("Monogamy", "Polygamy"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)), 
                     trait_b = sample(c("Monogamy", "Polygamy"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)),
                     bh = as.character(rep(NA, 100)))

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
    sample(c("Monogamy", "Polygamy"), sum(is.na(population$bh)&population$trait_c == "Monogamy"), prob = c(0.5, 1 - 0.5), replace = TRUE)
  population[is.na(population$bh)& population$trait_c == "Polygamy",]$bh <-
    sample(c("Polygamy", "Monogamy"), sum(is.na(population$bh)&population$trait_c == "Polygamy"), prob = c(0.5, 1 - 0.5), replace = TRUE)
}

p_bh <- sum(population$bh == "Monogamy") / 100

previous_population <- population


matingpool <- previous_population %>%
  add_column(bonded = "No", bond_ID = as.character("NA"), offspring = 0, ID = 1:100)

# Create empty offspring tibble
# This is the empty table of offspring that will constitute the next generation once filled by the mating loop

nextgen <- tibble(sex = sample(c("Female", "Male"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)),
                  trait_c = as.character(rep(NA, 100)), 
                  trait_b = as.character(rep(NA, 100)),
                  bh = as.character(rep(NA, 100)))

# The next loop is the mating process script that I wrote separately (this is not from the Acerbi book, though it used concepts an structures introduced there)

# Each iteration of the loop is a single mating encounter with a 100% chance of reproduction. This is an assumption that could eventually be relaxed.

# Mating is random - the only rule is that monogamous individuals mate only once. Polygamous individuals are not limited in their mating opportunities (for now)

# Because I wrote this myself, this is likelihood the least optimised bit of this entire script for now.

for (i in 1:100) {
  
  # Randomly draw the first participant of this mating encounter from the mating pool
  # NB. Has to be Monogamous and Unmated or Polygamous.
  
  mate1 <- slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$offspring < 3)|(matingpool$bh=="Polygamy"),], n=1)
  
  # Based on the sex of the first participant, an opposite sex partner is drawn from the same population (same rules)
  
  if (mate1$sex=="Female" & mate1$bonded=="No") {
    mate2 <- 
      slice_sample(matingpool[matingpool$sex=="Male" & matingpool$bonded == "No",], n=1)
  } 
  
  if (mate1$sex=="Male" & mate1$bonded=="No") {
    mate2 <- 
      slice_sample(matingpool[matingpool$sex=="Female" & matingpool$bonded == "No",], n=1)
  } 
  
  if (mate1$bonded=="Yes") {
    mate2 <- 
      matingpool[mate1$bond_ID,]
  }
  
  # Record mating encounter in matingpool tibble
  
  if (mate1$bh=="Monogamy"){
  matingpool[mate1$ID,]$bonded <- "Yes"
  matingpool[mate1$ID,]$bond_ID <- as.character(mate2$ID)
  }
  
  if (mate2$bh=="Monogamy") {
  matingpool[mate2$ID,]$bonded <- "Yes"
  matingpool[mate2$ID,]$bond_ID <- as.character(mate1$ID)
  }
  
  matingpool[mate1$ID,]$offspring <- matingpool[mate1$ID,]$offspring + 1
  matingpool[mate2$ID,]$offspring <- matingpool[mate2$ID,]$offspring + 1
  
  ## Next, we start to populate the data row for the offspring that this encounter produced, starting with the cultural trait
  
  # Both parents are culturally monogamous, thus child adopts cultural Monogamy
  if (mate1$trait_c == "Monogamy" & mate2$trait_c == "Monogamy") {
    nextgen$trait_c[i] <- "Monogamy"  
  }
  
  # Both parents are culturally polygamous, thus child adopts cultural Polygamy
  if (mate1$trait_c == "Polygamy" & mate2$trait_c == "Polygamy") {
    nextgen$trait_c[i] <- "Polygamy" 
  }
  
  # If they have mixed parental cultural traits (i.e. one Monogamy and one Polygamy parent)
  if (is.na(nextgen$trait_c[i])) {  
    # They adopt cultural Monogamy with probability b_c
    nextgen$trait_c[i] <- 
      sample(c("Monogamy", "Polygamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  }
  
  # Determine 'mutant' individuals (for cultural traits, this is roughly analogous to "individual learners")
  mutate_c <- sample(c(TRUE, FALSE), 1, prob = c(0, 1), replace = TRUE)
  
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
  
  # Rename to standardize
  
  if (nextgen$trait_c[i]=="Polygamy_m"){
    nextgen$trait_c[i] <- "Polygamy"
  }
  
  ## We repeat this with the biological traits
  
  # Both parents are biologically monogamous, thus child inherits biological Monogamy
  if (mate1$trait_b == "Monogamy" & mate2$trait_b == "Monogamy") {
    nextgen$trait_b[i] <- "Monogamy"  
  }
  
  # Both parents are culturally polygamous, thus child inherits biological Polygamy
  if (mate1$trait_b == "Polygamy" & mate2$trait_b == "Polygamy") {
    nextgen$trait_b[i] <- "Polygamy" 
  }
  
  # If they have mixed parental biological traits (i.e. one Monogamy and one Polygamy parent)
  if (is.na(nextgen$trait_b[i])) {  
    # They adopt biological Monogamy with probability b_b
    nextgen$trait_b[i] <- 
      sample(c("Monogamy", "Polygamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  }
  
  # Determine 'mutant' individuals (for cultural traits, this is roughly analogous to "individual learners")
  mutate_b <- sample(c(TRUE, FALSE), 1, prob = c(0, 1), replace = TRUE)
  
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
  
  # Rename to standardize
  
  if (nextgen$trait_b[i]=="Polygamy_m"){
    nextgen$trait_b[i] <- "Polygamy"
  }
  
  ## Determining behavioural tendencies based on cultural and biological trait combination
  
  # Both traits favour monogamy, thus behavioural preference is monogamous
  if (nextgen$trait_c[i] == "Monogamy" & nextgen$trait_b[i] == "Monogamy") {
    nextgen$bh[i] <- "Monogamy"  
  }
  
  # Both traits favour polygamy, thus behavioural preference is polygamous
  if (nextgen$trait_c[i] == "Polygamy" & nextgen$trait_b[i] == "Polygamy") {
    population$bh[i] <- "Polygamy" 
  }
  
  # If the have mixed cultural/biological traits
  if (is.na(nextgen$bh[i]) & nextgen$trait_c[i] == "Monogamy") {  
    # They will adopt the behavioural variant corresponding to their cultural trait with probability cb_bh
    nextgen$bh[i] <- 
      sample(c("Monogamy", "Polygamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  } else {
    nextgen$bh[i] <- 
      sample(c("Polygamy", "Monogamy"), 1, prob = c(0.5, 1 - 0.5), replace = TRUE)
  }
  
  
  
}

