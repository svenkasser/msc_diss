library(tidyverse)
library(dplyr)
library(ggplot2)



ABMmodel2 <- function(N, t_max, r_max, mu_c = 0, mu_b = 0, pc_0 = 0.5, pb_0 = 0.5, b_c = 0.5, b_b = 0.5, cb_bh = 0.5, pf_0 = 0.5, b_f = 0.5) {
  
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
        add_column(mated = "No", offspring = 0, ID = 1:N)
      
      # Create empty offspring tibble
      # This is the empty table of offspring that will constitute the next generation once filled by the mating loop
      
      nextgen <- tibble(sex = sample(c("Female", "Male"), N, replace = TRUE, prob = c(b_f, 1 - b_f)),
                        trait_c = as.character(rep(NA, N)), 
                        trait_b = as.character(rep(NA, N)),
                        bh = as.character(rep(NA, N)))
      
      matepairings <- tibble(mate1_ID = as.character(rep(NA, N),
                              mate1_sex = as.character(rep(NA, N)),
                              mate1_trait_c = as.character(rep(NA, N)),
                              mate1_trait_b = as.character(rep(NA, N)),
                              mate1_ID = as.character(rep(NA, N),
                              mate2_sex = as.character(rep(NA, N)),
                              mate2_trait_c = as.character(rep(NA, N)),
                              mate2_trait_b = as.character(rep(NA, N)),)
      
     
      
      for (i in 1:N) {
        
        # Randomly draw the first participant of this mating encounter from the mating pool
        # NB. Has to be Monogamous and Unmated or Polygamous.
        
        mate1 <- slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$mated=="No")|(matingpool$bh=="Polygamy"),], n=1)
        
        # Based on the sex of the first participant, an opposite sex partner is drawn from the same population (same rules)
        
        if (mate1$sex=="Female") {
          mate2 <- 
            slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$mated=="No" & matingpool$sex=="Male")|(matingpool$bh=="Polygamy" & matingpool$sex=="Male"),], n=1)
        }
        else if (mate1$sex=="Male") {
          mate2 <- 
            slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$mated=="No" & matingpool$sex=="Female")|(matingpool$bh=="Polygamy" & matingpool$sex=="Female"),], n=1)
        }
        
        # Record mating encounter in matingpool tibble
        
        matingpool[mate1$ID,]$mated <- "Yes"
        matingpool[mate2$ID,]$mated <- "Yes"
        
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
            sample(c("Monogamy", "Polygamy"), 1, prob = c(b_c, 1 - b_c), replace = TRUE)
        }
        
        # Determine 'mutant' individuals (for cultural traits, this is roughly analogous to "individual learners")
        mutate_c <- sample(c(TRUE, FALSE), 1, prob = c(mu_c, 1 - mu_c), replace = TRUE)
        
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
            sample(c("Monogamy", "Polygamy"), 1, prob = c(b_b, 1 - b_b), replace = TRUE)
        }
        
        # Determine 'mutant' individuals (for cultural traits, this is roughly analogous to "individual learners")
        mutate_b <- sample(c(TRUE, FALSE), 1, prob = c(mu_b, 1 - mu_b), replace = TRUE)
        
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
        sum(matingpool$mated == "Yes" & matingpool$offspring == 1) / N
      
      
    }
    
  }
  
  # Export data from function
  output
  
}



###########################################


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
  add_column(mated = "No", offspring = 0, ID = 1:100)



# Create empty offspring tibble
# This is the empty table of offspring that will constitute the next generation once filled by the mating loop

nextgen <- tibble(sex = sample(c("Female", "Male"), 100, replace = TRUE, prob = c(0.5, 1 - 0.5)),
                  trait_c = as.character(rep(NA, 100)), 
                  trait_b = as.character(rep(NA, 100)),
                  bh = as.character(rep(NA, 100)))

matepairings <- tibble(mate1_ID = sample(1:100, 100, replace = TRUE),
                       mate1_bh = as.character(rep(NA, 100)),
                       mate2_ID = as.character(rep(NA, 100)),
                       pairing_ID = as.character(rep(1:100,1)))


# Making sure that monogamists only mate once

matepairings$mate1_bh <- matingpool[matepairings$mate1_ID,]$bh



matepairings[duplicated(matepairings$mate1_ID),]$mate1_ID












matepairings[replace_monogamist1,]$mate1_ID <- sample(matingpool[matingpool$bh=="Polygamy"|(matingpool[matepairings[replace_monogamist1==FALSE,]$mate1_ID]),]$ID,
                                                      sum(replace_monogamist1),
                                                      replace = TRUE)














matingpool[matingpool$bh=="Polygamy"|(matingpool[matepairings[replace_monogamist1==FALSE,]$mate1_ID,]),]$ID

matingpool$bh=="Polygamy"

matingpool[matingpool$bh=="Polygamy",]$ID
matingpool[matingpool[matepairings[replace_monogamist1==FALSE,]$mate1_ID,]$ID,]

matingpool[matingpool$bh=="Polygamy",]

matingpool$bh=="Polygamy"

matepairings[replace_monogamist1==FALSE,]$mate1_ID


#sum(replace_monogamist1)


#while (sum(replace_monogamist1) > 0) {
#  matepairings[replace_monogamist1,]$mate1_ID <- sample(1:100, sum(replace_monogamist1), replace = TRUE)
  
#  replace_monogamist <- matingpool[matepairings$mate1_ID,]$bh == "Monogamy" & duplicated(matepairings$mate1_ID)
#}

# Find mate 2

matepairings[matingpool[matepairings$mate1_ID,]$sex == "Male",]$mate2_ID <- sample(as.character(matingpool[matingpool$sex == "Female",]$ID), 
                                  size=sum(matingpool[matepairings$mate1_ID,]$sex == "Male"),
                                  replace = TRUE)

matepairings[matingpool[matepairings$mate1_ID,]$sex == "Female",]$mate2_ID <- sample(as.character(matingpool[matingpool$sex == "Male",]$ID), 
                                                                                   size=sum(matingpool[matepairings$mate1_ID,]$sex == "Female"),
                                                                                   replace = TRUE)


sum(matingpool[matepairings$mate1_ID,]$sex == "Male")

matepairings[matingpool[matepairings$mate1_ID,]$sex == "Male",]


replace_monogamist2 <- matingpool[matepairings$mate1_ID,]$bh == "Monogamy" & duplicated(matepairings$mate1_ID)


polygnists <- matingpool[matepairings$mate1_ID,]$bh == "Polygamy" & duplicated(matepairings$mate1_ID)

sum(replace_monogamist1)
sum(polygnists)

sample(matingpool[matingpool$sex == "Female",]$ID, 
       size=37,
       replace = TRUE)
