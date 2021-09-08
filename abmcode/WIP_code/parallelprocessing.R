library(tidyverse)
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)

### Set up parallel processing

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

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

n.cores <- parallel::detectCores() - 1

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

foreach::getDoParWorkers()


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



ABMmodel5 <- function(N, t_max, r_max, pc_0 = 0.5, pb_0 = 0.5, b_c = 0.5, b_b = 0.5, cb_bh = 0.5, pf_0 = 0.5, b_f = 0.5, mono_max = 3, poly_max = 2) {
  
  # Create output file - this is he data what we want to save from each run and generation
  
  # Loop for each run (essentially each simulated population)
  
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
    
    output <- tibble(generation = 1, 
                     p_c = as.numeric(rep(NA, 1)),
                     p_b = as.numeric(rep(NA, 1)),
                     p_bh = as.numeric(rep(NA, 1)),
                     p_f = as.numeric(rep(NA, 1)),
                     p_mono = as.numeric(rep(NA, 1)),
                     p_polyg = as.numeric(rep(NA,1)),
                     p_polya = as.numeric(rep(NA,1)),
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
      
      
      matingpool <- previous_population %>%
        add_column(matesearched = "No", available = "Yes", bond_ID1 = as.character("NA"), bond_ID2 = as.character("NA"), partners = 0, offspring = 0, ID = 1:N)
      
      # Create empty offspring tibble
      # This is the empty table of offspring that will constitute the next generation once filled by the mating loop
      
      matepairings <- tibble(sex = sample(c("Female", "Male"), N, replace = TRUE, prob = c(b_f, 1 - b_f)),
                             ID_p1 = as.character(rep(NA, N)),
                             trait_c_p1 = as.character(rep(NA, N)),
                             trait_b_p1 = as.character(rep(NA, N)),
                             ID_p2= as.character(rep(NA, N)),
                             trait_c_p2 = as.character(rep(NA, N)),
                             trait_b_p2 = as.character(rep(NA, N)))
      
      
      matingmodel_output <- matingmarket_sim2(matingpool, matepairings, N = N, mono_maxRS = mono_max, poly_maxRS = poly_max)
      
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
      
      
      nextgen <- as_tibble(cbind(matingmodel_output[1]$output$sex, 
                                 matingmodel_output[1]$output$trait_c_p1, 
                                 matingmodel_output[1]$output$trait_b_p1, 
                                 matingmodel_output[1]$output$trait_c_p2, 
                                 matingmodel_output[1]$output$trait_b_p2))
      
      nextgen <- rename(nextgen, sex = V1, trait_c_p1 = V2, trait_b_p1 = V3, trait_c_p2 = V4 , trait_b_p2 = V5)
      
      population <- select(nextgen, -c(trait_c_p1, trait_c_p2, trait_b_p1, trait_b_p2)) %>%
        add_column(trait_c = as.character(rep(NA, N)), 
                   trait_b = as.character(rep(NA, N)),
                   bh = as.character(rep(NA, N)))
      
      population <- trait_determination(population, nextgen, b_c = b_c, b_b = b_b, cb_bh = cb_bh)
      
      new_output <- tibble(generation = t, 
                           p_c = as.numeric(rep(NA, 1)),
                           p_b = as.numeric(rep(NA, 1)),
                           p_bh = as.numeric(rep(NA, 1)),
                           p_f = as.numeric(rep(NA, 1)),
                           p_mono = as.numeric(rep(NA, 1)),
                           p_polyg = as.numeric(rep(NA, 1)),
                           p_polya = as.numeric(rep(NA, 1)),
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
      
      output[output$generation == t-1 & output$run == 1, ]$p_polyg <-
        sum(matingpool$partners == 2 & matingpool$sex == "Male") / sum(matingpool$sex=="Male")
      
      output[output$generation == t-1 & output$run == 1, ]$p_polya <-
        sum(matingpool$partners == 2 & matingpool$sex == "Female") / sum(matingpool$sex=="Female")
      
      
      output <- rbind(output,new_output)
      
    }
    return(output)
  }
  
  # Export data from function
  final_output
  
}

Before2 <- Sys.time()

tstst2 <- ABMmodel5(150, 5, 2, mono_max = 1)

After2 <- Sys.time()

Runtime2 <- After1 - Before1

Runtime2



ggplot(data = tstst2, aes(y = p_polygyny, x = generation)) +
  geom_line(aes(colour = run)) +
  #stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with Monogamy behaviour)", x = "Generation")




######## Testing

N = 100
t_max = 10
r_max = 3
pc_0 = 0.5
pb_0 = 0.5
b_c = 0.5
b_b = 0.5
cb_bh = 0.5
pf_0 = 0.5
b_f = 0.5
mono_max = 1
poly_max = 3

final_output <- foreach(r = 1:r_max, .packages = 'tidyverse', .combine = 'rbind') %dopar% {
  
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
  
  output <- tibble(generation = 1, 
                   p_c = as.numeric(rep(NA, 1)),
                   p_b = as.numeric(rep(NA, 1)),
                   p_bh = as.numeric(rep(NA, 1)),
                   p_f = as.numeric(rep(NA, 1)),
                   p_mono = as.numeric(rep(NA, 1)),
                   RS_poly = as.numeric(rep(NA, 1)),
                   RS_mono = as.numeric(rep(NA, 1)),
                   p_polygyny = as.numeric(rep(NA, 1)),
                   p_polyandry = as.numeric(rep(NA, 1)),
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
    
    
    matingpool <- previous_population %>%
      add_column(matesearched = "No", available = "Yes", bond_ID1 = as.character("NA"), bond_ID2 = as.character("NA"), partners = 0, offspring = 0, ID = 1:100)
    
    # Create empty offspring tibble
    # This is the empty table of offspring that will constitute the next generation once filled by the mating loop
    
    matepairings <- tibble(sex = sample(c("Female", "Male"), N, replace = TRUE, prob = c(b_f, 1 - b_f)),
                           ID_p1 = as.character(rep(NA, N)),
                           trait_c_p1 = as.character(rep(NA, N)),
                           trait_b_p1 = as.character(rep(NA, N)),
                           ID_p2= as.character(rep(NA, N)),
                           trait_c_p2 = as.character(rep(NA, N)),
                           trait_b_p2 = as.character(rep(NA, N)))
    
    
    matingmodel_output <- matingmarket_sim2(matingpool, matepairings, N = N, mono_maxRS = mono_max, poly_maxRS = poly_max)
    
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
    
    
    nextgen <- as_tibble(cbind(matingmodel_output[1]$output$sex, 
                               matingmodel_output[1]$output$trait_c_p1, 
                               matingmodel_output[1]$output$trait_b_p1, 
                               matingmodel_output[1]$output$trait_c_p2, 
                               matingmodel_output[1]$output$trait_b_p2))
    
    nextgen <- rename(nextgen, sex = V1, trait_c_p1 = V2, trait_b_p1 = V3, trait_c_p2 = V4 , trait_b_p2 = V5)
    
    population <- select(nextgen, -c(trait_c_p1, trait_c_p2, trait_b_p1, trait_b_p2)) %>%
      add_column(trait_c = as.character(rep(NA, N)), 
                 trait_b = as.character(rep(NA, N)),
                 bh = as.character(rep(NA, N)))
    
    population <- trait_determination(population, nextgen, b_c = b_c, b_b = b_b, cb_bh = cb_bh)
    
    new_output <- tibble(generation = t, 
                     p_c = as.numeric(rep(NA, 1)),
                     p_b = as.numeric(rep(NA, 1)),
                     p_bh = as.numeric(rep(NA, 1)),
                     p_f = as.numeric(rep(NA, 1)),
                     p_mono = as.numeric(rep(NA, 1)),
                     RS_poly = as.numeric(rep(NA, 1)),
                     RS_mono = as.numeric(rep(NA, 1)),
                     p_polygyny = as.numeric(rep(NA, 1)),
                     p_polyandry = as.numeric(rep(NA, 1)),           
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
    
    # Get average RS for polygamists
    output[output$generation == t-1 & output$run == r, ]$RS_poly <-
      mean(as.integer((matingpool[matingpool$bh=="Polygamy",]$offspring)))
    
    # Get average RS for monogamists
    output[output$generation == t-1 & output$run == r, ]$RS_mono <-
      mean(as.integer((matingpool[matingpool$bh=="Monogamy",]$offspring)))
   
    # Get average RS for polygamists
    output[output$generation == t-1 & output$run == r, ]$p_polygyny <-
      sum(matingpool$sex == "Male" & matingpool$partners > 0) / sum(population$sex == "Male") 
    
    # Get average RS for monogamists
    output[output$generation == t-1 & output$run == r, ]$p_polyandry <-
      sum(matingpool$sex == "Female" & matingpool$partners > 0) / sum(population$sex == "Female") 
    
    
    output <- rbind(output,new_output)
    
  }
  return(output)
}


final_output

mean(as.integer((matingpool$partners)))

mean(as.integer((matingpool$offspring)))


mean(as.integer((matingpool[matingpool$bh=="Polygamy",]$offspring)))

mean(as.integer((matingpool[matingpool$bh=="Monogamy",]$offspring)))




mean(as.integer((matingpool$offspring)))


final_output

final_output$p_polyandry

lm(final_output$p_polygny ~ final_output$p_polyandry, data=final_output)

ggplot(data = final_output, aes(y = p_polygyny, x = generation)) +
  geom_line(aes(colour = run)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "p (proportion of individuals with Monogamy behaviour)", x = "Generation")
