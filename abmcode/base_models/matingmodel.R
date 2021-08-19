library(tidyverse)
library(dplyr)

## Custom mating market function

matingmarket_sim <- function(input, output, N = N, mono_maxRS = mono_max, poly_maxRS = poly_max){
  
  matingpool <- input
  matepairings <- output
  
  while (sum(is.na(matepairings$trait_b_p1)) > 0) {
    
    
    # The following code is necessary to make sure that any given population, regardless of composition, can maintain its fixed population size
    # In rare cases, all available mates of either sex have been paired off and are withdrawn from the matingpool before there has been enough reproduction to maintain the population size
    # In that case, additional reproductive events are drawn from the existing population of unique pairings
    
    if(sum(matingpool[matingpool$sex=="Female",]$available=="Yes") == 0 | sum(matingpool[matingpool$sex=="Male",]$available=="Yes") == 0) {
      
      cellref_end <- cellref_end + 1
      
      uniq_pairings <- distinct(select(matepairings,-sex))
      
      matepairings[cellref_end,2:7] <- slice_sample(uniq_pairings[1:nrow(uniq_pairings)-1,],n=1)
      
      matingpool[matepairings[cellref_end,]$ID_p1,]$offspring <- matingpool[matepairings[cellref_end,]$ID_p1,]$offspring + 1
      matingpool[matepairings[cellref_end,]$ID_p2,]$offspring <- matingpool[matepairings[cellref_end,]$ID_p2,]$offspring + 1
      
      next
    }
    
    # Randomly draw the first participant of this mating encounter from the mating pool
    # NB. Has to be Monogamous and Unmated or Polygamous.
    
    mate1 <- slice_sample(matingpool[(matingpool$bh=="Monogamy" & matingpool$available == "Yes")|(matingpool$bh=="Polygamy" & matingpool$available == "Yes" & matingpool$matesearched == "No"),], n=1)
    
    # Based on the sex of the first participant, an opposite sex partner is drawn from the same population (same rules)
    
    if (mate1$sex=="Female") {
      mate2 <- 
        slice_sample(matingpool[matingpool$sex=="Male" & matingpool$available=="Yes",], n=1)
    } 
    
    if (mate1$sex=="Male") {
      mate2 <- 
        slice_sample(matingpool[matingpool$sex=="Female" & matingpool$available=="Yes",], n=1)
    } 
    
    # Record mating encounter in matingpool tibble
    
    matingpool[mate1$ID,]$matesearched <- "Yes"
    
    if (mate1$bh=="Monogamy"| (mate1$bh=="Polygamy" & mate1$partners==1)){
      matingpool[mate1$ID,]$available <- "No"
    }
    
    if (mate1$bond_ID1=="NA") {
      matingpool[mate1$ID,]$bond_ID1 <- as.character(mate2$ID)
    } else {
      matingpool[mate1$ID,]$bond_ID2 <- as.character(mate2$ID)
    }
    
    if (mate2$bh=="Monogamy"| (mate2$bh=="Polygamy" & mate2$partners==1)) {
      matingpool[mate2$ID,]$available <- "No"
    }
    
    if (mate1$bond_ID1=="NA") {
      matingpool[mate2$ID,]$bond_ID1 <- as.character(mate1$ID)
    } else {
      matingpool[mate2$ID,]$bond_ID2 <- as.character(mate1$ID)
    }
    
    if (mate1$bh=="Monogamy"){
      RS <- sample(1:mono_maxRS, 1)
    }
    
    if (mate1$bh=="Polygamy"){
      RS <- sample(1:poly_maxRS, 1)
    }
    
    cellref_start <- N + 1 - sum(is.na(matepairings$trait_b_p1))
    cellref_end <- cellref_start + RS - 1
    
    if (cellref_end > N) {
      RS_old <- RS
      RS <- RS_old - (cellref_end - N) 
      cellref_end <- N
    }
    
    
    for (i in cellref_start:cellref_end) {
      
      matepairings[i,]$trait_c_p1 <- mate1$trait_c
      matepairings[i,]$trait_b_p1 <- mate1$trait_b
      matepairings[i,]$trait_c_p2 <- mate2$trait_c
      matepairings[i,]$trait_b_p2 <- mate2$trait_b
      
    }
    
    matepairings[cellref_start:cellref_end,]$ID_p1 <- rep(as.character(mate1$ID), cellref_end - cellref_start + 1)
    matepairings[cellref_start:cellref_end,]$ID_p2 <- rep(as.character(mate2$ID), cellref_end - cellref_start + 1)
    
    matingpool[mate1$ID,]$offspring <- matingpool[mate1$ID,]$offspring + RS
    matingpool[mate1$ID,]$partners <- matingpool[mate1$ID,]$partners + 1
    matingpool[mate2$ID,]$offspring <- matingpool[mate2$ID,]$offspring + RS
    matingpool[mate2$ID,]$partners <- matingpool[mate2$ID,]$partners + 1
    
    if(mate1$bh=="Monogamy" | matingpool[mate1$ID,]$partners == 2 | sum(is.na(matepairings$trait_b_p1)==0)){
      next
    }
    
    if(mate1$sex == "Male"){
      OSR <- sum(matingpool$sex == "Female" & matingpool$available=="Yes")/sum(matingpool$available=="Yes")
    }
    
    if(mate1$sex == "Female"){
      OSR <- sum(matingpool$sex == "Male" & matingpool$available=="Yes")/sum(matingpool$available=="Yes")
    }
    
    poly_chance <- sample(c(TRUE,FALSE), prob = c(OSR, 1 - OSR),1)
    
    if(isFALSE(poly_chance)){
      next
    }
    
    if (mate1$sex=="Female") {
      mate2 <- 
        slice_sample(matingpool[matingpool$sex=="Male" & matingpool$available=="Yes",], n=1)
    } 
    
    if (mate1$sex=="Male") {
      mate2 <- 
        slice_sample(matingpool[matingpool$sex=="Female" & matingpool$available=="Yes",], n=1)
    } 
    
    # Record mating encounter in matingpool tibble
    
    if (matingpool[mate1$ID,]$partners==1){
      matingpool[mate1$ID,]$available <- "No"
    }
    
    matingpool[mate1$ID,]$bond_ID2 <- as.character(mate2$ID)
    
    if (mate2$bh=="Monogamy" | (mate2$bh=="Polygamy" & mate2$partners==1)) {
      matingpool[mate2$ID,]$available <- "No"
    }
    
    if (mate2$bond_ID1=="NA") {
      matingpool[mate2$ID,]$bond_ID1 <- as.character(mate1$ID)
    } else {
      matingpool[mate2$ID,]$bond_ID2 <- as.character(mate1$ID)
    }
    
    RS <- sample(1:poly_maxRS, 1)
    
    cellref_start <- N + 1 - sum(is.na(matepairings$trait_b_p1))
    cellref_end <- cellref_start + RS - 1
    
    for(i in cellref_start:cellref_end) {
      
      matepairings[i,]$trait_c_p1 <- mate1$trait_c
      matepairings[i,]$trait_b_p1 <- mate1$trait_b
      matepairings[i,]$trait_c_p2 <- mate2$trait_c
      matepairings[i,]$trait_b_p2 <- mate2$trait_b
      
    }
    
    matingpool[mate1$ID,]$offspring <- matingpool[mate1$ID,]$offspring + RS
    matingpool[mate1$ID,]$partners <- matingpool[mate1$ID,]$partners + 1
    matingpool[mate2$ID,]$offspring <- matingpool[mate2$ID,]$offspring + RS
    matingpool[mate2$ID,]$partners <- matingpool[mate2$ID,]$partners + 1
    
    
  }
  
  matepairings -> output
  matingpool -> input
  return(lst(output,input))
  
}
