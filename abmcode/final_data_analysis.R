# Full data analysis code

# Load required packages

## Key packages for visualisation / data handling

library(ggplot2)
library(tidyverse)

## Various packages to help arrange plots

library(gridExtra)
library(ggpubr)
library(lemon)
library(grid)
library(cowplot)
library(ggstatsplot)

## Package for robust analysis variant

library(robustbase)
library(lme4)

## Package for APA-friendly regression table

library(stargazer)

# Write helper function to extract summary statistics

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Load required data

## Load and transform data

Modelrun_1_a <- read.csv("Modelrun_1_a.csv")
Modelrun_2_a <- read.csv("Modelrun_2_a.csv")
Modelrun_3_a <- read.csv("Modelrun_3_a.csv")
Modelrun_4_a <- read.csv("Modelrun_4_a.csv")
Modelrun_5_a <- read.csv("Modelrun_5_a.csv")
Modelrun_6_a <- read.csv("Modelrun_6_a.csv")
Modelrun_7_a <- read.csv("Modelrun_7_a.csv")
Modelrun_8_a <- read.csv("Modelrun_8_a.csv")
Modelrun_9_a <- read.csv("Modelrun_9_a.csv")
Modelrun_10_a <- read.csv("Modelrun_10_a.csv")

Modelrun_1_a$run <- as.factor(Modelrun_1_a$run)
Modelrun_2_a$run <- as.factor(Modelrun_2_a$run)
Modelrun_3_a$run <- as.factor(Modelrun_3_a$run)
Modelrun_4_a$run <- as.factor(Modelrun_4_a$run)
Modelrun_5_a$run <- as.factor(Modelrun_5_a$run)
Modelrun_6_a$run <- as.factor(Modelrun_6_a$run)
Modelrun_7_a$run <- as.factor(Modelrun_7_a$run)
Modelrun_8_a$run <- as.factor(Modelrun_8_a$run)
Modelrun_9_a$run <- as.factor(Modelrun_9_a$run)
Modelrun_10_a$run <- as.factor(Modelrun_10_a$run)

## Corresponding parameterisations

# Run 1a - mono_max = 5, poly_max = 3, pf_0 = 0.5, b_f = 0.5
# Run 2a - mono_max = 5, poly_max = 3, pf_0 = 0.2, b_f = 0.5
# Run 3a - mono_max = 5, poly_max = 3, pf_0 = 0.4, b_f = 0.5
# Run 4a - mono_max = 5, poly_max = 3, pf_0 = 0.6, b_f = 0.5
# Run 5a - mono_max = 5, poly_max = 3, pf_0 = 0.8, b_f = 0.5

# Run 6a - mono_max = 5, poly_max = 3, pf_0 = 0.5, b_f = 0.5 # This is identical to scenario 1
# Run 7a - mono_max = 5, poly_max = 3, pf_0 = 0.2, b_f = 0.2
# Run 8a - mono_max = 5, poly_max = 3, pf_0 = 0.4, b_f = 0.4
# Run 9a - mono_max = 5, poly_max = 3, pf_0 = 0.6, b_f = 0.6
# Run 10a - mono_max = 5, poly_max = 3, pf_0 = 0.8, b_f = 0.8

# Create "Final state" datasets (for end-of-run analyses)

A1_final <- Modelrun_1_a[Modelrun_1_a$generation==500,]
A2_final <- Modelrun_2_a[Modelrun_2_a$generation==500,]
A3_final <- Modelrun_3_a[Modelrun_3_a$generation==500,]
A4_final <- Modelrun_4_a[Modelrun_4_a$generation==500,]
A5_final <- Modelrun_5_a[Modelrun_5_a$generation==500,]

final_a_startingSR <- bind_rows(A1_final, A2_final, A3_final, A4_final, A5_final, .id = "Parameter setting")
final_a_startingSR$`Parameter setting` <- as.factor(c(rep("0.5", 100), rep("0.2", 100), rep("0.4", 100), rep("0.6", 100), rep("0.8", 100)))

A6_final <- Modelrun_6_a[Modelrun_6_a$generation==500,]
A7_final <- Modelrun_7_a[Modelrun_7_a$generation==500,]
A8_final <- Modelrun_8_a[Modelrun_8_a$generation==500,]
A9_final <- Modelrun_9_a[Modelrun_9_a$generation==500,]
A10_final <- Modelrun_10_a[Modelrun_10_a$generation==500,]

final_a_contSR <- bind_rows(A6_final, A7_final, A8_final, A9_final, A10_final, .id = "Parameter setting")
final_a_contSR$`Parameter setting` <- as.factor(c(rep("0.5", 100), rep("0.2", 100), rep("0.4", 100), rep("0.6", 100), rep("0.8", 100)))

## Create simple overview table for equilibrium outcomes (outcomes where only one behaviour type still exists in the population)

final_a <- bind_rows(final_a_startingSR, final_a_contSR, .id = "Sex ratio variant")
final_a$`Sex ratio variant` <- as.factor(c(rep("Starting", 500), rep("Continuous", 500)))
final_a <- final_a %>%
  add_column(equilibrium = as.character("NA"))
final_a[final_a$p_bh==0,]$equilibrium <- c("Poly_complete")
final_a[final_a$p_bh==1,]$equilibrium <- c("Mono_complete")
final_a[final_a$p_bh<1 & final_a$p_bh>0,]$equilibrium <- c("Mixed")

table(final_a$equilibrium)

xtabs(~equilibrium + `Parameter setting` + `Sex ratio variant`, data=final_a)

## Compute summary statistic of trait frequencies at model termination (timestep 500)

final_a_summary_bh <- summarySE(final_a, measurevar="p_bh", groupvars=c("`Parameter setting`","`Sex ratio variant`"))
final_a_summary_c <- summarySE(final_a, measurevar="p_c", groupvars=c("`Parameter setting`","`Sex ratio variant`"))
final_a_summary_b <- summarySE(final_a, measurevar="p_b", groupvars=c("`Parameter setting`","`Sex ratio variant`"))

## Visualisation (not included in dissertation as not particularly informative)

ggplot(final_a_summary_bh, aes(x=`Parameter setting`, y=p_bh, fill=`Sex ratio variant`)) + 
  geom_bar(position=position_dodge(), stat="identity", alpha=0.5) +
  geom_errorbar(aes(ymin=p_bh-se, ymax=p_bh+se), width=.2,position=position_dodge(.9)) +
  theme_bw() +
  ylim(c(0,0.1)) +
  labs(y = "Proportion of individuals \nwith Monogamy behaviour", x = "Parameter value for sex ratio") +
  ggtitle("Effect of sex ratio parameter \non frequency of monogamous behaviour at end of simulation") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title.y=element_text(size=10))


# Timeseries visualisation - Plotting the dynamics of monogamy trait frequencies

## Plots for "Starting SR" conditions

### Phenotype behaviour plot

fixed1_a <- ggplot(data = Modelrun_1_a, aes(y = p_bh, x = generation)) +
  stat_summary(data = Modelrun_2_a, aes(colour="Strongly male-biased \nSR (0.2)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_3_a, aes(colour="Slightly male-biased \nSR (0.4)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_1_a, aes(colour="Even starting SR"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_4_a, aes(colour="Slightly female-biased \nSR (0.6)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_5_a, aes(colour="Strongly femmale-biased \nSR (0.8)"), fun = mean, geom = "line", size = 1) +
  #ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "Proportion Monogamy (behaviour)", x = "Generation", col = "Parameter setting") +
  #ggtitle("Timeseries of average trait frequency \nfor Monogamy across 100 runs (Starting SR Manipulation)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = "bottom")

### Cultural trait plot

fixed2_a <- ggplot(data = Modelrun_1_a, aes(y = p_c, x = generation)) +
  stat_summary(data = Modelrun_2_a, aes(colour="Strongly male-biased \nstarting SR (0.2)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_3_a, aes(colour="Slightly male-biased \nstarting SR (0.4)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_1_a, aes(colour="Even starting SR"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_4_a, aes(colour="Slightly female-biased \nstarting SR (0.6)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_5_a, aes(colour="Strongly femmale-biased \nstarting SR (0.8)"), fun = mean, geom = "line", size = 1) +
  #ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "Proportion Monogamy (cultural trait)", x = "Generation", col = "Parameter setting") +
  #ggtitle("Starting SR Manipulation") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = "right")

### Biological trait plot

fixed3_a <- ggplot(data = Modelrun_1_a, aes(y = p_b, x = generation)) +
  stat_summary(data = Modelrun_2_a, aes(colour="Strongly male-biased \nstarting SR (0.2)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_3_a, aes(colour="Slightly male-biased \nstarting SR (0.4)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_1_a, aes(colour="Even starting SR"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_4_a, aes(colour="Slightly female-biased \nstarting SR (0.6)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_5_a, aes(colour="Strongly femmale-biased \nstarting SR (0.8)"), fun = mean, geom = "line", size = 1) +
  #ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "Proportion Monogamy (biological trait)", x = "Generation", col = "Parameter setting") +
  #ggtitle("Timeseries of average trait frequency \nfor Monogamy across 100 runs (Starting SR Manipulation)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = "right")

## Plots for "Continuous SR" conditions

### Behaviour phenotype plot

cont1_a <- ggplot(data = Modelrun_6_a, aes(y = p_bh, x = generation)) +
  stat_summary(data = Modelrun_7_a, aes(colour="Strongly male-biased \ncontinuous SR (0.2)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_8_a, aes(colour="Slightly male-biased \ncontinuous SR (0.4)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_6_a, aes(colour="Even continuous SR"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_9_a, aes(colour="Slightly female-biased \ncontinuous SR (0.6)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_10_a, aes(colour="Strongly femmale-biased \ncontinuous SR (0.8)"), fun = mean, geom = "line", size = 1) +
  #ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "Proportion Monogamy (behaviour)", x = "Generation", col = "Parameter setting") +
  #ggtitle("Timeseries of average trait frequency \nfor Monogamy across 100 runs (Continuous SR manipulation)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = "right")

### Cultural trait plot

cont2_a <- ggplot(data = Modelrun_6_a, aes(y = p_c, x = generation)) +
  stat_summary(data = Modelrun_7_a, aes(colour="Strongly male-biased \ncontinuous SR (0.2)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_8_a, aes(colour="Slightly male-biased \ncontinuous SR (0.4)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_6_a, aes(colour="Even continuous SR"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_9_a, aes(colour="Slightly female-biased \ncontinuous SR (0.6)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_10_a, aes(colour="Strongly femmale-biased \ncontinuous SR (0.8)"), fun = mean, geom = "line", size = 1) +
  #ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "Proportion Monogamy (cultural trait)", x = "Generation", col = "Parameter setting") +
  #ggtitle("Timeseries of average trait frequency \nfor Monogamy across 100 runs (Continuous SR manipulation)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = "right")

### Biological trait plot

cont3_a <- ggplot(data = Modelrun_6_a, aes(y = p_b, x = generation)) +
  stat_summary(data = Modelrun_7_a, aes(colour="Strongly male-biased \ncontinuous SR (0.2)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_8_a, aes(colour="Slightly male-biased \ncontinuous SR (0.4)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_6_a, aes(colour="Even continuous SR"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_9_a, aes(colour="Slightly female-biased \ncontinuous SR (0.6)"), fun = mean, geom = "line", size = 1) +
  stat_summary(data = Modelrun_10_a, aes(colour="Strongly femmale-biased \ncontinuous SR (0.8)"), fun = mean, geom = "line", size = 1) +
  #ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "Proportion Monogamy (biological trait)", x = "Generation", col = "Parameter setting") +
  #ggtitle("Timeseries of average trait frequency \nfor Monogamy across 100 runs (Continuous SR manipulation)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = "right")


## Combine plots using "grid arrange"

grob1_a <- as_grob(fixed1_a + theme(legend.position ="none"))
grob2_a <- as_grob(fixed2_a + theme(legend.position ="none"))
grob3_a <- as_grob(fixed3_a + theme(legend.position ="none"))

combined_fixed <- grid.arrange(
  grob1_a,
  grob2_a,
  grob3_a,
  nrow = 1,
  top = grid::textGrob("Parameter 1. Starting Sex Ratio\n", x = 0, hjust = 0)
)

grob1c_a <- as_grob(cont1_a + theme(legend.position ="none"))
grob2c_a <- as_grob(cont2_a + theme(legend.position ="none"))
grob3c_a <- as_grob(cont3_a + theme(legend.position ="none"))


combined_cont <- grid.arrange(
  grob1c_a,
  grob2c_a,
  grob3c_a,
  nrow = 1,
  top = grid::textGrob("Parameter 2. Continuous Sex Ratio\n", x = 0, hjust = 0)
)

combined_total <- grid.arrange(
  combined_fixed,
  combined_cont,
  ncol = 1
)

## Alternative plot structure, including legend

combined_fixed <- grid_arrange_shared_legend(
  fixed1_a,
  fixed2_a,
  fixed3_a,
  cont1_a,
  cont2_a,
  cont3_a,
  nrow = 2,
  ncol = 3,
  position = "bottom",
  top = "Timeseries of average trait frequency \nfor Monogamy across 100 runs (Continuous SR manipulation)"
)

## Example regression analysis confirming null results of condition on percentage monogamy trait at model termination

### Prepare data and set the default scenaria (even starting and continuous SR) as the base category

final_a$parameters <- paste(final_a$`Sex ratio variant`,final_a$`Parameter setting`)
final_a_analysis <- final_a[final_a$parameters != "Continuous 0.5",]
final_a_analysis$parameters <- as.factor(final_a_analysis$parameters)

### Set the balanced scenario as the base comparison level

final_a_analysis$parameters <- relevel(final_a_analysis$parameters, ref="Starting 0.5")

### Compute standard linear regression model

nr_model <- lm(scale(p_bh)~ parameters, data=final_a_analysis)

summary(nr_model)

### Assumption checks

# Linearity

plot(nr_model, 1)

# Normality

plot(nr_model, 2)

#Homogeneity of variance

plot(nr_model, 3)

# Outliers

plot(nr_model, 4)
plot(nr_model, 5)


### Create regression output table in APA style using stargazer

regtable1 <- stargazer(nr_model, 
                       ci=TRUE, 
                       single.row = TRUE,
                       title = "Regression output with standardised coefficients. 95% confidence intervals reported in parantheses.",
                       column.labels = c("Proportion monogamy behaviour"), 
                       covariate.labels = c("Constant (SR 0.5)", "Continuous SR - 0.2", "Continuous SR - 0.4", "Continuous SR - 0.6",
                                            "Continuous SR - 0.8", "Starting SR - 0.2", "Starting SR - 0.4", "Starting SR - 0.6", "Starting SR - 0.6"),
                       intercept.bottom = FALSE,
                       report = ('vcsp'),
                       out="test.html")


# Additional analysis: Why does a biased sex ratio favour polygamy in the biased continuous sex ratio scenarios?


## Combine all "continuous" condition data runs into one dataset, adding an identifying column for the parameter setting under which the data was generated

full_continuous_data <- bind_rows(Modelrun_6_a, Modelrun_7_a, Modelrun_8_a, Modelrun_9_a, Modelrun_10_a, .id="sr_parameter")
full_continuous_data$sr_parameter <- as.factor(c(rep("0.5", 50000), rep("0.2", 50000), rep("0.4", 50000), rep("0.6", 50000), rep("0.8", 50000)))

## Subset data to only include first 100 timesteps and check summary stats

full_continuous_data_subset <- full_continuous_data[full_continuous_data$generation < 101,]

summary_stats <- full_continuous_data_subset %>% 
  group_by(sr_parameter) %>% 
  summarise_at(vars(-generation, -run), funs(mean(., na.rm=T), se = sd(., na.rm=T)/sqrt(sum(!is.na(.)))))

summary_stats


## Create two cross-sectional dataset, only including data from a single timestep (as noted in the dissertation, this prevents autocorrelation issues)
## N.B. in the dissertation, we only used the output data for timestep 1

full_continuous_data_subset_singlegen <- full_continuous_data[full_continuous_data$generation==1,]
full_continuous_data_subset_singlegen2 <- full_continuous_data[full_continuous_data$generation==20,]

## Combine data to plot in ggplot

subset_mono <- select(full_continuous_data_subset_singlegen, -RS_poly)  %>%
  add_column(Behaviour = "Monogamous")

subset_mono <- dplyr::rename(subset_mono, RS = RS_mono)

subset_poly <- select(full_continuous_data_subset_singlegen, -RS_mono)  %>%
  add_column(Behaviour = "Polygamous")

subset_poly <- dplyr::rename(subset_poly, RS = RS_poly)

subset_new <- bind_rows(subset_mono, subset_poly)

mono_mean <- mean(subset_mono$RS)
poly_mean <- mean(subset_poly$RS)

## Create plot

ggplot(subset_new, aes(x=sr_parameter, y=RS, fill=Behaviour)) + 
  geom_hline(aes(yintercept=mono_mean, linetype="Monogamous"), colour="red") +
  geom_hline(aes(yintercept=poly_mean, linetype="Polygamos"), colour="cyan") +
  geom_boxplot(alpha=0.7, position=position_dodge(.9)) +
  geom_violin(alpha=0.2, position=position_dodge(.9)) +
  stat_summary(aes(group=Behaviour), fun.y=mean, geom="point", shape=20, size=5, position=position_dodge(.9)) + 
  theme_bw() +
  labs(y = "Number of offspring (Reproductive success)", x = "Sex ratio (% female)") +
  #ggtitle("Reproductive success in timestep 1") +
  scale_linetype_manual(name = "Cross-parameter \naverage", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("red", "cyan")))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))


## Confirm visual inspection using a regression model

### Prepare data for regression analysis by computing fitness differential variable

subset_modeldata <- full_continuous_data_subset_singlegen %>%
  add_column(RS_diff = full_continuous_data_subset_singlegen$RS_poly - full_continuous_data_subset_singlegen$RS_mono)

### Set 0.5 as the reference level for linear model

subset_modeldata$sr_parameter <- relevel(subset_modeldata$sr_parameter, ref="0.5")

## Run base model 

model1 <- lm(RS_diff ~ sr_parameter, data = subset_modeldata)
aov1 <- aov(RS_diff ~ sr_parameter, data = subset_modeldata)
summary(model1)
summary(aov1)

## Assumption checks

# Linearity

plot(model1, 1)

# Normality

plot(model1, 2)

#Homogeneity of variance

plot(model1, 3)
car::ncvTest(model1)

# Outliers

plot(model1, 4)
plot(model1, 5)

## Based on the slightly unusual q-q plot, we elect to run a "Robust regression" from the robustbase package instead

## Robust regression model

robmodel1 <- lmrob(scale(RS_diff) ~ sr_parameter, data = subset_modeldata)
summary(robmodel1)

# To get F-statistic

robmodel0 <- lmrob(scale(RS_diff) ~ 1, data = subset_modeldata)
anova(robmodel0, robmodel1)


# Polygyny / polyandry analysis - to showcase the divergent patterns of polygyny and polyandry under different sex ratio conditions

## Create dataset for plotting in ggplot

full_continuous_data_subset_singlegen <- as_tibble(full_continuous_data_subset_singlegen)

subset_poly_m <- dplyr::select(full_continuous_data_subset_singlegen, -RS_poly_f, -p_polya)  %>%
  add_column(Sex = "Male")

subset_poly_m <- dplyr::rename(subset_poly_m, RS = RS_poly_m, p_poly = p_polyg)

subset_poly_f <- dplyr::select(full_continuous_data_subset_singlegen, -RS_poly_m, -p_polyg)  %>%
  add_column(Sex = "Female")

subset_poly_f <- dplyr::rename(subset_poly_f, RS = RS_poly_f, p_poly = p_polya)

subset_poly_new <- bind_rows(subset_poly_m, subset_poly_f)

poly_mean_m <- mean(subset_poly_m$RS)
poly_mean_f <- mean(subset_poly_f$RS)

## Create boxplots for both average reproductive success and average proportion of polygamy by sex

RS_plot <- ggplot(subset_poly_new, aes(x=sr_parameter, y=RS, fill=Sex)) + 
  geom_boxplot(alpha=0.7, position=position_dodge(.9), outlier.shape = NA) +
  stat_summary(aes(group=Sex), fun.y=mean, geom="point", shape=20, size=5, position=position_dodge(.9)) + 
  theme_bw() +
  ylim(c(0, 10)) +
  labs(y = "Number of offspring (Reproductive success)", x = "Sex ratio (% female)") +
  ggtitle("Reproductive success in timestep 1") +
  scale_linetype_manual(name = "Cross-parameter \naverage", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("red", "cyan")))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  scale_fill_brewer(palette="Pastel2")

RS_plot

proportion_plot <- ggplot(subset_poly_new, aes(x=sr_parameter, y=p_poly, fill=Sex)) + 
  geom_boxplot(alpha=0.7, position=position_dodge(.9), outlier.shape = NA) +
  stat_summary(aes(group=Sex), fun.y=mean, geom="point", shape=20, size=5, position=position_dodge(.9)) + 
  theme_bw() +
  ylim(c(0, 1)) +
  labs(y = "Proportion polygamous", x = "Sex ratio (% female)") +
  ggtitle("Proportion polygamous in timestep 1") +
  scale_linetype_manual(name = "Cross-parameter \naverage", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("red", "cyan")))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  scale_fill_brewer(palette="Pastel2")

proportion_plot

## Combine plots for use in dissertation

grob1_poly <- as_grob(RS_plot + theme(legend.position ="none"))
grob2_poly <- as_grob(proportion_plot)

combined_poly <- grid.arrange(
  grob1_poly,
  grob2_poly,
  nrow = 1,
  widths=c(0.44, 0.56),
  top = grid::textGrob("", x = 0, hjust = 0)
)


