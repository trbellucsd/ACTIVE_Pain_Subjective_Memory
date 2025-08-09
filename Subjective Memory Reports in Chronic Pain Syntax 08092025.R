# Load required packages
library(OpenMx)
library(psych)
library(devtools)
library(umx)
library(haven)
library(dplyr)
library(car)
library(lme4)
library(nlme)
library(lmerTest)
library(ggplot2)
library(tidyverse)
library(multilevelTools)
library(JWileymisc)
library(effects)
library(sjPlot)
library(sjmisc)
library(reghelper)
library(Rcmdr)
library(emmeans)
library(interactions)
#library(ReporteRs)
library(expss)
library(tidyr)
library(table1)

setwd("Z:/Tyler/Research/aa_Research Projects and Publications/ACTIVE/Data")

ACTIVE <- read_sav("ACTIVE_DATA_Long_23September2022.sav",encoding="Latin2")
names(ACTIVE)
names(ACTIVE)<- toupper(names(ACTIVE))
# Computing chronic pain status 

ACTIVE$BP_1 <- ifelse(ACTIVE$SF3621_1 >= 0 & ACTIVE$SF3621_1 <= 3, 0, 
                      ifelse(ACTIVE$SF3621_1 >= 4, 1, NA))
ACTIVE$BP_2 <- ifelse(ACTIVE$SF3621_2 >= 0 & ACTIVE$SF3621_2 <= 3, 0, 
                      ifelse(ACTIVE$SF3621_2 >= 4, 1, NA))

ACTIVE$BP_SUM_1 <- ACTIVE$BP_1 + ACTIVE$BP_2

ACTIVE$CP <- ifelse(ACTIVE$BP_SUM_1 == 2, 1, 0)

ACTIVE$CP <- factor(ACTIVE$CP, levels = c(0, 1), labels = c("CP-", "CP+"))

# Using complete.cases
ACTIVE_no_missing <- ACTIVE[complete.cases(ACTIVE[, c("CP", "MFQFF")]), ]



####

#Speed
SOPCOMP.mlm <-lmer(scale(SOPCOMP_R) ~   scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                     as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                     as_factor(CP)*scale(AGE_LONG_TIMESINCE)*scale(MFQ_FF_B) + 
                     scale(MFQFF_WP)*as_factor(CP) + (1| ID ), ACTIVE, REML=FALSE) 
summary(SOPCOMP.mlm)

SOPCOMP.mlm.tab<-tab_model(SOPCOMP.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)
SOPCOMP.mlm.tab

#Memory
MEM_COMP.mlm <-lmer(scale(MEM_COMP) ~   scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                      as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                      as_factor(CP)*scale(AGE_LONG_TIMESINCE)*scale(MFQ_FF_B) + 
                      scale(MFQFF_WP)*as_factor(CP) + (1| ID ), ACTIVE, REML=FALSE) 
summary(MEM_COMP.mlm)

MEM_COMP.mlm.tab<-tab_model(MEM_COMP.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)

MEM_COMP.mlm.tab

#Reasoning
REASCOMP.mlm <-lmer(scale(REASCOMP) ~ scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                      as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                      as_factor(CHRONICPAIN_RC)*scale(AGE_LONG_TIMESINCE)*scale(MFQ_FF_B) + 
                      scale(MFQFF_WP)*as_factor(CHRONICPAIN_RC) + (1| ID ), ACTIVE, REML=FALSE) 
summary(REASCOMP.mlm)
REASCOMP.mlm.tab<-tab_model(REASCOMP.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)

#CESDTOT
ACTIVE$CESDTOT
CESDTOT.mlm <-lmer(scale(CESDTOT) ~ scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                      as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                      as_factor(CP)*scale(AGE_LONG_TIMESINCE)*scale(MFQ_FF_B) + 
                      scale(MFQFF_WP)*as_factor(CP) + (1| ID ), ACTIVE, REML=FALSE) 
summary(CESDTOT.mlm)
CESDTOT.mlm.tab<-tab_model(CESDTOT.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)

#SF36SF_L
ACTIVE$SF36SF_L
SF36SF_L.mlm <-lmer(scale(SF36SF_L) ~ scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                      as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                      as_factor(CP)*scale(AGE_LONG_TIMESINCE)*scale(MFQ_FF_B) + 
                      scale(MFQFF_WP)*as_factor(CP) + (1| ID ), ACTIVE, REML=FALSE) 
summary(SF36SF_L.mlm)

SF36SF_L.mlm.tab<-tab_model(SF36SF_L.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)

#SF36PF_L
ACTIVE$SF36PF_L
SF36PF_L.mlm <-lmer(scale(SF36PF_L) ~ scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                     as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                     as_factor(CP)*scale(AGE_LONG_TIMESINCE)*scale(MFQ_FF_B) + 
                     scale(MFQFF_WP)*as_factor(CP) + (1| ID ), ACTIVE, REML=FALSE) 
summary(SF36PF_L.mlm)
SF36PF_L.mlm.tab<-tab_model(SF36PF_L.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)

#TURN360
ACTIVE$TURN360
TURN360.mlm <-lmer(scale(TURN360) ~ scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                      as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                      as_factor(CP)*scale(AGE_LONG_TIMESINCE)*scale(MFQ_FF_B) + 
                      scale(MFQFF_WP)*as_factor(CP) + (1| ID ), ACTIVE, REML=FALSE) 
summary(TURN360.mlm)
TURN360.mlm.tab<-tab_model(TURN360.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)



#DYN_AVG
ACTIVE$DYN_AVG
DYN_AVG.mlm <-lmer(scale(DYN_AVG_M) ~ scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                      as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                      as_factor(CP)*scale(AGE_LONG_TIMESINCE)*scale(MFQ_FF_B) + 
                      scale(MFQFF_WP)*as_factor(CP) + (1| ID ), ACTIVE, REML=FALSE) 
summary(DYN_AVG.mlm)
DYN_AVG.mlm.tab<-tab_model(DYN_AVG.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)


#IADL
ACTIVE$IADL
IADL.mlm <-lmer(scale(MDS_IADL_LIMITATIONS) ~ scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                     as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                     as_factor(CP)*scale(AGE_LONG_TIMESINCE)*scale(MFQ_FF_B) + 
                     scale(MFQFF_WP)*as_factor(CP) + (1| ID ), ACTIVE, REML=FALSE) 
summary(IADL.mlm)
IADL.mlm.tab<-tab_model(IADL.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)

ACTIVE$IADL
IADL.mlm <-lmer(scale(MFQ_FF_L_R) ~ scale(AGEB) + as_factor(MARSTAT_CAT) + as_factor(RACE_CAT) + scale(EDUCLEVL) + scale(CESDTOT) + 
                  as_factor(INTGRPR)*scale(AGE_LONG_TIMESINCE) + 
                  as_factor(CP)*scale(AGE_LONG_TIMESINCE) + (1| ID ), ACTIVE, REML=FALSE) 
summary(IADL.mlm)
IADL.mlm.tab<-tab_model(IADL.mlm, p.val = "kr", show.stat=TRUE,show.se = TRUE)

###FIGURES###############################################################################################################

ACTIVE$MFQ_FF_BR
# Create the median split variable
median_MFQ_FF_B <- median(ACTIVE$MFQ_FF_B, na.rm = TRUE)
ACTIVE$MFQ_FF_B_median_split <- ifelse(ACTIVE$MFQ_FF_B >= median_MFQ_FF_B, "High", "Low")
ACTIVE$GROUP <- ifelse(ACTIVE$MFQ_FF_B_median_split == "Low" & ACTIVE$CP == 0, 1,
                       ifelse(ACTIVE$MFQ_FF_B_median_split == "Low" & ACTIVE$CP == 1, 2,
                              ifelse(ACTIVE$MFQ_FF_B_median_split == "High" & ACTIVE$CP == 0, 3,
                                     ifelse(ACTIVE$MFQ_FF_B_median_split == "High" & ACTIVE$CP == 1, 4, NA))))
ACTIVE$GROUP <- factor(ACTIVE$GROUP,
                      levels = c(1, 2, 3, 4),
                      labels = c("CP-/SMP-HIGH", "CP+/SMP-HIGH", "CP-/SMP-LOW", "CP+/SMP-LOW"))


table(ACTIVE$AGE_LONG)      

median_MFQ_FF_B <- median(ACTIVE$MFQ_FF_BR, na.rm = TRUE)
ACTIVE$MFQ_FF_B_median_split_r <- ifelse(ACTIVE$MFQ_FF_BR >= median_MFQ_FF_B, "High", "Low")

#Memory Figure
## Subset the data to include only rows with complete cases for REASCOMP and GROUP
ACTIVE_MEM_sub <- subset(ACTIVE, complete.cases(MEM_COMP, GROUP))

## Fit the linear mixed-effects model
MEM_COMP.fig.mlm <- lmer(MEM_COMP ~ AGEB + MARSTAT_CAT + RACE_CAT + EDUCLEVL + CESDTOT + 
                             INTGRPR + INTGRPR * AGE_LONG_TIMESINCE + (1 | ID), 
                           data = ACTIVE_MEM_sub, REML = FALSE)

## Display the model summary
summary(MEM_COMP.fig.mlm)

## Add residuals to the dataset
ACTIVE_MEM_sub$MEM_COMP_RES <- residuals(MEM_COMP.fig.mlm)



figure.MEM_COMP.mlm <- ggplot(ACTIVE_MEM_sub, aes(x =AGE_LONG_TIMESINCE, y = MEM_COMP_RES, color = as.factor(GROUP), linetype = as.factor(GROUP))) +
  geom_smooth(aes(color = as.factor(GROUP), linetype = as.factor(GROUP)), method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Time Since Baseline", y = "Memory", color = "Group", linetype = "Group") +  # Corrected 'linetype' argument
  theme_minimal() +
  theme(legend.position = "bottom") 
figure.MEM_COMP.mlm

#Reasoning Figure
## Subset the data to include only rows with complete cases for REASCOMP and GROUP
ACTIVE_reas_sub <- subset(ACTIVE, complete.cases(REASCOMP, GROUP))

## Fit the linear mixed-effects model
reasoncomp.fig.mlm <- lmer(REASCOMP ~ AGEB + MARSTAT_CAT + RACE_CAT + EDUCLEVL + CESDTOT + 
                             INTGRPR + INTGRPR * AGE_LONG_TIMESINCE + (1 | ID), 
                           data = ACTIVE_reas_sub, REML = FALSE)

## Display the model summary
summary(reasoncomp.fig.mlm)

## Add residuals to the dataset
ACTIVE_reas_sub$REASCOMP_RES <- residuals(reasoncomp.fig.mlm)

##Create figure 
figure.reasoncomp.mlm <- ggplot(ACTIVE_reas_sub, aes(x =AGE_LONG_TIMESINCE, y = REASCOMP_RES, color = as.factor(GROUP), linetype = as.factor(GROUP))) +
  geom_smooth(aes(color = as.factor(GROUP), linetype = as.factor(GROUP)), method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Years Since Baseline", y = "Reasoning (z-score)", color = "", linetype = "") +  # Corrected 'linetype' argument
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 14),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 14))
figure.reasoncomp.mlm

#SOP Figure
# Subset the data to include only rows with complete cases for REASCOMP and GROUP
ACTIVE_SOPCOMP_sub <- subset(ACTIVE, complete.cases(SOPCOMP, GROUP,AGE_LONG_TIMESINCE,
                                                      CESDTOT,AGEB,MARSTAT_CAT,RACE_CAT,EDUCLEVL,INTGRPR ))
## Fit the linear mixed-effects model
SOPCOMPcomp.fig.mlm <- lmer(SOPCOMP ~ AGEB + MARSTAT_CAT + RACE_CAT + EDUCLEVL + CESDTOT + 
                             INTGRPR + INTGRPR * AGE_LONG_TIMESINCE + (1 | ID), 
                           data = ACTIVE_SOPCOMP_sub, REML = FALSE)

## Display the model summary
summary(SOPCOMPcomp.fig.mlm)

## Add residuals to the dataset
ACTIVE_SOPCOMP_sub$SOPCOMP_RES <- residuals(SOPCOMPcomp.fig.mlm)


figure.ACTIVE_SOPCOMP_R.mlm <- ggplot(ACTIVE_SOPCOMP_sub, aes(x =AGE_LONG_TIMESINCE, y = scale(SOPCOMP), color = as.factor(GROUP), linetype = as.factor(GROUP))) +
  geom_smooth(aes(color = as.factor(GROUP), linetype = as.factor(GROUP)), method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Years Since Baseline", y = "Speed of Processing (z-score)", color = "", linetype = "") +  # Corrected 'linetype' argument
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 14))
figure.ACTIVE_SOPCOMP_R.mlm


#SF36PF_L
## Subset the data ensuring complete cases for the specified variables
ACTIVE_SF36PF_L_sub <- subset(ACTIVE, complete.cases(SF36PF_L, GROUP, AGE_LONG_TIMESINCE,
                                                     CESDTOT, AGEB, MARSTAT_CAT, RACE_CAT, EDUCLEVL, INTGRPR))

## Fit the linear mixed-effects model
SF36PF_L.fig.mlm <- lmer(SF36PF_L ~ AGEB + MARSTAT_CAT + RACE_CAT + EDUCLEVL + CESDTOT + 
                           INTGRPR + INTGRPR * AGE_LONG_TIMESINCE + (1 | ID), 
                         data = ACTIVE_SF36PF_L_sub, REML = FALSE)

## Display the model summary
summary(SF36PF_L.fig.mlm)

## Add residuals to the dataset
ACTIVE_SF36PF_L_sub$SF36PF_L_RES <- residuals(SF36PF_L.fig.mlm)

## Create the figure using ggplot2
figure.ACTIVE_SF36PF_L.mlm <- ggplot(ACTIVE_SF36PF_L_sub, 
                                     aes(x = AGE_LONG_TIMESINCE, y = scale(SF36PF_L), 
                                         color = as.factor(GROUP), linetype = as.factor(GROUP))) +
  geom_smooth(aes(color = as.factor(GROUP), linetype = as.factor(GROUP)), 
              method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Years Since Baseline", y = "SF36 Physical Function (z-score)", color = "", linetype = "") +
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 14)) + theme(legend.position = "none") 

## Display the figure
figure.ACTIVE_SF36PF_L.mlm



#TURN360 - SMP
# Subset the data ensuring complete cases for the specified variables
ACTIVE_TURN360_sub <- subset(ACTIVE, complete.cases(TURN360, GROUP, AGE_LONG_TIMESINCE,
                                                    CESDTOT, AGEB, MARSTAT_CAT, RACE_CAT, EDUCLEVL, INTGRPR))

# Fit the linear mixed-effects model
TURN360.fig.mlm <- lmer(TURN360 ~ AGEB + MARSTAT_CAT + RACE_CAT + EDUCLEVL + CESDTOT + 
                          INTGRPR + INTGRPR * AGE_LONG_TIMESINCE + (1 | ID), 
                        data = ACTIVE_TURN360_sub, REML = FALSE)

# Display the model summary
summary(TURN360.fig.mlm)

# Add residuals to the dataset
ACTIVE_TURN360_sub$TURN360_RES <- residuals(TURN360.fig.mlm)

# Create the figure using ggplot2
figure.ACTIVE_TURN360.mlm <- ggplot(ACTIVE_TURN360_sub, 
                                    aes(x = AGE_LONG_TIMESINCE, y = scale(TURN360), 
                                        color = as.factor(MFQ_FF_B_median_split), linetype = as.factor(MFQ_FF_B_median_split))) +
  geom_smooth(aes(color = as.factor(MFQ_FF_B_median_split), linetype = as.factor(MFQ_FF_B_median_split)), 
              method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Years Since Baseline", y = "Turn360 (z-score)", color = "", linetype = "") +
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 14))

# Display the figure
figure.ACTIVE_TURN360.mlm


#TURN360 - CP
## Subset the data ensuring complete cases for the specified variables
ACTIVE_TURN360_sub <- subset(ACTIVE, complete.cases(TURN360, GROUP, AGE_LONG_TIMESINCE,
                                                    CESDTOT, AGEB, MARSTAT_CAT, RACE_CAT, EDUCLEVL, INTGRPR))

## Fit the linear mixed-effects model
TURN360.fig.mlm <- lmer(TURN360 ~ AGEB + MARSTAT_CAT + RACE_CAT + EDUCLEVL + CESDTOT + 
                          INTGRPR + INTGRPR * AGE_LONG_TIMESINCE + (1 | ID), 
                        data = ACTIVE_TURN360_sub, REML = FALSE)

## Display the model summary
summary(TURN360.fig.mlm)

## Add residuals to the dataset
ACTIVE_TURN360_sub$TURN360_RES <- predict(TURN360.fig.mlm)

## Create the figure using ggplot2
figure.ACTIVE_TURN360_CP.mlm <- ggplot(ACTIVE_TURN360_sub, 
                                    aes(x = AGE_LONG_TIMESINCE, y = scale(TURN360), 
                                        color = as.factor(CP), linetype = as.factor(CP))) +
  geom_smooth(aes(color = as.factor(CP), linetype = as.factor(CP)), 
              method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Time Since Baseline", y = "Turn360 (z-score)", color = "", linetype = "") +
  theme_minimal()  +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 14))

# Display the figure
figure.ACTIVE_TURN360_CP.mlm


#DYN_AVG_M
# Subset the data ensuring complete cases for the specified variables
ACTIVE_DYN_AVG_M_sub <- subset(ACTIVE, complete.cases(DYN_AVG_M, GROUP, AGE_LONG_TIMESINCE,
                                                      CESDTOT, AGEB, MARSTAT_CAT, RACE_CAT, EDUCLEVL, INTGRPR))

# Fit the linear mixed-effects model
DYN_AVG_M.fig.mlm <- lmer(DYN_AVG_M ~ AGEB + MARSTAT_CAT + RACE_CAT + EDUCLEVL + CESDTOT + 
                            INTGRPR + INTGRPR * AGE_LONG_TIMESINCE + (1 | ID), 
                          data = ACTIVE_DYN_AVG_M_sub, REML = FALSE)

# Display the model summary
summary(DYN_AVG_M.fig.mlm)

# Add residuals to the dataset
ACTIVE_DYN_AVG_M_sub$DYN_AVG_M_RES <- residuals(DYN_AVG_M.fig.mlm)

# Create the figure using ggplot2
figure.ACTIVE_DYN_AVG_M.mlm <- ggplot(ACTIVE_DYN_AVG_M_sub, 
                                      aes(x = AGE_LONG_TIMESINCE, y = scale(DYN_AVG_M), 
                                          color = as.factor(MFQ_FF_B_median_split_r), linetype = as.factor(MFQ_FF_B_median_split_r))) +
  geom_smooth(aes(color = as.factor(MFQ_FF_B_median_split_r), linetype = as.factor(MFQ_FF_B_median_split_r)), 
              method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Time Since Baseline", y = "Grip Strength (z-score)", color = "", linetype = "") +
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 14)) + theme(legend.position = "none") 

# Display the figure
figure.ACTIVE_DYN_AVG_M.mlm




#SF36SF_L
## Subset the data ensuring complete cases for the specified variables
ACTIVE_SF36SF_L_sub <- subset(ACTIVE, complete.cases(SF36SF_L, GROUP, AGE_LONG_TIMESINCE,
                                                     CESDTOT, AGEB, MARSTAT_CAT, RACE_CAT, EDUCLEVL, INTGRPR))

## Fit the linear mixed-effects model
SF36SF_L.fig.mlm <- lmer(SF36SF_L ~ AGEB + MARSTAT_CAT + RACE_CAT + EDUCLEVL + CESDTOT + 
                           INTGRPR + INTGRPR * AGE_LONG_TIMESINCE + (1 | ID), 
                         data = ACTIVE_SF36SF_L_sub, REML = FALSE)

## Display the model summary
summary(SF36SF_L.fig.mlm)

## Add residuals to the dataset
ACTIVE_SF36SF_L_sub$SF36SF_L_RES <- residuals(SF36SF_L.fig.mlm)

## Create the figure using ggplot2
figure.ACTIVE_SF36SF_L.mlm <- ggplot(ACTIVE_SF36SF_L_sub, 
                                     aes(x = AGE_LONG_TIMESINCE, y = scale(SF36SF_L_RES), 
                                         color = as.factor(MFQ_FF_B_median_split_r), linetype = as.factor(MFQ_FF_B_median_split_r))) +
  geom_smooth(aes(color = as.factor(MFQ_FF_B_median_split), linetype = as.factor(MFQ_FF_B_median_split_r)), 
              method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Time Since Baseline", y = "SF36 Social Function (z-score)", color = "", linetype = "") +
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 14)) + theme(legend.position = "none") 

## Display the figure
figure.ACTIVE_CP_SMC.mlm

figure.ACTIVE_CP_SMC.mlm <- ggplot(ACTIVE, 
                                     aes(x = AGE_LONG_TIMESINCE, y = scale(MFQ_FF_L_R), 
                                         color = as.factor(CP), linetype = as.factor(MFQ_FF_B_median_split_r))) +
  geom_smooth(aes(color = as.factor(CP), linetype = as.factor(CP)), 
              method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Time Since Baseline", y = "SF36 Social Function (z-score)", color = "", linetype = "") +
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 14)) + theme(legend.position = "none") 

## Display the figure
figure.ACTIVE_CP_SMC.mlm

# Extract the legend from your plot
figure.ACTIVE_DYN_AVG_M.mlm2 <- ggplot(ACTIVE_DYN_AVG_M_sub, 
                                      aes(x = AGE_LONG_TIMESINCE, y = scale(DYN_AVG_M_RES), 
                                          color = as.factor(MFQ_FF_B_median_split_r), linetype = as.factor(MFQ_FF_B_median_split_r))) +
  geom_smooth(aes(color = as.factor(MFQ_FF_B_median_split_r), linetype = as.factor(MFQ_FF_B_median_split_r)), 
              method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Time Since Baseline", y = "Grip Strength (z-score)", color = "Group", linetype = "Group") +
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 12)) 

# Display the figure
figure.ACTIVE_DYN_AVG_M.mlm2


legend_grob <- get_legend(figure.ACTIVE_DYN_AVG_M.mlm2)
plot(legend_grob)

figure.reasoncomp.mlm2 <- ggplot(ACTIVE_reas_sub, aes(x =AGE_LONG_TIMESINCE, y = REASCOMP_RES, color = as.factor(GROUP), linetype = as.factor(GROUP))) +
  geom_smooth(aes(color = as.factor(GROUP), linetype = as.factor(GROUP)), method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Years Since Baseline", y = "Reasoning (z-score)", color = "Group", linetype = "Group") +  # Corrected 'linetype' argument
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 12))
figure.reasoncomp.mlm2

legend_grob2 <- get_legend(figure.reasoncomp.mlm2)
plot(legend_grob2)

figure.ACTIVE_TURN360_CP.mlm2 <- ggplot(ACTIVE_TURN360_sub, 
                                       aes(x = AGE_LONG_TIMESINCE, y = scale(TURN360_RES), 
                                           color = as.factor(CP), linetype = as.factor(CP))) +
  geom_smooth(aes(color = as.factor(CP), linetype = as.factor(CP)), 
              method = "lm", se = TRUE, fullrange = TRUE, size = 1.5) +
  labs(x = "Time Since Baseline", y = "Turn360 (z-score)", color = "CP", linetype = "CP") +
  theme_minimal() + ylim(-.5,.5) +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 12),         
        axis.title = element_text(family = "serif", size = 14, face = "bold"),  
        axis.text = element_text(family = "serif", size = 12))

# Display the figure
figure.ACTIVE_TURN360_CP.mlm2
legend_grob3 <- get_legend(figure.ACTIVE_TURN360_CP.mlm2)
plot(legend_grob3)


