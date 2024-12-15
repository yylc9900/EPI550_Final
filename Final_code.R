##########################################################
##########################################################
## EPI550 Final 
##########################################################
##########################################################

##### Load package
library(haven)
library(gtsummary)
library(labelled)
library(ggplot2)
library(survival)
library(survminer)
library(ggcorrplot)
library(dplyr)
library(mgcv)
library(tidymv)
library(multcomp)
library(sandwich)


##### Load clean dataset
#getwd()

path <- "C:/Users/21055/OneDrive - Emory/EPI550_final"
Final550 <- read_sas(paste0(path,"/fham.sas7bdat"))


# Preview of data 
summary(Final550)# BMI, educ have NA's

###########################################################################################
##### Part A: Table1 (DEATH, SEX, BMI, AGE_BL, educ, CURSMOKE) stratified by HYPERTEN_BL
###########################################################################################

Final550_PERIOD1 <- subset(Final550, PERIOD==1)# baseline visit

#Labeling each level
Final550_PERIOD1$DEATH <-
  factor(Final550_PERIOD1$DEATH,
         levels=c(0,1),
         labels=c("Alive", "Died"))

Final550_PERIOD1$SEX <-
  factor(Final550_PERIOD1$SEX,
         levels=c(1,2),
         labels=c("Men", "Women"))

Final550_PERIOD1$educ <-
  factor(Final550_PERIOD1$educ,
         levels=c(1,2,3,4),
         labels=c("0-11 years", "High School Diploma, GED", "Some College, Vocational School","College (BS,BA) degree or more"))

Final550_PERIOD1$CURSMOKE <-
  factor(Final550_PERIOD1$CURSMOKE,
         levels=c(0,1),
         labels=c("Not current smoker", "Current smoker"))

Final550_PERIOD1$HYPERTEN_BL <-
  factor(Final550_PERIOD1$HYPERTEN_BL,
         levels=c(0,1),
         labels=c("Non-Hypertensive", "Hypertensive"))

#Create Table 1
Final550_PERIOD1 |> 
  dplyr::select(DEATH, SEX, BMI, AGE_BL, educ, CURSMOKE, HYPERTEN_BL) |>
  tbl_summary(
    by = HYPERTEN_BL,
    missing = "no",
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous()  ~ "{mean} ({sd})",
                     all_categorical() ~ "{n}    ({p}%)"),
    digits = list(all_continuous()  ~ c(1, 1),
                  all_categorical() ~ c(0, 0)),
    label  = list(DEATH ~ "Death from any cause (outcome)",
                  SEX ~ "Participant sex",
                  BMI ~ "BMI",
                  AGE_BL ~ "Age at baseline (years)",
                  educ ~ "Attained Education",
                  CURSMOKE ~ "Current cigarette smoking at exam")
  ) |>
  add_overall() |>
  modify_header(label = "Characteristics") |>
  modify_caption("Table1: Description of population characterisitics in the Framingham Heart Study, stratified by hypertension status at baseline") |>
  modify_footnote(
    all_stat_cols() ~ "*Values may not sum to the total due to missing data; *Mean (SD); *n (%)") |>
  bold_labels()


###########################################################################################
#### Part B: Crude Association
###########################################################################################
#B1: Kaplan-Meier curves - whole dataset
fit <- survfit(Surv(TIMEDTH,DEATH) ~ HYPERTEN_BL, cluster = RANDID, data=Final550)
ggsurvplot(fit) + ggtitle ("Hypertension at Enrollment Survival Estimates")

#B2:log-rank test
survdiff(formula= Surv(TIMEDTH,DEATH) ~ HYPERTEN_BL, data=Final550)

###########################################################################################
#### Part C: Proportional Hazards Assumption
###########################################################################################

######################################
#### Age at baseline 
######################################

####### Log-log survival curves
#Group by 2: Dichotomize at the median
median(Final550_PERIOD1$AGE_BL) #49
Final550$AGE_2cat <- ifelse(Final550$AGE_BL<=49, 1, 2)

#Generate plots
model_AGE2_log <- survfit(Surv(TIMEDTH,DEATH) ~ AGE_2cat, cluster = RANDID, data=Final550)
ggsurvplot(model_AGE2_log , fun = "cloglog") + 
  ggtitle("PH assumption for AGE(2 groups): Log-log survival curves")

#Group by 3
tertiles <- quantile(Final550_PERIOD1$AGE_BL, probs = c(1/3, 2/3), na.rm = TRUE)
tertiles # 45,54

Final550$AGE_3cat <- ifelse(Final550$AGE_BL <= 45, "1", 
                            ifelse(Final550$AGE_BL <= 54, "2", "3"))
#Generate plots
model_AGE3_log <- survfit(Surv(TIMEDTH,DEATH) ~ AGE_3cat, cluster = RANDID, data = Final550)
ggsurvplot(model_AGE3_log , fun = "cloglog") + 
  ggtitle("PH assumption for AGE(3 groups): Log-log survival curves")


######### Goodness-of-fit
model_AGE_GOF <- coxph(Surv(TIMEDTH,DEATH) ~ AGE_BL, cluster = RANDID, data = Final550, ties = "breslow")
cox.zph(model_AGE_GOF)

######### Time-dependent covariates 
#g(t) = t
coxph(Surv(TIMEDTH,DEATH) ~ AGE_BL + tt(AGE_BL), cluster = RANDID, data = Final550, ties = "breslow",
      tt = function(x,t,...) as.numeric(x)*t
)
#g(t) = ln(t)
coxph(Surv(TIMEDTH,DEATH) ~ AGE_BL + tt(AGE_BL),cluster = RANDID, data = Final550, ties = "breslow",
      tt = function(x,t,...) as.numeric(x)*log(t)
)


######################################
#### Education 
######################################

####### Log-log survival curves
model_educ_log <- survfit(Surv(TIMEDTH,DEATH) ~ educ, cluster = RANDID, data = Final550)
ggsurvplot(model_educ_log , fun = "cloglog") + 
  ggtitle("PH assumption for Education: Log-log survival curves")

######### Goodness-of-fit
model_educ_GOF <- coxph(Surv(TIMEDTH,DEATH) ~ educ, cluster = RANDID, data = Final550, ties = "breslow")
cox.zph(model_educ_GOF)

######### Time-dependent covariates
#g(t) = t
coxph(Surv(TIMEDTH,DEATH) ~ educ + tt(educ),cluster = RANDID, data = Final550, ties = "breslow",
      tt = function(x,t,...) as.numeric(x)*t
)
#g(t) = ln(t)
coxph(Surv(TIMEDTH,DEATH) ~ educ + tt(educ),cluster = RANDID, data = Final550,ties = "breslow",
      tt = function(x,t,...) as.numeric(x)*log(t)
)

######################################
## Prevalent hypertension at baseline 
######################################

####### Log-log survival curves
model_Hyper_log <- survfit(Surv(TIMEDTH,DEATH) ~ HYPERTEN_BL, cluster = RANDID, data = Final550)
ggsurvplot(model_Hyper_log , fun = "cloglog") + 
  ggtitle("PH assumption for hypertension at baseline: Log-log survival curves")

######### Goodness-of-fit
model_Hyper_GOF <- coxph(Surv(TIMEDTH,DEATH) ~ HYPERTEN_BL, cluster = RANDID, data = Final550, ties = "breslow")
cox.zph(model_Hyper_GOF)

######### Time-dependent covariates
#g(t) = t
coxph(Surv(TIMEDTH,DEATH) ~ HYPERTEN_BL + tt(HYPERTEN_BL),cluster = RANDID, data = Final550, ties = "breslow",
      tt = function(x,t,...) as.numeric(x)*t
)
#g(t) = ln(t)
coxph(Surv(TIMEDTH,DEATH) ~ HYPERTEN_BL + tt(HYPERTEN_BL),cluster = RANDID, data = Final550,ties = "breslow",
      tt = function(x,t,...) as.numeric(x)*log(t)
)


######################################
## Sex 
######################################

####### Log-log survival curves
model_SEX_log <- survfit(Surv(TIMEDTH,DEATH) ~ SEX, cluster = RANDID, data = Final550)
ggsurvplot(model_SEX_log , fun = "cloglog") + 
  ggtitle("PH assumption for SEX: Log-log survival curves")

######### Goodness-of-fit
model_SEX_GOF <- coxph(Surv(TIMEDTH,DEATH) ~ SEX, cluster = RANDID, data = Final550, ties = "breslow")
cox.zph(model_SEX_GOF)

######### Time-dependent covariates
#g(t) = t
coxph(Surv(TIMEDTH,DEATH) ~ SEX + tt(SEX),cluster = RANDID, data = Final550, ties = "breslow",
      tt = function(x,t,...) as.numeric(x)*t
)
#g(t) = ln(t)
coxph(Surv(TIMEDTH,DEATH) ~ SEX + tt(SEX),cluster = RANDID, data = Final550,ties = "breslow",
      tt = function(x,t,...) as.numeric(x)*log(t)
)


###########################################################################################
#### Part D: SURVIVAL ANALYSIS for Question A
###########################################################################################

fit_ModelA <- coxph(Surv(TIMEDTH,DEATH) ~ HYPERTEN_BL + AGE_BL + BMI + CURSMOKE + tt(HYPERTEN_BL) + 
                      strata(SEX) + strata(educ),
                    cluster = RANDID, ties="breslow", data=Final550,
                    tt = function(x,t, ...) {
                      gt1 <- ifelse(t >= 2000 & t < 4000, 1, 0)
                      gt2 <- ifelse(t >= 4000 & t < 6000, 1, 0)
                      gt3 <- ifelse(t >= 6000, 1, 0)
                      cbind(HYPERTEN_BL1=as.numeric(x) * gt1,
                            HYPERTEN_BL2=as.numeric(x) * gt2,
                            HYPERTEN_BL3=as.numeric(x) * gt3)
                    })

# HR (95% CI) for 0-2000 days - beta1
summary(fit_ModelA)

K_ModelA <- rbind("HR - 2000 ≤t <4000 days" = c(1,0,0,0,1,0,0),
                  "HR - 4000 ≤t <6000 days" = c(1,0,0,0,0,1,0),
                  "HR - t ≥6000 days" = c(1,0,0,0,0,0,1))

contrasts1 <- glht(fit_ModelA,linfct = K_ModelA)

# HRs (95% CI) for rest three time ranges
cbind(exp(coef(contrasts1)), exp(confint.default(contrasts1)))


###########################################################################################
#### Part E: OTHER REGRESSION TYPES for Question A
###########################################################################################

### Restrict this analysis to the based line visit - Final550_PERIOD1
table(Final550_PERIOD1$DEATH) # dichotomous outcome

## First, we use unconditional logistic regression to get OR

fit_log <- glm(data = Final550_PERIOD1,
               DEATH ~ HYPERTEN_BL + AGE_BL + BMI + CURSMOKE + AGE_BL + SEX + educ,
               family=binomial (link='logit'))
summary(fit_log)

cbind(exp(coef(fit_log)), exp(confint.default(fit_log)))

## Next, we use poisson regression to get RR
# Turn outcome variable back to numeric to run the model
if (is.factor(Final550_PERIOD1$DEATH)) {
  # Convert DEATH to numeric, assuming it has levels like "0" and "1"
  Final550_PERIOD1$DEATH <- as.numeric(as.character(Final550_PERIOD1$DEATH))
}

prop.table(table(Final550_PERIOD1$DEATH))
is.numeric(Final550_PERIOD1$DEATH)

# fit the model with robust variance estimates
fit_poi <- glm(data = Final550_PERIOD1,
               DEATH ~ HYPERTEN_BL + AGE_BL + BMI + CURSMOKE + AGE_BL + SEX + educ,
               family = 'poisson')

# Calculate 95% CI
sandwich_se <- diag(vcovHC(fit_poi, type = "HC"))^0.5

lowerlim <- coef(fit_poi) - 1.96*sandwich_se
upperlim <- coef(fit_poi) + 1.96*sandwich_se

PrevRatio <- coef(fit_poi)
cbind(exp(PrevRatio), exp(lowerlim), exp(upperlim))



###########################################################################################
##### Part F: Survival analysis for CHD diagnosis: BMI + CURSMOKE + BMI*SEX
###########################################################################################
Final550_PREVCHD0 <- Final550[Final550$PREVCHD==0,] # N = 10785

# Outcome is the incident CHD diagnosis
table(Final550_PREVCHD0$ANYCHD)

# BMI recoded to nominal variable
Final550_PREVCHD0$BMI_OW <- ifelse(is.na(Final550_PREVCHD0$BMI),NA,
                                   ifelse(Final550_PREVCHD0$BMI>=25 & Final550_PREVCHD0$BMI < 30,1,0))
Final550_PREVCHD0$BMI_OB <- ifelse(is.na(Final550_PREVCHD0$BMI),NA,
                                   ifelse(Final550_PREVCHD0$BMI>=30,1,0))

### Fit the model
fit_modelE <- coxph(Surv(TIMECHD,ANYCHD) ~ BMI_OW + BMI_OB + CURSMOKE + (BMI_OW + BMI_OB + CURSMOKE):SEX + 
                      strata(SEX),
                    data = Final550_PREVCHD0,
                    ties = "breslow")

k_modelE <- rbind('Overweight, males'= c(1,0,0,1,0,0),
                  'Overweight, females'= c(1,0,0,2,0,0),
                  'Obesity, males' = c(0,1,0,0,1,0),
                  'Obseity, females'= c(0,1,0,0,2,0))

contrasts2 <- glht(fit_modelE,linfct = k_modelE)

# HRs (95% CI) for each stratum
cbind(exp(coef(contrasts2)), exp(confint.default(contrasts2)))

### LRT - statistical interaction?
## Full model - fit_modelE
## Reduced model
fit_modelE_reduce <- coxph(Surv(TIMECHD,ANYCHD) ~ BMI_OW + BMI_OB + CURSMOKE + strata(SEX),
                    data = Final550_PREVCHD0,
                    ties = "breslow")
## Test statistic
anova(fit_modelE_reduce, fit_modelE, test = 'chisq')

### Confounding by current smoke?
## Gold standard model - fit_modelE
## Exclude CURSMOKE
fit_modelE_1 <- coxph(Surv(TIMECHD,ANYCHD) ~ BMI_OW + BMI_OB + (BMI_OW + BMI_OB):SEX + 
                      strata(SEX),
                    data = Final550_PREVCHD0,
                    ties = "breslow")
# New HRs
k_modelE_1 <- rbind('Overweight, males'= c(1,0,1,0),
                    'Overweight, females'= c(1,0,2,0),
                    'Obesity, males' = c(0,1,0,1),
                    'Obseity, females'= c(0,1,0,2))

contrasts3 <- glht(fit_modelE_1,linfct = k_modelE_1)
exp(coef(contrasts3))

