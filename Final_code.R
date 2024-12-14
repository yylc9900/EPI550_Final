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
library(table1)
library(survival)
library(survminer)
library(ggcorrplot)
library(dplyr)
library(mgcv)
library(tidymv)
library(haven)
library(multcomp)

##### Load clean dataset
getwd()
setwd("E:\\1A-Emory\\#3 FA 24\\RA\\Dataset")
Final550 <- read_sas("E:/1A-Emory/#3 FA 24/EPI550/550 Final/fham.sas7bdat")


###########################################################################################
##### Part A: Table1 (DEATH, SEX, BMI, AGE_BL, educ, CURSMOKE) stratified by HYPERTEN_BL
###########################################################################################

Final550_PERIOD1 <- subset(Final550, PERIOD==1)

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
  factor(Final550_PERIOD1$HYPERTEN,
         levels=c(0,1),
         labels=c("Non-Hypertensive", "Hypertensive"))

#Create Table 1
Final550_PERIOD1 |> 
  select(DEATH, SEX, BMI, AGE_BL, educ, CURSMOKE, HYPERTEN_BL) |>
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
  modify_header(label = "Demographic characteristics") |>
  modify_caption("Table1: Description of demographics in the Framingham Heart Study, stratified by hypertension status at baseline") |>
  modify_footnote(
    all_stat_cols() ~ "*Values may not sum to the total due to missing data; *Mean (SD); *n (%)") |>
  bold_labels()


###########################################################################################
#### Part B: Crude Association
###########################################################################################
#B1: Kaplan-Meier curves
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

#Group by 2: Dichotomize at the median
Final550$AGE_2cat <- ifelse(Final550$AGE<=49, 1, 2)

#Create dummy variables
Final550$group <- with(Final550, ifelse(
  SEX==1 & educ==1,1, ifelse(
    SEX==1 & educ==2,2, ifelse(
      SEX==1 & educ==3,3, ifelse(
        SEX==1 & educ==4,4, ifelse(
          SEX==2 & educ==1,5, ifelse(
            SEX==2 & educ==2,6, ifelse(
              SEX==2 & educ==3,7, ifelse(
                SEX==2 & educ==4,8, NA)))))))))

#Initialize dummy variables to zero
Final550[,c("Z1","Z2","Z3","Z4","Z5","Z6","Z7","Z8")] <- 0
Final550$Z1[Final550$group==1] <- 1
Final550$Z2[Final550$group==2] <- 1
Final550$Z3[Final550$group==3] <- 1
Final550$Z4[Final550$group==4] <- 1
Final550$Z5[Final550$group==5] <- 1
Final550$Z6[Final550$group==6] <- 1
Final550$Z7[Final550$group==7] <- 1
Final550$Z8[Final550$group==8] <- 1


# Full Model
fit_ModelA <- coxph(Surv(TIMEDTH,DEATH) ~ HYPERTEN_BL + BMI + CURSMOKE + tt(HYPERTEN_BL) + 
                      strata(AGE_2cat) + strata(Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8),
              cluster = RANDID, ties="breslow", data=Final550,
              tt = function(x,t, ...) {
                gt1 <- ifelse(t >= 2000 & t < 4000, 1, 0)
                gt2 <- ifelse(t >= 4000 & t < 6000, 1, 0)
                gt3 <- ifelse(t >= 6000, 1, 0)
                cbind(HYPERTEN_BL1=as.numeric(x) * gt1,
                      HYPERTEN_BL2=as.numeric(x) * gt2,
                      HYPERTEN_BL3=as.numeric(x) * gt3)
              })

summary(fit_ModelA)

K_ModelA <- rbind("HR - 2000 ≤t <4000 days" = c(1,0,0,1,0,0,0),
                  "HR - 4000 ≤t <6000 days" = c(1,0,0,0,1,0,0),
                  "HR - t ≥6000 days" = c(0,0,0,0,0,0,0))

contrasts1 <- glht(fit1,linfct = K1)
cbind(exp(coef(contrasts1)), exp(confint.default(contrasts1)))


#The period from 0 – 1,999 days of follow-up 
summary(fit_ModelA)
#The period from 2,000 – 3,999 days of follow-up 
#The period from 4,000 – 5,999 days of follow-up 
#The period from 6,000 days of follow-up through the end of the study period 





#WAY2 (Wrong, can be deleted)
#The period from 0 – 1,999 days of follow-up 
fit1 <- coxph(Surv(TIMEDTH,DEATH) ~ HYPERTEN_BL + educ + tt(BMI) + tt(CURSMOKE) + tt(HYPERTEN_BL) + strata(AGE_2cat) + strata(SEX),
              cluster = RANDID, ties="breslow", data=Final550,
              tt = function(x,t, ...) {
                gt <- ifelse (t < 1999,1,0)
                as.numeric(x) * gt
              })
K1 <- rbind("HR - 0 to 1,999 days" = c(1,0,0,0,1),
            "HR - After 1,999 day" = c(1,0,0,0,0))
contrasts1 <- glht(fit1,linfct = K1)
cbind(exp(coef(contrasts1)), exp(confint.default(contrasts1)))

summary(fit1)




