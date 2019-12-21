library(haven)
library(tidyverse)
library(survival)
library(survminer)
library(flexsurv)

setwd("C:/Users/hirsh/OneDrive/Desktop/NC State/02 Courses/02 Fall/AA502-Analytics Methods & Applications I/Fall 3/Survival Analysis/Homework")
hurricane <- read_sas("hurricane.sas7bdat")

# Remove H1-H48 variables
hurricane2 <- hurricane[, c(1:8,57:59)]
hurricane2$index <- row.names(hurricane2)

hurricane_surv <- Surv(time = hurricane2$hour, event = hurricane2$reason == 1)

# Plot Distributions
aft.mod.llog <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                            slope + elevation, data = hurricane2, dist = 'llogis')
plot(aft.mod.llog, type = "cumhaz", ci = TRUE, conf.int = FALSE, las = 1, bty = "n",
     xlab = "week", ylab = "Cumulative Hazard", main = "Log-Log Distribution")

aft.mod.gamma <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                            slope + elevation, data = hurricane2, dist = 'gamma')
plot(aft.mod.gamma, type = "cumhaz", ci = TRUE, conf.int = FALSE, las = 1, bty = "n",
     xlab = "week", ylab = "Cumulative Hazard", main = "Gamma Distribution")

aft.mod.weibull <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                            slope + elevation, data = hurricane2, dist = 'weibull')
plot(aft.mod.weibull, type = "cumhaz", ci = TRUE, conf.int = FALSE, las = 1, bty = "n",
     xlab = "week", ylab = "Cumulative Hazard", main = "Weibull Distribution")

aft.mod.lnorm <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                            slope + elevation, data = hurricane2, dist = 'lnorm')
plot(aft.mod.lnorm, type = "cumhaz", ci = TRUE, conf.int = FALSE, las = 1, bty = "n",
     xlab = "week", ylab = "Cumulative Hazard", main = "Log-Normal Distribution")

aft.mod.exp <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                            slope + elevation, data = hurricane2, dist = 'exp')
plot(aft.mod.exp, type = "cumhaz", ci = TRUE, conf.int = FALSE, las = 1, bty = "n",
     xlab = "week", ylab = "Cumulative Hazard", main = "Exponential Distribution")


# Test Distributions
loglik.llog <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                              slope + elevation, data = hurricane2, dist = 'llogis')$loglik

loglik.gamma <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                               slope + elevation, data = hurricane2, dist = 'gamma')$loglik

loglik.weibull <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                                 slope + elevation, data = hurricane2, dist = 'weibull')$loglik

loglik.lnorm <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                               slope + elevation, data = hurricane2, dist = 'lnorm')$loglik

loglik.exp <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                             slope + elevation, data = hurricane2, dist = 'exp')$loglik

loglik.genf <- flexsurvreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                            slope + elevation, data = hurricane2, dist = 'genf')$loglik
# gen f does not converge, exclude from analysis

pval.e_g <- 1 - pchisq((-2*(loglik.exp-loglik.gamma)), 2) #  gamma
pval.e_w <- 1 - pchisq((-2*(loglik.exp-loglik.weibull)), 1) #  weibull
pval.w_g <- 1 - pchisq((-2*(loglik.weibull-loglik.gamma)), 1) # weibull
pval.ln_g <- 1 - pchisq((-2*(loglik.lnorm-loglik.gamma)), 1) # gamma
# go with weibull

tests <- c('exp vs gamma', 'exp vs weibull', 'weibull vs gamma', 'log-normal vs gamma')
pvals <- c(pval.e_g, pval.e_w, pval.w_g, pval.ln_g)
dist.tests <- data.frame(tests, pvals)

aft.mod.gamma <- survreg(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
                               slope + elevation, data = hurricane2, dist = 'weibull')


# install.packages("pec")
library(pec)
selectCox(hurricane_surv ~ backup + age + bridgecrane + servo + gear + trashrack +
            slope + elevation, data = hurricane2, rule = 'aic') # backup, servo, slope significant @ p=0.03

aft.mod.gamma2 <- survreg(hurricane_surv ~ servo + slope + backup, data = hurricane2, dist = 'weibull')


# Predicted survival probabilities and time
survprob.actual <- 1 - psurvreg(hurricane2$hour,
                                mean = predict(aft.mod.gamma2, type = "lp"),
                                scale = aft.mod.gamma2$scale,
                                distribution = aft.mod.gamma2$dist)

survprob.10wk <- 1 - psurvreg(10,
                              mean = predict(aft.mod.gamma2, type = "lp"),
                              scale = aft.mod.gamma2$scale,
                              distribution = aft.mod.gamma2$dist)


new_time_servo <-  qsurvreg(1 - survprob.actual,
                      mean = predict(aft.mod.gamma2, type = "lp") + coef(aft.mod.gamma2)['servo'],
                      scale = aft.mod.gamma2$scale,
                      distribution = aft.mod.gamma2$dist)

new_time_backup <-  qsurvreg(1 - survprob.actual,
                            mean = predict(aft.mod.gamma2, type = "lp") + coef(aft.mod.gamma2)['backup'],
                            scale = aft.mod.gamma2$scale,
                            distribution = aft.mod.gamma2$dist)

mean(new_time_backup - hurricane2$hour)
mean(new_time_servo - hurricane2$hour)


# $2.5M to spend on pump upgrades
# $100k per backup upgrade
# $150k per servo upgrade

hurricane2$add_servo <- new_time_servo
hurricane2$add_backup <- new_time_backup
hurricane2$diff_servo <- new_time_servo - hurricane2$hour
hurricane2$diff_backup <- new_time_backup - hurricane2$hour

servo_only_pumps <- hurricane2[hurricane2$hour < 48 & hurricane2$servo == 0 & hurricane2$backup == 1 & hurricane2$reason == 1,]
servo_only_pumps$add_backup <- NULL
servo_only_pumps$diff_time <- servo_only_pumps$add_servo - servo_only_pumps$hour
servo_only_pumps <- servo_only_pumps[servo_only_pumps$add_servo > 48,]

backup_only_pumps <- hurricane2[hurricane2$hour < 48 & hurricane2$backup == 0 & hurricane2$servo == 1 & hurricane2$reason == 1,]
backup_only_pumps$add_servo <- NULL
backup_only_pumps$diff_time <- backup_only_pumps$add_backup - backup_only_pumps$hour
backup_only_pumps <- backup_only_pumps[backup_only_pumps$add_backup > 48,]

both_pumps <- hurricane2[hurricane2$hour < 48 & hurricane2$backup == 0 & hurricane2$servo == 0 & hurricane2$reason == 1 ,]
both_pumps <- both_pumps[both_pumps$add_backup > 48 | both_pumps$add_servo > 48,] # upgrade 1 servo/7 backups

(nrow(servo_pumps) * 150000) + (nrow(backup_pumps) * 100000) # upgrading all remaining pumps costs $4.15M
