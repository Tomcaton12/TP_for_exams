## ---------------------------
## R version : 4.2.2
## 
## Script name: 
##
## Purpose of script: Inspiration
##
## Author: X
##
## Date Created: 2022-11-15
##
## Copyright (c) X, 2022
## Email: 
##
## ---------------------------
##
## Notes:
##    Sans offset
##    Avec offset
##     
## ---------------------------
##
## load up the packages we will need:
pacman::p_load(tictoc, readr, readxl, tidyverse, glmmTMB, lme4, lmerTest, nnet, 
               broom, broom.mixed, leaps, DataEditR, VGAM, pROC, geeM, MuMIn, 
               purrr, forcats, ResourceSelection)

## ---------------------------
## load dataset
cig <- readRDS("D:/M2 - SMSDS/M2 SMSDS/Modèles linéaires généralisés, modèles mixtes/Cours/TP-20221111/cig.rds")

##%######################################################%##
#                                                          #
####                    Sans offset                     ####
#                                                          #
##%######################################################%##
## Data

cig <- tribble(
  ~n_cig, ~pa, ~n_cas,
  0, 1000,  0,
  5, 1000,  0,
  10, 1000,  2,
  15, 1000,  2,
  20, 1000,  9,
  30, 1000, 10
)

## Régression de Poisson

# Modèle log-linéaire : log(E(ncas_i))=β0+βcig ncig_i
  


fit1 <- glm(n_cas ~ n_cig, family = poisson, data = cig)


# Coefficients

betas <- coef(fit1)

# (Intercept)       n_cig 
# -0.5413352   0.1024930 

tidy(fit1, conf.int = TRUE)

# Statitistiques / adéquation


glance(fit1)

# Prédiction sur l'échelle du prédicteur linéaire


(x <- cbind(1, cig$n_cig))

x %*% betas # nombre de cig * pente de cig. Pour intercept osf donc colonne de 1

predict(fit1)

# Prédiction en nombre de cas

exp(predict(fit1))

predict(fit1, type = "response") # même chose que exp(predict())

augment(fit1) %>% 
  mutate(.pred = exp(.fitted))

## Variance des coefficients de régression

# Var(β^)ˆ=(X^T * W * X)^−1 avec W matrice diagonale des prédictions :

w <- diag(predict(fit1, type = "response")) # mettre en diagonal pour 1 ligne par prediction

(var_beta <- solve(t(x) %*% w %*% x))
(se_beta <- sqrt(diag(var_beta))) # erreur standard

tibble(
  parametre = c("beta_0", "beta_cig"),
  est = betas,
  se = se_beta,
  wald_ci_inf = betas + qnorm(.025) * se,
  wald_ci_sup = betas + qnorm(.975) * se,
  stat = betas / se,
  p = 2 * (pnorm(-abs(stat)))
)

## Variance du prédicteur linéaire

# Var(η^)ˆ= X Var(β^)ˆ XT=X(XTWX)−1XT

(var_eta <- x %*% var_beta %*% t(x))

# Variance des prédictions (prédicteur linéaire + nombre de cas)

cig %>% 
  bind_cols(
    tibble(
      eta = predict(fit1),
      se_eta = sqrt(diag(var_eta)),
      eta_lower = eta + qnorm(.025) * se_eta,
      eta_upper = eta + qnorm(.975) * se_eta,
      n_pred = exp(eta),
      n_pred_lower = exp(eta_lower),
      n_pred_upper = exp(eta_upper),
    )
  )
augment(fit1, se_fit = TRUE)


##%######################################################%##
#                                                          #
####                    Avec offset                     ####
#                                                          #
##%######################################################%##
## Data

cig2 <- readRDS("D:/M2 - SMSDS/M2 SMSDS/Modèles linéaires généralisés, modèles mixtes/Cours/TP-20221111/cig2.rds")


cig2 <- tribble(
  ~n_cig, ~pa, ~n_cas,
   0, 1400,  0,
   5,  900,  0,
  10, 1000,  2,
  15,  850,  2,
  20, 1500,  9,
  30, 1400, 10
)

## Régression de Poisson

# Modèle log-linéaire :

# log(E(ncasi))=β0+βcig ncigi+log(npai)

(fit2 <- glm(n_cas ~ n_cig, offset = log(pa), family = poisson, data = cig2))

# Syntaxe alternative

(fit2 <- glm(n_cas ~ n_cig + offset(log(pa)), family = poisson, data = cig2))

# Coefficients

tidy(fit2, conf.int = TRUE)

Statitistiques / adéquation

glance(fit2)

# Prédiction

augment(fit2) %>% 
  mutate(.pred = exp(.fitted))