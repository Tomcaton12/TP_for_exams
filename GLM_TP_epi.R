## ---------------------------
## R version : 4.2.2
## 
## Script name: 
##
## Purpose of script: Inspiration
##
## Author: X
##
## Date Created: 2022-11-16
##
## Copyright (c) X, 2022
## Email: 
##
## ---------------------------
##
## Notes:
##    Regression de Poisson (nombre total de crise Ysum)
##    Regression de Poisson AVEC offset = log(time) [nombre de crise par période]
##     
## ---------------------------
##
## load up the packages we will need:
pacman::p_load(tictoc, readr, readxl, tidyverse, glmmTMB, lme4, lmerTest, nnet, 
               broom, broom.mixed, leaps, DataEditR, VGAM, pROC, geeM, MuMIn, 
               purrr, forcats, ResourceSelection)

## ---------------------------
## load dataset
#  It is better to use the "long format" for storing data and use the
#  "wide format" at the very end of a data analysis process to reduce 
#  the data dimensionality. 

# Data
## Description variables
# Format "wide"
(epi <- robustbase::epilepsy)

# Format "long" pour `Y1:Y4` c a d 4 colonne distinctes pour le nombre de crise
# sur 4 périodes différentes d'étude
epi_long <- epi %>% 
  rename(Y0 = Base) %>% # change le nom de la colonne Base en Y0
  pivot_longer(cols = c(Y1:Y0), names_to = "period", values_to = "n", names_prefix = "Y") %>% 
  # transforme les colonnes Y1 à Y0 (4 colonne) en 1 : Ysum
  # 1assigner les valeurs du nombre de colonne prise a la colonne n,
  #  et mettre Y devant pour qu'il apparaissent pas en valeur
  select(ID, Trt, period, n) %>% # colonne a garder
  mutate(
    visite1_4 = ifelse(period != 0, 1, 0), # 1 visite dans les périodes Y1 a 4
    t = ifelse(period == 0, 8, 2), # si je rencontre un 0 je met 8 sinon un 2
  ) %>% 
  arrange(ID, period) # arrange le tableau par ID et period dans l'ordre croissant

epi_long

# Format "long" pour "Ysum"
epi_long2 <- epi %>% 
  pivot_longer(cols = c(Base, Ysum), names_to = "period", values_to = "n") %>% 
  select(ID, Trt, period, n)

epi_long2

## Description du nombre de crises
epi %>% 
  group_by(Trt) %>% # étudier les traitements sur le summarize séparement
  summarise(
    mean = mean(Ysum),
    sd = sd(Ysum),
    Q1 = quantile(Ysum, .25),
    Q3 = quantile(Ysum, .75),
    min = min(Ysum),
    max = max(Ysum),
    sum = sum(Ysum) # Mon rajout
  )

epi %>% 
  # filter(ID != 207) %>% # filter de l'outlier il y en a dautres, tous sauf le 207
  ggplot(aes(x = Trt, y = Ysum, col = Trt)) +
  geom_boxplot() +
  theme_bw()

# Outlier
epi %>% filter(Ysum == 302) # filtre de bcp de crise
epi_long2 %>% 
  filter(ID != 207) %>% # filtre de l'individu 207
  filter(period %in% c("Base", "Ysum")) %>% # %in% pour identifie si c'est un 
  # vecteur ou df et selectionner une colonne dans une df, créer une nouvelle 
  # variable en colonne, ou retirer une colonne
  ggplot(aes(x = period, y = n, col = Trt)) +
  geom_boxplot() +
  theme_bw()

## Corrélation entre `Base`et `Ysum`
cor.test(
  epi[epi$ID != 207, "Base"], 
  epi[epi$ID != 207, "Ysum"]
)

# Dans le groupe `progabide`, on veut le moins de corrélation possible entre base et ysum
# pour montrer l'efficacité du traitement
cor.test(
  epi[epi$Trt == "progabide" & epi$ID != 207, "Base"], 
  epi[epi$Trt == "progabide" & epi$ID != 207, "Ysum"]
)

# Dans le groupe `placebo`
cor.test(
  epi[epi$Trt == "placebo", "Base"], 
  epi[epi$Trt == "placebo", "Ysum"]
)

epi %>% 
  filter(ID != 207) %>% 
  ggplot(aes(Base, Ysum, col = Trt)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_bw()

# Régression de Poisson (nombre total de crises `Ysum`)

# Sans interaction

fit_glm1 <- glm(n ~ period + Trt, data = epi_long2 %>% filter(ID != 207), family = poisson)

fit_glm1 %>% tidy(conf.int = TRUE)
# Interprétation :
# Intercept : le nombre de crise moyen lorsque la periode est
# 0 est estimée à 3.49
# PeriodYsum : lorsque l'unité d'Ysum augmente d'une unité
# le nombre de crise diminue de -0.02576
# TrtProgabide : le traitement progabide, réduit le nombre de crise
# au cours du temps de - 0.25582

# Avec interaction
fit_glm2 <- glm(
  n ~ period * Trt, 
  data = epi_long2 %>% filter(ID != 207), 
  family = poisson
)
fit_glm2 %>% tidy(conf.int = TRUE)
# Interprétation:
# On estime que sur la périodeBase, le nombre moyen de crise est de 3.43 tout groupe
# confondu, lorsque la PeriodeYsum est de 0. Lorsqu'on augmente une unité de periodeYsum
# le nombre de crise augmente de 0.111, tout groupe confondu.
# Le coefficient de Trtprogabide estime la différence d'intercept entre
# la periodeBase et TrtProgabide, cette différence est de - 0.108 crises.
# Le coefficient PeriodYsum:Trtprogabide estime la différence de pente entre
# le groupe Progabide pendant la Base au groupe Progabide pendant la période Ysum.
# On estime donc le TrtProgabide, produit de 3.43 + (- 0.108) = 3.322 nombre
# de crise pour la periode Ysum a 0, et lorsqu'on augmente la période Ysum d'une
# unité, le nombre de crise pour le groupe TrtProgabide est le produit de
# 0.111 - 0.302 = - 0.191 (donc diminue de - 0 .191 crise)


# Régression de Poisson avec offset = log(time) (nombre de crises par période)
# On ajoute un offset pour rapporter le nombre de crises au temps de mesure.
# Par rapport à `fit_glm2`, l'intercept est donc modifié mais pas 
# les coefficients qui correspondent au log des risques relatifs pour les prédicteurs.

fit_glm3 <- glm(
  n ~ visite1_4 * Trt + offset(log(t)), 
  data = epi_long %>% filter(ID != 207), 
  family = poisson
)
fit_glm3 %>% tidy(conf.int = TRUE)

# Modèle log-linéaire :

# log(E(Yi,j))= β0 + β1xi,1+ β2xi,2 + β3xi,1 xi,2 + log(ti,j)

# Et :

# Yi,j : nombre de crises épileptiques sur l’intervalle j
# ti,j : durée de l’intervalle j
# xi,1 = 0 pour les semaines 0 à 8 (baseline) et 1 pour les semaines suivantes (traitement)
# xi,2 = 0 pour le groupe placebo et 1 pour le groupe progabide

# Interprétation des paramètres de régression
# 
# | Traitement | Visite   | log(E(Yi,j) / ti,j)
# :-----------:|:--------:|:------------------------------------------------:|
# | Placebo    | Baseline | β0                          |
# |            | 1 à 4    | β0+β1                           |
# | Progabide  | Baseline | β0+β2                           |
# |            | 1 à 4    | β0+β1+β2+β3         |


# NB : les coefficients de régression snt les mêmes avec et sans offset, seul l'intercept change.
# Dans le premier cas, il représente le log du nombre de crises sur 8 semaines (E(Yi,j)) et 
# dans le second cas sur 1 semaine  (E(Yi,j) / ti,j).
# On retrouve ce ratio (ou cette différence, sur une échelle log) entre les intercepts : exp(3.4270508 - 1.3476092) = 8
  
  
  
  ## Régression GEE
  # Résolution ci-dessous avec le package `geeM`. Possible également avec les packages `gee`, `yags`, `geepack` ou `wgee` notamment.

### Modèle de Poisson avec offset = log(time), corrélation échangeable
fit_gee_e <- geem(
  n ~ visite1_4 * Trt + offset(log(t)), 
  id = ID, 
  data = epi_long %>% filter(ID != 207) %>% as.data.frame(), 
  family = poisson, 
  corstr = "exch" # ECHANGEABLE
)

summary(fit_gee_e)

QIC(fit_gee_e)

# Matrice de travail
fit_gee_e$biggest.R.alpha


### Modèle de Poisson avec offset = log(time), corrélation non structurée
### GEE estime pour une population, pas l'échantillon
### https://rlbarter.github.io/Practical-Statistics/2017/05/10/generalized-estimating-equations-gee/

# Selon Hardin & Hilbe (2003), il faut analyser le critère QIC et choisir le modèle
# avec le plus petit QIC. Faut-il que ce soit en valeur absolue ?

fit_gee_u <- geem(
  n ~ visite1_4 * Trt + offset(log(t)), 
  id = ID, 
  data = epi_long %>% filter(ID != 207) %>% as.data.frame(), 
  family = poisson, 
  corstr = "unstr" # NON STRUCTUREE
)

summary(fit_gee_u)

QIC(fit_gee_u)

# Matrice de travail
fit_gee_u$biggest.R.alpha # diagonale ?

### Modèle de Poisson avec offset =log(time), corrélation autorégressive
fit_gee_a <- geem(
  n ~ visite1_4 * Trt + offset(log(t)), 
  id = ID, 
  data = epi_long %>% filter(ID != 207) %>% as.data.frame(), 
  family = poisson, 
  corstr = "ar1"
)

summary(fit_gee_a)

QIC(fit_gee_a)

# Matrice de travail
fit_gee_a$biggest.R.alpha

# NB : $0.149 = 0.621 * 0.239 = 0.621^2 * 0.386 = 0.621^3 * 0.621 = 0.621^4$
  
  ## Modèle de Poisson avec effet aléatoire (random intercept)
fit_glmm <- glmer(
  n ~ visite1_4 * Trt + (1 | ID) + offset(log(t)),
  data = epi_long %>% filter(ID != 207),
  family = poisson
)
summary(fit_glmm)

## Comparaison des modèles GEE et GLMM
epi_pred <- epi_long %>% filter(ID != 207) %>% 
  mutate(
    pred_glmm = predict(fit_glmm, type = "response"),
    pred_gee_e = fitted(fit_gee_e),
    pred_gee_u = fitted(fit_gee_u),
    pred_gee_a = fitted(fit_gee_a)
  )
# NB : la fonction predict ne fonctionne pas 
epi_pred

epi_pred %>% 
  pivot_longer(cols = starts_with("pred"), names_to = "model", values_to = "pred") %>% 
  mutate(model = str_remove(model, "pred_")) %>% 
  ggplot(aes(n, pred, color = model)) +
  geom_point() +
  theme_bw() +
  labs(x = "nombre de crises observé", y = "nombre de crises prédit") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

# Les prédictions des modèles GEE sont assez proches quelle que soit la matrice 
# de corrélation choisie. Pour un modèle donné, les prédictions 
# (rapportées au temps d'exposition) sont les mêmes pour une même combinaison de
# prédicteurs indépendamment du sujet concerné. Les modèles random intercept (GLMM) 
# donnent presque toujours de meilleures prédictions car elles prennent en
# compte les caractéristiques individuelles du sujet.

# NB : dans ces modèles, on n'a pas pris en compte l'évolution temporelle du 
# nombre de crises (variable `period`), ce qui aurait pu avoir un
# intérêt et explique que les prédictions soient les mêmes à tous les temps 
# d'observation (même combinaison de prédicteurs).
