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
##   Regression linéaire
##   Regression logistique
##   Regression polytomique / multinomial
##   
## ---------------------------
##
## load up the packages we will need:
pacman::p_load(tictoc, readr, readxl, tidyverse, glmmTMB, lme4, lmerTest, nnet, 
               broom, broom.mixed, leaps, DataEditR, VGAM)

## ---------------------------
## load dataset
dat <- read_delim("D:/M2 - SMSDS/M2 SMSDS/Modèles linéaires généralisés, modèles mixtes/Cours/TP-20221111/dataset1_smsds_import_raw.csv", ";", escape_double = FALSE, trim_ws = TRUE)


# Régression linéaire
# "elisa" = critère continu

## Échelle originelle
dat %>% 
  ggplot(aes(elisa)) +
  geom_histogram() +
  theme_bw()

(fit01 <- lm(elisa ~ I(age / 10), data = dat))

# Interprétation : lorsque l'âge augmente de 10 ans, le titre ELISA augmente en 
# moyenne de `r fit01$coefficients[[2]]` (ou diminue en moyenne de 
# fit01$coefficients[[2]]

# Échelle log10
dat %>% 
  ggplot(aes(log10(elisa))) +
  geom_histogram() +
  theme_bw()

(fit02 <- lm(log10(elisa) ~ I(age / 10), data = dat))
10 ^ fit02$coefficients[[2]]
# Interprétation :
# 
# - lorsque l'âge augmente de 10 ans, le log10 du titre ELISA augmente en 
# moyenne de fit02$coefficients[[2]]`(ou diminue en moyenne 
# de fit02$coefficients[[2]]
# - lorsque l'âge augmente de 10 ans, le titre ELISA varie en moyenne 
# d'un facteur 10 ^ fit02$coefficients[[2]] (ou diminue en moyenne de 
# (1 - 10 ^ fit02$coefficients[[2]]) * 100`%)


## Échelle log2
dat %>% 
  ggplot(aes(log2(elisa))) +
  geom_histogram() +
  theme_bw()

(fit03 <- lm(log2(elisa) ~ I(age / 10), data = dat))
2 ^ fit03$coefficients[[2]]

# Interprétation :
#   
# - lorsque l'âge augmente de 10 ans, le log2 du titre ELISA augmente en 
# moyenne de fit03$coefficients[[2]] (ou diminue en moyenne 
# de fit03$coefficients[[2]])
# - lorsque l'âge augmente de 10 ans, le titre ELISA varie en moyenne 
# d'un facteur 2 ^ fit03$coefficients[[2]] (même interprétation que pour 
# toute transformation logarithmique)


## Échelle log (naturel = népérien)
dat %>% 
  ggplot(aes(log(elisa))) +
  geom_histogram() +
  theme_bw()

(fit04 <- lm(log(elisa) ~ I(age / 10), data = dat))
exp(fit04$coefficients[[2]])
```
# Interprétation :
# 
# - lorsque l'âge augmente de 10 ans, le log du titre ELISA augmente en 
# moyenne de fit04$coefficients[[2]] (ou diminue en moyenne 
# de fit04$coefficients[[2]]
# - lorsque l'âge augmente de 10 ans, le titre ELISA varie en moyenne d'un 
# facteur exp(fit04$coefficients[[2]]) 
# (même interprétation que pour toute transformation logarithmique)
# - Remarque : la valeur du coefficient non transformé est très proche de 
# la variation relative (-7%). Cette approximation, issue du développement 
# limité de log(1 + x) lorsque x tend vers 0, est acceptable lorsque 
# le coefficient est "petit" et permet de visualiser directement 
# l'effet d'une variable sans transformation exponentielle.

# Régression logistique
# elisa2 = recodage en 2 classes

dat <- dat %>% 
  mutate(
    elisa2 = case_when(
      elisa >= 1.1 ~ 1,
      elisa < 1.1 ~ 0,
    ),
    age65 = case_when(
      age >= 65 ~ 1,
      age < 65 ~ 0,
    ),
  )

xtabs(~age65 + elisa2, data = dat)

(fit11 <- glm(elisa2 ~ age65, data = dat, family = binomial))
plogis(fit11$coefficients[[1]])
plogis(sum(fit11$coefficients))
exp(fit11$coefficients[[2]])

tidy(fit11, conf.int = TRUE)
tidy(fit11, conf.int = TRUE, exponentiate = TRUE)

# Interprétation :
#   
# - le risque absolu (probabilité de test ELISA positif) 
# chez les sujets < 65 ans est estimé à plogis(fit11$coefficients[[1]])
# - le risque absolu (probabilité de test ELISA positif) 
# chez les sujets ≥ 65 ans est estimé à plogis(sum(fit11$coefficients))
# l'OR associé à un âge ≥ 65an est estimé à exp(fit11$coefficients[[2]])`

# Modèles multivariés

(fit12 <- glm(elisa2 ~ age65 + covid, data = dat, family = binomial))

# Test de l'interaction

(fit13 <- glm(elisa2 ~ age65 + covid + age65:covid, data = dat, family = binomial))
(fit13b <- glm(elisa2 ~ age65 * covid, data = dat, family = binomial))

## Comparaison de modèles emboîtés
fit12$coefficients
fit13$coefficients
anova(fit12, fit13, test = "LRT")

# Interprétation : d'après le test du rapport de vraisemblance, le terme 
# d'interaction n'améliore pas significativement l'adéquation du modèle.

# Ce test apporte une information similaire au test du coefficient 
# dans la description du modèle :

tidy(fit13)

# Il serait d'un intérêt particulier si on avait une covariable 
# à $k > 2$ modalités : on aurait alors $k - 1$ tests de coefficients 
# dans la description du modèle mais un seul résultat de test en comparant 
# les modèles emboîtés résumant l'intérêt de la variable pour
#  l'ensemble de ses modalités.

(fit14 <- glm(elisa2 ~ covid, data = dat %>% mutate(age20 = cut(age, seq(20, 80, 20))), family = binomial))
(fit15 <- glm(elisa2 ~ age20 + covid, data = dat %>% mutate(age20 = cut(age, seq(20, 80, 20))), family = binomial))

tidy(fit15)
anova(fit14, fit15, test = "Chisq")

# La fonction `drop1` permet de tester simultanément l'effet marginal de toutes
# les covariables d'un modèle (modèle complet comparé aux modèles emboîtés 
# en retirant chacune de ces variables, une à la fois).

(fit16 <- glm(elisa2 ~ age65 + covid + sexe, data = dat, family = binomial))
drop1(fit16, test = "Chisq")

## Comparaison de modèles non emboîtés
# Le test du rapport de vraisemblance est ici impossible.

fit17 <- glm(formula = elisa2 ~ sexe + covid, family = binomial, data = dat)
fit12$coefficients
fit17$coefficients

anova(fit12, fit17, test = "Chisq")

# On peut alors avoir recours à des indices d'adéquation comme l'AIC ou la BIC :
glance(fit12)
glance(fit17)

## Écriture binomiale

dat_bin <- dat %>% 
    filter(!is.na(elisa2)) %>% 
  count(age65, covid, elisa2) %>% 
    pivot_wider(names_from = elisa2, values_from = n)

dat_outcomes <- tibble(
  pos = dat_bin[, 4],
  neg = dat_bin[, 3]
) %>% 
  as.matrix()

dat_covar <- tibble(
  age65 = dat_bin[, 1],
  covid = dat_bin[, 2]
) %>% 
  as.matrix()


(fit18 <- glm(dat_outcomes ~ dat_covar, family = binomial))


# Régression polytomique
# elisa3 = recodage en 3 classes

dat <- dat %>% 
  mutate(
    elisa3 = case_when(
      elisa >= 1.1 ~ "pos",
      elisa >= .8 ~ "dou",
      elisa < .8 ~ "neg"
    ) %>% 
      fct_relevel("neg") %>% # on choisit "neg" comme catégorie de référence
      ordered() # il s'agit d'une variable ordinale (l'ordre des niveaux a un sens)
)

# Package `nnet`
summary(multinom(elisa3 ~ age65 + covid, data = dat))

# Chaque coefficient estime l'effet par rapport à la catégorie `neg`.
# Les résultats sont proches de ceux qu'on aurait obtenus sur des modèles 
# comparant les classes 2 à 2 mais on a ici une seule vraisemblance pour 
# l'ensemble du modèle.

tidy(glm(elisa3 ~ age65 + covid, data = dat %>% filter(elisa3 != "pos"), family = binomial))
tidy(glm(elisa3 ~ age65 + covid, data = dat %>% filter(elisa3 != "dou"), family = binomial))

# Package `VGAM`

# Le choix de la famille permet d'indiquer comment les associations entre les 
# différents niveaux de la variable réponse doivent être définies.  
# Ici, `acat()` estime les log-OR définis entre des niveaux successifs 
# (1 = `dou` vs. `neg`, 2 = `pos` vs. `dou`).

vglm(elisa3 ~ age65 + covid, data = dat, family = acat()) %>% summary

# Le log-OR entre `pos` et `neg` peut en être déduit et correspond aux 
# estimations obtenues précédemment.
# -4.3924 + 0.8012 # (Intercept)
# 0.2429 + -1.0271 # age65
# 0.4201 + 0.9434 # covid

# D'autres définitions sont possibles, comme l'effet cumulatif 
# comparant 2 à 2 les valeurs agrégées des niveaux inférieurs ou supérieurs 
# à un seuil donné.

vglm(elisa3 ~ age65 + covid, data = dat, family = cumulative(reverse = TRUE)) %>% summary

# On peut également imposer que l'effet des différentes covariables soit 
# identique pour les différentes comparaisons de niveaux (dans ce cas, seul 
# l'intercept dépend des niveaux comparés et reflète les écarts relatifs 
# entre la fréquence des différents niveaux).

vglm(elisa3 ~ age65 + covid, data = dat, family = cumulative(parallel = TRUE, reverse = TRUE)) %>% summary