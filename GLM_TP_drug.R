
pacman::p_load(tictoc, readr, readxl, tidyverse, glmmTMB, lme4, lmerTest, nnet, 
               broom, broom.mixed, leaps, DataEditR, VGAM, pROC, geeM, MuMIn, 
               purrr, forcats, ResourceSelection)


# Data

# Format binomial

drug <- tribble(
  ~drug, ~x, ~r, ~n,
  "A",	.1,	1,	10,
  "A",	.23,	2,	12,
  "A",	.67,	1,	9,
  "B",	.2,	3,	13,
  "B",	.3,	4,	15,
  "B",	.45,	5,	16,
  "B",	.78,	5,	13,
  "C",	.04,	0,	10,
  "C",	.15,	0,	11,
  "C",	.56,	1,	12,
  "C",	.7,	2,	12,
  "D",	.34,	5,	10,
  "D",	.6,	5,	9,
  "D",	.7,	8,	10,
  "E",	.2,	12,	20,
  "E",	.34,	15,	20,
  "E",	.56,	13,	15,
  "E",	.8,	17,	20
)

drug <- readRDS("D:/M2 - SMSDS/M2 SMSDS/Modèles linéaires généralisés, modèles mixtes/Cours/TP-20221111/drug.rds")


# NB : code pour transformer les données en format "binaire" (une ligne par sujet)


drug2 <- drug %>% 
  mutate(id_exp = row_number()) %>% 
  nest(data = c(r, n)) %>% 
  mutate(
    r = map(
      data, 
      ~c(rep(1L, .x$r), rep(0L, .x$n - .x$r))
    )
  ) %>%
  select(id_exp, drug, x, r) %>% 
  unnest(cols = r)


# Régression logistique

Modèle :
  E(y)=E(log(p / 1−p))=β0+βx×x+βB×drugB+βC×drugC+βD×drugD+βE×drugE

  # Recodage en variables indicatrices

drug <- drug %>% 
  mutate(
    B = ifelse(drug == "B", 1, 0),
    C = ifelse(drug == "C", 1, 0),
    D = ifelse(drug == "D", 1, 0),
    E = ifelse(drug == "E", 1, 0)
  )

drug

drug %>% 
  count(drug, B, C, D, E)

## Modèle `fit1` avec covariables `x` et `drug`, variables indicatrices

# Modélisation des proportions de succès (pondération par le nombre d'essais)

(fit1 <- glm(r/n ~ x + B + C + D + E, data = drug, family = binomial, weights = n)) # n c'est nombre de personne et r le nombre de succès
# on estime une proporition influencé par le n de chaque groupe

# Modélisation des nombres de succès / échecs

(fit1 <- glm(cbind(r, n - r) ~ x + B + C + D + E, data = drug, family = binomial))

# Idem avec `drug` en variable catégorielle

(fit1 <- glm(cbind(r, n - r) ~ x + drug, data = drug, family = binomial))

# Alternative : catégorie de référence = `E`

(fit1b <- glm(
  cbind(r, n - r) ~ x + drug, 
  data = drug %>% mutate(drug = fct_relevel(drug, "E")), 
  family = binomial)
 )

# Coefficients

tidy(fit1, conf.int = TRUE)
tidy(fit1b, conf.int = TRUE)

# Statitistiques / adéquation

glance(fit1)
glance(fit1b)

# NB : on peut comparer les coefficients et les statistiques d'adéquation
# entre ce modèle "binomial" et celui utilisant les données au format
# binaire

(fit2 <- glm(r ~ x + drug, data = drug2, family = binomial))

tidy(fit2, conf.int = TRUE)
glance(fit2)

# Log-vraisemblance

(ll_fit1 <- logLik(fit1))

# AIC

-2 * ll_fit1 + 2 * attr(ll_fit1, "df")

# Probabilité de réponse :
#   pi^=1 / ('1+e−yi)
# Avec yi^=β0+βxxi
# Nombre de réponses prédites : ri^=ni∗pi^

drug_pred <- drug %>% 
  mutate(
    p_pred = predict(fit1, type = "response"),
    r_pred = n * p_pred
  )
drug_pred

# NB : NB : ∑iri^=∑iri

sum(drug_pred$r_pred)
sum(drug_pred$r)

# Calcul direct de yi^hat

predict(fit1, type = "link")

# Calcul direct de $\hat{p_i}$

predict(fit1, type = "response")

# Log-vraisemblance

drug_pred <- drug_pred %>% 
  mutate(
    log_lik = r * log(p_pred) + (n - r) * log(1 - p_pred),
    )

# Résidus de Pearson
 # RPi=ri−ri^/ √(nipi^(1−pi^)

drug_pred <- drug_pred %>% 
  mutate(
    resid_pearson = (r - r_pred) / sqrt(r_pred * (1 - p_pred))
  )
drug_pred


# Calcul direct de RPi

resid(fit1, type = "pearson")

#  Déviance résiduelle :
drug_pred <- drug_pred %>% 
  mutate(
    resid_dev = ifelse(
      r == 0,
      sign(r - r_pred) * sqrt(2 * (n - r) * log((n - r) / (n - r_pred))),
      sign(r - r_pred) * sqrt(2 * r * log(r / r_pred) + 2 * (n - r) * log((n - r) / (n - r_pred)))
    )
  )
drug_pred

# Calcul direct de di

resid(fit1, type = "deviance")

# Statistiques par observation

glance(fit1)

## Modèle `fit0` sans covariable

(fit0 <- glm(cbind(r, n - r) ~ 1, data = drug, family = binomial))

# Coefficients

tidy(fit0)

# Statistiques / adéquation

glance(fit0)

# Analyse de la déviance `fit1` / `fit0`

anova(fit1, fit0, test = "Chisq")

# Alternative (test de l'effet marg#inal de chaque variable)

drop1(fit1, test = "Chisq")

# Prédiction

drug %>% 
  bind_cols(
    augment(fit1) %>% 
      select(-(1:2))
  )
