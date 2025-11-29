## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6
)

## ----setup--------------------------------------------------------------------
library(crmPack)

## ----data-ordinal-1-----------------------------------------------------------
empty_ordinal_data <- DataOrdinal(
  doseGrid = c(seq(from = 10, to = 100, by = 10)),
  yCategories = c("No tox" = 0L, "Sub-tox AE" = 1L, "DLT" = 2L),
  placebo = FALSE
)

## ----data-ordinal-2, fig.alt = "A graph showing Patient ID on the x axis and dose administered on the y axis.  The shape and colour of the symbols indicate the toxicity status of the patient: red triangles for DLTs, orange circles for sub-toxic AEs and black triangles for no reported toxicities.  Patients 1 to 4 are dosed at 10, 20, 30 and 40, with no toxicitis reported.  patients 5 to 7 are dosed at 50, with patient 6 reporting a sub-toxic AE.  Patients 8 to 10 are treated at 60.  Patient 9 reports a sub-toxic AE and patient 10 a DLT."----
dose_grid_range(empty_ordinal_data)

ordinal_data <- update(empty_ordinal_data, x = 10, y = 0)
ordinal_data <- update(ordinal_data, x = 20, y = 0)
ordinal_data <- update(ordinal_data, x = 30, y = 0)
ordinal_data <- update(ordinal_data, x = 40, y = 0)
ordinal_data <- update(ordinal_data, x = 50, y = c(0, 1, 0))
ordinal_data <- update(ordinal_data, x = 60, y = c(0, 1, 2))

plot(ordinal_data)

## ----logisticlogordinal-------------------------------------------------------
ordinal_model <- LogisticLogNormalOrdinal(
  mean = c(3, 4, 0),
  cov = diag(c(4, 3, 1)),
  ref_dose = 55
)

## -----------------------------------------------------------------------------
opts <- .DefaultMcmcOptions()

samples <- mcmc(ordinal_data, ordinal_model, opts)

## ----samples-slot-names-------------------------------------------------------
names(samples@data)

## ----fit-1--------------------------------------------------------------------
fit(samples, ordinal_model, ordinal_data, grade = 1L)
fit(samples, ordinal_model, ordinal_data, grade = 2L)

## ----fit-2--------------------------------------------------------------------
fit(samples, ordinal_model, ordinal_data, grade = 1L, cumulative = FALSE)
fit(samples, ordinal_model, ordinal_data, grade = 2L, cumulative = FALSE)

## ----plot1, fig.alt = "A graph of the posterior probability of toxicity (DLT only) against dose.  The mean probability of toxicity is barely above 0% at a dose of zero and rises in a sigmoidal curve to around 65% at a dose of 100.  The confidence interval is relatively narrow for low doses but widens considerably for doses over 60, extending from around 15% to 100% for a dose of 100."----
plot(samples, ordinal_model, ordinal_data, grade = 2L)

## ----plot2, fig.alt = "A graph of the posterior cumulative probability of toxicity (sub-toxic AE or DLT) against dose.  The mean probability of toxicity is barely above 0% at a dose of zero and rises in a sigmoidal curve to around 75% at a dose of 100.  The confidence interval is relatively narrow for low doses but widens considerably for doses over 60, extending from around 30% to 100% for a dose of 100."----
plot(samples, ordinal_model, ordinal_data, grade = 1L)

## ----plot3, fig.alt = "A graph of the posterior probability of sub toxic AE against dose.  The mean probability of toxicity is barely above 0% at a dose of zero, rises to a peak of about 18% at a dose of 60 before falling to around 12% at a dose of 100.  The confidence interval is relatively narrow for low doses but widens considerably for doses over 60, extending from around 30% to 100% for a dose of 100."----
plot(samples, ordinal_model, ordinal_data, grade = 1L, cumulative = FALSE)

## ----rules-1------------------------------------------------------------------
dlt_rule <- CohortSizeDLT(intervals = 0:2, cohort_size = c(1, 3, 5))
ordinal_rule_1 <- CohortSizeOrdinal(grade = 1L, rule = dlt_rule)
ordinal_rule_2 <- CohortSizeOrdinal(grade = 2L, rule = dlt_rule)

size(ordinal_rule_1, 50, empty_ordinal_data)
size(ordinal_rule_2, 50, empty_ordinal_data)
size(ordinal_rule_1, 50, ordinal_data)
size(ordinal_rule_2, 50, ordinal_data)

## ----rules-2------------------------------------------------------------------
ordinal_rule_1 <- IncrementsOrdinal(
  grade = 1L,
  rule = IncrementsRelativeDLT(intervals = 0:2, increments = c(3, 1.5, 0.67))
)
maxDose(ordinal_rule_1, ordinal_data)
ordinal_rule_2 <- IncrementsOrdinal(
  grade = 2L,
  rule = IncrementsRelativeDLT(intervals = 0:1, increments = c(3, 0.5))
)
maxDose(ordinal_rule_2, ordinal_data)

## ----rules-3------------------------------------------------------------------
trial_rule <- IncrementsMin(list(ordinal_rule_1, ordinal_rule_2))
maxDose(trial_rule, ordinal_data)

## ----logisticlognormal--------------------------------------------------------
model <- LogisticLogNormal(
  mean = c(-3, 1),
  cov = matrix(c(4, -0.5, -0.5, 3), ncol = 2),
  ref_dose = 45
)

model@params@cov

## ----logisticlognormal-samples------------------------------------------------
data <- Data(doseGrid = seq(10, 100, 10))
options <- McmcOptions(
  samples = 30000,
  rng_kind = "Mersenne-Twister",
  rng_seed = 8191316
)
samples <- mcmc(data, model, options)

## ----logisticlognormal-covariance---------------------------------------------
d <- as.matrix(cbind(samples@data$alpha0, log(samples@data$alpha1)))
sigmaHat <- cov(d)
sigmaHat

## ----ordinal-with-covariance-2------------------------------------------------
ordinal_model_temp <- ordinal_model
ordinal_model_temp@params@cov <- matrix(c(4, -0.5, -0.5, -0.5, 3, -0.5, -0.5, -0.5, 1), ncol = 3)

ordinal_model_temp@params@cov

## ----ordinal-with-covariance-3------------------------------------------------
ordinal_data <- DataOrdinal(doseGrid = seq(10, 100, 10))
ordinal_samples <- mcmc(ordinal_data, ordinal_model_temp, options)

## ----ordinal-with-covariance-4------------------------------------------------
ordinalD <- as.matrix(
  cbind(
    ordinal_samples@data$alpha1,
    ordinal_samples@data$alpha2,
    log(ordinal_samples@data$beta)
  )
)
sigmaHat <- cov(ordinalD)
sigmaHat

## ----environment, echo = FALSE------------------------------------------------
sessionInfo()

