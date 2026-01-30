## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(crmPack)

# Define the dose-grid.
emptydata <- Data(
    doseGrid = c(0.1, 0.2, 0.5, 1, 3, 5, 10, 15, 20, 25, 40, 50, 60, 70, 80, 100)
)

# Define the dose-toxicity model.
model <- LogisticLogNormal(
    mean = c(-0.85, 1),
    cov = matrix(c(5, -0.5, -0.5, 5), nrow = 2),
    ref_dose = 56
)

# Choose the rule for selecting the next dose.
myNextBest <- NextBestNCRM(
    target = c(0.2, 0.35),
    overdose = c(0.35, 1),
    max_overdose_prob = 0.25
)

# Choose the rule for stopping.
myStopping1 <- StoppingMinCohorts(nCohorts = 3)
myStopping2 <- StoppingTargetProb(
    target = c(0.2, 0.35),
    prob = 0.5
)
myStopping3 <- StoppingMinPatients(nPatients = 40)
myStopping4 <- StoppingPatientsNearDose(nPatients = 10L, percentage = 30, include_backfill = FALSE)
myStopping <- (myStopping1 & myStopping2 & myStopping4) |
    myStopping3 |
    StoppingMissingDose()

# Choose the rule for dose increments.
myIncrements <- IncrementsRelative(
    intervals = c(0, 20, 50),
    increments = c(1, 0.67, 0.33)
)

## -----------------------------------------------------------------------------
design_no_backfill <- Design(
    model = model,
    nextBest = myNextBest,
    stopping = myStopping,
    increments = myIncrements,
    cohort_size = CohortSizeConst(3),
    data = emptydata,
    startingDose = 3
)
design_no_backfill@backfill

## -----------------------------------------------------------------------------
backfill_simple <- Backfill(
    cohort_size = CohortSizeConst(3),
    max_size = 12,
    opening = OpeningMinCohorts(min_cohorts = 1),
    recruitment = RecruitmentUnlimited(),
    priority = "lowest"
)
backfill_simple

## -----------------------------------------------------------------------------
design_simple_backfill <- design_no_backfill
design_simple_backfill@backfill <- backfill_simple

## -----------------------------------------------------------------------------
backfill_complex <- Backfill(
    cohort_size = CohortSizeRandom(min_size = 1, max_size = 6),
    opening = OpeningMinCohorts(min_cohorts = 3) &
      OpeningMinResponses(
        min_responses = 1, 
        include_lower_doses = TRUE
    ),
    recruitment = RecruitmentRatio(ratio = 1 / 2),
    priority = "highest",
    max_size = 20
)
backfill_complex

## -----------------------------------------------------------------------------
design_complex_backfill <- design_no_backfill
design_complex_backfill@backfill <- backfill_complex

## -----------------------------------------------------------------------------
# Assumed dose-response probability function.
mytruthResponse <- function(dose) {
    plogis(- 4 + 0.2 * dose) / 4
}
curve(mytruthResponse(x), from = 0, to = max(emptydata@doseGrid), 
      xlab = "Dose", ylab = "Probability of Response / Toxicity", 
      main = "Assumed Functions", ylim = c(0, 1))

myTruth <- probFunction(design_simple_backfill@model, alpha0 = 3, alpha1 = 3)
curve(myTruth(x), from = 0, to = max(emptydata@doseGrid), 
      add = TRUE, col = "red")

## -----------------------------------------------------------------------------
# For real applications, use e.g. McmcOptions() with defaults.
mcmcOptions <- McmcOptions(
    burnin = 10, 
    step = 1, 
    samples = 100, 
    rng_kind = "Mersenne-Twister", 
    rng_seed = 12345
)

# Simple backfill design simulation:
sims_simple <- simulate(
    design_simple_backfill,
    truth = myTruth,
    nsim = 10, # For real applications, increase to 1000 e.g.
    seed = 819,
    mcmcOptions = mcmcOptions,
    parallel = FALSE,
    firstSeparate = FALSE
)

# Complex backfill design simulation:
sims_complex <- simulate(
    design_complex_backfill,
    truth = myTruth,
    truthResponse = mytruthResponse,
    nsim = 10, # For real applications, increase to 1000 e.g
    seed = 819,
    mcmcOptions = mcmcOptions,
    parallel = FALSE,
    firstSeparate = FALSE
)

## -----------------------------------------------------------------------------
plot(sims_simple@data[[3]], mark_backfill = TRUE)

## -----------------------------------------------------------------------------
plot(sims_complex@data[[5]], mark_backfill = TRUE, mark_response = TRUE)

## -----------------------------------------------------------------------------
get_backfill_counts <- function(sims) {
    sapply(sims@data, \(d) sum(d@backfilled))
}
backfill_counts_simple <- get_backfill_counts(sims_simple)
backfill_counts_complex <- get_backfill_counts(sims_complex)
table(backfill_counts_simple)
table(backfill_counts_complex)

## -----------------------------------------------------------------------------
get_backfill_doses <- function(sims) {
    lapply(sims@data, \(d) d@x[d@backfilled])
}
backfill_doses_simple <- get_backfill_doses(sims_simple)
backfill_doses_complex <- get_backfill_doses(sims_complex)

## -----------------------------------------------------------------------------
head(backfill_doses_simple, 3)

## -----------------------------------------------------------------------------
all_backfill_doses_simple <- unlist(backfill_doses_simple)
table(all_backfill_doses_simple)

all_backfill_doses_complex <- unlist(backfill_doses_complex)
table(all_backfill_doses_complex)

## -----------------------------------------------------------------------------
summary(sims_simple, truth = myTruth)

