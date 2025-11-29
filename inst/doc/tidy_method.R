## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages({
  library(crmPack)
  library(knitr)
  library(kableExtra)
  library(tidyr)
  library(magrittr)
  library(dplyr)
})

## ----error = TRUE-------------------------------------------------------------
try({
CohortSizeConst(size = 3) %>% tidy()
})

## -----------------------------------------------------------------------------
IncrementsRelative(
  intervals = c(0, 20),
  increments = c(1, 0.33)
) %>%
  tidy()

## -----------------------------------------------------------------------------
cs_max <- maxSize(
  CohortSizeConst(3),
  CohortSizeDLT(intervals = 0:1, cohort_size = c(1, 3))
)
cs_max %>% tidy()

## -----------------------------------------------------------------------------
options <- McmcOptions(
  burnin = 100,
  step = 1,
  samples = 2000
)

emptydata <- Data(doseGrid = c(1, 3, 5, 10, 15, 20, 25, 40, 50, 80, 100))

model <- LogisticLogNormal(
  mean = c(-0.85, 1),
  cov =
    matrix(c(1, -0.5, -0.5, 1),
      nrow = 2
    ),
  ref_dose = 56
)
samples <- mcmc(emptydata, model, options)
tidySamples <- samples %>% tidy()
tidySamples %>% head()

## -----------------------------------------------------------------------------
CohortSizeRange(
  intervals = c(0, 50, 300),
  cohort_size = c(1, 3, 5)
) %>%
  tidy() %>%
  kable(
    col.names = c("Min", "Max", "Cohort size"),
    caption = "Rules for selecting the cohort size"
  ) %>%
  add_header_above(c("Dose" = 2, " " = 1))

## ----fig.width = 6, fig.height = 4--------------------------------------------
options <- McmcOptions(
  burnin = 5000,
  step = 1,
  samples = 40000
)

data <- Data(
  doseGrid = c(1, 3, 5, 10, 15, 20, 25, 40, 50, 80, 100),
  x = c(1, 3, 5, 10, 15, 15, 15),
  y = c(0, 0, 0, 0, 0, 1, 0),
  ID = 1L:7L,
  cohort = as.integer(c(1:4, 5, 5, 5))
)

model <- LogisticLogNormal(
  mean = c(-1, 0),
  cov =
    matrix(c(3, -0.1, -0.1, 4),
      nrow = 2
    ),
  ref_dose = 56
)
samples <- mcmc(data, model, options)
tidySamples <- samples %>% tidy()

# The magrittr pipe is necessary here
tidySamples$data %>%
  expand(
    nesting(!!!.[1:10]),
    Dose = data@doseGrid[2:11]
  ) %>%
  mutate(Prob = probFunction(model, alpha0 = alpha0, alpha1 = alpha1)(Dose)) %>%
  ggplot() +
  geom_density(aes(x = Prob, colour = as.factor(Dose)), adjust = 1.5) +
  labs(
    title = "Posterior dose-specific PDFs for p(Tox)",
    caption = "Dose 1 omitted as p(Tox) is essentially 0",
    x = "p(Tox)"
  ) +
  scale_colour_discrete("Dose") +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )

## -----------------------------------------------------------------------------
sessionInfo()

