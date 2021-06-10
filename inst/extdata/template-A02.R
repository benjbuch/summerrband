#' EMSA Band Analysis: Affinity                                     TEMPLATE A02
#'
#' analysis by:
#'          on: <<TODAY>>
#'
#' template by: Benjamin Buchmuller
#'          on: 210528
#'
#' <<RVERSION>>
#'
#' -----------------------------------------------------------------------------
#' BASIC USAGE
#' Modify from template as need be.
#' -----------------------------------------------------------------------------

library(summerrband)
library(tidyverse)

# create list with metadata for objects

gel_metadata <- list(
  list(file = "gel_001.txt",
       protein = "protA",
       probe = "p01",
       conc  = c(2^seq(10, 0), 0)),  ## use `seq` to generate series
  list(file = "gel_002.txt",
       protein = "protA",
       probe = rep(c("p03", "p04"), each = 6),  ## use `rep` to repeat elements
       exclude = c("lane_15"),
       conc  = rep(c(0, 2, 5, 10, 20, 50), 2)),
  # ... (add above this line; don't forget comma)
  NULL
)

data_emsa <- iqtl_import_all(gel_metadata, path = select_directory())

data_fits <- data_emsa %>%
  filter(band_id == "band_1") %>%
  group_by(protein, probe) %>%
  model_cleanly_groupwise(fit_Kd, formula = vol_frac ~ conc, R0 = 2,
                          newdata = data.frame(conc = 10^seq(0, 3, length.out = 100)))

data_fits %>%
  select(tidy) %>% unnest(tidy) %>%
  filter(term == "K_d")

data_fits %>%
  model_display(color = CPG_COMBINATIONS[probe]) +
  scale_x_log10() +
  # cale_color_manual(values = c("C/C" = "black", "C/mC" = "red", "X/Y" = "blue")) +
  labs(x = "concentration", y = "fraction bound", color = "probe") +
  facet_wrap(vars(protein)) +
  theme_linedraw(base_size = 12) + theme(
    legend.position = "bottom"
  )
