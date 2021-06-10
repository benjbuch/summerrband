#' EMSA Band Analysis: Selectivity                                  TEMPLATE A01
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

# define probe sets (if needed)

probe_set_1 <- CPG_COMBINATIONS[1:15]
probe_set_2 <- CPG_COMBINATIONS[c("p01", "q02", "q03", "q04", "q05", "p06")]

# create list with metadata for objects

gel_metadata <- list(
  list(file = "gel_001.txt",
       protein = "protA",
       probe = probe_set_1,
       conc  = 1000),
  list(file = "gel_002.txt",
       protein = "protA",
       probe = probe_set_2,
       conc  = 1000),
  list(file = "gel_003.txt",
       protein = "protA",
       probe = probe_set_2,
       exclude = c("lane_15"),
       conc  = 1000),
  # ... (add above this line; don't forget comma)
  NULL
)

data_emsa <- iqtl_import_all(gel_metadata, path = select_directory())

data_emsa %>%
  filter(band_id == "band_1") %>%
  group_by(protein, conc) %>%
  # mutate(rel_frac = vol_frac / vol_frac[probe == "mC/mC"]) %>%
  ggplot(aes(x = factor(probe, levels = CPG_COMBINATIONS),
             y = vol_frac, fill = as.factor(conc))) +
  geom_col(position = "dodge") +
  labs(x = NULL, y = "fraction bound", fill = "concentration") +
  facet_wrap(vars(protein)) +
  theme_linedraw(base_size = 12) + theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
