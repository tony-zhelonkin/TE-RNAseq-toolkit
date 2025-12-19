suppressPackageStartupMessages({
  library(data.table)
  library(tibble)
  library(dplyr)
  library(stringr)
})

# Expect TE ids like "HAL1:L1:LINE" (subfamily:family:class)
parse_te_label <- function(x) {
  p <- strsplit(x, ":", fixed = TRUE)[[1]]
  if (length(p) >= 3) list(subfamily = p[1], family = p[2], class = p[3]) else list(subfamily = NA, family = NA, class = NA)
}

# Vectorized parser → data.frame with columns Symbol, Ensembl, subfamily, family, class, type
build_te_annotation <- function(te_ids) {
  parsed <- lapply(te_ids, parse_te_label)
  subf   <- vapply(parsed, `[[`, "", "subfamily")
  fam    <- vapply(parsed, `[[`, "", "family")
  cls    <- vapply(parsed, `[[`, "", "class")

  # We’ll put TE label into Symbol column for convenience; Ensembl is NA
  tibble::tibble(
    Symbol    = te_ids,
    Ensembl   = NA_character_,
    subfamily = subf,
    family    = fam,
    class     = cls,
    type      = "TE"
  )
}