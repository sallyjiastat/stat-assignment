## ============================================================================
## Program : 03_create_tf101.R
## Purpose : TFL - Adverse Events Summary Table using
## Author  : Yiyang Jia
## ============================================================================

# Load library
library(gt)
library(gtsummary)
library(dplyr)
library(stringr)
library(pharmaverseadam)

# -- Log setup ----------------------------------------------------------------
dr      <- new.env()
logpath <- "/cloud/project"
rout    <- file(file.path(logpath, "tf101.log"), open = "wt", encoding = "UTF-8")
sink(rout, append = TRUE, split = TRUE)
sink(rout, append = TRUE, type = "message")

#-------------------------------#
#         Start of Main Program #
#-------------------------------#

# Load datasets
adae <- pharmaverseadam::adae
adsl <- pharmaverseadam::adsl

# Treatment group N counts (denominator)
trt_n <- adsl %>%
  filter(SAFFL == "Y") %>%
  count(ACTARM, name = "N")

# Filter TEAEs
teae <- adae %>%
  filter(TRTEMFL == "Y", SAFFL == "Y")

# Subject counts per SOC and PT
ae_counts <- teae %>%
  group_by(ACTARM, AESOC, AETERM) %>%
  summarise(n_subj = n_distinct(USUBJID), .groups = "drop")

# Subject counts at SOC level
soc_counts <- teae %>%
  group_by(ACTARM, AESOC) %>%
  summarise(n_subj = n_distinct(USUBJID), .groups = "drop") %>%
  mutate(AETERM = AESOC, level = "SOC")

# Subject counts overall TEAE row
overall_counts <- teae %>%
  group_by(ACTARM) %>%
  summarise(n_subj = n_distinct(USUBJID), .groups = "drop") %>%
  mutate(
    AESOC  = "Treatment Emergent AEs",
    AETERM = "Treatment Emergent AEs",
    level  = "TOTAL"
  )

# PT level
pt_counts <- ae_counts %>%
  mutate(level = "PT")

# Combine
all_counts <- bind_rows(overall_counts, soc_counts, pt_counts)

# Merge N
all_counts <- all_counts %>%
  left_join(trt_n, by = "ACTARM") %>%
  mutate(
    pct  = round(n_subj / N * 100, 1),
    cell = paste0(n_subj, " (", pct, "%)")
  )

# wide format: one column per treatment arm
wide <- all_counts %>%
  select(AESOC, AETERM, level, ACTARM, cell) %>%
  tidyr::pivot_wider(
    names_from  = ACTARM,
    values_from = cell,
    values_fill = "0 (0%)"
  )

# Sort:
soc_order <- wide %>%
  filter(level == "SOC") %>%
  arrange(AESOC) %>%
  pull(AESOC)

pt_sorted <- wide %>%
  filter(level == "PT") %>%
  left_join(
    all_counts %>%
      filter(level == "PT", grepl("High", ACTARM)) %>%
      select(AETERM, n_subj),
    by = "AETERM"
  ) %>%
  arrange(AESOC, desc(n_subj))

total_row <- wide %>% filter(level == "TOTAL")

soc_pt_rows <- bind_rows(
  lapply(soc_order, function(s) {
    bind_rows(
      wide %>% filter(level == "SOC", AESOC == s),
      pt_sorted %>% filter(AESOC == s)
    )
  })
)

final_tbl <- bind_rows(total_row, soc_pt_rows)

# label
final_tbl <- final_tbl %>%
  mutate(
    label = case_when(
      level == "TOTAL" ~ AETERM,
      level == "SOC"   ~ AESOC,
      level == "PT"    ~ paste0("  ", AETERM)
    )
  )

trt_cols <- trt_n$ACTARM

col_labels <- setNames(
  paste0(
    trt_cols,
    "\nN = ",
    trt_n$N[match(trt_cols, trt_n$ACTARM)]
  ),
  trt_cols
)

# Create gt table
tbl <- final_tbl %>%
  select(label, all_of(trt_cols)) %>%
  gt() %>%
  cols_label(
    label = gt::html(
      "<b>Primary System Organ Class</b><br><b>Reported Term for the Adverse Event</b>"
    ),
    .list = setNames(
      lapply(trt_cols, function(x)
        gt::html(
          paste0(
            "<b>", x, "</b><br>N = ",
            trt_n$N[trt_n$ACTARM == x],
            "<sup>t</sup>"
          )
        )
      ),
      trt_cols
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = final_tbl$level == "SOC"
    )
  ) %>%
  tab_style(
    style = cell_text(indent = px(20)),
    locations = cells_body(
      columns = label,
      rows = final_tbl$level == "PT"
    )
  ) %>%
  tab_header(
    title = "Adverse Event Summary by System Organ Class and Preferred Term",
    subtitle = "Safety Population"
  ) %>%
  tab_footnote(
    footnote = "N = number of subjects in treatment group",
    locations = cells_column_labels(
      columns = all_of(trt_cols[1])
    )
  ) %>%
  opt_align_table_header(align = "left") %>%
  cols_align(align = "center", columns = all_of(trt_cols)) %>%
  cols_align(align = "left", columns = label)

# QC
message(sprintf(
  "TEAE subjects: %d | SOCs: %d | PTs: %d",
  n_distinct(teae$USUBJID),
  n_distinct(teae$AESOC),
  n_distinct(teae$AETERM)
))

# Save outputs
gtsave(tbl, file.path(logpath, "ae_summary_table.html"))

#log-Check
ut_rlogcheck <- function(logfile = NULL, outfile = stdout()) {
  
  # Keywords to search for in the log
  searchfor <- c(
    "error:",
    "warning:",
    "Warning message",
    "not found",
    "does not exist",
    "invalid",
    "missing values were generated"
  )
  
  # Read log file
  if (is.null(logfile) || !file.exists(logfile)) {
    message("ERROR: log file not found or not specified")
    return(invisible(NULL))
  }
  
  log_lines <- readLines(logfile, warn = FALSE)
  
  # Search for issues (case-insensitive)
  hits <- log_lines[
    grepl(
      paste(searchfor, collapse = "|"),
      log_lines,
      ignore.case = TRUE
    )
  ]
  
  # Write report
  out <- if (is.character(outfile)) file(outfile, "w") else outfile
  
  writeLines(paste("Log file scanned :", logfile), con = out)
  writeLines(paste("Scanned at      :", Sys.time()), con = out)
  writeLines(paste("Total lines read :", length(log_lines)), con = out)
  writeLines(paste("Issues found     :", length(hits)), con = out)
  writeLines(paste(rep("-", 60), collapse = ""), con = out)
  
  if (length(hits) == 0) {
    writeLines("There're no errors and warnings.", con = out)
  } else {
    writeLines("Lines containing potential issues:", con = out)
    writeLines(hits, con = out)
  }
  
  if (is.character(outfile)) close(out)
  
  return(
    invisible(
      list(
        hits = hits,
        n_issues = length(hits)
      )
    )
  )
}

print(sessionInfo())

sink()
sink(type = "message")
close(rout)

dr$logname      <- file.path(logpath, "tf101.log")
dr$logcheckname <- file.path(logpath, "tf101_rlogcheck.lst")
dr$logcheck     <- ut_rlogcheck(logfile = dr$logname,outfile = dr$logcheckname
)