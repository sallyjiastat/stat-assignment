## ============================================================================
## Program : 03_create_tf102.R
## Purpose : TFL - Plots
## Author  : Yiyang Jia
## ============================================================================

# Continue on Q3 - 01 summary table
library(ggplot2)

# -- Log setup ----------------------------------------------------------------
dr      <- new.env()
logpath <- "/cloud/project"
rout    <- file(file.path(logpath, "tf102.log"), open = "wt", encoding = "UTF-8")
sink(rout, append = TRUE, split = TRUE)
sink(rout, append = TRUE, type = "message")

#-------------------------------#
#         Start of Main Program #
#-------------------------------#

# --- Plot 1: AE Severity Distribution by Treatment ---------------------------
sev_data <- teae %>%
  filter(!is.na(AESEV)) %>%
  count(ACTARM, AESEV) %>%
  mutate(AESEV = factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))

p1 <- ggplot(sev_data, aes(x = ACTARM, y = n, fill = AESEV)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      "MILD" = "#F8766D",
      "MODERATE" = "#00BA38",
      "SEVERE" = "#619CFF"
    ),
    name = "Severity/Intensity"
  ) +
  labs(
    title = "AE severity distribution by treatment",
    x = "Treatment Arm",
    y = "Count of AEs"
  ) +
  theme_gray() +
  theme(legend.position = "right")

ggsave(
  file.path(logpath, "ae_severity_bar.png"),
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

# --- Plot 2: Top 10 Most Frequent AEs with 95% Clopper-Pearson CI ----------
n_total <- n_distinct(teae$USUBJID)

top10 <- teae %>%
  group_by(AETERM) %>%
  summarise(
    n_subj = n_distinct(USUBJID),
    .groups = "drop"
  ) %>%
  arrange(desc(n_subj)) %>%
  slice_head(n = 10) %>%
  mutate(
    pct = n_subj / n_total,
    # Clopper-Pearson exact CI
    ci_lo = mapply(function(x, n) binom.test(x, n)$conf.int[1],
                   n_subj, n_total),
    ci_hi = mapply(function(x, n) binom.test(x, n)$conf.int[2],
                   n_subj, n_total),
    AETERM = factor(AETERM, levels = rev(AETERM))
  )

p2 <- ggplot(top10, aes(x = pct, y = AETERM)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Top 10 Most Frequent Adverse Events",
    subtitle = paste0("n = ", n_total, " subjects; 95% Clopper-Pearson CIs"),
    x = "Percentage of Patients (%)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(logpath, "ae_top10_ci.png"),
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

#log-Check
ut_rlogcheck <- function(logfile = NULL, outfile = stdout()) {
  
  # Keywords to search for in the log
  searchfor <- c(
    "error:",
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
  
  return(invisible(list(
    hits = hits,
    n_issues = length(hits)
  )))
}

print(sessionInfo())

sink()
sink(type = "message")
close(rout)

dr$logname      <- file.path(logpath, "tf102.log")
dr$logcheckname <- file.path(logpath, "tf102_rlogcheck.lst")
dr$logcheck     <- ut_rlogcheck(
  logfile = dr$logname, outfile = dr$logcheckname)

