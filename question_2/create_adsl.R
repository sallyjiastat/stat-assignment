## ============================================================================
## Program : 02_create_adsl.R
## Purpose : Build ADaM ADSL dataset
## Author  : Yiyang Jia
## ============================================================================

# Load library
library(admiral)
library(pharmaversesdtm)
library(dplyr)
library(lubridate)
library(stringr)

# -- Log file setup ------------------------------------------------------------
dr      <- new.env()
logpath <- "/cloud/project"
logname <- file.path(logpath, "adam_adsl.log")
rout <- file(logname, open = "wt", encoding = "UTF-8")
sink(rout, append = TRUE, split = TRUE)
sink(rout, append = TRUE, type = "message")


#--------------------------------------#
#         Start of Main Program        #
#--------------------------------------#

#Load SDTM datasets
dm <- pharmaversesdtm::dm
vs <- pharmaversesdtm::vs
ex <- pharmaversesdtm::ex
ds <- pharmaversesdtm::ds
ae <- pharmaversesdtm::ae

#ADSL from DM
adsl <- dm %>%
  filter(ARMCD != "Scrnfail" | is.na(ARMCD)) %>%
  select(STUDYID, USUBJID, SUBJID, SITEID, AGE, AGEU, SEX, RACE,
         ETHNIC, ARM, ARMCD, ACTARM, ACTARMCD, COUNTRY,
         RFSTDTC, RFENDTC, RFXSTDTC, RFXENDTC, RFICDTC, RFPENDTC,
         DTHDTC, DTHFL)

#AGEGR9 and AGEGR9N
adsl1 <- adsl %>%
  mutate(
    AGEGR9 = case_when(
      AGE < 18               ~ "<18",
      AGE >= 18 & AGE <= 50  ~ "18 - 50",
      AGE > 50               ~ ">50",
      TRUE                   ~ NA_character_
    ),
    AGEGR9N = case_when(
      AGEGR9 == "<18"        ~ 1,
      AGEGR9 == "18 - 50"    ~ 2,
      AGEGR9 == ">50"        ~ 3,
      TRUE                   ~ NA_real_
    )
  )

#ITTFL
adsl2 <- adsl1 %>%
  mutate(ITTFL = if_else(!is.na(ARM) & ARM != "", "Y", "N"))

#TRTSDTM and TRTSTMF
# TRTSDTM: first exposure per subject
ex_dt <- ex %>%
  filter(
    !is.na(EXSTDTC),
    nchar(EXSTDTC) >= 10,
    EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", toupper(EXTRT)))
  ) %>%
  mutate(
    TRTSDTM = as.POSIXct(paste(EXSTDTC, "00:00:00"), tz = "UTC"),
    TRTSTMF = "H:M:S"
  ) %>%
  group_by(USUBJID) %>%
  slice_min(TRTSDTM, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(USUBJID, TRTSDTM, TRTSTMF)

# TRTEDTM: last exposure per subject
ex_end <- ex %>%
  filter(
    !is.na(EXENDTC),
    nchar(EXENDTC) >= 10,
    EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", toupper(EXTRT)))
  ) %>%
  mutate(
    TRTEDTM = as.POSIXct(paste(EXENDTC, "00:00:00"), tz = "UTC")
  ) %>%
  group_by(USUBJID) %>%
  slice_max(TRTEDTM, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(USUBJID, TRTEDTM)

adsl3 <- adsl2 %>%
  left_join(ex_dt, by = "USUBJID") %>%
  left_join(ex_end, by = "USUBJID")


#LSTAVLDT
# LSTAVLDT using admiral derive_vars_extreme_event()
# Sources: VS, AE, DS, and TRTEDTM from adsl
vs_ext <- vs %>%
  filter(!(is.na(VSSTRESN) & is.na(VSSTRESC))) %>%
  mutate(VSDT = as.Date(substr(VSDTC, 1, 10))) %>%
  filter(!is.na(VSDT)) %>%
  arrange(USUBJID, VSDT, VSSEQ) %>%
  group_by(USUBJID, VSDT) %>%
  slice_max(VSSEQ, n = 1, with_ties = FALSE) %>%  # keep one record per subject/date
  ungroup()

ae_ext <- ae %>%
  mutate(AEDT = as.Date(substr(AESTDTC, 1, 10))) %>%
  filter(!is.na(AEDT)) %>%
  group_by(USUBJID, AEDT) %>%
  slice_max(AESEQ, n = 1, with_ties = FALSE) %>%
  ungroup()

ds_ext <- ds %>%
  mutate(DSDT = as.Date(substr(DSSTDTC, 1, 10))) %>%
  filter(!is.na(DSDT)) %>%
  group_by(USUBJID, DSDT) %>%
  slice_max(DSSEQ, n = 1, with_ties = FALSE) %>%
  ungroup()

adsl3 <- adsl3 %>% mutate(TRTEDT = as.Date(TRTEDTM))

adsl4 <- adsl3 %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "vs_ext",
        order        = exprs(VSDT, VSSEQ),
        condition    = !is.na(VSDT),
        set_values_to = exprs(LSTAVLDT = VSDT)
      ),
      event(
        dataset_name = "ae_ext",
        order        = exprs(AEDT, AESEQ),
        condition    = !is.na(AEDT),
        set_values_to = exprs(LSTAVLDT = AEDT)
      ),
      event(
        dataset_name = "ds_ext",
        order        = exprs(DSDT, DSSEQ),
        condition    = !is.na(DSDT),
        set_values_to = exprs(LSTAVLDT = DSDT)
      ),
      event(
        dataset_name = "adsl3",
        condition    = !is.na(TRTEDT),
        set_values_to = exprs(LSTAVLDT = TRTEDT)
      )
    ),
    source_datasets = list(
      vs_ext = vs_ext,
      ae_ext = ae_ext,
      ds_ext = ds_ext,
      adsl3  = adsl3     # now has TRTEDT
    ),
    tmp_event_nr_var = event_nr,
    order    = exprs(LSTAVLDT, event_nr),
    mode     = "last",
    new_vars = exprs(LSTAVLDT)
  ) %>%
  select(-TRTEDT, -TRTEDTM)

#labels
attr(adsl4$AGEGR9,  "label") <- "Analysis Age Group 9"
attr(adsl4$AGEGR9N, "label") <- "Analysis Age Group 9 (N)"
attr(adsl4$TRTSDTM, "label") <- "Datetime of First Exposure to Treatment"
attr(adsl4$TRTSTMF, "label") <- "Time Imputation Flag for TRTSDTM"
attr(adsl4$ITTFL,   "label") <- "Intent-to-Treat Population Flag"
attr(adsl4$LSTAVLDT,"label") <- "Last Known Alive Date"

#QC ----------------------------------------------------------------------
message(sprintf("ADSL records: %d; Subjects: %d", nrow(adsl4), n_distinct(adsl4$USUBJID)))

#output
saveRDS(adsl4, "/cloud/project/adam_adsl.rds")

#log-Check
ut_rlogcheck <- function(logfile = NULL, outfile = stdout()) {
  # Keywords to search for in the log
  searchfor <- c("error:", "warning:", "Warning message", "not found",
                 "does not exist", "invalid", "missing values were generated")
  # Read log file
  if (is.null(logfile) || !file.exists(logfile)) {
    message("ERROR: log file not found or not specified")
    return(invisible(NULL))
  }
  
  log_lines <- readLines(logfile, warn = FALSE)
  # Search for issues (case-insensitive)
  hits <- log_lines[grepl(paste(searchfor, collapse = "|"), log_lines, ignore.case = TRUE)]
  
  # Write report
  out <- if (is.character(outfile)) file(outfile, "w") else outfile
  
  writeLines(paste("Log file scanned :", logfile),           con = out)
  writeLines(paste("Scanned at      :", Sys.time()),         con = out)
  writeLines(paste("Total lines read :", length(log_lines)), con = out)
  writeLines(paste("Issues found    :", length(hits)),       con = out)
  writeLines(paste(rep("-", 60), collapse = ""),            con = out)
  
  if (length(hits) == 0) {
    writeLines("There're no errors and warnings.", con = out)
  } else {
    writeLines("Lines containing potential issues:", con = out)
    writeLines(hits, con = out)
  }
  
  if (is.character(outfile)) close(out)
  return(invisible(list(hits = hits, n_issues = length(hits))))
}

print(sessionInfo())
##----- Sink output back to console -----
sink()
sink(type = "message")
close(rout)
##----- Check the log for issues -----
dr$logname <- "/cloud/project/adam_adsl.log"
dr$logcheckname <- "/cloud/project/adam_adsl_rlogcheck.lst"
dr$logcheck <- ut_rlogcheck(logfile = dr$logname, outfile = dr$logcheckname)