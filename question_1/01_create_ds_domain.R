## ===============================================================
## Program : 01_create_sdtm_ds.R
## Purpose : Build SDTM DS domain from pharmaverseraw::ds_raw
## Author  : Yiyang Jia
## ===============================================================

# Load library
library(sdtm.oak)
library(pharmaverseraw)
library(pharmaversesdtm)
library(dplyr)
library(tibble)
library(stringr)

# Load the study_ct data
study_ct <- data.frame(
  stringsAsFactors = FALSE,
  codelist_code = rep("C66727", 10),
  term_code = c("C41331","C25250","C28554","C48226","C48227",
                "C48250","C142185","C49628","C49632","C49634"),
  term_value = c("ADVERSE EVENT","COMPLETED","DEATH","LACK OF EFFICACY",
                 "LOST TO FOLLOW-UP","PHYSICIAN DECISION","PROTOCOL VIOLATION",
                 "SCREEN FAILURE","STUDY TERMINATED BY SPONSOR",
                 "WITHDRAWAL BY SUBJECT"),
  collected_value = c("Adverse Event","Complete","Dead","Lack of Efficacy",
                      "Lost To Follow-Up","Physician Decision",
                      "Protocol Violation","Trial Screen Failure",
                      "Study Terminated By Sponsor",
                      "Withdrawal by Subject"),
  term_preferred_term = c("AE","Completed","Died",NA,NA,NA,
                          "Violation", "Failure to Meet Inclusion/Exclusion Criteria",
                          NA,"Dropout"),
  term_synonyms = c("ADVERSE EVENT","COMPLETE","Death",NA,NA,NA,NA,NA,NA,
                    "Discontinued Participation")
)

# -- Log file setup ------------------------------------------------------
dr      <- new.env()
logpath <- "/cloud/project"
logname <- file.path(logpath, "sdtm_ds.log")
rout    <- file(logname, open = "wt", encoding = "UTF-8")
sink(rout, append = TRUE, split = TRUE)
sink(rout, append = TRUE, type = "message")

#--------------------------------------#
#           Start of Main Program      #
#--------------------------------------#

#Load raw data and add trace keys
# generate_oak_id_vars() adds oak_id, raw_source, patient_number columns
# which are required join keys for all sdtm.oak assign_xxx functions
ds_raw <- pharmaverseraw::ds_raw %>%
  generate_oak_id_vars( pat_var = "PATNUM", raw_src = "ds_raw")

#Map topic variable first
ds <- assign_no_ct(
  raw_dat = ds_raw,
  raw_var = "IT.DSTERM",
  tgt_var = "DSTERM",
  id_vars = oak_id_vars()
)

#Map qualifiers onto ds
# DSDECOD: CT lookup against C66727
# Build a unified lookup: any of the 3 columns -> put into 'term_value'
ct_c66727 <- study_ct %>%
  filter(codelist_code == "C66727")

ct_unified <- bind_rows(
  ct_c66727 %>% select(match_val = collected_value,  term_value),
  ct_c66727 %>% select(match_val = term_synonyms,  term_value),
  ct_c66727 %>% select(match_val = term_preferred_term,  term_value),
  data.frame(match_val = "Randomized",  term_value = "RANDOMIZED"),
  data.frame(match_val = "Screen Failure", term_value = "SCREEN FAILURE")
) %>%
  filter(!is.na(match_val)) %>%
  distinct(match_val, .keep_all = TRUE)

# Direct lookup: raw IT.DSDECOD to DSDECOD using all 3 columns as keys
ds_raw <- ds_raw %>%
  mutate(
    DSDECOD_mapped = ct_unified$term_value[
      match(tolower(IT.DSDECOD),
            tolower(ct_unified$match_val))
    ]
  )

# Assign using the pre-mapped variable
ds <- assign_no_ct(
  tgt_dat = ds,
  raw_dat = ds_raw,
  raw_var = "DSDECOD_mapped",
  tgt_var = "DSDECOD",
  id_vars = oak_id_vars()
)

# DSDTCOL and IT.DSSTDAT
# convert to YYYY-MM-DD
ds_raw <- ds_raw %>%
  mutate(
    IT.DSSTDAT_iso = format(as.Date(IT.DSSTDAT, "%m-%d-%Y"), "%Y-%m-%d")
  )

# DSDTC
# Extract time from "Final Lab Visit" rows
time_lookup <- ds_raw %>%
  filter(!is.na(DSTMCOL) & is.na(IT.DSTERM)) %>%
  select(patient_number, DSDTCOL, DSTMCOL) %>%
  mutate(
    DSDTCOL_iso = format(as.Date(DSDTCOL, "%m-%d-%Y"), "%Y-%m-%d")
  )

# Join time onto rows that have DSTERM and same subject/date
ds_raw <- ds_raw %>%
  mutate(
    DSDTCOL_iso = format(as.Date(DSDTCOL, "%m-%d-%Y"), "%Y-%m-%d")
  ) %>%
  left_join(
    time_lookup %>%
      select(patient_number, DSDTCOL_iso, DSTMCOL_from_lab = DSTMCOL),
    by = c("patient_number", "DSDTCOL_iso")
  ) %>%
  mutate(
    DSTMCOL_final = case_when(
      !is.na(DSTMCOL) ~ DSTMCOL,
      !is.na(DSTMCOL_from_lab) ~ DSTMCOL_from_lab,
      TRUE ~ NA_character_
    )
  )

# Then build DSDTC using DSTMCOL_final instead of DSTMCOL
ds1 <- ds %>%
  left_join(
    ds_raw %>%
      select(oak_id, patient_number, DSDTCOL_iso, DSTMCOL_final), by = c("oak_id", "patient_number")
  ) %>%
  mutate(
    DSDTC = case_when(
      !is.na(DSDTCOL_iso) & !is.na(DSTMCOL_final)  ~ paste0(DSDTCOL_iso, "T", DSTMCOL_final),
      !is.na(DSDTCOL_iso)       ~ DSDTCOL_iso,
      TRUE          ~ NA_character_
    )
  ) %>%
  select(-DSDTCOL_iso, -DSTMCOL_final)

# DSSTDTC: date only
ds1 <- assign_datetime(
  tgt_dat = ds1,
  raw_dat = ds_raw,
  raw_var = c("IT.DSSTDAT_iso"),
  tgt_var = "DSSTDTC",
  raw_fmt = c("y-m-d"),
  id_vars = oak_id_vars()
)

# VISIT: direct copy from INSTANCE
ds2 <- assign_no_ct(
  tgt_dat = ds1,
  raw_dat = ds_raw,
  raw_var = "INSTANCE",
  tgt_var = "VISIT",
  id_vars = oak_id_vars()
)

#Merge STUDY/PATNUM back for USUBJID construction
ds3 <- ds2 %>%
  left_join(
    ds_raw %>% select(oak_id, patient_number, STUDY, PATNUM), by = c("oak_id", "patient_number")
  ) %>%
  mutate(
    STUDYID = STUDY, DOMAIN  = "DS", USUBJID = paste(STUDY, PATNUM, sep = "-")
  )

#DSCAT: per aCRF, randomisation records get "PROTOCOL MILESTONE"
ds4 <- ds3 %>%
  mutate(
    DSCAT = case_when(
      toupper(DSDECOD) == "RANDOMIZED"         ~ "PROTOCOL MILESTONE",
      !is.na(DSTERM)   ~ "DISPOSITION EVENT",
      TRUE   ~ NA_character_
    )
  )

# -- VISITNUM: numeric visit mapping from aCRF schedule ------------------
visit_map <- tribble(
  ~VISIT,               ~VISITNUM,
  "Screening 1",        1,
  "Baseline",           2,
  "Week 2",             3,
  "Week 4",             4,
  "Week 6",             5,
  "Week 8",             6,
  "Week 12",            7,
  "Week 16",            8,
  "Week 20",            9,
  "Week 24",           10,
  "Week 26",           11,
  "Retrieval",         12,
  "Ambul Ecg Removal", 13,
  "Unscheduled 1.1",  201,
  "Unscheduled 4.1",  204,
  "Unscheduled 5.1",  205,
  "Unscheduled 6.1",  206,
  "Unscheduled 8.2",  208,
  "Unscheduled 13.1", 213
)

ds5 <- ds4 %>% left_join(visit_map, by = "VISIT")

#DSSTDY: study day relative to RFSTDTC from DM
#day 1 = reference date, no day 0
dm_dates <- pharmaversesdtm::dm %>%
  select(USUBJID, RFSTDTC) %>%
  filter(!is.na(RFSTDTC)) %>%
  mutate( PATNUM = sub("^01-", "", USUBJID)
  )

ds6 <- ds5 %>%
  left_join(
    dm_dates %>% select(PATNUM, RFSTDTC), by = "PATNUM"
  ) %>%
  mutate(
    DSSTDY = case_when(
      is.na(DSSTDTC) | is.na(RFSTDTC)     ~ NA_integer_,
      as.Date(as.character(DSSTDTC)) >= as.Date(RFSTDTC) 
              ~ as.integer(as.Date(as.character(DSSTDTC)) - as.Date(RFSTDTC) + 1L),
      TRUE   ~ as.integer(as.Date(as.character(DSSTDTC)) - as.Date(RFSTDTC)
      )
    )
  ) %>%
  select(-RFSTDTC)

#DSSEQ and final variable selection
ds_final <- ds6 %>%
  filter(!is.na(DSTERM)) %>%
  arrange(USUBJID, DSSTDTC,     DSTERM ) %>%
  group_by(USUBJID) %>%
  mutate(DSSEQ = row_number()) %>%
  ungroup() %>%
  select(STUDYID,DOMAIN,USUBJID,DSSEQ,DSTERM,DSDECOD,DSCAT,VISITNUM,VISIT,DSDTC,DSSTDTC,DSSTDY
  )

#labels
attr(ds_final$STUDYID,  "label") <- "Study Identifier"
attr(ds_final$DOMAIN,   "label") <- "Domain Abbreviation"
attr(ds_final$USUBJID,  "label") <- "Unique Subject Identifier"
attr(ds_final$DSSEQ,    "label") <- "Sequence Number"
attr(ds_final$DSTERM,   "label") <- "Reported Term for the Disposition Event"
attr(ds_final$DSDECOD,  "label") <- "Standardized Disposition Term"
attr(ds_final$DSCAT,    "label") <- "Category for Disposition Event"
attr(ds_final$VISITNUM, "label") <- "Visit Number"
attr(ds_final$VISIT,    "label") <- "Visit Name"
attr(ds_final$DSDTC,    "label") <- "Date/Time of Collection"
attr(ds_final$DSSTDTC,  "label") <- "Start Date/Time of Disposition Event"
attr(ds_final$DSSTDY,   "label") <-"Study Day of Start of Disposition Event"


#QC and output ------------------------------------------------------------
message(sprintf("DS records: %d ; Subjects: %d",nrow(ds_final), n_distinct(ds_final$USUBJID)))
print(count(ds_final,DSDECOD,sort = TRUE))
colSums(is.na(ds_final))

# Output
saveRDS(ds_final,  "/cloud/project/sdtm_ds.rds")

#log-Check
ut_rlogcheck <- function( logfile = NULL, outfile = stdout()) {

  # Keywords to search for in the log
  searchfor <- c( "error:","warning:","Warning message","not found", "does not exist","invalid", "missing values were generated")

  # Read log file
  if (is.null(logfile) || !file.exists(logfile)) {
    message("ERROR: log file not found or not specified")
    return(invisible(NULL))
  }

  log_lines <- readLines(logfile,warn = FALSE  )

  # Search for issues (case-insensitive)
  hits <- log_lines[
    grepl(
      paste(searchfor, collapse = "|"),
      log_lines,
      ignore.case = TRUE
    )
  ]

  # Write report
  out <- if (is.character(outfile))  file(outfile, "w")   else  outfile
  writeLines(paste("Log file scanned :", logfile), con = out)
  writeLines(paste("Scanned at      :", Sys.time()), con = out)
  writeLines(paste("Total lines read :", length(log_lines)),con = out )
  writeLines(paste("Issues found    :", length(hits)),con = out )
  writeLines(paste(rep("-", 60), collapse = ""),con = out)

  if (length(hits) == 0) {
    writeLines("There're no errors and warnings.",con = out)
  } else {
    writeLines("Lines containing potential issues:",con = out)
    writeLines(hits,con = out)
  }
  if (is.character(outfile))
    close(out)
  return(invisible(list(hits = hits,n_issues = length(hits))))
}


#----- Display package and version information -----
print(sessionInfo())
##----- Sink output back to console -----
sink()
##----- Sink errors/warnings back to console -----
sink(type = "message")
##----- Close the log file -----
close(rout)

#----- Check the log for issues -----
dr$logname <-
  "/cloud/project/sdtm_ds.log"
dr$logcheckname <-
  "/cloud/project/sdtm_ds_rlogcheck.lst"
dr$logcheck <- ut_rlogcheck(logfile = dr$logname,outfile = dr$logcheckname)
