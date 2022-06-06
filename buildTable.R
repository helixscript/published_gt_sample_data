library(dplyr)
library(stringr)
library(RMySQL)
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))

intSiteDataPath <- '/media/lorax/data/export/projects/exchange'

# Requires specimen_management db entry in ~/.my.cnf file.
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp')

# Parse data files, retrieve anything that looks like a GTSP id.
r <- bind_rows(lapply(list.files('data', full.names = FALSE), function(x){
  o <- unlist(str_split(x, '_'))
  s <- tibble(GTSP = unique(unlist(str_extract_all(readLines(file.path('data', x)), 'GTSP\\d+'))))
  s <- left_join(s, select(samples, SpecimenAccNum, Trial, Patient, Timepoint, CellType), by = c('GTSP' = 'SpecimenAccNum'))
  s$SRA_ID   <- o[1]
  s$PubMedID <- o[2]
  s$Vector   <- o[3]
  s
}))


expandTimePoints <- function (tps) {
  d <- tibble::tibble(tp = sub("_", ".", tps))
  d$n <- 1:nrow(d)
  d$timePointType <- stringr::str_match(base::toupper(d$tp),"[DMY]")
  d$timePointType[which(is.na(d$timePointType))] <- "X"
  
  d <- dplyr::bind_rows(lapply(split(d, d$timePointType), function(x) {
    n <- as.numeric(stringr::str_match(x$tp, "[\\d\\.]+")) * ifelse(grepl("\\-", x$tp), -1, 1)
    
    if (x$timePointType[1] == "D") {
      x$timePointMonths <- base::round(n/30.4167, digits = 0)
      x$timePointDays <- base::round(n, digits = 0)
    }
    else if (x$timePointType[1] == "M") {
      x$timePointMonths <- base::round(n, digits = 0)
      x$timePointDays <- base::round(n * 30.4167, digits = 0)
    }
    else if (x$timePointType[1] == "Y") {
      x$timePointMonths <- base::round(n * 12, digits = 0)
      x$timePointDays <- base::round(n * 365, digits = 0)
    }
    else {
      message("Warning - could not determine date unit for: ", 
              paste0(unique(x$timePoint), collapse = ", "))
      x$timePointMonths <- n
      x$timePointDays <- n
    }
    x
  }))
  
  dplyr::arrange(d, n) %>% dplyr::select(timePointMonths, timePointDays)
}

# Expand time  points to days and months.
r <- bind_cols(r, expandTimePoints(r$Timepoint))

# Manual vector overrides.
r[grepl('SCID1', r$Trial),]$Vector <- 'gamma'


# Write output files.
openxlsx::write.xlsx(r, file = 'published_gt_samples.xlsx')
readr::write_tsv(r, 'published_gt_samples.tsv')

# Retrieve all site data, filter against published ids, write out.
sites <- bind_rows(lapply(system(paste('find', intSiteDataPath, '-name intSites.tsv'), intern = TRUE), function(x){
                          readr::read_tsv(x, col_types = cols(timePoint = "c", internalSampleID = "c", externalSampleID = "c"))
         })) %>% filter(internalSampleID %in% r$GTSP)

write_tsv(sites, 'published_gt_sites.tsv')
system('gzip published_gt_sites.tsv')

# Vector view.
o <- group_by(r, Trial, Vector) %>% summarise(publishedSamples = n()) %>% ungroup()
openxlsx::write.xlsx(o, file = 'trial_vectors.xlsx')