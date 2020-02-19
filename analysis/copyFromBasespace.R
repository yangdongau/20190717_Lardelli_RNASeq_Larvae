library(magrittr)
library(tidyverse)

# Where the files are
from <- "/data/basespace/Projects/Lardelli_10_July_2019/Samples"
to <- "/data/20190717_Lardelli_RNASeq_Larvae/0_rawData/fastq"
dirs <- list.files(from)
wSpace <- grepl(" ", dirs) #Whitespace in the path. Grrr!

stopifnot(file.exists(from))
stopifnot(file.exists(to))
stopifnot(length(dirs) > 0)

# Copy where there's no whitespace in the path
file.path(from, dirs[!wSpace], "Files") %>%
  sapply(list.files, pattern = "fastq.gz", full.names = TRUE, simplify = FALSE) %>%
  lapply(function(x){
    d <- unique(dirname(x)) 
    fl <- basename(x)
    smp <- unique(basename(dirname(d)))
    smp <- str_replace(smp, " ", "_") 
    smp <- str_replace_all(smp, "[\\(\\)]", "")
    smp <- paste0(smp, "_R1.fq.gz")
    message("Copying\n", x, "\nto\n", file.path(to, smp))
    system2("cat", paste(paste(x, collapse = " "), ">", file.path(to, smp)))
  })

# Copy where there's whitespace in the path
file.path(from, dirs[wSpace], "Files") %>%
  sapply(list.files, full.name = TRUE, simplify = FALSE) %>%
  lapply(function(x){
    outPath <- file.path(to, basename(x))
    message("Copying\n", paste(x, collapse = "\n"), "\nto\n", paste(outPath, collapse = "\n"))
    file.copy(x, outPath, overwrite = TRUE)
    smp <- unique(basename(dirname(dirname(x))))
    smp <- str_replace(smp, " ", "_")
    smp <- str_replace_all(smp, "[\\(\\)]", "")
    smp <- paste0(smp, "_R1.fq.gz")
    cat("Merging\n", paste(outPath, collapse = "\n"), "\ninto\n", file.path(to, smp), "\n")
    system2("cat",
            paste(paste(outPath, collapse = " "), ">", file.path(to, smp)))
    rmFiles <- list.files(to, pattern = "00[1-4]", full.names = TRUE)
    message("Deleting\n", paste(rmFiles, collapse = "\n"))  
    file.remove(rmFiles)
  })
