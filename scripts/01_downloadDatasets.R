#----------------------------------------------------#
# Download datasets needed for the ChIP-seq workshop #
#----------------------------------------------------#
# GEO ID: GSE112221
# DOI: 10.1002/hep.30211
# !! Careful, this step takes a long time !!
#----------------------------------------------------#
library(doParallel) # Allow parallel for loops
registerDoParallel(cores=2)

# Load samples info
samples <- read.csv("data/GSE112221_samples.csv",
                    stringsAsFactors = FALSE)

# Downloas SRA in FASTQ format
SRAToFastq <- function(SRA) {
  message(paste(">> Convert to fastq", SRA))
  line <- paste("fastq-dump",
                "--split-3", # send each pair to a different file
                "--gzip", # compress fastq output
                "-O data/ChIP-seq/fastq/",
                SRA)
  system(line)
}

foreach(i=samples$sra_id) %dopar% SRAToFastq(i)

# Rename files to match sample names
path <- "data/ChIP-seq/fastq/"

file.rename(paste0(path, samples$sra_id, "_1.fastq.gz"),
            paste0(path, samples$sample_names, ".fastq.gz"))

# Remove other files
files <- list.files(path,
                    pattern="SRR",
                    full.names=TRUE)
file.remove(files)
