#--------------------------------------#
# ChIP-seq Quality Control with FastQC #
#--------------------------------------#
# 1) List all fastq files
fastq <- list.files("data/ChIP-seq/fastq",
                    pattern=".fastq.gz",
                    full.names=TRUE)

# 2) Create directory for output FastQC report
dir.create("data/ChIP-seq/FastQC/")

# 3) Run FastQC
system(paste("fastqc",
             paste0(fastq, collapse=" "), # Print all filenames
             "-o data/ChIP-seq/FastQC", # Output folder for FastQC reports
             "-t 5")) # Use 5 threads
