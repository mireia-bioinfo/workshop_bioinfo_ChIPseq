#----------------------------------------#
# ChIP-seq Post-Processing with Samtools #
#----------------------------------------#
# 1) List aligned files
files <- list.files("data/ChIP-seq/BAM",
                    pattern="*.sam$",
                    full.names=TRUE)

# 2) Convert to BAM and sort
toBAMandSort <- function(sam) {
  message(paste(">> Convert to BAM and sort", sam))
  line <- paste("samtools view",
                sam,
                "-@ 6",
                "-b", # Convert to BAM
                "|", # Pipe
                "samtools sort",
                "-", # Input from pipe
                "-m 2G -@ 6",
                "-o", gsub(".sam", ".bam", sam))
  system(line)
}

sapply(files, toBAMandSort)

# 3) Remove duplicated reads
# Select raw bam files
files <- list.files("data/ChIP-seq/BAM",
                    pattern="*.raw.bam$",
                    full.names=TRUE)

rmDups <- function(bam) {
  message(paste(">> Remove duplicates", bam))
  line <- paste("samtools markdup",
                bam, # Input
                gsub(".raw.bam", ".rmdup.bam", bam), # Output
                "-r", # Remove duplicates
                "-s", # Report stats
                "-@ 6")
  system(line)
}

v <- sapply(files, rmDups)

# 4) Index raw and rmdup BAM
# Select all BAM files
files <- list.files("data/ChIP-seq/BAM",
                    pattern="*.bam",
                    full.names=TRUE)

indexBAM <- function(bam) {
  line <- paste("samtools index",
                bam,
                "-@ 6")
  system(line)
}

v <- sapply(files, indexBAM)