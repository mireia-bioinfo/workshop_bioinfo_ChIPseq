#----------------------------------#
# ChIP-seq Peak Calling with MACS2 #
#----------------------------------#
# 1) Select BAM files
bam <- list.files("data/ChIP-seq/BAM",
                  pattern="*rmdup.bam$",
                  full.names=TRUE)

# 2) Divide in inputs and samples
input_files <- bam[grepl("input", bam)]
sample_files <- bam[!grepl("input", bam)]

# 3) Run peak calling
dir.create("data/ChIP-seq/peaks/", F)
dir.create("tmp/", F)

peakCalling <- function(sample, input) {
  message(paste(">> Calling peaks for", sample))
  line <- paste("macs2 callpeak",
                "-f BAM -t", sample, # Control sample
                "-c", input, # Input sample
                "-g hs", # Genome size (hs = human)
                "-n", gsub("BAM/", "peaks/", gsub(".rmdup.bam", "", sample)), # Output name
                "--tempdir tmp/", # Custom tmp directory
                # "--bdg", # Output bedgraph
                "--broad --broad-cutoff 0.1 --nomodel") # Generate broadpeaks
  system(line)
}

v <- mapply(peakCalling, sample=sample_files, input=input_files)

system("rm -r tmp/")