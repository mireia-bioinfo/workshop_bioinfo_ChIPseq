#----------------------------------#
# ChIP-seq Alignment with Bowtie 2 #
#----------------------------------#

#---------------------------------------------------
# Pre-step: Download reference genome and
# create Bowtie2 index.
# !! You will only have to do this the first time
#---------------------------------------------------
# 1) Download hg19 genome
system("rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit data/reference_genome/.")

# 2) Uncompress with twoBitToFa (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
system("twoBitToFa data/reference_genome/hg19.2bit data/reference_genome/hg19.fa")

# 3) Generate index with Bowtie2
system("bowtie2-build data/reference_genome/hg19.fa data/reference_genome/hg19 --threads 5")

#---------------------------------------------------
# Align files with Bowtie2
#---------------------------------------------------
# 1) List fastq files
fastq <- list.files("data/ChIP-seq/fastq",
                    pattern=".fastq.gz",
                    full.names=TRUE)

# 2) Define alignment function
align <- function(fastq_file) {
  message(paste(">> Alignment for", fastq_file))
  line <- paste("bowtie2",
                "-t", # Print time used to load indexes
                "-x data/reference_genome/hg19", # Index to use (only preffix)
                "-U", fastq_file, # Single-end file
                "-S", gsub(".fastq.gz", ".raw.sam", fastq_file)) # Output sam file name
  print(line)
  system(line)
} 

# 3) Run alignment for each fastq file
sapply(fastq, align)
