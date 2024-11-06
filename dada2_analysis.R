# Load required library
library(dada2)

# Define paths
run_directory <- commandArgs(trailingOnly = TRUE)[1]
trimmed_path <- file.path(run_directory, "trimmed_files")
output_path <- file.path(run_directory, "dada2_output")
dir.create(output_path, showWarnings = FALSE)

# List only paired forward and reverse reads (excluding unpaired files)
fnFs <- sort(list.files(trimmed_path, pattern="_R1_.*\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(trimmed_path, pattern="_R2_.*\\.fastq\\.gz$", full.names = TRUE))

# Exclude unpaired files by checking for corresponding pairs
paired_fnFs <- fnFs[sapply(fnFs, function(f) {
  paired_r <- sub("_R1_", "_R2_", f)
  paired_r %in% fnRs
})]

paired_fnRs <- fnRs[sapply(fnRs, function(r) {
  paired_f <- sub("_R2_", "_R1_", r)
  paired_f %in% fnFs
})]

# Check if the number of paired forward and reverse reads match
if (length(paired_fnFs) != length(paired_fnRs)) {
  stop("Mismatched number of forward and reverse read files.")
}

# Filter out sequences with ambiguous bases (N) before denoising
filt_path <- file.path(run_directory, "filtered_files")
dir.create(filt_path, showWarnings = FALSE)
filtFs <- file.path(filt_path, basename(paired_fnFs))
filtRs <- file.path(filt_path, basename(paired_fnRs))

# Apply filterAndTrim to remove ambiguous bases (maxN=0 removes sequences with Ns)
filterAndTrim(paired_fnFs, filtFs, paired_fnRs, filtRs, maxN=0, truncQ=2, multithread=TRUE)

# Learn error rates from the filtered reads
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Perform DADA2 denoising
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Save output tables and sequences
saveRDS(dadaFs, file.path(output_path, "dadaFs.rds"))
saveRDS(dadaRs, file.path(output_path, "dadaRs.rds"))
