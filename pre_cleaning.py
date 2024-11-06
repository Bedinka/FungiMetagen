import os
import pandas as pd
import subprocess
from datetime import datetime
import logging

metadata_path = "/Users/baledi/Novobiom/FungiMetagen/LMS-16S-metadata.tsv"
fastq_directory = "/Users/baledi/Novobiom/FungiMetagen/paired-files"
input_path = "paired-files"
os.makedirs(input_path, exist_ok=True)

# Create a timestamped directory for this run
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
run_directory = f"/Users/baledi/Novobiom/FungiMetagen/runs/run_{timestamp}"
os.makedirs(run_directory, exist_ok=True)

# Define subdirectories for outputs
fastqc_output_directory = os.path.join(run_directory, "fastqc_output")
trimmed_output_directory = os.path.join(run_directory, "trimmed_files")
post_trim_fastqc_directory = os.path.join(run_directory, "post_trim_fastqc")
os.makedirs(fastqc_output_directory, exist_ok=True)
os.makedirs(trimmed_output_directory, exist_ok=True)
os.makedirs(post_trim_fastqc_directory, exist_ok=True)

# Define path for the log file
log_file = os.path.join(run_directory, "run_log.txt")

# Setup logging
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logging.info(f"Starting pipeline run at {timestamp}")

# Primer sequences (for adapter trimming)
forward_primer = "CCTACGGGNGGCWGCAG"
reverse_primer = "GACTACHVGGGTATCTAATCC"
min_length = 300

# Load metadata
metadata = pd.read_csv(metadata_path, sep='\t')
logging.info("Metadata loaded successfully.")
logging.info(f"Metadata head: \n{metadata.head()}")

# Iterate through FASTQ files, run FastQC, and trim with Trimmomatic
for fastq_file in os.listdir(fastq_directory):
    if fastq_file.endswith('.fastq.gz'):
        # Define input and output paths for FastQC
        input_path = os.path.join(fastq_directory, fastq_file)
        fastqc_output_path = os.path.join(fastqc_output_directory, fastq_file)
        trimmed_output_path = os.path.join(trimmed_output_directory, f"trimmed_{fastq_file}")
        unpaired_output_path = os.path.join(trimmed_output_directory, f"unpaired_{fastq_file}")

        logging.info(f"Started FastQC analysis for {fastq_file}.")
        # Run FastQC on the original file
        subprocess.run([
            'fastqc', 
            input_path, 
            '--outdir', fastqc_output_directory
        ])
        logging.info(f"FastQC completed for: {fastq_file}")

        # Extract sample ID from file name for matching with metadata
        sample_id = fastq_file.split('_L001')[0]
        logging.info(f"Extracted sample ID: {sample_id} for file: {fastq_file}")

        # Check if paired R2 file exists
        r2_file = input_path.replace("_R1_", "_R2_")
        if not os.path.exists(r2_file):
            logging.warning(f"Paired R2 file missing for {fastq_file}. Skipping this pair.")
            continue

        logging.info(f"Started trimming for {fastq_file} with Trimmomatic.")
        # Run Trimmomatic for trimming primers and reads
        subprocess.run([
            'trimmomatic',
            'PE',  # Specify paired-end mode
            '-phred33',  # Quality encoding (phred33 is most common)
            input_path, r2_file,
            trimmed_output_path,  # Output file for paired trimmed reads (R1)
            unpaired_output_path.replace("_R1_", "_R1_unpaired_"),  # Output file for unpaired reads (R1)
            trimmed_output_path.replace("_R1_", "_R2_"),  # Output file for paired trimmed reads (R2)
            unpaired_output_path.replace("_R1_", "_R2_unpaired_"),  # Output file for unpaired reads (R2)
            'ILLUMINACLIP:adapters.fa:2:30:10',  # Adapter clipping (you'll need an adapter file)
            'SLIDINGWINDOW:4:20',  # Sliding window trimming (4 bases with average quality below 20)
            f'MINLEN:{min_length}'  # Discard reads shorter than 300 bp
        ])
        logging.info(f"Trimming completed for: {fastq_file}, output saved to: {trimmed_output_path}")

        # Run FastQC on the trimmed output files (R1 and R2)
        trimmed_r1_path = trimmed_output_path
        trimmed_r2_path = trimmed_output_path.replace("_R1_", "_R2_")
        
        logging.info(f"Started FastQC analysis for trimmed files of {fastq_file}.")
        subprocess.run([
            'fastqc',
            trimmed_r1_path,
            '--outdir', post_trim_fastqc_directory
        ])
        subprocess.run([
            'fastqc',
            trimmed_r2_path,
            '--outdir', post_trim_fastqc_directory
        ])
        logging.info(f"FastQC completed for trimmed files of: {fastq_file}")

# Run MultiQC on the original FastQC reports after processing all files
logging.info("Running MultiQC on original FastQC reports.")
subprocess.run([
    'multiqc',
    fastqc_output_directory,
    '--outdir', run_directory
])
logging.info("MultiQC report for original files generated.")

# Run MultiQC on the post-trimming FastQC reports after processing all files
logging.info("Running MultiQC on post-trimming FastQC reports.")
subprocess.run([
    'multiqc',
    post_trim_fastqc_directory,
    '--outdir', run_directory
])
logging.info("MultiQC report for trimmed files generated.")

# R script for DADA2 analysis (assumes you have an R script set up)
dada2_script_path = "/Users/baledi/Novobiom/FungiMetagen/dada2_analysis.R"
logging.info(f"Running DADA2 denoising using R script at {dada2_script_path}")
subprocess.run(["Rscript", dada2_script_path, run_directory])
logging.info("DADA2 denoising completed.")

logging.info("Pipeline run completed.")

