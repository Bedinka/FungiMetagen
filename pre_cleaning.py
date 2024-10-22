import os
import pandas as pd
import subprocess
from datetime import datetime
import logging

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
run_directory = f"/Users/baledi/Novobiom/FungiMetagen/run_{timestamp}"
os.makedirs(run_directory, exist_ok=True)
log_file = os.path.join(run_directory, "run_log.txt")


metadata_path = "/Users/baledi/Novobiom/FungiMetagen/LMS-16S-metadata.tsv"
fastq_directory = "/Users/baledi/Novobiom/FungiMetagen/paired-files"

fastqc_output_directory = os.path.join(run_directory, "fastqc_output")
trimmed_output_directory = os.path.join(run_directory, "trimmed_files")
post_trim_fastqc_directory = os.path.join(run_directory, "post_trim_fastqc")
os.makedirs(fastqc_output_directory, exist_ok=True)
os.makedirs(trimmed_output_directory, exist_ok=True)
os.makedirs(post_trim_fastqc_directory, exist_ok=True)

forward_primer = "CCTACGGGNGGCWGCAG"
reverse_primer = "GACTACHVGGGTATCTAATCC"
min_length = 300

logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logging.info(f"Starting pipeline run at {timestamp}")

metadata = pd.read_csv(metadata_path, sep='\t')
print(metadata.head())

logging.info("Metadata loaded successfully.")
logging.info(f"Metadata head: \n{metadata.head()}")

for fastq_file in os.listdir(fastq_directory):
    if fastq_file.endswith('.fastq.gz'):
        
        input_path = os.path.join(fastq_directory, fastq_file)
        fastqc_output_path = os.path.join(fastqc_output_directory, fastq_file)
        trimmed_output_path = os.path.join(trimmed_output_directory, f"trimmed_{fastq_file}")
        unpaired_output_path = os.path.join(trimmed_output_directory, f"unpaired_{fastq_file}")
        sample_id = fastq_file.split('_L001')[0]  
        print(f"Extracted sample ID: {sample_id} for file: {fastq_file}")
        logging.info(f"Started FastQC analysis for {fastq_file}.")
        subprocess.run([
            'fastqc', 
            input_path, 
            '--outdir', fastqc_output_directory
        ])
        logging.info(f"FastQC completed for: {fastq_file}")

        if sample_id in metadata['sample-id'].values:
            logging.info(f"Matched {fastq_file} to sample ID: {sample_id}")
        else:
            logging.warning(f"Sample ID for {fastq_file} not found in metadata.")


        print(f"FastQC completed for: {fastq_file}")

        sample_id = fastq_file.split('_L001')[0]  
        print(f"Extracted sample ID: {sample_id} for file: {fastq_file}")
        if sample_id in metadata['sample-id'].values:
            logging.info(f"Matched {fastq_file} to sample ID: {sample_id}")
        else:
            logging.warning(f"Sample ID for {fastq_file} not found in metadata.")

        r2_file = input_path.replace("_R1_", "_R2_")
        if not os.path.exists(r2_file):
            logging.warning(f"Paired R2 file missing for {fastq_file}. Skipping this pair.")
            continue
        logging.info(f"Started trimming for {fastq_file} with Trimmomatic.")
        subprocess.run([
            'trimmomatic',
            'PE',  # Specify paired-end mode
            '-phred33',  # Quality encoding (phred33 is most common)
            input_path,
            input_path.replace("_R1_", "_R2_"),  # Adjust for paired-end naming convention
            trimmed_output_path,  # Output file for paired trimmed reads (R1)
            unpaired_output_path.replace("_R1_", "_R1_unpaired_"),  # Output file for unpaired reads (R1)
            trimmed_output_path.replace("_R1_", "_R2_"),  # Output file for paired trimmed reads (R2)
            unpaired_output_path.replace("_R1_", "_R2_unpaired_"),  # Output file for unpaired reads (R2)
            'ILLUMINACLIP:adapters.fa:2:30:10',  # Adapter clipping (you'll need an adapter file)
            'SLIDINGWINDOW:4:20',  # Sliding window trimming (4 bases with average quality below 20)
            f'MINLEN:{min_length}'  # Discard reads shorter than 300 bp
        ])
        print(f"Trimming completed for: {fastq_file}, output saved to: {trimmed_output_path}")
        logging.info(f"Trimming completed for: {fastq_file}, output saved to: {trimmed_output_path}")

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

        logging.info("Running MultiQC on original FastQC reports.")
        
subprocess.run([
    'multiqc',
    fastqc_output_directory,
    '--outdir', run_directory
])
logging.info("MultiQC report for original files generated.")

# Run MultiQC on the post-trimming FastQC reports
logging.info("Running MultiQC on post-trimming FastQC reports.")
subprocess.run([
    'multiqc',
    post_trim_fastqc_directory,
    '--outdir', run_directory
])
logging.info("MultiQC report for trimmed files generated.")

logging.info("Pipeline run completed.")
"""
2: Specifies that a match is considered valid if at least 2 bases match the primer sequence.
30: Specifies the mismatch threshold.
10: Sets the minimum adapter length."""