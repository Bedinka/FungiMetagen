import os
import pandas as pd
import subprocess

metadata_path = "/Users/baledi/Novobiom/FungiMetagen/LMS-16S-metadata.tsv"
fastq_directory = "/Users/baledi/Novobiom/FungiMetagen/paired-files"
fastqc_output_directory = "/Users/baledi/Novobiom/FungiMetagen/FastQC_output"
trimmed_output_directory = "/Users/baledi/Novobiom/FungiMetagen/Trimming_output"
forward_primer = "CCTACGGGNGGCWGCAG"
reverse_primer = "GACTACHVGGGTATCTAATCC"
min_length = 300


metadata = pd.read_csv(metadata_path, sep='\t')
print(metadata.head())

os.makedirs(fastqc_output_directory, exist_ok=True)
os.makedirs(trimmed_output_directory, exist_ok=True)

for fastq_file in os.listdir(fastq_directory):
    if fastq_file.endswith('.fastq.gz'):
        
        input_path = os.path.join(fastq_directory, fastq_file)
        output_path = os.path.join(fastqc_output_directory, fastq_file)
        trimmed_output_path = os.path.join(trimmed_output_directory, f"trimmed_{fastq_file}")
        unpaired_output_path = os.path.join(trimmed_output_directory, f"unpaired_{fastq_file}")

        subprocess.run([
            'fastqc', 
            input_path, 
            '--outdir', fastqc_output_directory
        ])

        print(f"FastQC completed for: {fastq_file}")

        sample_id = fastq_file.split('_')[0]  
        if sample_id in metadata['sample-id'].values:
            print(f"Matched {fastq_file} to sample ID: {sample_id}")
        else:
            print(f"Sample ID for {fastq_file} not found in metadata.")
        

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
"""
2: Specifies that a match is considered valid if at least 2 bases match the primer sequence.
30: Specifies the mismatch threshold.
10: Sets the minimum adapter length."""