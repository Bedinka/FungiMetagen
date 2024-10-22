import qiime2
import os

input_path = "paired-files"  
output_sequences = "paired-end-sequences.qza"


os.system(f"""
qiime tools import \\
  --type 'SampleData[PairedEndSequencesWithQuality]' \\
  --input-path {input_path} \\
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \\
  --output-path {output_sequences}
""")

output_trimmed_sequences = "trimmed-seqs.qza"
cores = 4  


adapter_f = "GACTACHVGGGTATCTAATCC"
front_f = "CCTACGGGNGGCWGCAG"

os.system(f"""
qiime cutadapt trim-paired \\
    --i-demultiplexed-sequences {output_sequences} \\
    --p-cores {cores} \\
    --p-front-f "{front_f}" \\
    --p-front-r "{adapter_f}" \\
    --o-trimmed-sequences {output_trimmed_sequences}
""")

trim_left_f = 0
trim_left_r = 0
trunc_len_f = 220
trunc_len_r = 180

table_output = "table.qza"
rep_seqs_output = "rep-seqs.qza"
denoising_stats_output = "denoising-stats.qza"

os.system(f"""
qiime dada2 denoise-paired \\
    --i-demultiplexed-seqs {output_trimmed_sequences} \\
    --p-trim-left-f {trim_left_f} \\
    --p-trim-left-r {trim_left_r} \\
    --p-trunc-len-f {trunc_len_f} \\
    --p-trunc-len-r {trunc_len_r} \\
    --p-n-threads {cores} \\
    --o-table {table_output} \\
    --o-representative-sequences {rep_seqs_output} \\
    --o-denoising-stats {denoising_stats_output}
""")

os.system(f"""
qiime metadata tabulate \\
    --m-input-file {denoising_stats_output} \\
    --o-visualization denoising-stats.qzv
""")

os.system(f"""
qiime feature-table summarize \\
    --i-table {table_output} \\
    --o-visualization table.qzv
""")

os.system(f"""
qiime feature-table tabulate-seqs \\
    --i-data {rep_seqs_output} \\
    --o-visualization rep-seqs.qzv
""")
