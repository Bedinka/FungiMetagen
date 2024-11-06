import qiime2
import os

input_path = "paired-files"
output_sequences = "paired-end-sequences.qza"
output_trimmed_sequences = "trimmed-seqs.qza"
cores = 24
adapter_f = "GACTACHVGGGTATCTAATCC"
front_f = "CCTACGGGNGGCWGCAG"
trim_left_f = 0
trim_left_r = 0
trunc_len_f = 240
trunc_len_r = 220
table_output = "table.qza"
rep_seqs_output = "rep-seqs.qza"
denoising_stats_output = "denoising-stats.qza"
taxonomy_output = "taxonomy.qza"
taxonomy_visualization = "taxonomy.qzv"
classifier_path = "/Users/baledi/Novobiom/FungiMetagen/silva-99-taxonomy.qza"

os.system(f"""
qiime tools import \\
  --type 'SampleData[PairedEndSequencesWithQuality]' \\
  --input-path {input_path} \\
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \\
  --output-path {output_sequences}
""")

os.system(f"""
qiime demux summarize \\
  --i-data {output_sequences} \\
  --o-visualization demux.qzv
""")

os.system(f"""
qiime cutadapt trim-paired \\
  --i-demultiplexed-sequences {output_sequences} \\
  --p-cores {cores} \\
  --p-front-f "{front_f}" \\
  --p-front-r "{adapter_f}" \\
  --o-trimmed-sequences {output_trimmed_sequences}
""")


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

os.makedirs("phyloseq", exist_ok=True)
os.system(f"""
qiime tools export \\
  --input-path {table_output} \\
  --output-path phyloseq
""")

os.system(f"""
biom convert \\
  -i phyloseq/feature-table.biom \\
  -o phyloseq/otu_table.tsv \\
  --to-tsv
""")

os.system("""
sed -i '' '1d' phyloseq/otu_table.tsv
sed -i '' 's/#OTU ID//' phyloseq/otu_table.tsv
""")

os.system(f"""
qiime tools export \\
  --input-path {rep_seqs_output} \\
  --output-path phyloseq
""")

if os.path.exists(rep_seqs_output):
    os.system(f"""
    qiime feature-classifier classify-sklearn \\
      --i-classifier {classifier_path} \\
      --i-reads {rep_seqs_output} \\
      --o-classification {taxonomy_output}
    """)
    if os.path.exists(taxonomy_output):
        os.system(f"""
        qiime metadata tabulate \\
          --m-input-file {taxonomy_output} \\
          --o-visualization {taxonomy_visualization}
        """)
    else:
        print("Taxonomy classification output not found.")
else:
    print("Representative sequences not found. Check DADA2 step.")
