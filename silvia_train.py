import os

version = "138.2"
identity_threshold = 0.99

os.system(f"""
qiime rescript get-silva-data \\
    --p-version {version} \\
    --p-target "SSURef_NR99" \\
    --p-include-species-labels \\
    --o-silva-sequences silva-sequences.qza \\
    --o-silva-taxonomy silva-taxonomy.qza
""")
os.system(f"""
qiime rescript cull-seqs \\
    --i-sequences silva-sequences.qza \\
    --o-clean-sequences silva-clean-seqs.qza
""")
# Dereplicate at 99% similarity to reduce redundancy
os.system(f"""
qiime rescript dereplicate \\
    --i-sequences silva-clean-seqs.qza \\
    --i-taxa silva-taxonomy.qza \\
    --p-rank-handles "silva" \\
    --o-dereplicated-sequences silva-99-seqs.qza \\
    --o-dereplicated-taxa silva-99-taxonomy.qza
""")
os.system(f"""
qiime feature-classifier fit-classifier-naive-bayes \\
    --i-reference-reads silva-99-seqs.qza \\
    --i-reference-taxonomy silva-99-taxonomy.qza \\
    --o-classifier silva-99-classifier.qza
""")
