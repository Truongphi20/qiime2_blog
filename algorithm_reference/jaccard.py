import biom
import pandas as pd
import numpy as np
import skbio
import sklearn.metrics
import skbio.diversity

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/beta.py:187
def jaccard(table: biom.Table, n_jobs: int = 1) -> skbio.DistanceMatrix:
    counts = table.matrix_data.toarray().T.copy()
    sample_ids = table.ids(axis='sample')
    jaccard_table = skbio.diversity.beta_diversity(
        metric='jaccard',
        counts=counts,
        ids=sample_ids,
        validate=False,
        pairwise_func=sklearn.metrics.pairwise_distances,
        n_jobs=n_jobs
    )
    return jaccard_table.to_data_frame()

if __name__ == "__main__":
    # Mini table (10 ASVs x 11 samples) extracted from the "table.qza" of DADA2 output 
    table_path = "/workspaces/qiime2_blog/support_data/mini-feature-table.tsv"

    # Read table
    with open(table_path, 'r') as f:
        table = biom.Table.from_tsv(f, None, None, lambda x: x)

    # Compute
    results = jaccard(table)
    print(results)