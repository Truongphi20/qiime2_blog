import biom
import pandas as pd
import numpy as np
import skbio
import sklearn.metrics
from skbio.stats.distance import DistanceMatrix

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_driver.py:367
def beta_diversity(metric, counts, ids=None, pairwise_func=None):
    counts = (counts > 0.0)
    distances = pairwise_func(counts, metric=metric)
    return DistanceMatrix(distances, ids)

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/beta.py:187
def jaccard(table: biom.Table, n_jobs: int = 1) -> skbio.DistanceMatrix:
    counts = table.matrix_data.toarray().T.copy()
    sample_ids = table.ids(axis='sample')
    jaccard_table = beta_diversity(
            metric='jaccard', 
            counts=counts, 
            ids=sample_ids,
            pairwise_func=sklearn.metrics.pairwise_distances
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

# Expected output:
#          L1S105    L1S140    L1S208    L1S257    L1S281     L1S57     L1S76      L1S8    L2S155    L2S175    L2S204
# L1S105  0.000000  0.600000  0.714286  0.666667  0.600000  0.250000  0.250000  0.500000  0.857143  0.875000  0.666667
# L1S140  0.600000  0.000000  0.400000  0.250000  0.000000  0.500000  0.500000  0.666667  1.000000  1.000000  1.000000
# L1S208  0.714286  0.400000  0.000000  0.200000  0.400000  0.666667  0.666667  0.571429  0.714286  0.750000  0.875000
# L1S257  0.666667  0.250000  0.200000  0.000000  0.250000  0.600000  0.600000  0.714286  0.857143  0.875000  0.857143
# L1S281  0.600000  0.000000  0.400000  0.250000  0.000000  0.500000  0.500000  0.666667  1.000000  1.000000  1.000000
# L1S57   0.250000  0.500000  0.666667  0.600000  0.500000  0.000000  0.000000  0.400000  0.833333  0.857143  0.833333
# L1S76   0.250000  0.500000  0.666667  0.600000  0.500000  0.000000  0.000000  0.400000  0.833333  0.857143  0.833333
# L1S8    0.500000  0.666667  0.571429  0.714286  0.666667  0.400000  0.400000  0.000000  0.714286  0.571429  0.875000
# L2S155  0.857143  1.000000  0.714286  0.857143  1.000000  0.833333  0.833333  0.714286  0.000000  0.200000  0.400000
# L2S175  0.875000  1.000000  0.750000  0.875000  1.000000  0.857143  0.857143  0.571429  0.200000  0.000000  0.500000
# L2S204  0.666667  1.000000  0.875000  0.857143  1.000000  0.833333  0.833333  0.875000  0.400000  0.500000  0.000000