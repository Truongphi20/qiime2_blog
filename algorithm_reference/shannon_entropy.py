import biom
import pandas as pd
import numpy as np

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/skbio/_methods.py:75
def _shannon(counts, base=2):
    freqs = counts / counts.sum()
    nonzero_freqs = freqs[freqs.nonzero()]
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/alpha.py:136
def shannon_entropy(table: biom.Table,
                    base: float = 2) -> pd.Series:
    if base == 'e':
        base = np.e

    results = []
    for v in table.iter_data(dense=True):
        v = np.reshape(v, (1, len(v)))
        results.extend([_shannon(c, base=base)for c in v])
    results = pd.Series(results, index=table.ids(), name='shannon_entropy')
    return results



if __name__ == "__main__":
    # Mini table (10 ASVs x 11 samples) extracted from the "table.qza" of DADA2 output 
    table_path = "/workspaces/qiime2_blog/support_data/mini-feature-table.tsv"

    # Read table
    with open(table_path, 'r') as f:
        table = biom.Table.from_tsv(f, None, None, lambda x: x)

    # Compute
    results = shannon_entropy(table)
    print(results)