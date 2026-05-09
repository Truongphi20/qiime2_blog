import biom
import pandas as pd
import numpy as np

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/skbio/_methods.py:75
def _shannon(counts, base=2):
    freqs = counts / counts.sum()
    nonzero_freqs = freqs[freqs.nonzero()]
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/skbio/_methods.py:69
def _p_evenness(counts):
    return _shannon(counts, base=np.e) / np.log(np.count_nonzero(counts))

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/alpha.py:114
def pielou_evenness(table: biom.Table) -> pd.Series:

    results = []
    for v in table.iter_data(dense=True):
        v = np.reshape(v, (1, len(v)))
        results.extend([_p_evenness(c)for c in v])
    results = pd.Series(results, index=table.ids(), name='pielou_evenness')
    return results


if __name__ == "__main__":
    # Mini table (10 ASVs x 11 samples) extracted from the "table.qza" of DADA2 output 
    table_path = "/workspaces/qiime2_blog/support_data/mini-feature-table.tsv"

    # Read table
    with open(table_path, 'r') as f:
        table = biom.Table.from_tsv(f, None, None, lambda x: x)

    # Compute
    results = pielou_evenness(table)
    print(results)