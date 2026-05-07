import biom
import pandas as pd
import numpy as np

# Reference: /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/alpha.py:92
def _skbio_alpha_diversity_from_1d(v):
    # alpha_diversity expects a 2d structure
    v = np.reshape(v, (1, len(v)))

    results = pd.Series([c.sum() for c in v], index=['placeholder', ])
    return results.iloc[0]

# Reference: /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/alpha.py:103
def observed_features(table: biom.Table) -> pd.Series:
    # Convert the table to presence/absence data
    presence_absence_table = table.pa(inplace=False)
    
    results = []
    # Loop through each sample
    for v in presence_absence_table.iter_data():          
        results.append(_skbio_alpha_diversity_from_1d(v.astype(int)))
    
    results = pd.Series(results, index=table.ids(), name='observed_features')
    return results





if __name__ == "__main__":
    # Mini table (10 ASVs x 11 samples) extracted from the "table.qza" of DADA2 output 
    table_path = "/workspaces/qiime2_blog/support_data/mini-feature-table.tsv"

    # Read table
    with open(table_path, 'r') as f:
        table = biom.Table.from_tsv(f, None, None, lambda x: x)

    # Compute
    results = observed_features(table)
    print(results)