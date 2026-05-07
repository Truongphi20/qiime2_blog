import biom
import pandas as pd
import skbio
import numpy as np

def _skbio_alpha_diversity_from_1d(v, metric):
    # alpha_diversity expects a 2d structure
    v = np.reshape(v, (1, len(v)))
    result = skbio.diversity.alpha_diversity(metric=metric,
                                             counts=v,
                                             ids=['placeholder', ],
                                             validate=False)
    return result.iloc[0]

def observed_features(table: biom.Table) -> pd.Series:
    # Convert the table to presence/absence data
    presence_absence_table = table.pa(inplace=False)
    
    results = []
    # Loop through each sample
    for v in presence_absence_table.iter_data():          
        results.append(_skbio_alpha_diversity_from_1d(v.astype(int),
                                                      'observed_otus'))
    
    results = pd.Series(results, index=table.ids(), name='observed_features')
    return results





if __name__ == "__main__":
    # Mini table (10 ASVs x 10 samples) extracted from the "table.qza" of DADA2 output 
    table_path = "/workspaces/qiime2_blog/support_data/mini-feature-table.tsv"

    # Read table
    with open(table_path, 'r') as f:
        table = biom.Table.from_tsv(f, None, None, lambda x: x)

    # Compute
    results = observed_features(table)
    print(results)