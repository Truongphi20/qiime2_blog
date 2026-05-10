import biom
import skbio
from skbio.diversity import alpha
import pandas as pd
import numpy as np

def calculate_faith_pd(biom_path, tree_path):
    # 1. Load the feature table
    table = biom.load_table(biom_path)
    # Convert to dense array: rows are ASVs, columns are Samples
    counts = table.matrix_data.toarray().T.astype(int) 
    sample_ids = table.ids(axis='sample')
    otu_ids = table.ids(axis='observation')

    # 2. Load the phylogenetic tree
    tree = skbio.TreeNode.read(tree_path)

    # 3. Calculate Faith's PD for each sample
    # skbio.diversity.alpha.faith_pd requires:
    # - counts: a 1D array of abundances
    # - otu_ids: a list of labels for those abundances
    # - tree: the skbio TreeNode object
    
    pd_results = []
    for i, sample_id in enumerate(sample_ids):
        sample_counts = counts[i]
        
        # Calculate the metric
        # Note: If an ASV in the table is missing from the tree, skbio will raise an error.
        result = alpha.faith_pd(sample_counts, otu_ids, tree)
        pd_results.append(result)

    # 4. Format as a Pandas Series
    return pd.Series(pd_results, index=sample_ids, name='faith_pd')

if __name__ == "__main__":
    BIOM_FILE = "/workspaces/qiime2_blog/support_data/feature-table-rarefied.biom"
    TREE_FILE = "/workspaces/qiime2_blog/support_data/tree.nwk"
    
    faith_pd_series = calculate_faith_pd(BIOM_FILE, TREE_FILE)
    print(faith_pd_series)