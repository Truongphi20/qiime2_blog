import biom
import skbio
from skbio.diversity import alpha
import pandas as pd
import numpy as np

def calculate_faith_pd(table: biom.Table, tree: skbio.TreeNode):
    
    # Convert to dense array: rows are ASVs, columns are Samples
    counts = table.matrix_data.toarray().T.astype(int) 
    sample_ids = table.ids(axis='sample')
    otu_ids = table.ids(axis='observation')

    # Calculate Faith's PD for each sample
    pd_results = []
    for i in range(len(sample_ids)):
        sample_counts = counts[i]
        
        # Calculate the metric
        result = alpha.faith_pd(sample_counts, otu_ids, tree)
        pd_results.append(result)

    return pd.Series(pd_results, index=sample_ids, name='faith_pd')

if __name__ == "__main__":
    BIOM_FILE = "/workspaces/qiime2_blog/support_data/feature-table-rarefied.biom"
    TREE_FILE = "/workspaces/qiime2_blog/support_data/tree.nwk"

    # Load the feature table
    table = biom.load_table(BIOM_FILE)

    # Load the phylogenetic tree
    tree = skbio.TreeNode.read(TREE_FILE)
    
    faith_pd_series = calculate_faith_pd(table, tree)
    print(faith_pd_series)

    # Expected output:
    # q2..L1S105     6.384621
    # q2..L1S140     5.900037
    # q2..L1S208     6.700060
    # q2..L1S257     7.362996
    # q2..L1S281     5.489152
    # q2..L1S57      5.572742
    # q2..L1S76      6.159797
    # q2..L1S8       4.973571
    # q2..L2S155    15.671622
    # q2..L2S175    13.210765
    # q2..L2S204    15.803438
    # q2..L2S222    18.029198
    # q2..L2S240     8.010033
    # q2..L2S309    13.595307
    # q2..L2S357     8.792001
    # q2..L2S382    12.776710
    # q2..L3S294     8.676022
    # q2..L3S313     7.818368
    # q2..L3S378     4.646372
    # q2..L4S112    15.468081
    # q2..L4S137    13.628922
    # q2..L4S63     21.770786
    # q2..L5S104     5.608632
    # q2..L5S155     4.361857
    # q2..L5S174     5.357451
    # q2..L5S203     5.903046
    # q2..L5S222     5.311506
    # q2..L5S240     5.060430
    # q2..L6S20      6.737120
    # q2..L6S68      6.374415
    # q2..L6S93      6.317089