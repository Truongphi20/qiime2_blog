import biom
import skbio
from skbio.diversity import beta_diversity
import pandas as pd
import numpy as np

def calculate_weighted_unifrac(table, tree):
    """
    Calculates the weighted UniFrac distance matrix.
    """
    # Prepare Table Data
    counts = table.matrix_data.toarray().T.astype(int)
    sample_ids = table.ids(axis='sample')
    otu_ids = table.ids(axis='observation')

    # Calculate Distance Matrix
    dm = beta_diversity(
        metric='weighted_unifrac',
        counts=counts,
        ids=sample_ids,
        taxa=otu_ids,
        tree=tree
    )

    # Convert to a readable DataFrame
    return dm.to_data_frame()

if __name__ == "__main__":
    BIOM_FILE = "/workspaces/qiime2_blog/support_data/feature-table-rarefied.biom"
    TREE_FILE = "/workspaces/qiime2_blog/support_data/tree.nwk"

    # Load data
    table = biom.load_table(BIOM_FILE)
    tree = skbio.TreeNode.read(TREE_FILE)
    
    # Calculate
    weighted_unifrac_df = calculate_weighted_unifrac(table, tree)
    
    print(weighted_unifrac_df.iloc[:5, :5])

# Expected output:
#             q2..L1S105  q2..L1S140  q2..L1S208  q2..L1S257  q2..L1S281
# q2..L1S105    0.000000    0.105481    0.101449    0.129749    0.137230
# q2..L1S140    0.105481    0.000000    0.103910    0.151287    0.156065
# q2..L1S208    0.101449    0.103910    0.000000    0.066404    0.083570
# q2..L1S257    0.129749    0.151287    0.066404    0.000000    0.080584
# q2..L1S281    0.137230    0.156065    0.083570    0.080584    0.000000