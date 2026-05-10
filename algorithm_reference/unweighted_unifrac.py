import biom
import skbio
from skbio.diversity import beta_diversity
import pandas as pd
import numpy as np

def calculate_unweighted_unifrac(table, tree):
    """
    Calculates the Unweighted UniFrac distance matrix.
    """
    # Prepare Table Data
    counts = table.matrix_data.toarray().T.astype(int)
    sample_ids = table.ids(axis='sample')
    otu_ids = table.ids(axis='observation')

    # Calculate Distance Matrix
    dm = beta_diversity(
        metric='unweighted_unifrac',
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
    unweighted_unifrac_df = calculate_unweighted_unifrac(table, tree)
    
    print(unweighted_unifrac_df.iloc[:5, :5])

# Expected output:
#             q2..L1S105  q2..L1S140  q2..L1S208  q2..L1S257  q2..L1S281
# q2..L1S105    0.000000    0.441940    0.545205    0.540621    0.575789
# q2..L1S140    0.441940    0.000000    0.523709    0.534961    0.461661
# q2..L1S208    0.545205    0.523709    0.000000    0.346487    0.320912
# q2..L1S257    0.540621    0.534961    0.346487    0.000000    0.393405
# q2..L1S281    0.575789    0.461661    0.320912    0.393405    0.000000