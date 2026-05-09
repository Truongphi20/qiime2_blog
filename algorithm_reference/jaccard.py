import biom
import pandas as pd
import numpy as np

# commands/scipy_7dcd8c5_src/distance_impl.h:706
def dist_to_squareform_from_vector_double(M_flat, X, d):
    v = 0
    
    for i in range(1, d):
        # The Horizontal Fill
        it1 = (i - 1) * d + i 
        length = d - i
        M_flat[it1 : it1 + length] = X[v : v + length]
        
        # The Vertical Fill (Symmetry)
        it2 = i * (d + 1) - 1
        
        for _ in range(i, d):
            M_flat[it2] = X[v]
            v += 1
            it2 += d 

    return M_flat.reshape(d, d)

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/scipy/spatial/distance.py:2196
def squareform(X, force="no", checks=True):
    s = X.shape

    # Grab the closest value to the square root of the number
    # of elements times 2 to see if the number of elements
    # is indeed a binomial coefficient.
    d = int(np.ceil(np.sqrt(s[0] * 2)))

    # Allocate memory for the distance matrix.
    M = np.zeros(d**2, dtype=X.dtype)

    # Fill in the values of the distance matrix.
    M = dist_to_squareform_from_vector_double(M, X, d)

    # Return the distance matrix.
    return M

def jaccard_distance(u, v):
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.jaccard.html
    
    num = np.sum((u != v) & (u | v))
    denom = np.sum(u | v)
    
    # Return 0 if both vectors are all zeros
    return num / denom if denom != 0 else 0.0

# commands/scipy_7dcd8c5_src/distance_pybind.cpp:485
def pdist_jaccard(X: np.array):
    X = np.asanyarray(X).astype(bool)
    n = X.shape[0]
    
    # Calculate the size of the condensed distance matrix: nC2
    out_size = n * (n - 1) // 2
    dm = np.zeros(out_size, dtype=np.double)
    
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            dm[k] = jaccard_distance(X[i], X[j])
            k += 1
            
    return dm

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/beta.py:187
def jaccard(table: biom.Table, n_jobs: int = 1) -> pd.DataFrame:
    counts = table.matrix_data.toarray().T.copy()
    sample_ids = table.ids(axis='sample')
    
    counts = (counts > 0.0)
    distances = squareform(pdist_jaccard(counts))

    jaccard_table = pd.DataFrame(distances, columns=sample_ids, index=sample_ids)
    return jaccard_table

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