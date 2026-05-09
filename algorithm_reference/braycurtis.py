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
def squareform(X):
    s = X.shape

    # Grab the closest value to the square root of the number
    # of elements times 2 to see if the number of elements
    # is indeed a binomial coefficient.
    d = int(np.ceil(np.sqrt(s[0] * 2)))

    # Allocate memory for the distance matrix.
    M = np.zeros(d**2, dtype=X.dtype)
    
    # Fill in the values of the distance matrix.
    M = dist_to_squareform_from_vector_double(M, X, d)

    return M

def braycurtis(u: np.array, v:np.array):
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.braycurtis.html

    num = np.abs(u-v).sum()
    denom = np.abs(u+v).sum()

    return num/denom

def pdist_braycurtis(X: np.array):

    n = X.shape[0]

    # Calculate the size of the condensed distance matrix: nC2
    out_size = n * (n - 1) // 2
    dm = np.zeros(out_size, dtype=np.double)
    
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            dm[k] = braycurtis(X[i], X[j])
            k += 1
            
    return dm

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/beta.py:172
def bray_curtis(table: biom.Table):
    counts = table.matrix_data.toarray().T.copy()
    sample_ids = table.ids(axis='sample')

    distances = squareform(pdist_braycurtis(counts))
    beta_diversity = pd.DataFrame(distances, columns=sample_ids, index=sample_ids)

    return beta_diversity


if __name__ == "__main__":
    # Mini table (10 ASVs x 11 samples) extracted from the "table.qza" of DADA2 output 
    table_path = "/workspaces/qiime2_blog/support_data/mini-feature-table.tsv"

    # Read table
    with open(table_path, 'r') as f:
        table = biom.Table.from_tsv(f, None, None, lambda x: x)

    # Compute
    results = bray_curtis(table)
    print(results)

# Expected output:
#           L1S105    L1S140    L1S208    L1S257    L1S281     L1S57     L1S76      L1S8    L2S155    L2S175    L2S204
# L1S105  0.000000  0.648485  0.583208  0.606228  0.768390  0.199953  0.243231  0.255174  0.995701  0.996196  0.987325
# L1S140  0.648485  0.000000  0.121081  0.321053  0.239909  0.642618  0.707509  0.668111  1.000000  1.000000  1.000000
# L1S208  0.583208  0.121081  0.000000  0.203825  0.204274  0.725427  0.753743  0.771498  0.994821  0.995467  0.997010
# L1S257  0.606228  0.321053  0.203825  0.000000  0.219063  0.758565  0.790724  0.813720  0.995288  0.996085  0.994885
# L1S281  0.768390  0.239909  0.204274  0.219063  0.000000  0.806679  0.837012  0.860712  1.000000  1.000000  1.000000
# L1S57   0.199953  0.642618  0.725427  0.758565  0.806679  0.000000  0.082817  0.106373  0.995811  0.996282  0.989868
# L1S76   0.243231  0.707509  0.753743  0.790724  0.837012  0.082817  0.000000  0.166556  0.996036  0.996460  0.990439
# L1S8    0.255174  0.668111  0.771498  0.813720  0.860712  0.106373  0.166556  0.000000  0.993928  0.993509  0.990227
# L2S155  0.995701  1.000000  0.994821  0.995288  1.000000  0.995811  0.996036  0.993928  0.000000  0.390575  0.473239
# L2S175  0.996196  1.000000  0.995467  0.996085  1.000000  0.996282  0.996460  0.993509  0.390575  0.000000  0.694297
# L2S204  0.987325  1.000000  0.997010  0.994885  1.000000  0.989868  0.990439  0.990227  0.473239  0.694297  0.000000