import biom
import pandas as pd
import numpy as np

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_util.py:18
def _validate_counts_vector(counts, cast_int=False):
    """Validate and convert input to an acceptable counts vector type.

    Parameters
    ----------
    counts : array_like of int or float of shape (n_taxa,)
        Vector of counts.
    cast_int : bool, optional
        Cast values into integers, if not already. ``False`` by default.

    Returns
    -------
    ndarray of int or float of shape (n_taxa,)
        Valid counts vector.

    Raises
    ------
    ValueError
        If input array has an invalid data type.
    ValueError
        If input array is not 1-D.
    ValueError
        If there are negative values.

    Notes
    -----
    This function will return the original ``counts`` if it is already a valid counts
    vector. Otherwise it will return an edited copy that is valid.

    The data type of counts must be any subtype of ``np.integer`` (integers) or
    ``np.floating`` (floating-point numbers; excluding complex numbers) [1]_.

    See Also
    --------
    _validate_counts_matrix

    References
    ----------
    .. [1] https://numpy.org/doc/stable/reference/arrays.scalars.html

    """
    counts = np.asarray(counts)

    # counts must be int or float
    if np.issubdtype(dtype := counts.dtype, np.floating):
        # cast values into integers
        if cast_int:
            counts = counts.astype(int)

    elif not np.issubdtype(dtype, np.integer) and dtype is not np.dtype("bool"):
        raise ValueError("Counts must be integers or floating-point numbers.")

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")

    if (counts < 0).any():
        raise ValueError("Counts vector cannot contain negative values.")

    return counts

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/skbio/_methods.py:75
def _shannon(counts, base=2):
    counts = _validate_counts_vector(counts)
    freqs = counts / counts.sum()
    nonzero_freqs = freqs[freqs.nonzero()]
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/alpha.py:136
def shannon_entropy(table: biom.Table,
                    drop_undefined_samples: bool = False,
                    base: float = 2) -> pd.Series:
    if base == 'e':
        base = np.e

    if drop_undefined_samples:
        table = table.remove_empty(inplace=False)

    results = []
    for v in table.iter_data(dense=True):
        # using in-house metrics temporarily
        # results.append(_skbio_alpha_diversity_from_1d(v, 'shannon'))
        v = np.reshape(v, (1, len(v)))
        results.extend([_shannon(c, base=base)for c in v])
    results = pd.Series(results, index=table.ids(), name='shannon_entropy')
    return results



if __name__ == "__main__":
    # Mini table (10 ASVs x 10 samples) extracted from the "table.qza" of DADA2 output 
    table_path = "/workspaces/qiime2_blog/support_data/mini-feature-table.tsv"

    # Read table
    with open(table_path, 'r') as f:
        table = biom.Table.from_tsv(f, None, None, lambda x: x)

    # Compute
    results = shannon_entropy(table)
    print(results)