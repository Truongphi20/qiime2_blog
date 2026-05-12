from biom import load_table

# 1. Load the table from DADA2 biom table
table = load_table('/workspaces/qiime2_blog/support_data/feature-table.biom')

# 4. Subset first 10 observes and 11 samples
with open('mini-feature-table.tsv', 'w') as f:
    f.write(table.head(10,11).to_tsv())