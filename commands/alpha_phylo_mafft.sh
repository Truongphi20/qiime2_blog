bashdb /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft \
                --preservecase \
                --inputorder --thread 1 /workspaces/qiime2_blog/support_data/dna-sequences.fasta

# (/opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2307)
# replaceu -i orig > infile 2>>"/dev/stderr"

# (/opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2673)
# disttbfast -q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 1-1 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100 < /workspaces/qiime2_blog/support_data/alpha_phy_mafft/infile   > pre 2>>"/dev/stderr"

# (/opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2699)
# restoreu -a pre -i orig > restored

# (/opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2804)
# f2cl -n -1 -f -l 60 < /workspaces/qiime2_blog/support_data/alpha_phy_mafft/pre 2>"/dev/null"