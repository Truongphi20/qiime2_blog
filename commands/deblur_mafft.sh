
## Interface of mafft to find the actually command
# bashdb /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft --quiet \
#         --preservecase \
#         --parttree \
#         --auto \
#         --thread 1 \
#         /workspaces/qiime2_blog/support_data/L6S68_30_L001_R1_001.fastq.gz.trim.derep.no_artifacts

## Real command: /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2804
## because gdb dont be compatible with bash operator "<", run executor fist, then add args to run 
gdb /workspaces/qiime2_blog/commands/mafft-7.525-with-extensions/build/core/f2cl

## ADD to run: (gdb) set args -n -1 -f -l 60 < /workspaces/qiime2_blog/support_data/mafft_tmp/pre