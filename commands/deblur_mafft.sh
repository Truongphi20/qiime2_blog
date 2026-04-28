
## Interface of mafft to find the actually command
bashdb /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft --quiet \
        --preservecase \
        --parttree \
        --auto \
        --thread 1 \
        /workspaces/qiime2_blog/support_data/L6S68_30_L001_R1_001.fastq.gz.trim.derep.no_artifacts

## Real command: 
## because gdb dont be compatible with bash operator "<", run executor fist, then add args to run

# /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2307
# gdb /opt/conda/envs/qiime2-amplicon-2026.1/libexec/mafft/replaceu
# (gdb) set args -i orig > infile 2>>/dev/null

# /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2673
# gdb /workspaces/qiime2_blog/commands/mafft-7.525-with-extensions/build/core/disttbfast
# (gdb) set args -q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 1-1 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100 < infile

# /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2691
# gdb /opt/conda/envs/qiime2-amplicon-2026.1/libexec/mafft/dndpre
# (gdb) set args -y hat2 -b 62 -M 2 -C 1 < pre     > /dev/null 2>>/dev/null

# /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2693
# gdb /opt/conda/envs/qiime2-amplicon-2026.1/libexec/mafft/dvtditr
# (gdb) set args -W 0.00001 -E 0.0 -s 0.0 -C 1 -t 0 -F -z 50 -b 62 -f -1.53 -Q 100.0 -h 0 -I 2 -X 0.1 -p BAATARI2 -K 0 < pre     > /dev/null 2>>/dev/null

# /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2699
# gdb /opt/conda/envs/qiime2-amplicon-2026.1/libexec/mafft/restoreu
# (gdb) set args -a pre -i orig > restored

# /opt/conda/envs/qiime2-amplicon-2026.1/bin/mafft:2804
# gdb /workspaces/qiime2_blog/commands/mafft-7.525-with-extensions/build/core/f2cl
# (gdb) set args -n -1 -f -l 60 < /workspaces/qiime2_blog/support_data/mafft_tmp/pre