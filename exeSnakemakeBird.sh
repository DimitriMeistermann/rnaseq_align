conda activate rnaseq_align
snakemake -rp -j 20 --cluster 'qsub'  --jobscript ./sge.sh --latency-wait 10 --max-jobs-per-second 1 --configfile configFileSave/

snakemake -rp -j 20 --cluster 'qsub'  --jobscript ./sgeAtlas.sh --latency-wait 10 --max-jobs-per-second 1