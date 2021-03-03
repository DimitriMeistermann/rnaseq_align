conda activate rnaseq_align
snakemake -rp -j 30 --cluster 'qsub'  --jobscript ./sge.sh --latency-wait 10 --max-jobs-per-second 1 --configfile configFileSave/