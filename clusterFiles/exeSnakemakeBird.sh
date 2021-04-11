conda activate rnaseq_align
snakemake -rp -j 20 --cluster 'qsub'  --jobscript clusterFiles/jobscriptBird.sh --latency-wait 10 --max-jobs-per-second 1 --configfile configFileSave/[yourCounfigFile]
