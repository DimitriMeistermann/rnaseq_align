conda activate rnaseq_align
snakemake -rp -j 20 --cluster 'qsub'  --jobscript clusterFiles/jobscriptAtlas.sh --latency-wait 60 --max-jobs-per-second 1 --configfile configFileSave/[yourCounfigFile]
