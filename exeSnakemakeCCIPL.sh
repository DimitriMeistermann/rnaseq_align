snakemake -rp -j 30 --cluster 'sbatch -o ./log/out -e ./log/err -p Bird -w bigmem001 -n {cluster.cpu} --mem={cluster.mem}' \
--cluster-config configFileSave/cluster.json  --latency-wait 10 --max-jobs-per-second 1 --configfile configFileSave/configTest.json
