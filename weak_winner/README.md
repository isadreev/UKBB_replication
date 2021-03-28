# Perform simulation analyssis in two parts:

1. Without replication
2. With replication

## To run everything on bc4

Use Snakemake to orchestrate the analysis. This will submit the jobs to slurm. 
So when you're in a screen session, to run the snakemake process:

```
module add languages/anaconda3/5.2.0-tflow-1.11
source ~/.bash_profile
snakemake -prk \
-j 100 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output}"
```

