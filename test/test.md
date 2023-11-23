## log
### 2023-11-21 13:47:54
```bash
conda activate metabolicv4

perl METABOLIC-G.pl -t $SLURM_NTASKS -in test -kofam-db full -o test/out
perl METABOLIC-G-test.pl -t $SLURM_NTASKS -in test -kofam-db small -o test/out_compare
```
