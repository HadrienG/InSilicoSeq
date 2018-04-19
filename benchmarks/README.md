# Benchmarks

## Introduction

20 random genomes, selected with

```bash
# only 10k reads because we only care about the genomes and mode basic (faster)
iss generate --ncbi bacteria viruses archaea -u 10 5 5 -n 10k --mode basic \
    --cpus 8 --output _reads
mv _reads_genomes.fasta test_genomes.fasta
```

## Software used in this benchmarks

### Simulators

See #49

### Other

* GNU time (1.9)

## Commands used and time info

### iss

```bash
env time -v iss generate --genomes test_genomes.fasta -n 0.5m --model miseq --cpus 4 --output iss/test_x
```

| run | time (m:s) | memory (kb) |
| --- | --- | --- |
| 1 | 8:19.88 | 335100 |
| 2 | 7:14.63 | 312708 |
| 3 | 8:38.22 | 439920 |
| 4 | 8:36.25 | 360900 |
| 5 | 9:05.43 | 700648 |
| 6 | 8:12.44 | 307652 |
| 7 | 7:44.63 | 394100 |
| 8 | 8:51.43 | 530440 |
| 9 | 8:31.71 | 305276 |
| 10 | 8:31.07 | 299068 |

### BEAR

```bash
env time -v bash run_bear.sh -g test_genomes.fasta -r reads -o bear/test_x
```
