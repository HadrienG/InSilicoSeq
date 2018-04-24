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

*Failed to run*


### FASTQSim

```bash

```

*Failed to run*

### GemSim

#### Error model generation

```bash
env time -v python GemErr.py -f ../../MiSeq_300bp/final.contigs.fa -s ../../MiSeq_300bp/45_S1_L001_R1_001.fastq.bam -n test -i 21 -r 300
```

```bash
env time -v
```

*Both failed to run*

### Grinder

```bash
env time -v grinder -reference_file ../test_genomes.fasta -abundance_model powerlaw 0.1 -total_reads 500000 -rd 300 -md poly4 3e-3 3.3e-8 -ql 30 10 -fq 1 -insert_dist 350 normal 20 -bn test_1
```

| run | time (h:m:s) | memory (kb) |
| --- | --- | --- |
| 1 | 13:39:23 | 164332 |
| 2 | 14:27:16 | 159760 |
| 3 | 13:11:57 | 164344 |
| 4 | 13:12:05 | 164332 |
| 5 |

### pIRS

```bash
env time -v pirs simulate -x 18 -l 300 -m 800 -t 4 ../test_genomes.fasta -B miseq_profile.count.matrix -I miseq_indel.InDel.matrix -G /opt/sw/pirs/2.0.2/Profiles/GC-depth_Profiles/humNew.gcdep_200.dat -s test_x
```

| run | time (m:s) | memory (kb) |
| --- | --- | --- |
| 1 | 0:36.08 | 36912 |
| 2 | 0:35.72 | 37040 |
| 3 | 0:36.06 | 35192 |
| 4 | 0:36.36 | 35192 |
| 5 | 0:36.43 | 35160 |
| 6 | 0:36.53 | 35276 |
| 7 | 0:35.68 | 36416 |
| 8 | 0:36.03 | 35208 |
| 9 | 0:36.70 | 36200 |
| 10 | 0:36.06 | 36724 |
