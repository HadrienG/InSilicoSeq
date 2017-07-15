# InSilicoSeq
## A sequencing simulator

[![Build Status](https://travis-ci.org/HadrienG/InSilicoSeq.svg?branch=master)](https://travis-ci.org/HadrienG/InSilicoSeq)
[![codecov](https://codecov.io/gh/HadrienG/InSilicoSeq/branch/master/graph/badge.svg)](https://codecov.io/gh/HadrienG/InSilicoSeq)
[![LICENSE](https://img.shields.io/badge/license-MIT-lightgrey.svg)](LICENSE)

InSilicoSeq (or iss) is a sequencing simulator producing (relatively) realistic
Illumina reads primarily intended for simulating metagenomic samples, although it can be used to produce sequencing data from a single genome.

InSilicoSeq is written in python, and use a kernel density estimation model to
model the read quality of real sequencing data.

InSilicoSeq support substitution, insertion and deletion errors. If you don't have
the use for insertion and deletion error a basic error model is provided.

## Installation

To install InSilicoSeq, simply type the following in your terminal:

`pip install git+https://github.com/HadrienG/InSilicoSeq.git@0.2.0`

Alternatively, with docker:

```shell
docker pull hadrieng/insilicoseq:0.2.0
docker run -it --rm hadrieng/insilicoseq:0.2.0 iss
```

## Usage

InSilicoSeq has two modes: one to generate Illumina reads, the other to create
an error model from which the reads will take their characteristics.

We provide pre-computed error models that should be sufficient for most use
cases.

### Generate reads with a pre-computed error model

for generating 1M reads using the MiSeq300bp error model:

```shell
iss generate --genomes genomes.fasta --abundance abundance_file.txt \
    --model_file MiSeq300bp --output HiSeq_reads
```

where `genomes.fasta` is a (multi-)fasta file containing the reference genome from which the simulated reads will be generated, and `abundance_file.txt` a tab-delimited file containing abundance information.

Currently InSilicoSeq comes with 3 error models: `HiSeq2500`, `MiSeq250bp` and
`MiSeq300bp`

### Example of genomes and abundance file

```
# multi-fasta file
>genome_A
ATGC...
>genome_B
CCGT...
...

# abundance file (total abundance must be 1!)
genome_A    0.2
genome_B    0.4
...
```

### Create your own error model

If you do not wish to use the pre-computed error models provided with
InSilicoSeq, it is possible to create your own.

Align you reads against the reference:

```shell
    bowtie2-build genomes.fasta genomes
    bowtie2 -x genomes -1 reads_R1.fastq.gz \
        -2 reads_R2.fastq.gz | samtools view -bS | samtools sort -o genomes.bam
    samtools index genomes.bam
```

then build the model:

```shell
    iss model -b genomes.bam -o genomes
```

which will create a `genome.npz` file containing your newly built model

## License

Code is under the [MIT](LICENSE) license.

## Issues

Found a bug or have a question? Please open an [issue](https://github.com/HadrienG/InSilicoSeq/issues)

## Contributing

*TODO: Code of conduct and instruction for Contributing coming soon!*

## Citation

A paper will be on its way. In the meantime if you use InSilicoSeq in your
research, please cite the poster

> Gourlé, Hadrien (2017): Simulating Illumina data with InSilicoSeq. figshare. https://doi.org/10.6084/m9.figshare.5053327.v1
