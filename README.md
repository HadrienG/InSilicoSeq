# InSilicoSeq
## A sequencing simulator

[![Build Status](https://travis-ci.org/HadrienG/InSilicoSeq.svg?branch=master)](https://travis-ci.org/HadrienG/InSilicoSeq)
[![Documentation Status](https://readthedocs.org/projects/insilicoseq/badge/?version=0.8.0)](http://insilicoseq.readthedocs.io/en/0.8.0/?badge=0.8.0)
[![codecov](https://codecov.io/gh/HadrienG/InSilicoSeq/branch/master/graph/badge.svg)](https://codecov.io/gh/HadrienG/InSilicoSeq)
[![LICENSE](https://img.shields.io/badge/license-MIT-lightgrey.svg)](LICENSE)

InSilicoSeq (or iss) is a sequencing simulator producing (relatively) realistic
Illumina reads primarily intended for simulating metagenomic samples, although
it can be used to produce sequencing data from a single genome.

InSilicoSeq is written in python, and use a kernel density estimation model to
model the read quality of real sequencing data.

InSilicoSeq support substitution, insertion and deletion errors. If you don't
have the use for insertion and deletion error a basic error model is provided.

## Installation

To install InSilicoSeq, type the following in your terminal:

`pip install InSilicoSeq`

Alternatively, with docker:

```shell
docker pull hadrieng/insilicoseq:0.8.1
```

## Usage

InSilicoSeq has two modes: one to generate Illumina reads, the other to create
an error model from which the reads will take their characteristics.

We provide pre-computed error models that should be sufficient for most use
cases.

### Generate reads with a pre-computed error model

for generating 1 million reads modelling a MiSeq instrument:

```shell
iss generate --genomes genomes.fasta --model_file MiSeq --output MiSeq_reads
```

where `genomes.fasta` is a (multi-)fasta file containing the reference genome
from which the simulated reads will be generated.

Currently InSilicoSeq comes with 2 error models: `HiSeq2500` and `MiSeq`

If you have built your own model, give the `.npz` file to the `--model_file`
argument to simulate reads from your own error model.

For 10 million reads and a custom error model:

```shell
iss generate --genomes genomes.fasta -n 10000000 --model_file my_model.npz \
--output my_model_reads
```

For more examples and a full list of options, please refer to the full
[documentation](http://insilicoseq.readthedocs.io)

### Generate reads without input genomes

We can download some for you! InSilicoSeq can download random genomes from the
ncbi using the infamous [eutils](https://www.ncbi.nlm.nih.gov/books/NBK25501/)

```shell
iss generate --ncbi bacteria -n_genomes 10 --model_file MiSeq \
--output ncbi_reads
```

For full usage and other option, refer to the full
[documentation](http://insilicoseq.readthedocs.io)

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

We welcome contributions from the community! See our
[Contributing](CONTRIBUTING.md) guidelines

## Citation

A paper will be on its way. In the meantime if you use InSilicoSeq in your
research, please cite the poster

> Gourlé, Hadrien (2017): Simulating Illumina data with InSilicoSeq. figshare. https://doi.org/10.6084/m9.figshare.5053327.v1
