# InSilicoSeq
## A sequencing simulator

[![Build Status](https://travis-ci.org/HadrienG/InSilicoSeq.svg?branch=master)](https://travis-ci.org/HadrienG/InSilicoSeq)
[![codecov](https://codecov.io/gh/HadrienG/InSilicoSeq/branch/master/graph/badge.svg)](https://codecov.io/gh/HadrienG/InSilicoSeq)
[![LICENSE](https://img.shields.io/badge/license-MIT-lightgrey.svg)](LICENSE)

InSilicoSeq (or iss) is a sequencing simulator. iss produces
(relatively) realistic Illumina reads, and was primarily intended for
simulating metagenomic samples, although it can be used to produce sequencing
data from a single genome.

InSilicoSeq is written in python, and use a kernel density estimation model to
model the read quality of real sequencing data. InSilicoSeq support
substitution, insertion and deletion errors. If you don't have the use for
insertion and deletion error a basic error model is provided.

## Installation

To install InSilicoSeq, simply type the following in your terminal:

`pip install git+https://github.com/HadrienG/InSilicoSeq.git`

Alternatively, with docker:

```shell
docker pull hadrieng/insilicoseq:0.1.0
docker run -it --rm hadrieng/insilicoseq:0.1.0 iss
```

## Usage

InSilicoSeq has two modes: one to generate Illumina reads, the other to create
an error model from which the reads will take their characteristics.

We provide pre-computed error models that should be sufficient for most use
cases.

### Generate Reads

```
usage: iss generate [-h] [--quiet] [--debug] --genomes <fasta>
                    [--abundance <txt>] [--n_reads <int>]
                    [--model ['cdf', 'kde', 'basic']] [--model_file <npz>]
                    --output <fastq>

simulate reads from an error model

arguments:
  -h, --help            show this help message and exit
  --quiet, -q           Disable info logging. (default: False).
  --debug, -d           Enable debug logging. (default: False).
  --genomes <fasta>, -g <fasta>
                        Input genome(s) from where the reads will originate
                        (Required)
  --abundance <txt>, -a <txt>
                        abundance file for coverage calculations (default:
                        None)
  --n_reads <int>, -n <int>
                        Number of reads to generate (default: 1000000)
  --model ['cdf', 'kde', 'basic'], -m ['cdf', 'kde', 'basic']
                        Error model. If not specified, using kernel density
                        estimation (default: kde). Can be 'kde', 'cdf' or
                        'basic'
  --model_file <npz>, -f <npz>
                        Error model file. If not specified, using a basic
                        error model instead (default: None)
  --output <fastq>, -o <fastq>
                        Output file prefix (Required)
```

### Create an Error Model

*TODO*

## License

Code is under the [MIT](LICENSE) license.

## Issues

Found a bug or have a question? Please open an [issue](https://github.com/HadrienG/InSilicoSeq/issues)

## Contributing

*TODO*
