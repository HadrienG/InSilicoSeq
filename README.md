# InSilicoSeq — Annotation-Aware Fork

This repository is an annotation-aware fork of [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq), a sequencing simulator for generating realistic Illumina reads.

The fork preserves the standard InSilicoSeq functionality and adds a separate annotation-aware workflow for generating reads from genomes with GTF/GFF coding-sequence annotations.

The main new command is:

```bash
iss generate-annotated
```

It can be used to generate paired-end reads while recording or restricting properties of R1 according to the source genome annotation.

## Main additions

The annotation-aware extension adds support for:

- loading CDS coordinates, strand, phase, and identifiers from GTF or GFF annotations;
- determining the CDS overlap of each generated R1 read;
- determining whether R1 is fully contained within a single CDS;
- calculating the reference-derived reading frame of R1;
- restricting generated reads to one of the six coding frames:
  `+1`, `+2`, `+3`, `-1`, `-2`, or `-3`;
- generating only reads fully contained within a CDS;
- optionally trimming simulated reads to a requested final read length;
- writing a separate annotation table containing genomic and CDS information for every generated R1 read;
- supporting chromosome, contig, and scaffold identifiers directly from FASTA/GTF/GFF files;
- efficiently sampling valid genomic positions when both a specific reading frame and full-CDS containment are requested.

The original `iss generate` workflow remains available and uses the standard InSilicoSeq generation path.

---

# Installation

## Clone this fork

```bash
git clone https://github.com/stavperez/InSilicoSeq.git
cd InSilicoSeq
```

Create an environment if desired:

```bash
conda create -n insilicoseq-annotated python=3.10
conda activate insilicoseq-annotated
```

Install the fork:

```bash
pip install -e .
```

After installation, verify that the new command is available:

```bash
iss generate-annotated --help
```

The standard InSilicoSeq commands remain available:

```bash
iss generate --help
iss model --help
```

---

# Quick start

The annotation-aware workflow can start from a CSV containing NCBI assembly accessions.

## 1. Create an organism CSV

Example:

```csv
organism_name,assembly_accession
Escherichia_coli_K12,GCF_000005845.2
Bacillus_subtilis_168,GCF_000009045.1
```

An example file is provided at:

```text
annotation_aware/example_organisms.csv
```

## 2. Download genomes and annotations

Run:

```bash
python annotation_aware/prepare_genome_csv.py \
    --input-csv annotation_aware/example_organisms.csv \
    --genomes-dir data/genomes \
    --output-csv data/organisms_with_paths.csv \
    --email your.email@example.com
```

This downloads the genomic FASTA and the corresponding GTF annotation when available, falling back to GFF when necessary.

The resulting CSV contains local paths such as:

```csv
organism_name,assembly_accession,fna_path,gtf_path,gff_path,download_source,download_status,download_note
Escherichia_coli_K12,GCF_000005845.2,/path/genomic.fna.gz,/path/genomic.gtf.gz,,RefSeq,ok,...
```

## 3. Generate annotation-aware reads

For example, generate 4,000 read pairs per organism using the MiSeq error model, trim the final reads to 150 bp, and require R1 to:

- originate from reading frame `+1`; and
- be fully contained within a single CDS.

```bash
iss generate-annotated \
    --csv data/organisms_with_paths.csv \
    --output-dir data/generated_reads \
    --read-pairs 4000 \
    --model miseq \
    --output-read-length 150 \
    --reading-frame 1 \
    --only-full-coding \
    --seed 42
```

---

# Input CSV

`iss generate-annotated` expects a CSV containing at least:

```text
organism_name
fna_path
```

When reading-frame filtering or full-CDS filtering is requested, the CSV must also contain a valid annotation path in either:

```text
gtf_path
```

or:

```text
gff_path
```

Example:

```csv
organism_name,fna_path,gtf_path,gff_path
Organism_A,/data/A_genomic.fna.gz,/data/A_genomic.gtf.gz,
Organism_B,/data/B_genomic.fna.gz,,/data/B_genomic.gff.gz
```

Both compressed (`.gz`) and uncompressed genome/annotation files can be used.

The FASTA sequence identifiers and GTF/GFF sequence identifiers must correspond.

For example, if a FASTA record has the identifier:

```text
NZ_CP012345.1
```

the corresponding CDS features in the GTF/GFF file must also use:

```text
NZ_CP012345.1
```

as their sequence identifier.

This is particularly important for draft assemblies containing multiple contigs or scaffolds.

---

# `iss generate-annotated`

## Basic usage

```bash
iss generate-annotated \
    --csv organisms.csv \
    --output-dir reads
```

By default, the command generates:

- 4,000 read pairs per organism;
- using the InSilicoSeq MiSeq KDE error model;
- with GC-bias filtering enabled;
- without reading-frame restriction;
- without requiring full-CDS containment;
- at the native read length of the selected error model.

## Number of read pairs

```bash
--read-pairs 4000
```

The value specifies the requested number of paired-end read pairs **per organism**.

## Error model

The default error model is:

```bash
--model miseq
```

A different InSilicoSeq model can be selected, for example:

```bash
--model hiseq
```

or:

```bash
--model novaseq
```

The error-generation mode can be selected with:

```bash
--mode kde
```

Valid modes are:

```text
kde
basic
perfect
```

For `kde`, `--model` can also point to a custom InSilicoSeq error-model file.

## Random seed

For reproducible simulation:

```bash
--seed 42
```

## Final read length

By default, generated reads retain the read length of the selected InSilicoSeq error model.

A shorter final read length can be requested with:

```bash
--output-read-length 150
```

Trimming is applied after error simulation.

The amount removed is distributed as evenly as possible between the two ends of each read, and the corresponding quality scores are retained.

For example, a 300-bp simulated read with:

```bash
--output-read-length 150
```

is trimmed by 75 bases from each end.

The requested output length cannot be greater than the native read length of the selected error model.

---

# CDS-aware generation

## Reading-frame filtering

Use:

```bash
--reading-frame FRAME
```

Valid values are:

```text
-3
-2
-1
1
2
3
```

For example:

```bash
--reading-frame 1
```

retains R1 reads assigned to reference-derived frame `+1`.

Reading-frame filtering requires a GTF/GFF annotation.

## Full-CDS filtering

Use:

```bash
--only-full-coding
```

to require the final R1 genomic interval to be completely contained within one CDS.

This is deliberately stricter than requiring 100% aggregate overlap with coding sequence: the complete read must fit within a single CDS feature.

A GTF/GFF annotation is required.

## Combining reading-frame and CDS filtering

For example:

```bash
iss generate-annotated \
    --csv organisms.csv \
    --output-dir reads \
    --reading-frame -2 \
    --only-full-coding
```

When both a reading frame and full-CDS containment are requested, valid genomic start positions are calculated directly and sampled from those positions.

This avoids relying entirely on repeated rejection sampling for this common use case.

Other annotation-aware modes use rejection sampling.

---

# Reading-frame definition

Reading frame is calculated from the genomic source of R1 using:

- the genomic coordinates of the read;
- the strand of the best-overlapping CDS;
- the CDS phase in the GTF/GFF annotation.

Possible coding frames are:

```text
+1
+2
+3
-1
-2
-3
```

Positive frames correspond to CDSs on the forward strand.

Negative frames correspond to CDSs on the reverse strand.

A value of:

```text
0
```

means that no CDS-derived reading frame was assigned.

The reported reading frame is **reference-derived**.

It represents the coding frame from which the read originated in the reference genome, rather than a frame inferred by translating the final error-containing read sequence.

When `--output-read-length` is used, CDS overlap, CDS containment, and reading frame are calculated using the final trimmed R1 genomic interval.

---

# GC-bias filtering

GC-bias filtering is enabled by default for `iss generate-annotated`.

Disable it with:

```bash
--no-gc-bias
```

It can be explicitly enabled with:

```bash
--gc-bias
```

GC-bias filtering preferentially retains reads with GC fractions between 0.40 and 0.60 while retaining other reads with 90% probability.

---

# Fragment-length parameters

Custom fragment-length parameters can be supplied with:

```bash
--fragment-length LENGTH
--fragment-length-sd SD
```

Both parameters must be supplied together.

If they are omitted, the fragment-length behavior of the selected InSilicoSeq error model is used.

---

# Existing output files

If all expected output files for an organism already exist, that organism is skipped.

Use:

```bash
--force
```

to regenerate and replace existing outputs.

---

# Output files

For each organism, the command generates paired FASTQ files, a mutation table, and an R1 annotation table.

When frame `+1` is requested:

```text
<organism>_frame_p1_R1.fastq
<organism>_frame_p1_R2.fastq
<organism>_frame_p1_mutations.tsv
<organism>_frame_p1_R1_annotations.tsv
```

For negative frames, `m` is used.

For example, frame `-2` produces filenames containing:

```text
frame_m2
```

When no reading-frame restriction is requested, the `_frame_*` component is omitted.

The output directory also contains:

```text
generation_status.csv
```

---

# R1 annotation table

The annotation table contains:

```text
read_id
chromosome
start
end
is_cds
fully_contained_in_cds
overlap_pct
coding_bases
read_length
reading_frame
cds_strand
cds_phase
best_cds_id
```

## Column definitions

### `read_id`

Identifier of the generated R1 read.

### `chromosome`

FASTA sequence, chromosome, contig, or scaffold identifier used to match the genome sequence to the annotation.

### `start`

Zero-based genomic start coordinate of the final R1 interval.

### `end`

Zero-based exclusive genomic end coordinate of the final R1 interval.

### `is_cds`

Whether R1 overlaps at least one CDS.

### `fully_contained_in_cds`

Whether the complete R1 interval is contained within a single CDS.

### `overlap_pct`

Percentage of R1 bases overlapping the union of CDS intervals.

### `coding_bases`

Number of R1 bases overlapping CDS sequence.

### `read_length`

Length of the genomic R1 interval used for annotation.

### `reading_frame`

Reference-derived reading frame:

```text
-3
-2
-1
0
1
2
3
```

`0` means that no CDS-derived frame was assigned.

### `cds_strand`

Strand of the best-overlapping CDS.

### `cds_phase`

CDS phase reported by the source GTF/GFF annotation.

### `best_cds_id`

Identifier of the CDS with the greatest overlap with R1.

---

# Generation status

`generation_status.csv` contains:

```text
organism_name
requested_read_pairs
generated_read_pairs
status
```

Possible successful statuses include:

```text
ok
ok_partial
skipped_existing
```

`ok_partial` means that generation stopped before the requested number of read pairs could be produced.

This can occur when the requested constraints leave too few acceptable genomic positions or when the simulation reaches its attempt limit.

---

# Preparing genomes from NCBI

The helper:

```text
annotation_aware/prepare_genome_csv.py
```

accepts a CSV containing:

```text
organism_name
assembly_accession
```

and resolves each assembly through NCBI.

For each assembly it:

1. resolves the NCBI Assembly accession;
2. prefers the RefSeq assembly path when available;
3. otherwise uses the GenBank assembly path;
4. downloads the genomic FASTA;
5. attempts to download the genomic GTF;
6. falls back to genomic GFF when GTF is unavailable;
7. records the downloaded paths and status in an output CSV.

Example:

```bash
python annotation_aware/prepare_genome_csv.py \
    --input-csv annotation_aware/example_organisms.csv \
    --genomes-dir data/genomes \
    --output-csv data/organisms_with_paths.csv \
    --email your.email@example.com
```

The email address is required for NCBI Entrez access.

---

# Generating all six frames

The six reading frames can be generated independently:

```bash
for FRAME in -3 -2 -1 1 2 3; do
    iss generate-annotated \
        --csv data/organisms_with_paths.csv \
        --output-dir data/generated_reads \
        --read-pairs 4000 \
        --model miseq \
        --output-read-length 150 \
        --reading-frame "${FRAME}" \
        --only-full-coding \
        --seed 42
done
```

Each frame receives a distinct filename suffix, so the outputs can coexist in the same directory.

---

# Standard InSilicoSeq functionality

This fork retains the standard InSilicoSeq commands.

For example:

```bash
iss generate \
    --genomes genomes.fasta \
    --model miseq \
    --output simulated_reads
```

and:

```bash
iss model \
    --bam aligned_reads.bam \
    --output custom_model
```

The annotation-aware functionality is intentionally provided through a separate command:

```bash
iss generate-annotated
```

so that the standard InSilicoSeq generation path remains distinct from annotation-aware simulation.

For documentation of the standard InSilicoSeq functionality, see:

- https://github.com/HadrienG/InSilicoSeq
- https://insilicoseq.readthedocs.io/

---

# Differences from upstream InSilicoSeq

The annotation-aware fork adds a separate generation path that can use GTF/GFF CDS annotations.

The main additional capabilities are:

| Capability                          | Standard `iss generate` | `iss generate-annotated` |
| ----------------------------------- | ----------------------: | -----------------------: |
| Standard InSilicoSeq error models   |                     Yes |                      Yes |
| Paired-end simulation               |                     Yes |                      Yes |
| CDS overlap annotation              |                      No |                      Yes |
| Full-CDS filtering                  |                      No |                      Yes |
| Reference-derived reading frame     |                      No |                      Yes |
| Reading-frame filtering             |                      No |                      Yes |
| Per-read CDS annotation table       |                      No |                      Yes |
| Optional final read-length trimming |                      No |                      Yes |
| CSV-based multi-organism workflow   |                      No |                      Yes |

The original generation functions remain available and are not replaced by the annotation-aware API.

---

# Important interpretation notes

## Annotation is based on genomic origin

CDS properties are calculated from the reference genomic coordinates from which a read was simulated.

They are not inferred from alignment or translation of the final simulated read.

## Sequencing errors do not redefine the source frame

Substitutions and indels introduced by the error model may alter the final read sequence.

The reported reading frame nevertheless describes the frame of the source CDS in the reference genome.

## R1 defines annotation-aware filtering

Reading-frame and CDS-containment filtering are applied to R1.

R2 belongs to the same simulated read pair but is not independently required to satisfy the R1 CDS/frame constraint.

---

# Requirements

The annotation-aware extension adds the following dependencies to the standard InSilicoSeq requirements:

```text
pandas
gffutils
```

The NCBI preparation helper additionally uses dependencies already used or available in the project environment, including:

```text
requests
biopython
```

---

# Upstream project

This repository is a fork of:

[HadrienG/InSilicoSeq](https://github.com/HadrienG/InSilicoSeq)

The original InSilicoSeq software and standard sequencing-simulation functionality were developed by Hadrien Gourlé and contributors.

This fork adds the annotation-aware generation functionality described above.

---

# Citation

If you use InSilicoSeq or this fork, please cite the original InSilicoSeq publication:

> Gourlé H, Karlsson-Lindsjö O, Hayer J, Bongcam-Rudloff E.  
> Simulating Illumina data with InSilicoSeq.  
> _Bioinformatics_. 2019;35(3):521–522.  
> https://doi.org/10.1093/bioinformatics/bty630

If a publication describing the annotation-aware extension becomes available, its citation can be added here.

---

# License

InSilicoSeq is distributed under the MIT License.

This fork retains the upstream license and copyright notices.

See:

```text
LICENSE
```

for the full license text.

---

# Issues

For problems specifically related to the annotation-aware functionality in this fork, open an issue at:

https://github.com/stavperez/InSilicoSeq/issues

For questions or issues involving the original InSilicoSeq implementation, consult the upstream project:

https://github.com/HadrienG/InSilicoSeq
