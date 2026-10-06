#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gc
import logging
import os
import random
import sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
import gffutils
from iss import abundance, download, util
from iss.error_models import basic, kde, perfect
from iss.util import load, rev_comp


# ============================================================
# Original InSilicoSeq generator
# ============================================================


def simulate_reads(
    record,
    error_model,
    n_pairs,
    cpu_number,
    forward_handle,
    reverse_handle,
    mutations_handle,
    sequence_type,
    gc_bias=False,
    mode="default",
):
    """Simulate reads from one genome (or sequence) according to an ErrorModel

    This function makes use of the `simulate_read` function to simulate reads
    and save them in a fastq file

    Args:
        record (SeqRecord): sequence or genome of reference
        error_model (ErrorModel): an ErrorModel
        n_pairs (int): the number of reads to generate
        cpu_number (int): an int indentifying the cpu that is used by the
            function. Is used for naming the output file
        forward_handle (file): a file handle to write the forward reads to
        reverse_handle (file): a file handle to write the reverse reads to
        mutations_handle (file): a file handle to write the mutations to
        sequencing_type (str): metagenomics or amplicon sequencing used
        gc_bias (bool): if set, the function may skip a read due to abnormal
            GC content

    """
    logger = logging.getLogger(__name__)

    # load the record from disk if mode is memmap
    if mode == "memmap":
        record_mmap = load(record)
        record = record_mmap

    logger.debug("Cpu #%s: Generating %s read pairs" % (cpu_number, n_pairs))

    for forward_record, reverse_record, mutations in reads_generator(
        n_pairs, record, error_model, cpu_number, gc_bias, sequence_type
    ):
        SeqIO.write(forward_record, forward_handle, "fastq-sanger")
        SeqIO.write(reverse_record, reverse_handle, "fastq-sanger")
        write_mutations(mutations, mutations_handle)


def reads_generator(n_pairs, record, error_model, cpu_number, gc_bias, sequence_type):
    logger = logging.getLogger(__name__)

    i = 0
    while i < n_pairs:
        try:
            forward, reverse, mutations = simulate_read(record, error_model, i, cpu_number, sequence_type)

        except AssertionError:
            logger.warning("%s shorter than read length for this ErrorModel" % record.id)
            logger.warning("Skipping %s. You will have less reads than specified" % record.id)
            break
        else:
            if gc_bias:
                stiched_seq = forward.seq + reverse.seq
                gc_content = gc_fraction(stiched_seq)
                if 40 < gc_content < 60:
                    yield (forward, reverse, mutations)
                    i += 1
                elif np.random.rand() < 0.90:
                    yield (forward, reverse, mutations)
                    i += 1
                else:
                    continue
            else:
                yield (forward, reverse, mutations)
                i += 1


def simulate_read(record, error_model, i, cpu_number, sequence_type):
    """From a read pair from one genome (or sequence) according to an
    ErrorModel

    Each read is a SeqRecord object
    returns a tuple containing the forward and reverse read.

    Args:
        record (SeqRecord): sequence or genome of reference
        error_model (ErrorModel): an ErrorModel class
        i (int): a number identifying the read
        cpu_number (int): cpu number. Is added to the read id.
        sequence_type (str): metagenomics or amplicon sequencing used

    Returns:
        tuple: tuple containg a forward read and a reverse read
    """
    logger = logging.getLogger(__name__)
    sequence = record.seq
    header = record.id

    read_length = error_model.read_length

    if error_model.fragment_length is not None and error_model.fragment_sd is not None:
        fragment_length = int(np.random.normal(error_model.fragment_length, error_model.fragment_sd))
        insert_size = fragment_length - (read_length * 2)
    else:
        insert_size = error_model.random_insert_size()
        fragment_length = insert_size + (read_length * 2)

    # generate the forward read
    try:  # a ref sequence has to be longer than 2 * read_length + i_size
        assert read_length < len(record.seq)
        # assign the start position of the forward read
        # if sequence_type == metagenomics, get a random start position
        # if sequence_type == amplicon, start position is the start of the read
        if sequence_type == "metagenomics":
            forward_start = random.randrange(0, len(record.seq) - fragment_length)
        elif sequence_type == "amplicon":
            forward_start = 0
        else:
            raise RuntimeError(f"sequence type '{sequence_type}' is not supported")
    except AssertionError:
        raise
    except ValueError as e:
        logger.debug("%s shorter than template length for this ErrorModel:%s" % (record.id, e))
        forward_start = max(0, random.randrange(0, len(record.seq) - read_length))

    forward_end = forward_start + read_length
    bounds = (forward_start, forward_end)
    # create a perfect read
    forward = SeqRecord(
        Seq(str(sequence[forward_start:forward_end])), id="%s_%s_%s/1" % (header, i, cpu_number), description=""
    )
    forward.annotations["mutations"] = []
    forward.annotations["original"] = str(forward.seq)

    # add the indels, the qual scores and modify the record accordingly
    forward = error_model.introduce_indels(forward, "forward", sequence, bounds)
    forward = error_model.introduce_error_scores(forward, "forward")
    forward = error_model.mut_sequence(forward, "forward")

    # generate the reverse read
    # assign start position reverse read
    # if sequence_type == metagenomics, get a start position based on insert_size
    # if sequence_type == amplicon, start position is the end of the read
    if sequence_type == "metagenomics":
        reverse_start = forward_end + insert_size
        reverse_end = reverse_start + read_length
    elif sequence_type == "amplicon":
        reverse_start = len(record.seq) - read_length
        reverse_end = reverse_start + read_length
    else:
        raise ValueError(f"Sequence type {sequence_type} not known")
    if reverse_end > len(record.seq):
        # we use random insert when the modelled template length distribution
        # is too large
        reverse_end = random.randrange(read_length, len(record.seq))
        reverse_start = reverse_end - read_length
    bounds = (reverse_start, reverse_end)
    # create a perfect read
    reverse = SeqRecord(
        Seq(rev_comp(str(sequence[reverse_start:reverse_end]))),
        id="%s_%s_%s/2" % (header, i, cpu_number),
        description="",
    )
    reverse.annotations["mutations"] = []
    reverse.annotations["original"] = str(reverse.seq)

    # add the indels, the qual scores and modify the record accordingly
    reverse = error_model.introduce_indels(reverse, "reverse", sequence, bounds)
    reverse = error_model.introduce_error_scores(reverse, "reverse")
    reverse = error_model.mut_sequence(reverse, "reverse")

    return (forward, reverse, forward.annotations["mutations"] + reverse.annotations["mutations"])


def to_fastq(generator, output):
    """Write reads to a fastq file

    Take a generator or a list containing read pairs (tuples) and write them
        in two fastq files: output_R1.fastq and output_R2.fastq

    Args:
        generator (generator): a read generator (or list)
        output (string): the output files prefix
    """
    logger = logging.getLogger(__name__)
    # define name of output files
    output_forward = output + "_R1.fastq"
    output_reverse = output + "_R2.fastq"

    try:
        f = open(output_forward, "a")
        r = open(output_reverse, "a")
    except PermissionError as e:
        logger.error("Failed to open output file(s): %s" % e)
        sys.exit(1)
    else:
        with f, r:
            for read_tuple in generator:
                SeqIO.write(read_tuple[0], f, "fastq-sanger")
                SeqIO.write(read_tuple[1], r, "fastq-sanger")


def worker_iterator(work, error_model, cpu_number, worker_prefix, seed, sequence_type, gc_bias):
    """A utility function to run the reads simulation of each record in a loop for a specific cpu"""
    logger = logging.getLogger(__name__)
    try:
        forward_handle = open(f"{worker_prefix}_R1.fastq", "w")
        reverse_handle = open(f"{worker_prefix}_R2.fastq", "w")
        mutation_handle = open(f"{worker_prefix}.vcf", "w")
    except PermissionError as e:
        logger.error("Failed to write temporary output file(s): %s" % e)
        sys.exit(1)

    if seed is not None:
        random.seed(seed + cpu_number)
        np.random.seed(seed + cpu_number)

    with forward_handle, reverse_handle, mutation_handle:
        for record, n_pairs, mode in work:
            simulate_reads(
                record=record,
                n_pairs=n_pairs,
                mode=mode,
                error_model=error_model,
                cpu_number=cpu_number,
                forward_handle=forward_handle,
                reverse_handle=reverse_handle,
                mutations_handle=mutation_handle,
                sequence_type=sequence_type,
                gc_bias=gc_bias,
            )


def generate_work_divider(
    fasta_file, readcount_dic, abundance_dic, n_reads, coverage, coverage_file, error_model, output, chunk_size
):
    """Yields a list of tuples containing the records and the number of reads to generate for each record

    Yields:
        list[tuple[SeqIO.Record, int, str]]: a list of tuples containing the records,
            the number of reads to generate for each record and the mode
    """
    logger = logging.getLogger(__name__)

    current_chunk = 0
    total_reads_generated = 0
    total_reads_generated_unrounded = 0

    chunk_work = []

    for record in fasta_file:
        # generate reads for records
        if readcount_dic is not None:
            if record.id not in readcount_dic:
                logger.warning(f"Record {record.id} not found in readcount file")
                continue
            n_pairs_unrounded = readcount_dic[record.id] / 2
        elif abundance_dic is not None:
            if record.id not in abundance_dic:
                logger.warning(f"Record {record.id} not found in abundance file")
                continue
            record_abundance = abundance_dic[record.id]

            genome_size = len(record.seq)

            if coverage or coverage_file:
                record_coverage = record_abundance
            else:
                record_coverage = abundance.to_coverage(
                    n_reads,
                    record_abundance,
                    error_model.read_length,
                    genome_size,
                )
            n_pairs_unrounded = ((record_coverage * len(record.seq)) / error_model.read_length) / 2
        else:
            raise RuntimeError("No readcount or abundance file provided")

        # check that the rounding does not cause to drop read pairs
        n_pairs = round(n_pairs_unrounded)
        total_reads_generated_unrounded += n_pairs_unrounded
        total_reads_generated += n_pairs
        if round(total_reads_generated_unrounded) > total_reads_generated:
            logger.debug("Adding a pair to correct rounding error")
            n_pairs += 1
            total_reads_generated += 1

        logger.debug("Will generate %s read pairs for %s" % (n_pairs, record.id))

        if n_pairs == 0:
            continue

        # due to a bug in multiprocessing
        # https://bugs.python.org/issue17560
        # we can't send records taking more than 2**31 bytes
        # through serialisation.
        # In those cases we use memmapping
        if sys.getsizeof(str(record.seq)) >= 2**31 - 1:
            logger.warning("record %s unusually big." % record.id)
            logger.warning("Using a memory map.")
            mode = "memmap"

            record_mmap = "%s.memmap" % output
            if os.path.exists(record_mmap):
                os.unlink(record_mmap)
            util.dump(record, record_mmap)
            del record
            record = record_mmap
            gc.collect()
        else:
            mode = "default"

        n_pairs_remaining = n_pairs
        while n_pairs_remaining > 0:
            chunk_remaining = chunk_size - current_chunk

            if n_pairs_remaining <= chunk_remaining:
                # Record fits in the current chunk
                chunk_work.append((record, n_pairs_remaining, mode))
                n_pairs_added = n_pairs_remaining
            else:
                # Record does not fit in the current chunk
                chunk_work.append((record, chunk_remaining, mode))
                n_pairs_added = chunk_remaining

            n_pairs_remaining -= n_pairs_added
            current_chunk += n_pairs_added

            if current_chunk == chunk_size:
                yield chunk_work
                chunk_work = []
                current_chunk = 0

    if chunk_work:
        # Yield the last (not full) chunk
        yield chunk_work


def load_error_model(mode, seed, model, fragment_length, fragment_length_sd, store_mutations):
    """
    Load the error model based on the specified mode and parameters.

    Args:
        mode (str): The mode of the error model. Possible values are 'kde', 'basic', and 'perfect'.
        seed (int): The random seed to use for generating random numbers.
        model (str): The model to use for the error model. Only applicable for the 'kde' mode.
        fragment_length (float): The mean fragment length for the error model.
        fragment_length_sd (float): The standard deviation of the fragment length for the error model.

    Returns:
        err_mod: The loaded error model based on the specified mode and parameters.
    """
    logger = logging.getLogger(__name__)

    logger.info("Using %s ErrorModel" % mode)

    precomputed_error_models = {
        "hiseq": os.path.join(os.path.dirname(__file__), "profiles/HiSeq"),
        "novaseq": os.path.join(os.path.dirname(__file__), "profiles/NovaSeq"),
        "miseq": os.path.join(os.path.dirname(__file__), "profiles/miSeq_0.npz"),
        "miseq-20": os.path.join(os.path.dirname(__file__), "profiles/miSeq_20.npz"),
        "miseq-24": os.path.join(os.path.dirname(__file__), "profiles/miSeq_24.npz"),
        "miseq-28": os.path.join(os.path.dirname(__file__), "profiles/miSeq_28.npz"),
        "miseq-32": os.path.join(os.path.dirname(__file__), "profiles/miSeq_32.npz"),
        "miseq-36": os.path.join(os.path.dirname(__file__), "profiles/miSeq_36.npz"),
        "nextseq": os.path.join(os.path.dirname(__file__), "profiles/nextSeq.npz"),
    }

    if fragment_length is not None and fragment_length_sd is not None:
        logger.info(
            f"Using custom fragment length {fragment_length} and default fragment length sd {fragment_length_sd}"
        )
    elif bool(fragment_length) ^ bool(fragment_length_sd):
        logger.error("fragment_length and fragment_length_sd must be specified together")
        sys.exit(1)

    if seed:
        logger.info("Setting random seed to %i" % seed)
        random.seed(seed)
        np.random.seed(seed)
    if mode == "kde":
        if model is None:
            logger.error("--model is required in --mode kde")
            sys.exit(1)
        elif model.lower() in precomputed_error_models:
            npz = precomputed_error_models[model.lower()]
        else:
            npz = model
        err_mod = kde.KDErrorModel(npz, fragment_length, fragment_length_sd, store_mutations)
    elif mode == "basic":
        if model is not None:
            logger.warning("--model %s will be ignored in --mode %s" % (model, mode))

        err_mod = basic.BasicErrorModel(fragment_length, fragment_length_sd, store_mutations)
    elif mode == "perfect":
        if model is not None:
            logger.warning("--model %s will be ignored in --mode %s" % (model, mode))

        err_mod = perfect.PerfectErrorModel(fragment_length, fragment_length_sd)

    return err_mod


def load_genomes(genomes, draft, ncbi, n_genomes_ncbi, output, n_genomes):
    """Load genomes from different sources and concatenate them into a single file.

    Args:
        genomes (list): List of paths to genome files.
        draft (list): List of paths to draft genome files.
        ncbi (bool): Flag indicating whether to download genomes from NCBI.
        n_genomes_ncbi (list): List of tuples specifying the NCBI genome IDs and the number of genomes to download.
        output (str): Path to the output file.
        n_genomes (int): Number of genomes to select randomly from the concatenated genome file.

    Returns:
        tuple: A tuple containing the number of genomes and the path to the concatenated genome file.
    """

    logger = logging.getLogger(__name__)
    if not (genomes or draft or ncbi):
        logger.error("One of --genomes/-g, --draft, --ncbi/-k is required")
        sys.exit(1)

    genome_files = []
    if genomes:
        genome_files.extend(genomes)

    if draft:
        genome_files.extend(draft)

    if ncbi and n_genomes_ncbi:
        util.genome_file_exists(output + "_ncbi_genomes.fasta")
        if len(*ncbi) != len(*n_genomes_ncbi):
            logger.error(
                "--ncbi and --n_genomes_ncbi of unequal lengths. \
                Aborting"
            )
            sys.exit(1)

        for g, n in zip(*ncbi, *n_genomes_ncbi):
            genomes_ncbi = download.ncbi(g, n, output + "_ncbi_genomes.fasta")
        genome_files.append(genomes_ncbi)

    if ncbi and not n_genomes_ncbi:
        logger.error("--ncbi/-k requires --n_genomes_ncbi/-U. Aborting.")
        sys.exit(1)

    genome_file = output + ".iss.tmp.genomes.fasta"
    util.concatenate(genome_files, output=genome_file)

    # for n_genomes we use reservoir sampling to draw random genomes
    # from the concatenated genome file. We then override the file.
    if n_genomes and not draft and not ncbi:
        genome_count = util.count_records(genome_file)
        genome_files = [genome for genome in util.reservoir(SeqIO.parse(genome_file, "fasta"), genome_count, n_genomes)]
        SeqIO.write(genome_files, genome_file, "fasta")

    if os.stat(genome_file).st_size == 0:
        logger.error("Genome(s) file seems empty: %s" % genome_file)
        sys.exit(1)

    try:
        f = open(genome_file, "r")
        with f:  # count the number of records
            genome_list = util.count_records(f)
    except IOError as e:
        logger.error("Failed to open genome(s) file:%s" % e)
        sys.exit(1)

    return genome_list, genome_file


def load_readcount_or_abundance(
    readcount_file,
    abundance_file,
    coverage_file,
    coverage,
    abundance_distribution,
    draft,
    genome_list,
    genome_file,
    n_reads,
    output,
    error_model,
):
    """Load readcount or abundance information based on the provided input parameters.

    Args:
        readcount_file (str): Path to the readcount file.
        abundance_file (str): Path to the abundance file.
        coverage_file (str): Path to the coverage file.
        coverage (str): Coverage distribution type.
        abundance_distribution (str): Abundance distribution type.
        draft (str): Draft mode.
        genome_list (list): List of genomes.
        genome_file (str): Path to the genome file.
        n_reads (int): Number of reads.
        output (str): Output file path.
        error_model (object): Error model object.

    Returns:
        tuple: A tuple containing the readcount dictionary and abundance dictionary.
    """
    logger = logging.getLogger(__name__)
    abundance_dispatch = {
        "uniform": abundance.uniform,
        "halfnormal": abundance.halfnormal,
        "exponential": abundance.exponential,
        "lognormal": abundance.lognormal,
        "zero_inflated_lognormal": abundance.zero_inflated_lognormal,
    }
    readcount_dic = None
    abundance_dic = None
    if readcount_file:
        logger.info("Using readcount file:%s" % readcount_file)
        logger.warning("--readcount_file disables --n_reads, n_reads will be calculated from the readcount file")
        if draft:
            raise RuntimeError("readcount_file is only supported using --genomes, not --draft")
        readcount_dic = abundance.parse_readcount_file(readcount_file)

    elif abundance_file:
        logger.info("Using abundance file:%s" % abundance_file)
        if draft:
            abundance_dic_short = abundance.parse_abundance_file(abundance_file)
            complete_genomes_dic = {k: v for k, v in abundance_dic_short.items() if k not in draft}
            draft_dic = abundance.expand_draft_abundance(abundance_dic_short, draft)
            abundance_dic = {**complete_genomes_dic, **draft_dic}
        else:
            abundance_dic = abundance.parse_abundance_file(abundance_file)
    elif coverage_file:
        logger.warning("--coverage_file is an experimental feature")
        logger.warning("--coverage_file disables --n_reads")
        logger.info("Using coverage file:%s" % coverage_file)
        if draft:
            coverage_dic = abundance.parse_abundance_file(coverage_file)
            complete_genomes_dic = {k: v for k, v in coverage_dic.items() if k not in draft}
            draft_dic = abundance.expand_draft_abundance(coverage_dic, draft, mode="coverage")
            abundance_dic = {**complete_genomes_dic, **draft_dic}
        else:
            abundance_dic = abundance.parse_abundance_file(coverage_file)
    elif coverage in abundance_dispatch:
        # todo coverage distribution with --draft
        logger.info("Using %s coverage distribution" % coverage)
        if draft:
            abundance_dic = abundance.draft(
                genome_list,
                draft,
                abundance_dispatch[abundance_distribution],
                output,
                mode="coverage",
            )
        else:
            abundance_dic = abundance_dispatch[coverage](genome_list)
        if n_reads:
            n_reads = util.convert_n_reads(n_reads)
            logger.info("scaling coverage to %s reads" % n_reads)
            abundance_dic = abundance.coverage_scaling(n_reads, abundance_dic, genome_file, error_model.read_length)
        abundance.to_file(abundance_dic, output, mode="coverage")
    elif abundance_distribution in abundance_dispatch:
        logger.info("Using %s abundance distribution" % abundance_distribution)
        if draft:
            abundance_dic = abundance.draft(
                genome_list,
                draft,
                abundance_dispatch[abundance_distribution],
                output,
            )
        else:
            abundance_dic = abundance_dispatch[abundance_distribution](genome_list)
            abundance.to_file(abundance_dic, output)
    else:
        logger.error("Could not get abundance, or coverage or readcount information")
        sys.exit(1)

    return readcount_dic, abundance_dic


def write_mutations(mutations, mutations_handle):
    """Write mutations to a file

    Args:
        mutations (list): List of mutations.
        mutations_handle (file): File handle to write the mutations to.
    """
    for vcf_dict in mutations:
        mutations_handle.write(
            "\t".join(
                [
                    str(vcf_dict["id"]),
                    str(vcf_dict["position"] + 1),  # vcf files have 1-based index
                    ".",
                    vcf_dict["ref"],
                    str(vcf_dict["alt"]),
                    str(vcf_dict["quality"]),
                    "",
                    "",
                ]
            )
            + "\n"
        )


# ============================================================
# Annotation-aware extension
# ============================================================

# Annotation-aware public API

def simulate_and_add_info_to_reads(
    annotation_file,
    record,
    error_model,
    n_pairs,
    cpu_number,
    forward_handle,
    reverse_handle,
    mutations_handle,
    sequence_type="metagenomics",
    gc_bias=False,
    target_reading_frame=None,
    only_full_coding=False,
    annotation_handle=None,
    output_read_length=None,
):
    """
    Generate annotation-aware paired reads from a FASTA genome.

    CDS coordinates, strand, phase, and identifiers are loaded from a
    GTF/GFF annotation. R1 can optionally be restricted to a requested
    reading frame and to reads fully contained within a single CDS.

    For the common combination of ``only_full_coding=True`` and a specified
    reading frame, valid genomic start positions are calculated directly.
    Other annotation-aware modes use rejection sampling.

    Returns:
        int: Number of accepted read pairs written to the output handles.
    """
    valid_frames = {-3, -2, -1, 1, 2, 3}

    if (target_reading_frame is not None and target_reading_frame not in valid_frames):
        raise ValueError("target_reading_frame must be one of -3, -2, -1, 1, 2, or 3")

    if (target_reading_frame is not None or only_full_coding) and not annotation_file:
        raise ValueError("An annotation file is required for reading-frame or full-CDS filtering")

    if annotation_file and not os.path.isfile(annotation_file):
        raise FileNotFoundError(f"Annotation file not found: {annotation_file}")
    
    logger = logging.getLogger(__name__)
    trim_left, trim_right = _get_trim_lengths(
        error_model.read_length,
        output_read_length,
    )

    cds_regions = {}
    if annotation_file and os.path.exists(annotation_file):
        db = create_annotation_db(annotation_file)
        cds_regions = query_cds_regions(db)
        logger.info("Loaded CDS regions for %d chrom keys", len(cds_regions))
        logger.debug("Example CDS keys: %s", list(cds_regions.keys())[:5])

    genome_records = [
        r for r in SeqIO.parse(record, "fasta")
        if len(r.seq) > error_model.read_length
    ]
    logger.debug("Example FASTA record ids: %s", [r.id for r in genome_records[:5]])
    if not genome_records:
        logger.warning("No usable FASTA records found in %s", record)
        return 0

    if (sequence_type == "metagenomics" and only_full_coding and target_reading_frame is not None):
        eligible = []

        for rec in genome_records:
            ranges = _build_valid_start_ranges_full_coding_frame(
                rec,
                cds_regions,
                read_length=error_model.read_length,
                target_reading_frame=target_reading_frame,
                trim_left=trim_left,
                trim_right=trim_right,
            )

            n_valid = _count_start_ranges(ranges)
            if n_valid > 0:
                eligible.append((rec, ranges, n_valid))
                
        if not eligible:
            logger.warning(
                "No eligible contigs with full-coding frame=%s in %s",
                target_reading_frame,
                record,
            )
            return 0

        accepted_total = 0
        attempts_total = 0
        max_attempts_total = max(n_pairs * 20, 2000)

        while accepted_total < n_pairs and attempts_total < max_attempts_total:
            attempts_total += 1

            # Weight contigs by their number of valid start positions so that
            # valid genomic positions are sampled uniformly across the genome.
            rec, ranges, _ = random.choices(eligible, weights=[n_valid for _, _, n_valid in eligible], k=1)[0]
            forward_start = _sample_from_start_ranges(ranges)

            try:
                forward_record, reverse_record, mutations = simulate_annotated_read(
                    rec,
                    error_model,
                    accepted_total,
                    cpu_number,
                    sequence_type,
                    cds_regions,
                    forced_forward_start=forward_start,
                    output_read_length=output_read_length,
                )
            except AssertionError:
                continue

            if gc_bias:
                stitched_seq = forward_record.seq + reverse_record.seq
                gc_content = gc_fraction(stitched_seq)
                if not (0.4 < gc_content < 0.6 or np.random.rand() < 0.90):
                    continue

            # safety checks
            if not forward_record.annotations.get("fully_contained_in_cds", False):
                continue
            if int(forward_record.annotations.get("reading_frame", 0)) != int(target_reading_frame):
                continue

            _write_read_annotation(forward_record, annotation_handle)

            SeqIO.write(forward_record, forward_handle, "fastq-sanger")
            SeqIO.write(reverse_record, reverse_handle, "fastq-sanger")
            write_mutations(mutations, mutations_handle)

            accepted_total += 1

        if accepted_total < n_pairs:
            logger.warning(
                "Stopped early: accepted %d/%d pairs after %d total attempts",
                accepted_total,
                n_pairs,
                attempts_total,
            )
        return accepted_total

    # Generic rejection-sampling path for other annotation-aware modes.
    accepted_total = 0
    attempts_total = 0
    max_attempts_total = max(n_pairs * 200, 10000)
    per_record_attempts = 10

    while accepted_total < n_pairs and attempts_total < max_attempts_total:
        rec = random.choice(genome_records)

        for _ in range(per_record_attempts):
            attempts_total += 1
            if attempts_total >= max_attempts_total:
                break

            gen = annotated_reads_generator(
                n_pairs=1,
                record=rec,
                error_model=error_model,
                cpu_number=cpu_number,
                gc_bias=gc_bias,
                sequence_type=sequence_type,
                cds_regions=cds_regions,
                target_reading_frame=target_reading_frame,
                only_full_coding=only_full_coding,
                max_attempt_factor=1,
                output_read_length=output_read_length,
            )

            got_pair = False
            for forward_record, reverse_record, mutations in gen:
                got_pair = True

                _write_read_annotation(forward_record, annotation_handle)

                SeqIO.write(forward_record, forward_handle, "fastq-sanger")
                SeqIO.write(reverse_record, reverse_handle, "fastq-sanger")
                write_mutations(mutations, mutations_handle)

                accepted_total += 1
                break

            if got_pair:
                break

    if accepted_total < n_pairs:
        logger.warning(
            "Stopped early: accepted %d/%d pairs after %d total attempts",
            accepted_total, n_pairs, attempts_total
        )
    return accepted_total


def annotated_reads_generator(
    n_pairs,
    record,
    error_model,
    cpu_number,
    gc_bias,
    sequence_type,
    cds_regions=None,
    target_reading_frame=None,
    only_full_coding=False,
    max_attempt_factor=50,
    output_read_length=None,
):
    """
    Generate annotation-aware read pairs from one reference sequence.

    Read pairs are simulated using the selected InSilicoSeq error model.
    R1 may optionally be restricted to reads fully contained within a CDS
    and/or to a requested BLASTX-like reading frame.

    Reads that fail GC-bias, CDS-containment, or frame filters are rejected
    until the requested number of pairs is reached or the attempt limit is
    exceeded.
    """
    logger = logging.getLogger(__name__)

    accepted = 0
    attempts = 0
    max_attempts = max(n_pairs * max_attempt_factor, 1000)

    while accepted < n_pairs:
        attempts += 1
        if attempts > max_attempts:
            logger.warning(
                "Stopped early for %s after %d attempts. "
                "Accepted %d/%d pairs for target_reading_frame=%s",
                record.id, attempts, accepted, n_pairs, target_reading_frame
            )
            break

        try:
            forward, reverse, mutations = simulate_annotated_read(
                record,
                error_model,
                accepted,
                cpu_number,
                sequence_type,
                cds_regions or {},
                output_read_length=output_read_length,
            )
        except AssertionError:
            continue

        if gc_bias:
            stitched_seq = forward.seq + reverse.seq
            gc_content = gc_fraction(stitched_seq)
            if not (0.4 < gc_content < 0.6 or np.random.rand() < 0.90):
                continue
            
        if only_full_coding:
            if not forward.annotations.get("fully_contained_in_cds", False):
                continue

        if target_reading_frame is not None:
            if forward.annotations.get("reading_frame", 0) != target_reading_frame:
                continue

        yield (forward, reverse, mutations)
        accepted += 1


def simulate_annotated_read(
    record,
    error_model,
    i,
    cpu_number,
    sequence_type,
    cds_regions,
    forced_forward_start=None,
    output_read_length=None,
):
    """
    Simulate one paired-end read from a reference sequence and attach
    CDS-derived metadata to R1.

    R1 annotations include genomic coordinates, CDS overlap, CDS identity,
    strand, phase, and BLASTX-like reading frame. When
    ``output_read_length`` is specified, both mates are symmetrically
    trimmed after error simulation.

    Args:
        record: Reference SeqRecord.
        error_model: InSilicoSeq error model.
        i: Read-pair index.
        cpu_number: Worker identifier used in read names.
        sequence_type: ``metagenomics`` or ``amplicon``.
        cds_regions: CDS intervals indexed by reference sequence ID.
        forced_forward_start: Optional genomic start coordinate for R1.
        output_read_length: Optional final read length after trimming.

    Returns:
        Tuple containing R1, R2, and the combined mutation records.
    """
    logger = logging.getLogger(__name__)
    sequence = record.seq
    seq_len = len(sequence)
    header = record.id

    read_length = error_model.read_length
    trim_left, trim_right = _get_trim_lengths(read_length, output_read_length,)

    if error_model.fragment_length is not None and error_model.fragment_sd is not None:
        fragment_length = int(np.random.normal(error_model.fragment_length, 
                                               error_model.fragment_sd))
        insert_size = fragment_length - (read_length * 2)
    else:
        insert_size = error_model.random_insert_size()
        fragment_length = insert_size + (read_length * 2)

    # generate the forward read
    if forced_forward_start is None:
        try:
            assert read_length < len(record.seq)
            if sequence_type == "metagenomics":
                forward_start = random.randrange(0, seq_len - fragment_length)
            elif sequence_type == "amplicon":
                forward_start = 0
            else:
                raise RuntimeError(f"sequence type '{sequence_type}' is not supported")
        except AssertionError:
            raise
        except ValueError as e:
            logger.debug("%s shorter than template length for this ErrorModel:%s" % (record.id, e))
            forward_start = max(0, random.randrange(0, len(record.seq) - read_length))
    else:
        forward_start = int(forced_forward_start)
        if forward_start < 0 or (forward_start + read_length) > seq_len:
            raise AssertionError("forced_forward_start is outside the contig")

    forward_end = forward_start + read_length
    bounds = (forward_start, forward_end)

    forward_read_info = {
        "start": forward_start,
        "end": forward_end,
        "strand": "+",
        "chrom": str(header),
    }

    forward_read = Seq(str(sequence[forward_start:forward_end]))
    forward_ann = get_read_annotation(forward_read_info, cds_regions)
    forward_info = build_read_description(
        forward_read_info,
        cds_regions,
        forced_reading_frame=forward_ann["reading_frame"],
    )

    # create a perfect read
    forward = SeqRecord(
        forward_read, id="%s_%s_%s/1" % (header, i, cpu_number), description=forward_info)
    forward.annotations["mutations"] = []
    forward.annotations["original"] = str(forward.seq)
    _apply_read_annotation(forward, forward_ann)

    # add the indels, the qual scores and modify the record accordingly
    forward = error_model.introduce_indels(forward, "forward", sequence, bounds)
    forward = error_model.introduce_error_scores(forward, "forward")
    forward = error_model.mut_sequence(forward, "forward")

    if trim_left > 0 or trim_right > 0:
        trim_end = (len(forward.seq) - trim_right if trim_right > 0 else len(forward.seq))

        if "phred_quality" in forward.letter_annotations:
            trimmed_phred_quality = forward.letter_annotations["phred_quality"][trim_left:trim_end]
        else:
            trimmed_phred_quality = None

        forward.letter_annotations = {}
        forward.seq = Seq(str(forward.seq[trim_left:trim_end]))

        if trimmed_phred_quality is not None:
            forward.letter_annotations["phred_quality"] = trimmed_phred_quality

        forward_read_info["start"] += trim_left
        forward_read_info["end"] -= trim_right

        forward_ann = get_read_annotation(forward_read_info, cds_regions,)

        forward.description = build_read_description(
            forward_read_info,
            cds_regions,
            forced_reading_frame=forward_ann["reading_frame"],
        )
        _apply_read_annotation(forward, forward_ann)

    # generate the reverse read
    # assign start position reverse read
    # if sequence_type == metagenomics, get a start position based on insert_size
    # if sequence_type == amplicon, start position is the end of the read
    if sequence_type == "metagenomics":
        reverse_start = forward_end + insert_size
        reverse_end = reverse_start + read_length
    elif sequence_type == "amplicon":
        reverse_start = seq_len - read_length
        reverse_end = reverse_start + read_length
    else:
        raise ValueError(f"Sequence type {sequence_type} not known")
    if reverse_end > len(record.seq):
        # we use random insert when the modelled template length distribution
        # is too large
        reverse_end = random.randrange(read_length, len(record.seq))
        reverse_start = reverse_end - read_length

    bounds = (reverse_start, reverse_end)

    reverse_read_info = {"start": reverse_start, "end": reverse_end, 
                         "strand": "-", "chrom": str(header)}

    reverse_substring = sequence[reverse_start:reverse_end]
    rev_comp_seq = rev_comp(str(reverse_substring))

    backward_info = build_read_description(
        reverse_read_info,
        cds_regions,
        forced_reading_frame=forward.annotations["reading_frame"],
    )
    # create a perfect read
    reverse = SeqRecord(Seq(rev_comp_seq), id="%s_%s_%s/2" % (header, i, cpu_number), 
                        description=backward_info)
    reverse.annotations["reading_frame"] = forward.annotations["reading_frame"]
    reverse.annotations["is_cds"] = forward.annotations["is_cds"]
    reverse.annotations["chromosome"] = forward.annotations["chromosome"]
    reverse.annotations["start"] = reverse_read_info["start"]
    reverse.annotations["end"] = reverse_read_info["end"]
    reverse.annotations["pair_best_cds_id"] = forward.annotations["best_cds_id"]
    reverse.annotations["mutations"] = []
    reverse.annotations["original"] = str(reverse.seq)

    # add the indels, the qual scores and modify the record accordingly
    reverse = error_model.introduce_indels(reverse, "reverse", sequence, bounds)
    reverse = error_model.introduce_error_scores(reverse, "reverse")
    reverse = error_model.mut_sequence(reverse, "reverse")

    if trim_left > 0 or trim_right > 0:
        trim_end = (len(reverse.seq) - trim_right if trim_right > 0 else len(reverse.seq))

        if "phred_quality" in reverse.letter_annotations:
            trimmed_phred_quality = reverse.letter_annotations["phred_quality"][trim_left:trim_end]
        else:
            trimmed_phred_quality = None

        reverse.letter_annotations = {}
        reverse.seq = Seq(str(reverse.seq[trim_left:trim_end]))

        if trimmed_phred_quality is not None:
            reverse.letter_annotations["phred_quality"] = trimmed_phred_quality

        reverse_read_info["start"] += trim_left
        reverse_read_info["end"] -= trim_right

        reverse.description = build_read_description(
            reverse_read_info,
            cds_regions,
            forced_reading_frame=forward.annotations["reading_frame"],
        )
    combined_mutations = forward.annotations["mutations"] + reverse.annotations["mutations"]

    return (forward, reverse, combined_mutations)

# Annotation parsing and read annotation helpers

def create_annotation_db(annotation_file):
    """Create an in-memory gffutils database for a GFF/GTF annotation file."""
    return gffutils.create_db(
        annotation_file,
        dbfn=":memory:",
        force=True,
        merge_strategy="merge",
        keep_order=False,
        sort_attribute_values=False,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )


def _normalize_chrom_key(raw_chrom):
    return str(raw_chrom).strip()


def query_cds_regions(db):
    """
    Return:
      cds_regions[chrom] = list of
      (start0, end0, strand, phase, cds_id)

    start0 is 0-based inclusive
    end0 is 0-based exclusive
    """
    cds_regions = {}

    for feature in db.features_of_type("CDS"):
        chrom = _normalize_chrom_key(feature.chrom)
        start0 = max(0, int(feature.start) - 1)
        end0 = int(feature.end)
        strand = feature.strand
        phase = _safe_phase(getattr(feature, "frame", 0))

        attrs = getattr(feature, "attributes", {}) or {}
        cds_id = (
            (attrs.get("protein_id", [""])[0] if "protein_id" in attrs else "")
            or (attrs.get("locus_tag", [""])[0] if "locus_tag" in attrs else "")
            or (attrs.get("gene_id", [""])[0] if "gene_id" in attrs else "")
            or getattr(feature, "id", None)
            or ""
        )

        cds_regions.setdefault(chrom, []).append((start0, end0, strand, phase, cds_id))

    for chrom in cds_regions:
        cds_regions[chrom].sort(key=lambda x: x[0])

    return cds_regions


def _safe_phase(val):
    if val in (None, ".", ""):
        return 0
    try:
        return int(val)
    except Exception:
        return 0


def _merged_overlap_len(read_start, read_end, cds_list):
    """
    Return the union overlap length between [read_start, read_end) and CDS intervals.
    Prevents counting the same read bases multiple times when CDS features overlap.
    """
    overlaps = []

    for cds_item in cds_list:
        cds_start, cds_end = cds_item[0], cds_item[1]

        if cds_end <= read_start:
            continue
        if cds_start >= read_end:
            break

        ov_start = max(read_start, cds_start)
        ov_end = min(read_end, cds_end)
        if ov_start < ov_end:
            overlaps.append((ov_start, ov_end))

    if not overlaps:
        return 0

    overlaps.sort()
    merged_total = 0
    cur_s, cur_e = overlaps[0]

    for s, e in overlaps[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged_total += (cur_e - cur_s)
            cur_s, cur_e = s, e

    merged_total += (cur_e - cur_s)
    return merged_total


def get_cds_overlap_info(read_info, cds_regions):
    read_start = int(read_info["start"])
    read_end = int(read_info["end"])
    read_length = max(0, read_end - read_start)
    strand = read_info["strand"]
    chrom_key = _normalize_chrom_key(read_info["chrom"])

    cds_list = cds_regions.get(chrom_key, [])
    coding_bases = _merged_overlap_len(read_start, read_end, cds_list)

    if read_length > 0:
        overlap_percentage = min(100.0, (coding_bases / read_length) * 100.0)
    else:
        overlap_percentage = 0.0

    is_cds = coding_bases > 0

    return {
        "start": read_start,
        "end": read_end,
        "strand": strand,
        "chromosome": chrom_key,
        "coding_bases": coding_bases,
        "read_length": read_length,
        "overlap_pct": round(overlap_percentage, 3),
        "is_cds": is_cds,
    }


def find_best_cds_hit(chrom_key, read_start, read_end, cds_regions):
    best_overlap = 0
    best_hit = None

    for cds_start, cds_end, cds_strand, cds_phase, cds_id in cds_regions.get(
        chrom_key, []
    ):
        if cds_end <= read_start:
            continue
        if cds_start >= read_end:
            break

        overlap = min(read_end, cds_end) - max(read_start, cds_start)

        if overlap > best_overlap:
            best_overlap = overlap
            best_hit = (cds_start, cds_end, cds_strand, cds_phase, cds_id,)

    return best_overlap, best_hit


def compute_r1_reading_frame(read_start, read_end, best_hit):
    """
    Return the reference-derived BLASTX-like query frame of R1.

    The frame is determined from the source genomic coordinates and CDS phase,
    before considering sequencing errors introduced into the simulated read.

    Positive values (+1, +2, +3) indicate translation on the forward
    strand; negative values (-1, -2, -3) indicate translation on the
    reverse-complement strand. Returns 0 when no CDS overlaps the read.
    """
    if best_hit is None:
        return 0

    cds_start, cds_end, cds_strand, cds_phase, _ = best_hit

    if cds_strand == "+":
        offset = (read_start - cds_start - cds_phase) % 3
        plus_map = {0: 1, 1: 3, 2: 2}
        return plus_map[offset]

    if cds_strand == "-":
        offset = (cds_end - read_end - cds_phase) % 3
        minus_map = {0: -1, 1: -3, 2: -2}
        return minus_map[offset]

    return 0

def get_read_annotation(read_info, cds_regions, forced_reading_frame=None):
    read_start = int(read_info["start"])
    read_end = int(read_info["end"])
    strand = str(read_info["strand"])
    chrom_key = _normalize_chrom_key(read_info["chrom"])

    overlap_info = get_cds_overlap_info(read_info, cds_regions)
    best_overlap, best_hit = find_best_cds_hit(chrom_key, read_start, read_end, cds_regions)
    fully_contained_in_cds = (best_hit is not None and best_overlap == overlap_info["read_length"])

    cds_strand = "."
    cds_phase = "."
    best_cds_id = ""

    if best_hit is not None:
        _, _, cds_strand, cds_phase, best_cds_id = best_hit

    reading_frame = (
        int(forced_reading_frame)
        if forced_reading_frame is not None
        else compute_r1_reading_frame(read_start, read_end, best_hit)
    )

    return {
        "start": read_start,
        "end": read_end,
        "strand": strand,
        "chromosome": chrom_key,
        "is_cds": overlap_info["is_cds"],
        "overlap_pct": overlap_info["overlap_pct"],
        "coding_bases": overlap_info["coding_bases"],
        "read_length": overlap_info["read_length"],
        "fully_contained_in_cds": fully_contained_in_cds,
        "cds_strand": cds_strand,
        "cds_phase": cds_phase,
        "best_cds_id": best_cds_id,
        "reading_frame": int(reading_frame),
    }


def build_read_description(read_info, cds_regions, forced_reading_frame=None,):
    ann = get_read_annotation(
        read_info,
        cds_regions,
        forced_reading_frame=forced_reading_frame,
    )

    cds_info = "CDS" if ann["is_cds"] else "Non-CDS"

    return (
        f"Start:{ann['start']}\n"
        f"End:{ann['end']}\n"
        f"Strand:{ann['strand']}\n"
        f"Chromosome:{ann['chromosome']}\n"
        f"CDS_info:{cds_info}_{ann['overlap_pct']}%_coding\n"
        f"CDS_strand:{ann['cds_strand']}\n"
        f"CDS_phase:{ann['cds_phase']}\n"
        f"Reading_frame:{ann['reading_frame']}\n"
        f"Best_CDS:{ann['best_cds_id']}"
    )


def _apply_read_annotation(record, annotation):
    record.annotations["reading_frame"] = annotation["reading_frame"]
    record.annotations["is_cds"] = annotation["is_cds"]
    record.annotations["chromosome"] = annotation["chromosome"]
    record.annotations["start"] = annotation["start"]
    record.annotations["end"] = annotation["end"]
    record.annotations["overlap_pct"] = annotation["overlap_pct"]
    record.annotations["coding_bases"] = annotation["coding_bases"]
    record.annotations["read_length"] = annotation["read_length"]
    record.annotations["fully_contained_in_cds"] = annotation["fully_contained_in_cds"]
    record.annotations["best_cds_id"] = annotation["best_cds_id"]
    record.annotations["cds_strand"] = annotation["cds_strand"]
    record.annotations["cds_phase"] = annotation["cds_phase"]


# Frame-specific sampling helpers

def _record_chrom_key(record):
    """
    Use the raw FASTA seqid as the lookup key.

    For NCBI genomic FASTA/GTF/GFF, the stable join key is the sequence id
    itself (record.id / feature.chrom), not a chromosome number parsed from
    the human-readable description line.

    This is required for draft / scaffold / shotgun assemblies such as GCA
    where descriptions do not contain a real chromosome number and appending
    '_0' breaks the join to CDS annotations.
    """
    return _normalize_chrom_key(record.id)


def _get_trim_lengths(model_read_length, output_read_length):
    """
    Return the number of bases to trim from the left and right sides
    of a simulated read.
    If output_read_length is None, no trimming is applied.
    """
    model_read_length = int(model_read_length)

    if output_read_length is None:
        return 0, 0

    output_read_length = int(output_read_length)

    if output_read_length <= 0:
        raise ValueError("output_read_length must be a positive integer")

    if output_read_length > model_read_length:
        raise ValueError(
            f"output_read_length ({output_read_length}) cannot exceed "
            f"the error-model read length ({model_read_length})"
        )

    total_trim = model_read_length - output_read_length
    trim_left = total_trim // 2
    trim_right = total_trim - trim_left

    return trim_left, trim_right


def _count_start_ranges(ranges):
    total = 0
    for start, stop, step in ranges:
        if stop > start:
            total += (stop - start + step - 1) // step
    return total


def _sample_from_start_ranges(ranges):
    total = _count_start_ranges(ranges)
    if total <= 0:
        raise ValueError("No valid start positions to sample from")

    k = random.randrange(total)

    for start, stop, step in ranges:
        n = 0
        if stop > start:
            n = (stop - start + step - 1) // step
        if k < n:
            return start + k * step
        k -= n

    raise RuntimeError("Failed to sample start from ranges")


def _build_valid_start_ranges_full_coding_frame(
    record,
    cds_regions,
    read_length,
    target_reading_frame,
    trim_left=0,
    trim_right=0,
):
    """
    Build forward-start ranges that guarantee:
      1. trimmed R1 is fully inside one CDS
      2. trimmed R1 has the requested BLASTX-like reading frame

    Returned ranges are tuples: (start, stop_exclusive, step)
    """
    chrom_key = _record_chrom_key(record)
    cds_list = cds_regions.get(chrom_key, [])

    if not cds_list:
        return []

    max_forward_start = len(record.seq) - read_length
    if max_forward_start < 0:
        return []

    ranges = []

    # compute_r1_reading_frame() expresses the genomic displacement as the
    # negative modulo-3 equivalent of the query-relative translation offset.
    # These mappings convert between that internal offset and BLASTX-like
    # query frames.
    plus_inv = {1: 0, 3: 1, 2: 2}
    minus_inv = {-1: 0, -3: 1, -2: 2}

    for cds_start, cds_end, cds_strand, cds_phase, cds_id in cds_list:
        # trimmed interval is:
        # [forward_start + trim_left, forward_start + read_length - trim_right)
        lo = max(0, cds_start - trim_left)
        hi = min(max_forward_start, cds_end - (read_length - trim_right))

        if lo > hi:
            continue

        tf = int(target_reading_frame)

        if cds_strand == "+":
            if tf not in plus_inv:
                continue
            wanted_offset = plus_inv[tf]
            residue = (wanted_offset - trim_left + cds_start + cds_phase) % 3
        elif cds_strand == "-":
            if tf not in minus_inv:
                continue
            wanted_offset = minus_inv[tf]
            residue = (cds_end - cds_phase - (read_length - trim_right) - wanted_offset) % 3

        else:
            continue

        first = lo + ((residue - lo) % 3)
        if first <= hi:
            ranges.append((first, hi + 1, 3))

    return ranges

# Output helper

def _write_read_annotation(record, handle):
    if handle is None:
        return

    handle.write(
        f"{record.id}\t"
        f"{record.annotations.get('chromosome', '')}\t"
        f"{record.annotations.get('start', '')}\t"
        f"{record.annotations.get('end', '')}\t"
        f"{record.annotations.get('is_cds', '')}\t"
        f"{record.annotations.get('fully_contained_in_cds', '')}\t"
        f"{record.annotations.get('overlap_pct', '')}\t"
        f"{record.annotations.get('coding_bases', '')}\t"
        f"{record.annotations.get('read_length', '')}\t"
        f"{record.annotations.get('reading_frame', '')}\t"
        f"{record.annotations.get('cds_strand', '')}\t"
        f"{record.annotations.get('cds_phase', '')}\t"
        f"{record.annotations.get('best_cds_id', '')}\n"
    )
