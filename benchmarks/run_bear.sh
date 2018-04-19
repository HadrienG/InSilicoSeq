#!/usr/bin/env bash

usage="$(basename "${0}") [-h] [-g fasta] [-r reads] [-o output] --
        Wrapper script around BEAR

        options:
    		-g <multifasta> reference genomes in fasta format
            -r <reads> read file in fastq format
            -o <str> prefix for output files"


while getopts "hg:r:o:" option
    do
        case "${option}" in
            h) echo "${usage}"
            exit 0
                ;;
            g) reference=$(readlink -f "${OPTARG}")
                ;;
            r) reads=$(readlink -f "${OPTARG}")
                ;;
            o) output="${OPTARG}"
                ;;
            :) printf "missing argument for -%s\n" "${OPTARG}" >&2
            echo "${usage}" >&2
            exit 1
                ;;
            \?) printf "illegal option -%s\n" "${OPTARG}" >&2
            echo "${usage}">&2
            exit 1
                ;;
        esac
    done
shift $((OPTIND-1))

if [ -z "${reference+x}" ]
    then
        printf "%s\n\nEmply value for -g\n" "${usage}" >&2
        exit 1
fi
if [ -z "${reads+x}" ]
    then
        printf "%s\n\nEmply value for -r\n" "${usage}" >&2
        exit 1
fi
if [ -z "${output+x}" ]
    then
        printf "%s\n\Bad value for -o\n" "${usage}" >&2
        exit 1
fi

module load bear

parametric_abundance.pl "${reference}" med "${output}.txt"
drisee.py --percent -n 3 -b 10000 -x 10000 -t fastq -s 500000 "${reads}" \
    "${output}_err"
error quality.pl "${reads}" "${output}_err.uc"
generate_reads.py -r "${reference}" -a "${output}.txt" -o "${output}.fastq" \
    -t 500000 -l 300 -i 350 -s 10
trim_reads.py -i "${reads}" -f "${reference}" -o "${output}.fastq" -r \
    "${output}_err.per" -q "${output}_err.err.qual" -m "${output}_err.err.matr"
