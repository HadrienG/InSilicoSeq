#!/usr/bin/env bash

usage="$(basename "${0}") [-h] [-g fasta] [-r reads] [-o output] --
        Wrapper script around FASTQSim

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


FASTQcharacterize.sh -input "${reads}" -reference 
FASTQspike.sh
