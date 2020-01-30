#!/bin/bash

ORIGINAL_DIR=$(pwd)
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

## if no arguments, show manual
if [[ $# -eq 0 ]]; then
    man -l ${SCRIPT_DIR}/MANUAL_get_seqs.1
    exit 1
fi


# BED_DEFAULT='/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed'
FEATURE_DEFAULT='gene'
REFERENCE_DEFAULT="${SCRIPT_DIR}/data/1135acc.csv"
DOMAIN_DEFAULT=''
DOMAIN_F_DEFAULT=''

while (( "$#" )); do
    case "$1" in
        -g|--gene) GENE="${2}";;
        -a|--acc) ACC="${2}";;
        -i|--input) ACCS_F="${2}";;
        -f|--feature) FEATURE="${2}";;
        -d|--dir) DIR="${2}";;
        -c|--chrom) CHROM="${2}";;
        -s|--start) START="${2}";; ## 1-indexed
        -e|--end) END="${2}";; ## 1-indexed
        -b|--bed) BED="${2}";;
        -r|--reference) REFERENCE="${2}";;
        -o|--out) FOUT="${2}";;
        --domain) DOMAIN="${2}";;
        --domain-file) DOMAIN_F="${2}";;
        --domain-db) RPS_DB="${2}";;
        --input-dir) INPUT_DIR="${2}";;
        --bed-out) FOUT_BED="${2}";;
        --merge) MERGE='True';;
        --translate) TRANSLATE='True';;
        --no-header) HEADER='False';;
        --complete) COMPLETE='True';;
        --adjust-dir) ADJ_DIR='True';;
        --by-gene) BY_GENE='True';;
        --qname-dname) QDNAME="${2}";;
        --qstart-qend) QSTARTEND="${2}";;
        --no-bed) NO_BED='True';;
        -h|--help) man -l ${SCRIPT_DIR}/MANUAL_get_seqs.1; exit 0;;
        --readme) cat ${SCRIPT_DIR}/README; exit 0;;
    esac
    shift
done

DIR="${DIR:-${DIR_DEFAULT}}"
BED="${BED:-${BED_DEFAULT}}"
FEATURE="${FEATURE:-${FEATURE_DEFAULT}}"
REFERENCE="${REFERENCE:-${REFERENCE_DEFAULT}}"
HEADER="${HEADER:-True}"
COMPLETE="${COMPLETE:-False}"
MERGE="${MERGE:-False}"
TRANSLATE="${TRANSLATE:-False}"
SINGLE_FILE="${SINGLE_FILE:-False}"
ADJ_DIR="${ADJ_DIR:-False}"
BY_GENE="${BY_GENE:-False}"
DOMAIN="${DOMAIN:-${DOMAIN_DEFAULT}}"
DOMAIN_F="${DOMAIN_F:-${DOMAIN_F_DEFAULT}}"

## throw error if directory not provided
if [ -z "${DIR}" ]; then
    echo "Directory required. Please provide directory using '-d <path to directory>'"
    exit 1
## throw error if ACCS_F, ACC, and INPUT_DIR are not provided
elif [ -z "${ACCS_F}" ] && [ -z "${ACC}" ] && [ -z "${INPUT_DIR}" ]; then
    echo "Accession name(s) or number(s) required. Please provide a file of accession name(s)/number(s) using '-i <path to file>' or provide a single accession name or number using '-a <accession name/number>', or a directory of files of accession numbers (named according to gene ID) using '--input-dir <path to directory>'"
    exit 1
## else if any combination of 2 of ACCS_F, ACC, and INPUT_DIR are provided for some reason
elif (! [ -z "${ACCS_F}" ] && ! [ -z "${ACC}" ]) || \
         (! [ -z "${ACCS_F}" ] && ! [ -z "${INPUT_DIR}" ]) || \
         (! [ -z "${ACC}" ] && ! [ -z "${INPUT_DIR}" ]); then
    echo "Please only use either '-a <accession name/number>', '-i <path to file>', or '--input-dir <path to directory>', and not 2 or more at the same time. These parameters are mutually exclusive."
    exit 1
## throw error if GENE or INPUT_DIR or CHROM+START+END are not provided
elif [ -z "${GENE}" ] && [ -z "${INPUT_DIR}" ] && \
         ([ -z "${CHROM}" ] || [ -z "${START}" ] || [ -z "${END}" ]); then
    echo "Gene ID(s) or genome range required. Please provide gene ID(s) using '-g <gene ID(s)>', a directory of files of accession numbers (name according to gene ID) using '--input-dir <path to directory>', or a genome range using '-c <chromosome> -s <start position (inclusive)> -e <end position (inclusive)>'"
    exit 1
## else if any combination of 2 of GENE or INPUT_DIR or CHROM+START+END are provided for some reason
elif (! [ -z "${GENE}" ] && ! [ -z "${INPUT_DIR}" ]) || \
         (! [ -z "${GENE}" ] && (! [ -z "${CHROM}" ] || ! [ -z "${START}" ] || ! [ -z "${END}" ])) || \
         (! [ -z "${INPUT_DIR}" ] && (! [ -z "${CHROM}" ] || ! [ -z "${START}" ] || ! [ -z "${END}" ])); then
    echo "Please only use either '-g <gene ID(s)>', '--input-dir <path to directory>', or '-c <chromosome> -s <start position (inclusive)> -e <end position (inclusive)>', and not 2 or more at the same time. These parameters are mutually exclusive."
    exit 1
## else if genomic range is provided but information is incomplete
elif [ -z "${GENE}" ] && [ -z "${INPUT_DIR}" ] && \
         ([ -z "${CHROM}" ] || [ -z "${START}" ] || [ -z "${END}" ]) && \
         (! [ -z "${CHROM}" ] || ! [ -z "${START}" ] || ! [ -z "${END}" ]); then
     echo "Please provide all of the required parameters to define a genomic range using '-c <chromosome>', '-s <start position(inclusive)> -e <end position (inclusive)>'"
fi

## convert values of variables storing file and directory paths to absolute paths
path_vars_to_convert=( "GENE" "ACCS_F" "DIR" "BED" "REFERENCE" "FOUT" "DOMAIN_F" "RPS_DB" "INPUT_DIR" "FOUT_BED" )
for varname in ${path_vars_to_convert[@]}; do
    if ! [ -z "${!varname}" ] && ([ -f "${!varname}" ] || [ -d "${!varname}" ]); then
        eval ${varname}="$(realpath ${!varname})"
    fi
done

## move to output dir, create if doesn't exist
mkdir -p ${DIR}
cd ${DIR}
DIR=$(pwd)
echo "Output files will be generated in $(pwd)"

start_t=$(date '+%Y-%m-%d %T') ## for combining files
## ## start parsing data
extra=''
if ! [ -z "${QDNAME}" ]; then
    extra+=", qname_dname=${QDNAME}"
fi
if ! [ -z "${QSTARTEND}" ]; then
    extra+=", qstart_qend=${QSTARTEND}"
fi
## if not using --input-dir
if [ -z "${INPUT_DIR}" ]; then
    ## if accessions are provided (i.e. not requesting ref Col-0), reformat them into new variable $ACCS_INPUT
    if ! [ -z "${ACCS_F}" ] || (! [ -z "${ACC}" ] && ! [ "${ACC}" == 'ref' ]); then
        if [ -z "${ACCS_F}" ] && ! [ -z "${ACC}" ]; then
            ACCS_INPUT="('$(echo ${ACC} | sed 's/,/\x27,\x27/g')',)"
            ## else if a file of accessions is provided
        elif ! [ -z "${ACCS_F}" ] && [ -z "${ACC}" ]; then
            sed -i 's/\r//g' ${ACCS_F}
            ACCS_INPUT="('$(awk 'NF > 0' ${ACCS_F} | awk 'NR>1{print PREV} {PREV=$0} END{printf("%s",$0)}' | cat | tr '\n' ',' | sed 's/,/\x27,\x27/g')',)"
        fi
        echo "ACCS_INPUT: ${ACCS_INPUT}"
    fi
    ## if genomic range provided
    if ! [ -z "${CHROM}" ]; then
        ## if 'ref' is provided as $ACC
        if [ -z "${ACCS_F}" ] && ! [ -z "${ACC}" ] && [ "${ACC}" == 'ref' ]; then
            python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_ref_by_range('${CHROM}', $(( ${START}-1 )), ${END}, '${DIR}')"
        else
            python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_1001pseudogenome_by_range(${ACCS_INPUT}, '${CHROM}', ${START}, ${END}, '${DIR}')"
        fi
    ## else if genes provided
    else
        if [ -f "${GENE}" ]; then
            readarray -d ',' -t parsed_genes <<< "$(tr '\n' ',' < ${GENE} | sed 's/,\+$//g')"
        else
            readarray -d , -t parsed_genes <<< "$(echo ${GENE} | sed 's/\n//g')"
        fi
        for parsed_gene in "${parsed_genes[@]}"; do
            ## if 'ref' is provided as $ACC
            if [ -z "${ACCS_F}" ] && ! [ -z "${ACC}" ] && [ "${ACC}" == 'ref' ]; then
                python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_ref_by_gene(\"\"\"${parsed_gene}\"\"\".replace('\n',''), '${FEATURE}', '${DIR}', bed='${BED}', complete=(${COMPLETE}), domain='${DOMAIN}', domain_f='${DOMAIN_F}', merge=(${MERGE}), translate=(${TRANSLATE}), adj_dir=(${ADJ_DIR}), by_gene=(${BY_GENE})${extra})"
            else
                python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_1001pseudogenome_by_gene(${ACCS_INPUT}, \"\"\"${parsed_gene}\"\"\".replace('\n',''), '${FEATURE}', '${DIR}', bed='${BED}', ref_file='${REFERENCE}', header=(${HEADER}), complete=(${COMPLETE}), domain='${DOMAIN}', domain_f='${DOMAIN_F}', adj_dir=(${ADJ_DIR})${extra})"
            fi
        done
    fi
else
    for file in ${INPUT_DIR}/AT?G*; do
        fname=$(basename -- "${file}")
        GENE="${fname%.*}"
        sed -i 's/\r//g' ${file}
        ACCS_INPUT="('$(awk 'NF > 0' ${file} | awk 'NR>1{print PREV} {PREV=$0} END{printf("%s",$0)}' | cat | tr '\n' ',' | sed 's/,/\x27,\x27/g')')"
        python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_1001pseudogenome_by_gene(${ACCS_INPUT}, '${GENE}', '${FEATURE}', '${DIR}', bed='${BED}', ref_file='${REFERENCE}', header=(${HEADER}), complete=(${COMPLETE}), domain='${DOMAIN}', domain_f='${DOMAIN_F}', adj_dir=(${ADJ_DIR})${extra})"
    done
fi

## cat .fasta files if output file path provided
if ! [ -z "${FOUT}" ]; then
    cat $(find ${DIR}/*.fasta -maxdepth 1 -type f -newermt "${start_t}") > ${DIR}/tmp.fa
    rm $(find ${DIR}/*.fasta -maxdepth 1 -type f -newermt "${start_t}")
    if [ $(dirname "${FOUT}") == '.' ]; then
        fout=${DIR}/${FOUT}
    else
        fout=${FOUT}
    fi
    mv ${DIR}/tmp.fa ${fout}
    echo "FASTA files successfully concatenated into ${fout}"
fi


## delete .bed files if no-bed raised
if ! [ -z "${NO_BED}" ]; then
    rm $(find ${DIR}/bed/*.bed -maxdepth 1 -type f -newermt "${start_t}")
    find ${DIR}/bed -type d -empty -delete ## remove ./bed if empty
elif ! [ -z "${FOUT_BED}" ]; then ## cat .bed files if output .bed path provided
    cat $(find ${DIR}/bed/*.bed -maxdepth 1 -type f -newermt "${start_t}") > ${DIR}/bed/tmp.txt
    rm $(find ${DIR}/bed/*.bed -maxdepth 1 -type f -newermt "${start_t}")
    if [ $(dirname "${FOUT_BED}") == '.' ]; then
        fout=${DIR}/${FOUT_BED}
    else
        fout=${FOUT_BED}
    fi
    mv ${DIR}/bed/tmp.txt ${fout}
    find ${DIR}/bed -type d -empty -delete ## remove ./bed if empty
    echo "BED files successfully concatenated into ${fout}"
fi

exit 0
