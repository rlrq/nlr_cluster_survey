#!/bin/bash

run_blast6 () {
    params=${@}
    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            -outfmt) local executable="$(echo ${params} | sed "s/${2}/\x27${2}\x27/")";
                     local header="$(echo ${2} | sed 's/^[[:digit:]]\+ //' | tr ' ' '\t')";;
            -out) local fout="${2}";
                  local tmp="${fout}_tmp.txt";;
        esac
        shift
    done
    eval "${executable}" ## execute
    printf "${header}\n$(grep -v '^#' ${fout})" > $tmp; mv $tmp $fout
}
