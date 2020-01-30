import os
import re
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from data_manip import *

# cdd_versions_fname = "/mnt/chaelab/rachelle/data/cdd/cdd.versions.tsv"
## cdd_versions_fname should point to a tsv with headers (important columns are: pssmid, shortn)

def map_pssmid(fname, pssm_col, fout, cdd_versions_fname, header = False, output_header = False, sep = '\t',
               new_col = "name", pattern = "\d+"):
    with open(fname, 'r') as f:
        data = [x[:-1].split(sep) for x in f.readlines() if x]
    ## parse header (if present)
    if header and \
       re.search(pattern, data[0][pssm_col]) and \
       re.search(pattern, data[0][pssm_col]).group(0).isnumeric():
        print("Number detected in header row, PSSM-Id column. Parsing first row as data.")
        header = False
    elif header:
        data_header = data[0]
        data = data[1:]
    ## read cdd_versions
    with open(cdd_versions_fname, 'r') as f:
        cdd_data = [x[:-1].split('\t') for x in f.readlines() if x]
        cdd_get = make_custom_get(cdd_data[0])
        cdd_data = {cdd_get(x, "pssmid"): cdd_get(x, "shortname") for x in cdd_data[1:]}
    ## map shortname to pssmid
    data = [x + [cdd_data.get(int(re.search(pattern, x[pssm_col]).group(0)), "NA")] for x in data]
    ## write out
    if header and output_header:
        data = [data_header + [new_col]] + data
    with open(fout, "w+") as f:
        f.write('\n'.join(['\t'.join([str(y) for y in x]) for x in data]))
    return
