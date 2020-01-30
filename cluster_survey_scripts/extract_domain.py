import os
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from data_manip import *

def get_dat(make_fname, cluster, *d):
    with open(make_fname(cluster, d), 'r') as f:
        data = [x[:-1].split('\t') for x in f.readlines()]
    header = data[0]
    data = data[1:]
    return (header, data)

def make_custom_get(header, parse_num = True):
    def get_col(colname, data):
        return [get_col_in_row(x, colname) for x in data]
    def get_col_in_row(row, colname):
        output = row[header.index(colname)]
        if isinstance(output, (list, tuple)):
            if [str(x).isdigit() for x in output].count(True) == len(output):
                output = [int(x) for x in output]
            elif [str(x).replace('.','',1).replace('-','',1).isdigit() \
                  for x in output].count(True) == len(output):
                output = [float(x) for x in output]
            return output
        else:
            return output if not str(output).replace('.','',1).replace('-','',1).isdigit() else \
                float(output) if not str(output).isdigit() else int(output)
    def helper(data, *colnames):
        if isinstance(data[0], (list, tuple)):
            output = [get_col(colname, data) for colname in colnames]
            return output[0] if len(output) == 1 else [[output[r][c] for r in range(len(output))] \
                                                       for c in range(len(output[0]))]
        else:
            output = [get_col_in_row(data, colname) for colname in colnames]
            return output[0] if (len(output) == 1) else output
    return helper

def merge_domain_and_write_and_get_seq_multi(combos, fasta, header = [], colnames = [], adj_dir = False,
                                             **for_merge_domain_and_write_multi):
    if not header:
        with open(make_fname(combos[0][0], combos[0][1:]), 'r') as f:
            header = f.readlines()[0].split('\t')
    if not colnames:
        colnames = header
    fin_d = merge_domain_and_write_multi(combos, header = header, colnames = colnames, adj_dir = adj_dir,
                                         **{k: v for k, v in for_merge_domain_and_write_multi.items() if
                                            k != "make_domain_seq_fout"})
    print("Writing sequences to file")
    make_domain_seq_fout = for_merge_domain_and_write_multi["make_domain_seq_fout"]
    for k, fin in fin_d.items():
        fout = make_domain_seq_fout(*k)
        print("Writing " + fout)
        get_domain_seqs(fin, fasta, fout, header = colnames + (["strand"] if adj_dir else []), domain = k[1], adj_dir = adj_dir)

def merge_domain_and_write_multi(combos, **for_merge_domain_and_write):
    fout_d = {}
    for combo in combos:
        print("Processing " + str(combo))
        fout_d = {**fout_d, **merge_domain_and_write(*combo, **for_merge_domain_and_write)}
    return fout_d

def merge_domain_and_write(cluster, *domains, min_len = {}, min_cds_len = {}, adj_dir = False,
                           colnames = [], header = [], **for_merge_domain):
    make_fname = for_merge_domain["make_fname"]
    header, dat = get_dat(make_fname, cluster, *domains)
    fout_d = {}
    make_merged_fout = for_merge_domain["make_merged_fout"]
    for d in domains:
        merged = merge_domain(dat, d, min_len = min_len.get(d, 0), min_cds_len = min_cds_len.get(d, 0),
                              header = header, colnames = colnames, adj_dir = adj_dir,
                              **{k: v for k, v in for_merge_domain.items()
                                 if k not in ("make_merged_fout", "make_fname")})
        fout = make_merged_fout(cluster, d)
        write_table(merged, fout, header = list(colnames) + (["strand"] if adj_dir else []) + \
                                           (["cds.len"] if min_cds_len.get(d, 0) else [0]))
        fout_d[(cluster, d)] = fout
    return fout_d

## merge overlapping ranges if same domain
def merge_domain(data, domain, buffer_range = 100, min_id = 90, match_containing = False, min_len = 100,
                 header = [], adj_dir = False, min_cds_len = 0, equivalents = {},
                 colnames = ("contig", "hit.start", "hit.end", "X..identity", "accID", "domain",
                             "q..frame", "s..frame")):
    get = make_custom_get(header)
    dat = [get(x, *colnames) for x in data]
    get = make_custom_get(colnames)
    dat_complete = [x for x in dat if ((get(x, "domain")==domain)or(match_containing and (domain in get(x,"domain"))))]
    dat_complete.sort(key = lambda x: (get(x, "contig"), min(get(x, "hit.start", "hit.end"))))
    if len(dat_complete) == 0:
        return []
    
    def merge(dat, min_len):
        output = []
        last = {equivalents.get(x, x): get(dat[0], x) for x in colnames}
        s, e = last["hit.start"], last["hit.end"]
        last["hit.start"], last["hit.end"] = min(s,e), max(s,e)
        plus_minus = [0, 0]
        for i, entry in enumerate(dat[1:]):
            curr = {equivalents.get(x, x): get(entry, x) for x in colnames}
            s, e = curr["hit.start"], curr["hit.end"]
            curr["hit.start"], curr["hit.end"] = min(s,e), max(s,e)
            ## check if current entry overlaps with stored merged entries (same contig + overlapping range)
            if curr["contig"] != last["contig"] or curr["hit.start"] > (last["hit.end"] + buffer_range):
                ## if no overlap, check if last merged entries meet minimum len and id requirement
                if last["X..identity"] >= min_id and last["alignment.length"] >= min_len:
                    last["domain"] = domain
                    if adj_dir:
                        last["strand"] = '-' if plus_minus[-1] > plus_minus[0] else '+'
                    output.append([last[equivalents.get(x, x)]
                                   for x in colnames + (["strand"] if adj_dir else [])])
                ## update last merged entries
                last, to_revcomp, plus_minus = curr, [0,0], [0,0] ## reset plus_minus
            ## update stored merged entries with data from current entry
            else:
                last["hit.end"] = max(curr["hit.end"], last["hit.end"])
                last["X..identity"] = max(curr["X..identity"], last["X..identity"])
                last["alignment.length"] = last["hit.end"] - last["hit.start"] + 1
                ## store whether +/- or +/+ is majority; use this to decide whether to rev_comp at final step
                if adj_dir:
                    if (int(curr["q..frame"]) > 0 and int(curr["s..frame"]) < 0): ## check if -/+ alignment
                        plus_minus[-1] += curr["alignment.length"]
                    else:
                        plus_minus[0] += curr["alignment.length"]
            ## if is last entry and meets minimum id and len rquirement, add to output
            if i == len(dat)-2 and last["X..identity"] >= min_id and last["alignment.length"] >= min_len:
                last["domain"] = domain
                if adj_dir:
                    last["strand"] = '-' if (plus_minus[-1] > plus_minus[0]) else '+'
                output.append([last[equivalents.get(x, x)]
                               for x in colnames + (["strand"] if adj_dir else [])])
        return output
    
    complete_output = merge(dat_complete, min_len)
    
    if min_cds_len <= 0:
        return [x + [NaN] for x in complete_output]
    else:
        dat_cds = [x for x in dat if (domain in get(x, "domain") and "CDS" in get(x, "domain"))]
        if len(dat_cds) == 0:
            print("No CDS detected. Returning empty data.")
            return []
        dat_cds.sort(key = lambda x: (get(x, "contig"), min(get(x, "hit.start", "hit.end"))))
        cds_output = sorted(merge(dat_cds, 0),
                            key = lambda x: (get(x, "contig"), min(get(x, "hit.start", "hit.end"))))
        output = []
        def overlap_size(a, b):
            return max(0, min(max(a), max(b)) - max(min(a), min(b)) + 1)
        ## get CDS-complete overlap for each merged complete range
        cds_i_start = 0
        complete_output.sort(key = lambda x: (get(x, "contig"), min(get(x, "hit.start", "hit.end"))))
        for entry in complete_output:
            overlap, contig, start, end = 0, *get(entry, "contig", "hit.start", "hit.end")
            start, end = min(start, end), max(start, end)
            cds_i = cds_i_start
            ## find first overlapping cds entry
            while cds_i < len(cds_output) - 1 and (get(cds_output[cds_i], "contig") != contig or \
                  (get(cds_output[cds_i], "contig") == contig and \
                   max(get(cds_output[cds_i], "hit.end", "hit.start")) < start)):
                cds_i += 1
            ## update cds_last_start so we don't keep searching earlier cds entries
            cds_i_start = cds_i_start if cds_i >= len(cds_output) - 1 else cds_i
            ## while cds entries, overlap, get total overlap size
            while cds_i < len(cds_output) and \
                  get(cds_output[cds_i], "contig") == contig and \
                  min(get(cds_output[cds_i], "hit.start", "hit.end")) <= end:
                overlap += overlap_size(get(cds_output[cds_i], "hit.start", "hit.end"), (start, end))
                cds_i += 1
            if overlap >= min_cds_len:
                output.append(entry + [overlap])
        return output

def write_table(data, fout, header = [], sep = '\t'):
    if header:
        data = [header] + data
    to_write = '\n'.join(['\t'.join([str(y) for y in x]) for x in data])
    f = open(fout, "w+")
    f.write(to_write)
    f.close()

## read table in and get sequences
def get_domain_seqs(domain_f, fasta, fout, header = [], domain = '', adj_dir = False):
    ## get domain ranges
    with open(domain_f, 'r') as f:
        domain_dat = [x[:-1].split('\t') for x in f.readlines()]
    get = make_custom_get(domain_dat[0])
    domain_dat = domain_dat[1:]
    ## parse fasta file
    from fasta_manip import fasta_to_dict
    seqs = fasta_to_dict(fasta)
    ## get sequences
    output = {}
    num = 1
    for i, entry in enumerate(domain_dat):
        ## note: domain_dat is 1-indexed, start and end inclusive, but Python splicing is 0-indexed
        seq = seqs[get(entry, "contig")][get(entry, "hit.start") - 1:get(entry, "hit.end")]
        if adj_dir and get(entry, "strand") == '-':
            seq = seq.reverse_complement()
        output['|'.join([str(x) for x in \
                         # (get(entry, "contig", "accID", "group") + [domain, num] + \
                         (get(entry, "contig", "accID") + [domain, num] + \
                          ['-'.join([str(x) for x in get(entry, "hit.start", "hit.end")])] + \
                          (["revcomp"] if adj_dir and get(entry, "strand") == '-' else []))])] = seq
        num += 1
    from fasta_manip import dict_to_fasta
    dict_to_fasta(output, fout)

def extract_domains_multi(coi_combos, vdw_contigs, cols,
                          make_fname, equivalents = {}, **kwargs):
    with open(make_fname(coi_combos[0][0], coi_combos[0][1:]), 'r') as f:
        header = f.readlines()[0][:-1].split('\t')
        get = make_custom_get(header)
    merge_domain_and_write_and_get_seq_multi(coi_combos, vdw_contigs, match_containing = False,
                                             buffer_range = 150, min_id = 85, colnames = cols,
                                             min_cds_len = {"TIR": 240, "RX-CC_like": 150, "NB-ARC":200},
                                             adj_dir = True, header = header, equivalents = equivalents,
                                             make_fname = make_fname, **kwargs)
    
