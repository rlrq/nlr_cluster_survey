import os
import Bio
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from itertools import groupby
from fasta_manip import fasta_to_dict, dict_to_fasta, remove_empty_pos, remove_incomplete_pos
from data_manip import *

IGNORE_GAPS = True
GAP_CHAR = '-'


def calculate_pi_multi(domains = ("TIR", "NB-ARC"), **kwargs):
    for d in domains:
        for ignore_gaps in (True, False):
            for cds_only in (True, False):
                print(f"Ignore gaps: {ignore_gaps} || CDS-only: {cds_only}")
                calculate_pi(domain = d, IGNORE_GAPS = ignore_gaps, cds_only = cds_only, **kwargs)
    return

def calculate_pi(domain, IGNORE_GAPS, GAP_CHAR = '-', preview = False, cds_only = False, in_all_ref = False,
                 prefix = "nlr165_AL70", fout_dir = '',
                 make_seq_fname = lambda domain: f"/mnt/chaelab/rachelle/hap_tags/results/nlr165/alignment/nlr165_col0-AL70-Alyrata_{domain}_mafft.fa",
                 make_pred_fname = lambda domain: f"/mnt/chaelab/rachelle/hap_tags/results/nlr165/predicted_identity/nlr165_AL70_{domain}_predictedIdentity.txt"):
    
    dat_seqs = fasta_to_dict(make_seq_fname(domain))
    with open(make_pred_fname(domain), 'r') as f:
        id_dat = [x[:-1].split('\t') for x in f.readlines()]
    
    get = make_custom_get(id_dat[0])
    ids = id_dat[1:]
    import itertools
    gene_cluster = set(list(itertools.chain.from_iterable([list(zip(get(x, "gene").split(';'), get(x, "cluster").split(';'))) for x in ids])))
    gene_cluster = {x[0]: x for x in gene_cluster}
    genes = set(gene_cluster.keys())
    
    pi_raw = {}
    pi_totalAlnLen = {}
    pi_pairAlnLen = {}
    
    for i, gene in enumerate(genes):
        print(i + 1, gene)
        
        ## get sequences assigned to the gene
        seq_ids = set(get(entry, "contig") for entry in ids
                      if gene in get(entry, "gene"))
        seqs_raw = {seq_id: seq for seq_id, seq in dat_seqs.items() if seq_id in seq_ids}
        
        ## remove empty positions
        if cds_only:
            ref_names = tuple(seq_id for seq_id in dat_seqs.keys() if
                              f"Col-0_ref|{gene}" in seq_id and "complete" not in seq_id)
            if not ref_names:
                continue
            from fasta_manip import trim_alignment_to_seqs
            seqs = trim_alignment_to_seqs(seqs = {**seqs_raw,
                                                  **{seq_id: seq for seq_id, seq in dat_seqs.items()
                                                     if seq_id in ref_names}},
                                          gap_char = GAP_CHAR, write = False, ref_names = ref_names,
                                          in_all_ref = in_all_ref)
            seqs = [seq for seq_id, seq in seqs.items() if "Col-0_ref" not in seq_id]
        else:
            seqs = remove_empty_pos(list(seqs_raw.values()), empty_char = GAP_CHAR)
    
        ## calculate pi while accounting for every pairwise alignment length
        n_seqs = len(seqs)
        freq_constant = 1/(n_seqs * n_seqs)
        pi_raw[gene] = 0
        pi_pairAlnLen[gene] = 0
        for i in range(n_seqs - 1):
            for j in range(i + 1, n_seqs):
                seq_i, seq_j = ( remove_incomplete_pos([seqs[i], seqs[j]], empty_char = GAP_CHAR) if IGNORE_GAPS
                                 else [seqs[i], seqs[j]] )
                if len(seq_i) == 0:
                    continue
                curr_pi = ( freq_constant *
                            len( tuple( 1 for c_i in range(len(seq_i))
                                        if ( seq_i[c_i] != seq_j[c_i] ) ) ) )
                pi_raw[gene] += curr_pi
                pi_pairAlnLen[gene] += curr_pi/len(seq_i)
        
        ## throw in the constants
        pi_raw[gene] *= 2
        pi_totalAlnLen[gene] = pi_raw[gene]/len(seqs[0])
        pi_pairAlnLen[gene] *= 2
    
    pi_header = ["gene", "cluster", "pi", "pi/alnlen", "pi/pairalnlen"]
    pi_get = make_custom_get(pi_header)
    genes = set(get(ids, "gene"))
    pi_dat = [gene_cluster[gene] + (pi_raw[gene], pi_totalAlnLen[gene], pi_pairAlnLen[gene])
              for gene in pi_raw]
    
    ## print different things?
    if preview:
        for cluster in set(get(ids, "cluster")):
            entries = sorted(entry for entry in pi_dat
                             if pi_get(entry, "cluster") == cluster)
            for entry in entries:
                print('\t'.join(str(x) for x in entry))
        
        for entry in sorted(pi_dat, key = lambda x: pi_get(x, "pi/alnlen")):
            print('\t'.join(str(x) for x in entry))
    
    ## write
    with open(f"{fout_dir}/{prefix}_{domain}_{'gapsExc' if IGNORE_GAPS else 'gapsInc'}_{'CDSonly' if cds_only else 'CDScomplete'}_pi.tsv", "w+") as f:
        f.write('\n'.join('\t'.join(str(y) for y in x)
                          for x in [pi_header] +
                          sorted(pi_dat, key = lambda x: pi_get(x, "gene"))))
    return

# calculate_pi_multi()
