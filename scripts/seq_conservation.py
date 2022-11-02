import re
import os
import Bio
import sys
import copy
import dendropy
import itertools

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from Bio import Phylo
from itertools import groupby
from dendropy.calculate import popgenstat
from fasta_manip import fasta_to_dict, dict_to_fasta, remove_empty_pos, remove_incomplete_pos
from data_manip import *

IGNORE_GAPS = True
GAP_CHAR = '-'

# dstats, pgstats, et, em = popgenome_pi("NB-ARC", cds_only = True, fout_dir = "/mnt/chaelab/rachelle/hap_tags/popg_test", keep_fasta = True)
dstats, pgstats, et, em = popgenome_pi("NB-ARC", cds_only = True, fout_dir = "/mnt/chaelab/rachelle/hap_tags/popg_test", keep_fasta = True, models = ["M0", "M1", "M2", "M3"])
dstats_mk, pgstats_mk = popgenome_pi("NB-ARC", cds_only = True, fout_dir = "/mnt/chaelab/rachelle/hap_tags/popg_test", keep_fasta = True, models = []) ## set models to empty list if not running dN/dS analysis
# dstats, pgstats, et, em = popgenome_pi("NB-ARC", cds_only = True, fout_dir = "/mnt/chaelab/rachelle/hap_tags/popg_test", keep_fasta = True, models = ["M3"])
dstats_mk, pgstats_mk = popgenome_pi("NB-ARC", cds_only = False, fill_empty = True, fout_dir = "/mnt/chaelab/rachelle/hap_tags/popg_test", keep_fasta = True, models = []) ## set models to empty list if not running dN/dS analysis


## write output
def write_pgstats_mk(stats, pgstats, fout_dir, prefix, domain, cds_only, intron_only):
    ## flatten mkstats
    def flatten_mk(mkstats):
        sn_fp = ["syn_fix", "nonsyn_fix", "syn_poly", "nonsyn_poly"]
        snfp = ["syn", "non", "fix", "poly"]
        header = ["mk"] + ['_'.join(["site_counts", x]) for x in sn_fp] + \
                 ["tot"] + ['_'.join(["tot", x]) for x in snfp] + \
                 ['_'.join(["exp", x]) for x in sn_fp]
        if isinstance(mkstats, str):
            return header, [mkstats] * 14
        else:
            output = ["NA" if not mkstats["p"] else mkstats["p"]]
            output += [mkstats["site_counts"][x] for x in sn_fp]
            output += ([mkstats["tot"]["tot"]] + [mkstats["tot"]['_'.join(["tot", x])] for x in snfp])
            output += (["NA"] * 4) if not mkstats["exp"] else \
                      [mkstats["exp"][x] for x in sn_fp]
            return header, output
    ## headers
    stats_col = ["gene", "numSeqRaw", "lenSeqTrim", "uniqSeqTrim", "numSeqTrim"]
    pgstats_col = ["nucleotide_diversity", "tajimas_d", "wattersons_theta"]
    mkstats_col = list(list(flatten_mk(dat["mktest"])[0] for dat in pgstats.values())[0])
    ## stuff to write out
    to_write = [["domainId"] + stats_col + pgstats_col + mkstats_col] + \
               [[domainId] + \
                ["NA" if not x in stats[domainId] else stats[domainId][x]
                 for x in stats_col] + \
                (["NA", "NA", "NA"] if (("lenSeqTrim" not in stats[domainId]) or
                                        stats[domainId]["lenSeqTrim"] == "NA" or
                                        stats[domainId]["uniqSeqTrim"] <= 1) else
                 [pgstats[domainId][x] for x in pgstats_col]) + \
                (flatten_mk("NA" if not domainId in pgstats else
                            pgstats[domainId]["mktest"])[1]) for
                domainId n sorted(list(stats.keys()))]
    ## open file and write .tsv
    # with open(f"{fout_dir}/{prefix}_{domain}_{'CDSonly' if cds_only else 'intrononly' if intron_only else 'CDScomplete'}_stats_pg_mk.tsv", "w+") as f:
    with open(f"{fout_dir}/{prefix}_{domain}_{'CDSonly' if cds_only else 'intrononly' if intron_only else 'CDScomplete'}_stats.tsv", "w+") as f:
        f.write('\n'.join('\t'.join(str(y) for y in x) for x in to_write))
    return

# ## write :)
# write_pgstats_mk(dstats_mk, pgstats_mk, "/mnt/chaelab/rachelle/hap_tags/popg_test",
#                  "nlr164", "NB-ARC", cds_only = True, intron_only = False)

def make_combo(main, *args):
    max_len = len(main)
    output = [main]
    output += [max_len * [x] for x in args]
    return output

def calc_dnds(domainId, fout_dir, models, prefix, domain, cds_only, intron_only):
    print(domainId, models)
    f_fasta = f"{fout_dir}/{prefix}_{domain}_{'CDSonly' if cds_only else 'intrononly' if intron_only else 'CDScomplete'}_{domainId}.fasta"
    x = fasta_to_dict(f_fasta)
    if len(x) <= 1:
        print("Sequences for {domainId} <= 1")
        return {domainId: None}, {f"none.{domainId}": None}
    from ete3 import EvolTree
    t = EvolTree(newick = f"{fout_dir}/{prefix}_{domain}_{'CDSonly' if cds_only else 'intrononly' if intron_only else 'CDScomplete'}_{domainId}.nwk")
    t.link_to_alignment(f_fasta)
    t.workdir = fout_dir
    evol_models = {}
    for model in models:
        model_name = f"{model}.{domainId}"
        t.run_model(model_name)
        evol_models[model_name] = t.get_evol_model(model_name)
    print(f"{domainId} completed")
    return {domainId: t}, evol_models

## trim sequences if required
def trim_intron_exon(seqs_raw, ref_seqs, cds_only = False, intron_only = False,
                     gap_char = GAP_CHAR, in_all_ref = False):
    ## if no trimming required
    if not cds_only and not intron_only:
        seqs = seqs_raw
    ## if references provided
    elif ref_seqs:
        kwargs = {"seqs": {**seqs_raw, **ref_seqs}, "ref_names": set(ref_seqs.keys()),
                 "gap_char": GAP_CHAR, "write": False, "in_all_ref": in_all_ref}
        if cds_only:
            from fasta_manip import trim_alignment_to_seqs
            seqs = trim_alignment_to_seqs(**kwargs)
        elif intron_only:
            from fasta_manip import trim_alignment_to_seqs_inv
            seqs = trim_alignment_to_seqs_inv(**kwargs)
        ## remove reference sequence
        seqs = {seq_id: seq for seq_id, seq in seqs.items() if "Col-0_ref" not in seq_id}
    ## if no reference provided, return empty dictionary
    else:
        seqs = dict()
    ## return final sequences
    return seqs

def popgenome_pi(domain, GAP_CHAR = '-', preview = False, cds_only = False, in_all_ref = False,
                 prefix = "nlr164_AL70", fout_dir = '', intron_only = False, keep_fasta = False,
                 tree_exists = False, models = ["M0"], fill_empty = False,
                 make_seq_fname = lambda domain: f"/mnt/chaelab/rachelle/hap_tags/results/nlr164/alignment/nlr164_col0-AL70-Alyrata_{domain}_mafft.fa",
                 make_pred_fname = lambda domain: f"/mnt/chaelab/rachelle/hap_tags/results/nlr164/predicted_identity/nlr164_AL70_{domain}_predictedIdentity.txt",
                 make_tree_fname = lambda domain: f"/mnt/chaelab/rachelle/hap_tags/results/nlr164/tree/nlr164_col0-AL70-Alyrata_{domain}_mafft_ML.nwk",
                 make_pred_araly_fname = lambda domain: f"/mnt/chaelab/rachelle/hap_tags/results/nlr164/predicted_identity/nlr164_Alyrata_{domain}_predictedIdentity.txt"):
    
    dat_seqs = fasta_to_dict(make_seq_fname(domain))
    with open(make_pred_fname(domain), 'r') as f:
        id_dat = [x[:-1].split('\t') for x in f.readlines()]
    with open(make_pred_araly_fname(domain), 'r') as f:
        id_dat_araly = [x[:-1].split('\t') for x in f.readlines()]
    tree = Phylo.read(make_tree_fname(domain), "newick")
    
    get = make_custom_get(id_dat[0])
    ## organise A. thaliana data
    ids = id_dat[1:]
    domain_cluster = set(list(itertools.chain.from_iterable([list(zip(get(x, "domain").split(';'), get(x, "cluster").split(';'))) for x in ids])))
    domain_cluster = {x[0]: x for x in domain_cluster}
    domains = set(domain_cluster.keys())
    # ## organise A. lyrata data
    ids_araly = id_dat_araly[1:]
    # domain_cluster_araly = set(list(itertools.chain.from_iterable([list(zip(get(x, "domain").split(';'), get(x, "cluster").split(';'))) for x in ids_araly])))
    # domain_cluster_araly = {x[0]: x for x in domain_cluster_araly}
    
    ## dictionary to summarise each domain (whether it has valid sequences or not, and how many there are)
    stats = dict()
    pgstats = {} ## dict for stats
    fouts = []
    
    for i, domainId in enumerate(domains):
        print(i + 1, domainId)
        gene, domainNum = domainId.split('.')
        
        ## get sequences assigned to the gene
        seq_ids = set(get(entry, "contig") for entry in ids + ids_araly
                      if domainId in get(entry, "domain"))
        seqs_raw = {seq_id: seq for seq_id, seq in dat_seqs.items() if seq_id in seq_ids}
        
        ## trim sequences if required
        ref_names = tuple(seq_id for seq_id in dat_seqs.keys() if
                                  re.search("Col-0_ref\|" + gene + "(\|[^|]+){3}\|" + str(domainNum), seq_id)
                                  and "complete" not in seq_id)
        ref_seqs = {seq_id: dat_seqs[seq_id] for seq_id in ref_names}
        seqs_all = trim_intron_exon(seqs_raw, ref_seqs, cds_only = cds_only, intron_only = intron_only,
                                    gap_char = GAP_CHAR, in_all_ref = in_all_ref)
        ## write some general stats
        stats[domainId] = {"gene": gene,
                           "numSeqRaw": len([seq_id for seq_id in seqs_all if "tig" in seq_id]),
                           "numSeqRawAraly": len([seq_id for seq_id in seqs_all if "tig" not in seq_id])}
        
        ## if no A. thaliana sequences, skip
        if not stats[domainId]["numSeqRaw"]:
            continue
        
        ## remove all empty positions and sequences and update stats
        seqs_all = {seq_id: seq.upper()
                    for seq_id, seq in remove_empty_pos(seqs_all, empty_char = GAP_CHAR).items()
                    if len(seq) > 0}
        seqs = {seq_id: seq for seq_id, seq in seqs_all.items() if "tig" in seq_id}
        stats[domainId]["numSeqTrim"] = len(seqs)
        stats[domainId]["numSeqTrimAraly"] = len(seqs_all) - len(seqs)
        ## the following stats are only for A. thaliana
        if len(seqs) == 0:
            stats[domainId]["uniqSeqTrim"] = "NA"
            stats[domainId]["lenSeqTrim"] = "NA"
        else:
            stats[domainId]["uniqSeqTrim"] = len(set(str(seq) for seq in seqs.values()))
            stats[domainId]["lenSeqTrim"] = len(list(seqs.values())[0])
        
        # ## write to file so we can pass it into R later (or, uh, use dendropy's popgenstat)
        fasta_out = f"{fout_dir}/{prefix}_{domain}_{'CDSonly' if cds_only else 'intrononly' if intron_only else 'CDScomplete'}_{domainId}.fasta"
        # fouts.append(fasta_out)
        # dict_to_fasta(seqs, fasta_out)
        
        ## popgenstat
        if "NA" != stats[domainId]["uniqSeqTrim"] > 1:
            cmat = dendropy.DnaCharacterMatrix.get(file = open(fasta_out), schema = "fasta")
            pgstats[domainId] = dict()
            try:
                pgstats[domainId]["nucleotide_diversity"] = popgenstat.nucleotide_diversity(cmat)
            except ZeroDivisionError:
                pgstats[domainId]["nucleotide_diversity"] = "NA"
            try:
                pgstats[domainId]["tajimas_d"] = popgenstat.tajimas_d(cmat)
            except ZeroDivisionError:
                pgstats[domainId]["tajimas_d"] = "NA"
            try:
                pgstats[domainId]["wattersons_theta"] = popgenstat.wattersons_theta(cmat)
            except ZeroDivisionError:
                pgstats[domainId]["wattersons_theta"] = "NA"
            ## try McDonald-Kreitman if Araly sequences exist
            if stats[domainId]["numSeqTrimAraly"] >= 1:
                from Bio.codonalign.codonalignment import CodonSeq, default_codon_alphabet
                from Bio.Alphabet import generic_dna
                from Bio.SeqRecord import SeqRecord
                from Bio.Alphabet import IUPAC, Gapped
                from fasta_manip import remove_incomplete_pos, unknown_codon_to_gap
                from codonalignment_modified import CodonAlignment, mktest ## modified _G_test to return stats
                ungapped_codons = {seq_id: seq for seq_id, seq
                                   in remove_incomplete_pos(unknown_codon_to_gap(seqs_all, fill_empty = fill_empty)).items()
                                   if len(seq.replace('-', '')) > 0}
                ungapped_codons_thaliana = [SeqRecord(CodonSeq(unknown_codon_to_gap(seq, fill_empty = fill_empty)), id = seq_id)
                                            for seq_id, seq in ungapped_codons.items()
                                            if seq_id in seqs]
                ungapped_codons_araly = [SeqRecord(CodonSeq(unknown_codon_to_gap(seq, fill_empty = fill_empty)), id = seq_id)
                                         for seq_id, seq in ungapped_codons.items()
                                         if not seq_id in seqs]
                if ungapped_codons_thaliana and ungapped_codons_araly:
                    # try:
                        codonseqthaliana = CodonAlignment(ungapped_codons_thaliana)
                        codonseqaraly = CodonAlignment(ungapped_codons_araly)
                        pgstats[domainId]["mktest"] = mktest([codonseqthaliana, codonseqaraly])
                    # except ValueError:
                    #     pgstats[domainId]["mktest"] = "mktest_ValueError"
                    # except ZeroDivisionError:
                    #     pgstats[domainId]["mktest"] = "mktest_ZeroDivisionError"
                else:
                    pgstats[domainId]["mktest"] = "no_seq"
            else:
                pgstats[domainId]["mktest"] = "NA"
        
        ## get relevant subtree/branches + do some other calculations
        if stats[domainId]["lenSeqTrim"] != "NA":
            ## count number of sites or codons where >= 1 base differs across all sequences
            stats[domainId]["nSegCodons"] = sum([1 if len(set(seq[i:i+3] for seq in seqs.values())) > 1 else 0
                                                   for i in range(0, stats[domainId]["lenSeqTrim"], 3)])
            stats[domainId]["nSegSites"] = sum([1 if len(set(seq[i] for seq in seqs.values())) > 1 else 0
                                                  for i in range(0, stats[domainId]["lenSeqTrim"])])
            # ## get subtree/branches
            # mrca = copy.deepcopy(tree.common_ancestor(tree.find_elements(name = '|'.join(seq_id.replace('|', '\|')for seq_id in seqs.keys()))))
            # to_prune = set(x.name for x in mrca.get_terminals()) - set(seqs.keys())
            # for leaf in to_prune:
            #     x = mrca.prune(leaf)
            # Phylo.write(mrca, f"{fout_dir}/{prefix}_{domain}_{'CDSonly' if cds_only else 'intrononly' if intron_only else 'CDScomplete'}_{domainId}.nwk", "newick")
        
    # ## pass into R or do something with pgstats??
    
    # ## write pgstats
    # stats_col = ["gene", "numSeqRaw", "lenSeqTrim", "uniqSeqTrim", "numSeqTrim"]
    # pgstats_col = ["nucleotide_diversity", "tajimas_d", "wattersons_theta"]
    # to_write = [["domainId"] + stats_col + pgstats_col] + \
    #            [[domainId] + [stats[domainId][x] for x in stats_col] +
    #             (["NA", "NA", "NA"] if (stats[domainId]["lenSeqTrim"] == "NA" or
    #                                     stats[domainId]["uniqSeqTrim"] <= 1) else
    #              [pgstats[domainId][x] for x in pgstats_col]) for
    #             domainId in sorted(list(stats.keys()))]
    # with open(f"{fout_dir}/{prefix}_{domain}_{'CDSonly' if cds_only else 'intrononly' if intron_only else 'CDScomplete'}_stats.tsv", "w+") as f:
    #     f.write('\n'.join('\t'.join(str(y) for y in x) for x in to_write))
    
    # ## calculate dN/dS for each domain
    # if models:
    #     evol_trees = {}
    #     evol_models = {}
    #     from concurrent import futures
    #     with futures.ProcessPoolExecutor() as pool:
    #         for et, em in pool.map(calc_dnds,
    #                                *make_combo([domainId for domainId in stats.keys()
    #                                             if "NA" != stats[domainId]["uniqSeqTrim"] > 1],
    #                                            fout_dir, models, prefix, domain, cds_only, intron_only)):
    #             # for et, em in pool.map(calc_dnds, domains, [fout_dir]*ndomains, [models]*ndomains, [prefix]*ndomains, [domain]*ndomains, [cds_only]*ndomains, [intron_only]*ndomains):
    #             evol_trees = {**evol_trees, **et}
    #             evol_models = {**evol_models, **em}
    
    # ## delete fasta files if "keep_fasta" is not raised
    # if not keep_fasta:
    #     for fasta_out in fouts:
    #         os.remove(fasta_out)
    
    return stats, pgstats#, evol_trees, evol_models
    # return stats, pgstats



def calculate_pi_multi(domains = ("TIR", "NB-ARC"), **kwargs):
    for d in domains:
        for ignore_gaps in (True, False):
            for cds_only in (True, False):
                print(f"Ignore gaps: {ignore_gaps} || CDS-only: {cds_only} || Intron-only: False")
                calculate_pi(domain = d, IGNORE_GAPS = ignore_gaps, cds_only = cds_only,
                             intron_only = False, **kwargs)
            print(f"Ignore gaps: {ignore_gaps} || CDS-only: False || Intron-only: True")
            calculate_pi(domain = d, IGNORE_GAPS = ignore_gaps, cds_only = False, intron_only = True, **kwargs)
    return

def calculate_pi(domain, IGNORE_GAPS, GAP_CHAR = '-', preview = False, cds_only = False, in_all_ref = False,
                 prefix = "nlr164_AL70", fout_dir = '', intron_only = False,
                 make_seq_fname = lambda domain: f"/mnt/chaelab/rachelle/hap_tags/results/nlr164/alignment/nlr164_col0-AL70-Alyrata_{domain}_mafft.fa",
                 make_pred_fname = lambda domain: f"/mnt/chaelab/rachelle/hap_tags/results/nlr164/predicted_identity/nlr164_AL70_{domain}_predictedIdentity.txt"):
    
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
        elif intron_only:
            ref_names = tuple(seq_id for seq_id in dat_seqs.keys() if
                              f"Col-0_ref|{gene}" in seq_id and "complete" not in seq_id)
            if not ref_names:
                continue
            from fasta_manip import trim_alignment_to_seqs
            seqs = trim_alignment_to_seqs_inv(seqs = {**seqs_raw,
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
    with open(f"{fout_dir}/{prefix}_{domain}_{'gapsExc' if IGNORE_GAPS else 'gapsInc'}_{'CDSonly' if cds_only else 'intrononly' if intron_only else 'CDScomplete'}_pi.tsv", "w+") as f:
        f.write('\n'.join('\t'.join(str(y) for y in x)
                          for x in [pi_header] +
                          sorted(pi_dat, key = lambda x: pi_get(x, "gene"))))
    return

# calculate_pi_multi()
