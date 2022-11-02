import os
import re
import sys
import copy
import numpy
import scipy
import dendropy
import itertools

dir_results = '' ## replace with path to results directory
dir_scripts = '' ## replace with path to scripts directory
sys.path.append(dir_scripts)

from basics import *
from data_manip import *
from fasta_manip import fasta_to_dict, dict_to_fasta
from Bio import Phylo
from dendropy.calculate import popgenstat

## TODO:
# - get number of sequences in each clade

## parse data
t_nwk = f"{dir_results}/tree/nlr164_col0-AL70-Alyrata_NB-ARC_mafft_ML.nwk"
predicted_identity_f = f"{dir_results}/predicted_identity/nlr164_AL70_NB-ARC_predictedIdentity.txt"
## some sequences to exclude if so desired
to_exclude = ["5993.tig00000076|5993|NB-ARC|716|4515-5174|revcomp"]
## get root and gene clades (grouped by mono- and polyphyly)
root, gene_clades = get_gene_clades(t_nwk, "placeholder", predicted_identity_f, to_exclude = to_exclude,
                                    exclude_ref = False)

fout_dir = f"{dir_results}/popg_test/clades"
f_aln_raw = f"{dir_results}/alignment/nlr164_col0-AL70-Alyrata_NB-ARC_mafft.fa"

## these are confidence values, starting from closer to the root to further from the root
## (take, these, steps, sequentially)
## [take, these, at, the, same, level]
cluster_groups = {"DM2": [(0.998, 0.448, 0.526, 0.658, 0.922, 1.00, 0.531, 1.000, 1.000, 0.838, 0.944,
                           0.800, 0.927, 0.653, 0.998, 0.999, [0.898, 1.000])],
                  "DM4": [(0.998, 0.922, 0.890, 1.00, 0.989, 0.973,
                           0.911, 1.000, [(0.299, [1.000, 1.000]), 1.000])],
                  "DM8": [(0.998, 0.448, 0.526, 0.658, 0.922, 1.00, 0.531, 0.581, 0.999,
                           1.000, [(0.804, 0.992), (0.945, 0.919, 0.996, [0.871, 1.000])])]}
radiation_groups = {"DM2": [(0.998, 0.448, 0.526, 0.658, 0.922, 1.00, 0.531, 1.000, 1.000, 0.838, 0.944,
                             0.800, 0.927, 0.653, 0.998, 0.999, 0.898)],
                    "DM4": [(0.998, 0.922, 0.890, 1.00, 0.989, 0.973,
                             0.911, 1.000, 0.299, 1.000, {0.831})],
                    "DM8": [(0.998, 0.448, 0.526, 0.658, 0.922, 1.00, 0.531, 0.581, 0.999,
                             1.000, 0.945, 0.919, 0.996, 0.871)],
                    "B3": [(0.998, 0.922, 0.999, 0.813, 0.585, 0.986, 0.857, 0.836,
                            0.971, 0.992, 0.947, 1.000)]}
hifi_groups = {"DM2": [(0.998, 0.448, 0.526, 0.658, 0.922, 1.00, 0.531, 1.000, 1.000, 0.838, 0.944,
                        0.800, 0.927, 0.653, 0.998, 0.999, 1.000)],
               "DM4": [(0.998, 0.922, 0.890, 1.00, 0.989, 0.973,
                        0.911, 1.000, [(0.299, 1.000, {0.365}), 1.000])],
               "DM8": [(0.998, 0.448, 0.526, 0.658, 0.922, 1.00, 0.531, 0.581, 0.999,
                        1.000, [(0.804, 0.992), (0.945, [(0.255, 1.000, 0.00), (0.919, 0.996, 1.000)])])]}
hifi_multi_groups = {"RPP13": [(0.998, 0.922, 0.890, 1.000, 0.996)],
                     "B5_around": [(0.998, 0.448, 0.526, 0.658, 0.922, 0.960, 0.631)]}

## get clades in cluster (manually curated)
def get_clades_from_conf_seqs_dict(conf_seqs, root = root):
    return {k: get_clades_from_confidence_seq_multi(root, seqs) for k, seqs in conf_seqs.items()}
cluster_clades = get_clades_from_conf_seqs_dict(cluster_groups, root)
radiation_clades = get_clades_from_conf_seqs_dict(radiation_groups, root)
hifi_clades = get_clades_from_conf_seqs_dict(hifi_groups, root)
hifi_multi_clades = get_clades_from_conf_seqs_dict(hifi_multi_groups, root)

pgstats, et, em = calculate_popg_multi({"cluster": cluster_clades, "radiation": radiation_clades,
                                        "hifi": hifi_clades, "hifi_multi": hifi_multi_clades},
                                       fout_dir = f"{dir_results}/popg_test/clades",
                                       f_aln_raw = f"{dir_results}/alignment/nlr164_col0-AL70-Alyrata_NB-ARC_mafft.fa")

pgstats = calculate_popg_multi({"cluster": cluster_clades, "radiation": radiation_clades,
                                "hifi": hifi_clades, "hifi_multi": hifi_multi_clades},
                               fout_dir = f"{dir_results}/popg_test/clades",
                               f_aln_raw = f"{dir_results}/alignment/nlr164_col0-AL70-Alyrata_NB-ARC_mafft.fa")
f = open(f"{dir_results}/popg_test/clades/nlr164_NB-ARC_CDSonly_clade_stats.tsv", "w+")
header_stats = ["nucleotide_diversity", "tajimas_d", "wattersons_theta"]
header = ["type", "group", "num", "tgid"] + header_stats
f.write('\n'.join(['\t'.join([str(x) for x in line])
                   for line in [header] + \
                   [tgid.split(',') + [tgid] + [tgid_stats[stat] for stat in header_stats]
                    for tgid, tgid_stats in pgstats.items()]]))
f.close()

## whole sequence (without distinguishing CDS from intron)
pgstats = calculate_popg_multi({"cluster": cluster_clades, "radiation": radiation_clades,
                                "hifi": hifi_clades, "hifi_multi": hifi_multi_clades},
                               fout_dir = f"{dir_results}/popg_test/clades",
                               f_aln_raw = f"{dir_results}/alignment/nlr164_col0-AL70-Alyrata_NB-ARC_mafft.fa",
                               cds_only = False, models = [])
f = open(f"{dir_results}/popg_test/clades/nlr164_NB-ARC_CDScomplete_clade_stats.tsv", "w+")
header_stats = ["nucleotide_diversity", "tajimas_d", "wattersons_theta"]
header = ["type", "group", "num", "tgid"] + header_stats
f.write('\n'.join(['\t'.join([str(x) for x in line])
                   for line in [header] + \
                   [tgid.split(',') + [tgid] + [tgid_stats[stat] for stat in header_stats]
                    for tgid, tgid_stats in pgstats[0].items()]]))
f.close()


# ## get branch lengths within each cluster clade
# cluster_branch_lengths = {cluster: {clade: get_branch_lengths(clade) for clade in clades}
#                           for cluster, clades in cluster_clades.items()}
## calculate some branch length stats
cluster_branch_length_stats = {cluster: {clade: get_clade_branch_length_stats(clade = clade) for clade in clades}
                               for cluster, clades in cluster_clades.items()}
## calculate branch length stats for mono- and polyphyletic clades
gene_branch_length_stats = {group: {clade: get_clade_branch_length_stats(clade = clade)
                                    for clade, genes in clade_genes}
                            for group, clade_genes in gene_clades.items()}
## print branch length stats
print_branch_length_stats(cluster_branch_length_stats)
print_branch_length_stats(gene_branch_length_stats)

## write stats to file
stats_to_include = ("mean", "sd", "median", "max", "min")
fout = f"{dir_results}/branch_lengths/clade_stats.tsv"
write_branch_length_stats_multi(fout, {"gene": gene_branch_length_stats,
                                       "cluster": cluster_branch_length_stats},
                                stats_to_include = stats_to_include)

# try:
#     x = split_and_calculate_bl_stats(cluster_clades["DM2"][0])
# except RecursionError:
#     print("oops")
# # sys.setrecursionlimit(5000)

stats_to_include = ("mean", "sd", "median", "max", "min")

## calculate some stats while splitting clades by max branch length and write to file
fout_bl_split = f"{dir_results}/branch_lengths/clade_split_longest_stats_mmsd2.tsv"
fout_bl_split_tmp = f"{dir_results}/branch_lengths/tmp.tsv"
layer_names = ["type", "group"]
split_calculate_write_bl_stats_multi(fout_bl_split,
                                     {"gene": {k: [x[0] for x in v] for k,v in gene_clades.items()},
                                      "cluster": cluster_clades, "radiation": radiation_clades,
                                      "hifi": hifi_clades, "hifi_multi": hifi_multi_clades},
                                     layer_names, monophyly_only = False, max_depth = 40,
                                     max_mean_sd_diff = 1)

# ## function to flatten branch length stats output from split_and_calculate_bl_stats
# def flatten_split_bl_stats(master_split_bl_stats, path = [],
#                            stats_to_include = ("mean", "sd", "median", "max", "min", "branch_length",
#                                                "status", "size")):
#     def helper(split_bl_stats, path):
#         ## if terminal
#         if not split_bl_stats[1]:
#             return [ [','.join(str(p) for p in path)] + 
#                      [split_bl_stats[0].get(stat, None) for stat in stats_to_include] ]
#         ## else assign path to each feature
#         else:
#             output = helper( [split_bl_stats[0], []], path ) ## get current level's stats
#             ## get next levels' stats
#             for i, bl_stats in enumerate(split_bl_stats[1]):
#                 next_path = path + [i + 1]
#                 output.extend( flatten_split_bl_stats(bl_stats, next_path) )
#             return output
#     return helper(master_split_bl_stats, path)


## trim sequences if required
def trim_intron_exon(seqs_raw, ref_seqs, cds_only = False, intron_only = False,
                     gap_char = '-', in_all_ref = False):
    ## if no trimming required
    if not cds_only and not intron_only:
        seqs = seqs_raw
    ## if references provided
    elif ref_seqs:
        kwargs = {"seqs": {**seqs_raw, **ref_seqs}, "ref_names": set(ref_seqs.keys()),
                 "gap_char": gap_char, "write": False, "in_all_ref": in_all_ref}
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

## tag job with name
def tagged_job(f):
    def output(n, *args, **kwargs):
        return [n, f(*args, **kwargs)]
    return output

## run model packaged into convenient function for parallelisation
def run_get_model(t, m):
    t.run_model(m)
    return t, t.get_evol_model(m)

def tagged_run_get_model(n, t, m):
    return tagged_job(run_get_model)(n, t, m)

## calculate dN/dS
def calc_dnds(models, job_name, f_fasta, f_nwk):
    print(job_name, models)
    x = fasta_to_dict(f_fasta)
    ## return empty dictionaries if insufficient sequences
    if len(x) <= 1:
        print("Sequences for {job_name} <= 1")
        return dict(), dict()
    ## make EvolTree object
    from ete3 import EvolTree
    t = EvolTree(newick = f_nwk)
    t.link_to_alignment(f_fasta)
    t.workdir = fout_dir
    ## parallelise model running
    evol_trees = {}
    evol_models = {}
    from concurrent import futures
    model_names = [f"{model}.{job_name}" for model in models]
    with futures.ProcessPoolExecutor() as pool:
        for tag, output in pool.map(tagged_run_get_model, model_names, [t]*len(model_names), model_names):
            evol_trees[tag] = output[0]
            evol_models[tag] = output[1]
    print(f"{job_name} completed")
    return evol_trees, evol_models

## tree: Biopython tree (clade) object
## aln: dictionary of sequences {seq_id: seq}
## to_keep, to_exclude: iterable of seq_id of sequences to keep or exclude respectively
def calculate_popg(clade_name, tree, aln, fout_dir, to_keep = [], to_exclude = [], tmp_pref = 'tmp',
                   models = ["M0"], cds_only = True, intron_only = False):
    ## get the desired sequences
    to_keep = set(to_keep) if to_keep else set(aln.keys())
    to_exclude = set(to_exclude) if to_exclude else set()
    seqs_raw = {seq_id: seq for seq_id, seq in aln.items() if seq_id in to_keep and not seq_id in to_exclude}
    ## trim and remove reference from alignment
    seqs = trim_intron_exon({seq_id: seq for seq_id, seq in seqs_raw.items() if not "Col-0" in seq_id},
                            {seq_id: seq for seq_id, seq in seqs_raw.items()
                             if "Col-0" in seq_id and "complete" not in seq_id},
                            cds_only = cds_only, intron_only = intron_only,
                            gap_char = '-', in_all_ref = False)
    output = dict()
    if len(seqs) < 1:
        output["nucleotide_diversity"] = "no_seqs"
        output["tajimas_d"] = "no_seqs"
        output["wattersons_theta"] = "no_seqs"
        em, et = {}, {}
    else:
        ## write to file so that ete3 can work with it
        fasta_tmp = f"{fout_dir}/{tmp_pref}.fasta"
        nwk_tmp = f"{fout_dir}/{tmp_pref}.nwk"
        dict_to_fasta(seqs, fasta_tmp)
        Phylo.write(tree, nwk_tmp, "newick")
        ## calculate stats
        cmat = dendropy.DnaCharacterMatrix.get(file = open(fasta_tmp), schema = "fasta")
        output = dict()
        try:
            output["nucleotide_diversity"] = popgenstat.nucleotide_diversity(cmat)
        except ZeroDivisionError:
            output["nucleotide_diversity"] = "NA"
        try:
            output["tajimas_d"] = popgenstat.tajimas_d(cmat)
        except ZeroDivisionError:
            output["tajimas_d"] = "NA"
        try:
            output["wattersons_theta"] = popgenstat.wattersons_theta(cmat)
        except ZeroDivisionError:
            output["wattersons_theta"] = "NA"
        # ## calculate model stuff
        # from concurrent import futures
        # nmodels = len(models)
        # em, et = calc_dnds(models, clade_name, fasta_tmp, nwk_tmp)
        # print("calculate_popg", {clade_name: em}, {clade_name: et})
    return output, {}, {}# , {clade_name: em}, {clade_name: et}

## takes dictionary of clades (number of levels to reach clade must be same for all clades), output directory, and alignment fasta file containing all relevant sequences
def calculate_popg_multi(clades_d, fout_dir, f_aln_raw, **kwargs):    
    def helper(d, layers):
        dnds_et, dnds_em, pgstats = dict(), dict(), dict()
        ## if not terminal layer, call helper recursively for each dictionary
        if set(type(x) for x in d.values()) == {dict}:
            for k, v in d.items():
                tmp_et, tmp_em, tmp_pgstats = helper(v, layers = [k])
                dnds_et = {**dnds_et, **tmp_et}
                dnds_em = {**dnds_em, **tmp_em}
                pgstats = {**pgstats, **tmp_pgstats}
        ## if terminal layer, calculate stats and flatten them
        else:
            for k, clades in d.items():
                for i, clade in enumerate(clades):
                    job_name = ','.join(layers + [k, str(i + 1)])
                    print(f"Processing {job_name}")
                    ## get split by branch length stats for the clade
                    from fasta_manip import remove_empty_pos
                    leaves = set(leaf.name for leaf in clade.get_terminals() if not "Col-0" in leaf.name)
                    ref_seqs_pattern = set(''.join((re.search("^.+?\|CDS\|", leaf.name).group(0), "[^|]+",
                                                    re.search("\|NB-ARC\|\d", leaf.name).group(0))).replace('|',
                                                                                                            '\|')
                                           for leaf in clade.get_terminals()
                                           if "Col-0" in leaf.name)
                    seqs = {seq_id: seq for seq_id, seq in fasta_to_dict(f_aln_raw).items() if seq_id in leaves}
                    ref_seqs = {seq_id: seq for seq_id, seq in fasta_to_dict(f_aln_raw).items()
                                if True in [bool(re.search(pattern, seq_id)) for pattern in ref_seqs_pattern]}
                    seqs = remove_empty_pos({**seqs, **ref_seqs})
                    ## check if a reference exists; if not, skip
                    if not ref_seqs:
                        print("No reference found. Skipping.")
                        continue
                    ## make a copy of the tree, remove the reference, and pass into calculate_popg
                    clade_copy = copy.deepcopy(clade)
                    new_root = prune_tree(clade_copy, to_keep = leaves)
                    tmp_et, tmp_em, tmp_pgstats = calculate_popg(job_name, clade_copy, seqs, fout_dir, **kwargs)
                    tmp_et = {job_name: tmp_et}
                    tmp_em = {job_name: tmp_em}
                    tmp_pgstats = {job_name: tmp_pgstats}
                    dnds_et = {**dnds_et, **tmp_et}
                    dnds_em = {**dnds_em, **tmp_em}
                    pgstats = {**pgstats, **tmp_pgstats}
        return dnds_et, dnds_em, pgstats
    return helper(clades_d, [])

def flatten_iterable(master_l):
    output = []
    def helper(l):
        for element in l:
            if not isinstance(element, list) and not isinstance(element, tuple):
                output.append(element)
            else:
                helper(element)
        return
    helper(master_l)
    return output

## quick stats for nested branch length stats (output of split_and_calculate_bl_stats)
def count_nested_layers(master_l, item_wrapper_is_list = True):
    max_depth = [1]
    num_splits = [0]
    def is_list_or_tuple(x):
        return isinstance(x, list) or isinstance(x, tuple)
    def helper(l, depth):
        ## check if there are more layers. If this is last layer, return
        if not True in [is_list_or_tuple(x) for x in l]:
            return
        ## check if items in layer are of the same type (i.e. result of a split). If yes, add to split count
        if is_list_or_tuple(l) and len(set(type(x) for x in l)) == 1:
            num_splits[0] += max(len(l) - 1, 0)
        ## iterate through all items. If is layer, run helper recursively
        for next_l in l:
            if isinstance(next_l, list) or isinstance(next_l, tuple):
                if depth + 1 > max_depth[0]:
                    max_depth[0] = depth + 1
                helper(next_l, depth + 1)
        return
    helper(master_l, 1)
    denom = (2 if item_wrapper_is_list else 1)
    return {"max_depth": max_depth[0] / denom, "num_splits": num_splits[0]}

# ## split by longest branch length and calculate stats after each split
# def split_and_calculate_bl_stats(clade, max_depth = 5, stats_last = {}, max_mean_sd_diff = 5,
#                                  stop_when = (lambda stats_last, stats_now:
#                                               stats_now["max"] < (stats_last["mean"] + 5 * stats_last["sd"])),
#                                  **kwargs):
#     def sd_stop_condition(stats_last, stats_now):
#         return stats_now["max"] < (stats_last["mean"] + max_mean_sd_diff * stats_last["sd"])
#     def helper(clade, depth, stats_last):
#         ## stop recursion if clade is terminal
#         if clade.is_terminal():
#             return [{**get_clade_branch_length_stats(branch_lengths = [clade.branch_length, clade.branch_length]),
#                      "size": 1, "status": "terminal"}, []]
#         else:
#             ## get branch lengths + branch length stats
#             branch_lengths = get_branch_lengths(clade)
#             clade_bl_stats = {**get_clade_branch_length_stats(branch_lengths = branch_lengths),
#                               "size": len(branch_lengths), "branch_length": clade.branch_length,
#                               "confidence": clade.confidence}
#             ## check termination conditions
#             if depth == 1:
#                 return [{**clade_bl_stats, "status": f"max_depth_reached;max_depth={max_depth}"}, []]
#             elif stats_last and sd_stop_condition(stats_last, clade_bl_stats):
#                 return [{**clade_bl_stats, "status": f"max_mean_sd_diff_fulfilled={max_mean_sd_diff}"}, []]
#             ## if no termination condition met, continue
#             else:
#                 clades_split = split_by_longest_branch(clade, **kwargs)
#                 return [{**clade_bl_stats, "status": "continued"},
#                         [helper(c, stats_last = clade_bl_stats, depth = depth - 1)
#                          for c in clades_split]]
#     return helper(clade, max_depth, {})

# ## function to split a clade by longest brach
# def split_by_longest_branch(clade, monophyly_only = True, **kwargs):
#     ## helper functions
#     def get_clade_with_longest_branch(c):
#         return max([x for x in c.get_terminals() + c.get_nonterminals() if not x is c],
#                    key = lambda subclade: subclade.branch_length)
#     def get_children_along_path(path):
#         ## note that path must be in this direction: [most basal, least basal]
#         return list(subclade for clade in path for subclade in clade.clades if not subclade in path)
#     ## get clade with longest branch
#     subclade = get_clade_with_longest_branch(clade)
#     if subclade in clade.clades: ## if subclade is direct descendent of clade
#         return list(clade.clades)
#     elif monophyly_only: ## break into multiple monophyletic clades if para-/polyphyly not allowed
#         path_to_clade = clade.get_path(subclade)[:-1] ## does not include clade and subclade
#         return list(c for c in clade.clades if not c is subclade) + \
#             list(get_children_along_path(path_to_clade)) + [subclade]
#     else: ## make deep copy and prune if para-/poly phyly allowed
#         clade_copy = copy.deepcopy(clade) ## make deep copy of clade -- this will be pruned
#         subclade = get_clade_with_longest_branch(clade_copy) ## re-identify clade w/ longest branch in copied tree
#         subclade_copy = copy.deepcopy(subclade) ## make deep copy of subclade
#         subclade_terminals = subclade.get_terminals() ## get terminals of original subclade to be pruned
#         for subclade_terminal in subclade_terminals: ## remove each terminal from the deep copy of clade
#             pruned = clade_copy.prune(subclade_terminal)
#         return [clade_copy, subclade_copy] ## return pruned para-/polyphyletic clade + longest branch clade

# ## write multiple branch length stats dictionaries to file
# def write_branch_length_stats_multi(fout, stats_d, stats_to_include = ("mean", "sd", "median", "max", "min")):
#     f = open(fout, "w+")
#     def join_stats(bl_stats, clade_type, cols = () ):
#         return '\n'.join(['\t'.join([clade_type, group, ','.join(clade.genes)] +
#                                     [str(stats[col]) for col in cols])
#                           for group, clade_stats in bl_stats.items() for clade, stats in clade_stats.items()])
#     f.write('\n'.join(['\t'.join(["type", "group", "genes"] + list(stats_to_include))] +
#                       [join_stats(bl_stats, bl_stats_name, cols = stats_to_include)
#                        for bl_stats_name, bl_stats in stats_d.items()]))
#                        join_stats(gene_branch_length_stats, "gene", cols = stats_to_include),
#                        join_stats(cluster_branch_length_stats, "cluster", cols = stats_to_include)]))
#     f.close()
#     return



# def print_branch_length_stats(all_bl_stats):
#     for cluster, clade_bl_stats in all_bl_stats.items():
#         print(f"--{cluster}--")
#         for clade, bl_stats in clade_bl_stats.items():
#             print(f"mean: {bl_stats['mean']}\tsd: {bl_stats['sd']}\tmedian: {bl_stats['median']}\tmax: {bl_stats['max']}\tmin: {bl_stats['min']}")
#     return

# ## calculate some common stats values
# def get_clade_branch_length_stats(clade = None, branch_lengths = [],
#                                   stats = ("mean", "sd", "median", "max", "min")):
#     if not branch_lengths and not clade:
#         return {}
#     elif not branch_lengths:
#         branch_lengths = get_branch_lengths(clade)
#     results = {"mean": numpy.mean(branch_lengths), "sd": numpy.std(branch_lengths),
#                "median": numpy.median(branch_lengths), "max": max(branch_lengths), "min": min(branch_lengths)}
#     stats = set(stats)
#     return {k: v for k, v in results.items() if k in results}

## flatten nested confidence value sequences
## tuples contain sequential nodes, lists contain simultaneous nodes (x-furcation, where you take both paths), and sets contain "look-ahead" nodes
def flatten_confidence_seqs(conf_seq):
    if isinstance(conf_seq, float) or isinstance(conf_seq, set):
        return [(conf_seq,)]
    elif len(conf_seq) == 0:
        return [()]
    elif isinstance(conf_seq[0], float) or isinstance(conf_seq[0], set):
        return [(conf_seq[0],) + x for x in flatten_confidence_seqs(conf_seq[1:])]
    elif isinstance(conf_seq[0], list):
        output = []
        for nested_conf_seq in conf_seq[0]:
            output.extend(flatten_confidence_seqs(nested_conf_seq))
        return output

def set_conf_seqs(seqs):
    output = []
    for seq in seqs:
        if not seq in output:
            output.append(seq)
    return output

## get clade from a single sequence of confidence values
def get_clades_from_confidence_seq(first_c, master_conf_seq):
    conf_seqs = set_conf_seqs(flatten_confidence_seqs(master_conf_seq))
    def helper(c, conf_seq):
        if len(conf_seq) == 0:
            return [c]
        elif (isinstance(conf_seq[0], float) and
              not conf_seq[0] in [next_c.confidence for next_c in c.clades]) or \
             (isinstance(conf_seq[0], set) and
              not conf_seq[0].issubset(set(next_c.confidence for next_c in c.clades))):
            return []
        if isinstance(conf_seq[0], set):
            if helper(c, list(conf_seq[0])):
                return [c]
            else:
                return []
        else:
            output = []
            for next_c in c.clades:
                if next_c.confidence == conf_seq[0]:
                    output.extend(helper(next_c, conf_seq[1:]))
            return output
    output = set(itertools.chain(*[helper(first_c, conf_seq) for conf_seq in conf_seqs]))
    return list(output)

## get clade from multiple confidence values
def get_clades_from_confidence_seq_multi(first_c, conf_seqs):
    return list(itertools.chain(*[get_clades_from_confidence_seq(first_c, conf_seq) for conf_seq in conf_seqs]))

## calculate branch lengths within radiation clads vs non-radiation clades
def get_gene_clades(t_nwk, fout, predicted_identity_f, to_exclude = [],
                    cluster = '', genes = [], seq_names_f = '', seq_names = [],
                    exclude_araly = True, exclude_ref = True):
    with open(predicted_identity_f, 'r') as f:
        dat_pred = [x[:-1].split('\t') for x in f.readlines()]
        get_pred = make_custom_get(dat_pred[0])
        dat_pred = dat_pred[1:]
    if cluster:
        dat_pred = [get_pred(x, "contig", "gene", "cluster", "domain") for x in dat_pred
                    if get_pred(x, "cluster") == cluster]
        seq_names = set(get_pred(x, "contig") for x in dat_pred)
    elif genes:
        dat_pred = [get_pred(x, "contig", "gene", "cluster", "domain") for x in dat_pred
                    if get_pred(x, "gene") in genes]
        seq_names = set(get_pred(x, "contig") for x in dat_pred)
    elif not (seq_names_f or seq_names):
        seq_names = set(get_pred(dat_pred, "contig"))
        dat_pred = [get_pred(x, "contig", "gene", "cluster", "domain") for x in dat_pred]
    else:
        if seq_names_f:
            seq_names = set(x[:-1] for x in open(seq_names_f, 'r').readlines())
        elif seq_names:
            seq_names = set(seq_names)
        dat_pred = [get_pred(x, "contig", "gene", "cluster", "domain") for x in dat_pred
                    if get_pred(x, "contig") in seq_names]
    get_pred_reduced = make_custom_get(["name", "gene", "cluster", "domain"])
    ## assume rooted tree
    t = Phylo.read(t_nwk, "newick")
    # root = prune_tree(t, seq_names, to_prune = to_exclude)
    to_keep = seq_names
    if not exclude_ref:
        to_keep |= set(leaf.name for leaf in t.get_terminals() if "Col-0" in leaf.name)
    root = prune_tree(t, to_keep = to_keep)
    ## start assigning genes to terminal and internal nodes
    assign_leaf_genes(t, dat_pred, get_pred_reduced, unit = "domain")
    assign_inner_genes(t, root = root) ## assumes that all terminal nodes have been assigned a gene
    return (root, define_gene_clades(t, root = root)) ## return poly- and monophyletic gene clades

## requires a tree where all genes have been assigned to the root
## ## ideally, one gene per clade and one clade per gene, but multiple genes will be grouped into a clade in the presence of poly-/paraphyly
def define_gene_clades(t, root = None, sd = 2, i_interval = 20, sd_threshold = 5): ## sd and sd_threshold not used
    ## get root
    if not root:
        print("Getting root")
        root = t.common_ancestor(t.get_terminals())
    ## get mrca of each gene
    genes = root.genes
    print("Getting gene clades")
    gene_clades = []
    for i, gene in enumerate(genes):
        if i % i_interval == 0:
            print(i)
        node = None
        next_node = root
        while next_node == None or not (next_node is node or is_mono(next_node)):
            node = next_node
            next_node = get_containing_child(node, gene)
        gene_clades.append([next_node, next_node.genes])
    ## combine clades that overlap (i.e. non-monophyletic wrt gene)
    print("Merging overlapping clades")
    merged_clades = []
    for clade, genes in gene_clades:
        clade_genes_old = set()
        clade_genes_new = genes
        ## get all clades that overlap with gene set--begin with the original clade's genes as a seed
        while clade_genes_old != clade_genes_new:
            clade_genes_old = clade_genes_new
            clade_genes_new = set(itertools.chain(*[q_genes for q_genes in [x[1] for x in gene_clades]
                                                    if bool(q_genes & clade_genes_old)]))
        # for gene in genes:
        node = None
        next_node = root
        while next_node == None or not (next_node is node or is_mono(next_node, clade_genes_new)):
            node = next_node
            next_node = get_containing_child(node, clade_genes_new)
        merged_clades.append( (next_node, clade_genes_new) )
    ## get unique clade:genes combo
    print("Finding final clades")
    unique_clades = []
    for clade, genes in merged_clades:
        if True in [(clade == clade_o and genes == genes_o) for clade_o, genes_o in unique_clades]:
            continue
        else:
            unique_clades.append( (clade, genes) )
    ## for non-monophyletic genes, re-assign clades
    clades_poly = [clade_gene for clade_gene in unique_clades if len(clade_gene[1]) > 1]
    clades_mono = [clade_gene for clade_gene in unique_clades if len(clade_gene[1]) == 1]
    # # branch_length_stats_mono = {clade: {"genes": genes,
    # #                                     "max": max(get_branch_lengths(clade)),
    # #                                     "sd": numpy.std(get_branch_lengths(clade)),
    # #                                     "median": numpy.median(get_branch_lengths(clade))}
    # #                             for clade, genes in clades_mono}
    # # branch_length_stats_poly = {clade: {"genes": genes,
    # #                                     "max": max(get_branch_lengths(clade)),
    # #                                     "sd": numpy.std(get_branch_lengths(clade)),
    # #                                     "median": numpy.median(get_branch_lengths(clade))}
    # #                             for clade, genes in clades_poly}
    # ## find long branches in monophyletic clades
    # longbranch_mono = []
    # for clade, genes in clades_mono:
    #     branch_lengths = get_branch_lengths(clade)
    #     longbranch_mono.append((clade, genes, numpy.mean(branch_lengths), numpy.std(branch_lengths),
    #                             [(branch_lengths[i], x) for i, x in enumerate(stats.zscore(branch_lengths))
    #                              if abs(x) >= sd_threshold]))
    # longbranch_mono = [len([x > numpy.mean(get_branch_lengths(c)) + y.std(
    # for clade in clades_poly:
    #     ## EHH
    # ## okay I give up. I can't find a way to automate dividing up the tree into sensible monophyletic clades without a lot of manual curation :(
    # ## just return polyphyletic & monophyletic gene clades
    return {"mono": clades_mono, "poly": clades_poly}



## get all branch lengths in clade (excluding current branch)
def get_branch_lengths(c, include_current = False):
    return [x.branch_length for x in c.clades] + \
        list(itertools.chain(*[get_branch_lengths(x) for x in c.clades]))

## requires a clade where genes have be assigned to all inner nodes
## ## The given node, 'c', should also have been assigned to the 'gene'
def get_containing_child(c, genes):
    if c.is_terminal():
        return c
    ## coerce 'gene' into a set
    if isinstance(genes, str):
        genes = {genes}
    else:
        genes = set(genes)
    containing_children = [child for child in c.clades if genes.issubset(child.genes)]
    if len(containing_children) > 1:
        return c
    else:
        return containing_children[0]

def is_para_clade(c):
    child1, child2 = c.clades
    for gene in c.genes:
        if gene in child1 and gene in child2:
            return (child1.genes != child2.genes)
    return True
            

def is_mono(c, genes = []):
    if not genes:
        return len(c.genes) == 1
    else:
        return c.genes == set(genes)

## prune tree and return root
def prune_tree(t, to_keep = [], to_prune = [], i_interval = 500):
    if to_keep:
        to_prune = set(to_prune)
        to_prune |= set(x.name for x in t.get_terminals()) - set(to_keep)
    ## prune unwanted sequences
    print("Pruning")
    for i, terminal_name in enumerate(to_prune):
        if i % i_interval == 0:
            print(i)
        tmp = t.prune(terminal_name)
    print("Getting MRCA")
    return t.common_ancestor(t.get_terminals())

## assign genes and clusters to terminal nodes
def assign_leaf_genes(t, data, get, unit = "gene", i_interval = 1000):
    ## make dictionary of data
    data = {get(dat, "name"): get(dat, unit, "cluster") for dat in data}
    get = make_custom_get([unit, "cluster"])
    for i, leaf in enumerate(t.get_terminals()):
        if i % i_interval == 0:
            print(i)
        if not "Col-0" in leaf.name:
            dat = data[leaf.name]
            leaf.genes = set(get(dat, unit).split(';'))
            leaf.cluster = set(get(dat, "cluster").split(';'))
        else:
            leaf.genes = {re.search("AT\dG\d{5}", leaf.name).group(0)}
            leaf.cluster = {"None"}
    return

## assign genes to inner nodes
## all terminal nodes must have been assigned a gene for it to work properly
## ## if no gene assigned, leaf will be ignored
def assign_inner_genes(t, i_interval = 500, root = None):
    if not root:
        print("Getting root")
        root = t.common_ancestor(t.get_terminals())
    print("Resetting non-terminal genes")
    for nonterm in t.get_nonterminals():
        nonterm.genes = set()
    print("Assigning genes to inner nodes")
    for i, leaf in enumerate(t.get_terminals()):
        if i % i_interval == 0:
            print(i)
        ## get gene assigned to leaf
        try:
            leaf_genes = leaf.genes
        except AttributeError:
            continue
        parent_nodes = [root] + t.get_path(leaf.name)
        ## add gene to parent node if parent node doesn't contain gene
        for parent_node in parent_nodes[-2::-1]:
            if leaf_genes.issubset(parent_node.genes):
                break
            else:
                parent_node.genes |= leaf_genes
    ## PROBLEM: genes from CNLs keep getting assigned to genes from TNLs or vice versa (FIXED!!)
    return
