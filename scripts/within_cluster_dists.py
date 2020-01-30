import os
import re
import sys
from Bio import Phylo

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from data_manip import *

def assign_dist_to_ref(t, q_nodes):
    root = t.root
    for i, q_node in enumerate(q_nodes):
        if i % 500 == 0:
            print("Query", i)
        parent_nodes = [root] + t.get_path(q_node)
        distance = q_node.branch_length
        q_node.ref_distance = {}
        q_node.ref_sisters = {}
        for parent_node in parent_nodes[-2::-1]:
            try:
                q_node.ref_distance = {**q_node.ref_distance,
                                       **{k: v + distance for k, v in parent_node.ref_distance.items()
                                          if not k in q_node.ref_distance}}
                if not q_node.ref_sisters:
                    q_node.ref_sisters = {k: v + distance for k, v in parent_node.ref_distance.items()}
            except AttributeError:
                continue
            try:
                distance += parent_node.branch_length
            except TypeError:
                continue
    return

def assign_inner_dist_to_ref(t, ref_nodes):
    root = t.root
    for ref_node in ref_nodes:
        parent_nodes = [root] + t.get_path(ref_node)
        distance = ref_node.branch_length
        for parent_node in parent_nodes[-2::-1]:
            try:
                parent_node.ref_distance = {**parent_node.ref_distance,
                                            **{ref_node: distance}}
            except AttributeError:
                parent_node.ref_distance = {ref_node: distance}
            try:
                distance += parent_node.branch_length
            except TypeError:
                continue
    return

def get_all_pairwise_distance(t_nwk, fout, predicted_identity_f, cluster = '', genes = [], seq_names_f = '', seq_names = [], all_mrca_desc = False):
    with open(predicted_identity_f, 'r') as f:
        dat_pred = [x[:-1].split('\t') for x in f.readlines()]
        get_pred = make_custom_get(dat_pred[0])
        dat_pred = dat_pred[1:]
    if cluster:
        dat_pred = [get_pred(x, "contig", "gene", "cluster") for x in dat_pred if get_pred(x, "cluster") == cluster]
        seq_names = set(get_pred(x, "contig") for x in dat_pred)
    elif genes:
        dat_pred = [get_pred(x, "contig", "gene", "cluster") for x in dat_pred if get_pred(x, "gene") in genes]
        seq_names = set(get_pred(x, "contig") for x in dat_pred)
    else:
        if seq_names_f:
            seq_names = set(x[:-1] for x in open(seq_names_f, 'r').readlines())
        elif seq_names:
            seq_names = set(seq_names)
        dat_pred = [get_pred(x, "contig", "gene", "cluster") for x in dat_pred if get_pred(x, "contig") in seq_names]
    ref_genes_raw = set(get_pred(x, "gene") for x in dat_pred)
    ref_genes = set()
    for ref_gene in ref_genes_raw:
        ref_genes |= set(ref_gene.split(';'))
    ## assume rooted tree
    t = Phylo.read(t_nwk, "newick")
    nodes = tuple(t.find_elements(name='|'.join(seq_name.replace('|', '\|') for seq_name in seq_names) + "|Col-0_ref.(" + '|'.join(ref_genes) + ")\|.+"))
    mrca = t.common_ancestor(nodes)
    if all_mrca_desc:
        nodes = mrca.get_terminals()
    assign_inner_dist_to_ref(mrca, nodes)
    assign_dist_to_ref(mrca, nodes)
    to_write = ['\n'.join(['\t'.join([node.name, other.name, str(distance)]) for other, distance in node.ref_distance.items()]) for node in nodes]
    f = open(fout, "w+")
    f.write('\n'.join(['\t'.join(["seqa", "seqb", "distance"])] + to_write))
    f.close()
    return

def get_all_pairwise_distance_multi(t_dir, predicted_identity_dir, cluster_gene_f, out_dir,
                                    domain = "NB-ARC", clusters = [], genes = {}, exclude_cluster = ["singleton"],
                                    *args):
    t_nwk = f"{t_dir}/nlr164_col0-AL70-Alyrata_{domain}_mafft_ML.nwk"
    with open(cluster_gene_f, 'r') as f:
        dat_genes = [x[:-1].split('\t') for x in f.readlines()]
        get_genes = make_custom_get(dat_genes[0])
        dat_genes = dat_genes[1:]
    clusters = set(clusters) if clusters else set(get_genes(dat_genes, "cluster")) if not genes else set(genes.keys())
    clusters = clusters - set(exclude_cluster)
    predicted_identity_al70 = f"{predicted_identity_dir}/nlr164_AL70_{domain}_predictedIdentity.txt"
    predicted_identity_araly = f"{predicted_identity_dir}/nlr164_Alyrata_{domain}_predictedIdentity.txt"
    with open(predicted_identity_al70, 'r') as f:
        dat_pred = [x[:-1].split('\t') for x in f.readlines()]
        get_pred = make_custom_get(dat_pred[0])
        dat_pred = dat_pred[1:]
    with open(predicted_identity_araly, 'r') as f:
        dat_pred.extend([x[:-1].split('\t') for x in f.readlines()][1:])
    for cluster in clusters:
        print(f"Working on {cluster}")
        fout = f"{out_dir}/nlr164_{domain}_{cluster}_distances.tsv"
        genes = [get_genes(x, "gene") for x in dat_genes if get_genes(x, "cluster") == cluster]
        seq_names = [get_pred(x, "contig") for x in dat_pred if get_pred(x, "gene") in genes]
        get_all_pairwise_distance(t_nwk, fout, predicted_identity_al70, seq_names = seq_names, all_mrca_desc = False)
    return
