import os
import re
import sys
from Bio import Phylo

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from basics import *
from data_manip import *


def write_closest_ref(t_nwk, fout_al70, fout_alyrata, domain, coi_f, root_pattern = '',
                      ref_pattern = "Col-0_ref.+?complete.*", al70_pattern = "(_R_)?\d+?\.tig\d+?.+",
                      alyrata_pattern = "X[MR]_.+", rooted = False):
    t = Phylo.read(t_nwk, "newick")
    if not rooted:
        if not root_pattern:
            print(f"Rooting {domain} tree at midpoint")
            t.root_at_midpoint()
        elif not rooted:
            print(f"Rooting {domain} tree with outgroup")
            root = tuple(t.find_elements(name=root_pattern))
            t.root_with_outgroup(root)
    ref_nodes_complete = tuple(t.find_elements(name=ref_pattern))
    al70_nodes_predicted = tuple(t.find_elements(name=al70_pattern))
    alyrata_nodes = tuple(t.find_elements(name=alyrata_pattern))
    collapse_iter = lambda iterable, char: [char.join([str(y) for y in x]) for x in iterable]
    def generate_to_write (nodes):
        to_write = '\n'.join(['\t'.join(["contig", "gene", "cluster", "distance", "ref_name", "sisters"])] +
                             ['\t'.join(collapse_iter([(node.name,), node.gene, node.cluster,(node.closest_dist,),
                                                       tuple(r.name for r in node.closest_ref),
                                                       tuple(collapse_iter([(sis.name, d) for sis, d in
                                                                            node.ref_sisters.items()], ','))],
                                                      ';'))
                              for node in nodes])
        return to_write
    print("Working on AL70")
    assign_closest_gene_and_cluster(t, al70_nodes_predicted, ref_nodes_complete, coi_f, domain,
                                    uncategorised = "singletons")
    f = open(fout_al70, "w+")
    f.write(generate_to_write(al70_nodes_predicted))
    f.close()
    print("Working on A. lyrata")
    assign_closest_gene_and_cluster(t, alyrata_nodes, ref_nodes_complete, coi_f, domain,
                                    uncategorised = "singletons")
    f = open(fout_alyrata, "w+")
    f.write(generate_to_write(alyrata_nodes))
    f.close()
    return

def assign_closest_gene_and_cluster(t, q_nodes, ref_nodes, coi_f, domain,
                                    uncategorised = "singletons"):
    clusters = read_clusters(coi_f)
    print("Assigning genes and clusters to reference nodes")
    assign_gene_cluster_to_ref(ref_nodes, clusters, uncategorised = "singletons")
    print("Finding inner distance to references")
    assign_inner_dist_to_ref(t, ref_nodes)
    print("Finding distances from references to queries")
    assign_dist_to_ref(t, q_nodes)
    print("Finding closest reference to queries")
    for i, q_node in enumerate(q_nodes):
        closest_ref = min(q_node.ref_distance, key = lambda x: q_node.ref_distance[x])
        q_node.closest_dist = q_node.ref_distance[closest_ref]
        q_node.closest_ref = tuple(ref_node for ref_node, dist in q_node.ref_distance.items() if
                                   dist == q_node.closest_dist)
        q_node.gene = tuple(ref_node.gene for ref_node in q_node.closest_ref)
        q_node.cluster = tuple(ref_node.cluster for ref_node in q_node.closest_ref)
        q_node.domain = domain
        # q_node.isoform_domain = tuple(ref_node.isoform_domain for ref_node in q_node.closest_ref)
    return

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

def read_clusters(coi_f):
    with open(coi_f, 'r') as f:
       coi = [x[:-1].split('\t') for x in f.readlines()]
    get = make_custom_get(coi[0])
    coi = coi[1:]
    output = {cluster: {"genes": set()} for cluster in set(get(coi, "cluster"))}
    for entry in coi:
        output[get(entry, "cluster")]["genes"] |= {get(entry, "gene")}
    return output

def assign_gene_cluster_to_ref(ref_nodes, clusters, uncategorised = "singletons"):
    for ref_node in ref_nodes:
        ref_node.gene = re.search("AT.G\d+", ref_node.name).group(0)
        # ref_node.isoform_domain = re.search("(?<=\|)AT.G\d+?.\d+?\|.+?\|\d+", ref_node.name).group(0)
        possible_clusters = tuple(k for k, v in clusters.items() if ref_node.gene in v["genes"])
        ref_node.cluster = uncategorised if not possible_clusters else possible_clusters[0]
    return

def get_closest_ref(dist_to_ref):
    output = {q_node: min(v, key = lambda x: v[x]) for q_node, v in dist_to_ref.items()}
    return output    

def get_dist_to_ref(t, q_nodes, ref_nodes):
    output = {q_node: {ref_node: t.distance(q_node, ref_node)
                       for ref_node in ref_nodes}
              for q_node in q_nodes}
    return output

## read tree.
# coi_f = "/mnt/chaelab/rachelle/hap_tags/results/nlr165/nlr165_final_wPred.tsv"

def get_closest_ref_nlr164(domain, prefix = 'nlr164', root_pattern = '', ref_pattern = "Col-0_ref.+",
                           al70_pattern = "(_R_)?\d+?.+?\.tig\d+?.+",
                           alyrata_pattern = "\d+?\|.+", t_nwk = '', out_dir = '',
                           **kwargs):
    if not t_nwk:
        t_nwk = "/mnt/chaelab/rachelle/hap_tags/results/nlr165/tree/nlr165_col0-AL70-Alyrata_" + domain + "_mafft_ML.nwk"
    if not out_dir:
        fout_al70 = "/mnt/chaelab/rachelle/hap_tags/results/nlr165/predicted_identity/nlr165_AL70_" + domain + "_predictedIdentity.txt"
        fout_alyrata = "/mnt/chaelab/rachelle/hap_tags/results/nlr165/predicted_identity/nlr165_Alyrata_" + domain + "_predictedIdentity.txt"
    else:
        fout_al70 = f"{out_dir}/{prefix}_AL70_{domain}_predictedIdentity.txt"
        fout_alyrata = f"{out_dir}/{prefix}_Alyrata_{domain}_predictedIdentity.txt"
    
    write_closest_ref(t_nwk, fout_al70, fout_alyrata, domain, coi_f = kwargs["coi_f"],
                      root_pattern = root_pattern, ref_pattern = ref_pattern, al70_pattern = al70_pattern,
                      alyrata_pattern = alyrata_pattern, **{k: v for k, v in kwargs.items() if k != "coi_f"})
    return
