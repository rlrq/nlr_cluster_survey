import os
import sys
import subprocess

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__))))
from file_manip import *

def fasta_to_dict(fname):
    """
    returns dictionary of sequences indexed by sequence name
    """
    from Bio import SeqIO
    seqs = {}
    for seq_record in SeqIO.parse(fname, "fasta"):
        seqs[seq_record.id] = seq_record.seq
    return seqs

def dict_to_SeqRecordList(d, description = '', seq_id_func = lambda x:x):
    out_l = []
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC, Gapped
    for seq_id,seq in d.items():
        out_l.append(SeqRecord(seq if isinstance(seq, Seq) else
                               Seq(seq, Gapped(IUPAC.unambiguous_dna)),
                               id = seq_id_func(seq_id), description = description))
    return out_l

def dict_to_fasta(d, fout):
    from Bio import SeqIO
    SeqIO.write(dict_to_SeqRecordList(d), fout, "fasta")

# fasta="/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta"
# fout="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0_10kbflank/anna_lena70.Col-0_10kbflank.blastn_topContigs_AT5G58120.fasta"
# titles="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0_10kbflank/anna_lena70.Col-0_10kbflank.blastn_topContigs_AT5G58120.txt"
def filter_title(fasta, fout, titles, match_exact = False, delim = '\t'):
    """
    keep only sequences with titles that (partially) match a set of specified titles
    """
    seq_dict = fasta_to_dict(fasta)
    titles_list = single_col_to_list(titles)
    output = []
    for title in titles_list:
        output.extend([(k,v) for k,v in seq_dict.items() if title == k] if match_exact else \
                      [(k,v) for k,v in seq_dict.items() if title in k])
    dict_to_fasta(dict(output), fout)
    
def trim_alignment(fasta, ref_title, fout):
    """
    trims alignment to start and end of reference sequence (specified by title using 'ref_title' variable; partial match allowed; user must provide unique pattern, this function doesn't check for uniqueness)
    """
    seq_dict = fasta_to_dict(fasta)
    ref_seq = [v for k,v in seq_dict.items() if ref_title in k][0]
    first_base, last_base = 0, len(ref_seq)
    for i,v in enumerate(ref_seq):
        if v != '-':
            first_base = i
            break
    for i,v in reversed(list(enumerate(ref_seq))):
        if v != '-':
            last_base = i
            break
    output = {}
    for k,v in seq_dict.items():
        output[k] = v[first_base:last_base]
    dict_to_fasta(output, fout)

    
def trim_seq(fasta, start, end, fout):
    """
    Trims sequences in a file. Uses 0-indexing.
    """
    seq_dict = fasta_to_dict(fasta)
    output = {}
    for k,v in seq_dict.items():
        if start < 0:
            tstart, tend = len(v)+start+1, len(v)+end+1
            start, end = min(tstart, tend), max(tstart, tend)
        output[k + "_pos" + str(start) + '-' + str(end)] = v[start:end]
    dict_to_fasta(output, fout)

    
def split_seq(fasta, fout, size = None, num = None):
    """
    Splits sequences in file according to final seq length or number of sequences
    """
    seq_dict = fasta_to_dict(fasta)
    output = {}
    import math
    for k, v in seq_dict.items():
        curr_window = size if size else math.ceil(len(v)/num)
        i = 0
        while i < len(v):
            output[k + "_pos" + str(i) + '-' + str(min(i+curr_window, len(v)))] = v[i:i+curr_window]
            i += curr_window
    dict_to_fasta(output, fout)
    
    
def extract_feature(fasta, gff_bed, ref_title, geneID, start_pos, fout, feature_type,
                    geneID_col = -1, feature_col = 7, phase_col = None, encoding = "utf-8"):
    """
    based on a single reference sequence, all sequences in alignment are trimmed to only retain what aligns with a specific type of feature of the specified reference sequence
    """
    # with open(gff_bed, 'r') as f:
    #     gff_gene_raw = [x[:-1].split('\t') for x in f.readlines()]
    # gff_gene = [x for x in gff_gene_raw if (geneID in x[geneID_col] and\
    #                                         x[feature_col] == feature_type)]
    # cds_ranges = [(int(x[1]) - start_pos, int(x[2]) - start_pos) for x in gff_gene]
    grep = subprocess.Popen(("grep", geneID, gff_bed), stdout=subprocess.PIPE)
    awk = subprocess.Popen(("awk", "$8 == \"" + feature_type + "\""), stdin=grep.stdout, stdout=subprocess.PIPE)
    merged = [entry.split('\t') for entry in
              subprocess.check_output(("bedtools", "merge"),
                                      stdin=awk.stdout).decode(encoding).split('\n') if entry]
    cds_ranges = [(int(x[1]) - start_pos, int(x[2]) - start_pos) for x in merged]
    seq_dict = fasta_to_dict(fasta)
    ref_seq = [v for k,v in seq_dict.items() if ref_title in k][0]
    adj_ranges = [(adjusted_pos(ref_seq,x[0]), adjusted_pos(ref_seq, x[1])) for x in cds_ranges]
    output = {k:extract_ranges(v, adj_ranges) for k,v in seq_dict.items()}
    dict_to_fasta(output, fout)

# fasta = "/mnt/chaelab/rachelle/TE_SV/results/alignment_20190529/AT5G58120_10kbflank_80acc_Col-0_mafft.fa"
# gff_bed = "/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed"
# geneID = "AT5G58120"
geneID_col, feature_col, phase_col = -1, 7, 8
# start_pos = 23507256-1 ## make sure it's 0-indexed
# ref_title = "Col0"
# fout = "/mnt/chaelab/rachelle/TE_SV/results/alignment_20190529/AT5G58120_10kbflank_80acc_Col-0_mafft_CDSonly.fa"
def extract_CDS(fasta, gff_bed, ref_title, geneID, start_pos, fout,
                geneID_col = -1, feature_col = 7, phase_col = 8):
    """
    based on a single reference sequence, all sequences in alignment are trimmed to only retain what aligns with the CDS of the specified reference sequence
    """
    extract_feature(fasta, gff_bed,ref_title, geneID, start_pos, fout, "CDS",
                    geneID_col = geneID_col, feature_col = feature_col, phase_col = phase_col)

def extract_gene(fasta, gff_bed, ref_title, geneID, start_pos, fout,
                geneID_col = -1, feature_col = 7):
    """
    based on a single reference sequence, all sequences in alignment are trimmed to only retain what aligns with the gene body of the specified reference sequence
    """
    extract_feature(fasta, gff_bed,ref_title, geneID, start_pos, fout, "gene",
                    geneID_col = geneID_col, feature_col = feature_col, phase_col = phase_col)

def extract_CDS_seqs(fasta, gff_bed, fout, fgeneIDs,
                    geneID_col = -1, feature_col = 7, phase_col = None):
    """
    Extracts CDS sequences from genome fasta file (cat together by parent protein) based on gene titles
    and a gff3 file converted to bed format
    (ONLY WORKS FOR ARABIDOPSIS GENOME)
    """
    with open(gff_bed, 'r') as f:
        gff_gene_raw = [x[:-1].split('\t') for x in f.readlines()]
    geneIDs = '(' + ')|('.join(single_col_to_list(fgeneIDs)) + ')'
    import re
    gff_gene = [x for x in gff_gene_raw if (re.findall(geneIDs, x[geneID_col]) and
                                            x[feature_col] == "CDS")]
    genome = fasta_to_dict(fasta)
    cds_seqs = {}
    for parent in sorted(list(set([x[geneID_col] for x in gff_gene]))):
        cds_ranges = [x[:6] for x in gff_gene if x[geneID_col] == parent]
        cds_seqs[parent] = extract_ranges(genome[cds_ranges[0][0]],
                                          [(int(x[1]),int(x[2])) for x in cds_ranges],
                                          strand = cds_ranges[0][5])
    dict_to_fasta(cds_seqs, fout)
    
def translate_aligned_CDS(fasta, ref_title, fout, stop_symbol = '*'):
    seq_dict = fasta_to_dict(fasta)
    ref_seq = [v for k,v in seq_dict.items() if ref_title in k][0]
    for k,seq in seq_dict.items():
        tmp = seq.ungap('-')
        for i,v in enumerate(seq):
            if v != '-':
                tmp = 'n' * ((i - ref_seq[:i].count('-')) % 3) + tmp
                break
        tmp += '' if (len(tmp) % 3 == 0) else 'n' if (len(tmp) % 3 == 2) else 'nn'
        seq_dict[k] = tmp.translate(stop_symbol=stop_symbol)
    dict_to_fasta(seq_dict, fout)

def adjusted_pos(seq, pos):
    """
    returns position, adjusted for gaps in alignment
    """
    last_pos = 0
    while True:
        curr_gaps = seq[last_pos:pos].count('-')
        if curr_gaps == 0:
            return pos
        last_pos = pos
        pos += curr_gaps

def extract_ranges(seq, ranges, strand = '+'):
    ranges_sorted = sorted(ranges, key = lambda x: int(x[0]), reverse = (strand == '-'))
    output = seq[:0]
    for start, end in ranges:
        output += seq[int(start):int(end)] if strand == '+' else \
                  seq[int(end)-1:int(start)-1:-1]
    return output
    
def revcomp_seqs(fname, seq_names):
    """
    returns dictionary of sequences indexed by name
    """
    seq_dict = fasta_to_dict(fname)
    fin_dict = {}
    for k,v in seq_dict.items():
        if k in seq_names:
            fin_dict["R " + k] = v.reverse_complement()
        else:
            fin_dict[k] = v
    return fin_dict
