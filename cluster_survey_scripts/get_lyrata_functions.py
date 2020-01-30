import os
import sys
import subprocess


sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from fasta_manip import fasta_to_dict, dict_to_fasta, extract_ranges
from data_manip import make_custom_get
from str_manip import split_extract


def has_overlap(r1, r2):
    return (r2[0] <= r1[0] <=r2[1]) or (r1[0] <= r2[0] <= r1[1])

def grep_bedmerge(bed, fasta, feature, out_dir = '', gene = '', encoding = "utf-8", merge = False):
    ## get chrom, start, end of gene's features
    if gene:
        first_pass = subprocess.Popen(("grep", gene, bed), stdout=subprocess.PIPE)
    else:
        first_pass = subprocess.Popen(("cat", bed), stdout=subprocess.PIPE) 
    if merge:
        awk = subprocess.Popen(("awk", "$5 == \""+feature+"\""),
                               stdin=first_pass.stdout, stdout=subprocess.PIPE)
        data = [entry.split('\t') for entry in
                subprocess.check_output(("bedtools","merge","-c","7,8,9","-o","collapse,collapse,collapse"),
                                        stdin=awk.stdout).decode(encoding).split('\n') if entry]
    else:
        data = [entry.split('\t') for entry in
                subprocess.check_output(("awk", "$5 == \""+feature+"\""),
                                        stdin=first_pass.stdout).decode(encoding).split('\n') if entry]
        data = [[x[i] for i in (0,1,2,6,7,8)] for x in data]
    return data

## output is...0-indexed, I think...it seems like it follows .bed conventions...
def get_domain_in_genome_coords(query, domain, domain_f, out_dir, fasta, bed,
                                pos_type="aa", encoding="utf-8", start_inc = True, end_inc = True,
                                qname_dname=("name", "domain name"), qstart_qend=("start", "end"),
                                restrict_pid_exact = False):
    qname, dname = qname_dname
    qstart, qend = qstart_qend
    with open(domain_f, 'r') as f:
        domain_raw = [x[:-1].split('\t') for x in f.readlines()]
    domain_header = domain_raw[0]
    domain_data = [x for x in domain_raw[1:] if len(domain_raw) > 1 and len(x) == len(domain_header) and \
                   (((not domain) or domain == x[domain_header.index(dname)]) and \
                    query in x[domain_header.index(qname)])]
    output = []
    for domain_dat in domain_data:
        ## convert to 0-index, start inclusive, stop exclusive + adjust for unit (e.g. aa or nt)
        def adj_pos(start, end):
            account_unit = lambda x: x * (3 if pos_type == "aa" else 1)
            return (account_unit(int(start) - (1 if start_inc else 0)),
                    account_unit(int(end) - (0 if end_inc else 1)))
        domain_start, domain_end = adj_pos(int(domain_dat[domain_header.index(qstart)]),
                                           int(domain_dat[domain_header.index(qend)]))
        ## get CDS to define boundaries of translated nucleotides
        cds_dat = [x for x in grep_bedmerge(bed, fasta, "CDS", out_dir = out_dir, gene = query)
                   if ( (not restrict_pid_exact) or (restrict_pid_exact and '=' + query + ';' in x[-1]) )]
        ## extract boundaries
        chrom = cds_dat[0][0]
        plus_strand = (cds_dat[0][3] == '+')
        bounds = sorted([(int(x[1]), int(x[2])) for x in cds_dat],
                        key = lambda x: x[0], reverse = (not plus_strand))
        last_end = 0
        genome_start, genome_end = None, None
        for i, v in enumerate(bounds):
            curr_cds_start = last_end
            curr_cds_end = last_end + (v[1] - v[0])
            get_g_pos = lambda x: ((v[0] + (x - curr_cds_start)) if plus_strand else\
                                   (v[1] - (x - curr_cds_start)))
            if curr_cds_start <= domain_start < curr_cds_end:
                genome_start = get_g_pos(domain_start)
            if curr_cds_start < domain_end <= curr_cds_end: ## recall that end is exclusive
                genome_end = get_g_pos(domain_end)
            last_end = curr_cds_end
        output.append((min(genome_start, genome_end), max(genome_start, genome_end)))
    return output

def get_cds(fout, fasta, bed, domain_f = '', domain = '',
            complete = False, adjust_dir = False, translate = False, protein_id_field = "protein_id",
            domain_pid_f = lambda x: x, **kwargs):
    
    ## extract isoform ranges
    data = grep_bedmerge(bed, fasta, "CDS")
    data.sort(key = lambda x: x[-1])
    get = make_custom_get(["chrom", "start", "end", "strand", "phase", "name"])
    ## extract domain ranges
    if domain and domain_f:
        with open(domain_f) as f:
            domain_data = [x[:-1].split('\t') for x in f.readlines()]
        domain_get = make_custom_get(domain_data[0])
        domain_data = [x for x in domain_data[1:] if domain_get(x, "domain") == domain]
        domain_seqs = set(domain_pid_f(x) for x in domain_get(domain_data, "qseqid"))
        # domain_data = {seq_name: [x for x in domain_data if domain_get(x, "qseqid") == seq_name]
        #                for seq_name in domain_seqs}
        data = [x for x in data if split_extract(split_extract(get(x, "name"), protein_id_field + '=', 1), ';', 0) in domain_seqs]
    isoforms_dat = {}
    isoform = split_extract(split_extract(get(data[0], "name"), protein_id_field + '=', 1), ';', 0)
    isoform_dat = []
    for i, entry in enumerate(data):
        curr_isoform = split_extract(split_extract(get(entry, "name"), protein_id_field + '=', 1), ';', 0)
        if curr_isoform == isoform:
            isoform_dat.append(get(entry, "chrom", "start", "end", "strand", "phase"))
        else:
            ## if complete, get min, max of combined CDS in the isoform
            if complete:
                tmp_isoform_dat = isoform_dat[0]
                tmp_isoform_dat[1] = min(get(isoform_dat, "start"))
                tmp_isoform_dat[2] = max(get(isoform_dat, "end"))
                isoforms_dat[isoform] = [tmp_isoform_dat]
            else:
                isoforms_dat[isoform] = isoform_dat
            isoform_dat = [get(entry, "chrom", "start", "end", "strand", "phase")]
            isoform = curr_isoform
        ## if last entry, write to isoforms_dat
        if i == len(data) - 1:
            isoforms_dat[isoform] = isoform_dat
    
    ## domain magic
    if domain_f and domain:
        updated_isoforms_dat = {}
        for isoform in domain_seqs:
            if not isoform in isoforms_dat:
                continue
            isoform_dat = isoforms_dat[isoform]
            chrom, strand = get(isoform_dat[0], "chrom", "strand")
            domain_ranges = get_domain_in_genome_coords(isoform, domain, domain_f,
                                                        os.path.dirname(fout),
                                                        qname_dname = ("qseqid", "domain"),
                                                        qstart_qend = ("qstart", "qend"),
                                                        bed = bed, fasta = fasta,
                                                        **kwargs)
            for i, domain_dat in enumerate(domain_ranges):
                d_start, d_end = domain_dat
                domain_cds = [[chrom, max(min(x), d_start), min(max(x), d_end), strand]
                              for x in get(isoform_dat, "start", "end") if
                              has_overlap(x, (d_start, d_end))]
                updated_isoforms_dat[isoform + f"|{domain}|{i+1}"] = domain_cds
        isoforms_dat = updated_isoforms_dat
    
    ## extract sequences
    ref_seqs = {split_extract(k, ' ', 0): v for k, v in fasta_to_dict(fasta).items()}
    isoforms_seq = {}
    for isoform, isoform_dat in isoforms_dat.items():
        if not isoform_dat:
            continue
        chrom, strand = get(isoform_dat[0], "chrom", "strand")
        isoform_seq = extract_ranges(ref_seqs[str(chrom)],
                                     get(isoform_dat, "start", "end"))
        if (adjust_dir or (translate and not complete)) and strand == '-':
            isoform_seq = isoform_seq.reverse_complement()
        if translate and not complete:
            isoform_seq = isoform_seq.translate()
        isoforms_seq[isoform +
                     ("|revcomp" if strand == '-' and (adjust_dir or (translate and not complete)) else '') +
                     ("|complete" if complete else '')] = isoform_seq
    dict_to_fasta(isoforms_seq, fout)
    return
