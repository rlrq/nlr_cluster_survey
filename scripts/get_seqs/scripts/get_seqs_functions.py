import re
import os
import sys
import itertools
import unidecode
import subprocess
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__))))

chrom_pref="Chr"
data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "data", "1135acc.csv")
# bed_path = "/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed"
bed_path = ''
ref_fasta_pref = "/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_"
ref_fasta = {chrom_pref + '1': ref_fasta_pref + "chr01.fasta",
             chrom_pref + '2': ref_fasta_pref + "chr02.fasta",
             chrom_pref + '3': ref_fasta_pref + "chr03.fasta",
             chrom_pref + '4': ref_fasta_pref + "chr04.fasta",
             chrom_pref + '5': ref_fasta_pref + "chr05.fasta",
             chrom_pref + 'M': ref_fasta_pref + "mt.fasta",
             chrom_pref + 'C': ref_fasta_pref + "pltd.fasta"}

def has_overlap(r1, r2):
    r1 = sorted(r1)
    r2 = sorted(r2)
    return (r2[0] <= r1[0] <=r2[1]) or (r1[0] <= r2[0] <= r1[1])

def has_any_overlap(l1, l2):
    any_overlap = [has_overlap(r1, r2) for r1 in l1 for r2 in l2]
    return True in any_overlap

def merge_ranges(*l):
    all_ranges = list(sorted(itertools.chain(*l)))
    if len(all_ranges) <= 1:
        return(all_ranges)
    final_ranges = [tuple(all_ranges.pop(0))]
    if has_overlap(final_ranges[-1], all_ranges[0]):
        r1, r2 = tuple(final_ranges.pop(-1)), tuple(all_ranges.pop(0))
        final_ranges.append(tuple(min(*r1, *r2), max(*r1, *r2)))
    else:
        final_ranges.append(all_ranges.pop(0))
    return(final_ranges)

def grep_bedmerge(gene, bed, feature, encoding, out_dir, merge = False):
    ## get chrom, start, end of gene's features
    grep = subprocess.Popen(("grep", gene, bed), stdout=subprocess.PIPE)
    if merge:
        awk = subprocess.Popen(("awk", "$8 == \""+feature+"\""), stdin=grep.stdout, stdout=subprocess.PIPE)
        data = [entry.split('\t') for entry in
                subprocess.check_output(("bedtools","merge","-c","6,8,10","-o","collapse,collapse,collapse"),
                                        stdin=awk.stdout).decode(encoding).split('\n') if entry]
    else:
        data = [entry.split('\t') for entry in
                subprocess.check_output(("awk", "$8 == \""+feature+"\""),
                                        stdin=grep.stdout).decode(encoding).split('\n') if entry]
        data = [[x[i] for i in (0,1,2,5,7,9)] for x in data]
    ## write bed file of regions used (0-indexed)
    fout = os.path.join(out_dir, "bed", gene + '_' + feature + ".bed")
    os.makedirs(os.path.dirname(fout), exist_ok=True)
    f = open(fout, "w+")
    f.write('\n'.join(['\t'.join(x) for x in data]))
    f.close()
    return {"fout": fout, "data": data}

## only accepts .tsv and .csv file if 'delim' is not specified.
## accepts a tuple (or list) of accession names
## returns [(acc_num:acc_name)] list, all strings
def get_acc_num(acc_names, ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
    ## read data
    if not delim:
        delim = ',' if ref_file[-4:] == ".csv" else '\t' if ref_file[-4:] == ".tsv" else ''
    if not delim:
        print("Unable to detect delimiter. Please specify delimiter with 'delim' keyword.")
        sys.exit(1)
    with open(ref_file, 'r') as f:
        ## convert all names to lower case, only keep alphanumeric characters (data)
        data = {(re.sub(r'[^a-zA-Z0-9]', '',
                        unidecode.unidecode(x.split(delim)[name_col]))).lower():x.split(delim)[num_col]\
                for x in f.readlines()[0 if not header else 1:]}
    ## convert all names to lower case, only keep alphanumeric characters (query), and get acc num
    return [(data[(re.sub(r'[^a-zA-Z0-9]', '', unidecode.unidecode(x))).lower()], x) for x in acc_names]

## only accepts .tsv and .csv file if 'delim' is not specified.
## accepts a tuple (or list) of accession names
## returns [(acc_num:acc_name)] list, all strings
def get_acc_name(acc_nums, ref_file=data_path, num_col=1, name_col=0, delim='', header=True):
    ## read data
    with open(ref_file, 'r') as f:
        if not delim:
            delim = ',' if (ref_file[-4:] == ".csv" or ref_file[-4:] == ".txt") else \
                    '\t' if ref_file[-4:] == ".tsv" else ''
        if not delim:
            print("Unable to detect delimiter. Please specify delimiter with 'delim' keyword.")
            sys.exit(1)
        data = {x.split(delim)[num_col]:x.split(delim)[name_col]
                for x in f.readlines()[0 if not header else 1:]}
    return [(str(x), data[str(x)]) for x in acc_nums]

## returns address from which to retrieve sequence data from 1001genomes
def make_1001pseudogenome_api(chrom, start, end, accs, add_chr="Chr", remove_chr=True):
    chrom = add_chr + (chrom if not remove_chr else re.search("[Cc]hr(.+)", chrom).group(1))
    return "http://tools.1001genomes.org/api/v1/pseudogenomes/strains/" + ','.join(map(str, accs)) +\
        "/regions/" + chrom + ':' + str(start) + ".." + str(end)

## parse accession num/name
def raw_accs_to_id(accs, ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
    accs = [acc for acc in accs if acc]
    if str(accs[0]).isdigit():
        accs = dict(get_acc_name(accs, ref_file=ref_file, num_col=num_col, name_col=name_col,
                                 delim=delim, header=header))
    else:
        accs = dict(get_acc_num(accs, ref_file=ref_file, num_col=num_col, name_col=name_col,
                                delim=delim, header=header))
    return accs

def get_1001pseudogenome_raw(accs, chrom, start, end, fasta_out,
                             ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
    accs = raw_accs_to_id(accs, ref_file=data_path, num_col=0, name_col=1, delim='', header=True)
    subprocess.run(("curl", "-o", fasta_out,
                    make_1001pseudogenome_api(chrom, start, end, accs.keys())))
    return fasta_out

## get sequence info using genomic range
def get_1001pseudogenome_by_range(accs, chrom, start, end, out_dir,
                                  ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
    fasta_out = os.path.join(out_dir, "chr{}_p{}-{}.fasta".format(chrom, start, end))
    get_1001pseudogenome_raw(accs, chrom_pref + chrom, start, end, fasta_out,
                             ref_file=ref_file, num_col=num_col, name_col=name_col, delim=delim,
                             header=header)
    ## ensure nice line width
    from fasta_manip import fasta_to_dict
    tmp = fasta_to_dict(fasta_out)
    from fasta_manip import dict_to_fasta
    dict_to_fasta(tmp, fasta_out)
    print("Sequences were successfully written to: {}".format(fasta_out))
    
## get sequence info of genes
def get_1001pseudogenome_by_gene(accs, gene, feature, out_dir, bed=bed_path, encoding="utf-8",
                                 ref_file=data_path, num_col=0, name_col=1, delim='', header=True,
                                 complete=False, domain="", domain_f="", merge=False, adj_dir=False,
                                 **for_get_domain_in_genome_coords):
    data = grep_bedmerge(gene, bed, feature, encoding, out_dir, merge=merge)["data"]
    ## extract sequences
    chrom = data[0][0]
    start = min([int(x[1]) for x in data])
    end = max([int(x[2]) for x in data])
    strand = data[0][3]
    print(data[0], chrom, start, end)
    if feature == "gene":
        isoforms = {gene: data}
    else:
        isoforms = {re.search("=(.+\.\d+)(?=[,;]|$)", isoform).group(1): [x for x in data if x[-1] == isoform]
                    for isoform in set([x[-1] for x in data])}
    fasta_out_l = []
    ## iterate through isoforms
    for isoform, isoform_dat in isoforms.items():
        ## if extracting only domain-specific range
        if domain_f:
            d_start, d_end = get_domain_in_genome_coords(gene, domain, domain_f, out_dir,
                                                         bed=bed, encoding=encoding,
                                                         isoform = isoform,
                                                         start_inc = start_inc, end_inc = end_inc,
                                                         **{k: v for k, v
                                                            in for_get_domain_in_genome_coords.items()
                                                            if k in ["qname_dname", "qstart_qend"]})
            if (not d_start) or (not d_end):
                continue
        ## get fasta file of sequence data
        fasta_out = os.path.join(out_dir, isoform + '_' + feature + ("_complete" if complete else '') + \
                                 (('_' + ("domain" if not domain else domain))if domain_f else '') + ".fasta")
        get_1001pseudogenome_raw(accs, chrom, start + 1, end, fasta_out)
        accs = raw_accs_to_id(accs, ref_file=data_path, num_col=0, name_col=1, delim='', header=True)
        ## rename seqs + trim if required
        from fasta_manip import fasta_to_dict
        tmp_seqs = fasta_to_dict(fasta_out)
        seqs = {"{}|{}|{}|{}|{}".format(k.split('|')[3], accs[k.split('|')[3]], gene, feature, isoform) +\
                (('|' + ("domain" if not domain else domain)) if domain_f else '' +
                 ("|complete" if complete else '')): seq \
                for k, seq in tmp_seqs.items()}
        if (not complete) or domain_f:
            if complete and domain_f and d_start and d_end:
                ranges = [(max(start, d_start) - start, min(end, d_end) - start)]
            elif domain_f and d_start and d_end:
                ranges = [(max(int(x[1]), d_start) - start, min(int(x[2]), d_end) - start) \
                          for x in isoform_dat if has_overlap((int(x[1]), int(x[2])), (d_start, d_end))]
            else:
                ranges = [(int(x[1])-start, int(x[2])-start) for x in isoform_dat]
            from fasta_manip import extract_ranges
            seqs = {k: extract_ranges(seq, ranges) for k, seq in seqs.items()}
        if adj_dir and strand == '-':
            seqs = {k + "|revcomp": seq.reverse_complement() for k, seq in seqs.items()}
        from fasta_manip import dict_to_fasta
        dict_to_fasta(seqs, fasta_out)
        fasta_out_l.append(fasta_out)
    fasta_out_final = os.path.join(out_dir, gene + '_' + feature + ("_complete" if complete else '') + \
                                   (('_' + ("domain" if not domain else domain)) if domain_f else '') + ".fasta")
    if fasta_out_l and fasta_out_l[0] != fasta_out_final:
        from file_manip import cat_files
        cat_files(sorted(fasta_out_l), fasta_out_final)
        for fasta_out in fasta_out_l:
            os.remove(fasta_out)
    print("Sequences were successfully written to: {}".format(fasta_out_final))

###################
## get reference ##
###################

def get_ref_raw(chrom, start, end, fasta_out, encoding="utf-8",
            ref_fasta_files=ref_fasta):
    from fasta_manip import fasta_to_dict
    ref_seq = list(fasta_to_dict(ref_fasta_files[chrom]).values())[0][start:end] ## 0-indexed
    from fasta_manip import dict_to_fasta
    dict_to_fasta({"Col-0_ref|{}:{}..{}".format(chrom, start + 1, end): ref_seq}, fasta_out)

## note: setting by_gene to True will collapse identical entries from all isoforms
def get_ref_by_gene(gene, feature, out_dir, bed=bed_path, encoding="utf-8",
                    ref_fasta_files=ref_fasta, complete=False, domain="", domain_f="",
                    start_inc=True, end_inc=True, merge=False, translate=False, adj_dir=False,
                    by_gene=False, **for_get_domain_in_genome_coords):
    if domain_f:
        feature = "CDS"
    data = grep_bedmerge(gene, bed, feature, encoding, out_dir, merge=merge)["data"]
    ## extract sequences from fasta file
    if not data:
        return
    chrom = data[0][0]
    start = min([int(x[1]) for x in data])
    end = max([int(x[2]) for x in data])
    strand = data[0][3]
    if feature == "gene":
        isoforms = {gene: data}
    else:
        isoforms = {re.search("=(.+\.\d+)(?=[,;]|$)", isoform).group(1): [x for x in data if x[-1] == isoform]
                    for isoform in set([x[-1] for x in data])}
    fasta_out_l = []
    seq_ranges = {}
    ## iterate through isoforms
    for isoform, isoform_dat in isoforms.items():
        ## if extracting only domain-specific range
        if domain_f:
            domain_data = get_domain_in_genome_coords(gene, domain, domain_f, out_dir,
                                                      bed=bed, encoding=encoding,
                                                      isoform = isoform,
                                                      start_inc = start_inc, end_inc = end_inc,
                                                      **{k: v for k, v
                                                         in for_get_domain_in_genome_coords.items()
                                                         if k in ["qname_dname", "qstart_qend"]})
            if (not domain_data):
                continue
        else:
            domain_data = [(start, end)]
        ## get fasta file of sequence data
        fasta_out = os.path.join(out_dir, isoform + "_ref_" + feature + ("_complete" if complete else '') + \
                                 (('_' + ("domain" if not domain else domain)) if domain_f else '') + \
                                 ("_protein" if (translate and (feature=="CDS")) else '') + ".fasta")
        get_ref_raw(chrom, start, end, fasta_out, encoding=encoding, ref_fasta_files=ref_fasta_files)
        seqs_to_write = {}
        for i, domain_range in enumerate(domain_data):
            from fasta_manip import fasta_to_dict
            ref_seq = list(fasta_to_dict(fasta_out).values())[0]
            d_start, d_end = domain_range
            ranges = [(d_start, d_end)]
            ## trim sequence if complete flag not raised or if domain required
            if (not complete) or domain_f:
                if complete and domain_f and d_start and d_end:
                    ranges = [(max(start, d_start) - start, min(end, d_end) - start)]
                elif domain_f and d_start and d_end:
                    ranges = [(max(int(x[1]), d_start) - start, min(int(x[2]), d_end) - start) \
                              for x in isoform_dat if has_overlap((int(x[1]), int(x[2])), (d_start, d_end))]
                else:
                    ranges = [(int(x[1])-start, int(x[2])-start) for x in isoform_dat]
                from fasta_manip import extract_ranges
                ref_seq = extract_ranges(ref_seq, ranges)
            if (adj_dir or translate) and strand == '-':
                ref_seq = ref_seq.reverse_complement()
            ## translate sequence if translate flag raised AND feature is CDS
            if translate:
                if feature == "CDS" and not complete:
                    ref_seq = ref_seq.translate(to_stop = True)
                else:
                    print("Translation is only possible when the selected feature is 'CDS' and the flag 'complete' is not raised.")
            seq_name = "Col-0_ref|{}|{}|{}".format(gene, feature, isoform) + \
                       (('|' + ("domain" if not domain else domain) + f"|{i+1}") if domain_f else '') + \
                       ("|complete" if complete else '') + \
                       ("|revcomp" if adj_dir and strand == '-' else '')
            seqs_to_write[seq_name] = ref_seq
            ## for by_gene
            if by_gene:
                overlap_ranges = []
                overlap_seq_names = []
                for logged_ranges, logged_seq_names in seq_ranges.items():
                    if has_any_overlap(ranges, logged_ranges):
                        overlap_ranges.append(logged_ranges)
                        overlap_seq_names.extend(logged_seq_names)
                if overlap_ranges:
                    for logged_ranges in overlap_ranges:
                        del(seq_ranges[logged_ranges])
                    ranges = merge_ranges(*overlap_ranges)
                seq_ranges[tuple(sorted(ranges))] = seq_ranges.get(tuple(sorted(ranges)), []) + [seq_name] + overlap_seq_names
            else:
                seq_ranges[tuple(sorted(ranges))] = seq_ranges.get(tuple(sorted(ranges)), []) + [seq_name]
        from fasta_manip import dict_to_fasta
        dict_to_fasta(seqs_to_write, fasta_out)
        fasta_out_l.append(fasta_out)
    fasta_out_final = os.path.join(out_dir, gene + "_ref_" + feature + \
                                   ("_complete" if complete else '') + \
                                   (('_' + ("domain" if not domain else domain)) if domain_f else '') + \
                                   ("_protein" if (translate and (feature=="CDS") and not complete) else '')+\
                                   ".fasta")
    if fasta_out_l:
        if fasta_out_l[0] != fasta_out_final:
            from file_manip import cat_files
            cat_files(sorted(fasta_out_l), fasta_out_final)
            for fasta_out in fasta_out_l:
                os.remove(fasta_out)
        if by_gene:
            isoform_seqs = fasta_to_dict(fasta_out_final)
            final_seqs = {}
            i = 0
            for ranges, seq_names in sorted(seq_ranges.items()):
                seq_name_l = seq_names[0].split('|')
                seq_name_l[3] = ','.join(f'{r[0]}-{r[1]}' for r in ranges)
                if domain_f:
                    seq_name_l[5] = str(i + 1)
                seq_name = '|'.join(seq_name_l)
                final_seqs[seq_name] = isoform_seqs[seq_names[0]]
                i += 1
            dict_to_fasta(final_seqs, fasta_out_final)
        print("Sequences were successfully written to: {}".format(fasta_out_final))
    elif not fasta_out_l:
        f = open(fasta_out_final, "w+")
        f.write('')
        f.close()
        print("{} is an empty file".format(fasta_out_final))

def get_ref_by_range(chrom, start, end, out_dir, encoding="utf-8",
                     ref_fasta_files=ref_fasta):
    fasta_out = os.path.join(out_dir, "chr{}_p{}-{}_ref.fasta".format(chrom, start + 1, end))
    get_ref_raw(chrom_pref + chrom, start, end, fasta_out, encoding=encoding, ref_fasta_files=ref_fasta_files)
    print("Sequences were successfully written to: {}".format(fasta_out))

## output is...0-indexed, I think...it seems like it follows .bed conventions...
def get_domain_in_genome_coords(gene, domain, domain_f, out_dir, pos_type="aa", isoform='',
                                bed=bed_path, encoding="utf-8", start_inc = True, end_inc = True,
                                qname_dname=("name", "domain name"), qstart_qend=("start", "end")):
    qname, dname = qname_dname
    qstart, qend = qstart_qend
    with open(domain_f, 'r') as f:
        domain_raw = [x[:-1].split('\t') for x in f.readlines()]
    domain_header = domain_raw[0]
    domain_data = [x for x in domain_raw[1:] if len(domain_raw) > 1 and len(x) == len(domain_header) and \
                   (((not domain) or domain == x[domain_header.index(dname)]) and \
                    (gene if not isoform else isoform) in x[domain_header.index(qname)])]
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
        cds_dat = grep_bedmerge((gene if not isoform else isoform), bed, "CDS", encoding, out_dir)["data"]
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

# ## get sequence info of genes
# def get_1001pseudogenome(accs, gene, feature, out_dir, bed=bed_path, encoding="utf-8",
#                          ref_file=data_path, num_col=0, name_col=1, delim='', header=True):
#     accs = [acc for acc in accs if acc]
#     if str(accs[0]).isdigit():
#         accs = dict(get_acc_name(accs, ref_file=ref_file, num_col=num_col, name_col=name_col,
#                                  delim=delim, header=header))
#     else:
#         accs = dict(get_acc_num(accs, ref_file=ref_file, num_col=num_col, name_col=name_col,
#                                 delim=delim, header=header))
#     ## get chrom, start, end of gene's features
#     grep = subprocess.Popen(("grep", gene, bed), stdout=subprocess.PIPE)
#     awk = subprocess.Popen(("awk", "$8 == \"" + feature + "\""), stdin=grep.stdout, stdout=subprocess.PIPE)
#     merged = [entry.split('\t') for entry in
#               subprocess.check_output(("bedtools", "merge"),
#                                       stdin=awk.stdout).decode(encoding).split('\n') if entry]
#     ## write bed file of regions used (0-indexed)
#     f = open(os.path.join(out_dir, gene + '_' + feature + ".bed"), "w+")
#     f.write('\n'.join(['\t'.join(x) for x in merged]))
#     f.close()
#     ## extract sequences
#     chrom = merged[0][0]
#     start = min([int(x[1]) for x in merged]) + 1
#     end = max([int(x[2]) for x in merged])
#     ## get fasta file of sequence data
#     fasta_out = os.path.join(out_dir, gene + '_' + feature + ".fasta")
#     subprocess.run(("curl", "-o", fasta_out,
#                     make_1001pseudogenome_api(chrom, start, end, accs.keys())))
#     ## trim sequences in file
#     from fasta_manip import extract_feature
#     seqs_trimmed = extract_feature(fasta_out, bed, tuple(accs.keys())[0], gene, start - 1, fasta_out, feature,
#                                    encoding = encoding)
#     from fasta_manip import fasta_to_dict
#     seqs = fasta_to_dict(fasta_out)
#     seqs_fin = {k.split('|')[3] + '|' + accs[k.split('|')[3]] + '|' + gene + '|' + feature:v \
#                 for k,v in seqs.items()}
#     from fasta_manip import dict_to_fasta
#     dict_to_fasta(seqs_fin, fasta_out)
#     print("Sequences were successfully written to: {}".format(fasta_out))

    
# def get_ref(gene, feature, out_dir, bed=bed_path, encoding="utf-8",
#             ref_fasta_files=ref_fasta, include_within_features=False):
#     ## get chrom, start, end of gene's features
#     grep = subprocess.Popen(("grep", gene, bed), stdout=subprocess.PIPE)
#     awk = subprocess.Popen(("awk", "$8 == \"" + feature + "\""), stdin=grep.stdout, stdout=subprocess.PIPE)
#     merged = [entry.split('\t') for entry in
#               subprocess.check_output(("bedtools", "merge"),
#                                       stdin=awk.stdout).decode(encoding).split('\n') if entry]
#     ## write bed file of regions used (0-indexed)
#     f = open(os.path.join(out_dir, gene + '_' + feature + ".bed"), "w+")
#     f.write('\n'.join(['\t'.join(x) for x in merged]))
#     f.close()
#     ## extract sequences from fasta file
#     chrom = merged[0][0]
#     start = min([int(x[1]) for x in merged])
#     end = max([int(x[2]) for x in merged])
#     ## get fasta file of sequence data
#     fasta_out = os.path.join(out_dir, gene + '_' + feature + ("_includewithinfeatures" if include_within_features else '') + ".fasta")
#     from fasta_manip import fasta_to_dict
#     ref_seq = list(fasta_to_dict(ref_fasta_files[chrom]).values())[0][start:end]
#     if not include_within_features:
#         ranges = [(int(x[1])-start, int(x[2])-start) for x in merged]
#         from fasta_manip import extract_ranges
#         ref_seq = extract_ranges(ref_seq, ranges)
#     from fasta_manip import dict_to_fasta
#     dict_to_fasta({"Col-0_ref|{}|{}".format(gene, feature): ref_seq}, fasta_out)
#     print("Sequences were successfully written to: {}".format(fasta_out))    
