import os
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from map_pssmid import *

fname = sys.argv[1]
fout_isoform = sys.argv[2]
fout_gene = sys.argv[3]
fout_domain = sys.argv[4]
cdd_versions_fname = sys.argv[5]
gene_pattern = "AT.G\d{5}" if len(sys.argv) == 6 else sys.argv[6]
isoform_pattern = "AT.G\d{5}.\d+" if len(sys.argv) < 8 else sys.argv[7]

print(sys.argv)

map_pssmid(fname = fname, pssm_col = 1, fout = fname, cdd_versions_fname = cdd_versions_fname,
           header = True, output_header = True, new_col = 'domain')
with open(fname, 'r') as f:
    data = [x[:-1].split('\t') for x in f.readlines()]
header = data[0]
get = make_custom_get(header)
data = data[1:]

all_domains = [tuple(str(y) for y in x) for x in get(data, 'domain', 'sseqid')]
domains = sorted(set(all_domains))
domains = {x: all_domains.count(x) for x in domains}
domains = sorted(list(domains.items()), key = lambda x: x[1], reverse = True)
# print("Top 30: \n", head(domains, n = 30), "...")
f = open(fout_domain, "w+")
f.write('\n'.join(['\t'.join(["domain", "pssmid", "count"])] + ['\t'.join(dp + (str(c),)) for dp, c in domains]))
f.close()

data = sorted(data)
output_isoform = []
qsseq = get(data[0], "qseqid", "sseqid", "domain")
count = 0
for i, entry in enumerate(data):
    curr_qsseq = get(entry, "qseqid", "sseqid", "domain")
    if qsseq == curr_qsseq:
        count += 1
    else:
        output_isoform.append([re.search(gene_pattern, qsseq[0]).group(0),
                               re.search(isoform_pattern, qsseq[0]).group(0)] + qsseq[1:] + [str(count)])
        qsseq = curr_qsseq
        count = 1
    if i == len(data) - 1:
        output_isoform.append([re.search(gene_pattern, qsseq[0]).group(0),
                               re.search(isoform_pattern, qsseq[0]).group(0)] + qsseq[1:] + [str(count)])

## summarise by isoform
output_isoform_header = ["gene", "isoform", "pssmid", "domain", "count"]
f = open(fout_isoform, "w+")
f.write('\n'.join(['\t'.join([str(y) for y in x]) for x in [output_isoform_header] + output_isoform]))
f.close()

## summarise by max domain count among all isoforms of each gene
output_isoform_get = make_custom_get(output_isoform_header)
data = sorted(output_isoform, key = lambda x: output_isoform_get(x, "gene", "pssmid", "domain"))
output_gene = []
qsseq = output_isoform_get(data[0], "gene", "pssmid", "domain")
count = output_isoform_get(data[0], "count")
for i, entry in enumerate(data):
    curr_qsseq = output_isoform_get(entry, "gene", "pssmid", "domain")
    if qsseq == curr_qsseq:
        count = max(count, output_isoform_get(entry, "count"))
    else:
        output_gene.append(qsseq + [str(count)])
        qsseq = curr_qsseq
        count = output_isoform_get(entry, "count")
    if i == len(data) - 1:
        output_gene.append(qsseq + [str(count)])

output_gene_header = ["gene", "pssmid", "domain", "count"]
f = open(fout_gene, "w+")
f.write('\n'.join(['\t'.join([str(y) for y in x]) for x in [output_gene_header] + output_gene]))
f.close()
