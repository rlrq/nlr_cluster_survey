#!/bin/bash

###################
##  COMMON VARS  ##
##    + SETUP    ##
###################

## directories
this_script_dir=$(realpath $(dirname "$0"))
scripts_dir="${this_script_dir}/scripts" ## make sure scripts directory is in same folder as this script, otherwise this should be updated to reflect location of scripts
data_dir="${this_script_dir}/data" ## data directory; to be updated accordingly
db_dir="${this_script_dir}/db" ## blastdb/rpsblastdb director; to be updated accordingly or removed if not needed
results_dir="${this_script_dir}/results" ## output directory; to be modified if so desired
novikova_dir="${data_dir}/Novikova_etal" ## directory containing Novikova et al. fastq files; to be updated accordingly

mkdir -p "${results_dir}"

## temporary files
tmp_f="${results_dir}/tmp.txt"
tmp_f2="${results_dir}/tmp2.txt"

## VdW contig sequence data
vdw_contigs="${data_dir}/VdW.fasta" ## assembly FASTA files downloaded from 2blades, cat into a single file

## A. thaliana data (to be updated accordingly)
arath_all="${data_dir}/a_thaliana_all.fasta" ## all A. thaliana chromosomes, cat into a single file (chromosome prefix 'Chr')
arath_gff="${data_dir}/TAIR10_GFF3_genes.gff" ## downloaded from TAIR
arath_bed="${data_dir}/TAIR10_GFF3_genes.bed" ## to be generated using the following:
## generate BED file from GFF3 annotations
sed 's/; /;/g' ${arath_gff} | tr ' ' '=' > ${tmp_f}
${scripts_dir}/gff_to_bed.sh ${tmp_f} ${tmp_f2}
grep -v '^#' ${tmp_f2} | awk 'NR > 1 {print}' > ${arath_bed}

## A. lyrata data (to be updated accordingly)
araly_aa="${data_dir}/Araly1_GeneModels_FilteredModels6_aa.fasta" ## downloaded from genome.jgi.doe.gov
araly_fa="${data_dir}/Araly1_assembly_scaffolds.fasta" ## downloaded from genome.jgi.doe.gov
araly_gff="${data_dir}/Araly1_GeneModels_FilteredModels6.gff" ## downloaded from genome.jgi.doe.gov
araly_bed="${data_dir}/Araly1_GeneModels_FilteredModels6.gff.bed" ## to be generated using the following:
## generate BED file from GFF3 annotations
sed 's/; /;/g' ${araly_gff} | tr ' ' '=' > ${tmp_f}
${scripts_dir}/gff_to_bed.sh ${tmp_f} ${tmp_f2}
grep -v '^#' ${tmp_f2} | awk 'NR > 1 {print}' > ${araly_bed}
rm ${tmp_f} ${tmp_f2}
## positions of A. lyrata NLR
araly_nlr_txt=/path/to/Araly1_NLR_full.txt ## download from github repo and update accordingly

## A. lyrata read data from Novikova et al.
novikova_araly=/path/to/ploidy2_alyrata.tsv ## download from the github repo and update accordingly
sras=( $( cut -f17 ${novikova_araly} | head -18 | tail -17 | tr '\n' ' ' ) )
ncbi_sra_download_dir=/path/to/ncbi/public/sra/ ## to be updated accordingly
## If you haven't yet downloaded Novikova et al.'s fastq files, execute the following (files will be downloaded into ${novikova_dir} (defined at start of this script file); this requires 'fastq-dump' from the SRA toolkit
for sra in ${sras[@]}; do
    prefetch ${sra}
    mv ${ncbi_sra_download_dir}/${sra}.sra ${novikova_dir} ## update accordigly
    fastq-dump --split-files --gzip -O ${novikova_dir} ${novikova_dir}/${sra}.sra
done

## accession/gene/protein/cluster/domain names
accs_map="${data_dir}/accID_AL70.tsv" ## supplementary table 1, in tsv format (w/ header)
genes_f="${data_dir}/nlr164_extended_wPred.tsv" ## supplementary table 2, in tsv format (w/ header)
nlr164="${data_dir}/nlr164_NLR.txt" ## fourth column of supplementary table 1, only rows with 'y' in 'original.bait' column (w/o header)
nlraraly="${data_dir}/Araly1_NLR.txt" ## first column of supplementary table 7 (w/o header)
cdd_v_f="${data_dir}/cdd.versions.tsv" ## location of file mapping pssmid to pssm names, downloadable from same source as CDD database

## databases
db_cdd="${db_dir}/Cdd" ## location of installed CDD database; to be updated accordingly (please set up according to instructions from CDD's README file)
db_al="${db_dir}/VdW/VdW" ## location of VdW contig database; to be updated accordingly
## NOTE: please set up db_al using the following line:
makeblastdb -in "${vdw_contigs}" -dbtype nucl -out "${db_al}"

## executables
get_seq="${scripts_dir}/get_seqs/get_seqs.sh -b ${arath_bed}" ## get_seqs execution
## NOTE: ${scripts_dir}/get_seqs/scripts/get_seqs_functions.py needs to be updated with path to FASTA file for each A. thaliana chromosome
fast_tree="${scripts_dir}/FastTree" ## FastTree binary; to be updated accordingly
yeppp_base=/path/to/yeppp/yeppp-1.0.0/ ## Yeppp! mathematical library, used for CNVnator; to be updated accordingly
root_dir=/path/to/root/dir/ ## ROOT (https://root.cern/install/), used for CNVnator; to be updated accordingly

## other variables
declare -A domains_a=( ['NB-ARC']='gnl|CDD|307194' ['TIR']='gnl|CDD|214587' ) ## domains + pssmid

## source functions
source ${scripts_dir}/run_blast6.sh ## has function that adds header row to outfmt 6 files


######################
##  GET TIR/CC/NBS  ##
##     SEQUENCES    ##
######################

## variables and directory
cdd_dir="${results_dir}/cdd"
domain_seq="${results_dir}/domain_seq"
mkdir -p "${cdd_dir}" "${domain_seq}"
nlr164_aa="${results_dir}/nlr164_protein.fasta"
nlr164_cdd_pref="${cdd_dir}/nlr164_cdd"

## get A. thaliana NLR protein sequences
$get_seq -g "${nlr164}" -a ref -d "${results_dir}" --out "${nlr164_aa}" --translate -f CDS --no-bed

## find domains in A. thaliana NLR proteins
echo "Running RPS-BLAST"
run_blast6 rpsblast+ -db "${db_cdd}" -query "${nlr164_aa}" -out "${nlr164_cdd_pref}.tsv" -outfmt "6 qseqid sseqid pident length qstart qend"

## view types of domains found and map pssm id to pssm name
python3 ${scripts_dir}/summarise_domains.py "${nlr164_cdd_pref}.tsv" "${nlr164_cdd_pref}_isoform.tsv" "${nlr164_cdd_pref}_gene.tsv" "${nlr164_cdd_pref}_domain.tsv" "${cdd_v_f}"

## extract domains from A. thaliana NLRs
for domain in ${!domains_a[@]}; do
    pssm=${domains_a["${domain}"]}
    echo ${domain} ${pssm}
    printf "$(head -1 ${nlr164_cdd_pref}.tsv)\n$(grep ${pssm} ${nlr164_cdd_pref}.tsv)\n" > ${tmp_f}
    cds_complete="${domain_seq}/nlr164_${domain}_ref_CDS_complete.fa"
    cds_only="${domain_seq}/nlr164_${domain}_ref_CDS.fa"
    $get_seq -g "${nlr164}" -a ref -d "${domain_seq}" --qname-dname "('qseqid', 'domain')" --qstart-qend "('qstart', 'qend')" -f CDS --adjust-dir --domain "${domain}" --domain-file "${tmp_f}" --no-bed --complete --out "${cds_complete}" --by-gene
    $get_seq -g "${nlr164}" -a ref -d "${domain_seq}" --qname-dname "('qseqid', 'domain')" --qstart-qend "('qstart', 'qend')" -f CDS --adjust-dir --domain "${domain}" --domain-file "${tmp_f}" --no-bed --out "${cds_only}" --by-gene
    rm "${tmp_f}"
done

## rename extension
for f in ${domain_seq}/*.fa; do
    mv ${f} ${f}sta
done

## same as above, but for newly predicted genes (to be executed only after "completeness check")
nlr_pred='AT4G16930,AT5G38344,AT5G45440,AT5G48775' ## update accordingly if so desired
nlr_pred_aa="${results_seq}/nlr164_6909pred_protein.fasta"
nlr_pred_cdd_pref="${cdd_dir}/nlr164_6909pred_cdd"
$get_seq -g "${nlr_pred}" -a ref -d "${results_dir}" --out "${nlr_pred_aa}" --translate -f CDS --no-bed
## find domains
echo "Running RPS-BLAST"
run_blast6 rpsblast+ -db "${db_cdd}" -query "${nlr_pred_aa}" -out "${nlr_pred_cdd_pref}.tsv" -outfmt "6 qseqid sseqid pident length qstart qend"
## view types of domains found (map pssm id to pssm name)
python3 ${scripts_dir}/summarise_domains.py "${nlr_pred_cdd_pref}.tsv" "${nlr_pred_cdd_pref}_isoform.tsv" "${nlr_pred_cdd_pref}_gene.tsv" "${nlr164_cdd_pref}_domain.tsv" "${cdd_v_f}"
## extract domain
for domain in ${!domains_a[@]}; do
    pssm=${domains_a["${domain}"]}
    echo ${domain} ${pssm}
    printf "$(head -1 ${nlr_pred_cdd_pref}.tsv)\n$(grep ${pssm} ${nlr_pred_cdd_pref}.tsv)\n" > ${tmp_f}
    cds_complete="${domain_seq}/nlr164_${domain}_6909pred_ref_CDS_complete.fa"
    cds_only="${domain_seq}/nlr164_${domain}_6909pred_ref_CDS.fa"
    $get_seq -g "${nlr_pred}" -a ref -d "${domain_seq}" --qname-dname "('qseqid', 'domain')" --qstart-qend "('qstart', 'qend')" -f CDS --adjust-dir --domain "${domain}" --domain-file "${tmp_f}" --no-bed --complete --out "${cds_complete}" --by-gene
    $get_seq -g "${nlr_pred}" -a ref -d "${domain_seq}" --qname-dname "('qseqid', 'domain')" --qstart-qend "('qstart', 'qend')" -f CDS --adjust-dir --domain "${domain}" --domain-file "${tmp_f}" --no-bed --out "${cds_only}" --by-gene
    rm "${tmp_f}"
done
## rename extension
for f in ${domain_seq}/*.fa; do
    mv ${f} ${f}sta
done


###################
##  GET DOMAINS  ##
##   IN AL70 &   ##
##   A. lyrata   ##
###################

domain_vdw="${results_dir}/domain_vdw"
domain_alyrata="${results_dir}/domain_alyrata"
mkdir "${domain_vdw}" "${domain_alyrata}"

## Araly based on Ya-Long Guo's 2011 paper's A. lyrata NLR protein IDs
lyrata_ylg_aa="${domain_alyrata}/araly1_NLR_protein.fasta"
lyrata_ylg_cdd_pref="${domain_alyrata}/araly1_NLR_protein.cdd"
## get lyrata protein sequences from YLG2011
python3 -c "import sys; sys.path.append('${scripts_dir}'); from fasta_manip import *; f = open('${nlraraly}', 'r'); nlr_id = [x[:-1] for x in f.readlines()]; f.close(); proteins = fasta_to_dict('${araly_aa}'); proteins_nlr = {seq_id: seq for seq_id, seq in proteins.items() if seq_id.split('|')[2] in nlr_id}; dict_to_fasta(proteins_nlr, '${lyrata_ylg_aa}')"
## CDD search using YLG2011
run_blast6 rpsblast+ -db "${db_cdd}" -query "${lyrata_ylg_aa}" -out "${lyrata_ylg_cdd_pref}.tsv" -outfmt "6 qseqid sseqid pident length qstart qend"
python3 -c "import sys; sys.path.append('${scripts_dir}'); from map_pssmid import *; map_pssmid('${lyrata_ylg_cdd_pref}.tsv', 1, '${lyrata_ylg_cdd_pref}_reformat.tsv', header = True, output_header = True, new_col = 'domain', cdd_versions_fname = '${cdd_v_f}')"
mv ${lyrata_ylg_cdd_pref}_reformat.tsv ${lyrata_ylg_cdd_pref}.tsv


## BLAST against VdW and superimpose over whole gene alignments to count genes :D
## -- of course, do some post-processing to filter out VdW contigs not identified by full gene recip_blast
for domain in ${!domains_a[@]}; do
    
    echo ${domain}
    
    echo "Processing nlr164"
    pref="nlr164_${domain}_ref_CDS"
    suff_vdw=vdw.blastn-dc_summary
    suff_alyrata=alyrata.blastn-dc_summary
    cols="6 qseqid sseqid pident length mismatch gaps qstart qend sstart send qlen slen evalue bitscore qframe sframe"
    run_blast6 blastn -outfmt "${cols}" -db "${db_al}" -query "${domain_seq}/${pref}_complete.fasta" -out "${domain_vdw}/${pref}_complete.${suff_vdw}.tsv"
    run_blast6 blastn -outfmt "${cols}" -db "${db_al}" -query "${domain_seq}/${pref}.fasta" -out "${domain_vdw}/${pref}.${suff_vdw}.tsv"
    
    echo "Processing YLG A. lyrata"
    ## get corresponding domains from A. lyrata (from genome and predicted mRNA)
    ## YLG2011, use earlier CDD search results to define domain boundaries of YLG2011 NLRs
    lyrata_ylg_domain_pref="${domain_alyrata}/araly1_${domain}_CDS"
    python3 -c "import sys; sys.path.append('${scripts_dir}'); fasta = '${araly_fa}'; gff_bed = '${araly_bed}'; pid_field = 'proteinId'; domain_pid_f = lambda x: x.split('|')[2]; from get_lyrata_functions import *; get_cds('${lyrata_ylg_domain_pref}.fasta', bed = gff_bed, fasta = fasta, domain_f = '${lyrata_ylg_cdd_pref}.tsv', domain = '${domain}', adjust_dir = True, protein_id_field = pid_field, domain_pid_f = domain_pid_f, restrict_pid_exact = True); get_cds('${lyrata_ylg_domain_pref}_complete.fasta', domain_f = '${lyrata_ylg_cdd_pref}.tsv', bed = gff_bed, fasta = fasta, domain = '${domain}', adjust_dir = True, complete = True, protein_id_field = pid_field, domain_pid_f = domain_pid_f, restrict_pid_exact = True)"
done

## reformat blast results to genome
Rscript -e "source('${scripts_dir}/reformat_for_domain_extraction.R'); reformat('nlr164', genes_f = '${genes_f}', vdw_dir = '${domain_vdw}')"

## output in ${domain_vdw}

####################
##     EXTRACT    ##
##  DOMAINS FROM  ##
##   ANNALENA70   ##
####################

python3 -c "import sys; sys.path.append('${scripts_dir}'); cols = ['qseqid', 'sseqid', 'pident', 'length', 'qframe', 'sframe',  'gene', 'accID', 'tigID', 'contig', 'contig.len', 'hit.start', 'hit.end', 'geneID', 'domain', 'group']; from extract_domain import *; extract_domains_multi(coi_combos = [['nlr164', 'TIR'], ['nlr164', 'NB-ARC']], vdw_contigs = '${vdw_contigs}', cols = cols, make_fname = (lambda cluster, domains: '${domain_vdw}/{}_{}_ref.annaLena70.blastn-dc_summary_reformat.tsv'.format(cluster, '-'.join(domains))), equivalents = {'pident': 'X..identity', 'length': 'alignment.length', 'qframe': 'q..frame', 'sframe': 's..frame'}, make_merged_fout = (lambda cluster, d: '${domain_vdw}/{}_col0-{}_merged.tsv'.format(cluster, d)), make_domain_seq_fout = (lambda cluster, d: '${domain_vdw}/{}_col0-{}_predicted.fasta'.format(cluster, d)))"


##########################
##  COMPLETENESS CHECK  ##
##   CHECK 6909 GENES   ##
##      ANNALENA70      ##
##########################

completeness_dir="${results_dir}/completeness"
mkdir -p "${completeness_dir}"
arath_pred_bed="${completeness_dir}/6909.predicted.bed"
arath_pred_gff_bed="${completeness_dir}/6909.predicted.GFF3.bed"
arath_pred_gff_bed_cds="${completeness_dir}/6909.predicted.GFF3_CDS.bed"
arath_pred_tsv="${completeness_dir}/6909.predicted.best2_hits.tsv"

## To check completeness of AL70 predictions (sensitivity check)
for domain in ${!domains_a[@]}; do
    echo ${domain}
    pred_6909="${domain_vdw}/nlr164_col0-${domain}_predicted_6909.fasta"
    pred_6909_b6="${completeness_dir}/6909.${domain}_predicted.tsv"
    python3 -c "import sys; sys.path.append('${scripts_dir}'); from fasta_manip import *; seqs = fasta_to_dict('${domain_vdw}/nlr164_col0-${domain}_predicted.fasta'); output = {seq_id: seq for seq_id, seq in seqs.items() if '6909.tig' in seq_id}; dict_to_fasta(output, '${pred_6909}')"
    run_blast6 blastn -query "${pred_6909}" -subject "${arath_all}" -out "${pred_6909_b6}" -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send qlen slen evalue bitscore qframe sframe'
done

## manipulate a little with completeness_nlr164.R to get bed format (see line with alcolb6) (alt, reproduced here in a single line lol)
Rscript -e "source('${scripts_dir}/completeness.R'); b6_to_bed('${completeness_dir}', '${arath_pred_bed}', pattern = '_predicted.tsv', cds_len_src = '${domain_vdw}', cds_len_f_pattern = '_merged.tsv')"

## intersect predicted 6909 with respective GFF3
bedtools intersect -wo -a ${arath_pred_bed} -b ${arath_bed} > ${arath_pred_gff_bed}
## from intersected bed file, get only CDS and merge
awk -F'\t' '{if ($19==CDS) {split($20,a,"."); split(a[1],g,"="); print $1"\t"$2"\t"$3"\t"$4"\t"g[2]}}' ${arath_pred_gff_bed} | bedtools sort | bedtools merge -c 4,5 -o distinct | head > ${arath_pred_gff_bed_cds}

## manipulate some more with completeness_nlr164.R, results would be easier to see with excel or libreoffice calc
Rscript -e "source('${scripts_dir}/completeness.R'); summarise_pred_genes('${arath_pred_gff_bed}', '${arath_pred_gff_bed_cds}', '${arath_pred_tsv}', '${nlr164}')"
## the output is ${arath_pred_tsv}

## go back up and repeat getting CDS and CDS complete for these new genes before proceeding


#######################
##        BUILD      ##
##    DOMAIN TREES   ##
##  W A.LYRATA SEQS  ##
#######################

# accs_map=/mnt/chaelab/rachelle/data/anna_lena/accID_AL70.tsv
alignment_in_dir="${results_dir}/alignment_input"
alignment_dir="${results_dir}/alignment"
tree_dir="${results_dir}/tree"
threads=10
mkdir -p "${alignment_in_dir}" "${alignment_dir}" "${tree_dir}"
## MAFFT MY OLD FRIEND (into ${results_dir}/alignment and ${results_dir}/tree)
for domain in "${!domains_a[@]}"; do
    col0_araly_cds="${alignment_in_dir}/nlr164_col0-Alyrata_${domain}_CDS.fasta"
    col0_araly_cds_mafft="${alignment_dir}/nlr164_col0-Alyrata_${domain}_CDS_mafft.fa"
    col0_araly_cds_complete="${alignment_in_dir}/nlr164_col0-Alyrata_${domain}_CDS_complete.fasta"
    col0_araly_cds_complete_mafft="${alignment_dir}/nlr164_col0-Alyrata_${domain}_CDS_complete_mafft.fa"
    aligned_fa="${alignment_dir}/nlr164_col0-AL70-Alyrata_${domain}_mafft.fa"
    aligned_fa_pruned="${alignment_dir}/nlr164_col0-AL70-Alyrata_${domain}_mafft_pruned.fa"
    tree_nwk="${tree_dir}/nlr164_col0-AL70-Alyrata_${domain}_mafft_ML.nwk"
    ## align CDS-only from Col-0 and Alyrata
    echo "Aligning ${domain} from CDS-only Col-0 and A. lyrata"
    cat "${domain_seq}/nlr164_${domain}_ref_CDS.fasta" "${domain_seq}/nlr164_${domain}_6909pred_ref_CDS.fasta" "${domain_alyrata}/araly1_${domain}_CDS.fasta" > "${col0_araly_cds}"
    mafft --thread ${threads} "${col0_araly_cds}" > "${col0_araly_cds_mafft}"
    ## align CDS complete from Col-0 and Alyrata to CDS-only
    echo "Aligning complete CDS from Col-0 and A. lyrata to CDS-only from Col-0 and A. lyrata"
    cat "${domain_seq}/nlr164_${domain}_ref_CDS_complete.fasta" "${domain_seq}/nlr164_${domain}_6909pred_ref_CDS_complete.fasta" "${domain_alyrata}/araly1_${domain}_CDS_complete.fasta" > "${col0_araly_cds_complete}"
    mafft --thread ${threads} --add "${col0_araly_cds_complete}" "${col0_araly_cds_mafft}" > "${col0_araly_cds_complete_mafft}"
    ## align all other sequences from VdW and Alyrata
    echo "Aligning AL70"
    mafft --thread ${threads} --adjustdirectionaccurately --large --add "${domain_vdw}/nlr164_col0-${domain}_predicted.fasta" "${col0_araly_cds_complete_mafft}" > "${aligned_fa}"
    ## duplicate the newly discovered 6909 genes/pseudogenes and rename them for easier gene assignment, remove Col-0 and A. lyrata sequences not marked as 'complete'
    python3 -c "import re; import sys; sys.path.append('${scripts_dir}'); from fasta_manip import *; from data_manip import *; accs = [x.split('\t')[0] for x in open('${accs_map}', 'r').readlines()][1:]; alignment = {seq_id: seq for seq_id, seq in fasta_to_dict('${aligned_fa}').items() if re.search('^(_R_)?(' + '|'.join(accs) + ')\.tig|complete(\|revcomp)?$', seq_id)}; covered_genes = [x.split('|')[1] for x in alignment.keys() if 'Col-0_ref' in x]; f = open('${completeness_dir}/6909.predicted.best2_hits.tsv', 'r'); to_dupl_dat = [x[:-1].split('\t') for x in f.readlines() if len(x) > 1]; f.close(); get = make_custom_get(to_dupl_dat[0]); to_dupl_dat = to_dupl_dat[1:]; to_dupl_seqs = {get(entry, 'seqid'): get(entry, 'gene') for entry in to_dupl_dat if get(entry, 'nlr164') == 'FALSE' and get(entry, 'pident') == 100 and get(entry, 'gene') not in covered_genes}; to_dupl_seqs = {**to_dupl_seqs, **{'_R_' + get(entry, 'seqid'): get(entry, 'gene') for entry in to_dupl_dat}}; dupl_seqs = {f'Col-0_ref|{to_dupl_seqs[seqid]}|CDS|{to_dupl_seqs[seqid]}.1|${domain}|1|complete': seq for seqid, seq in alignment.items() if seqid in to_dupl_seqs}; alignment = {**alignment, **dupl_seqs}; dict_to_fasta(alignment, '${aligned_fa_pruned}')"
    ## make tree
    echo "Making tree for ${domain}"
    ${fast_tree} -gtr -nt < "${aligned_fa_pruned}" > "${tree_nwk}"
done



#######################
##  ASSIGN GENES IN  ##
##  VdW & A. lyrata  ##
#######################
## Assign closest reference gene in tree as tentative identity of each predicted domain sequence
## ## manually root the NB-ARC tree using TIR/CC genes as the most basal divide
##   see find_dist_to_ref.py
pred_dir="${results_dir}/predicted_identity"
mkdir -p "${pred_dir}"
python3 -c "import sys; sys.path.append('${scripts_dir}'); from find_dist_to_ref import *; get_closest_ref_nlr164('NB-ARC', rooted = True, coi_f = '${genes_f}', out_dir = '${pred_dir}', t_nwk = '${tree_dir}/nlr164_col0-AL70-Alyrata_NB-ARC_mafft_ML.nwk', prefix = 'nlr164'); get_closest_ref_nlr164('TIR', coi_f = '${genes_f}', out_dir = '${pred_dir}', t_nwk = '${tree_dir}/nlr164_col0-AL70-Alyrata_TIR_mafft_ML.nwk', prefix = 'nlr164')"


##########################
##  SEQUENCE DIVERSITY  ##
##########################

## remap clusters
# genes_f=${results_dir}/nlr164_final_wPred_newClust.tsv
diversity_dir="${results_dir}/seq_diversity"
mkdir -p "${diversity_dir}"
python3 -c "import sys; sys.path.append('${scripts_dir}'); from seq_conservation import *; calculate_pi_multi(fout_dir = '${diversity_dir}', make_seq_fname = (lambda domain: f'${alignment_dir}/nlr164_col0-AL70-Alyrata_{domain}_mafft.fa'), prefix = 'nlr164_AL70', make_pred_fname = (lambda domain: f'${pred_dir}/nlr164_AL70_{domain}_predictedIdentity.txt'))"
Rscript -e "library(dplyr); genes <- read.table('${genes_f}', sep = '\t', stringsAsFactors = FALSE, header = TRUE) %>% select(gene, cluster); for (fname in list.files('${diversity_dir}', full.names = TRUE)){print(fname); dat <- read.table(fname, sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>% rename(cluster_old = cluster) %>% left_join(genes, by = 'gene'); write.table(dat, fname, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)}"


#####################
##  INTRA & INTER  ##
##     CLUSTER     ##
##    DISTANCES    ##
#####################

# genes_f=${results_dir}/nlr164_final_wPred_newClust.tsv
distances_dir="${results_dir}/seq_distances"
mkdir -p "${distances_dir}"
python3 -c "import sys; sys.path.append('${scripts_dir}'); from within_cluster_dists import *; get_all_pairwise_distance_multi(t_dir = '${tree_dir}', predicted_identity_dir = '${pred_dir}', cluster_gene_f = '${genes_f}', out_dir = '${distances_dir}')"



# #########################
# ##  INTRA- and INTER-  ##
# ##  CLUSTER DISTANCE   ##
# #########################

# nlr_bed=${results_dir}/nlr164_final_wPred_GFF3.bed
# grep -f ${results_dir}/nlr164/nlr164_final_wPred.txt ${arath_bed} | awk '$8=="gene"||$8=="pseudogene"' > ${nlr_bed}

# ## see physical_distance.R for analysis


#################
##  A. lyrata  ##
##     CNV     ##
#################

## align to reference A. lyrata genome
araly_cnv_dir=${results_dir}/araly_cnv
mkdir -p ${araly_cnv_dir}
mkdir ${araly_cnv_dir}/sam ${araly_cnv_dir}/bam
bwa index ${araly_fa} ## index reference
samtools faidx ${araly_fa}
threads=20
for sra in ${sras[@]}; do
    echo "Working on ${sra}"
    fq_pref=${novikova_dir}/${sra}
    sam_f=${araly_cnv_dir}/sam/${sra}.sam
    bam_f=${araly_cnv_dir}/bam/${sra}.bam
    bam_sorted_f=${araly_cnv_dir}/bam/${sra}.sorted.bam
    echo "Running bwa"
    bwa mem -v 2 -t ${threads} ${araly_fa} ${fq_pref}_1.fastq.gz ${fq_pref}_2.fastq.gz > ${sam_f}
    echo "Converting SAM file to BAM file"
    samtools view -S -b ${sam_f} > ${bam_f}
    echo "Sorting BAM file"
    samtools sort ${bam_f} -o ${bam_sorted_f}
    echo "Indexing BAM file"
    samtools index ${bam_sorted_f}
    rm ${sam_f} ${bam_f}
done

## run CNVnator
## get the environment ready...
export YEPPPLIBDIR=${yeppp_base}/binaries/linux/x86_64
export YEPPPINCLUDEDIR=${yeppp_base}/library/headers
export LD_LIBRARY_PATH=$YEPPPLIBDIR:$LD_LIBRARY_PATH
cd ${root_dir}
source ./config/thisroot.sh
cnvnator_dir=${araly_cnv_dir}/cnvnator
cd ${cnvnator_dir}
ln -s /path/to/samtools ${cnvnator_dir}/samtools ## update /path/to/samtools appropriately
ln -s /path/to/htslib ${cnvnator_dir}/htslib ## update /path/to/samtools appropriately

## run cnvnator
araly_nlr_chrom=$(cut -f1 ${araly_nlr_txt} | sort | uniq | tr '\n' ' ')
## extract individual chromosome sequences into .fa files (because cnvnator -his ... requires it :()
araly_nlr_chrom_a=( ${araly_nlr_chrom} )
for chrom in ${araly_nlr_chrom_a[@]}; do
    python3 -c "import sys; sys.path.append('${scripts_dir}'); from fasta_manip import fasta_to_dict, dict_to_fasta; dict_to_fasta({k: v for k, v in fasta_to_dict('${araly_fa}').items() if k == '${chrom}'}, '${cnvnator_dir}/${chrom}.fa')"
done
## seems like this has got to be done for each sample separately :(
binsize=500
for sra in ${sras[@]}; do
    echo "Processing ${sra}"
    cnvnator_root_chrom=$( echo "-root ${cnvnator_dir}/aralyCNV_${sra}.root -chrom ${araly_nlr_chrom}" )
    ## EXTRACT from BAM
    echo "Extracting from BAM file"
    cnvnator ${cnvnator_root_chrom} -tree ${araly_cnv_dir}/bam/${sra}.sorted.bam
    ## HISTOGRAM
    echo "Generating histogram"
    cnvnator ${cnvnator_root_chrom} -his ${binsize} -d ${cnvnator_dir}
    ## CALCULATE STATS
    echo "Calculating stats"
    cnvnator ${cnvnator_root_chrom} -stat ${binsize}
    ## RD SIGNAL PARTITIONING
    echo "Partitioning RD signal"
    cnvnator ${cnvnator_root_chrom} -partition ${binsize}
    ## CNV calling
    echo "Calling CNV"
    cnvnator ${cnvnator_root_chrom} -call ${binsize} > ${cnvnator_dir}/aralyCNV_called_${sra}.tsv
    ## convert output file of CNV calling to bed file (well...technically only the first 3 columns will adhere to the bed fomat...); also, append SRA number to end of each row
    fname=${cnvnator_dir}/aralyCNV_called_${sra}
    python3 -c "f = open('${fname}.tsv', 'r'); data = [x[:-1].split('\t') for x in f.readlines()]; f.close(); data = [ [x[1].split(':')[0]] + [str(int(x[1].split(':')[1].split('-')[0]) - 1)] + [x[1].split(':')[1].split('-')[1]] + [x[0]] + x[2:] + ['${sra}'] for x in data]; f = open('${fname}.bed', 'w+'); f.write('\n'.join(['\t'.join(x) for x in data])); f.close()"
done

## bedtools intersect output with positions of canonical A. lyrata NLRs
## BED1: all gff_bed entries for A. lyrata NLRs
araly_nlr_map=${araly_cnv_dir}/Araly1_NLR_full.map
araly_nlr_bed=${araly_cnv_dir}/Araly1_NLR_GFF3.bed
araly_nbs_complete_bed=${araly_cnv_dir}/araly_NB-ARC_complete.bed
araly_nbs_bed=${araly_cnv_dir}/araly_NB-ARC.bed
python3 -c "import re; gff_bed = [x[:-1].split('\t') for x in open('${araly_bed}', 'r').readlines()]; nlr_pid = set(x[:-1].split('\t')[-1] for x in open('${araly_nlr_txt}', 'r').readlines()); nlr_gid = set(re.search('(?<=name=)[^;]+', x[-1]).group(0) for x in gff_bed if re.search('proteinId', x[-1]) and re.search('(?<=proteinId=)\d+', x[-1]).group(0) in nlr_pid); nlr_entries = [x for x in gff_bed if re.search('(?<=name=)[^;]+', x[-1]).group(0) in nlr_gid]; f = open('${araly_nlr_bed}', 'w+'); f.write('\n'.join(['\t'.join(x) for x in nlr_entries])); f.close(); mapped = sorted(tuple(set(tuple([re.search('(?<=proteinId=)\d+', x[-1]).group(0), re.search('(?<=name=)[^;]+', x[-1]).group(0)[1:-1]]) for x in nlr_entries if re.search('proteinId', x[-1])))); f = open('${araly_nlr_map}', 'w+'); f.write('\n'.join(['\t'.join(x) for x in mapped])); f.close()"
## BED2: NB-ARC positions (complete)
python3 -c "import sys; sys.path.append('${scripts_dir}'); fasta = '${araly_fa}'; gff_bed = '${araly_bed}'; pid_field = 'proteinId'; domain_pid_f = lambda x: x.split('|')[2]; from get_lyrata_functions import *; get_cds('${tmp_f}', domain_f = '${lyrata_ylg_cdd_pref}.tsv', bed = gff_bed, fasta = fasta, domain = 'NB-ARC', adjust_dir = True, complete = True, protein_id_field = pid_field, domain_pid_f = domain_pid_f, restrict_pid_exact = True, bed_out = '${araly_nbs_bed_complete}', write_fasta = False)"
## BED3: NB-ARC positions (CDS only)
python3 -c "import sys; sys.path.append('${scripts_dir}'); fasta = '${araly_fa}'; gff_bed = '${araly_bed}'; pid_field = 'proteinId'; domain_pid_f = lambda x: x.split('|')[2]; from get_lyrata_functions import *; get_cds('${tmp_f}', domain_f = '${lyrata_ylg_cdd_pref}.tsv', bed = gff_bed, fasta = fasta, domain = 'NB-ARC', adjust_dir = True, complete = False, protein_id_field = pid_field, domain_pid_f = domain_pid_f, restrict_pid_exact = True, bed_out = '${cnvnator_dir}/../araly_NB-ARC.bed', write_fasta = False)"

cnv_nlr_bed=${cnvnator_dir}/aralyCNV_called_novikova14_NLR.bed
cnv_nbs_bed=${cnvnator_dir}/aralyCNV_called_novikova14_NB-ARC.bed
cnv_nbs_complete_bed=${cnvnator_dir}/aralyCNV_called_novikova14_NB-ARC_complete.bed
for sra in ${sras[@]}; do
    cnvnator_bed=${cnvnator_dir}/aralyCNV_called_${sra}.bed
    bedtools intersect -wa -wb -a ${cnvnator_bed} -b ${araly_nlr_bed} | bedtools sort | uniq >> ${cnv_nlr_bed}
    bedtools intersect -wa -wb -a ${cnvnator_bed} -b ${araly_nbs_bed} | bedtools sort | uniq >> ${cnv_nbs_bed}
    bedtools intersect -wa -wb -a ${cnvnator_bed} -b ${araly_nbs_complete_bed} | bedtools sort | uniq >> ${cnv_nbs_complete_bed}
done

## for parsing cnv_nbs_bed and cnv_nbs_complete_bed
cnv_nbs_summary=${cnvnator_dir}/aralyCNV_called_novikova14_NB-ARC_summary.tsv
## this one separates NB-ARC domains even when in same gene
python3 -c "import re; samples = ['$(echo ${sras[@]} | sed 's/ /\x27,\x27/g')']; gff_bed = [x[:-1].split('\t') for x in open('${cnv_nbs_bed}', 'r').readlines()]; stats_d = {domain: {cnv: set(tuple(x) for x in gff_bed if domain == x[-1] and x[3] == cnv) for cnv in ['deletion', 'duplication']} for domain in [x[-1] for x in gff_bed]}; [print(k, '\tdeletion:', len(set(x[11] for x in v['deletion'])), '\tduplication:', len(set(x[11] for x in v['duplication']))) for k, v in stats_d.items()]; output = [['proteinId', 'domainId', 'deletion', 'duplication']] + [[re.search('^\d+', k).group(0), k, ';'.join(set(','.join([x[11], x[5]]) for x in v['deletion'])), ';'.join(set(','.join([x[11], x[5]]) for x in v['duplication']))] for k, v in stats_d.items()]; f = open('${cnv_nbs_summary}', 'w+'); f.write('\n'.join(['\t'.join(x) for x in output])); f.close()"
