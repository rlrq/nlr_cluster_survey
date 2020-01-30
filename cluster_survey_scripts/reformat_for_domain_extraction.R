library(stringr)
library(tidyverse)

thisDir <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    # Rscript
    if (length(match) > 0) {
        fname <- normalizePath(sub(needle, "", cmdArgs[match]))
    }
    # 'source'd via R console
    else {
        fname <- normalizePath(sys.frames()[[1]]$ofile)
    }
    return(str_extract(fname, "^.+?(?=[^/]+?$)"))
}

source(paste0(thisDir(), "str_manip.R"))


reformat <- function(name, genes_f, vdw_dir){
    genes <- read.table(genes_f, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    make_fname_al70 <- function(domain, feature) {return(paste0(vdw_dir, "/", name, "_", domain, "_ref_", feature, ".annaLena70.blastn-dc_summary.tsv"))}
    make_fout_al70 <- function(domain) {return(paste0(vdw_dir, "/", name, "_", domain, "_ref.annaLena70.blastn-dc_summary_reformat.tsv"))}
    reformat_for_domain_extraction_by_domain(genes, c("TIR", "NB-ARC"), make_fname_al70, make_fout_al70)
}


reformat_for_domain_extraction_by_domain <- function(genes, coi_domains, make_fname, make_fout,
                                                     uncategorised="singletons"){
    g <- genes$gene
    genes <- genes %>% mutate(cluster = as.character(cluster))
    rownames(genes) <- g
    for (domain in coi_domains){
        dat_cds <- read.table(make_fname(domain, "CDS"), sep = '\t',
                              header = TRUE, stringsAsFactors = FALSE) %>%
            mutate(contig = sseqid, contig.len = slen,
                   hit.start = sstart, hit.end = send,
                   gene = sapply(qseqid, gen_split_extract('\\|', 2)),
                   accID = sapply(sseqid, gen_split_extract('\\.', 1)),
                   tigID = sapply(sseqid, gen_split_extract('\\.', 2)),
                   geneID = gene, domain = c(paste0(domain, "_CDS")),
                   group = sapply(gene, function(x){if (! x %in% rownames(genes)) {uncategorised}
                                                    else {genes[x, "cluster"]}}))
        dat_complete <- read.table(make_fname(domain, "CDS_complete"), sep = '\t',
                                   header = TRUE, stringsAsFactors = FALSE) %>%
            mutate(contig = sseqid, contig.len = slen,
                   hit.start = sstart, hit.end = send,
                   gene = sapply(qseqid, gen_split_extract('\\|', 2)),
                   accID = sapply(sseqid, gen_split_extract('\\.', 1)),
                   tigID = sapply(sseqid, gen_split_extract('\\.', 2)),
                   geneID = gene, domain = c(domain),
                   group = sapply(gene, function(x){if (! x %in% rownames(genes)) {uncategorised}
                                                    else {genes[x, "cluster"]}}))
        write.table(rbind(dat_cds, dat_complete), make_fout(domain), sep = '\t', quote = FALSE, row.names = FALSE)
    }
}
