library(bedr)
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)

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

## check how many 6909 domains are true positives and what the false positives look like
b6_to_bed <- function(dirname, fout, pattern = "_predicted.tsv", cds_len_src = '', cds_len_f_pattern = "_merged.tsv"){
    alcolb6 <- as.data.frame(rbindlist(lapply(list.files(path=dirname, pattern=pattern, full.names=TRUE), fread, sep='\t')), stringsAsFactors=FALSE) %>%
        mutate(domain = sapply(qseqid, gen_split_extract('\\|', 3))) %>%
        group_by(qseqid) %>% top_n(bitscore, n = 2) %>% ungroup()
    
    for (i in 1:nrow(alcolb6)){
        start <- alcolb6[i,"sstart"]
        end <- alcolb6[i,"send"]
        alcolb6[i,"sstart"] <- min(start, end) - 1
        alcolb6[i,"send"] <- max(start, end)
    }
    if (cds_len_src != ''){
        pred_dat <- as.data.frame(rbindlist(lapply(list.files(cds_len_src, pattern = cds_len_f_pattern, full.names = TRUE), fread, sep = '\t')), stringsAsFactors = FALSE) %>%
            filter(accID == 6909) %>% select(c(contig, accID, hit.start, hit.end, strand, cds.len))
        alcolb6 <- alcolb6 %>% 
            mutate(contig = sapply(qseqid, gen_split_extract('\\|', 1)),
                   accID = sapply(qseqid, gen_split_extract('\\|', 2)),
                   ranges = sapply(qseqid, gen_split_extract('\\|', 5))) %>%
            separate(ranges, into = c("hit.start", "hit.end"), sep = '-', convert = TRUE) %>%
            left_join(pred_dat, by = c("contig", "accID", "hit.start", "hit.end"))
        write.table(alcolb6 %>% select(c(sseqid, sstart, send, qseqid, length, pident, bitscore, qlen, domain, cds.len)),
                    fout, quote = FALSE, row.names = FALSE, sep = '\t', col.names = FALSE)
    } else {
        alcolb6$cds.len <- c(NA)
        write.table(alcolb6 %>% select(c(sseqid, sstart, send, qseqid, length, pident, bitscore, qlen, domain, cds.len)),
                    fout, quote = FALSE, row.names = FALSE, sep = '\t', col.names = FALSE)
    }
}

## summarise discovered genes
summarise_pred_genes <- function(fname, fname_cds, fout, genes_f){
    genes <- read.table(genes_f, sep = '\t', header = FALSE, stringsAsFactors = FALSE)[,1]
    arath_pred_raw <- read.table(fname, header = FALSE, sep = '\t', stringsAsFactors = FALSE,
                                 col.names = c("chrom", "start", "end", "contig", "length", "pident", "bitscore",
                                               "qlen", "domain", "predcdslen", "bchrom", "bstart", "bend", "name", "score",
                                               "frame", "source", "feature", "phase", "attributes", "overlap")) %>%
        select(-c(bchrom, source, score))
    arath_pred_gpsg <- arath_pred_raw %>% filter(feature %in%  c("gene", "pseudogene")) %>%
        mutate(gene = sapply(attributes, function(x){split_extract(x, "Name=", 2)}),
               canon = gene %in% genes)
    ## arath_pred_cds <- read.table(fname_cds, header = FALSE, sep = '\t', stringsAsFactors = FALSE,
    ##                             col.names = c("chrom", "start", "end", "contig", "gene")) %>%
    ##    separate_rows(contig, sep = ',')
    tmp_arath_pred_cds <- arath_pred_raw %>% filter(feature %in% c("CDS")) %>%
        mutate(gene = sapply(attributes, function(x){split_extract(split_extract(x, "arent=", 2), "\\.\\d,", 1)})) %>%
        select(c(chrom, start, end, bstart, bend, contig, gene))
    tmp_arath_pred_cds <- tmp_arath_pred_cds %>%
        mutate(gcontig = sapply(1:nrow(tmp_arath_pred_cds),
                                function(i){paste(tmp_arath_pred_cds[i, "gene"],
                                                  tmp_arath_pred_cds[i, "contig"],
                                                  sep = ';')}))
    arath_pred_cds <- data.frame(gcontig = as.character(), cds_overlap = as.integer())
    for (gc in unique(tmp_arath_pred_cds$gcontig)){
        curr_arath_pred_cds <- tmp_arath_pred_cds %>% filter(gcontig == gc)
        curr_arath_pred_cds.cds <- curr_arath_pred_cds %>% select(c(chrom, bstart, bend, gcontig)) %>% 
            bedr.sort.region(check.chr = FALSE, verbose = FALSE) %>% bedr.merge.region(check.chr = FALSE, verbose = FALSE)
        curr_arath_pred_cds.intersect <- bedr(input = list(a = curr_arath_pred_cds %>% select(c(chrom, start, end)) %>% distinct() %>% 
                                                               bedr.sort.region(check.chr = FALSE, verbose = FALSE),
                                                           b = curr_arath_pred_cds.cds),
                                              method = "intersect", params = "-wo", check.chr = FALSE, verbose = FALSE) %>%
            mutate(cds_overlap_each = as.integer(V8))
        curr_arath_pred_cds.summary <- curr_arath_pred_cds.intersect %>% rename(gcontig = names) %>%
            group_by(gcontig) %>% summarise(cds_overlap = sum(cds_overlap_each))
        arath_pred_cds <- rbind(arath_pred_cds, curr_arath_pred_cds.summary)
    }
    arath_pred <- arath_pred_gpsg %>%
        mutate(gcontig = sapply(1:nrow(arath_pred_gpsg),
                                function(i){paste(arath_pred_gpsg[i, "gene"],
                                                  arath_pred_gpsg[i, "contig"],
                                                  sep = ';')})) %>%
        left_join(arath_pred_cds, by = c("gcontig")) %>% mutate_if(is.numeric, replace_na, replace = 0) %>%
        mutate(coding = cds_overlap > 0, realcdslen = cds_overlap)
    write.table(arath_pred %>% mutate(seqid = contig, nlr164 = canon, seqlen = qlen) %>%
                select(c(seqid, seqlen, length, pident, bitscore, feature, gene, domain, coding, nlr164, predcdslen, realcdslen)),
                fout, quote = FALSE, sep = '\t', row.names = FALSE)
}

