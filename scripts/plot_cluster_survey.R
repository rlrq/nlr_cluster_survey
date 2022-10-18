library(tidyverse)

library(grid)
library(gridExtra)
library(ggpubr)
library(wesanderson)
library(RColorBrewer)

library(Biostrings)

library(data.table)

dir_output <- "" ## update accordingly
dir_results <- "" ## update accordingly. Should be the same path as stored in 'results_dir' variable in pipeline_cluster_survey.sh

## files
fname_table_s11 <- "/path/to/table/s11.tsv" ## extract table S11 to tab-separated file and update path accordingly
fname_acc_legacy <- "/path/to/383Legacy.tsv" ## download 383Legacy.tsv from the github repo and update accordingly
fname_acc_1135 <- "/path/to/1135acc.csv" ## download 1135acc.csv from the github repo and update accordingly
fname_vdw_arch <- "/path/to/vdw_arch.tsv" ## Table S2f from Van de Weyer et al. (2019). Download vdw_arch.tsv from the github repo and update accordingly
fname_members <- "/path/to/members.tsv" ## Table S3a from Van de Weyer et al. (2019), filtered for transcripts starting with ALYR or accession IDs, and retaining columns Transcript_ID to Arch_Type. Download members.tsv from the github repo and update accordingly
fname_orthogroups_members <- "/path/to/orthogroups_members.csv" ## transcripts grouped by orthogroups as sorted by Van de Weyer et al. (2019). Download orthogroups_members.csv from the github repo and update accordingly
fname_ploidy2_alyrata <- "/path/to/ploidy2_alyrata.tsv" ## download ploidy2_alyrata.tsv from the github repo and update accordingly


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


## collapse clusters
collapse_clusters <- function(df, cluster_colname = "cluster", output_colname = "cluster",
                              rename_old_col = "cluster_old", collapse_name = "minor_cluster",
                              collapse_pattern = "^cAT"){
    df[, rename_old_col] <- df[[cluster_colname]]
    df[, output_colname] <- sapply(df[[cluster_colname]],
                                   function(x){if (str_detect(x, collapse_pattern)) {collapse_name}
                                               else {x}})
    return(df)
}

## extract legend
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

## make plot location
make_plot_loc <- function(s){
    return(paste0(dir_output, "/", s))
}

## get sequences that are identical across a contiguous stretch
get_identical_genes <- function(seqs_dat){
    genes_misassigned <- list()
    for (i in 1:(length(seqs_dat)-1)){
        s1 <- seqs_dat[[i]]
        g1 <- split_extract(names(seqs_dat)[[i]], '\\|', 2)
        for (j in (i+1):length(seqs_dat)){
            s2 <- seqs_dat[[j]]
            g2 <- split_extract(names(seqs_dat)[[j]], '\\|', 2)
            if (g1==g2 || Position(function(x) identical(x, sort(c(g1, g2))), genes_misassigned, nomatch=0) > 0){
                next
            } else if (s1 == s2) {
                genes_misassigned[[length(genes_misassigned) + 1]] <- sort(c(g1, g2))
            }
        }
    }
    return(genes_misassigned)
}


#############
##  TREES  ##
#############
## for making trees plots

## clusters
## 4B 4C 4D: B3, DM2, B5
## 5C 5D 5E 5F: DM8/RPP4/RPP5, DM4/RPP8, RPP13, RPS6

library(ape)
library(geiger)
library(ggtree)

## borrowed and adapted function from either the ape or geiger package
get_mrca <- function(t, tip){
    Ntip <- length(t$tip.label)
    rootnd <- Ntip + 1L
    pars <- integer(t$Nnode)
    tnd <- if (is.character(tip)) {match(tip, t$tip.label)} else {tip}
    done_v <- logical(Ntip + t$Nnode)
    pvec <- integer(Ntip + t$Nnode)
    pvec[t$edge[,2]] <- t$edge[,1]
    nd <-tnd[1]
    pars <- integer(t$Nnode)
    for (k in 1:t$Nnode){nd <- pvec[nd]; pars[k] <- nd; if(nd == rootnd) break}
    pars <- pars[1:k]
    mrcind <- integer(max(pars))
    mrcand <- pars[1]
    for (i in 2:length(tnd)){
        cnd <- tnd[i]; done <- done_v[cnd]
        if (is.na(done)){next}
        while(!done){done_v[cnd] <- TRUE; cpar <- pvec[cnd]; done <- done_v[cpar];
            if (cpar %in% pars){
                if (cpar == rootnd) {return(rootnd)}
                if (mrcind[cpar] > mrcind[mrcand]) mrcand <- cpar; done_v[cpar] <- TRUE; done <- TRUE}
            cnd <- cpar
        }
    }
    mrcand
}

make_tree_nolab <- function(d, c, exclude = c(), g = c(), group_name = '', dist = NA, only_arath = FALSE,
                            bg_col = "darkgrey", h1_col = "black", h2_col = h1_col, otu_gene = c(),
                            otu_domainid = c(), otu_seq = list(), otu_colours = c(), terminal = FALSE,
                            show_ref = TRUE, tip_colours = c(), tip_seq = list(),
                            colour_bootstrap = FALSE, bootstrap_low = "darkgreen", bootstrap_high = "red",
                            bootstrap_legend_pos = "bottom", ...){ ## h1/2 for highlights
    arbitrary_args_sink <- list(...)
    ## get tree
    t_fname <- paste0(dir_results, "/tree/nlr164_col0-AL70-Alyrata_",
                      d,"_mafft_ML.nwk")
    t <- read.tree(file = t_fname)
    ## get reference names
    all_ref_leaves <- c(names(readDNAStringSet(paste0(dir_results, "/domain_seq/nlr164_", d, "_ref_CDS_complete.fasta"))), names(readDNAStringSet(paste0(dir_results, "/domain_seq/nlr164_", d, "_6909pred_ref_CDS_complete.fasta"))))
    ## get data of predicted gene ids
    tmp_dat_pred <- dat_pred %>% dplyr::filter(domain == d)
    tmp_dat_pred_araly <- dat_pred_araly %>% filter(domain == d)
    ## get genes in desired cluster
    c_genes <- genes %>% dplyr::filter(str_detect(cluster, c)) %>% pull(gene)
    c_genes <- c_genes[! c_genes %in% exclude]
    ## get names of sequences assigned to genes in c_genes
    al_leaves <- tmp_dat_pred %>% dplyr::filter(gene %in% c_genes) %>% pull(contig)
    araly_leaves <- tmp_dat_pred_araly %>% dplyr::filter(gene %in% c_genes) %>% pull(contig)
    ref_leaves <- data.frame(seq_id = all_ref_leaves) %>%
        mutate(gene = sapply(seq_id, gen_split_extract('\\|', 2))) %>%
        dplyr::filter(gene %in% c_genes) %>% pull(seq_id) %>% as.character()
    ## if only A. thaliana sequences requested, filter
    if (only_arath){
        c_leaves <- c(al_leaves, ref_leaves)
    } else {
        c_leaves <- c(al_leaves, araly_leaves, ref_leaves)
    }
    ## prune tree to only keep relevant branches
    c_t <- keep.tip(t, c_leaves)
    ## if colouring bootstrap
    if (colour_bootstrap){
        p <- ggtree(c_t, aes(colour = as.numeric(label))) +
            scale_colour_continuous(low = bootstrap_low, high = bootstrap_high) +
            theme(legend.pos = bootstrap_legend_pos) +
            labs(colour = "bootstrap confidence")
        if (show_ref){
            p <- p +
                geom_tippoint(aes(subset=(node %in% which(c_t$tip %in% ref_leaves))),fill="cyan",shape=21)
        }
        return(p)
    ## if using groupOTU
    } else if (length(otu_gene) + length(otu_domainid) + length(otu_seq) > 0){
        if (length(otu_gene) > 0){ ## if otu names are gene names, get sequences automatically
            tmp_get_seqs <- function(g){
                if (!only_arath){
                    araly_leaves <- tmp_dat_pred_araly %>% dplyr::filter(gene == g) %>% pull(contig)
                } else {
                    araly_leaves <- c()
                }
                vdw_leaves <- tmp_dat_pred %>% dplyr::filter(gene == g) %>% pull(contig)
                ref_leaves <- data.frame(seq_id = all_ref_leaves) %>%
                    mutate(gene = sapply(seq_id, gen_split_extract('\\|', 2))) %>%
                    dplyr::filter(gene == g) %>% pull(seq_id) %>% as.character()
                return(c(vdw_leaves, araly_leaves, ref_leaves))
            }
            otu <- lapply(otu_gene, tmp_get_seqs)
        } else if (length(otu_domainid) > 0){ ## else if otu is for domain id
            tmp_get_seqs <- function(did){
                if (!only_arath){
                    araly_leaves <- tmp_dat_pred_araly %>% dplyr::filter(domainId == did) %>% pull(contig)
                } else {
                    araly_leaves <- c()
                }
                vdw_leaves <- tmp_dat_pred %>% dplyr::filter(domainId == did) %>% pull(contig)
                ref_leaves <- data.frame(seq_id = all_ref_leaves) %>%
                    mutate(domainId = sapply(seq_id, paste0(gen_split_extract('\\|', 2),
                                                            gen_split_extract('\\|', 6)))) %>%
                    dplyr::filter(domainId == did) %>% pull(seq_id) %>% as.character()
                return(c(vdw_leaves, araly_leaves, ref_leaves))
            }
            otu <- lapply(otu_domainid, tmp_get_seqs)
        } else if (length(otu_seq) > 0){ ## else if otu is of sequences names
            otu <- otu_seq
        }
        if (length(otu_colours > 0)){
            colours <- otu_colours
        } else {
            colours <- c(bg_col, palette(brewer.pal(n = max(length(otu), 3), name = "Dark2"))[1:length(otu)])
        }
        if (terminal){
            get_colour_num <- function(tip_lab, otu, cols){
                if (class(otu) != "list"){ return(cols[[1]]) }
                else {
                    for (i in 1:length(otu)){
                        if (tip_lab %in% otu[[i]]){ return (cols[[i]]) }
                    }
                }
            }
            tmp <- data.frame(node = 1:(Nnode(c_t)+length(c_t$tip.label)),
                              colour = c(sapply(c_t$tip.label, function(x){get_colour_num(x, otu, colours)}),
                                         rep("black", Nnode(c_t))))
            p <- ggtree(c_t) %<+% tmp + aes(colour=I(colour))
        } else {
            p <- ggtree(groupOTU(c_t, otu), aes(colour = group)) +
                geom_tippoint(aes(subset=(node %in% which(c_t$tip %in% araly_leaves))),fill="red",shape=21) +
                scale_colour_manual(values = colours)
        }
        if (length(tip_colours) > 0 && length(tip_seq) > 0){
            get_colour_num <- function(tip_lab, tip_seq, cols){
                for (i in 1:length(tip_seq)){
                    if (tip_lab %in% tip_seq[[i]]){ return (cols[[i]]) }
                }
            }
            tmp <- data.frame(node = 1:length(c_t$tip.label),
                              colour = c(sapply(c_t$tip.label, function(x){get_colour_num(x, tip_seq,
                                                                                          tip_colours)})))
            p <- ggtree(c_t) %<+% tmp + aes(colour=I(colour))
        }
        if (show_ref){
            p <- p +
                geom_tippoint(aes(subset=(node %in% which(c_t$tip %in% ref_leaves))),fill="cyan",shape=21)
        }
        return(p)
    } else if (length(g) == 0){
        ## return(ggtree(c_t, aes(colour = I(c(col, "black")))) +
        return(ggtree(c_t, colour = h1_col) +
               geom_tippoint(aes(subset=(node %in% which(c_t$tip %in% ref_leaves))),fill="cyan",shape=21) +
               geom_tippoint(aes(subset=(node %in% which(c_t$tip %in% araly_leaves))),fill="red",shape=21))
    } else { ## if specific gene is to be highlighted
        ## rename some variables to keep cluster sequences stored
        c_ref_leaves <- ref_leaves
        c_al_leaves <- al_leaves
        c_araly_leaves <- araly_leaves
        ## get sequences of gene (g) to be included
        al_leaves <- tmp_dat_pred %>% dplyr::filter(gene %in%  g) %>% pull(contig)
        araly_leaves <- tmp_dat_pred_araly %>% dplyr::filter(gene %in% g) %>% pull(contig)
        ref_leaves <- data.frame(seq_id = all_ref_leaves) %>%
            mutate(gene = sapply(seq_id, gen_split_extract('\\|', 2))) %>%
            dplyr::filter(gene %in% g) %>% pull(seq_id) %>% as.character()
        g_leaves <- c(al_leaves, araly_leaves, ref_leaves)
        ## get inter-leaf distance data
        fname <- paste0(dir_results, "/seq_distances/nlr164_NB-ARC_",
                        c, "_distances.tsv")
        dat_dists <- read.table(fname, sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
            left_join(dat_pred %>% dplyr::select(contig, gene), by = c("seqa" = "contig")) %>%
            dplyr::rename(genea = gene) %>%
            left_join(dat_pred %>% dplyr::select(contig, gene), by = c("seqb" = "contig")) %>%
            dplyr::rename(geneb = gene) %>% 
            filter(seqa %in% ref_leaves & geneb %in% g)
        ## identify distinct clades for the gene in (g) by splitting according to distance from reference
        g1_leaves <- dat_dists %>% dplyr::filter(distance <= dist & seqb %in% g_leaves) %>% pull(seqb)
        g2_leaves <- dat_dists %>% dplyr::filter(distance > dist & seqb %in% g_leaves) %>% pull(seqb)
        ## plot and return
        return(ggtree(groupClade(c_t, c(getMRCA(c_t, g1_leaves), getMRCA(c_t, g2_leaves))),aes(colour=group)) +
               theme(legend.position = 'none') + scale_colour_manual(values = c(bg_col, h1_col, h2_col)) +
               geom_tippoint(aes(subset=(node %in% which(c_t$tip %in% c_ref_leaves))),fill="cyan",shape=21) +
               geom_tippoint(aes(subset=(node %in% which(c_t$tip %in% c_araly_leaves))),fill="red",shape=21))
    }
}

## fig saves
letter_height <- 279
letter_col <- 85
ggsave_1 <- function(p, fname, h, ...){ggsave(make_plot_loc(fname), p, height = h, width = 85,
                                              units = "mm", ...)}
ggsave_1_5 <- function(p, fname, h, ...){ggsave(make_plot_loc(fname), p, height = h, width = 114,
                                                units = "mm", ...)}
ggsave_2 <- function(p, fname, h, ...){ggsave(make_plot_loc(fname), p, height = h, width = 174,
                                              units = "mm", ...)}
ggsave_1s <- function(p, fname, h, s, ...){ggsave(make_plot_loc(fname), p, height = h, width = 85,
                                                  units = "mm", scale = s, ...)}
ggsave_2s <- function(p, fname, h, s, ...){ggsave(make_plot_loc(fname), p, height = h, width = 174,
                                                  units = "mm", scale = s, ...)}


################
##  METADATA  ##
##########################
##  ACCESSIONS + GENES  ##
##########################

## accession information
acc_relict <- read.table("/path/to/accID_relict.tsv", sep = '\t', header = FALSE)[,1] ## update accordingly with path to accID_relict.tsv clone from the repo
acc_map <- read.table(fname_table_s11,
                      sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(accID = as.character(accID), accID_old = as.character(accID_old),
           relict = sapply(accID, function(x){if(x %in% acc_relict){'y'}else{'n'}}),
           latitude = str_split(latitude, '/') %>% lapply(as.numeric) %>% lapply(mean) %>% unlist,
           longitude = str_split(longitude, '/') %>% lapply(as.numeric) %>% lapply(mean) %>% unlist,
           altitude = str_split(altitude, '/') %>% lapply(as.numeric) %>% lapply(mean) %>% unlist)
## coordinates
acc1135 <- read.table(fname_acc_1135, sep=',', header=T) %>%
    mutate(accID = as.character(tg_ecotypeid)) %>% dplyr::select(c(accID, latitude, longitude))
accLegacy <- read.table(fname_acc_legacy, sep='\t', header=T) %>%
    mutate(accID = as.character(Accession.ID), latitude = Lat, longitude = NA) %>%
    dplyr::select(accID, longitude, latitude)
accCoord <- rbind(acc1135, accLegacy %>% filter(! accID %in% acc1135$accID)) ## only acc1135 has longitude data

## genes
dat_vdw_arch <- read.table(fname_vdw_arch,
                           sep = '\t', stringsAsFactors = FALSE, header = TRUE) %>%
    tidyr::separate_rows(Transcript, sep = ';') %>%
    dplyr::filter(str_detect(Transcript, "AT\\dG\\d")) %>%
    dplyr::select(-c(pair_putpair_transcripts)) %>%
    dplyr::mutate(gene = str_extract(Transcript, "AT\\dG\\d{5}"),
                  TIR = str_detect(Subclass, 'T'),
                  CC = str_detect(Subclass, 'C'),
                  RPW8 = str_detect(Subclass, 'R')) ## integrate architecture info with genes
genes <- read.table(paste0(dir_results, "/nlr164_final_wPred_newClust.tsv"),
                    sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
    left_join(dat_vdw_arch %>% dplyr::select(gene, TIR, CC, RPW8), by = "gene")

## review genes supposedly without TIR/CC/R and manually add domains where necessary
for (i in 1:nrow(genes)){
    g <- genes[i, "gene"]
    ## no Nterm (explicitly defined to show we didn't forget these genes)
    if (g %in% c("AT1G61300", "AT3G41160", "AT4G09360", "AT5G38350", "AT5G45440", "AT5G47280", "AT1G10920",
                 "AT1G52660", ## this row and the one above: no coiled coil or rx_cc-like
                 "AT5G47250", "AT4G27220", "AT4G26090", "AT5G45490", "AT3G15700")) { ## this row: has coiled coil but not rx_cc-like
        next
    }
    ## TIR
    if (g %in% c("AT3G14470", "AT4G12020", "AT5G45260", "AT5G45050", "AT5G17890")) {
        genes[i, "TIR"] <- TRUE
    }
    ## CC
    if (g %in% c("AT3G46530", "AT3G46710", "AT3G46730", "AT1G58848", "AT1G59620", "AT3G07040", "AT1G59218")) { ## these genes: have rx_cc-like
        genes[i, "CC"] <- TRUE
    }
    ## R
    if (g %in% c()) {
        genes[i, "RPW8"] <- TRUE
    }
    ## no CC
    if (g %in% c("AT5G17890", "AT5G66900")){
        genes[i, "CC"] <- FALSE
    }
}

genes <- genes %>%
    left_join(genes %>% gather(Nterm, presence, TIR:RPW8) %>% dplyr::filter(presence) %>%
              dplyr::group_by(cluster) %>%
              dplyr::summarise(Nterms = paste0(Nterm %>% unique %>% sort, collapse = ';')),
              by = "cluster") %>%
    replace_na(list(Nterms = "None"))
pseudogenes <- (genes %>% filter(feature == "pseudogene"))$gene
clust_inc <- (genes %>% filter(! str_detect(cluster, "singleton")))$cluster %>% unique()

## domain counts per gene + genes w/ identical domains
library(rlist)
domains <- c("TIR", "NB-ARC")
domain_count <- list()
genes_misassigned <- list()
for (d in domains){
    ## tmp_seqs_f <- paste0(dir_results, "/nlr164_col0_", d, ".fasta")
    tmp_seqs_f <- paste0(dir_results, "/nlr164/domain_seq/nlr164_", d, "_ref_CDS_complete.fasta")
    tmp_seqs_dat <- readDNAStringSet(tmp_seqs_f)
    ## identify genes with anomalous NB-ARC count (!= 1)
    tmp_domain_count <- data.frame(seqid = names(tmp_seqs_dat)) %>%
        mutate(gene = sapply(seqid, gen_split_extract('\\|', 2)),
               domainNum = as.numeric(sapply(seqid, gen_split_extract('\\|', 6)))) %>% dplyr::select(-c(seqid))
    tmp_domain_count <- rbind(tmp_domain_count,
                              data.frame(gene = (genes %>%
                                                 filter(original.bait == 'n' & str_detect(domains, d)))$gene,
                                         domainNum = c(1)))
    domain_count[[d]] <- tmp_domain_count %>% group_by(gene) %>% summarise(domainNum = max(domainNum))
    genes_misassigned[[d]] <- get_identical_genes(tmp_seqs_dat)
}
## manually add the new genes/pseudogenes that are identical to other genes
genes_misassigned[["NB-ARC"]] <- list.append(genes_misassigned[["NB-ARC"]], c("AT1G72850", "AT1G72852"), c("AT5G43725", "AT5G43730"), c("AT5G48770", "AT5G48775"), c("AT4G26090", "AT4G26095"))


############
##  DATA  ##
##############################
##  predicted AL70 domains  ##
##############################
##  TABLES S4, S...genes  ##
############################

## get predicted gene identity data
dat_pred_raw <- as.data.frame(rbindlist(lapply(list.files(path=paste0(dir_results, "/predicted_identity"), pattern="nlr164_AL70", full.names=TRUE), fread, sep='\t')), stringsAsFactors=FALSE) %>%
    mutate(domainId = domain,
           accID = sapply(contig, function(x){unlist(strsplit(x, split="\\|"))[[2]]}),
           domain = sapply(contig, function(x){unlist(strsplit(x, split="\\|"))[[3]]})) %>%
    filter(accID %in% acc_map$accID) %>%
    mutate(relict = sapply(accID, function(x){if(x %in% acc_relict){"y"} else{"n"}}))
## drop sisters, expand columns about closest references, then get unique rows
dat_pred_all <- dat_pred_raw %>% separate_rows(gene, ref_name, cluster, sep = ';') %>% dplyr::select(-c(cluster)) %>% left_join(genes %>% dplyr::select(gene, cluster), by = "gene")
dat_pred_all <- dat_pred_all %>% left_join(dat_pred_all %>% group_by(contig) %>% summarise(count = n()), by = "contig")

## completely identical domains from different genes were doubly assigned to closest homologues
##   reduce dat_pred_all into dat_pred by only retaining one of each pair of doubly assigned genes
## this list and the following loop arbitrarily picks one of the identical genes to keep
dat_pred <- dat_pred_all
for (d in unique(dat_pred_all$domain)){
    throw_genes <- dat_pred_all %>% filter(count > 1 & domain == d) %>% group_split(contig) %>% lapply(function(l) {sort(l$gene)}) %>% unique() %>% lapply(function(l) {l[2:length(l)]}) %>% unlist()
    dat_pred <- dat_pred %>% filter(domain != d | (domain == d & ! gene %in% throw_genes))
}

## clean up data
dat_pred_genes <- rbind(dat_pred %>% dplyr::select(accID, domain, gene) %>% filter(domain == "NB-ARC") %>%
                        group_by(accID, gene) %>% summarise(count = n()) %>% ungroup() %>%
                        tidyr::complete(accID, gene) %>% mutate(domain = "NB-ARC"),
                        dat_pred %>% dplyr::select(accID, domain, gene) %>% filter(domain == "TIR") %>%
                        group_by(accID, gene) %>% summarise(count = n()) %>% ungroup() %>%
                        tidyr::complete(accID, gene) %>% mutate(domain = "TIR")) %>%
    replace_na(list(count = 0)) %>% left_join(acc_map %>% dplyr::select(accID, relict), by = c("accID")) %>%
    left_join(genes %>% dplyr::select(cluster, cluster_type, gene, original.bait, Nterms), by = c("gene")) %>%
    distinct() %>% dplyr::select(accID, relict, domain, cluster_type, cluster, gene, original.bait, count)

## ## TABLES S8 S9 (PREDICTED_NB-ARC_GENE) (PREDICTED_TIR_GENE)
## write.table(dat_pred_genes, paste0(dir_results, "/tables/nlr164_genes_accs.tsv"),
##             sep = '\t', row.names = FALSE, quote = FALSE)

## SUMMARISE predicted AL70 genes
## find accessions with minimum and maximum number of each domain (total, not limited to clusters)
dat_pred_total <- dat_pred_genes %>% group_by(accID, domain) %>% summarise(total=sum(count))
dat_pred_min_max <- dat_pred_total %>% group_by(domain) %>%
    summarise(min=min(total), max=max(total)) %>% column_to_rownames("domain")
dat_pred_total[,"total_sum"] <- sapply(1:nrow(dat_pred_total),
                                  function(i){domain <- unlist(dat_pred_total[i, "domain"]);
                                      count <- unlist(dat_pred_total[i, "total"]);
                                      if (count == dat_pred_min_max[domain, "max"]){"max"}
                                      else if (count == dat_pred_min_max[domain, "min"]){"min"}
                                      else {"other"}})

dat_pred_sum <- full_join(full_join(dat_pred_total,
                                    dat_pred %>% group_by(accID, domain, cluster) %>% summarise(count = n()),
                                    by = c("accID", "domain")),
                          do.call(expand.grid, list(domain = unique(dat_pred$domain),
                                                    cluster = unique(dat_pred$cluster),
                                                    accID = unique(dat_pred$accID))),
                          by=c("accID", "domain", "cluster")) %>%
    arrange(domain, total) %>% ungroup() %>% replace_na(list(count = 0))

dat_pred_sum <- dat_pred_sum %>%
    full_join(genes %>% dplyr::select(cluster, cluster_type, Nterms) %>% distinct(), by = "cluster") %>%
    left_join(acc_map %>% dplyr::select(accID, relict), by = "accID") %>% distinct()
dat_pred_sum[is.na(dat_pred_sum)] <- 0
dat_pred_sum[,"total_sum"] <- sapply(dat_pred_sum$total_sum, function(x){if(x == 0){"other"} else{x}})

dat_pred_sum <- dat_pred_sum %>% dplyr::select(-c(total)) %>%
    left_join(dat_pred_sum %>% group_by(accID, domain) %>% summarise(total = sum(count)),
              by = c("accID", "domain")) %>%
    left_join(dat_pred_sum %>% filter(! str_detect(cluster, "singleton")) %>%
              group_by(accID, domain) %>% summarise(total_cluster = sum(count)),
              by = c("accID", "domain")) %>%
    left_join(dat_pred_sum %>% filter(!str_detect(cluster, "^(cAT|singleton)")) %>%
              group_by(accID, domain) %>% summarise(total_major_cluster = sum(count)),
              by = c("accID", "domain"))

## ## TABLE S3 (PREDICTED_REPERTOIRE_SIZE)
## write.table(dat_pred_sum %>%
##             select(accID,relict,domain,cluster,cluster_type,count,total,total_cluster,total_major_cluster),
##             paste0(dir_results, "/tables/predicted_AL70_summary.tsv"),
##             sep = '\t', quote = FALSE, row.names = FALSE)


############
##  DATA  ##
###################################
##  predicted A. lyrata domains  ##
###################################

## ylg-only dataset
dat_pred_araly_raw <- as.data.frame(rbindlist(lapply(list.files(path=paste0(dir_results, "/predicted_identity"), pattern="nlr164_Alyrata", full.names=TRUE), fread, sep='\t')), stringsAsFactors=FALSE) %>%
    mutate(domainId = domain,
           domain = sapply(ref_name, function(x){split_extract(x, '\\|', 5)}),
           source = sapply(contig, function(x){if(str_detect(x, "^_R_")){x <- split_extract(x, "R_", 2)}
               return(split_extract(x, '\\|', 1))}))
dat_pred_araly_raw <- dat_pred_araly_raw %>%
    mutate(num = sapply(1:nrow(dat_pred_araly_raw),
                        function(i){split_extract(dat_pred_araly_raw[i, "contig"], '\\|', 3)}))

## drop sisters, expand columns about closest references, then get unique rows
dat_pred_araly_all <- dat_pred_araly_raw %>%  separate_rows(gene, ref_name, cluster, sep = ';') %>% dplyr::select(-c(sisters)) %>% dplyr::select(-c(cluster)) %>% left_join(genes %>% dplyr::select(gene, cluster), by = "gene")
dat_pred_araly_all <- dat_pred_araly_all %>% full_join(dat_pred_araly_all %>% group_by(contig) %>% summarise(count = n())) %>% ungroup()

## completely identical domains from different genes were doubly assigned to closest homologues
##   reduce dat_pred_araly_all into dat_pred_araly by only retaining one of each pair of doubly assigned genes
## this list and the following loop arbitrarily picks one of the identical genes to keep
dat_pred_araly <- dat_pred_araly_all
for (d in unique(dat_pred_araly_all$domain)){
    throw_genes <- dat_pred_araly_all %>% filter(count > 1 & domain == d) %>% group_split(contig) %>% lapply(function(l) {sort(l$gene)}) %>% unique() %>% lapply(function(l) {l[2:length(l)]}) %>% unlist()
    dat_pred_araly <- dat_pred_araly %>% filter(domain != d | (domain == d & ! gene %in% throw_genes))
}

dat_pred_araly_genes <- rbind(dat_pred_araly %>% dplyr::select(domain, gene) %>% filter(domain == "NB-ARC") %>%
                              group_by(gene) %>% summarise(count = n()) %>% ungroup() %>%
                              mutate(domain = "NB-ARC"),
                              dat_pred_araly %>% dplyr::select(domain, gene) %>% filter(domain == "TIR") %>%
                              group_by(gene) %>% summarise(count = n()) %>% ungroup() %>%
                              mutate(domain = "TIR")) %>%
    replace_na(list(count = 0)) %>%
    full_join(genes %>% dplyr::select(cluster, cluster_type, gene, original.bait), by = c("gene")) %>%
    distinct() %>% dplyr::select(domain, cluster_type, cluster, gene, original.bait, count)


## ## TABLE S6 (PREDICTED_ARALY_GENE)
## write.table(dat_pred_araly %>% select(-c(cluster)) %>%
##             left_join(genes %>% select(gene, cluster, cluster_type), by = "gene") %>%
##             select(source, domain, num, cluster_type, cluster, gene) %>%
##             dplyr::rename(araly_protein_ID = source, domain_order = num, Col0_gene = gene),
##             paste0(dir_results, "/tables/nlr164_genes_Araly.tsv"),
##             sep = '\t', row.names = FALSE, quote = FALSE)


############
##  DATA  ##
##############################
##  AL70 domains diversity  ##
##############################

## pi: nucleotide diversity
pi_f_nbs <- paste0(dir_results, "/seq_diversity/nlr164_AL70_NB-ARC_gapsExc_CDSonly_pi.tsv")
pi_f_tir <- paste0(dir_results, "/seq_diversity/nlr164_AL70_TIR_gapsExc_CDSonly_pi.tsv"
dat_pi <- rbind(read.table(pi_f_nbs, header = TRUE, stringsAsFactors = FALSE, sep = '\t') %>%
                dplyr::select(-c(cluster)) %>%
                left_join(genes %>% dplyr::select(gene, cluster), by = "gene") %>%
                collapse_clusters() %>% mutate(domain = "NB-ARC"),
                read.table(pi_f_tir, header = TRUE, stringsAsFactors = FALSE, sep = '\t') %>%
                dplyr::select(-c(cluster)) %>%
                left_join(genes %>% dplyr::select(gene, cluster), by = "gene") %>%
                collapse_clusters() %>% mutate(domain = "TIR"))

############
##  DATA  ##
##############################
##  INTRA-CLUSTER DISTANCE  ##
##############################

f_src_tmp <- function(s){
    if (str_detect(s, 'Col-0_ref')){'ref'} else if (str_detect(s, '.tig')){'vdw'} else {'araly'}}
tmp_arathly_pred <- rbind(dat_pred %>% dplyr::select(contig, gene, domainId, cluster),
                          dat_pred_araly %>% dplyr::select(contig, gene, domainId, cluster))
dat_dists <- as.data.frame(rbindlist(lapply(list.files(path=paste0(dir_results, "/seq_distances"), pattern="nlr164_NB-ARC.+distances.tsv", full.names=TRUE), fread, sep = '\t')),
                           stringsAsFactors = FALSE, fill = TRUE) %>%
    left_join(tmp_arathly_pred, by = c("seqa" = "contig")) %>%
    dplyr::rename(genea = gene, domainIda = domainId, clustera = cluster) %>%
    left_join(tmp_arathly_pred, by = c("seqb" = "contig")) %>%
    dplyr::rename(geneb = gene, domainIdb = domainId, clusterb = cluster) %>%
    dplyr::mutate(srca = sapply(seqa, f_src_tmp), srcb = sapply(seqb, f_src_tmp))


############
##  DATA  ##
################################
##  AL70 domains popgenstats  ##
################################

dat_pops <- rbind(read.table(paste0(dir_results, "/popg_test/nlr164_AL70_NB-ARC_CDSonly_stats.tsv"),
                             sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
                  mutate(region = "CDS"),
                  read.table(paste0(dir_results, "/popg_test/nlr164_AL70_NB-ARC_intrononly_stats.tsv"),
                             sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
                  mutate(region = "noncoding"),
                  read.table(paste0(dir_results, "/popg_test/nlr164_AL70_NB-ARC_CDScomplete_stats.tsv"),
                             sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
                  mutate(region = "complete")) %>%
    left_join(genes, by = c("gene"))
dat_pops_gathered <- dat_pops %>%
    gather(test, value, nucleotide_diversity:wattersons_theta)

dat_pops_clade <- read.table(paste0(dir_results, "/popg_test/nlr164_NB-ARC_CDSonly_clade_stats.tsv"),
                             sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(region = "CDS", cluster = group)
dat_pops_clade_gathered <- dat_pops_clade %>%
    gather(test, value, nucleotide_diversity:wattersons_theta) %>%
    mutate(value = as.numeric(value)) %>%
    drop_na()
dat_pops_clade <- dat_pops_clade_gathered %>%
    spread(test, value)


############
##  DATA  ##
###########################
##  Van de Weyer et al.  ##
###########################

dat_vdw <- read.table(fname_members,
                      header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(accID_old = sapply(Gene_ID, gen_split_extract('\\|', 1)),
           gene = sapply(Gene_ID, gen_split_extract('\\|', 2))) %>%
    filter(accID_old != 7063) %>%
    left_join(read.table(fname_table_s11,
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
              mutate(accID_old = as.character(accID_old)), by = c("accID_old")) %>%
    mutate(accID = as.character(accID))
dat_vdw_nb <- dat_vdw[sapply(dat_vdw$Architecture, function(x){"NB" %in% unlist(str_split(x, ','))}),]

dat_vdw_orth_raw <- read.table(fname_orthogroups_members,
                               header = FALSE, stringsAsFactors = FALSE)[,1]
dat_vdw_orth <- data.frame(seqid = as.character(), accID_old = as.character(), ogNum = as.character())
for (i in 1:length(dat_vdw_orth_raw)){
    tmp_og <- dat_vdw_orth_raw[[i]]
    tmp_genes <- unlist(str_split(tmp_og, ','))
    tmp <- data.frame(seqid = tmp_genes, accID_old = sapply(tmp_genes, gen_split_extract('\\|', 1)),
                      seq = sapply(tmp_genes, gen_split_extract('\\|', 2)), ogNum = i)
    dat_vdw_orth <- rbind(dat_vdw_orth, tmp)
}
dat_vdw_orth <- dat_vdw_orth %>% mutate(seq = sapply(seqid, gen_split_extract('\\|', 2))) %>%
    left_join(read.table(fname_table_s11,
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
              mutate(accID_old = as.character(accID_old)), by = c("accID_old")) %>%
    mutate(accID = as.character(accID), gene = sapply(seq, gen_split_extract('\\.', 1)))


######################
##  FIGURES. 1, S1  ##
##     VALIDITY,    ##
##     TOTAL_NLR    ##
######################

acc_relict_raw <- acc_relict
acc_relict <- data.frame(accID = as.character(acc_relict_raw)) %>%
    left_join(acc_map %>% dplyr::select(accID, accessions)) %>%
    pull(accessions)

## VDW QC plot data
dat_vdw_nb_sum_all <- dat_vdw_nb[sapply(dat_vdw_nb$Gene_ID, function(x){!str_detect(x, "6909\\|AT")}),] %>%
    group_by(accID) %>% summarise(count = n()) %>% filter(accID %in% acc_map$accID) %>%
    mutate(rank = rank(dplyr::desc(count),ties.method = "first"))
## summarise AL70 predicted data
dat_nb <- dat_pred_sum %>% filter(domain == "NB-ARC") %>% select(c(accID, total)) %>% distinct() %>%
    mutate(count = total, rank = rank(dplyr::desc(count), ties.method = "first")) %>% select(-total)
dat_vdw_nb_sum_canonOG <- dat_vdw_orth %>% filter(! str_detect(seqid, "6909\\AT")) %>% drop_na() %>%
    filter(ogNum %in% filter(dat_vdw_orth,gene %in% genes$gene)$ogNum & seqid %in% dat_vdw_nb$Transcript_ID
           & accID %in% acc_map$accID) %>%
    group_by(accID) %>% summarise(count = n()) %>% mutate(rank = rank(dplyr::desc(count),ties.method = "first"))
## combine both data to be plotted together    
dat_combined <- rbind(cbind(dat_nb, source = c("BLAST pipeline")),
                      cbind(dat_vdw_nb_sum_all, source = c("VdW")),
                      cbind(dat_vdw_nb_sum_canonOG, source=c("VdW nlr164OG")))%>%
    mutate(accID = as.character(accID)) %>%
    left_join(acc_map %>% dplyr::select(accID, accessions))
tmp_src <- dat_combined$source
    

## Figs. 1B, S1
## composition bar plot
bar_plots_all <- list()
bar_plots_cluster <- list()
bar_plots_major <- list()
for (d in unique(dat_pred_sum$domain)){
    tmp_dat <- dat_pred_sum %>% filter(domain %in% d) %>%
        left_join(acc_map %>% dplyr::filter(accID_old != 7063) %>% dplyr::select(accID, accessions)) %>%
        mutate(accID = accessions)
    tmp_dat_totals <- tmp_dat %>% group_by(accID) %>% summarise(total = sum(count))
    accID_order <- tmp_dat_totals %>% arrange(-total) %>% pull(accID)
    tmp_dat$accID <- factor(tmp_dat$accID, levels = accID_order)
    tmp_median <- tmp_dat_totals %>% pull(total) %>% median
    tmp_dat <- tmp_dat %>% collapse_clusters()
    bar_plots_all[[d]] <- ggplot(tmp_dat, aes(x = accID, y = count, fill = cluster)) +
        geom_hline(yintercept = tmp_median, lty = 2) +
        geom_bar(stat = "identity", colour="black", size=0.5) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle=90, colour=ifelse(accID_order %in% acc_relict, "red", "black"),
                                         size = 8, hjust = 1, vjust = 0.5),
              legend.text = element_text(size = 8)) +
        ## xlab("accID") +
        xlab("accession") +
        scale_fill_discrete(guide = guide_legend(ncol = 1)) +
        ylab(paste0("Number of ", d, " domains"))
    tmp_dat <- tmp_dat %>% filter(! str_detect(cluster_type, "singleton"))
    tmp_dat_totals <- tmp_dat %>% group_by(accID) %>% summarise(total = sum(count))
    accID_order <- tmp_dat_totals %>% arrange(-total) %>% pull(accID)
    tmp_dat$accID <- factor(tmp_dat$accID, levels = accID_order)
    tmp_median <- tmp_dat_totals %>% pull(total) %>% median
    bar_plots_cluster[[d]] <- ggplot(tmp_dat, aes(x = accID, y = count, fill = cluster)) +
        geom_hline(yintercept = tmp_median, lty = 2) +
        geom_bar(stat = "identity", colour="black", size=0.5) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle=90, colour=ifelse(accID_order %in% acc_relict, "red", "black"),
                                         size = 8, hjust = 1, vjust = 0.5),
              legend.text = element_text(size = 8)) +
        ## xlab("accID") +
        xlab("accession") +
        scale_fill_discrete(guide = guide_legend(ncol = 1)) +
        ylab(paste0("Number of ", d, " domains"))
    tmp_dat <- tmp_dat %>% filter(! str_detect(cluster_type, "minor"))
    tmp_dat_totals <- tmp_dat %>% group_by(accID) %>% summarise(total = sum(count))
    accID_order <- tmp_dat_totals %>% arrange(-total) %>% pull(accID)
    tmp_dat$accID <- factor(tmp_dat$accID, levels = accID_order)
    tmp_median <- tmp_dat_totals %>% pull(total) %>% median
    bar_plots_major[[d]] <- ggplot(tmp_dat, aes(x = accID, y = count, fill = cluster)) + 
        geom_hline(yintercept = tmp_median, lty = 2) +
        geom_bar(stat = "identity", colour="black", size=0.5) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle=90, colour=ifelse(accID_order %in% acc_relict, "red", "black"),
                                         size = 8, hjust = 1, vjust = 0.5),
              legend.text = element_text(size = 8)) +
        ## xlab("accID") +
        xlab("accession") +
        scale_fill_discrete(guide = guide_legend(ncol = 1)) +
        ylab(paste0("Number of ", d, " domains"))
}

## Fig. 1A
accID_order_al70 <- (dat_vdw_nb_sum_all %>% filter(accID %in% dat_nb$accID) %>% arrange(rank) %>%
                     left_join(acc_map %>% dplyr::filter(accID_old != 7063) %>%
                               dplyr::select(accID, accessions)) %>%
                     mutate(accID = accessions))$accID
p_qc <- ggplot(dat_combined %>% left_join(dat_vdw_nb_sum_all, by = c("accID")) %>%
               left_join(acc_map %>% dplyr::filter(accID_old != 7063) %>%
                         dplyr::select(accID, accessions)) %>%
               dplyr::mutate(accID = accessions) %>%
               dplyr::filter(!str_detect(source, "OG")) %>%
               mutate(count = count.x, rank = rank.y) %>% drop_na()) +
    geom_line(aes(x = reorder(accID, rank), y = count, group = source, colour = source),
              size = 1.3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,
                                     colour = ifelse(accID_order_al70 %in% acc_relict, "red", "black")),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 12), legend.pos = "bottom",
          legend.title = element_blank()) +
    ylab("Number of NB-ARC domains or NB-ARC-containing genes") + ## xlab("accID") +
    xlab("accession") +
    ## scale_colour_manual(values = wes_palette("Cavalcanti1"))##  +
    scale_colour_brewer(palette = "Dark2")

## Fig. 1
p1 <- arrangeGrob(arrangeGrob(p_qc + guides(colour = guide_legend(nrow = 1, byrow = TRUE)),
                             left = textGrob(expression(bold('A')), vjust = 2,
                                             x = unit(1,"npc"), y = unit(1, "npc"))),
                 arrangeGrob(bar_plots_major[["NB-ARC"]] +
                             guides(fill = guide_legend(nrow=3,byrow=FALSE)) +
                             theme(legend.pos = "bottom", legend.title = element_blank(),
                                   legend.spacing.y = unit(0.05, "cm"), legend.spacing.x = unit(0.1, "cm")),
                             left = textGrob(expression(bold('B')), vjust = 2,
                                             x = unit(1,"npc"), y = unit(1, "npc"))),
                 heights = c(1, 1), padding = unit(1, "npc"))
ggsave_2(p1, "Fig1.pdf", letter_height)
ggsave_2(p1, "Fig1.png", letter_height, dpi = 300)

## Fig. S1
ps1 <- arrangeGrob(arrangeGrob(arrangeGrob(bar_plots_all[["NB-ARC"]] + theme(legend.pos = "None"),
                                             left = textGrob(expression(bold('A')), vjust = 2,
                                                             x = unit(1, "npc"), y = unit(1, "npc"))),
                                 arrangeGrob(bar_plots_all[["TIR"]] + theme(legend.pos = "None"),
                                             left = textGrob(expression(bold('B')), vjust = 2,
                                                             x = unit(1, "npc"), y = unit(1, "npc"))),
                                 arrangeGrob(bar_plots_major[["TIR"]] + theme(legend.pos = "None"),
                                             left = textGrob(expression(bold('C')), vjust = 2,
                                                             x = unit(1, "npc"), y = unit(1, "npc"))),
                                 heights = c(1,1,1)),
                   arrangeGrob(g_legend(bar_plots_all[["NB-ARC"]] +
                                        guides(fill = guide_legend(nrow = 4)))),
                   heights = c(5,1))
ggsave_2(ps1, "FigS1.pdf", letter_height + 20)
ggsave_2(ps1, "FigS1.png", letter_height + 20, dpi = 300)

acc_relict <- acc_relict_raw

#####################
##  Figures 2, S2  ##
##    GEOGRAPHY    ##
#####################

## note that only the 57 of Anna-Lena's accessions in acc1135 have longitude data
d <- "NB-ARC"

## Figs. 2, S2 data
tmp_cluster_ranked <- dat_pred_sum %>% filter(domain == d & cluster_type == "major") %>%
    dplyr::select(c(accID, cluster, count, Nterms)) %>% group_by(cluster, Nterms) %>%
    mutate(rank = rank(dplyr::desc(count), ties.method = "min"))
clust_order <- (tmp_cluster_ranked %>% group_by(cluster, Nterms) %>% summarise(var = var(count)) %>%
                ungroup() %>% arrange(var))$cluster
long_dat_ll_ordered <- tmp_cluster_ranked %>% left_join(accCoord, by = "accID") %>% ungroup() %>%
    mutate(cluster = factor(cluster, levels = clust_order))

## Figs. 2, S2 data (all)
tmp_cluster_ranked_all <- dat_pred_sum %>% filter(domain == d) %>%
    dplyr::select(c(cluster_type, accID, cluster, count, Nterms)) %>% group_by(cluster_type, cluster, Nterms) %>%
    mutate(rank = rank(dplyr::desc(count), ties.method = "min"))
clust_order_all <- (tmp_cluster_ranked_all %>% group_by(cluster_type, cluster, Nterms) %>%
                    summarise(var = var(count)) %>% ungroup() %>% arrange(var))$cluster
long_dat_ll_ordered_all <- tmp_cluster_ranked_all %>% left_join(accCoord, by = "accID") %>% ungroup() %>%
    mutate(cluster = factor(cluster, levels = clust_order_all)) %>%
    dplyr::filter(!str_detect(cluster, "singleton")) %>%
    dplyr::mutate(Nterms = factor(Nterms, levels = c("TIR", "CC", "RPW8", "None")))
clust_order <- clust_order_all
long_dat_ll_ordered <- long_dat_ll_ordered_all

## Fig. 2
p2_hist_all <- ggplot(long_dat_ll_ordered_all, aes(x = count)) +
    geom_histogram(aes(fill = Nterms), alpha = 0.7, binwidth = 1) +
    geom_vline(data = long_dat_ll_ordered_all %>% filter(accID == 6909),
               aes(xintercept = count), colour = "red", lty = 2) +
    facet_wrap(~cluster, scales = "free") + theme_bw() +
    theme(legend.pos = "bottom", panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          strip.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 6)) +
    scale_x_continuous(breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]),
                                                by = (if (x[2] - x[1] <= 6){1}
                                                      else{tmp <- ceiling(min(floor(sqrt(x[2]) * 1.2), x[2]));
                                                    if(tmp %% 2 != 0){max(1, tmp - 1)}else{max(1, tmp)}}))) +
    xlab("number of NB-ARC domains in cluster") + ylab("number of accessions") +
    scale_fill_brewer(palette = "Dark2")
p2 <- p2_hist_all
ggsave_2s(p2_hist_all, "Fig2.pdf", 150, 1)
ggsave_2s(p2_hist_all, "Fig2.png", 150, 1, dpi = 300)

## Fig. S2
## coordinates, by rank
p2s<- ggplot(long_dat_ll_ordered %>% drop_na(), aes(x=longitude, y=latitude, size = -rank+1, fill = -rank+1)) +
    geom_point(pch = 21, colour = "black", alpha = 0.3, stroke = 0.5) +
    facet_wrap("cluster") + labs(size = "rank") +
    guides(fill = guide_legend(reverse = TRUE), size = guide_legend(reverse = TRUE)) +
    scale_size_continuous(name = "rank", labels = c(60, 40, 20, 1)) +
    scale_fill_gradient(name = "rank", high = "red", low = "black", labels = c(60, 40, 20, 1))
ggsave_1s(p2s, "FigS2.pdf", 70, 3)
ggsave_1s(p2s, "FigS2.png", 70, 3, dpi = 300)



#################
##  Table S5   ##
#################
## Table S5 (MORANI)

## spatial correlation
library(ape)
morani.geography <- list()
for (d in unique(dat_pred_sum$domain)){
    d_dat <- data.frame(cluster = as.character(), observed = as.numeric(),
                        expected = as.numeric(), sd = as.numeric(), p.value = as.numeric())
    long_dat_ll <- dat_pred_sum %>% filter(domain == d) %>% select(c(accID, cluster, count)) %>%
        group_by(cluster) %>% left_join(accCoord) %>% filter(!str_detect(cluster, "singleton")) %>% drop_na()
    for (c in unique(long_dat_ll$cluster)){
        dat <- long_dat_ll %>% filter(cluster == c)
        if (nrow(dat) == 0 || length(unique(dat$count)) == 1) {next}
        dat.dists <- as.matrix(dist(cbind(dat$longitude, dat$latitude)))
        dat.dists.inv <- 1/dat.dists
        dat.dists.inv[is.infinite(dat.dists.inv)] <- 0
        mi <- Moran.I(dat$count, dat.dists.inv)
        d_dat <- rbind(d_dat, data.frame(cluster = c, observed = mi$observed,
                                         expected = mi$expected, sd = mi$sd, p.value = mi$p.value))        
    }
    morani.geography[[d]] <- d_dat
    write.table(d_dat, paste0(dir_results, "/tables/nlr164_AL70_MoranI_geography_",
                              d, ".txt"), quote = FALSE, sep = '\t', row.names = FALSE)
}

geno <- read.table(paste0(dir_results, "/hierarchy/al70_acc_new_pairwiseSNPs_cAcc_upperTri.tsv"),
                   header = FALSE, fill = TRUE, sep = '\t', row.names = 1)
genos <- rownames(geno)
geno <- geno %>% arrange(-row_number())
m <- Matrix::forceSymmetric(as.matrix(geno), uplo="L")
geno <- as.data.frame(as.matrix(m), row.names = rev(genos))
colnames(geno) <- rev(genos)

morani.hierarchy <- list()
for (d in unique(dat_pred_sum$domain)){
    d_dat <- data.frame(cluster = as.character(), observed = as.numeric(),
                        expected = as.numeric(), sd = as.numeric(), p.value = as.numeric())
    long_dat_ll <- dat_pred_sum %>% filter(domain == d) %>% select(c(accID, cluster, count)) %>%
        group_by(cluster) %>% left_join(accCoord) %>% filter(!str_detect(cluster, "singleton")) %>%
        mutate(accID = sapply(accID, function(x){paste0('a', x)})) %>% arrange(accID) %>%
        filter(accID %in% genos)
    not_in_dat <- genos[! genos %in% long_dat_ll$accID]
    tmp_geno <- geno %>% filter(! rownames(geno) %in% not_in_dat) %>% select(-not_in_dat)
    for (c in unique(long_dat_ll$cluster)){
        dat <- long_dat_ll %>% drop_na() %>% filter(cluster == c)
        if (nrow(dat) == 0 || length(unique(dat$count)) == 1) {next}
        dat.dists <- as.matrix(tmp_geno)
        dat.dists.inv <- 1/dat.dists
        dat.dists.inv[row(dat.dists.inv) == col(dat.dists.inv)] <- 0
        dat.dists.inv[is.infinite(dat.dists.inv)] <- 0
        mi <- Moran.I(dat$count, dat.dists.inv)
        d_dat <- rbind(d_dat, data.frame(cluster = c, observed = mi$observed,
                                         expected = mi$expected, sd = mi$sd, p.value = mi$p.value))        
    }
    morani.hierarchy[[d]] <- d_dat
    write.table(d_dat, paste0(dir_results, "/tables/nlr164_AL70_MoranI_hierarchy_",
                              d, ".txt"), quote = FALSE, sep = '\t', row.names = FALSE)
}


################
##  Figure 3  ##
#####################
##  A. lyrata CNV  ##
#####################

acc_novikova <- read.table(fname_ploidy2_alyrata,
                           sep = '\t', stringsAsFactors = FALSE, header = TRUE) %>%
    mutate(subspecies = str_extract(sample, "(?<=Alyrata)[^\\d]+"))

## Note that araly domainIds that were not found to be deleted or duplicated in other accessions (and thus did not have their nRD stats output by CNVnator) were assigned nRD of 1
dat_araly_cnv <- read.table(paste0(dir_results, "/araly_cnv/aralyCNV_called_novikova14_NB-ARC_summary.tsv"), header = TRUE, sep = '\t', stringsAsFactors = FALSE) %>%
    dplyr::mutate(proteinId = as.character(proteinId), domainId = as.character(domainId)) %>%
    dplyr::mutate_all(~replace(., . == '', NA)) %>%
    tidyr::gather("cnv_type", "to_split", deletion:duplication) %>%
    rbind(data.frame(domainId = c(NA), proteinId = c(NA), cnv_type = c("unchanged"), to_split = c(NA))) %>%
    tidyr::separate_rows(to_split, sep = ';') %>%
    tidyr::separate(to_split, into = c("acc", "nRD"), sep = ',') %>%
    dplyr::full_join(dat_pred_araly %>% dplyr::select(c(contig, gene, cluster, source, domainId)) %>%
                     dplyr::mutate(col0domainId = domainId,
                                   domainId = stringr::str_extract(contig, "^\\d+\\|NB-ARC\\|\\d+")) %>%
                     dplyr::select(-c(contig)),
                     by = c("domainId", "proteinId" = "source")) %>%
    tidyr::complete(nesting(proteinId, domainId, cluster, gene, col0domainId), acc,
                    fill = list(cnv_type = "unchanged", nRD = 1)) %>%
    dplyr::left_join(genes %>% dplyr::select(c(cluster, cluster_type)), by = "cluster") %>%
    tidyr::drop_na() %>%
    dplyr::distinct() %>%
    dplyr::mutate(nRD = as.numeric(nRD))

## Table S7 (NOVIKOVA_CNVnator)
## ## write nRD to table (mean of nRD taken for domains with multiple nRD)
## ## ## domainId-acc combinations with multiple nRD: 883372|NB-ARC|1-SRR2040795 (deletion (0.145)/deletion (0.347)); 916988|NB-ARC|1-SRR2040792 (deletion (0.330)/duplication (1.87))
## write.table(dat_araly_cnv %>% dplyr::select(-c(cnv_type)) %>%
##             group_by(proteinId, domainId, cluster_type, cluster, gene, col0domainId, acc) %>%
##             summarise(meannRD = mean(nRD)) %>% ungroup() %>% spread(acc, meannRD),
##             paste0(dir_results, "/tables/nlr164_domain_Araly_CNVnator_meannRD.tsv"),
##             sep = '\t', quote = FALSE, row.names = FALSE)

## combine with Fig3
d <- "NB-ARC"
tmp_dat <- dat_pred_sum %>% filter(cluster_type != "singleton" & domain == d)
tmp_total <- tmp_dat %>% group_by(cluster, cluster_type) %>%
    summarise(total = sum(count)) %>% filter(total != 0) %>% ungroup()
tmp_dat <- tmp_dat %>% dplyr::filter(cluster %in% tmp_total$cluster) %>%
    dplyr::select(accID, cluster, cluster_type, count) %>% mutate(species = "A. thaliana") %>%
    rbind(dat_araly_cnv %>% group_by(acc, cluster_type, cluster) %>% dplyr::summarise(count = sum(nRD)) %>%
          ungroup() %>% filter(!str_detect(cluster, "singleton")) %>% dplyr::rename(accID = acc) %>%
          mutate(species = "A. lyrata")) %>%
    tidyr::complete(nesting(species, accID), nesting(cluster_type, cluster), fill = list(count = 0)) %>%
    dplyr::left_join(tmp_dat %>% dplyr::group_by(cluster) %>% dplyr::summarise(athaliana_count = sum(count)), by = "cluster") %>%
    dplyr::mutate(species = factor(species, levels = c("A. thaliana", "A. lyrata")))

## Fig. 3
library(ggsignif)
library(ggpubr)
scale <- 20
p3extcol <- ggplot(tmp_dat %>% drop_na(),
                   aes(x = reorder(cluster, -athaliana_count, group = factor(species)), y = count)) +
    geom_boxplot(fill = NA, outlier.alpha=0, width=0.7, alpha = 0, colour = NA) +
    geom_bar(data = tmp_total, aes(x = cluster, y = total / scale),
             stat = "identity", position = "dodge", fill = "light grey", colour = "black", width = 1) +
    geom_boxplot(aes(x = cluster, y = count, colour = species),
                 outlier.alpha=0.5, width=0.7) +
    facet_grid(.~cluster_type, space = "free_x", scales = "free_x") +
    ylab("number of NB-ARC per genome") + xlab("cluster") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=60, hjust=1), strip.text.x = element_text(size = 10),
          legend.pos = "bottom", panel.grid.major.x = element_blank()) +
    scale_y_continuous(sec.axis = sec_axis(~.*scale, name = "total NB-ARC in 64 accessions")) +
    coord_cartesian(ylim = c(-2.5, 30)) +
    scale_colour_brewer(palette = "Dark2") +
    stat_compare_means(aes(group = species), method = "wilcox.test", paired = FALSE, label = "p.signif",
                       label.y = -2, angle = 90)
p3 <- p3extcol
ggsave_2(p3extcol, "Fig3.pdf", 150)
ggsave_2(p3extcol, "Fig3.png", 150, dpi = 300)

## Fig. S3
## groups nRD by the col0 domainId assigned to their araly domainId
p_araly_cnv_col0domain <- ggplot(dat_araly_cnv %>%
       dplyr::filter(! str_detect(cluster, "singleton")) %>%
       dplyr::group_by(cluster_type, cluster, gene, col0domainId, acc) %>%
       dplyr::summarise(sumnRD = sum(nRD))) +
    geom_hline(aes(yintercept = 1), lty = 2) +
    geom_point(aes(x = acc, y = sumnRD, colour = acc), alpha = 0.5, size = 2) +
    theme_bw() +
    theme(legend.pos = "none",
          axis.text.x = element_text(angle = 90, hjust = 1,
                                     colour = ifelse((dat_araly_cnv %>% pull(acc) %>% unique %>% sort) %in%
                                                   (acc_novikova %>% filter(subspecies == "lyrata") %>%
                                                       pull(SRR)), "red", "black"))) +
    facet_wrap(~cluster) +
    geom_text(data = dat_araly_cnv %>% filter(! str_detect(cluster, "singleton")) %>%
                  dplyr::select(c(cluster_type, cluster, gene, col0domainId)) %>% dplyr::distinct() %>%
                  dplyr::group_by(cluster_type, cluster) %>% dplyr::summarise(count = n()),
              aes(label = paste0("Ndomains = ", count)), x = 1, y = 12, hjust = 0) +
    scale_y_continuous(limits = c(-1, 13)) +
    xlab("sample") + ylab("sum of normalised read depth")
p3s <- p_araly_cnv_col0domain
ggsave_2s(p3s, "FigS3.pdf", 150, 1.5)
ggsave_2s(p3s, "FigS3.png", 150, 1.5, dpi = 300)


#############################
##    GENE_COUNT_NB-ARC    ##
##       (bar, ecdf)       ##
##  GENE_COUNT_BOX_NB-ARC  ##
##    (bar, box, point)    ##
#############################
library(ggsignif)

make_summary_plots <- function(domain, reduced_dat, gene_count_dat, genes_misassigned = list(), scale = 64,
                               gene_acc_count_dat = data.frame(), gene_count_araly = data.frame(),
                               coloured_clusters = c(), show_1n = TRUE, show_2n = TRUE){
    gene_order_asc <- (reduced_dat %>% arrange(count))$gene
    gene_order_pos <- (reduced_dat %>% arrange(gene))$gene
    p_bar_asc <- ggplot(reduced_dat, aes(x = gene, y = count)) +
        geom_bar(data = reduced_dat %>% filter(cluster %in% coloured_clusters),
                 aes(x=reorder(gene, count), y = count, fill = cluster), stat="identity") +
        geom_bar(data = reduced_dat %>% filter(!cluster %in% coloured_clusters),
                 aes(x=reorder(gene, count), y = count), stat = 'identity', fill = 'grey') +
        geom_hline(yintercept = 64, lty = 2) +
        geom_hline(yintercept = 128, lty = 3) +
        geom_text(aes(label = per_gene),
                  position = position_dodge(width = 0.9), vjust = -0.25, colour = "red") +
        theme(axis.text.x = element_text(angle=90, size=7,
                                         colour=ifelse(gene_order_asc %in% pseudogenes, "red", "black")),
              legend.position = "bottom", axis.title = element_blank()) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        scale_x_discrete(limits = gene_order_asc)
    p_bar_pos <- ggplot(reduced_dat, aes(x = gene, y = count)) +
        geom_rect(data = reduced_dat[1,], fill = "red", alpha = 0.1, ymin = -Inf, ymax = Inf,
                  xmin = 0.5, xmax = nrow(filter(reduced_dat, chrom <= 1)) + 0.5) +
        geom_rect(data = reduced_dat[1,], fill = "orange", alpha = 0.1, ymin = -Inf, ymax = Inf,
                  xmin = nrow(filter(reduced_dat,chrom<=1))+0.5,xmax = nrow(filter(reduced_dat,chrom<=2))+0.5) +
        geom_rect(data = reduced_dat[1,], fill = "green", alpha = 0.1, ymin = -Inf, ymax = Inf,
                  xmin = nrow(filter(reduced_dat,chrom<=2))+0.5,xmax = nrow(filter(reduced_dat,chrom<=3))+0.5) +
        geom_rect(data = reduced_dat[1,], fill = "cyan", alpha = 0.1, ymin = -Inf, ymax = Inf,
                  xmin = nrow(filter(reduced_dat,chrom<=3))+0.5,xmax = nrow(filter(reduced_dat,chrom<=4))+0.5) +
        geom_rect(data = reduced_dat[1,], fill = "purple", alpha = 0.1, ymin = -Inf, ymax = Inf,
                  xmin = nrow(filter(reduced_dat,chrom<=4))+0.5,xmax = nrow(filter(reduced_dat,chrom<=5))+0.5) +
        geom_bar(data = reduced_dat %>% filter(cluster %in% coloured_clusters),
                 aes(x=gene, y = count, fill = cluster), stat="identity") +
        geom_bar(data = reduced_dat %>% filter(! cluster %in% coloured_clusters),
                 aes(x=gene, y = count), stat = 'identity', fill = "darkgrey") +
        geom_hline(yintercept=64, lty=2) + 
        geom_hline(yintercept=128, lty=3) + 
        geom_text(aes(label = per_gene),
                  position = position_dodge(width = 0.9), vjust = -0.25, colour = "red") +
        theme(axis.text.x = element_text(angle=90, size=7,
                                         colour=ifelse(gene_order_pos %in% pseudogenes, "red", "black")),
              legend.position = "bottom",
              axis.title = element_blank()) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        scale_x_discrete(limits = gene_order_pos)
    if (length(genes_misassigned) > 0){
        p_bar_asc <- p_bar_asc +
            geom_signif(comparisons = genes_misassigned, annotation = "", tip_length = -0.005,
                        y_position = 5 * seq(-1, -length(genes_misassigned), -1))
        p_bar_pos <- p_bar_pos +
        geom_signif(comparisons = genes_misassigned, annotation = "", tip_length = 0.005,
                    y_position = sapply(genes_misassigned,
                                        function(x) {max(filter(reduced_dat, gene %in% x)$count) + 10}))
    }
    p_bar_pos_box <- p_bar_pos +
        geom_boxplot(data = gene_count_dat %>% filter(cluster %in% coloured_clusters),
                     aes(x = gene, y = count * scale), outlier.alpha=0.4, fill = NA, width=0.5) +
        geom_boxplot(data = gene_count_dat %>% filter(! cluster %in% coloured_clusters),
                     aes(x = gene, y = count * scale), outlier.alpha=0.4, fill = NA, width=0.5) +
        scale_y_continuous(sec.axis = sec_axis(~./scale, name = paste0(domain, " domains per genome"),
                                               breaks = c(0, 1, 2, 5, 10)))
    if (nrow(gene_acc_count_dat) > 0){
        p_bar_pos_box <- p_bar_pos_box +
            geom_point(data = gene_acc_count_dat %>% filter(count > 0),
                       aes(x = gene, y = count), colour = "cyan", pch = 1)
    }
    if (nrow(gene_count_araly) > 0){
        p_bar_pos_box <- p_bar_pos_box +
            geom_point(data = gene_count_araly %>% full_join(data.frame(gene = reduced_dat$gene)) %>%
                           replace_na(list(count = 0)),
                       aes(x = gene, y = count * scale), colour = "red", pch = 1, size = 2)
    }
    p_bar_pos_l <- g_legend(p_bar_pos)
    p_ecdf <- ggplot() +
        stat_ecdf(data = reduced_dat, aes(count)) +
        geom_vline(xintercept = 64, lty = 2) + 
        ylab("cumulative genes frequency") + xlab(paste0("num assigned ", domain, " domains")) +
        ggtitle(paste0("ECDF of num. assigned ", domain, " domains per gene"))
    if (show_1n){
        p_bar_asc <- p_bar_asc +
            geom_text(aes(0, 64, label = "1N = 64", vjust=-1, hjust=-0.3))
        p_bar_pos <- p_bar_pos +
            geom_text(aes(0, 64, label = "1N = 64", vjust=-1, hjust=-0.3))
        p_ecdf <- p_ecdf +
            geom_text(aes(64, 0, label = "1N = 64", hjust = -0.3))
    }
    if (show_2n){
        p_bar_asc <- p_bar_asc +
            geom_text(aes(0, 128, label = "2N = 128", vjust=-1,hjust=-0.3))
        p_bar_pos <- p_bar_pos +
            geom_text(aes(0, 128, label = "2N = 128", vjust=-1,hjust=-0.3))
    }
    p_bar <- arrangeGrob(arrangeGrob(p_bar_pos + theme(legend.position="none"),
                                     left=textGrob('A', just = c("left", "top"),
                                                   x = unit(1,"npc"), y = unit(1, "npc"))),
                         arrangeGrob(p_bar_asc + theme(legend.position = "none"),
                                     left=textGrob('B', just = c("left", "top"),
                                                   x = unit(1,"npc"), y = unit(1, "npc"))),
                         heights = c(1, 1), bottom=textGrob("gene", gp=gpar(fontsize = 10)),
                         left=textGrob(paste0("num. assigned ",domain," domains"),gp=gpar(fontsize=10), rot=90))
    p_bar_l <- arrangeGrob(p_bar, p_bar_pos_l, heights = c(10, 1),
                           top = paste0("Number of ", domain, " domains assigned to each canonical NLR"))
    p <- arrangeGrob(p_bar_l, arrangeGrob(p_ecdf, left=textGrob('C', just = c("left", "top"),
                                                                x = unit(1,"npc"), y = unit(1, "npc"))),
                     top =textGrob(paste0("Anna-Lena 64 Gene Summary; ",domain," domain"),gp=gpar(fontsize=15)),
                     widths = c(3, 1))
    return(list(main = p, bar_l = p_bar_l, bar = p_bar, bar_asc = p_bar_asc, bar_pos = p_bar_pos,
                bar_pos_box = p_bar_pos_box))
}

d <- "TIR"
make_p4 <- function(d, show_per_gene = TRUE, ...){
    tmp_domain_count <- domain_count[[d]] %>% full_join(genes %>% dplyr::select(gene), by = "gene") %>%
        replace_na(list(domainNum = 0))
    tmp_genes_misassigned <- genes_misassigned[[d]]
    major_clusters <- (genes %>% dplyr::filter(gene %in% tmp_domain_count$gene &
                                               !str_detect(cluster, "^(cAT|singleton)")))$cluster %>% unique()
    minor_clusters <- (genes %>% dplyr::filter(gene %in% tmp_domain_count$gene &
                                               str_detect(cluster, "^cAT")))$cluster %>% unique()
    ## coerce data into nice, simple summaries
    tmp_dat_pred <- dat_pred %>% dplyr::filter(domain == d) %>% dplyr::select(-c(cluster)) %>%
        left_join(genes %>% dplyr::select(gene, cluster), by = "gene")
    tmp_dat_pred_sum <- dat_pred_genes %>% dplyr::filter(domain == d)
    tmp_reduced_dat <- tmp_dat_pred_sum %>% group_by(gene, cluster, cluster_type) %>%
        summarise(count = sum(count)) %>%
        full_join(tmp_domain_count, by = c("gene")) %>% mutate_all(~replace(., is.na(.), 0)) %>%
        mutate(chrom = as.numeric(sapply(gene, function(x){unlist(str_extract(x,"(?<=AT)\\d(?=G)"))[[1]]})))
    if (show_per_gene){
        tmp_reduced_dat <- tmp_reduced_dat %>%
            mutate(per_gene = sapply(domainNum, function(x) {if (x == 1) {''} else {as.character(x)}}))
    } else {
        tmp_reduced_dat <- tmp_reduced_dat %>%
            mutate(per_gene = '')
    }
    tmp_gene_acc_count_dat <- tmp_dat_pred %>% dplyr::select(-c(contig, distance)) %>%
        distinct() %>% group_by(gene) %>% summarise(count = n()) %>% full_join(genes) %>%
        mutate_all(~replace(., is.na(.), 0))
    tmp_gene_count_araly <- dat_pred_araly %>%
        dplyr::filter(domain == "NB-ARC" & str_detect(contig, "complete")) %>%
        ## filter(src_type == "ylg" & domain == "NB-ARC" & str_detect(contig, "complete")) %>%
        group_by(gene) %>% summarise(count = n())
    ## make plots: major clusters
    plots <- make_summary_plots(d, tmp_reduced_dat, tmp_dat_pred_sum, genes_misassigned = tmp_genes_misassigned,
                                gene_acc_count_dat = tmp_gene_acc_count_dat,
                                gene_count_araly = tmp_gene_count_araly,
                                coloured_clusters = major_clusters, ...)
    return(plots)
}

## save various plots
## Fig. 4A
## note that geom_segment used to be extended 0.3 units towards each other
p4a <- make_p4("NB-ARC")$bar_pos + theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
                                         legend.text = element_text(size = 8), legend.spacing.x = unit(1, "mm"),
                                         axis.title = element_text()) +
    guides(fill = guide_legend(nrow = 3, byrow = FALSE)) + ylab("number of NB-ARC domains") +
    geom_segment(aes(x = 51 - 0.5, xend = 55 + 0, y = -2, yend = -2)) +
    geom_segment(aes(x = 56 - 0, xend = 62 + 0.5, y = -2, yend = -2)) +
    geom_text(aes(x = 53, y = -4, label = "P+", vjust=1)) +
    geom_text(aes(x = 59, y = -4, label = "P-", vjust=1))

## Fig. 4B 4C 4D
d <- "NB-ARC"
b3 <- make_tree_nolab(d, "B3"); b3
d2 <- make_tree_nolab(d, "DM2"); d2
b5 <- make_tree_nolab(d, "B5", exclude = c("AT1G72920")); b5

## Fig. 4
p4a2 <- arrangeGrob(p4a + theme(legend.pos = "None"), g_legend(p4a), heights = c(4, 1),
                   left = textGrob(expression(bold('A')), vjust = 2, x = unit(1,"npc"), y = unit(1, "npc")))
p4bcd <- arrangeGrob(arrangeGrob(b3, left = textGrob(expression(bold('B')), vjust = 2,
                                                     x = unit(1,"npc"), y = unit(1, "npc"))),
                     arrangeGrob(d2, left = textGrob(expression(bold('C')), vjust = 2,
                                                     x = unit(1,"npc"), y = unit(1, "npc"))),
                     arrangeGrob(b5, left = textGrob(expression(bold('D')), vjust = 2,
                                                     x = unit(1,"npc"), y = unit(1, "npc"))),
                     widths = c(1, 1, 1))
p4 <- arrangeGrob(p4a2, p4bcd, heights = c(2, 1))
ggsave_2(p4, "Fig4.pdf", 150)
ggsave_2(p4, "Fig4.png", 150, dpi = 300)


## Fig. S5A
p5sn <- make_p4("NB-ARC", show_per_gene = FALSE, show_1n = FALSE, show_2n = FALSE)$bar_pos_box +
    theme(legend.text = element_text(size = 8), legend.spacing.x = unit(1, "mm"), axis.title = element_text()) +
    guides(fill = guide_legend(nrow = 3, byrow = FALSE)) + ylab("number of NB-ARC domains") +
    geom_segment(aes(x = 51 - 0.5, xend = 55 + 0, y = -20, yend = -20)) +
    geom_segment(aes(x = 56 - 0, xend = 62 + 0.5, y = -20, yend = -20)) +
    geom_text(aes(x = 53, y = -30, label = "P+", vjust=1)) +
    geom_text(aes(x = 59, y = -30, label = "P-", vjust=1))

## Fig. S5B
p5st <- make_p4("TIR", show_per_gene = FALSE, show_1n = FALSE, show_2n = FALSE)$bar_pos_box +
    theme(legend.text = element_text(size = 8), legend.spacing.x = unit(1, "mm"), axis.title = element_text()) +
    guides(fill = guide_legend(nrow = 3, byrow = FALSE)) + ylab("number of TIR domains") +
    geom_segment(aes(x = 51 - 0.5, xend = 55 + 0, y = -20, yend = -20)) +
    geom_segment(aes(x = 56 - 0, xend = 62 + 0.5, y = -20, yend = -20)) +
    geom_text(aes(x = 53, y = -30, label = "P+", vjust=1)) +
    geom_text(aes(x = 59, y = -30, label = "P-", vjust=1))

## Fig. S4
p5s <- arrangeGrob(arrangeGrob(p5sn + theme(legend.pos = "None"),
                               left = textGrob(expression(bold('A')), vjust = 2,
                                               x = unit(1,"npc"), y = unit(1, "npc"))),
                   arrangeGrob(p5st + theme(legend.pos = "None"),
                               left = textGrob(expression(bold('B')), vjust = 2,
                                               x = unit(1,"npc"), y = unit(1, "npc"))),
                   g_legend(p5sn + guides(fill = guide_legend(nrow = 2, byrow = FALSE))),
                   heights = c(4, 4, 1))
ggsave(make_plot_loc("FigS5.pdf"), p5s, height = 215, width = 297, units = "mm")
ggsave(make_plot_loc("FigS5.png"), p5s, height = 215, width = 297, units = "mm", dpi = 300)


####################
##  FIGURES 5 S6  ##
#########################
##    PI_BOX_NB-ARC    ##
##  (bar, box, point)  ##
#########################
library(ggpmisc)

plot_ordered_pi_bar <- function(dat, cluster_order, cluster_count, title, scale,
                                show_max = FALSE, d = "", ggproto = NULL){
    p <- ggplot() + ggproto +
        geom_bar(data = cluster_count, aes(x = reorder(cluster, -m), y = count/scale),
                 stat = "identity", fill = "grey") +
        geom_point(data = dat, aes(x = cluster, y = pi.pairalnlen, size = count),
                   ## colour = ifelse(dat$domainNum == 1, "black", "red"), alpha = 0.3) +
                   colour = "#085da4", alpha = 0.3) +
        geom_boxplot(data = dat %>% drop_na(),
                     aes(x = cluster, y = pi.pairalnlen), outlier.alpha = 0, width = 0.5, fill=NA) +
        ylab("pi") + xlab("cluster") +
        scale_y_continuous(sec.axis = sec_axis(~.*scale, name = paste("total number of",d,"domains",sep = ' ')),
                           expand = expand_scale(mult = c(0.01, 0.2))) +
        theme(axis.text.x = element_text(angle=60, hjust=1))
    if (show_max){
        p <- p +
            geom_text(data = dat %>% group_by(cluster, cluster_type) %>%
                          top_n(1, pi.pairalnlen) %>% ungroup() %>%
                          group_by(cluster, cluster_type, pi.pairalnlen) %>%
                          summarise(gene = paste(gene, collapse='/')),
                      aes(x = cluster, y = pi.pairalnlen, label = gene),
                      angle = 90, nudge_y = 0.01, size = 3, hjust = 0)
    }
    return(p)
}


## For Fig. 5A S6
scale <- 10000
domains <- c("NB-ARC", "TIR")
d <- "NB-ARC"

## pi, calculated for each gene (domains combined if same gene)
## some minor clusters expanded
tmp_dat_pred <- dat_pred %>% dplyr::filter(domain == d) %>% dplyr::select(-c(cluster)) %>% drop_na() %>%
    left_join(genes %>% dplyr::select(gene, cluster), by = "gene")
tmp_dat_pi <- dat_pi %>% filter(domain == d) %>% select(-c(cluster)) %>% drop_na() %>%
    left_join(genes %>% select(gene, cluster), by = "gene")
tmp_cluster_count <- tmp_dat_pred %>% filter(domain == d) %>%
    group_by(cluster) %>% summarise(count = sum(count)) %>% ungroup() %>%
    left_join(genes %>% select(cluster, cluster_type) %>% unique, by = "cluster")
tmp_gene_count <- tmp_dat_pred %>% group_by(gene) %>% summarise(count = sum(count))
tmp_dat_pi <- tmp_dat_pi %>% left_join(domain_count[[d]], by = c("gene")) %>%
    left_join(tmp_gene_count) %>% mutate(domainNum = as.factor(domainNum)) %>% select(-c(cluster)) %>%
    left_join(genes %>% select(gene, cluster, cluster_type), by = "gene")
tmp_dat_pi <- tmp_dat_pi %>% left_join(tmp_dat_pi %>% group_by(cluster) %>%
                                       summarise(m = median(pi.pairalnlen)), by = "cluster") %>%
    mutate(cluster_type = sapply(cluster_type, function(x){if(str_detect(x, "singleton")){""}else{x}}))
tmp_cluster_count <- tmp_cluster_count %>% left_join(tmp_dat_pi %>% group_by(cluster) %>%
                                                     summarise(m = median(pi.pairalnlen)), by = "cluster") %>%
    mutate(cluster_type = sapply(cluster_type, function(x){if(str_detect(x, "singleton")){""}else{x}}))

## pi, calculated by biopython, for each domain (no combining domains in same gene)
tmp_dat_pi_py <- dat_pops %>% dplyr::select(gene, cluster, cluster_type, nucleotide_diversity, domainId, region,
                                            numSeqRaw, numSeqTrim, tajimas_d, wattersons_theta) %>%
    dplyr::mutate(domain = "NB-ARC", domainNum = 1,
                  pi.pairalnlen = nucleotide_diversity, count = numSeqTrim) %>%
    dplyr::filter(region == "CDS")
tmp_dat_pi_py <- tmp_dat_pi_py %>% left_join(tmp_dat_pi_py %>% drop_na() %>% group_by(cluster) %>%
                                             summarise(m = median(pi.pairalnlen)), by = "cluster") %>%
    mutate(cluster_type = sapply(cluster_type, function(x){if(str_detect(x, "singleton")){""}else{x}}),
           gene = domainId)
tmp_cluster_count_py <- tmp_dat_pred %>% dplyr::filter(domain == d) %>%
    group_by(cluster) %>% summarise(count = sum(count)) %>% ungroup() %>%
    left_join(genes %>% dplyr::select(cluster, cluster_type) %>% unique, by = "cluster")
tmp_cluster_count_py <- tmp_cluster_count_py %>%
    left_join(tmp_dat_pi_py %>% drop_na() %>% group_by(cluster) %>%
              summarise(m = median(pi.pairalnlen)),by = "cluster") %>%
    mutate(cluster_type = sapply(cluster_type, function(x){if(str_detect(x, "singleton")){""}else{x}}))

## For Fig. 5B
## some other data manip
dat_split_genes <- dat_split %>%
    dplyr::filter(type %in% c("radiation", "hifi", "hifi_multi")) %>%
    ## dplyr::filter(type != "hifi_multi") %>% 
    dplyr::select(type, group, genes) %>%
    tidyr::separate_rows(genes, sep = ',') %>%
    dplyr::distinct() %>%
    dplyr::rename(domainId = genes)
dat_pops_gathered_clade <- dat_pops_gathered %>%
    left_join(dat_split_genes, by = c("domainId")) %>%
    replace_na(list(type = "unclassified")) %>%
    mutate(type = sapply(type, function(x){if(x == "hifi_multi"){"hifi"}else{x}}))
## reclassify DM2 and DM8 domains that are involved in both hifi and radiation as 'unclassified'
tmp_messy_domains <- c("AT3G44630.1", "AT4G16930.1")
dat_pops_gathered_clade <- dat_pops_gathered_clade %>%
    dplyr::filter(! (domainId %in% tmp_messy_domains & type %in% c("radiation")))
dat_pops_gathered_clade <- dat_pops_gathered_clade %>%
    dplyr::mutate(type = sapply(1:nrow(dat_pops_gathered_clade),
                                function(i){if(dat_pops_gathered_clade[i,"domainId"] %in% tmp_messy_domains)
                                            {return("unclassified")}else
                                                                   {dat_pops_gathered_clade[i,"type"]}})) %>%
    mutate(type = factor(type, levels = c("unclassified", "radiation", "hifi_multi", "hifi")))
tmp_conflated_domains <- c("AT5G35450.1", "AT4G16900.1", "AT3G46530.1", "AT1G72840.1", "AT1G72870.1")

## Fig. 5B
library(ggridges)
p_ntw <- ggplot(dat_pops_gathered_clade %>% dplyr::filter(region == "CDS")) +
    geom_vline(data = dat_pops_clade_gathered %>%
                   dplyr::filter(type != "hifi_multi") %>%
                   dplyr::filter(type != "cluster"),
               aes(xintercept = value, colour = type), size = 1.5, alpha = 0.5) +
    geom_vline(data = data.frame(test = c("tajimas_d"), vline = c(-2, 2)),
               aes(xintercept = vline), colour = "red", lty = 2) +
    geom_histogram(aes(x = value, fill = type), colour = NA, alpha = 1) +
    scale_fill_manual(values = wes_palette("Royal1")) +
    scale_colour_brewer(palette = "Dark2") +
    geom_histogram(data = dat_pops_gathered %>%
                       dplyr::filter(region != "complete" & region == "CDS" &
                                     domainId %in% tmp_conflated_domains),
                   aes(x = value), colour = "black", fill = NA) + ## show conflated genes
    facet_wrap(~test, ncol = 1, scales = "free") +
    labs(colour = "clade", fill = "domain")
p_5b <- p_ntw

## Fig. 5A
p5a_domainsep <- 
    plot_ordered_pi_bar(dat = tmp_dat_pi_py, cluster_order = c(), scale = scale,
                        cluster_count = tmp_cluster_count_py, show_max = TRUE, d = "NB-ARC",
                        title = paste0("Nucleotide diversity; Gaps excluded; CDS only; ",
                                       d, " domain"),
                        ggproto = geom_hline(data = tmp_dat_pi %>% dplyr::filter(domainNum == 2) %>%
                                                 complete(cluster_type, pi.pairalnlen),
                                             aes(yintercept = pi.pairalnlen),
                                             lty = 2, colour = "red", size = 0.6, alpha = 0.7)) +
    facet_grid(~cluster_type, scales = "free_x", space = "free_x") + theme(legend.position = "None")
p5a <- p5a_domainsep

## Fig. 5C 5D 5E 5F
rpp13 <- make_tree_nolab(d, "RPP13", g = c("AT3G46530"), dist = 0.1); rpp13
rps6 <- make_tree_nolab(d, "RPS6", g = c("AT5G46470"), dist = 0.05); rps6
d4 <- make_tree_nolab(d, "DM4_RPP8", g = c("AT5G35450"), dist = 0.05); d4
d8 <- make_tree_nolab(d, "DM8"); d8
p5cd <- arrangeGrob(arrangeGrob(d8, left = textGrob(expression(bold('C')), vjust = 2,
                                                    x = unit(1,"npc"), y = unit(1, "npc"))),
                    arrangeGrob(d4, left = textGrob(expression(bold('D')), vjust = 2,
                                                    x = unit(1,"npc"), y = unit(1, "npc"))),
                    heights = c(1,1))
p5ef <- arrangeGrob(arrangeGrob(rpp13, left = textGrob(expression(bold('E')), vjust = 2,
                                                       x = unit(1,"npc"), y = unit(1, "npc"))),
                    arrangeGrob(rps6, left = textGrob(expression(bold('F')), vjust = 2,
                                                      x = unit(1,"npc"), y = unit(1, "npc"))),
                    heights = c(1, 1))

## Fig. 5
p5 <- arrangeGrob(arrangeGrob(arrangeGrob(p5a_domainsep,
                                          left = textGrob(expression(bold('A')), vjust = 2,
                                                          x = unit(1,"npc"), y = unit(1, "npc"))),
                              arrangeGrob(p_ntw + facet_wrap(~test, nrow = 1, scales = "free") +
                                          theme(legend.pos = "bottom", legend.box = "horizontal") +
                                          guides(fill = guide_legend(nrow = 1)),
                                          left = textGrob(expression(bold('B')), vjust = 2,
                                                          x = unit(1,"npc"), y = unit(1, "npc"))),
                              heights = c(1.7,1)),
                  arrangeGrob(p5cd, p5ef, heights = c(1, 1)), widths = c(3, 1))
ggsave_2s(p5, "Fig5.pdf", 150, 1.2)
ggsave_2s(p5, "Fig5.png", 150, 1.2, dpi = 300)

## Fig. S6A
p6sa <- ggplot(tmp_dat_pi_py, aes(x = count, y = pi.pairalnlen)) +
    geom_point(alpha = 0.7) +
    stat_smooth(method = "lm", formula = y ~ x) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "right", label.y.npc = 0.1, formula = y~x, parse = TRUE, size = 3) +
    geom_vline(xintercept = 64, lty = 2) +
    geom_text(aes(64, 0.2, label = "1N = 64", hjust=-0.1, vjust=-1)) + ylab("pi")

## Fig. S6B
## scatter pi against sd
tmp_dat_pi_sd <- tmp_dat_pi_py %>%
    left_join(dat_pred %>% dplyr::filter(domain == "NB-ARC") %>%
              group_by(accID, domainId) %>% summarise(count = n()) %>% ungroup() %>%
              group_by(domainId) %>% summarise(sd = sd(count), total = sum(count)),
              by = "domainId")
p6sb <- ggplot(tmp_dat_pi_sd, aes(x = sd, y = pi.pairalnlen)) +
    geom_point(aes(size = count), alpha = 0.7) +
    stat_smooth(method = "lm", formula = y~x) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 label.x.npc = "right", label.y.npc = 0.9, formula = y~x, parse = TRUE, size = 3)

## Fig. S6
p6s <- arrangeGrob(arrangeGrob(p6sa, left = textGrob(expression(bold('A')), vjust = 2,
                                                    x = unit(1,"npc"), y = unit(1, "npc"))),
                   arrangeGrob(p6sb + ylab("pi") + guides(size = guide_legend(title = "homologues")),
                               left = textGrob(expression(bold('B')), vjust = 2,
                                               x = unit(1,"npc"), y = unit(1, "npc"))),
                   heights = c(1,1))
ggsave_1s(p6s, "FigS6.pdf", 60, 2)
ggsave_1s(p6s, "FigS6.png", 60, 2, dpi = 300)


################
##  FIGURE 6  ##
######################
##  BRANCH LENGTHS  ##
##   FOR RADIATION  ##
######################

dat_bl <- read.table(paste0(dir_results, "/branch_lengths/clade_stats.tsv"),
                     header = TRUE, stringsAsFactors = FALSE, sep = '\t')

library(ggrepel)
## sd vs median

## plot how different stats change after each longest branch is sequentially split
library(data.table)
## get data
dat_split_premerge <- rbind(read.table(paste0(dir_results, "/branch_lengths/clade_split_longest_stats_mmsdNA_all.tsv"), sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>% mutate(genes = c("root")),
                            as.data.frame(rbindlist(lapply(list.files(path=paste0(dir_results, "/branch_lengths"), pattern="clade_split_longest_stats_mmsdNA_[tw]", full.names=TRUE), fread, sep='\t')),
                                          stringsAsFactors=FALSE)) %>%
    mutate(abs_path = paste(type, group, genes, path, sep = ';'),
           last_path = str_extract(abs_path, "^.+(?=[;,][^;,]+$)"),
           root = str_extract(abs_path, "^.+(?=[;][^;]+$)"),
           level = str_count(path, ','),
           num_genes = as.factor(str_count(genes, ',') + 1))
dat_split <- dat_split_premerge %>%
    left_join(dat_split_premerge %>%
              dplyr::select(c(mean, median, max, min, branch_length, level, abs_path, sd, tgenes)) %>%
              dplyr::rename(last_mean = mean, last_median = median, last_max = max, last_min = min,
                            last_level = level, last_branch_length = branch_length, last_sd = sd,
                            last_tgenes = tgenes),
              by = c("last_path" = "abs_path")) %>%
    mutate(group = sapply(group, function(x){if (x == "B5_around") {"B5*"} else {x}}))

## coerce data for plotting
coerce_dat_split <- function(dat_split, y, x = "level", add_suffix = FALSE, scale = 0.008,
                             colour_type = "type", fill_type = "status", alpha_type = "size",
                             other_cols = c()){
    dat_split_last <- dat_split[, c("abs_path", paste0("last_", x), x, paste0("last_", y), y,
                                    colour_type, fill_type, alpha_type, other_cols)]
    colnames(dat_split_last) <- c("abs_path", "last_x", "curr_x", "last_y", "curr_y",
                                  "colour_type", "fill_type", "alpha_type", other_cols)
    dat_split_to_plot <- left_join(dat_split, dat_split_last, by = "abs_path")
    dat_split_to_plot <- dat_split_to_plot %>%
        left_join(dat_split_to_plot %>% group_by(root) %>%
                  summarise(alpha_max = max(alpha_type, na.rm = TRUE)) %>% ungroup(),
                  by = "root") %>%
        mutate(alpha_type = alpha_type / alpha_max)
    dat_split_to_plot_density_ish <- dat_split_to_plot %>%
        dplyr::select(c(root, last_x, curr_x, colour_type)) %>%
        distinct() %>% group_by(colour_type, curr_x) %>% summarise(count = n()) %>%
        mutate(count = scale * count / max(count, na.rm = TRUE),
               density_type = as.factor(colour_type)) %>%
        ungroup() %>% tidyr::complete(density_type, curr_x) %>% replace_na(list(count = 0)) %>%
        mutate(density_type = as.factor(density_type))
    dat_split_to_plot_labs <- dat_split_to_plot %>% dplyr::filter(level == 0) %>%
        group_by(group, colour_type) %>% mutate(rank = rank(dplyr::desc(curr_y))) %>% ungroup()
    if (add_suffix){
        dat_split_to_plot_labs <- dat_split_to_plot_labs %>%
            mutate(suf = sapply(1:nrow(dat_split_to_plot_labs),
                                function(i){c_type <- dat_split_to_plot_labs[i, "colour_type"]; if (c_type == "radiation"){"r"} else if (c_type == "hifi") {"hf"} else {""}}),
                   suf_n = sapply(1:nrow(dat_split_to_plot_labs),
                                  function(i){if (dat_split_to_plot_labs[i, "colour_type"] == "radiation") {dat_split_to_plot_labs[i, "rank"]} else {""}}))
    }
    return(list(dat = dat_split_to_plot, labs = dat_split_to_plot_labs, dens = dat_split_to_plot_density_ish))
}

## convert coerced data into list (one metric)
make_dat_list <- function(..., add_suffix = FALSE){
    dat_split_tmp <- coerce_dat_split(..., add_suffix = add_suffix)
    if (add_suffix){
        dat_split_tmp[["labs"]] <- dat_split_tmp[["labs"]] %>% mutate(lab = paste0(group, suf, suf_n))
    } else {
        dat_split_tmp[["labs"]] <- dat_split_tmp[["labs"]] %>% mutate(lab = colour_type)
    }
    return(dat_split_tmp)
}

## convert coerced data into list (multiple metrics)
make_dat_list_multi <- function(dat, stat, args = list(), ...){
    output <- list()
    for (s in stat){
        if (s %in% names(args)){s_args <- args[[s]]} else {s_args <- list()}
        output[[s]] <- do.call(make_dat_list, c(list(dat_split = dat, y = s), s_args, list(...)))
    }
    return(output)
}

## plot data of one metric after coercion + list conversion
plot_dat_split <- function(dat_split_to_plot, dat_split_to_plot_density_ish, dat_split_to_plot_labs, stat, ...,
                           clade_label = TRUE, label_member_change = FALSE, mark_member_change = FALSE,
                           mark_member_new = FALSE, rep_plot = FALSE, minimal_legend = TRUE,
                           show_termination_type = FALSE, ylab = c(), show_survival = FALSE,
                           colour_title = "clade type", breaks = seq(0,40,10), rev = FALSE,
                           size_title = "number of\nbranches in clade",
                           lty_title = "clade type\n(survival plot)", fill_title = "termination condition",
                           alpha_title = "number of branches in clade as \nfraction of initial clade size"){
    arbitrary_args_sink <- list(...)
    library(ggrepel)
    if (rep_plot){clade_label <- FALSE}
    p <- ggplot(dat_split_to_plot) +
        geom_segment(aes(x = last_x, xend = curr_x, y = last_y, yend = curr_y, colour = colour_type,
                         size = size, alpha = alpha_type),
                     lineend = "round") +
        xlab("Number of iterations") + ylab(stat) +
        ## ggtitle(paste0(stat, " as a function of number of iterations\na clade was split by the branch with the maximum length")) +
        scale_colour_brewer(palette = "Dark2") +
        scale_x_continuous(breaks = breaks) +
        theme_bw() + 
        theme(legend.pos = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if (clade_label){
        p <- p +
            geom_text_repel(data = dat_split_to_plot_labs, ## commented out for representative plot
                            aes(label = lab, x = 0, y = curr_y, colour = colour_type),
                            hjust = 1, direction = 'y', nudge_x = -5, segment.size = 0.5, fontface = "bold",
                            show.legend = FALSE)
    }
    if (show_termination_type){
        p <- p +
           geom_point(data = dat_split_to_plot %>% filter(! str_detect(status, "max_depth|continued|terminal")),
                      aes(x = curr_x, y = curr_y, fill = fill_type), alpha = 0.7, size = 2, pch = 21) +
           scale_fill_manual(values = c("cyan", "yellow", "red")) +
           scale_shape_manual(values = c(1, 2, 4))
    }
    if (show_survival){ ## basically a survival plot, but commented out for representative plot
        p <- p +
            geom_line(data = dat_split_to_plot_density_ish,
                     aes(x = curr_x, y = count, lty = density_type, colour = density_type)) +
           scale_colour_manual(values = wes_palette("Cavalcanti1")) +
           scale_y_continuous(sec.axis = sec_axis(~ ./ scale, name = "unresolved clades"))
    }
    if (label_member_change){
        p <- p +
            geom_text_repel(data = dat_split_to_plot %>% dplyr::filter(tgenes != last_tgenes),
                            aes(label = tgenes, x = curr_x, y = curr_y, colour = colour_type))
    }
    if (mark_member_change){
        p <- p +
            geom_point(data = dat_split_to_plot %>% dplyr::filter(tgenes != last_tgenes),
                       aes(size = size, x = last_x, y = last_y), pch = 1)
    }
    if (mark_member_new){
        p <- p +
            geom_point(data = dat_split_to_plot %>% dplyr::filter(tgenes != last_tgenes),
                       aes(size = size, x = curr_x, y = curr_y), pch = 1, colour = "red")
    }
    if (!rep_plot){p <- p + expand_limits(x = c(-10,0))} ## commented out for representative plot
    if (minimal_legend){
        p <- p + guides(colour = guide_legend(ncol = 1, title = colour_title),
               lty = guide_legend(ncol = 1, lty_title),
               size = FALSE, fill = FALSE, alpha = FALSE)
    } else {
        p <- p +
            guides(colour = guide_legend(ncol = 1, title = colour_title),
                   lty = guide_legend(ncol = 1, lty_title),
                   size = guide_legend(ncol = 1, title = size_title),
                   fill = guide_legend(ncol = 1, title = fill_title),
                   alpha = guide_legend(ncol = 1, title = alpha_title),
                   size = FALSE, fill = FALSE, alpha = FALSE)
    }
    if (rev){p + aes(group = rev(colour_type))}
    return(p)
}

## plot data of two metrics after coercion + list conversion
plot_dat_split_full <- function(dat, ..., margins = unit(c(0,0,0,0), "mm"), show_legend = TRUE, ylab = c(),
                                colour_title = "clade type", legend_nrow = 2, legend_by_row = TRUE,
                                plot_legend_ratio = c(8, 1), labs = c('A', 'B'), legend_rev = FALSE,
                                legend_margin = rep(-10,4), legend_key_size = 2, legend_text_size = 9,
                                hide_legend_title = FALSE, common_plot_lab = FALSE, legend_by_gene = FALSE){
    arbitrary_args_sink <- list(...)
    stat <- names(dat)
    dat_a <- dat[[stat[[1]]]]; dat_b <- dat[[stat[[2]]]]
    if (legend_by_gene) {
        tmp_strip_suf <- function(d){d %>% mutate(colour_type = str_replace_all(colour_type, "\\.\\d+", ""))}
        dat_a <- lapply(dat_a, tmp_strip_suf); dat_b <- lapply(dat_b, tmp_strip_suf)
    }
    p1 <- plot_dat_split(dat_a[["dat"]], dat_a[["dens"]], dat_a[["labs"]],
                         stat = (if ('1' %in% names(ylab)){ylab[['1']]}else{stat[[1]]}),
                         colour_title = colour_title, ...)
    p2 <- plot_dat_split(dat_b[["dat"]], dat_b[["dens"]], dat_b[["labs"]],
                         stat = (if ('2' %in% names(ylab)){ylab[['2']]}else{stat[[2]]}),
                         colour_title = colour_title, ...)
    p_for_legend <- p1 +
        guides(colour = guide_legend(nrow = legend_nrow, title = colour_title,
                                     byrow = legend_by_row, reverse = legend_rev,
                                     override.aes = list(size = legend_key_size))) +
        theme(legend.spacing.y = unit(-2, 'mm'), legend.box.margin = margin(legend_margin),
              legend.text = element_text(size = legend_text_size),
              legend.title = if(hide_legend_title){element_blank()}else{element_text()})
    common_legend <- g_legend(p_for_legend)
    if (common_plot_lab){
        p <- arrangeGrob(arrangeGrob(p1 + theme(legend.pos = "none", plot.margin = margins),
                                     p2 + theme(legend.pos = "none", plot.margin = margins),
                                     widths = c(1, 1), padding = unit(1, "npc")),
                         left = textGrob(bquote(bold(.(labs[[1]]))), vjust = 2,
                                         x = unit(0.5, "npc"), y = unit(1, "npc")))
    } else {
        p <- arrangeGrob(arrangeGrob(p1 + theme(legend.pos = "none", plot.margin = margins),
                                     left = textGrob(bquote(bold(.(labs[[1]]))), vjust = 2,
                                                     x = unit(0.5, "npc"), y = unit(1, "npc"))),
                         arrangeGrob(p2 + theme(legend.pos = "none", plot.margin = margins),
                                     left = textGrob(bquote(bold(.(labs[[2]]))), vjust = 2,
                                                     x = unit(0.5, "npc"), y = unit(1, "npc"))),
                         widths = c(1, 1), padding = unit(1, "npc"))
    }
    if (show_legend){p <- arrangeGrob(p, common_legend, heights = plot_legend_ratio, padding = unit(1, "npc"))}
    return(p)
}

## plot data with two metrics after coercion + list conversion with tree
plot_dat_split_full_tree <- function(dat, ..., domain = '', cluster = c(), t_plot = list(),
                                     plot_tree_ratio = c(3, 1), labs = c('A', 'B', 'C'),
                                     common_plot_lab = FALSE){## pre-made tree plot can be provided using t_plot
    labs_ab <- labs[c(1,2)]
    f_tree <- function(domain, cluster, lab, t_plot, ...){
        if (is.na(t_plot)){t_plot <- make_tree_nolab(domain, cluster, ...)}
        p_tmp <- arrangeGrob(t_plot,
                             left = textGrob(bquote(bold(.(lab))), vjust = 2,
                                             x = unit(1, "npc"), y = unit(1, "npc")))
        return(p_tmp)
    }
    if (common_plot_lab) {shift <- 1} else {shift <- 2}
    p_decay <- plot_dat_split_full(dat, labs = labs_ab, common_plot_lab = common_plot_lab, ...)
    if (length(t_plot) == 0){
        glist <- sapply(1:length(cluster),
                        function(n){output <- list();
                            output[[labs[[n+shift]]]] <- f_tree(domain, cluster[[n]],
                                                                labs[[n+shift]], NA, ...);
                            return(output)})
    } else {
        glist <- sapply(1:length(t_plot),
                        function(n){output <- list();
                            output[[labs[[n+shift]]]] <- f_tree(domain, NA,
                                                                labs[[n+shift]], t_plot[[n]], ...);
                            return(output)})
    }
    p_tree <- arrangeGrob(grobs = glist,
                          heights = rep(1, max(length(cluster), length(t_plot))))
    p <- arrangeGrob(p_decay, p_tree, widths = plot_tree_ratio, padding = unit(1, "npc"))
    return(p)
}


## data for both metrics
dat_split_to_plot_rep <- dat_split %>% filter(level >= 0 & status != "terminal" &
                                              type %in% c("hifi", "radiation", "hifi_multi") &
                                              group %in% c("RPP13", "DM2"))
dat_split_to_plot_dm4 <- dat_split %>% filter(level >= 1 & status != "terminal" &
                                              type == "cluster_whole" & group == "DM4")

## coerce data for both metrics
dat_split_rep_ratio <- make_dat_list_multi(dat_split_to_plot_rep %>%
                                           dplyr::mutate(r = sd/mean, last_r = last_sd/last_mean),
                                           c("mean", "r"))
dat_split_dm4_ratio <- make_dat_list_multi(dat_split_to_plot_dm4 %>%
                                           dplyr::mutate(r = sd/mean, last_r = last_sd/last_mean),
                                           c("mean", "r"), colour_type = "tgenes")

## plot + save functions
margin1 <- unit(c(5,2,1,-5), "mm")
margin2 <- unit(c(1,2,0,0), "mm") ## for plots with common labels
plot_decay_nbs_preset1 <- function(dat, ..., cluster = c(), t_plot = list(), only_arath = TRUE, rep_plot= TRUE,
                                   ggsave_f = ggsave_1s, save = TRUE, h = 30, s = 2.2, suf = "test",
                                   margins = margin1, labs_start = 1, mark_member_change = TRUE){
    labs <- LETTERS[labs_start:(labs_start + max(length(cluster), length(t_plot)) + 1)]
    p <- plot_dat_split_full_tree(dat, domain = "NB-ARC", cluster = cluster, t_plot = t_plot,
                                  only_arath = only_arath, mark_member_change = mark_member_change,
                                  rep_plot = rep_plot, margins = margins, labs = labs, ...)
    if (save){ggsave_f(p, paste0("hifi_rad_col_", suf, ".pdf"), h, s = s)}
    return(p)
}
plot_decay_nbs_preset2 <- function(...){ ## for legend title = members
    return(plot_decay_nbs_preset1(..., colour_title = "members"))
}
plot_decay_nbs_preset3 <- function(..., plot_legend_ratio = c(4,1)){ ## for 2 row legends+legend title = members
    return(plot_decay_nbs_preset2(..., legend_nrow = 2, plot_legend_ratio = plot_legend_ratio))
}
plot_decay_nbs_preset4 <- function(...){ ## for 1 row legends + legend title = members
    return(plot_decay_nbs_preset2(..., legend_nrow = 1))
}

## Figs. 6B 6C, S7A, S7B
p_rpp13_col <- make_tree_nolab(d = "NB-ARC", c = "RPP13", only_arath = TRUE,
                               h1_col = palette(brewer.pal(n=2,name = "Dark2"))[[2]])
p_dm2_col <- make_tree_nolab(d = "NB-ARC", c = "DM2_RPP1", only_arath = TRUE,
                             g = "AT3G44630", dist = 0.1,
                             h1_col = palette(brewer.pal(n = 4, name = "Dark2"))[[3]],
                             h2_col = palette(brewer.pal(n = 4, name = "Dark2"))[[1]])
p_rpp13_bs <- make_tree_nolab(d = "NB-ARC", c = "RPP13", only_arath = TRUE,
                              colour_bootstrap = TRUE)
p_dm2_bs <- make_tree_nolab(d = "NB-ARC", c = "DM2_RPP1", only_arath = TRUE,
                            colour_bootstrap = TRUE)

## Fig. 6E (manual colouring of tips + AT5G48620 clades to accommodate paraphyly in radiation), S7C
d4_20_tips <- as.data.frame(rbindlist(lapply(list.files(path=paste0(dir_results, "/plots_final/Fig6_tips/AT5G48620"), pattern="AT", full.names=TRUE), fread, sep='\t', header = FALSE), idcol = "origin"), stringsAsFactors=FALSE) %>% group_split(origin) %>% lapply(function(x){as.vector(x[,2]) %>% pull(V1)})
d4_50_tips <- as.data.frame(rbindlist(lapply(list.files(path=paste0(dir_results, "/plots_final/Fig6_tips"), pattern="AT5G35450", full.names=TRUE), fread, sep='\t', header = FALSE)), stringsAsFactors=FALSE)[,1]
d4_70_tips <- as.data.frame(rbindlist(lapply(list.files(path=paste0(dir_results, "/plots_final/Fig6_tips"), pattern="AT5G43470", full.names=TRUE), fread, sep='\t', header = FALSE)), stringsAsFactors=FALSE)[,1]
d4_20_single <- d4_20_tips[lapply(d4_20_tips, length) == 1]
d4_20_multi <- d4_20_tips[lapply(d4_20_tips, length) > 1]
d4_tips <- list(unlist(d4_20_tips), d4_50_tips, d4_70_tips) %>% unlist
p_dm4_col_tmp <- make_tree_nolab(d = "NB-ARC", c = "DM4_RPP8", only_arath = TRUE,
                                 otu_seq = c(list("AT5G35450" = d4_50_tips, "AT5G43470" = d4_70_tips),
                                             d4_20_multi),
                                 otu_colours = c(as.vector(palette(brewer.pal(n = 4, name ="Dark2"))[c(4,3,1)]),
                                                 rep(as.vector(palette(brewer.pal(n = 4, name = "Dark2"))[[2]]),
                                                     length(d4_20_multi))), show_ref = FALSE)
cvec <- list("2" = "#E7298A", "7" = "#D95F02", "5" = "#1B9E77")
tmp <- data.frame(taxa = c(unlist(d4_20_tips), d4_50_tips, d4_70_tips),
                  did = rep(c("2", "5", "7"),
                            unlist(lapply(list(unlist(d4_20_tips), d4_50_tips, d4_70_tips), length)))) %>%
    mutate(col = sapply(did, function(x){cvec[[as.character(x)]]})) %>%
    distinct()
tmp2 <- p_dm4_col_tmp$data %>% dplyr::filter(isTip) %>% left_join(tmp, by = c("label" = "taxa"))
p_dm4_col <- p_dm4_col_tmp + geom_tippoint(colour = I(tmp2$col), size = 0.2) +
    geom_tippoint(aes(subset=(str_detect(label, "Col-0_ref"))),fill="cyan",shape=21)
p_dm4_bs <- make_tree_nolab(d = "NB-ARC", c = "DM4_RPP8", only_arath = TRUE,
                            colour_bootstrap = TRUE)

## Fig. S7
p7s_legend <- g_legend(p_rpp13_bs)
p7s <- arrangeGrob(p_dm2_bs + theme(legend.pos = "none"),
                   p_rpp13_bs + theme(legend.pos = "none"),
                   p_dm4_bs + theme(legend.pos = "none"),
                   p7s_legend,
                   heights = c(3, 3, 3, 1))
ggsave_1s(p7s, "FigS7.pdf", 200, 1)
ggsave_1s(p7s, "FigS7.png", 200, 1, dpi = 300)

## Fig. 6ABCDE
save <- FALSE
p_decay_dm4_ratio <- arrangeGrob(plot_decay_nbs_preset1(dat_split_rep_ratio, save = save,
                                                        t_plot = list("DM2" = p_dm2_col,
                                                                      "RPP13" = p_rpp13_col),
                                                        common_plot_lab = TRUE, margins = margin2,
                                                        plot_legend_ratio = c(5,1), ylab = c('2' = "sd/mean"),
                                                        legend_nrow = 1, mark_member_change = FALSE),
                                 plot_decay_nbs_preset3(dat_split_dm4_ratio %>%
                                                        lapply(function(x){list("labs"=x[["labs"]],
                                                                                "dens"=x[["dens"]],
                                                                                "dat"=x[["dat"]] %>%
                                                                                    arrange(-size))}),
                                                        save = save,
                                                        t_plot = list("DM4" = p_dm4_col), rev = TRUE,
                                                        common_plot_lab = TRUE, margins = margin2,
                                                        plot_legend_ratio = c(3.5,1), labs_start = 4,
                                                        legend_by_gene = TRUE, ylab = c('2' = "sd/mean")),
                                 heights = c(1.2,1))
p6abcde <- p_decay_dm4_ratio

## p6f
library(ggridges)
p_ridge_dm4 <- ggplot(dat_dists %>%
                      dplyr::filter(clustera == "DM4_RPP8" & clusterb == clustera &
                                    srca == "vdw" & srcb == srca) %>%
                      drop_na() %>% mutate(geneb = as.factor(geneb)),
    ##                   aes(x = distance, height = ..count..)) + ## count
    ## geom_density_ridges(aes(y = reorder(geneb, dplyr::desc(geneb)), fill = geneb),
    ##                     alpha = 0.8, scale = 2, size = 0.2, stat = "density") +
                      aes(x = distance)) + ## density
    geom_density_ridges(aes(y = reorder(geneb, dplyr::desc(geneb)), fill = geneb),
                        alpha = 0.8, scale = 2, size = 0.2) +
    geom_density(aes(y = ..density../2.5)) + facet_wrap(~genea, nrow = 1) +
    scale_fill_manual(values = as.vector(palette(brewer.pal(n = 4, name = "Dark2"))[c(1,2,4)])) +
    theme(legend.pos = "none", axis.title.y = element_blank())
p6f <- p_ridge_dm4

p6_new_r <- arrangeGrob(p6abcde,
                        arrangeGrob(p_ridge_dm4,
                                    left = textGrob('F', vjust = 2,
                                                    x = unit(0.5, "npc"), y = unit(1, "npc"))),
                        heights = c(7,2))

p6 <- p6_new_r
ggsave_1s(p6_new_r, "Fig6.pdf", h = 80, s = 2)
ggsave_1s(p6_new_r, "Fig6.png", h = 80, s = 2, dpi = 300)



################
##  TABLE S4  ##
##################################
##  MULTIPLE LINEAR REGRESSION  ##
##     (now with altitude)      ##
##################################
## TABLE S4 (GLM)

## fit model and coerce summary data into to dataframe
glm_c_df <- data.frame()
glm_c_models <- list()
for (c in (dat_pred_sum$cluster %>% unique)){
    print(c)
    tmp_dat <- dat_pred_sum %>% filter(cluster == c & accID != 6909) %>%
        left_join(acc_map, by = "accID") %>% drop_na()
    all_interactions <- count ~ longitude*latitude*altitude
    glm_c_models[[c]] <- step(glm(all_interactions, data = tmp_dat))
    tmp_df <- summary(glm_c_models[[c]])$coefficients %>% as.data.frame() %>% rownames_to_column("predictor")
    colnames(tmp_df) <- c("predictor", "estimate", "stderror", "tvalue", "Pabstvalue")
    tmp_df <- cbind(cluster = c(c), tmp_df)
    glm_c_df <- rbind(glm_c_df, tmp_df)
}
## write table
write.table(glm_c_df, paste0(dir_results, "/tables/cluster_size_longlatalt.tsv"),
            sep = '\t', row.names = FALSE, quote = FALSE)


#################
##  FIGURE S4  ##
##########################
##  PAIRED CLADES TREE  ##
##########################

## RPS4 (paired NLRs) tree
library(ape)
d <- "NB-ARC"
t_fname <- paste0(dir_results, "/tree/nlr164_col0-AL70-Alyrata_",d,"_mafft_ML.nwk")
t <- read.tree(file = t_fname)
all_ref_leaves <- c(names(readDNAStringSet(paste0(dir_results, "/domain_seq/nlr164_", d, "_ref_CDS_complete.fasta"))), names(readDNAStringSet(paste0(dir_results, "/domain_seq/nlr164_", d, "_6909pred_ref_CDS_complete.fasta"))))
tmp_dat_pred <- dat_pred %>% filter(domain == d)
tmp_dat_pred_araly <- dat_pred_araly %>% filter(domain == d)
paired_c <- c("RPS4", "cAT2G17050", "cAT4G12010", "cAT4G36140", "CHS3", "LCD9", "RPP2", "TTR1")
c_genes <- genes %>% filter(cluster %in% paired_c) %>% pull(gene)
## get leaves
al_leaves <- tmp_dat_pred %>% filter(gene %in% c_genes) %>% pull(contig)
araly_leaves <- tmp_dat_pred_araly %>% filter(gene %in% c_genes) %>% pull(contig)
ref_leaves <- data.frame(seq_id = all_ref_leaves) %>%
    mutate(gene = sapply(seq_id, gen_split_extract('\\|', 2))) %>%
    filter(gene %in% c_genes) %>% pull(seq_id) %>% as.character()
## get paired tree (includes non-paired clades as well)
paired_mrca <- getMRCA(t, c(al_leaves, araly_leaves, ref_leaves))
paired_t <- extract.clade(t, paired_mrca)
## remove non-paired clade
library(geiger)
non_paired_tips <- data.frame(seq_id = all_ref_leaves) %>%
    mutate(gene = sapply(seq_id, gen_split_extract('\\|', 2))) %>%
    filter(gene %in% c("AT5G48770", "AT5G41740")) %>% pull(seq_id) %>% as.character()
paired_t <- drop.tip(paired_t, tips(paired_t, getMRCA(paired_t, non_paired_tips)))
## regenerate al_leaves, araly_leavs, and ref_leaves to include non-paired genes in clade
all_leaves <- tips(paired_t, getMRCA(paired_t, c(al_leaves, araly_leaves, ref_leaves)))
al_leaves <- tmp_dat_pred %>% filter(contig %in% all_leaves) %>% pull(contig)
araly_leaves <- tmp_dat_pred_araly %>% filter(contig %in% all_leaves) %>% pull(contig)
ref_leaves <- all_ref_leaves[all_ref_leaves %in% all_leaves]

## organise clades by gene
groups <- list()
genes_in_clades <- tmp_dat_pred %>% filter(contig %in% all_leaves) %>% pull(gene) %>% unique()

for (g in genes_in_clades){
    g_al_leaves <- tmp_dat_pred %>% filter(gene == g & contig %in% paired_t$tip) %>% pull(contig)
    g_ref_leaves <- data.frame(seq_id = all_ref_leaves) %>%
        mutate(gene = sapply(seq_id, gen_split_extract('\\|', 2))) %>%
        filter(gene == g) %>% pull(seq_id) %>% as.character()
    groups[[g]] <- c(g_al_leaves, g_ref_leaves)
}

## assign colours
library(scales)
clusters_in_clades <- genes %>% filter(gene %in% names(groups)) %>% pull(cluster) %>% unique
cluster_non_singletons <- clusters_in_clades[! str_detect(clusters_in_clades, "singleton")]
cluster_colours <- list("singleton" = "darkgrey")
for (i in 1:length(cluster_non_singletons)){
    cluster_colours[[cluster_non_singletons[[i]]]] <- hue_pal()(length(cluster_non_singletons))[[i]]
}

gene_colours <- list()
for (g in names(groups)){
    gene_colours[[g]] <- cluster_colours[[(genes %>% filter(gene == g) %>% pull(cluster))[[1]]]]
}
gene_colours[['0']] <- "black"

## plot
library(ggtree)
library(scales)
library(ggnetwork)

paired_p <- ggtree(groupOTU(paired_t, groups, "gene")) + aes(colour = gene) +
    scale_colour_manual(values = unlist(gene_colours)) +
    geom_tippoint(aes(subset=(node %in% which(paired_t$tip %in% ref_leaves))),fill="cyan",shape=21) +
    geom_tippoint(aes(subset=(node %in% which(paired_t$tip %in% araly_leaves))),fill="red",shape=21)

for (g in names(groups)){
    g_leaves <- groups[[g]][groups[[g]] %in% all_leaves]
    paired_p <- paired_p +
        geom_cladelabel(node = getMRCA(paired_t, g_leaves), align = TRUE,
                        label = (genes %>% filter(gene == g) %>% pull(cluster))[[1]],
                        colour = gene_colours[[g]], barsize = 2)
}

## show bootstrap values
ttr_rps4_mrca1 <- getMRCA(paired_t, ref_leaves[str_detect(ref_leaves, "AT5G45260|AT5G45050")])
ttr_rps4_mrca2 <- getMRCA(paired_t, ref_leaves[str_detect(ref_leaves, "AT5G45250|AT5G45060")])
paired_root <- getMRCA(paired_t, c(ttr_rps4_mrca1, ttr_rps4_mrca2))
nodes_to_root <- c(nodepath(paired_t, paired_root, ttr_rps4_mrca1),
                   nodepath(paired_t, paired_root, ttr_rps4_mrca2)) %>% unique()

paired_p <- paired_p +
    geom_nodelab(aes(subset = node %in% nodes_to_root), vjust = -.5, hjust = 1.2, colour = "grey45") +
    expand_limits(x = c(-0.05, 1))

ggsave_2(paired_p, "FigS4.pdf", 200)
ggsave_2(paired_p, "FigS4.png", 200, dpi = 300)

