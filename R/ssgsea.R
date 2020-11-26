#' Function to run fGSEA using the mSIGDB database from Broad
#'
#' @param res output table from 'DESeq2_module()' (all genes not just sig)
#' @param sig_res output table from 'DESeq2/limma/edger_module()' or combination thereof (used for table output)
#' @param msigdb_species one of msigdbr::msigdbr_show_species(), default:"Homo sapiens"
#' @param msigdb_cat one of 'c("H", paste0("C", c(1:7)))', see: gsea-msigdb.org/gsea/msigdb/collections.jsp
#' @param gene_col the name of the column in tibble with gene names found in pathways, set to "rownames" if they are the rownames (default)
#' @param rank_col the name of the column in tibble which ranks genes for fgsea (default: "stat" for DESeq2 results; limma - use "t", edgeR - unsure)
#' @param padj the significance threshold
#' @param output_dir path to where output goes
#' @param tag string used to prefix output
#' @param contrast string to define the contrast being made, tags output
#' @return msigdb_fgsea object
#' @export

ssgsea_run <- function(res, sig_res = NULL, msigdb_species = "Homo sapiens", msigdb_cat = "H", gene_col = NULL, rank_col = NULL, padj = 0.01, output_dir, tag, contrast) {

  print("Running: fgsea_plot()")

  ##output dir
  out_dir <- paste0(output_dir, "/fgsea")
  dir.create(paste0(out_dir, "/plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(out_dir, "/de_pathways"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(out_dir, "/all"), recursive = TRUE, showWarnings = FALSE)

  ##use res as significant set also
  if(is.null(sig_res)){
      sig_res <- res
  }

  ##name of gene column to use
  if(is.null(gene_col)){
      gene_col <- "external_gene_name"
  }

  ##determine which col to user as rank (based on colnames)
  sig_col <- ""
  if(length(grep("_padj", colnames(sig_res))) > 0){
    if("limma_padj" %in% colnames(sig_res)){
      sig_col <- "limma_padj"
    } else {
      sig_col <- "DESeq2_padj"
    }
  } else {
    ##if no '_padj' we are in single caller, it's padj (renamed by script)
    sig_col <- "padj"
  }

  ##this to guess at correct rank_col if NULL
  ##NB can't use edgeR, so stop and return warning if matches only edgeR
  if(is.null(rank_col)){
    if(length(grep("AveExpr", colnames(res))) > 0){
      if(!"t" %in% colnames(res)){
        rank_col <- "limma_t"
      } else {
        rank_col <- "t"
      }
    } else {
      if(length(grep("stat", colnames(res))) > 0){
        if(!"stat" %in% colnames(res)){
          rank_col <- "DESeq2_stat"
        } else {
          rank_col <- "stat"
        }
      } else {
        stop("Couldn't find limma or DESeq2 colnames, is this edgeR result data? IF so please use a different source, or specify 'rank_col' paramater input")
      }
    }
  }
  print(paste0("Found rank_col: ", rank_col))

  ##select/arrange by rank_col
  ##from https://stephenturner.github.io/deseq-to-fgsea/ (thanks!)
  res_rank <- dplyr::select(.data = res, !!gene_col, !!rank_col) %>%
              na.omit() %>%
              dplyr::distinct() %>%
              dplyr::group_by(!!as.symbol(gene_col)) %>%
              dplyr::summarize(rank = mean(!!as.symbol(rank_col)))
  rank_vec <- tibble::deframe(res_rank)

  ##MsigDB pathways
  msigdb_pathway <- msigdbr::msigdbr(species = msigdb_species, category = msigdb_cat)

  ##create list
  msigdb_pathlist <- lapply(unique(msigdb_pathway$gs_name), function(f){
    fsig <- dplyr::filter(.data = msigdb_pathway, gs_name %in% !!f)
    return(as.vector(unlist(dplyr::select(.data = fsig, gene_symbol))))
  })
  names(msigdb_pathlist) <- unique(msigdb_pathway$gs_name)

  ##run, order fgsea
  fgsea_res <- fgsea::fgsea(pathways = msigdb_pathlist,
                            stats = rank_vec,
                            nperm = 1000000)

  ##make tibble and add FDR, filters
  fgsea_res_tb <- dplyr::mutate(.data = tibble::as_tibble(fgsea_res))
  fgsea_res_tb <- dplyr::select(.data = fgsea_res_tb, -leadingEdge)
  readr::write_tsv(fgsea_res_tb, file = paste0(out_dir, "/all/", tag , ".fgsea_all.tsv"))

  fgsea_res_sig_tb <- dplyr::filter(.data = fgsea_res_tb, padj < !!padj)
  fgsea_res_sig_tb <- dplyr::arrange(.data = fgsea_res_sig_tb, desc(NES))

  ##plotting
  gg_fgsea <- ggplot2::ggplot(fgsea_res_sig_tb, ggplot2::aes(reorder(pathway, NES), NES)) +
              ggplot2::geom_col(ggplot2::aes(fill = padj)) +
              ggplot2::coord_flip() +
              ggplot2::labs(x = "Pathway",
                            y = "Normalized Enrichment Score",
                            title = paste0(msigdb_cat, "MsigDB pathways NES")) +
              ggplot2::theme_minimal()
  ggplot2::ggsave(gg_fgsea, file = paste0(out_dir, "/plots/", tag , ".fgsea_sig.ggplot2.pdf"))

  ##output results per gene
  pathways_sig_res_tb <- msigdb_pathlist[names(msigdb_pathlist) %in% fgsea_res_sig_tb$pathway] %>%
                     tibble::enframe("pathway", dplyr::all_of(gene_col)) %>%
                     tidyr::unnest(cols = gene_col) %>%
                     dplyr::inner_join(sig_res, by = gene_col) %>%
                     dplyr::filter(!!as.symbol(sig_col) < !!padj) %>%
                     dplyr::distinct() %>%
                     dplyr::mutate(contrast = !!contrast) %>%
                     dplyr::inner_join(fgsea_res_sig_tb[,c("pathway", "NES", "size")], by = "pathway") %>%
                     dplyr::select(1, NES, size, contrast, tidyselect::everything())
  readr::write_tsv(pathways_sig_res_tb, file = paste0(out_dir, "/de_pathways/", tag , ".fgsea_pathway.tsv"))
  return(pathways_sig_res_tb)
}

##run fgsea on each MR set individually, returning list same length as MRs, table of FGSEA
fgseamsViperIndiv <- function(mrs, mrsstat, pathway, qval, TAG, OUTDIR, gene_col=NULL) {

  ##mrs is output from msViper, with an es object with specified mrstat as name
  ##mrsstat is statistic used to rank data, comes from msViper
  ##pathway is pathway (from MSigDB in GMT format if *gmt, else named list)
  ##qval is the adjusted p value for all results herein
  ##gene_col is the name of the column in tibble with gene names found in pathways

  if(is.null(gene_col)){
    gene_col <- "external_gene_name"
  }

  if(is.null(geneset)){
    geneset <- names(mrs$es[[mrsstat]])
  }

  if(length(grep(".gmt$", pathway, perl=TRUE))==1){
    pathways_msig <- gmtPathways(pathway)
  }

  ##create ranks
  ranks <- as_tibble(data.frame(stat=mrs$es[[mrsstat]]), rownames="external_gene_name") %>%
           dplyr::filter(external_gene_name %in% geneset) %>%
           dplyr::select(gene_col, stat) %>%
           na.omit() %>%
           tibble::deframe()

  ##run, order fgsea
  fgsea_res <- fgsea(pathways=pathways_msig, stats=ranks, nperm=1000000)
  fgsea_res_le <- fgsea_res[lapply(fgsea_res$leadingEdge,length)>25,]
  fgsea_resTidy <- as_tibble(fgsea_res_le) %>%
                  dplyr::mutate(FDR = p.adjust(pval,method="BH")) %>%
                  dplyr::filter(padj < 0.1) %>%
                  dplyr::arrange(desc(NES))

  ##
  pathways_msig.DEsig.absFC2 <- pathways_msig %>% enframe("pathway", gene_col) %>%
                    unnest() %>%
                    inner_join(DESeqResults.t, by=gene_col) %>%
                    dplyr::mutate(absFC = logratio2foldchange(log2FoldChange)) %>%
                    dplyr::filter(padj < qval) %>%
                    dplyr::filter(abs(absFC) > 2) %>%
                    left_join(., fgsea_resTidy, by="pathway") %>%
                    dplyr::select(pathway, gene_col, padj.y, NES, padj.x, absFC) %>%
                    dplyr::rename(padj_pathway = "padj.y", padj_gene = "padj.x") %>%
                    dplyr::arrange(NES)

  ##plot
  pathways_msig.DEsig.absFC2.plt <- pathways_msig.DEsig.absFC2 %>%
                                    dplyr::select(pathway, NES, padj_pathway) %>%
                                    dplyr::rename(padj = "padj_pathway") %>%
                                    dplyr::filter(padj < qval) %>%
                                    distinct()

  ggplot(pathways_msig.DEsig.absFC2.plt, ggplot2::aes(reorder(pathway, NES), NES)) +
    geom_col(ggplot2::aes(fill=padj<qval)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(TAG, " pathways NES from GSEA")) +
    theme_minimal() +
    theme(axis.text.y=element_text(size=5))
  ggsave(paste0(OUTDIR, "/fgseaDESeq.", TAG, ".sig.absFC2.plt.pdf"))

  return(pathways_msig.DEsig.absFC2)
}
