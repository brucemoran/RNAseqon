#' Function to run fGSEA using the mSIGDB database from Broad
#'
#' @param res data frame (all genes not just sig)
#' @param sig_res data frame of significant results
#' @param msigdb_species one of msigdbr::msigdbr_show_species(), default:"Homo sapiens"
#' @param msigdb_cat one of 'c("H", paste0("C", c(1:7)))', see: gsea-msigdb.org/gsea/msigdb/collections.jsp
#' @param gene_col the name of the column in data frame with gene names found in pathways, set to "rownames" if they are the rownames (default)
#' @param rank_col the name of the column in tibble which ranks genes for fgsea (default: "stat" for DESeq2 results; limma - use "t", edgeR - unsure)
#' @param padj the significance threshold
#' @param output_dir path to where output goes
#' @param tag string used to prefix output
#' @param contrast string to define the contrast being made, tags output
#' @param plot_out logical plot output to screen (for knitr; default TRUE)
#' @return msigdb_fgsea object
#' @export

fgsea_plot <- function(res, sig_res = NULL, msigdb_species = "Homo sapiens", msigdb_cat = "H", gene_col = NULL, rank_col = NULL, padj = 0.01, output_dir, tag, contrast, plot_out = TRUE) {

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

  ##MsigDB pathways genelists
  msigdb_pathlist <- msigdb_pathways_to_list(msigdb_species, msigdb_cat)

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
  gg_fgsea <- fgsea_plotting(fgsea_res_tb = fgsea_res_sig_tb, msigdb_cat = msigdb_cat)
  if(!is.null(plot_out)){
    gg_fgsea
  }
  ggplot2::ggsave(gg_fgsea, file = paste0(out_dir, "/plots/", tag , ".fgsea_sig.ggplot2.pdf"))

  ##output results per gene
  size_list <- lapply(seq_along(msigdb_pathlist), function(f){
      fo <- data.frame(pathway = names(msigdb_pathlist)[f],
                       size = length(msigdb_pathlist[[f]]))
      return(fo)
  })
  size_tb <- tibble::as_tibble(do.call(rbind, size_list))

  pathways_sig_res_tb <- msigdb_pathlist[names(msigdb_pathlist) %in% fgsea_res_sig_tb$pathway] %>%
                         tibble::enframe("pathway", dplyr::all_of(gene_col)) %>%
                         tidyr::unnest(cols = gene_col) %>%
                         dplyr::inner_join(sig_res, by = gene_col) %>%
                         dplyr::left_join(size_tb) %>%
                         dplyr::filter(!!as.symbol(sig_col) < !!padj) %>%
                         dplyr::distinct() %>%
                         dplyr::mutate(contrast = !!contrast)
  readr::write_tsv(pathways_sig_res_tb, file = paste0(out_dir, "/de_pathways/", tag , ".fgsea_pathway.tsv"))
  return(pathways_sig_res_tb)
}

#' Function to take output table from multiple contrasts
#' take each level of contrasts and combine those levels into master tables
#' filter based on min genes per master table
#'
#' @param fgsea_contrast_list named list returned from lapply on fgsea_plot()
#' @param occupancy % of genes in the pathway that allow use of the pathway for ssgsea
#' @param output_dir path to where output goes
#' @param tag string used to prefix output
#' @return table of pathways for each contrast level
#' @export

per_contrast_fgsea_de <- function(fgsea_contrast_list, occupancy = 5, output_dir, tag) {

  pc_fg_pway_list <- lapply(fgsea_contrast_list, function(f){

    ##levels of contrast
    cont_levels <- unique(unlist(strsplit(unique(f$contrast), "-")))

    ##count gene occupancy per pathway and filter
    pc_fg_n <- dplyr::group_by(.data = f, pathway) %>%
               dplyr::tally() %>%
               dplyr::left_join(f) %>%
               dplyr::mutate(occ = 100*(n/size)) %>%
               dplyr::filter(occ > !!occupancy)

    ##second list of genesets per pathway per list
    pways <- unique(pc_fg_n$pathway)
    pways_list <- lapply(pways, function(p){
      dplyr::filter(.data = pc_fg_n, pathway %in% p) %>%
      dplyr::select(external_gene_name) %>%
      unlist() %>% as.vector() %>% unique()
    })
    names(pways_list) <- pways

    return(list(pc_fg_n, pways_list))
  })
  return(pc_fg_pway_list)
}


#' msigdb pathways in a nice list
#'
#' @param msigdb_species one of msigdbr::msigdbr_show_species()
#' @param msigdb_cat one of 'c("H", paste0("C", c(1:7)))',
#'        see: gsea-msigdb.org/gsea/msigdb/collections.jsp
#' @return msigdb_pathlist list object
#' @export

msigdb_pathways_to_list <- function(msigdb_species, msigdb_cat){

  msigdb_pathway <- msigdbr::msigdbr(species = msigdb_species, category = msigdb_cat)

  ##create list
  msigdb_pathlist <- lapply(unique(msigdb_pathway$gs_name), function(f){
    fsig <- dplyr::filter(.data = msigdb_pathway, gs_name %in% !!f)
    return(as.vector(unlist(dplyr::select(.data = fsig, gene_symbol))))
  })
  names(msigdb_pathlist) <- unique(msigdb_pathway$gs_name)
  return(msigdb_pathlist)
}

#' fgsea plotting
#'
#' @param fgsea_res_tb tibble of results from fgsea::fgsea
#' @param msigdb_cat one of 'c("H", paste0("C", c(1:7)))',
#'        see: gsea-msigdb.org/gsea/msigdb/collections.jsp
#' @return msigdb_pathlist list object
#' @export

fgsea_plotting <- function(fgsea_res_tb, msigdb_cat){
  gg_fgsea <- ggplot2::ggplot(fgsea_res_tb, ggplot2::aes(reorder(pathway, NES), NES)) +
              ggplot2::geom_col(ggplot2::aes(fill = padj)) +
              ggplot2::coord_flip() +
              ggplot2::labs(x = "Pathway",
                            y = "Normalized Enrichment Score",
                            title = paste0(msigdb_cat, "MsigDB pathways NES")) +
              ggplot2::theme_minimal()
  return(gg_fgsea)
}
