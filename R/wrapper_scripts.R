#! R

#' Wrap prep and module for data from brucemoran RNAseq_kallisto Nextflow pipeline
#'
#' @param metadata_csv path to file holding metadata
#' @param metadata_design string of design using colnames of metadata_csv
#' @param tag string used to prefix output
#' @param output_dir path to where output goes
#' @param data_dir path to where data recurse search begins
#' @param control_reference string indicating the 'control' intgroup for last element of design
#' @param genome_prefix string of a genome in biomart$dataset, suffixed with '_gene_ensembl'
#' @param msigdb_species one of msigdbr::msigdbr_show_species(), default:"Homo sapiens"
#' @param msigdb_cat one of 'c("H", paste0("C", c(1:7)))', see: gsea-msigdb.org/gsea/msigdb/collections.jsp
#' @param server_threads set to anything but NULL to allow threading
#' @return none rds and tsv files printed
#' @importFrom magrittr '%>%'
#' @export

run_prep_modules_bm <- function(metadata_csv, metadata_design, tag, output_dir = NULL, data_dir = NULL, control_reference = NULL, genome_prefix = "hsapiens", msigdb_species = "Homo sapiens", msigdb_cat = "H", server_threads = NULL) {

  ##allow "NULL" to become NULL for input from Nextflow
  if(is.null(control_reference) || control_reference == "NULL"){
    control_reference <- as.null(control_reference)
  }

  ##create an output called RNAseqR in current dir if no output_dir defined
  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/", tag, "/RNAseqon")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ##prepare data and save
  sot <- RNAseqon::brucemoran_rnaseq_kallisto_parser(metadata_csv,
                                           data_dir = data_dir,
                                           genome_prefix = genome_prefix)

  ##create required inputs to modules
  count_data <- so_to_raw_counts(sot[[1]])
  tpm_tb <- sot[[1]]$obs_raw_tpm$wide

  ##make log2tpm
  log2tpm <- dplyr::mutate(.data = tpm_tb, dplyr::across(where(is.numeric), log2))

  ##aggregate ens IDs to get single ext ID for TPM
  ##NB that ens IDs are unique, but map to multiple ext IDs
  agg_log2tpm_tb <- agg_log2tpm <- RNAseqon::group_agg_two(log2tpm, pattern = "_gene")
  agg_log2tpm_df <- as.data.frame(agg_log2tpm_tb)
  rownames(agg_log2tpm_df) <- agg_log2tpm_df$external_gene_name
  agg_log2tpm_df <- agg_log2tpm_df[, ! colnames(agg_log2tpm_df) %in% c("external_gene_name", "ensembl_gene_id")]
  agg_log2tpm_ext_mat <- as.matrix(agg_log2tpm_df)

  anno_tb <- tibble::as_tibble(sot[[2]])
  metadata_tb <- RNAseqon::get_metadata(metadata_csv, data_dir)
  outdir <- paste0(output_dir, "/RData")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  save(sot, count_data, anno_tb, tag, tpm_tb, metadata_csv, metadata_design, output_dir, control_reference, agg_log2tpm_tb, agg_log2tpm_ext_mat, anno_tb, metadata_tb, file = paste0(outdir, "/", tag, ".DE_ready.RData"))

  ##run modules
  RNAseqon::DESeq2_module(count_data = count_data,
                anno_tb = anno_tb,
                tag = tag,
                metadata_csv = metadata_csv,
                metadata_design = metadata_design,
                output_dir = output_dir,
                control_reference = control_reference,
                delim_samples = "\\.")
  RNAseqon::edgeR_module(count_data = count_data,
                anno_tb = anno_tb,
                tag = tag,
                metadata_csv = metadata_csv,
                metadata_design = metadata_design,
                output_dir = output_dir,
                control_reference = control_reference,
                delim_samples = "\\.")
  RNAseqon::limma_module(count_data = count_data,
                anno_tb = anno_tb,
                tag = tag,
                metadata_csv = metadata_csv,
                metadata_design = metadata_design,
                output_dir = output_dir,
                control_reference = control_reference,
                delim_samples = "\\.",
                run_voom = TRUE)

  ##create master list from module RDS files
  master_list <- RNAseqon::master_parse_join(input_dir = output_dir)

  ##find overlaps
  fitwo <- RNAseqon::found_in_two(master_list)
  fithree <- RNAseqon::found_in_three(master_list)

  ##Venn results WIP
  RNAseqon::venn_3(master_list = master_list,
                  tag = tag,
                  output_dir = output_dir)

  ##fgsea
  ##run on limma, DESeq2 output per contrast
  fgsea_list_limma_rank_fithree <- lapply(names(master_list[["limma"]]), function(f){
    RNAseqon::fgsea_plot(res = master_list[["limma"]][[f]],
               sig_res = fithree[[f]],
               msigdb_species = msigdb_species,
               msigdb_cat = "H",
               gene_col = NULL,
               padj = 0.01,
               output_dir = output_dir,
               tag = f,
               contrast = f)
  })
  fgsea_list_deseq2_stat_fithree <- lapply(names(master_list[["DESeq2"]]), function(f){
    RNAseqon::fgsea_plot(res = master_list[["DESeq2"]][[f]],
               sig_res = fithree[[f]],
               msigdb_species = msigdb_species,
               msigdb_cat = "H",
               gene_col = NULL,
               padj = 0.01,
               output_dir = output_dir,
               tag = f,
               contrast = f)
  })
  names(fgsea_list_limma_rank_fithree) <- names(fgsea_list_deseq2_stat_fithree)<- names(master_list[["limma"]])

  ##bind those into a single table
  fgsea_list_limma_rank_fithree_master <- do.call(rbind, fgsea_list_limma_rank_fithree)
  fgsea_list_deseq2_stat_fithree_master <- do.call(rbind, fgsea_list_deseq2_stat_fithree)

  ##save outputs
  save(master_list, fitwo, fithree,
       fgsea_list_limma_rank_fithree,
       fgsea_list_deseq2_stat_fithree,
       fgsea_list_limma_rank_fithree_master,
       fgsea_list_deseq2_stat_fithree_master,
       file = paste0(outdir, "/", tag, ".full_results.RData"))

  ##per contrast DE overlap with pathways, and gene sets in lists
  pc_fgsea_limma_de_list <- per_contrast_fgsea_de(fgsea_list_limma_rank_fithree, occupancy = 5)
  names(pc_fgsea_limma_de_list) <- names(fgsea_list_limma_rank_fithree)

  ##use these as input to ssGSEA in GSVA
  ##iterate over contrasts, using DE genesets in each pathway found associated
  ##first set up log2TPM for all genes (NB all master_list[[x]] have same geneset)
  total_geneset <- unique(as.vector(unlist(lapply(pc_fgsea_limma_de_list, function(f){
    return(unlist(f[[2]]))
  }))))
  limma_log2tpm_mat <- agg_log2tpm_ext_mat[rownames(agg_log2tpm_ext_mat) %in% total_geneset,]

  ##rotation_PCA plots use ssGSEA from genesets specified in pc_fgsea
  ##this is per level of the condition used for contrasts
  metadata_cov <- rev(gsub(" ", "", strsplit(metadata_design, "\\+")[[1]]))[1]
  metadata_pca <- dplyr::select(.data = metadata_tb, sample, !!metadata_cov)

  pc_ssgsea_list <- lapply(names(pc_fgsea_limma_de_list), function(f){
                      ssgsea_pca_list <- ssgsea_pca(pways = pc_fgsea_limma_de_list[[f]][[2]],
                                                    log2tpm_mat = agg_log2tpm_ext_mat,
                                                    msigdb_cat = "H",
                                                    output_dir = output_dir,
                                                    contrast = f,
                                                    metadata = metadata_pca)
                    })
  names(pc_ssgsea_list) <- names(pc_fgsea_limma_de_list)
}
