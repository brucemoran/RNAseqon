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
    output_dir <- paste0(getwd(), "/", tag, "/RNAseqR")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ##prepare data and save
  sot <- RNAseqR::brucemoran_rnaseq_kallisto_parser(metadata_csv,
                                           data_dir = data_dir,
                                           genome_prefix = genome_prefix)

  ##create required inputs to modules
  count_data <- so_to_raw_counts(sot[[1]])
  tpm_tb <- tpm_tb <- sot[[1]]$obs_raw_tpm$wide
  anno_tb <- tibble::as_tibble(sot[[2]])
  metadata_tb <- RNAseqR::get_metadata(metadata_csv, data_dir)
  outdir <- paste0(output_dir, "/RData")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  save(count_data, tpm_tb, anno_tb, metadata_tb, file = paste0(outdir, "/", tag, ".count_tpm_anno_metadata.RData"))

  ##run modules
  RNAseqR::DESeq2_module(count_data = count_data,
                anno_tb = anno_tb,
                tpm_tb = tpm_tb,
                tag = tag,
                metadata_csv = metadata_csv,
                metadata_design = metadata_design,
                output_dir = output_dir,
                control_reference = control_reference,
                delim_samples = "\\.")
  RNAseqR::edgeR_module(count_data = count_data,
                anno_tb = anno_tb,
                tpm_tb = tpm_tb,
                tag = tag,
                metadata_csv = metadata_csv,
                metadata_design = metadata_design,
                output_dir = output_dir,
                control_reference = control_reference,
                delim_samples = "\\.")
  RNAseqR::limma_module(count_data = count_data,
                anno_tb = anno_tb,
                tpm_tb = tpm_tb,
                tag = tag,
                metadata_csv = metadata_csv,
                metadata_design = metadata_design,
                output_dir = output_dir,
                control_reference = control_reference,
                delim_samples = "\\.",
                run_voom = TRUE)

  ##create master list from module RDS files
  master_list <- RNAseqR::master_parse_join(input_dir = output_dir)

  ##find overlaps
  fitwo <- RNAseqR::found_in_two(master_list)
  fithree <- RNAseqR::found_in_three(master_list)

  ##Venn results WIP
  RNAseqR::venn_3(master_list = master_list,
                  tag = tag,
                  output_dir = output_dir)

  ##fgsea
  ##run on limma, DESeq2 output per contrast
  fgsea_list_limma_rank_fithree <- lapply(names(master_list[["limma"]]), function(f){
    RNAseqR::fgsea_plot(res = master_list[["limma"]][[f]],
                       sig_res = fithree[[f]],
                       msigdb_species = msigdb_species,
                       msigdb_cat = msigdb_cat,
                       gene_col = NULL,
                       padj = 0.01,
                       output_dir = output_dir,
                       tag = f,
                       contrast = f)
  })
  fgsea_list_deseq2_stat_fithree <- lapply(names(master_list[["DESeq2"]]), function(f){
    RNAseqR::fgsea_plot(res = master_list[["DESeq2"]][[f]],
                       sig_res = fithree[[f]],
                       msigdb_species = msigdb_species,
                       msigdb_cat = msigdb_cat,
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
       fgsea_list_limma_rank_fithree_master,
       fgsea_list_deseq2_stat_fithree_master,
       file = paste0(outdir, "/", tag, ".full_results.RData"))

  ##per contrast DE overlap with pathways, and gene sets in lists
  pc_fgsea_limma_de_list <- per_contrast_fgsea_de(fgsea_list_limma_rank_fithree_master, occupancy = 5)
  names(pc_fgsea_limma_de_list) <- c("pc_fgsea_limma_de_list", "pc_genesets_limma_de_list")

  ##use these as input to ssGSEA in GSVA
  ##iterate over contrasts, using DE genesets in each pathway found associated
  ##first set up log2TPM for all genes (NB all master_list[[x]] have same geneset)
  limma_log2tpm_mat <- dplyr::select(.data = master_list[["limma"]][[1]],
                                 external_gene_name,
                                 tidyselect::ends_with("_tpm")) %>%
                   dplyr::mutate(dplyr::across(is.numeric, log2)) %>%
                   as.data.frame() %>%
                   tibble::column_to_rownames("external_gene_name") %>%
                   as.matrix()

  ##rotation_PCA plots use ssGSEA from genesets specified in pc_fgsea
  ##this is per level of the condition used for contrasts
  pc_ssgsea_list <- lapply(names(pc_fgsea_limma_de_list[[2]]), function(f){
                      ssgsea_pca_list <- ssgsea_pca(cont_pways = pc_fgsea_limma_de_list[[2]][[f]],
                                                    log2tpm_mat = limma_log2tpm_mat,
                                                    msigdb_cat = "H",
                                                    output_dir = output_dir,
                                                    hallmark_tb = hallmarks,
                                                    contrast = f)
                    })
}
