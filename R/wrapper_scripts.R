#! R

#' Wrap prep and module for data from brucemoran RNAseq_kallisto Nextflow pipeline
#'
#' @param metadata_csv path to file holding metadata
#' @param metadata_design string of design using colnames of metadata_csv
#' @param tag string used to prefix output
#' @param output_dir path to where output goes
#' @param data_dir path to where data recurse search begins
#' @param control_reference string indicating the 'control' intgroup for last element of design
#' @param sleuth_object pre-parsed sleuth_object via brucemoran_rnaseq_kallisto_parser
#' @param anno_tb string indicating the 'control' intgroup for last element of design
#' @param server_threads set to anything but NULL to allow threading
#' @return none rds and tsv files printed
#' @importFrom magrittr '%>%'
#' @export

run_prep_modules_bm <- function(metadata_csv, metadata_design, tag, output_dir, data_dir = NULL, control_reference = NULL, sleuth_object = NULL, anno_tb = NULL, server_threads = NULL) {

  ##prepare data and save
  if(is.null(sleuth_object)){
    sot <- brucemoran_rnaseq_kallisto_parser(metadata_csv, data_dir = data_dir)
  } else {
    sot <- list(sleuth_object, anno_tb)
  }

  ##check obs_raw_count, obs_raw_tpm exist
  if(!exists("sot[[1]]$obs_raw_tpm$wide")){
    sot[[1]] <- so_obs_raw_out(sot[[1]])
  }

  ##create required inputs to modules
  count_data <- so_to_raw_counts(sot[[1]])
  tpm_tb <- tpm_tb <- sot[[1]]$obs_raw_tpm$wide
  anno_tb <- tibble::as_tibble(sot[[2]])

  ##run modules
  DESeq2_module(count_data = count_data,
                anno_tb = anno_tb,
                tpm_tb = tpm_tb,
                tag = tag,
                metadata_csv = metadata_csv,
                metadata_design = metadata_design,
                output_dir = output_dir,
                control_reference = control_reference,
                delim_samples = "\\.")
  edgeR_module(count_data = count_data, anno_tb = anno_tb, tpm_tb = tpm_tb, tag = tag, metadata_csv = metadata_csv,  metadata_design = metadata_design, output_dir = output_dir, control_reference = control_reference, delim_samples = "\\.")
  limma_module(count_data = count_data, anno_tb = anno_tb, tpm_tb = tpm_tb, tag = tag, metadata_csv = metadata_csv,  metadata_design = metadata_design, control_reference = control_reference, output_dir = output_dir, delim_samples = "\\.")

  ##create master list from module RDS files
  master_list <- master_parse_join(output_dir)

  ##find overlaps
  # fitwo_list <- found_in_two(master_list)
  # fithree <- found_in_three(master_list)

  ##Venn results WIP
  venn_3(master_list = master_list, tag = tag, output_dir = output_dir)

  ##fgsea
  ##run on DESeq2 output per contrast
  fgsea_list <- lapply(names(master_list[["DESeq2"]]), function(f){
    fgsea_plot(res = master_list[["DESeq2"]][[f]], msigdb_species = "Homo sapiens", msigdb_cat = "H", gene_col = NULL, padj = 0.01)
  })
}
