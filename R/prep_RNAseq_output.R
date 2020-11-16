#! R

#' Parse featurecounts input as a follow-on from standard nf-core/rnaseq pipeline (v1.4 currently)
#'
#' @param tag string to prefix output
#' @return nz_fc_co raw count object with no lines summing to zeros
#' @importFrom magrittr '%>%'
#' @export

nf_core_rnaseq_featco_parser <- function(tag = NULL){

  if(is.null(tag)){
    tag <- "nf-core_rnaseq"
  }

  ##read in merged featureCounts results
  ##parse names using "." for nicer colnames
  fc_merge <- read_tsv("results/featureCounts/merged_gene_counts.txt")
  coln_samples <- colnames(fc_merge[3:length(colnames(fc_merge))])
  n_coln_samples <- unlist(lapply(coln_samples, function(f){
    strsplit(f, delim_samples)[[1]][1]
  }))
  colnames(fc_merge) <- c("ensembl_gene_id",
                          "external_gene_name",
                          n_coln_samples)
  on_coln_samples <- n_coln_samples[order(n_coln_samples)]
  fc_merge_so <- fc_merge %>% dplyr::select(1, 2, on_coln_samples) %>%
           arrange(ensembl_gene_id)
  write_tsv(fc_merge_so, paste0(output_dir, "/", tag, ".merged_gene_counts.fc.tsv"))

  fc_co <- fc_merge_so %>%
           dplyr::select(1, on_coln_samples) %>%
           as.data.frame() %>%
           column_to_rownames("ensembl_gene_id")

  ##remove zero-count genes
  nzero <- function(x){
    x[apply(x,1,sum) > 0,]
  }
  nz_fc_co <- nzero(fc_co)
  return(nz_fc_co)
}

#' Take metadata CSV, parse Kallisto output and return a Sleuth object with raw counts
#'
#' @param metadata_csv CSV format file with columns:
#'        'sample' that has exact dir matching which contains 'adundance.h5' (recursively)
#' @param agg_col column on which to aggregate gene_mode
#' @param genome_prefix string of a genome in biomart$dataset, suffixed with '_gene_ensembl'
#' @param server set to anything but NULL to allow threading
#' @return list object with sleuth object 'so' [[1]] and tx2gene annotation used to create it [[2]]
#' @importFrom magrittr '%>%'
#' @export

brucemoran_rnaseq_kallisto_parser <- function(metadata_csv, agg_col = "ensembl_gene_id", genome_prefix = "hsapiens", server_threads = NULL){

  ##how much comp resources to allocate
  if(!is.null(server_threads)){
    mem_free <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE))
    ram_per_core_gb <- 8
    thread_alloc <- round((memfree/1000000)/ram_per_core_gb, 0)-1
  } else {
    thread_alloc <- 1
  }

  ##parse in kallisto output to metadata
  metadata <- get_metadata(metadata_csv)

  ##annotation to use
  print("Getting tx2gene object")
  tx2gene <- get_tx2gene(genome_prefix)

  ##close open h5 files
  rhdf5::h5closeAll()
  so <- sleuth::sleuth_prep(metadata,
                            gene_mode = TRUE,
                            target_mapping = tx2gene,
                            aggregation_column = agg_col,
                            num_cores = thread_alloc,
                            max_bootstrap = 50)
  so <- so_obs_raw_out(so)
  rhdf5::h5closeAll()

  return(list(so, tx2gene))
}

#' Get annotation in tx2gene format; returns hsapiens_homolog when non-hsapiens input
#'
#' @param genome_prefix string of a genome in biomart$dataset, suffixed with '_gene_ensembl'
#' @return tx2gene format for sleuth_prep() in get_abundance_tsv()
#' @export

get_tx2gene <- function(genome_prefix){

  datasets <- biomaRt::listDatasets(biomaRt::useMart("ensembl"))

  if(! paste0(genome_prefix, "_gene_ensembl") %in% datasets$dataset){
    stop(paste0("'genome_prefix = ", genome_prefix, "' was not available, please use one of: ", datasets$dataset))
  } else {
    mart <- biomaRt::useMart(biomart = "ensembl", dataset = paste0(genome_prefix, "_gene_ensembl"))

    ##in case non-human, get ortholog
    if(genome_prefix != "hsapiens"){
      tx2gene <- tibble::as_tibble(biomaRt::getBM(attributes=c("ensembl_transcript_id",
                                             "ensembl_gene_id",
                                             "external_gene_name",
                                             "hsapiens_homolog_ensembl_gene",
                                             "hsapiens_homolog_associated_gene_name"),
                                   mart = mart))
    } else {
      tx2gene <- biomaRt::getBM(attributes=c("ensembl_transcript_id",
                                             "ensembl_gene_id",
                                             "external_gene_name"),
                                mart = mart)
    }
    colnames(tx2gene)[1] <- "target_id"
    return(tx2gene)
  }
}

#' Make metadata from CSV input, including generation of 'path' for abundance.h5 kallisto ouput
#' main job is recursive dir() to find abundance.h5 files
#' searches for dir with sampleID, then recursive-searches that for abundance.h5
#'
#' @param metadata_csv CSV format file with sample and ontehr columns
#' @export

get_metadata <- function(metadata_csv){

  ##find sampleID dirs, if none prompt to ensure same naming
  metadat <- readr::read_csv(metadata_csv)

  if(length(metadat$sample) == 0){
    stop("Please name a column 'sample'")
  } else {
    ab_list <- lapply(metadat$sample, function(f){
      files <- dir(recursive = TRUE, full.names = TRUE)
      abh5 <- grep("abundance.h5",
                   grep(paste0("/", f, "/"), files, value = TRUE),
                   value = TRUE)
      if(length(abh5) == 1){
        return(abh5)
      } else {
        print("Found multiple or no files: ")
        print(abh5)
        stop(paste0("Please rename sample contents to reflect single samples, or ensure 'abundance.h5' files are in the directory from which you are working (", getwd(), ")"))
      }
    })
    names(ab_list) <- metadat$sample
    ab_df <- data.frame(sample = names(ab_list), path = unlist(ab_list))
    metadata <- dplyr::left_join(metadat, ab_df)
    metadata$path <- as.character(metadata$path)
    return(metadata)
  }
}

#' Add wide and long raw count tibbles to sleuth object
#' sums transcripts from a gene-mode sleuth object into counts of those genes
#'
#' @param so sleuth object from sleuth_prep in gene_mode
#' @param which_value which 'values_from' to use for generating output (default: 'est_counts', else 'tpm')
#' @param group_col string of column to use in aggregating ("ensembl_gene_id" if NULL)
#' @return so which has $obs_raw_count$wide and $long
#' @importFrom magrittr '%>%'
#' @export

so_obs_raw_out <- function(so, group_col = NULL, which_value = "est_counts"){

    ##use first of the target_mapping if group_col NULL
    if(is.null(group_col)){
      group_col <- colnames(so$target_mapping)[!colnames(so$target_mapping) %in% "target_id"][1]
    }
    print(paste0("group_col -> ", group_col))

    ##cols to include as samples
    long_col <- so$sample_to_covariates$sample

    if(class(so) != "sleuth"){
        stop("Require a sleuth class object as input")
    }
    else{
      if(so$gene_mode != "TRUE"){
        stop("Require a sleuth object in gene_mode")
      }
      else{

        #est_counts
        so$obs_raw_count$wide <- tibble::as_tibble(obs_raw_widen(so, which_value = "est_counts", group_col))
        so$obs_raw_count$long <- so$obs_raw_count$wide %>%
                                 tidyr::pivot_longer(cols = dplyr::all_of(long_col),
                                                     names_to = "sample",
                                                     values_to = which_value) %>%
                                 dplyr::select(1,2,4,3)

        #tpm
        so$obs_raw_tpm$wide <- tibble::as_tibble(obs_raw_widen(so, which_value = "tpm", group_col))
        so$obs_raw_tpm$long <- so$obs_raw_tpm$wide %>%
                               tidyr::pivot_longer(cols = dplyr::all_of(long_col),
                                                   names_to = "sample",
                                                   values_to = which_value) %>%
                               dplyr::select(1,2,4,3)
        return(so)
    }
  }
}

#' Create wide tibble from sleuth object
#' sums transcripts from a gene-mode sleuth object into counts of those genes
#'
#' @param so sleuth object from sleuth_prep in gene_mode
#' @param which_value which 'values_from' to use for generating output (default: 'est_counts', else 'tpm')
#' @param group_col string of column to use in aggregating ("ensembl_gene_id" if NULL)
#' @return so which has $obs_raw_count$wide and $long
#' @importFrom magrittr '%>%'
#' @export

obs_raw_widen <- function(so, which_value, group_col) {

    ##cols to include as samples
    long_col <- so$sample_to_covariates$sample

   obs_raw_widened <- tidyr::pivot_wider(data = so$obs_raw,
                                         id_cols = target_id,
                                         names_from = sample,
                                         values_from = !!which_value) %>%
                      dplyr::inner_join(so$target_mapping, .) %>%
                      dplyr::group_by_at(group_col) %>%
                      dplyr::summarise_if(is.numeric, .funs = c(sum="sum")) %>%
                      dplyr::ungroup() %>%
                      dplyr::inner_join(so$target_mapping, .) %>%
                      dplyr::rename_at(unlist(purrr::map(long_col, dplyr::starts_with, vars = colnames(.))), list(~gsub("_sum", "", .)))  %>%
                      dplyr::select(-target_id) %>%
                      dplyr::distinct() %>%
                      dplyr::arrange(ensembl_gene_id)

  if(which_value == "est_counts"){
    obs_raw_widened <- obs_raw_widened %>%
                       dplyr::mutate_if(is.numeric, round)
  }
  return(obs_raw_widened)
}

#' so to raw_counts_df
#' @param so from sleuth object treated with so_obs_raw_out
#' @param output one of 'count' for counts, 'tpm' for tpm
#' @export

so_to_raw_counts <- function(so){

  ##remove zero count rows
  nzero <- function(x){
    x[apply(x,1,sum) > 0,]
  }

  ##df and header
  raw_df <- tibble::column_to_rownames(.data = as.data.frame(so$obs_raw_count$wide),
                                       var = "ensembl_gene_id")
  raw_df <- nzero(raw_df[,-1])
  return(raw_df)
}

#' Get StringTie TPM results
#' @importFrom magrittr '%>%'
#' @export

stringTie_res_join <- function(anno_data = NULL){

  ##FPKM/TPM from stringTie results
  st_res_tpm_list <- lapply(dir("results/stringtieFPKM", pattern="gene_abund.txt", full.names=TRUE), function(f){
    samp <- paste0(strsplit(strsplit(f,"\\/")[[1]][3], "\\.")[[1]][1], ".TPM")
    readr::read_tsv(f) %>%
    dplyr::rename(ensembl_gene_id = 1) %>%
    dplyr::select(1, !!samp := TPM) %>%
    arrange(ensembl_gene_id)
  })
  st_res_tpm <- tibble::as_tibble(join_all(st_res_tpm_list, type="left"))

  if(is.null(anno_data)){
    return(st_res_tpm)
  }
  if(!is.null(anno_data)){
    coln <- colnames(anno_data)[3]
    t1 <- st_res_tpm %>% dplyr::rename(!!coln := ensembl_gene_id)
    t2 <- dplyr::left_join(anno_data, t1)
    return(t2)
  }
}
