#! R

#' Limma module taking metadata_csv and raw count data frame to run DE analysis
#'
#' @param count_data object in matrix format (integers) with rownames as genes and colnames as sample IDs
#' @param anno_tb data frame of ensembl_gene_id - external_gene_id mappings for annotation
#' @param tpm_tb data frame of tpms to append to results for output
#' @param tag string used to prefix output
#' @param metadata_csv path to file holding metadata
#' @param metadata_design string of design using colnames of metadata_csv BN that final entry is the
#' @param control_reference string indicating the 'control' intgroup for last element of design
#' @param output_dir path to where output goes
#' @param delim_samples character to delimit sample IDs (default: ".")
#' @return nz_fc_co raw count object with no lines summing to zeros
#' @importFrom magrittr '%>%'
#' @export

limma_module <- function(count_data, anno_tb = NULL, tpm_tb = NULL, tag = NULL, metadata_csv = NULL, metadata_design = NULL, control_reference = NULL, output_dir = NULL, delim_samples = "\\."){

  print("Running: limma_module()")

  ##warnigns for undefined param
  if(is.null(tag)){
    stop("Please provide a tag to name your outputs")
  }

  if(is.null(metadata_csv)){
    stop("Please specify a metadata_csv file in CSV format")
  }

  if(is.null(metadata_design)){
    print("Please specify a metadata_design for DESeq2")
    print("NB this should be format: first + second + ... + last")
    print("NBB that 'last' in design is condition of interest in DE analysis")
    stop("NBBB that all elements of the design must be names of columns in 'metadata' file")
  }

  # if(is.null(control_reference)){
  #   print("Please specify a control_reference from a column in metadata_csv for DESeq2")
  #   print("NB that this should be one of the 'intgroup' levels")
  # }

  ##read in condition data
  cond <- readr::read_csv(metadata_csv)
  cond_df <- as.data.frame(cond)
  rownames(cond_df) <- cond_df[,1]
  # cond_df <- cond_df[,-1]

  ##break design into components
  design_vec <- unlist(lapply(metadata_design, function(f){gsub(" ", "", strsplit(f, "\\+")[[1]])}))
  CONDITION <- rev(design_vec)[1]
  cond_df[,CONDITION] <- factor(cond_df[,CONDITION])
  if(!is.null(control_reference)){
    cond_df[,CONDITION] <- relevel(cond_df[,CONDITION], ref = control_reference)
  }

  ##facotrise any other components of design_vec
  if(length(design_vec) > 1){
  for(x in 2:dim(cond_df)[2]){
    cond_df[,rev(design_vec)[x]] <- factor(cond_df[,rev(design_vec)[x]])
  }}

  ##ExpressionSet object
  # eset <- Biobase::ExpressionSet(assayData = as.matrix(count_data),
  #                                phenoData = Biobase::AnnotatedDataFrame(cond_df))
  #
  ##design
  limma_design <- model.matrix(formula(paste0("~ 0 + ", metadata_design)), data = cond_df)
  #
  # ##VOOM - Removing heteroscedascity from count data
  # pdf(paste0(output_dir, "/", tag, ".limma_voom.pdf"))
  #   limma_voom <- limma::voom(eset, limma_design, plot=TRUE)
  # dev.off()
  #
  # ##fitting design with eset
  # limma_fit <- limma::lmFit(limma_voom, limma_design)
  #
  ##make all contrasts of CONDITION, then set into named list
  contrasts <- apply(t(combn(levels(cond_df[,CONDITION]),2)), 1, function(f){
    paste(paste0(CONDITION, f), collapse="-")
  })
  limma_contrasts <- limma::makeContrasts(contrasts = contrasts,
                                          levels = paste0(CONDITION, levels(cond_df[,CONDITION])))
  # limma_cfit <- limma::contrasts.fit(limma_fit, limma_contrasts)
  # limma_fite <- limma::eBayes(limma_cfit)

  limma_dgelist <- edgeR::DGEList(counts = count_data,
                                  samples = cond_df,
                                  group = cond_df[,CONDITION])
  keep <- edgeR::filterByExpr(limma_dgelist)
  limma_dgelist <- limma_dgelist[keep,,keep.lib.sizes=FALSE]
  limma_dgelist <- edgeR::calcNormFactors(limma_dgelist)
  limma_logCPM <- edgeR::cpm(limma_dgelist, log=TRUE, prior.count=3)

  limma_fit <- limma::lmFit(limma_logCPM, limma_design)
  limma_fitc <- limma::contrasts.fit(limma_fit, limma_contrasts)
  limma_fite <- limma::eBayes(limma_fitc, trend=TRUE)

  limma_res_list <- lapply(seq_along(contrasts), function(f){
    ress_tb <- tibble::as_tibble(limma::topTable(fit = limma_fite,
                                      coef = f,
                                      number = Inf,
                                      adjust.method = "fdr"), rownames = "ensembl_gene_id")
    if(!is.null(anno_tb)){
      if(class(anno_tb)[1] != "tbl_df"){
        anno_tb <-  tibble::as_tibble(anno_tb)
      }
      if(any(colnames(anno_tb) %in% "target_id")==TRUE){
        anno_tb <- anno_tb %>% dplyr::select(-target_id)
      }
      ress_tb <- dplyr::left_join(anno_tb, ress_tb) %>%
                na.omit() %>%
                dplyr::arrange(adj.P.Val) %>%
                tibble::as_tibble() %>%
                dplyr::distinct()
    }
    if(!is.null(tpm_tb)){
      num_cols <- unlist(lapply(tpm_tb[1,], is.numeric))
      colnames(tpm_tb)[num_cols] <- paste0(colnames(tpm_tb)[num_cols], "_tpm")
      ress_tb <- dplyr::left_join(ress_tb, tpm_tb) %>%
                na.omit() %>%
                tibble::as_tibble() %>%
                dplyr::distinct() %>%
                dplyr::arrange(adj.P.Val)
    }
    ress_tb <- ress_tb %>% dplyr::rename("padj" = adj.P.Val)
    readr::write_tsv(ress_tb, paste0(output_dir, "/", tag, ".res.", contrasts[f], ".limma.tsv"))
  })
  names(limma_res_list) <- contrasts
  saveRDS(limma_res_list, file = paste0(output_dir, "/", tag, ".limma_res_list.rds"))

}
