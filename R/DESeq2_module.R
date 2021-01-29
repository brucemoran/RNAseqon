#! R

#' DESeq2 module taking metadata_csv nad raw count data frame to run DE analysis
#'
#' @param count_data object in matrix format (integers) with rownames as genes and colnames as sample IDs
#' @param anno_tb tibble of ensembl_gene_id - external_gene_id mappings for annotation
#' @param tpm_tb tibble of tpms to append to results for output
#' @param tag string used to prefix output
#' @param metadata_csv path to file holding metadata
#' @param metadata_design string of design using colnames of metadata_csv
#' @param control_reference string indicating the 'control' intgroup for last element of design
#' @param output_dir path to where output goes
#' @param ens_version ENSEMBL version to use
#' @param org_prefix organism prefix from ENSEMBL
#' @param delim_samples character to delimit sample IDs (default: ".")
#' @return nz_fc_co raw count object with no lines summing to zeros
#' @importFrom magrittr '%>%'
#' @export

DESeq2_module <- function(count_data, anno_tb = NULL, tpm_tb = NULL, tag = NULL, metadata_csv = NULL, metadata_design = NULL, control_reference = NULL, output_dir = NULL, delim_samples = "\\."){

  print("Running: DESeq2_module()")

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
  #   stop("NB that this should be one of the 'intgroup' levels")
  # }

  if(is.null(output_dir)){
    print("No directory specified for output, the current dir will be used")
    output_dir <- "./DEseq2"
  } else {
    output_dir <- paste0(output_dir, "/DEseq2")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ##read in condition data
  cond <- readr::read_csv(metadata_csv)
  cond_df <- as.data.frame(cond)
  rownames(cond_df) <- cond_df[,1]

  ##overlap w/count_data and sort
  cond_df <- cond_df[rownames(cond_df) %in% colnames(count_data),]
  cond_df <- cond_df[sort(rownames(cond_df)),]

  ##break design into components
  design_vec <- unlist(lapply(metadata_design, function(f){gsub(" ", "", strsplit(f, "\\+")[[1]])}))

  ##take only cols of interest and define major condition
  cond_df <- cond_df[,colnames(cond_df) %in% c("sample", design_vec)]
  CONDITION <- rev(design_vec)[1]

  ##create factors
  for(x in 1:length(colnames(cond_df))){
    cond_df[,x] <- factor(cond_df[,x])
  }

  ##change reference level if specified
  if(!is.null(control_reference)){
    cond_df[,CONDITION] <- relevel(cond_df[,CONDITION], ref = control_reference)
  }

  ##DESeq2DataSet object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data,
                                colData = cond_df,
                                design = formula(paste0("~ ", metadata_design)))

  ##run DESeq2, make results
  ddseq <- DESeq2::DESeq(dds)

  ##make all contrasts of CONDITION, then set into named list
  combn_mat <- t(combn(levels(cond_df[,CONDITION]),2))
  combns <- apply(combn_mat, 1, function(f){paste(f, collapse="-")})

  DESeq2_res_list <- lapply(combns, function(f){
    ress <- na.omit(DESeq2::results(ddseq, contrast = c(CONDITION, strsplit(f, "-")[[1]])))
    ress_tb <- tibble::as_tibble(as.data.frame(ress), rownames="ensembl_gene_id")

    if(!is.null(anno_tb)){
      if(class(anno_tb)[1] != "tbl_df"){
        anno_tb <-  tibble::as_tibble(anno_tb)
      }
      if(any(colnames(anno_tb) %in% "target_id")==TRUE){
        anno_tb <- anno_tb %>% dplyr::select(-target_id)
      }
      ress_tb <- dplyr::left_join(anno_tb, ress_tb) %>%
                na.omit() %>%
                dplyr::arrange(padj) %>%
                tibble::as_tibble() %>%
                dplyr::distinct()
    }

    if(!is.null(tpm_tb)){
      num_cols <- unlist(lapply(tpm_tb[1,], is.numeric))
      colnames(tpm_tb)[num_cols] <- paste0(colnames(tpm_tb)[num_cols], "_tpm")
      ress_tb <- dplyr::left_join(ress_tb, tpm_tb) %>%
                na.omit() %>%
                dplyr::arrange(padj) %>%
                tibble::as_tibble() %>%
                dplyr::distinct()
    }
    readr::write_tsv(ress_tb, paste0(output_dir, "/", tag, ".res.", f, ".DESeq2.tsv"))

    return(ress_tb)
  })
  names(DESeq2_res_list) <- unlist(lapply(combns, function(f){
      paste(paste0(CONDITION, strsplit(f, "-")[[1]]), collapse = "-")
    }))
  saveRDS(DESeq2_res_list, file = paste0(output_dir, "/", tag, ".DESeq2_res_list.rds"))

  ##plots
  vsd <- DESeq2::vst(dds, blind = TRUE)
  sampleDists <- dist(t(SummarizedExperiment::assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- grDevices::colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(256)
  hc <- hclust(sampleDists)

  pdf(paste0(output_dir, "/", tag, ".heatmap.pdf"))
    heatmap(sampleDistMatrix, Rowv=as.dendrogram(hc),
            symm=TRUE, col=colors,
            margins=c(2,10), labCol=FALSE, cexRow = 0.8 )
  dev.off()

  ##PCA
  lapply(design_vec, function(f){
    bmpcaplot <- BMplotPCA(vsd, intgroup=c(f), anno_tb = anno_tb, pc_limit = 5)
      pdf(paste0(output_dir, "/", tag, ".", f, ".PCA_PC_loadings.pdf"), onefile = TRUE)
      lapply(bmpcaplot[[1]], print)
      print(bmpcaplot[[2]])
      print(bmpcaplot[[3]])
    dev.off()
  })
}

#' Plotting PCA function
#' @param x variance stabilized DESeq2 object
#' @param intgroup which colname from metadata_csv to be output on plot
#' @param ntop top n genes to use by variance
#' @param anno_tb tibble of ensembl_gene_id - external_gene_id mappings for annotation
#' @param pc_limit integer percent variance accounted by PC for inclusion in plots  (default: 10%)
#' @return list of ggplot2 objects for printing (PCA, PCs, loadings)
#' @export

BMplotPCA <- function(x, intgroup = NULL, ntop = 15000, anno_tb, pc_limit = 10) {
    rv <- matrixStats::rowVars(SummarizedExperiment::assay(x))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
        length(rv)))]
    pca <- prcomp(t(SummarizedExperiment::assay(x)[select, ]))
    sdPc <- apply(pca$x, 2, sd)
    percentVar <- sdPc^2/sum(sdPc^2)
    sdRot <- apply(pca$rotation, 2, sd)
    percentVarRot <- sdRot^2/sum(sdRot^2)
    if (!all(intgroup %in% names(SummarizedExperiment::colData(x)))) {
        stop("The argument 'intgroup' should specify columns of colData(xx)")
    }
    intgroup.df <- as.data.frame(SummarizedExperiment::colData(x)[, intgroup, drop = FALSE])
    group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))

    pca_plot <- function(d, group, intgroup.df, pcv){

      if(nlevels(group)>6){
        ggp <- ggplot2::ggplot(data = d,
                               ggplot2::aes_string(x = "PC1",
                                                   y = colnames(d)[2],
                                                   group = intgroup)) +
               ggplot2::scale_shape_manual(values = 1:dim(unique(intgroup.df))[1]) +
               ggplot2::labs(paste0("PCA plot using ", colnames(intgroup.df)), x = paste0("PC1: ", round(pcv[1] *  100), "% variance"),y = paste0(colnames(d)[2], ": ", round(pcv[2] * 100), "% variance")) +
               ggrepel::geom_text_repel(label = rownames(d),
                                        colour = "black",
                                        size = 2,
                                        fontface = "bold") +
               ggplot2::geom_point(ggplot2::aes_string(shape = intgroup,
                                                       colour = intgroup,
                                                       fill = intgroup),
                                   size = 3) +
               ggplot2::ggtitle(paste0("PCA plot using ", colnames(intgroup.df)),
                                subtitle = paste0(colnames(d)[1], " vs. ", colnames(d)[2]))
        }

        if(nlevels(group)<=6){
        ggp <- ggplot2::ggplot(data = d, ggplot2::aes_string(x = "PC1",
                            y = colnames(d)[2],
                            group = intgroup)) +
                            # ggplot2::aes(x = PC1, y = PC2, group = group, shape = group, colour = group)) +
               ggplot2::geom_point(ggplot2::aes_string(shape = intgroup,
                                                       colour = intgroup,
                                                       fill = intgroup),
                                   size = 3) +
               ggplot2::scale_shape_discrete(solid = T) +
               ggplot2::xlab(paste0("PC1: ", round(pcv[1] *  100), "% variance")) +
               ggplot2::ylab(paste0(colnames(d)[2], ": ", round(pcv[2] * 100), "% variance")) +
               ggrepel::geom_text_repel(label = rownames(d),
                                        colour = "black",
                                        size = 2,
                                        fontface = "bold") +
                                        # ggplot2::annotate("text",x=pca$x[,1], y = pca$x[,2]-0.4, label = colnames(x), cex = 1.6) +
               ggplot2::ggtitle(paste0("PCA plot using ", colnames(intgroup.df)),
                                subtitle = paste0(colnames(d)[1], " vs. ", colnames(d)[2]))
      }
      return(ggp)
    }

    ##lapply to make all PC > pc_limit included
    pcv_u <- percentVar[percentVar > pc_limit/100]
    ggps <- lapply(2:length(pcv_u), function(f){

        d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, f], intgroup.df, names = colnames(x))
        colnames(d)[colnames(d) == "PC2"] <- paste0("PC", f)
        ggp <- pca_plot(d, group, intgroup.df, pcv = percentVar[c(1,f)])
        return(ggp)

    })


    ##scree plot showing contributions of PCs
    pcv <- data.frame(PC = colnames(pca$x), percent_variance = percentVar)
    pcv <- pcv[order(pcv$percent_variance, decreasing = TRUE),]
    levels(pcv$PC) <- pcv$PC
    ggs <- ggplot2::ggplot(data = pcv, ggplot2::aes(x = PC, y = percent_variance )) +
           ggplot2::geom_col() +
           ggplot2::ggtitle("Proprotion of Variances of Principle Components") +
           ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

    ##loading top 5, bottom 5 from 10 top PCs
    top_pc <- pcv[1:10,][!is.na(pcv[1:10, "PC"]),]
    top_pc$percent_variance <- paste0(round(top_pc$percent_variance, digits = 3)*100, "%")
    pca_rot <- pca$rotation[,top_pc$PC]
    fives_list <- lapply(colnames(pca_rot), function(f){
      fc <- c(tail(sort(pca_rot[,f]),5), head(sort(pca_rot[,f]),5))
      nfc <- dplyr::filter(.data = anno_tb, ensembl_gene_id %in% names(fc))
      names(fc) <- unlist(lapply(names(fc), function(ff){
        dist_nfc <- dplyr::filter(.data = nfc, ensembl_gene_id %in% ff) %>%
        dplyr::select(external_gene_name) %>% dplyr::distinct() %>% unlist()
        paste(dist_nfc, collapse = ",")
      }))
      return(data.frame(gene = names(fc), loading = fc, PC = f))
    })

    loadings_plot <- do.call(rbind,fives_list) %>%
                     dplyr::left_join(top_pc) %>%
                     dplyr::mutate(PC = paste0(PC, "\n(", percent_variance, ")"))

    ggl <- ggplot2::ggplot(data = loadings_plot,
                           ggplot2::aes_string(x = "PC", y = "loading", label = "gene")) +
           ggplot2::geom_point(ggplot2::aes_string(colour = "loading"),
                               size = 3) +
           ggplot2::scale_shape_discrete(solid = T) +
           ggrepel::geom_text_repel(colour = "black",
                                    size = 2,
                                    fontface = "bold") +
           ggplot2::ggtitle("Loadings plot, top 10 PCs by variance, top/bottom 5 genes") +
           ggplot2::xlab("PC (% variance accounted for)")
    return(list(ggps, ggs, ggl))
}

#' Parse information from STAR run to find used GTF file
#' @importFrom magrittr '%>%'
#' @export

genesGTF <- function(output_dir = OUTDR){

  ##get genes.gtf input
  genesgtf_file <- system("find work/stage -name genes.gtf | grep -v STAR", intern=T)
  system(paste0("cat ", genesgtf_file, " | perl -ane 'if(($F[1] eq \"ensembl\") && ($F[2] eq \"CDS\")){print $_;}' > ", output_dir, "/gene.ensembl.CDS.gtf"))
  cds_gtf_idnm <- read_table2(paste0(output_dir, "/gene.ensembl.CDS.gtf")) %>%
             dplyr::select(14,16)

  ##parse CDS from genes.gtf
  splitFun <- function(x){unlist(lapply(x, function(f){strsplit(f, '\\"')[[1]][2]}))}
  cds_gtf_nm <- tibble(ensembl_gene_id = splitFun(unlist(cds_gtf_idnm[,1])),
                       external_gene_name = splitFun(unlist(cds_gtf_idnm[,2]))) %>%
                distinct() %>%
                arrange(ensembl_gene_id)

  return(cds_gtf_nm)
}

#' Take two inputs of datasets in biomaRt::useMart and return the tibble of the matched orthologs
#'
#' @param human_genome is the genome against which we want orthologs (set as human but changeble with unknown results)
#' @param org_prefix string indicating the name of the organism which is suffixed with '_gene_ensembl' for biomaRt
#' @param ens_version string indicating version of ENSEMBL to use
#' @return ensid2gene2orth a tibble of orthologs in ensembl - gene name - ortholog
#' @importFrom magrittr '%>%'
#' @export

biomaRt_anno_orth <- function(human_genome = "hsapiens_gene_ensembl", org_prefix = NULL, ens_version = NULL){

  ##need a way to access versions
  ##if you know you need a specific version for acces to a specific genome version
  ##then specify version; otehrwise defaults to latest by below
  if(is.null(ens_version)){
    ens_version <- rev(strsplit(biomaRt::listMarts()[1,2], " ")[[1]])[1]
  }

  ##specifies the HOST connect (URL)
  HOST <- as_tibble(listEnsemblArchives()) %>%
          dplyr::filter(version %in% ens_version) %>%
          dplyr::select(url) %>% unlist()

  ##access the two genomes
  mart1 <- biomaRt::useMart(biomart = "ensembl", dataset = human_genome, host = HOST)
  mart2 <- biomaRt::useMart(biomart = "ensembl", dataset = paste0(org_prefix, "_gene_ensembl"), host = HOST)

  ##first genome genes with homologs for second genome
  ensid2gene1 <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name", paste0(org_prefix, "_homolog_ensembl_gene")), mart = mart1))

  ##second genome, with renaming to allow join with above
  renames <- c(paste0(org_prefix,"_homolog_ensembl_gene"),
               paste0(org_prefix,"_homolog_external_name"))
  ensid2gene2 <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = mart2))
  ensid2gene2 <- dplyr::select(.data = ensid2gene2, !!renames[1] := ensembl_gene_id, !!renames[2] := external_gene_name)

  ##join to get human and ortholog IDs
  ensid2gene2orth <- left_join(ensid2gene1, ensid2gene2, by=renames[1]) %>%
                     na.omit()
  return(ensid2gene2orth)
}
