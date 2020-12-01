#' Function to run ssGSEA using the GSVA package and return PCA of results
#'
#' @param cont_pways list of per-contrast pathway genesets to use in ssGSEA
#' @param log2tpm_tb tibble of log2tpm values per sample across genes in gene lists
#' @param msigdb_cat one of 'c("H", paste0("C", c(1:7)))', see: gsea-msigdb.org/gsea/msigdb/collections.jsp; if H, loads Process_Category data from Liberzon 2015 paper to colour PCA
#' @param hallmark_tb hallmark data with Process_Category column to colour by process
#' @param output_dir path to where output goes
#' @param contrast string used to prefix output naming the contrast being investigated
#' @return ssgsea results table of each geneset and sample
#' @export

ssgsea_pca <- function(cont_pways, log2tpm_mat, msigdb_cat = "H", hallmark_tb = NULL, output_dir, contrast) {

  print(paste0("Working on: ", contrast))

  ##output dir
  out_dir <- paste0(output_dir, "/ssgsea")
  dir.create(paste0(out_dir, "/plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(out_dir, "/data"), recursive = TRUE, showWarnings = FALSE)

  ##use res as significant set also
  if(msigdb_cat == "H"){
      # hallmarks <- system.file("extdata",
      #                          "Liberzon_2015_Table1_Hallmark_categories.csv",
      #                          package = "RNAseqR")
      hallmarks <- readr::read_csv("https://raw.githubusercontent.com/brucemoran/RNAseqR/master/extdata/Liberzon_2015_Table1_Hallmark_categories.csv")
      rns <- "Hallmark_Name"
  } else {
    print(paste0("No colouring will be added to PCA baecause you are using: ", msigdb_cat))
    rns <- msigdb_cat
  }

  ##GSVA ssGSEA
  ssgsea_res <- GSVA::gsva(log2tpm_mat,
                           cont_pways,
                           method='ssgsea',
                           min.sz=0,
                           max.sz=1000,
                           ssgsea.norm=T)

  ggp <- ssgsea_rotationPCA(x = ssgsea_res,
                            hallmark_tb = hallmarks,
                            contrast = contrast,
                            returnData = FALSE)

  readr::write_tsv(tibble::as_tibble(ssgsea_res, rownames = rns), file = paste0(out_dir, "/data/", contrast , ".ssgsea_res.tsv"))
  ggplot2::ggsave(ggp[[1]], file = paste0(out_dir, "/plots/", contrast , ".ssgsea_rotationPCA.pdf"))
  return(list(ssgsea_res, ggp))
}

#' Plotting PCA function
#' @param ssgsea output table
#' @param hallmark_tb hallmark data with Process_Category column to colour by process
#' @param contrast tag to title plot
#' @param returnData flag to allow data to be returned
#' @return ggp ggplot2 object for printing
#' @export

ssgsea_rotationPCA <- function(x, hallmark_tb = NULL, contrast, returnData = FALSE) {

    pca <- prcomp(t(x))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)

    d <- data.frame(PC1 = pca$rotation[, 1], PC2 = pca$rotation[, 2], Hallmark_Name = rownames( pca$rotation))

    if(!is.null(hallmark_tb)){
      hallmark_tb$Hallmark_Name <- paste0("HALLMARK_", hallmark_tb$Hallmark_Name)
      d <- dplyr::left_join(d, hallmark_tb, by = "Hallmark_Name")
      d$Process_Category[is.na(d$Process_Category)] <- "other"
    }
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }

    if("Process_Category" %in% colnames(d)){
      ggp <- ggplot2::ggplot(data = d, ggplot2::aes(x = PC1, y = PC2, group = Hallmark_Name, colour = Process_Category)) +
             ggplot2::labs(title = contrast,
                           subtitle = "PCA plot using MsigDB Hallmark with Process_Category",
                           x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
             ggplot2::geom_segment(data = d,
                                   ggplot2::aes(x = 0, y = 0, xend = (PC1*3),
                                   yend = (PC2*3), color = Process_Category),
                                   arrow = ggplot2::arrow(length = ggplot2::unit(1/2, "picas"),
                                                         ends="last", type = "closed"), size = 1) +
             ggrepel::geom_text_repel(data = d,
                                ggplot2::aes(x = (PC1*3), y = (PC2*3), label = gsub("HALLMARK_", "", Hallmark_Name), colour = Process_Category),
                                size = 2, fontface = "bold")
    } else {
      stop("Unfinished, suggest using Hallmarks...")
      # ggp <- ggplot2::ggplot(data = d, ggplot2::aes(x = PC1, y = PC2, group = intgroup, shape = intgroup, colour = intgroup)) +
      #        ggplot2::geom_point(size = 3) +
      #        ggplot2::scale_shape_discrete(solid = T) +
      #        ggplot2::xlab(paste0("PC1: ", round(percentVar[1] *  100), "% variance")) +
      #        ggplot2::ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      #        ggplot2::annotate("text",x=pca$x[,1], y = pca$x[,2]-0.4, label = colnames(x), cex = 1.6) +
      #        ggplot2::ggtitle(paste0("PCA plot using ", intgroup))
    }
    return(list(ggp, pca))
}
