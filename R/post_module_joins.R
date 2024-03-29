#! R

#' Master list, join to give found-in-2 and all-3 tables
#'
#' @param input_dir path where DESeq2, limma, edgeR RDS results files are held
#' @param tpm_tb tpm tibble for adding to master_ist_sig ([[x]][[2]])
#' @return list object of 'found-in-2' and 'all-3' results, including unique cols from each module
#' @importFrom magrittr '%>%'
#' @export

master_parse_join <- function(input_dir, tpm_tb = NULL){

  print("Running: master_parse_join()")

  ##read RDS results
  rdss <- dir(input_dir, pattern = "_res_list.rds", recursive = TRUE, full.names = TRUE)
  master_names <- unlist(lapply(rdss, function(f){
    filn <- rev(strsplit(f, "/")[[1]])[1]
    gsub("_res_list", "", strsplit(filn, "\\.")[[1]][2])
  }))
  master_list <- lapply(rdss, function(f){
    readRDS(f)
  })
  names(master_list) <- master_names

  ##check if multiple annotations arise and fix by aggregating if so
  ##esp for homology tables (i.e. where !genome_prefix == hsapiens)
  master_list_a <- lapply(master_list, function(ff){
    mlist <- lapply(ff, function(f){
      extg <- f$external_gene_name
      ensg <- f$ensembl_gene_id
      if(length(table(extg[table(extg)>1])) > 1 | length(table(ensg[table(ensg)>1])) > 1){
        print("Found multiple annotations for identifiers, aggregating results")

        colns <- grep("_gene", colnames(f), value = TRUE)

        if(length(colns) == 4){
          ##have 4 cols, first 2 to group_by, second 2 to aggregate
          f <- group_agg_multi(f, colns)
        }

        ##then aggregate ENS IDs so only single external_gene_name is found
        group_agg_two(f, colns)
      } else {
        return(f)
      }
    })
    return(mlist)
  })

  ##list within names of significant results (p < 0.01)
  master_list_sig <- lapply(master_list_a, function(f){

    ff <- dplyr:::filter(.data = f[[1]], padj < 0.01)

    ##add TPM to genes
    if(!is.null(tpm_tb)){

      ff <- dplyr::left_join(ff, tpm_tb) %>%
            dplyr::arrange(external_gene_name)
    }
    return(ff)
  })
  names(master_list_sig) <- names(master_list_a)

  save(master_list_sig, file = paste0(input_dir, "/RData/", tag, ".master_list_sig.RData"))

  return(master_list_a)
}

#' Aggregate columns 3, 4 based on grouping on columns 1, 2
#' NB that table can have multiple annotations in RNAseqon based on contrasts
#' this function removes this 'per-caller'
#'
#' @param f table object with colnames specified in colns
#' @param colns colnames of f on which to aggregate
#' @param pattern string to grep for colnames on which to operate
#' @return agg_f aggregated gene ID, names table
#' @importFrom magrittr '%>%'
#' @export

group_agg_multi <- function(f, colns = NULL, pattern = NULL){

  if(is.null(colns) && is.null(pattern)){
    stop("Please specify one of 'colns' or 'patterns'")
  }

  if(is.null(colns)){
    colns <- grep(pattern, colnames(f), value = TRUE)
  }

  agg_f <- f %>% dplyr:::rowwise() %>%
                 dplyr::distinct(dplyr::across(-!!as.symbol(colns[3])), .keep_all = TRUE) %>%
                 dplyr::ungroup()
  return(agg_f)
}

#' Aggregate column 1 based on grouping by column 2
#'
#' @param f table object with colnames specified in colns
#' @param colns colnames of f on which to aggregate, group in that order!
#' @param pattern string to grep for colnames on which to operate
#' @return agg_f aggregated gene ID, names table
#' @importFrom magrittr '%>%'
#' @export

group_agg_two <- function(f, colns = NULL, pattern = NULL){

  if(is.null(colns) && is.null(pattern)){
    stop("Please specify one of 'colns' or 'patterns'")
  }

  if(is.null(colns)){
    colns <- grep(pattern, colnames(f), value = TRUE)
  }

  agg_f <- f %>% dplyr:::rowwise() %>%
                 dplyr::distinct(dplyr::across(-!!as.symbol(colns[1])), .keep_all = TRUE) %>%
                 dplyr::ungroup()
  return(agg_f)
}

#' Found-in-2
#'
#' @param master_list list of 3 DE results
#' @param padj significance threshold below which to determine significance (adjusted by FDR)
#' @return fitwo_list list of tibble of results including unique cols from each module
#' @importFrom magrittr '%>%'
#' @export

found_in_two <- function(master_list, padj = 0.01){

  ##pairwise joins
  ##operate over multiple contrasts
  combns <- apply(t(combn(m = 2, names(master_list))), 1, function(f){
    paste(f, collapse = "-")
  })
  fitwo_list <- apply(t(combn(m = 2, names(master_list))), 1, function(f){
    per_contrast_list <- lapply(names(master_list[[1]]), function(nm){
      print(paste0("Working on found-in-two, contrast: ", nm))
      ml1 <- master_list

      ##rename colnames with method ID
      cn1 <- colnames(ml1[[f[1]]][[nm]])
      cn2 <- colnames(ml1[[f[2]]][[nm]])
      cn1i <- grep("_gene", grep("tpm", cn1, invert = TRUE, value = TRUE), invert = TRUE)
      cn2i <- grep("_gene", grep("tpm", cn2, invert = TRUE, value = TRUE), invert = TRUE)
      colnames(ml1[[f[1]]][[nm]])[cn1i] <- paste0(f[1], "_", cn1[cn1i])
      colnames(ml1[[f[2]]][[nm]])[cn2i] <- paste0(f[2], "_", cn2[cn2i])
      f1 <- paste0(f[1], "_padj")
      f2 <- paste0(f[2], "_padj")

      base::suppressMessages(dplyr::left_join(ml1[[f[1]]][[nm]], ml1[[f[2]]][[nm]])) %>%
      dplyr::filter(!!as.symbol(f1) < !!padj & !!as.symbol(f2) < !!padj)
    })

    names(per_contrast_list) <- names(master_list[[1]])
    return(per_contrast_list)
  })
  names(fitwo_list) <- combns
  return(fitwo_list)
}

#' Found-in-3
#'
#' @param master_list list of 3 DE results
#' @param padj significance threshold below which to determine significance (adjusted by FDR)
#' @return per_contrast_list list of tibble of results including unique cols from each module
#' @importFrom magrittr '%>%'
#' @export

found_in_three <- function(master_list, padj = 0.01){

  ##as per find-in-two
  per_contrast_list <- lapply(names(master_list[[1]]), function(nm){
    ml1 <- master_list
    f <- names(ml1)
    print(paste0("Working on found-in-three, contrast: ", nm))

    ##rename colnames with method ID
    cn1 <- colnames(ml1[[f[1]]][[nm]])
    cn2 <- colnames(ml1[[f[2]]][[nm]])
    cn3 <- colnames(ml1[[f[3]]][[nm]])
    cn1i <- grep("_gene", grep("tpm", cn1, invert = TRUE, value = TRUE), invert = TRUE)
    cn2i <- grep("_gene", grep("tpm", cn2, invert = TRUE, value = TRUE), invert = TRUE)
    cn3i <- grep("_gene", grep("tpm", cn3, invert = TRUE, value = TRUE), invert = TRUE)
    colnames(ml1[[f[1]]][[nm]])[cn1i] <- paste0(f[1], "_", cn1[cn1i])
    colnames(ml1[[f[2]]][[nm]])[cn2i] <- paste0(f[2], "_", cn2[cn2i])
    colnames(ml1[[f[3]]][[nm]])[cn3i] <- paste0(f[3], "_", cn3[cn3i])
    f1 <- paste0(f[1], "_padj")
    f2 <- paste0(f[2], "_padj")
    f3 <- paste0(f[3], "_padj")

    mlj <- base::suppressMessages(dplyr::left_join(ml1[[f[1]]][[nm]], ml1[[f[2]]][[nm]]))
    base::suppressMessages(dplyr::left_join(mlj, ml1[[f[3]]][[nm]])) %>%
    dplyr::filter(!!as.symbol(f1) < !!padj & !!as.symbol(f2) < !!padj & !!as.symbol(f3) < !!padj)

  })
  names(per_contrast_list) <- names(master_list[[1]])
  return(per_contrast_list)
}

#' Make venn diagram
#'
#' @param master_list list of 3 DE results
#' @param padj significance threshold below which to determine significance (adjusted by FDR)
#' @param output_dir path to where output goes
#' @param tag string used to prefix output
#' @return none, venn diagram of 3 DE callers printed
#' @importFrom magrittr '%>%'
#' @export

venn_3 <- function(master_list, tag, output_dir, padj = 0.01){

  print("Running: venn_3()")

  ##run over master_list to get number of overlapped genes
  print("Working on found-in-two")
  fitwo <- suppressMessages(RNAseqon::found_in_two(master_list))
  print("Working on found-in-three")
  fithree <- suppressMessages(RNAseqon::found_in_three(master_list))

  ##each master_list for contrast
  unc <- unique(t(combn(x = c(0,0,0,1,1,1,0,0,0,1,1,1), m = 3)))
  unc <- unc[order(rowSums(unc)),]
  unct <- unc == 1

  ##create list of ensembl_gene_ids per set defined by unc
  lapply(names(master_list[[1]]), function(contrast){
    print(paste0("Working on contrast: ", contrast))
    compf_list <- lapply(seq(1:dim(unct)[1]), function(f){
      compf <- names(master_list)[unct[f,]]
      print(compf)
      if(length(compf) == 0){
        return()
      } else if(length(compf) == 1){
        dat <- master_list[[compf]][[contrast]]
        if(dim(dat)[1] == 0){
          return()
        } else {
          dat1 <- dplyr::filter(.data = dat, padj < !!padj)
          dat2 <- dplyr::select(.data = dat1, ensembl_gene_id)
          as.vector(unlist(dat2))
        }
      } else if(length(compf) == 2){
        ##aready padjusted
        dat <- fitwo[[paste(compf, collapse = "-")]][[contrast]]
        if(dim(dat)[1] == 0){
          return()
        } else {
          dat1 <- dplyr::select(.data = dat, ensembl_gene_id)
          as.vector(unlist(dat1))
        }
      } else {
        dat <- fithree[[contrast]]
        if(dim(dat)[1] == 0){
          return()
        } else {
          dat1 <- dplyr::select(.data = dat, ensembl_gene_id)
          as.vector(unlist(dat1))
        }
      }
    })

    ##compare those lists
    vdf <- cbind(unc, unlist(lapply(compf_list, length)))
    colnames(vdf) <- c(names(master_list), "Counts")
    dir.create(paste0(output_dir, "/venn_3"), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(as.data.frame(vdf), file = paste0(output_dir, "/venn_3/", tag, ".", contrast, ".venn_3.csv"))

    ##write venn
    trip_venn <- VennDiagram::draw.triple.venn(area1 = vdf[4,4],
                                           area2 = vdf[3,4],
                                           area3 = vdf[2,4],
                                           n12 = vdf[6,4],
                                           n13 = vdf[7,4],
                                           n23 = vdf[5,4],
                                           n123 = vdf[8,4],
                                           category = colnames(vdf)[1:3],
                                           fill = c("dodgerblue", "orange", "forestgreen"),
                                           cat.cex = 2, cex = 2)
    pdf(paste0(output_dir, "/venn_3/", tag, ".", contrast, ".venn_3.pdf"))
      title <- grid::textGrob(contrast,
                              x = 0.2, y = 0.01,
                              vjust = 0,
                              gp = grid::gpar(fontsize = 12))
      gt <- grid::gTree(children = grid::gList(title, trip_venn))
      grid::grid.draw(gt)
    dev.off()
  })
}
