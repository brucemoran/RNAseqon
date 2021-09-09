#' Function to run fGSEA using the mSIGDB database from Broad
#'
#' @param gtf path/to/GTF.gtf
#' @param paired boolean of whether BAM data is paired
#' @param stranded string '0', '1', '2', indicating 'unstranded', 'forward', reverse' respectively
#' @param threads total threads used to run
#'
#' @return RData object of dupRadar outputs in list
#' @export

dupradar_run <- function(gtf, paired, stranded = '0', threads) {

  #define BAMs
  in_bams <- dir(pattern = "bam")
  dup_out <- as.list(in_bams)

  dir.create("dupradar")

  for (x in 1:length(in_bams)){

    dup_out <- lapply(seq_along(in_bams), function(x){
      dupo <- dupRadar::analyzeDuprates(in_bams[x], gtf, stranded, paired, threads)
      pdf(paste0("dupradar/", in_bams[[x]], ".duprate_exp_densplot.pdf"))
        dupRadar::duprateExpDensPlot(DupMat = dupo)
      dev.off()

      pdf(paste0("dupradar/", in_bams[[x]], ".duprate_exp_boxplot.pdf"))
        dupRadar::duprateExpBoxplot(DupMat = dupo)
      dev.off()
      return(dupo)
    })
    names(dup_out) <- in_bams
    save(dup_out, file="dupradar/dupradar.analyzeDuprates.RData")
  }
}
