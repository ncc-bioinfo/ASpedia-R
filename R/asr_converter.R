#rMATS result file, AS type, P-value cutoff, DPSI cutoff
#SUPPA result file, gtf file, ioe file, gene model, P-value cutoff, DPSI cutoff
#spliceR result file, gene model, genome version, pvalue cutoff

#' Generate ASpedia-R input format from DAS analysis tools result.
#' 
#' result of DAS analysis tools convert to ASpedia-R input format.
#'
#' @param das.analysis.result
#' name of DAS analysis tools result file.
#' @param program
#' name of DAS analysis tool. one of rMATS, SUPPA, or spliceR
#' @param pvalue.cutoff
#' value of pvalue cutoff. default value is 0.05
#' @param dpsi.cutoff
#' value of dPSI cutoff. default value is 0.1
#' @param gene.model
#' gene model of reference. One of Refseq, Ensembl, or GENCODE. (spliceR only)
#' @param genome.version
#' genome version of reference. One of hg18, GRCh19, or GRCh38. (spliceR only)
#' @param as.type
#' AS event type. One of A3SS, A5SS, SE, MXE, or RI. (rMATS only)
#' @param gtf.file.name
#' a GTF format file of reference. (SUPPA only)
#' @param ioe.file.name
#' name of ioe file generated by SUPPA generateEvents command. (SUPPA)
#'
#' @return converter.result
#' @export
#'
#' @examples 
#' ## rMATS
#' das.file.name <- system.file("extdata" "rMATS_test.txt", package="ASpediaR")
#' rmats.convert.result <- asr_converter(das.file.name, program="rMATS", as.type="SE")
#' 
#' ##SUPPA
#' das.file.name <- system.file("extdata", "SUPPA_test.txt", package="ASpediaR")
#' gtf.file.name <- system.file("extdata", "test_gtf.gtf", package="ASpediaR")
#' ioe.file.name <- system.file("extdata", "SUPPA_test.ioe", package="ASpediaR")
#' suppa.convert.result <- asr_converter(das.file.name, program="SUPPA",
#'                                         gtf.file=gtf.file.name, ioe.file=ioe.file.name)
#' 
#' ##spliceR
#' das.file.name <- system.file("extdata", "spliceR_test.txt", package="ASpediaR")
#' splilcer.convert.result <- asr_converter(das.file.name, program="spliceR",
#'                                           gene.model="Ensembl", genome.version="GRCh38")

asr_converter <- function(das.analysis.result="", program="", pvalue.cutoff=0.05, dpsi.cutoff=0.1, gene.model="Ensembl", genome.version="GRCh38", as.type="", gtf.file.name="", ioe.file.name="") {
  if(tolower(program) == "rmats") {
    converter.result <- rMATS.converter(das.analysis.result, as.type, pvalue.cutoff, dpsi.cutoff)

    return(converter.result)
  } else if(tolower(program) == "suppa") {
    converter.result <- SUPPA.converter(das.analysis.result, pvalue.cutoff, dpsi.cutoff, gtf.file.name, ioe.file.name)

    return(converter.result)
  } else if(tolower(program) == "splicer") {
    converter.result <- spliceR.converter(das.analysis.result, pvalue.cutoff, gene.model, genome.version)

    return(converter.result)
  } else {
    message("ERROR : Input program name is wrong. Please check program option. We currently support rMATS, SUPPA, or spliceR.")
    return()
  }
}
