#rMATS result file, AS type, P-value cutoff, DPSI cutoff
#' Generate ASpedia-R input format from rMATS result.
#' 
#' rMSTS result convert to ASpedia-R input format.
#'
#' @param rMATS.result
#' name of rMATS result file.
#' @param as.type
#' AS event type. One of A3SS, A5SS, SE, MXE, or RI.
#' @param pvalue.cutoff
#' value of pvalue cutoff. default value is 0.05
#' @param dpsi.cutoff
#' value of dPSI cutoff. default value is 0.1
#'
#' @return converting.result
#' 
#' @export
#'
#' @examples
#' rmats.result.file.name <- system.file(“extdata”, “rMATS_test.txt”, package=“ASpediaR”)
#' rmats.convert.result <- rMATS_converter(rmats.result.file.name, program=“rMATS”,
#'                                             as.type=“SE”)

rMATS.converter <- function(rMATS.result="", as.type="", pvalue.cutoff=0.05, dpsi.cutoff=0.1)
{
  if(rMATS.result == "" || file.exists(rMATS.result) == FALSE) {
    message(paste0("*** ERROR MESSAGE: No such input file. ", rMATS.result))
    return()
  }

  rMATS.input <- read.table(rMATS.result, header=TRUE, sep="\t", stringsAsFactors=FALSE)

  if(nrow(rMATS.input) == 0)
  {
    message("*** ERROR MESSAGE: Input file is empty. Please check input file.")
    return()
  }
  
  if(as.type == "" || (tolower(as.type) != "a3ss" && tolower(as.type) != "a5ss" && tolower(as.type) != "se" && tolower(as.type) != "mxe" && tolower(as.type) != "ri")) {
    message("*** ERROR MESSAGE: Input AS type is wrong. Please check AS.type option. We only support A5SS, A3SS, SE, MXE, RI.")
    return()
  }
  
  loaded.packages <- tolower((.packages()))
  
  if(("stringr" %in% loaded.packages) == FALSE) {
    library(stringr)
  }
  
  #"chr", "strand", "as_type", "upstream_exon", "splicing exon", "downstream exon", "gene_id", "gene_name", "P-value", "dPSI"
  if(as.type == "A3SS") {
    tmp.line <- rMATS.input[1, ]

    if(length(tmp.line) != 23) {
      message("*** ERROR MESSAGE: Input file format is wrong. Please check input file.")
      return()
    }

    rMATS.input.filtering <- rMATS.input[abs(rMATS.input$IncLevelDifference) > dpsi.cutoff & rMATS.input$FDR < pvalue.cutoff, ]
    
    if(str_length(tmp.line$chr) <= 2) {
      rMATS.input.filtering$chr <- paste0("chr", rMATS.input.filtering$chr)
    }

    rMATS.plus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "+", ]
    rMATS.minus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "-", ]

    plus.as.id <- paste0(rMATS.plus.strand$chr, ":", (as.numeric(rMATS.plus.strand$flankingES) + 1), ":", rMATS.plus.strand$flankingEE, ":", (as.numeric(rMATS.plus.strand$longExonStart_0base) + 1), ":", (as.numeric(rMATS.plus.strand$shortES) + 1), ":", rMATS.plus.strand$shortEE)
    plus.upstream.exon <- paste0((as.numeric(rMATS.plus.strand$flankingES) + 1), "-", rMATS.plus.strand$flankingEE)
    plus.splicing.exon <- paste0((as.numeric(rMATS.plus.strand$longExonStart_0base) + 1), "-", rMATS.plus.strand$longExonEnd)
    plus.downstream.exon <- paste0((as.numeric(rMATS.plus.strand$shortES) + 1), "-", rMATS.plus.strand$shortEE)

    plus.strand.result <- data.frame(chr=rMATS.plus.strand$chr, strand=rMATS.plus.strand$strand, as_type=as.type, upstream_exon=plus.upstream.exon, splicing_exon=plus.splicing.exon, downstream_exon=plus.downstream.exon, gene_name=rMATS.plus.strand$geneSymbol, pvalue=rMATS.plus.strand$FDR, dPSI=rMATS.plus.strand$IncLevelDifference, as_id=plus.as.id, gene_id=rMATS.plus.strand$GeneID)

    minus.as.id <- paste0(rMATS.minus.strand$chr, ":", rMATS.minus.strand$flankingEE, ":", (as.numeric(rMATS.minus.strand$flankingES) + 1), ":", rMATS.minus.strand$longExonEnd, ":", rMATS.minus.strand$shortEE, ":", (as.numeric(rMATS.minus.strand$longExonStart_0base) + 1))
    minus.upstream.exon <- paste0(rMATS.minus.strand$flankingEE, "-", (as.numeric(rMATS.minus.strand$flankingES) + 1))
    minus.splicing.exon <- paste0(rMATS.minus.strand$longExonEnd, "-", (as.numeric(rMATS.minus.strand$longExonStart_0base) + 1))
    minus.downstream.exon <- paste0(rMATS.minus.strand$shortEE, "-", (as.numeric(rMATS.minus.strand$shortES) + 1))

    minus.strand.result <- data.frame(chr=rMATS.minus.strand$chr, strand=rMATS.minus.strand$strand, as_type=as.type, upstream_exon=minus.upstream.exon, splicing_exon=minus.splicing.exon, downstream_exon=minus.downstream.exon, gene_name=rMATS.minus.strand$geneSymbol, pvalue=rMATS.minus.strand$FDR, dPSI=rMATS.minus.strand$IncLevelDifference, as_id=minus.as.id, gene_id=rMATS.minus.strand$GeneID)

    converting.result <- rbind(plus.strand.result, minus.strand.result)
  } else if(as.type == "A5SS") {
    tmp.line <- rMATS.input[1, ]

    if(length(tmp.line) != 23) {
      message("*** ERROR MESSAGE: Input file format is wrong. Please check input file.")
      return()
    }

    rMATS.input.filtering <- rMATS.input[abs(rMATS.input$IncLevelDifference) > dpsi.cutoff & rMATS.input$FDR < pvalue.cutoff, ]
    
    if(str_length(tmp.line$chr) <= 2) {
      rMATS.input.filtering$chr <- paste0("chr", rMATS.input.filtering$chr)
    }

    rMATS.plus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "+", ]
    rMATS.minus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "-", ]

    plus.as.id <- paste0(rMATS.plus.strand$chr, ":", (as.numeric(rMATS.plus.strand$longExonStart_0base) + 1), ":", rMATS.plus.strand$shortEE, ":", rMATS.plus.strand$longExonEnd, ":", (as.numeric(rMATS.plus.strand$flankingES) + 1), ":", rMATS.plus.strand$flankingEE)
    plus.upstream.exon <- paste0((as.numeric(rMATS.plus.strand$shortES) + 1), "-", rMATS.plus.strand$shortEE)
    plus.splicing.exon <- paste0((as.numeric(rMATS.plus.strand$longExonStart_0base) + 1), "-", rMATS.plus.strand$longExonEnd)
    plus.downstream.exon <- paste0((as.numeric(rMATS.plus.strand$flankingES) + 1), "-", rMATS.plus.strand$flankingEE)

    plus.strand.result <- data.frame(chr=rMATS.plus.strand$chr, strand=rMATS.plus.strand$strand, as_type=as.type, upstream_exon=plus.upstream.exon, splicing_exon=plus.splicing.exon, downstream_exon=plus.downstream.exon, gene_name=rMATS.plus.strand$geneSymbol, pvalue=rMATS.plus.strand$FDR, dPSI=rMATS.plus.strand$IncLevelDifference, as_id=plus.as.id, gene_id=rMATS.plus.strand$GeneID)

    minus.as.id <- paste0(rMATS.minus.strand$chr, ":", rMATS.minus.strand$longExonEnd, ":", (as.numeric(rMATS.minus.strand$shortES) + 1), ":", (as.numeric(rMATS.minus.strand$longExonStart_0base) + 1), ":", rMATS.minus.strand$flankingEE, ":", (as.numeric(rMATS.minus.strand$flankingES) + 1))
    minus.upstream.exon <- paste0(rMATS.minus.strand$shortEE, "-", (as.numeric(rMATS.minus.strand$shortES) + 1))
    minus.splicing.exon <- paste0(rMATS.minus.strand$longExonEnd, "-", (as.numeric(rMATS.minus.strand$longExonStart_0base) + 1))
    minus.downstream.exon <- paste0(rMATS.minus.strand$flankingEE, "-", (as.numeric(rMATS.minus.strand$flankingES) + 1))

    minus.strand.result <- data.frame(chr=rMATS.minus.strand$chr, strand=rMATS.minus.strand$strand, as_type=as.type, upstream_exon=minus.upstream.exon, splicing_exon=minus.splicing.exon, downstream_exon=minus.downstream.exon, gene_name=rMATS.minus.strand$geneSymbol, pvalue=rMATS.minus.strand$FDR, dPSI=rMATS.minus.strand$IncLevelDifference, as_id=minus.as.id, gene_id=rMATS.minus.strand$GeneID)

    converting.result <- rbind(plus.strand.result, minus.strand.result)
  } else if(as.type == "SE") {
    tmp.line <- rMATS.input[1, ]

    if(length(tmp.line) != 23) {
      message("*** ERROR MESSAGE: Input file format is wrong. Please check input file.")
      return()
    }

    rMATS.input.filtering <- rMATS.input[abs(rMATS.input$IncLevelDifference) > dpsi.cutoff & rMATS.input$FDR < pvalue.cutoff, ]
    
    if(str_length(tmp.line$chr) <= 2) {
      rMATS.input.filtering$chr <- paste0("chr", rMATS.input.filtering$chr)
    }

    rMATS.plus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "+", ]
    rMATS.minus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "-", ]

    plus.as.id <- paste0(rMATS.plus.strand$chr, ":", (as.numeric(rMATS.plus.strand$upstreamES) + 1), ":", rMATS.plus.strand$upstreamEE, ":", (as.numeric(rMATS.plus.strand$exonStart_0base) + 1), ":", rMATS.plus.strand$exonEnd, ":", (as.numeric(rMATS.plus.strand$downstreamES) + 1), ":", rMATS.plus.strand$downstreamEE)
    plus.upstream.exon <- paste0((as.numeric(rMATS.plus.strand$upstreamES) + 1), "-", rMATS.plus.strand$upstreamEE)
    plus.splicing.exon <- paste0((as.numeric(rMATS.plus.strand$exonStart_0base) + 1), "-", rMATS.plus.strand$exonEnd)
    plus.downstream.exon <- paste0((as.numeric(rMATS.plus.strand$downstreamES) + 1), "-", rMATS.plus.strand$downstreamEE)

    plus.strand.result <- data.frame(chr=rMATS.plus.strand$chr, strand=rMATS.plus.strand$strand, as_type=as.type, upstream_exon=plus.upstream.exon, splicing_exon=plus.splicing.exon, downstream_exon=plus.downstream.exon, gene_name=rMATS.plus.strand$geneSymbol, pvalue=rMATS.plus.strand$FDR, dPSI=rMATS.plus.strand$IncLevelDifference, as_id=plus.as.id, gene_id=rMATS.plus.strand$GeneID)

    minus.as.id <- paste0(rMATS.minus.strand$chr, ":", rMATS.minus.strand$downstreamEE, ":", (as.numeric(rMATS.minus.strand$downstreamES) + 1), ":", rMATS.minus.strand$exonEnd, ":", (as.numeric(rMATS.minus.strand$exonStart_0base) + 1), ":", rMATS.minus.strand$upstreamEE, ":", (as.numeric(rMATS.minus.strand$upstreamES) + 1))
    minus.upstream.exon <- paste0(rMATS.minus.strand$downstreamEE, "-", (as.numeric(rMATS.minus.strand$downstreamES) + 1))
    minus.splicing.exon <- paste0(rMATS.minus.strand$exonEnd, "-", (as.numeric(rMATS.minus.strand$exonStart_0base) + 1))
    minus.downstream.exon <- paste0(rMATS.minus.strand$upstreamEE, "-", (as.numeric(rMATS.minus.strand$upstreamES) + 1))

    minus.strand.result <- data.frame(chr=rMATS.minus.strand$chr, strand=rMATS.minus.strand$strand, as_type=as.type, upstream_exon=minus.upstream.exon, splicing_exon=minus.splicing.exon, downstream_exon=minus.downstream.exon, gene_name=rMATS.minus.strand$geneSymbol, pvalue=rMATS.minus.strand$FDR, dPSI=rMATS.minus.strand$IncLevelDifference, as_id=minus.as.id, gene_id=rMATS.minus.strand$GeneID)

    converting.result <- rbind(plus.strand.result, minus.strand.result)
  } else if(as.type == "MXE") {
    tmp.line <- rMATS.input[1, ]

    if(length(tmp.line) != 25) {
      message("*** ERROR MESSAGE: Input file format is wrong. Please check input file.")
      return()
    }

    rMATS.input.filtering <- rMATS.input[abs(rMATS.input$IncLevelDifference) > dpsi.cutoff & rMATS.input$FDR < pvalue.cutoff, ]
  
    if(str_length(tmp.line$chr) <= 2) {
      rMATS.input.filtering$chr <- paste0("chr", rMATS.input.filtering$chr)
    }

    rMATS.plus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "+", ]
    rMATS.minus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "-", ]

    plus.as.id <- paste0(rMATS.plus.strand$chr, ":", (as.numeric(rMATS.plus.strand$upstreamES) + 1), ":", rMATS.plus.strand$upstreamEE, ":", (as.numeric(rMATS.plus.strand$X1stExonStart_0base) + 1), ":", rMATS.plus.strand$X1stExonEnd, ":", (as.numeric(rMATS.plus.strand$X2ndExonStart_0base) + 1), ":", rMATS.plus.strand$X2ndExonEnd, ":", (as.numeric(rMATS.plus.strand$downstreamES) + 1), ":", rMATS.plus.strand$downstreamEE)
    plus.upstream.exon <- paste0((as.numeric(rMATS.plus.strand$upstreamES) + 1), "-", rMATS.plus.strand$upstreamEE)
    plus.splicing.exon <- paste0((as.numeric(rMATS.plus.strand$X1stExonStart_0base) + 1), "-", rMATS.plus.strand$X1stExonEnd, ";", (as.numeric(rMATS.plus.strand$X2ndExonStart_0base) + 1), "-", rMATS.plus.strand$X2ndExonEnd)
    plus.downstream.exon <- paste0((as.numeric(rMATS.plus.strand$downstreamES) + 1), "-", rMATS.plus.strand$downstreamEE)

    plus.strand.result <- data.frame(chr=rMATS.plus.strand$chr, strand=rMATS.plus.strand$strand, as_type=as.type, upstream_exon=plus.upstream.exon, splicing_exon=plus.splicing.exon, downstream_exon=plus.downstream.exon, gene_name=rMATS.plus.strand$geneSymbol, pvalue=rMATS.plus.strand$FDR, dPSI=rMATS.plus.strand$IncLevelDifference, as_id=plus.as.id, gene_id=rMATS.plus.strand$GeneID)

    minus.as.id <- paste0(rMATS.minus.strand$chr, ":", rMATS.minus.strand$downstreamEE, ":", (as.numeric(rMATS.minus.strand$downstreamES) + 1), ":", rMATS.minus.strand$X2ndExonEnd, ":", (as.numeric(rMATS.minus.strand$X2ndExonStart_0base) + 1), ":", rMATS.minus.strand$X1stExonEnd, ":", (as.numeric(rMATS.minus.strand$X1stExonStart_0base) + 1), ":", rMATS.minus.strand$upstreamEE, ":", (as.numeric(rMATS.minus.strand$upstreamES) + 1))
    minus.upstream.exon <- paste0(rMATS.minus.strand$downstreamEE, "-", (as.numeric(rMATS.minus.strand$downstreamES) + 1))
    minus.splicing.exon <- paste0(rMATS.minus.strand$X2ndExonEnd, "-", (as.numeric(rMATS.minus.strand$X2ndExonStart_0base) + 1), ";", rMATS.minus.strand$X1stExonEnd, "-", (as.numeric(rMATS.minus.strand$X1stExonStart_0base) + 1))
    minus.downstream.exon <- paste0(rMATS.minus.strand$upstreamEE, "-", (as.numeric(rMATS.minus.strand$upstreamES) + 1))

    minus.strand.result <- data.frame(chr=rMATS.minus.strand$chr, strand=rMATS.minus.strand$strand, as_type=as.type, upstream_exon=minus.upstream.exon, splicing_exon=minus.splicing.exon, downstream_exon=minus.downstream.exon, gene_name=rMATS.minus.strand$geneSymbol, pvalue=rMATS.minus.strand$FDR, dPSI=rMATS.minus.strand$IncLevelDifference, as_id=minus.as.id, gene_id=rMATS.minus.strand$GeneID)

    converting.result <- rbind(plus.strand.result, minus.strand.result)
  } else if(as.type == "RI") {
    tmp.line <- rMATS.input[1, ]

    if(length(tmp.line) != 23) {
      message("*** ERROR MESSAGE: Input file format is wrong. Please check input file.")
      return()
    }

    rMATS.input.filtering <- rMATS.input[abs(rMATS.input$IncLevelDifference) > dpsi.cutoff & rMATS.input$FDR < pvalue.cutoff, ]
    
    if(str_length(tmp.line$chr) <= 2) {
      rMATS.input.filtering$chr <- paste0("chr", rMATS.input.filtering$chr)
    }

    rMATS.plus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "+", ]
    rMATS.minus.strand <- rMATS.input.filtering[rMATS.input.filtering$strand == "-", ]

    plus.as.id <- paste0(rMATS.plus.strand$chr, ":", (as.numeric(rMATS.plus.strand$upstreamES) + 1), ":", rMATS.plus.strand$upstreamEE, ":", (as.numeric(rMATS.plus.strand$downstreamES) + 1), ":", rMATS.plus.strand$downstreamEE)
    plus.upstream.exon <- paste0((as.numeric(rMATS.plus.strand$upstreamES) + 1), "-", rMATS.plus.strand$upstreamEE)
    plus.splicing.exon <- paste0(rMATS.plus.strand$upstreamEE, "-", rMATS.plus.strand$downstreamES)
    plus.downstream.exon <- paste0((as.numeric(rMATS.plus.strand$downstreamES) + 1), "-", rMATS.plus.strand$downstreamEE)

    plus.strand.result <- data.frame(chr=rMATS.plus.strand$chr, strand=rMATS.plus.strand$strand, as_type=as.type, upstream_exon=plus.upstream.exon, splicing_exon=plus.splicing.exon, downstream_exon=plus.downstream.exon, gene_name=rMATS.plus.strand$geneSymbol, pvalue=rMATS.plus.strand$FDR, dPSI=rMATS.plus.strand$IncLevelDifference, as_id=plus.as.id, gene_id=rMATS.plus.strand$GeneID)

    minus.as.id <- paste0(rMATS.minus.strand$chr, ":", rMATS.minus.strand$downstreamEE, ":", (as.numeric(rMATS.minus.strand$downstreamES) + 1), ":", rMATS.minus.strand$upstreamEE, ":", (as.numeric(rMATS.minus.strand$upstreamES) + 1))
    minus.upstream.exon <- paste0(rMATS.minus.strand$downstreamEE, "-", (as.numeric(rMATS.minus.strand$downstreamES) + 1))
    minus.splicing.exon <- paste0(rMATS.minus.strand$downstreamES, "-", rMATS.minus.strand$upstreamEE)
    minus.downstream.exon <- paste0(rMATS.minus.strand$upstreamEE, "-", (as.numeric(rMATS.minus.strand$upstreamES) + 1))

    minus.strand.result <- data.frame(chr=rMATS.minus.strand$chr, strand=rMATS.minus.strand$strand, as_type=as.type, upstream_exon=minus.upstream.exon, splicing_exon=minus.splicing.exon, downstream_exon=minus.downstream.exon, gene_name=rMATS.minus.strand$geneSymbol, pvalue=rMATS.minus.strand$FDR, dPSI=rMATS.minus.strand$IncLevelDifference, as_id=minus.as.id, gene_id=rMATS.minus.strand$GeneID)

    converting.result <- rbind(plus.strand.result, minus.strand.result)
  } else {
    message("*** ERROR MESSAGE: Input AS type is wrong. Please check AS_type option. We only support A5SS, A3SS, SE, MXE, or RI.")
    return()
  }

  colnames(converting.result) <- c("chr", "strand", "as_type", "upstream_exon", "splicing_exon", "downstream_exon", "gene_name", "P-value", "dPSI", "as_id", "gene_id")
  
  return(converting.result)
}
