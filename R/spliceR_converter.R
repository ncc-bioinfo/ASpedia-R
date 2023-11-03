#spliceR result file, gene model, genome version, pvalue cutoff, dIF(delta Isoform fraction (IF) value) cutoff
#' Generate ASpedia-R input format from spliceR result.
#' 
#' spliceR result convert to ASpedia-R input format.
#'
#' @param spliceR.result
#' name of spliceR result file.
#' @param pvalue.cutoff
#' value of pvalue cutoff. default value is 0.05
#' @param gene.model
#' gene model of reference. One of Refseq, Ensembl, or GENCODE.
#' @param genome.version 
#' genome version of reference. One of hg18, GRCh19, or GRCh38.
#' 
#' @return
#' converting.result
#' @export
#'
#' @examples
#' spliceR.result.file <- system.file(“extdata”, “spliceR_test.txt”, package=“ASpediaR”)
#' splilcer.converter.result <- spliceR_converter(spliceR.result.file, program=“spliceR”,
#'                                          gene.model=“Ensembl”, genome.version=“GRCh38”)

spliceR.converter <- function(spliceR.result="", pvalue.cutoff=0.05, gene.model="Ensembl", genome.version="GRCh38") {
  if(file.exists(spliceR.result) == FALSE) {
    print(paste0("*** ERROR MESSAGE: No such input file. ", spliceR.result))
    return()
  }
  
  loaded.packages <- tolower((.packages()))
  
  if(("valr" %in% loaded.packages) == FALSE) {
    library(valr)
  }
  
  if(("stringr" %in% loaded.packages) == FALSE) {
    library(stringr)
  }
  
  data.dir <- paste0(.libPaths()[1], "/ASpediaR/data")
  
  if(file.exists(data.dir) == FALSE) {
    dir.create(data.dir)
  }

  asdb.prefix <- paste0(data.dir, "/", gene.model, ".", genome.version, ".AS")

  spliceR.input <- read.table(spliceR.result, header=TRUE, stringsAsFactor=FALSE, sep="\t")
  spliceR.input <- spliceR.input[spliceR.input$spliceR.iso_p_value < pvalue.cutoff, ]
  spliceR.input[is.na(spliceR.input)] <- ""
  #spliceR.input <- spliceR.input[spliceR.input$spliceR.iso_p_value < pvalue.cutoff & abs(spliceR.input$spliceR.dIF) > dif.cutoff, ]

  total.a3ss.count <- 0
  total.a5ss.count <- 0
  total.se.count <- 0
  total.mxe.count <- 0
  total.ri.count <- 0
  total.af.count <- 0
  total.al.count <- 0

  a3ss.exon.list <- data.frame()
  a5ss.exon.list <- data.frame()
  se.exon.list <- data.frame()
  mxe.exon.list <- data.frame()
  ri.exon.list <- data.frame()
  af.exon.list <- data.frame()
  al.exon.list <- data.frame()

  for (i in 1:nrow(spliceR.input)) {
    tmp.result <- spliceR.input[i, ]

    a3ss.count <- tmp.result$spliceR.A3
    a5ss.count <- tmp.result$spliceR.A5
    se.count <- tmp.result$spliceR.ESI
    mxe.count <- tmp.result$spliceR.MEE
    ri.count <- tmp.result$spliceR.ISI
    af.count <- tmp.result$spliceR.ATSS
    al.count <- tmp.result$spliceR.ATTS

    if(a3ss.count >= 1) {
      as.exon.start.list <- strsplit(as.character(tmp.result$spliceR.A3.start), split="[;,]")[[1]]
      as.exon.end.list <- strsplit(as.character(tmp.result$spliceR.A3.end), split="[;,]")[[1]]

      #total.a3ss.count <- total.a3ss.count + a3ss.count
      total.a3ss.count <- total.a3ss.count + length(as.exon.start.list)

      for(j in 1:length(as.exon.start.list)) {
        as.exon.start <- as.exon.start.list[j]
        as.exon.end <- as.exon.end.list[j]

        tmp.a3ss <- data.frame(tmp.result$seqname, as.exon.start, as.exon.end, tmp.result$spliceR.nearest_ref, 0, tmp.result$spliceR.strand, "A3SS", tmp.result$spliceR.gene_name, tmp.result$spliceR.iso_p_value, tmp.result$spliceR.dIF)
        a3ss.exon.list <- rbind(a3ss.exon.list, tmp.a3ss)
      }
    }

    if(a5ss.count >= 1) {
      as.exon.start.list <- strsplit(as.character(tmp.result$spliceR.A5.start), split="[;,]")[[1]]
      as.exon.end.list <- strsplit(as.character(tmp.result$spliceR.A5.end), split="[;,]")[[1]]

      #total.a5ss.count <- total.a5ss.count + a5ss.count
      total.a5ss.count <- total.a5ss.count + length(as.exon.start.list)

      for(j in 1:length(as.exon.start.list)) {
        as.exon.start <- as.exon.start.list[j]
        as.exon.end <- as.exon.end.list[j]

        tmp.a5ss <- data.frame(tmp.result$seqname, as.exon.start, as.exon.end, tmp.result$spliceR.nearest_ref, 0, tmp.result$spliceR.strand, "A5SS", tmp.result$spliceR.gene_name, tmp.result$spliceR.iso_p_value, tmp.result$spliceR.dIF)
        a5ss.exon.list <- rbind(a5ss.exon.list, tmp.a5ss)
      }
    }

    if(se.count >= 1) {
      as.exon.start.list <- strsplit(as.character(tmp.result$spliceR.ESI.start), split="[;,]")[[1]]
      as.exon.end.list <- strsplit(as.character(tmp.result$spliceR.ESI.end), split="[;,]")[[1]]

      #total.se.count <- total.se.count + se.count
      total.se.count <- total.se.count + length(as.exon.start.list)

      for(j in 1:length(as.exon.start.list)) {
        as.exon.start <- as.exon.start.list[j]
        as.exon.end <- as.exon.end.list[j]

        tmp.se <- data.frame(tmp.result$seqname, as.exon.start, as.exon.end, tmp.result$spliceR.nearest_ref, 0, tmp.result$spliceR.strand, "SE", tmp.result$spliceR.gene_name, tmp.result$spliceR.iso_p_value, tmp.result$spliceR.dIF)
        se.exon.list <- rbind(se.exon.list, tmp.se)
      }
    }

    if(mxe.count >= 1) {
      as.exon.start.list <- strsplit(as.character(tmp.result$spliceR.MEE.start), split="[;,]")[[1]]
      as.exon.end.list <- strsplit(as.character(tmp.result$spliceR.MEE.end), split="[;,]")[[1]]

      total.mxe.count <- total.mxe.count + length(as.exon.start.list)

      for(j in 1:length(as.exon.start.list)) {
        as.exon.start <- as.exon.start.list[j]
        as.exon.end <- as.exon.end.list[j]

        tmp.mxe <- data.frame(tmp.result$seqname, as.exon.start, as.exon.end, tmp.result$spliceR.nearest_ref, 0, tmp.result$spliceR.strand, "MXE", tmp.result$spliceR.gene_name, tmp.result$spliceR.iso_p_value, tmp.result$spliceR.dIF)
        mxe.exon.list <- rbind(mxe.exon.list, tmp.mxe)
      }
    }

    if(ri.count >= 1) {
      as.exon.start.list <- strsplit(as.character(tmp.result$spliceR.ISI.start), split="[;,]")[[1]]
      as.exon.end.list <- strsplit(as.character(tmp.result$spliceR.ISI.end), split="[;,]")[[1]]

      total.ri.count <- total.ri.count + length(as.exon.start.list)

      for(j in 1:length(as.exon.start.list)) {
        as.exon.start <- as.exon.start.list[j]
        as.exon.end <- as.exon.end.list[j]

        tmp.ri <- data.frame(tmp.result$seqname, as.exon.start, as.exon.end, tmp.result$spliceR.nearest_ref, 0, tmp.result$spliceR.strand, "RI", tmp.result$spliceR.gene_name, tmp.result$spliceR.iso_p_value, tmp.result$spliceR.dIF)
        ri.exon.list <- rbind(ri.exon.list, tmp.ri)
      }

    }

    if(af.count >= 1) {
      as.exon.start.list <- strsplit(as.character(tmp.result$spliceR.ATSS.start), split="[;,]")[[1]]
      as.exon.end.list <- strsplit(as.character(tmp.result$spliceR.ATSS.end), split="[;,]")[[1]]

      total.af.count <- total.af.count + length(as.exon.start.list)

      for(j in 1:length(as.exon.start.list)) {
        as.exon.start <- as.exon.start.list[j]
        as.exon.end <- as.exon.end.list[j]

        tmp.af <- data.frame(tmp.result$seqname, as.exon.start, as.exon.end, tmp.result$spliceR.nearest_ref, 0, tmp.result$spliceR.strand, "AF", tmp.result$spliceR.gene_name, tmp.result$spliceR.iso_p_value, tmp.result$spliceR.dIF)
        af.exon.list <- rbind(af.exon.list, tmp.af)
      }
    }

    if(al.count >= 1) {
      as.exon.start.list <- strsplit(as.character(tmp.result$spliceR.ATTS.start), split="[;,]")[[1]]
      as.exon.end.list <- strsplit(as.character(tmp.result$spliceR.ATTS.end), split="[;,]")[[1]]

      total.al.count <- total.al.count + length(as.exon.start.list)

      for(j in 1:length(as.exon.start.list)) {
        as.exon.start <- as.exon.start.list[j]
        as.exon.end <- as.exon.end.list[j]

        tmp.al <- data.frame(tmp.result$seqname, as.exon.start, as.exon.end, tmp.result$spliceR.nearest_ref, 0, tmp.result$spliceR.strand, "AL", tmp.result$spliceR.gene_name, tmp.result$spliceR.iso_p_value, tmp.result$spliceR.dIF)
        al.exon.list <- rbind(al.exon.list, tmp.al)
      }
    }
  }

  spliceR.inter.colname <- c("chrom", "start", "end", "transcript_id", "tmp", "strand", "as_type", "gene", "iso.pvalue", "dIF")

  asdb.inter.colname <- c("chrom", "start", "end", "asid", "tmp", "strand", "as_type", "gene_model", "genome_version")

  a3ss.as.file.name <- paste0(asdb.prefix, ".A3SS.exon.region.sorted.bed")
  a5ss.as.file.name <- paste0(asdb.prefix, ".A5SS.exon.region.sorted.bed")
  se.as.file.name <- paste0(asdb.prefix, ".SE.exon.region.sorted.bed")
  mxe.as.file.name <- paste0(asdb.prefix, ".MXE.exon.region.sorted.bed")
  ri.as.file.name <- paste0(asdb.prefix, ".RI.exon.region.sorted.bed")
  af.as.file.name <- paste0(asdb.prefix, ".AF.exon.region.sorted.bed")
  al.as.file.name <- paste0(asdb.prefix, ".AL.exon.region.sorted.bed")

  if(!(file.exists(a3ss.as.file.name))) {
    url <- paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/split_by_as_type/", gene.model, ".", genome.version, ".AS.A3SS.exon.region.sorted.bed")
    download.file(url, a3ss.as.file.name, method="auto")
  }

  if(!(file.exists(a5ss.as.file.name))) {
    url <- paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/split_by_as_type/", gene.model, ".", genome.version, ".AS.A5SS.exon.region.sorted.bed")
    download.file(url, a5ss.as.file.name, method="auto")
  }

  if(!(file.exists(se.as.file.name))) {
    url <- paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/split_by_as_type/", gene.model, ".", genome.version, ".AS.SE.exon.region.sorted.bed")
    download.file(url, se.as.file.name, method="auto")
  }

  if(!(file.exists(mxe.as.file.name))) {
    url <- paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/split_by_as_type/", gene.model, ".", genome.version, ".AS.MXE.exon.region.sorted.bed")
    download.file(url, mxe.as.file.name, method="auto")
  }

  if(!(file.exists(ri.as.file.name))) {
    url <- paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/split_by_as_type/", gene.model, ".", genome.version, ".AS.RI.exon.region.sorted.bed")
    download.file(url, ri.as.file.name, method="auto")
  }

  if(!(file.exists(af.as.file.name))) {
    url <- paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/split_by_as_type/", gene.model, ".", genome.version, ".AS.AF.exon.region.sorted.bed")
    download.file(url, af.as.file.name, method="auto")
  }

  if(!(file.exists(al.as.file.name))) {
    url <- paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/split_by_as_type/", gene.model, ".", genome.version, ".AS.AL.exon.region.sorted.bed")
    download.file(url, al.as.file.name, method="auto")
  }

  a3ss.as <- read.table(a3ss.as.file.name, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  a5ss.as <- read.table(a5ss.as.file.name, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  se.as <- read.table(se.as.file.name, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  mxe.as <- read.table(mxe.as.file.name, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  ri.as <- read.table(ri.as.file.name, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  af.as <- read.table(af.as.file.name, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  al.as <- read.table(al.as.file.name, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  #colnames(a3ss.as) <- asdb.colname

  inter.colname <- c("asdb.chr", "asdb.start", "asdb.end", "asdb.asid", "asdb.tmp", "asdb.strand", "asdb.as_type", "asdb.gene_model", "asdb.genome_version", "spliceR.start", "spliceR.end", "spliceR.transcript_id", "spliceR.tmp", "spliceR.strand", "spliceR.as_type", "spliceR.gene", "spliceR.iso.pvalue", "spliceR.dIF", "overlap")

  #"chr", "strand", "as_type", "upstream_exon", "splicing exon", "downstream exon", "gene_id", "gene_name", "P-value", "dPSI", "as_id"
  #A3SS
  colnames(a3ss.as) <- asdb.inter.colname
  a3ss.as$start <- as.numeric(a3ss.as$start)
  a3ss.as$end <- as.numeric(a3ss.as$end)

  colnames(a3ss.exon.list) <- spliceR.inter.colname
  a3ss.exon.list$start <- as.numeric(a3ss.exon.list$start)
  a3ss.exon.list$end <- as.numeric(a3ss.exon.list$end)

  a3ss.inter.result <- as.data.frame(bed_intersect(a3ss.as, a3ss.exon.list))
  colnames(a3ss.inter.result) <- inter.colname

  a3ss.as.result <- a3ss.inter.result[a3ss.inter.result$asdb.start <= a3ss.inter.result$spliceR.start & a3ss.inter.result$asdb.end >= a3ss.inter.result$spliceR.end, ]

  a3ss.as.result <- a3ss.as.result[, c("asdb.chr", "asdb.strand", "asdb.as_type", "spliceR.gene", "spliceR.iso.pvalue", "spliceR.dIF", "asdb.asid")]
  a3ss.asid <- a3ss.as.result$asdb.asid
  a3ss.asid.list <- as.data.frame(t(as.data.frame((str_split(a3ss.asid, ":")))))
  a3ss.upstream.exon <- paste0(a3ss.asid.list$V2, "-", a3ss.asid.list$V3)
  a3ss.splicing.exon <- paste0(a3ss.asid.list$V4, "-", a3ss.asid.list$V6)
  a3ss.downstream.exon <- paste0(a3ss.asid.list$V5, "-", a3ss.asid.list$V6)

  a3ss.result <- data.frame(chr=a3ss.as.result$asdb.chr, strand=a3ss.as.result$asdb.strand, as_type=a3ss.as.result$asdb.as_type, upstream_exon=a3ss.upstream.exon, splicing_exon=a3ss.splicing.exon, downstream_exon=a3ss.downstream.exon, gene_name=a3ss.as.result$spliceR.gene, pvalue=a3ss.as.result$spliceR.iso.pvalue, dIF=a3ss.as.result$spliceR.dIF, as_id=a3ss.as.result$asdb.asid)


  #A5SS
  colnames(a5ss.as) <- asdb.inter.colname
  a5ss.as$start <- as.numeric(a5ss.as$start)
  a5ss.as$end <- as.numeric(a5ss.as$end)

  colnames(a5ss.exon.list) <- spliceR.inter.colname
  a5ss.exon.list$start <- as.numeric(a5ss.exon.list$start)
  a5ss.exon.list$end <- as.numeric(a5ss.exon.list$end)

  a5ss.inter.result <- as.data.frame(bed_intersect(a5ss.as, a5ss.exon.list))
  colnames(a5ss.inter.result) <- inter.colname

  a5ss.as.result <- a5ss.inter.result[a5ss.inter.result$asdb.start <= a5ss.inter.result$spliceR.start & a5ss.inter.result$asdb.end >= a5ss.inter.result$spliceR.end,]

  a5ss.as.result <- a5ss.as.result[, c("asdb.chr", "asdb.strand", "asdb.as_type", "spliceR.gene", "spliceR.iso.pvalue", "spliceR.dIF", "asdb.asid")]
  a5ss.asid <- a5ss.as.result$asdb.asid
  a5ss.asid.list <- as.data.frame(t(as.data.frame((str_split(a5ss.asid, ":")))))
  a5ss.upstream.exon <- paste0(a5ss.asid.list$V2, "-", a5ss.asid.list$V3)
  a5ss.splicing.exon <- paste0(a5ss.asid.list$V2, "-", a5ss.asid.list$V4)
  a5ss.downstream.exon <- paste0(a5ss.asid.list$V5, "-", a5ss.asid.list$V6)

  a5ss.result <- data.frame(chr=a5ss.as.result$asdb.chr, strand=a5ss.as.result$asdb.strand, as_type=a5ss.as.result$asdb.as_type, upstream_exon=a5ss.upstream.exon, splicing_exon=a5ss.splicing.exon, downstream_exon=a5ss.downstream.exon, gene_name=a5ss.as.result$spliceR.gene, pvalue=a5ss.as.result$spliceR.iso.pvalue, dIF=a5ss.as.result$spliceR.dIF, as_id=a5ss.as.result$asdb.asid)

  #SE
  colnames(se.as) <- asdb.inter.colname
  se.as$start <- as.numeric(se.as$start)
  se.as$end <- as.numeric(se.as$end)

  colnames(se.exon.list) <- spliceR.inter.colname
  se.exon.list$start <- as.numeric(se.exon.list$start)
  se.exon.list$end <- as.numeric(se.exon.list$end)

  se.inter.result <- as.data.frame(bed_intersect(se.as, se.exon.list))
  colnames(se.inter.result) <- inter.colname

  se.as.result <- se.inter.result[se.inter.result$asdb.start == se.inter.result$spliceR.start & se.inter.result$asdb.end == se.inter.result$spliceR.end,]

  se.as.result <- se.as.result[, c("asdb.chr", "asdb.strand", "asdb.as_type", "spliceR.gene", "spliceR.iso.pvalue", "spliceR.dIF", "asdb.asid")]
  se.asid <- se.as.result$asdb.asid
  se.asid.list <- as.data.frame(t(as.data.frame((str_split(se.asid, ":")))))
  se.upstream.exon <- paste0(se.asid.list$V2, "-", se.asid.list$V3)
  se.splicing.exon <- paste0(se.asid.list$V4, "-", se.asid.list$V5)
  se.downstream.exon <- paste0(se.asid.list$V6, "-", se.asid.list$V7)

  se.result <- data.frame(chr=se.as.result$asdb.chr, strand=se.as.result$asdb.strand, as_type=se.as.result$asdb.as_type, upstream_exon=se.upstream.exon, splicing_exon=se.splicing.exon, downstream_exon=se.downstream.exon, gene_name=se.as.result$spliceR.gene, pvalue=se.as.result$spliceR.iso.pvalue, dIF=se.as.result$spliceR.dIF, as_id=se.as.result$asdb.asid)

  #MXE
  colnames(mxe.as) <- asdb.inter.colname
  mxe.as$start <- as.numeric(mxe.as$start)
  mxe.as$end <- as.numeric(mxe.as$end)

  colnames(mxe.exon.list) <- spliceR.inter.colname
  mxe.exon.list$start <- as.numeric(mxe.exon.list$start)
  mxe.exon.list$end <- as.numeric(mxe.exon.list$end)

  mxe.inter.result <- as.data.frame(bed_intersect(mxe.as, mxe.exon.list))
  colnames(mxe.inter.result) <- inter.colname

  mxe.as.result <- mxe.inter.result[mxe.inter.result$asdb.start == mxe.inter.result$spliceR.start & mxe.inter.result$asdb.end == mxe.inter.result$spliceR.end,]

  mxe.as.result <- mxe.as.result[, c("asdb.chr", "asdb.strand", "asdb.as_type", "spliceR.gene", "spliceR.iso.pvalue", "spliceR.dIF", "asdb.asid")]
  mxe.asid <- mxe.as.result$asdb.asid
  mxe.asid.list <- as.data.frame(t(as.data.frame((str_split(mxe.asid, ":")))))
  mxe.upstream.exon <- paste0(mxe.asid.list$V2, "-", mxe.asid.list$V3)
  mxe.splicing.exon <- paste0(mxe.asid.list$V4, "-", mxe.asid.list$V5, ";", mxe.asid.list$V6, "-", mxe.asid.list$V7)
  mxe.downstream.exon <- paste0(mxe.asid.list$V8, "-", mxe.asid.list$V9)

  mxe.result <- data.frame(chr=mxe.as.result$asdb.chr, strand=mxe.as.result$asdb.strand, as_type=mxe.as.result$asdb.as_type, upstream_exon=mxe.upstream.exon, splicing_exon=mxe.splicing.exon, downstream_exon=mxe.downstream.exon, gene_name=mxe.as.result$spliceR.gene, pvalue=mxe.as.result$spliceR.iso.pvalue, dIF=mxe.as.result$spliceR.dIF, as_id=mxe.as.result$asdb.asid)

  #RI
  colnames(ri.as) <- asdb.inter.colname
  ri.as$start <- as.numeric(ri.as$start)
  ri.as$end <- as.numeric(ri.as$end)

  colnames(ri.exon.list) <- spliceR.inter.colname
  ri.exon.list$start <- as.numeric(ri.exon.list$start)
  ri.exon.list$end <- as.numeric(ri.exon.list$end)

  ri.inter.result <- as.data.frame(bed_intersect(ri.as, ri.exon.list))
  colnames(ri.inter.result) <- inter.colname

  ri.as.result <- ri.inter.result[ri.inter.result$asdb.start <= ri.inter.result$spliceR.start & ri.inter.result$asdb.end >= ri.inter.result$spliceR.end,]

  ri.as.result <- ri.as.result[, c("asdb.chr", "asdb.strand", "asdb.as_type", "spliceR.gene", "spliceR.iso.pvalue", "spliceR.dIF", "asdb.asid")]
  ri.asid <- ri.as.result$asdb.asid
  ri.asid.list <- as.data.frame(t(as.data.frame((str_split(ri.asid, ":")))))
  ri.upstream.exon <- paste0(ri.asid.list$V2, "-", ri.asid.list$V3)
  ri.splicing.exon <- paste0(ri.asid.list$V3, "-", ri.asid.list$V4)
  ri.downstream.exon <- paste0(ri.asid.list$V4, "-", ri.asid.list$V5)

  ri.result <- data.frame(chr=ri.as.result$asdb.chr, strand=ri.as.result$asdb.strand, as_type=ri.as.result$asdb.as_type, upstream_exon=ri.upstream.exon, splicing_exon=ri.splicing.exon, downstream_exon=ri.downstream.exon, gene_name=ri.as.result$spliceR.gene, pvalue=ri.as.result$spliceR.iso.pvalue, dIF=ri.as.result$spliceR.dIF, as_id=ri.as.result$asdb.asid)

  #AF
  colnames(af.as) <- asdb.inter.colname
  af.as$start <- as.numeric(af.as$start)
  af.as$end <- as.numeric(af.as$end)

  colnames(af.exon.list) <- spliceR.inter.colname
  af.exon.list$start <- as.numeric(af.exon.list$start)
  af.exon.list$end <- as.numeric(af.exon.list$end)

  af.inter.result <- as.data.frame(bed_intersect(af.as, af.exon.list))
  colnames(af.inter.result) <- inter.colname

  af.as.result <- af.inter.result[af.inter.result$asdb.start == af.inter.result$spliceR.start & af.inter.result$asdb.end == af.inter.result$spliceR.end,]

  af.as.result <- af.as.result[, c("asdb.chr", "asdb.strand", "asdb.as_type", "spliceR.gene", "spliceR.iso.pvalue", "spliceR.dIF", "asdb.asid")]
  af.asid <- af.as.result$asdb.asid
  af.asid.list <- as.data.frame(t(as.data.frame((str_split(af.asid, ":")))))
  af.upstream.exon <- paste0(af.asid.list$V2, "-", af.asid.list$V3)
  af.splicing.exon <- paste0(af.asid.list$V4, "-", af.asid.list$V5)
  af.downstream.exon <- paste0(af.asid.list$V6, "-", af.asid.list$V7)

  af.result <- data.frame(chr=af.as.result$asdb.chr, strand=af.as.result$asdb.strand, as_type=af.as.result$asdb.as_type, upstream_exon=af.upstream.exon, splicing_exon=af.splicing.exon, downstream_exon=af.downstream.exon, gene_name=af.as.result$spliceR.gene, pvalue=af.as.result$spliceR.iso.pvalue, dIF=af.as.result$spliceR.dIF, as_id=af.as.result$asdb.asid)

  #AL
  colnames(al.as) <- asdb.inter.colname
  al.as$start <- as.numeric(al.as$start)
  al.as$end <- as.numeric(al.as$end)

  colnames(al.exon.list) <- spliceR.inter.colname
  al.exon.list$start <- as.numeric(al.exon.list$start)
  al.exon.list$end <- as.numeric(al.exon.list$end)

  al.inter.result <- as.data.frame(bed_intersect(al.as, al.exon.list))
  colnames(al.inter.result) <- inter.colname

  al.as.result <- al.inter.result[al.inter.result$asdb.start == al.inter.result$spliceR.start & al.inter.result$asdb.end == al.inter.result$spliceR.end,]

  al.as.result <- al.as.result[, c("asdb.chr", "asdb.strand", "asdb.as_type", "spliceR.gene", "spliceR.iso.pvalue", "spliceR.dIF", "asdb.asid")]
  al.asid <- al.as.result$asdb.asid
  al.asid.list <- as.data.frame(t(as.data.frame((str_split(al.asid, ":")))))
  al.upstream.exon <- paste0(al.asid.list$V2, "-", al.asid.list$V3)
  al.splicing.exon <- paste0(al.asid.list$V4, "-", al.asid.list$V5)
  al.downstream.exon <- paste0(al.asid.list$V6, "-", al.asid.list$V7)

  al.result <- data.frame(chr=al.as.result$asdb.chr, strand=al.as.result$asdb.strand, as_type=al.as.result$asdb.as_type, upstream_exon=al.upstream.exon, splicing_exon=al.splicing.exon, downstream_exon=al.downstream.exon, gene_name=al.as.result$spliceR.gene, pvalue=al.as.result$spliceR.iso.pvalue, dIF=al.as.result$spliceR.dIF, as_id=al.as.result$asdb.asid)

  result.data <- rbind(a3ss.result, a5ss.result, se.result, mxe.result, ri.result, af.result, al.result)

  result.data$spliceR.strand <- result.data$asdb.strand

  return(result.data)
}
