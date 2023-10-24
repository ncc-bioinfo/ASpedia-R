#SUPPA result file, gtf file, ioe file, gene model, P-value cutoff, DPSI cutoff

#input.file.name <- "/home/aspedia/project/01.Util.Package/01.test.data/SUPPA/origin/ESRP_dPSI_empirical.dpsi"
#input.file.name <- "/home/aspedia/project/01.Util.Package/01.test.data/SUPPA/origin/RBM47_dPSI_empirical.dpsi"
#gtf.file.name <-  "/home/aspedia/project/01.Util.Package/00.reference/Homo_sapiens.GRCh38.99.rm.mt.add.chr.gtf"
#ioe.file.name <- "/home/aspedia/project/01.Util.Package/00.reference/SUPPA/Homo_sapiens.GRCh38.99.rm.mt.add.chr.ioe"
#gene.model <- "Ensembl"
#pvalue.cutoff <- 0.05
#dpsi.cutoff <- 0.1
#result.file.name <- "/home/aspedia/project/01.Util.Package/01.test.data/SUPPA/origin/ESRP_dPSI_empirical_dPSI_0.1_pvalue_0.05_result.txt"
#result.file.name <- "/home/aspedia/project/01.Util.Package/01.test.data/SUPPA/origin/RBM47_dPSI_empirical_dPSI_0.1_pvalue_0.05_result.txt"
library(rtracklayer)

#' Title
#'
#' @param input.file.name 
#' @param gtf.file.name 
#' @param ioe.file.name 
#' @param pvalue.cutoff 
#' @param dpsi.cutoff 
#'
#' @return
#' @export
#'
#' @examples
SUPPA.converter <- function(input.file.name, gtf.file.name, ioe.file.name, pvalue.cutoff, dpsi.cutoff) {
  if(file.exists(input.file.name) == FALSE) {
    print(paste0("*** ERROR MESSAGE: No such input file. ", input.file.name))
    return()
  }

  if(file.exists(gtf.file.name) == FALSE) {
    print(paste0("*** ERROR MESSAGE: No such gtf file. ", gtf.file.name))
    return()
  }

  if(file.exists(ioe.file.name) == FALSE) {
    print(paste0("*** ERROR MESSAGE: No such .ioe file. ", ioe.file.name))
    return()
  }

  #load gtf
  #gtf.data <- data.frame()
  #gtf.data <- import(gtf.file.name)
  reference.start.time <- Sys.time()
  gtf.data <- readGFF(gtf.file.name)
  reference.end.time <- Sys.time()
  print("SUPPA reference load...")
  print(reference.end.time - reference.start.time)

  #if (gene.model == "GENCODE" | gene.model == "Refseq" | gene.model == "UCSC") {
  #  gtf.data <- readGFF(gtf.file.name)
  #}else {
  #  library(refGenome)
  #  genome.gtf <- ensemblGenome()
  #  read.gtf(genome.gtf, gtf.file.name, useBasedir=FALSE, sep="\t")
  #  gtf.data <- getGtf(genome.gtf)
  #}

  gtf.gene.id <- sapply(gtf.data$gene_id, function(x) {
    strsplit(x, split="[.]")[[1]][1] })

  gtf.data$gene_id <- gtf.gene.id

  gtf.transcript.id <- sapply(gtf.data$transcript_id, function(x) {
    strsplit(x, split="[.]")[[1]][1] })

  gtf.data$transcript_id <- gtf.transcript.id

  #head(gtf.data)

  #load suppa .ioe file
  suppa.ioe <- read.table(ioe.file.name, header=TRUE, stringsAsFactor=FALSE, sep="\t")

  #load suppa result file
  suppa.result <- read.table(input.file.name, header=TRUE, stringsAsFactor=FALSE, sep="\t")

  colnames(suppa.result) <- c("dPSI", "pvalue")
  suppa.result <- suppa.result[(abs(suppa.result$dPSI) > dpsi.cutoff) & (suppa.result$pvalue < pvalue.cutoff), ]
  print(nrow(suppa.result))
  suppa.event.id.list <- rownames(suppa.result)
  rownames(suppa.result) <- NULL

  as.event.list <- sapply(suppa.event.id.list, function(x) {
    tmp.event.id.list <- strsplit(x, split="[;]")[[1]]
    tmp.gene.id.list <- strsplit(tmp.event.id.list[1], split="[:]")[[1]]

    as.gene.id <- ""
    as.gene.name <- ""

    if(length(tmp.gene.id.list) == 2)
    {
      as.gene.id <- strsplit(tmp.gene.id.list[1], split="[.]")[[1]][1]
      as.gene.name <- tmp.gene.id.list[2]
    }else if(length(tmp.gene.id.list) == 1)
    {
      as.gene.id <- strsplit(tmp.gene.id.list[1], split="[.]")[[1]][1]
      as.gene.name <- unique(gtf.data[gtf.data$gene_id == as.gene.id, "gene_name"])
    }else
    {
      print(tmp.gene.id.list)
    }

    event.id <- tmp.event.id.list[2]
    c(as.gene.id, as.gene.name, event.id) })

  as.event.list <- as.data.frame(t(as.event.list), stringsAsFactors=FALSE)
  rownames(as.event.list) <- NULL

  suppa.event.id.list <- data.frame(suppa_id=suppa.event.id.list, stringsAsFactors=FALSE)

  suppa.result <- cbind(suppa.event.id.list, as.event.list, suppa.result)
  colnames(suppa.result) <- c("suppa_event_id", "gene_id", "gene_name", "as_id", "dPSI", "pvalue")
  rownames(suppa.result) <- NULL

  as.gene.id.list <- suppa.result$gene_id

  as.gene.id.list <- unique(as.gene.id.list)
  #head(as.result)

  #extract as related genes
  gtf.as.gene <- data.frame()
  gtf.as.gene <- gtf.data[(gtf.data$gene_id %in% as.gene.id.list & gtf.data$type == "exon"), ]

  #if (gene.model == "GENCODE" | gene.model == "Refseq" | gene.model == "UCSC") {
  #  gtf.as.gene <- gtf.data[(gtf.data$gene_id %in% as.gene.id.list & gtf.data$type == "exon"), ]
  #} else {
  #  gtf.as.gene <- gtf.data[(gtf.data$gene_id %in% as.gene.id.list & gtf.data$feature == "exon"), ]
  #}

  as.dump.result <- data.frame()

  for (i in 1:nrow(suppa.result)) {
    tmp.result <- suppa.result[i, ]
    suppa.event.id <- tmp.result$as_id
    gene.id <- tmp.result$gene_id
    gene.name <- tmp.result$gene_name
    tmp.pvalue <- tmp.result$pvalue
    tmp.dpsi <- tmp.result$dPSI

    tmp.ioe <- suppa.ioe[suppa.ioe$event_id == tmp.result$suppa_event_id, ]

    if(nrow(tmp.ioe) == 0)
    {
      #print(tmp.ioe)
      next
    }

    #extract transcript list
    as.transcript <- strsplit(tmp.ioe$alternative_transcripts, ",")[[1]]
    as.transcript <- as.vector(sapply(as.transcript, function(x) {
      strsplit(x, split="[.]")[[1]][1] }))

    total.transcript <- strsplit(tmp.ioe$total_transcripts, ",")[[1]]
    total.transcript <- as.vector(sapply(total.transcript, function(x) {
      strsplit(x, split="[.]")[[1]][1] }))

    diff.transcript <- setdiff(total.transcript, as.transcript)
    tid.set <- ""

    as.transcript.gtf <- gtf.as.gene[gtf.as.gene$transcript_id %in% as.transcript, ]
    diff.transcript.gtf <- gtf.as.gene[gtf.as.gene$transcript_id %in% diff.transcript, ]

    #make as id
    #"chr", "strand", "as_type", "upstream_exon", "splicing exon", "downstream exon", "gene_id", "gene_name", "P-value", "dPSI"
    event.id.list <- strsplit(suppa.event.id, split="[:]")[[1]]
    as.type <- event.id.list[1]
    chrom <- event.id.list[2]
    strand <- event.id.list[length(event.id.list)]

    result.as.type <- ""

    first.start <- ""
    last.end <- ""
    as.id <- ""
    as.id.list <- data.frame()
    as.id.result <- data.frame()

    #print(as.type)

    if(as.type == "A3") {
      #event id => e1-s2:e1-s3(+) , e2-s3:e1-s3(-)
      #as id => s1:e1:s2:s3:e2 or s1:e1:s2:s3:e3(+) , e3:s3:e2:e1:s2 or e3:s3:e2:e1:s1(-)

      result.as.type <- "A3SS"

      exon.list.1 <- strsplit(event.id.list[3], "-")[[1]]
      exon.list.2 <- strsplit(event.id.list[4], "-")[[1]]

      if(strand == "+") {
        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.1[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[1], exon.list.1[2], exon.list.2[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[2], "-", last.end), downstream_exon=paste0(exon.list.2[2], "-", last.end)))
        }

        for(tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.1[1]), "start"])
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.2[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[1], exon.list.1[2], exon.list.2[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[2], "-", last.end), downstream_exon=paste0(exon.list.2[2], "-", last.end)))
        }
      }else if(strand == "-") {
        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.1[2]), "end"])
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[2], exon.list.1[1], exon.list.2[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.1[2]), splicing_exon=paste0(exon.list.1[1], "-", last.end), downstream_exon=paste0(exon.list.2[1], "-", last.end)))
        }

        for(tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.1[2]), "end"])
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.2[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[2], exon.list.1[1], exon.list.2[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.1[2]), splicing_exon=paste0(exon.list.1[1], "-", last.end), downstream_exon=paste0(exon.list.2[1], "-", last.end)))
        }
      }
    }else if(as.type == "A5") {
      #event id => e2-s3:e1-s3(+) , e1-s2:e1-s3(-)
      #as id => s1:e1:e2:s3:e3 or s2:e1:e2:s3:e3(+) , e2:s3:s2:e1:s1 or e3:s3:s2:e1:s1(-)

      result.as.type <- "A5SS"

      exon.list.1 <- strsplit(event.id.list[3], "-")[[1]]
      exon.list.2 <- strsplit(event.id.list[4], "-")[[1]]

      if(strand == "+") {
        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.1[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.2[1], exon.list.1[1], exon.list.1[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.2[1]), splicing_exon=paste0(first.start, "-", exon.list.1[1]), downstream_exon=paste0(exon.list.1[2], "-", last.end)))
        }

        for(tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.2[1]), "start"])
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.1[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.2[1], exon.list.1[1], exon.list.1[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.2[1]), splicing_exon=paste0(first.start, "-", exon.list.1[1]), downstream_exon=paste0(exon.list.1[2], "-", last.end)))
        }
      }else if(strand == "-") {
        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.1[2]), "end"])
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.2[2], exon.list.1[2], exon.list.1[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.2[2]), splicing_exon=paste0(first.start, "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[1], "-", last.end)))
        }

        for(tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.2[2]), "end"])
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.2[2], exon.list.1[2], exon.list.1[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.2[2]), splicing_exon=paste0(first.start, "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[1], "-", last.end)))
        }
      }
    }else if(as.type =="SE") {
      #event id => e1-s2:e2-s3(+) , e1-s2:e2-s3(-)
      #as id => s1:e1:s2:e2:s3:e3(+) , e3:s3:e2:s2:e1:s1(-)

      result.as.type <- "SE"

      exon.list.1 <- strsplit(event.id.list[3], "-")[[1]]
      exon.list.2 <- strsplit(event.id.list[4], "-")[[1]]

      if(strand == "+") {
        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.2[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[1], exon.list.1[2], exon.list.2[1], exon.list.2[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[2], "-", exon.list.2[1]), downstream_exon=paste0(exon.list.2[2], "-", last.end)))
        }

        for (tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.1[1]), "start"])
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.2[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[1], exon.list.1[2], exon.list.2[1], exon.list.2[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[2], "-", exon.list.2[1]), downstream_exon=paste0(exon.list.2[2], "-", last.end)))
        }
      }else if(strand == "-") {
        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.2[2]), "end"])
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.2[2], exon.list.2[1], exon.list.1[2], exon.list.1[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.2[2]), splicing_exon=paste0(exon.list.2[1], "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[1], "-", last.end)))
        }

        for (tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.2[2]), "end"])
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.2[2], exon.list.2[1], exon.list.1[2], exon.list.1[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.2[2]), splicing_exon=paste0(exon.list.2[1], "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[1], "-", last.end)))
        }
      }
    }else if(as.type == "MX") {
      #event id => e1-s2:e2-s4:e1-s3:e3-s4(+) , e1-s2:e2-s4:e1-s3:e3-s4(-)
      #as id => s1:e2:s2:e2:s3:e3:s4:e4(+) , e4:s4:e3:s3:e2:s2:e1:s1(-)

      result.as.type <- "MXE"

      exon.list.1 <- strsplit(event.id.list[3], "-")[[1]]
      exon.list.2 <- strsplit(event.id.list[4], "-")[[1]]
      exon.list.3 <- strsplit(event.id.list[5], "-")[[1]]
      exon.list.4 <- strsplit(event.id.list[6], "-")[[1]]

      if(exon.list.1[1] != exon.list.3[1]) {
        next
      }

      if(exon.list.2[2] != exon.list.4[2]) {
        next
      }

      if(strand == "+") {
        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.4[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[1], exon.list.1[2], exon.list.2[1], exon.list.3[2], exon.list.4[1], exon.list.4[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[2], "-", exon.list.2[1], ";", exon.list.3[2], "-", exon.list.4[1]), downstream_exon=paste0(exon.list.4[2], "-", last.end)))
        }

        for (tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.1[1]), "start"])
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.4[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[1], exon.list.1[2], exon.list.2[1], exon.list.3[2], exon.list.4[1], exon.list.4[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[2], "-", exon.list.2[1], ";", exon.list.3[2], "-", exon.list.4[1]), downstream_exon=paste0(exon.list.4[2], "-", last.end)))
        }
      }else if(strand == "-") {
        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.4[2]), "end"])
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.4[2], exon.list.4[1], exon.list.3[2], exon.list.2[1], exon.list.1[2], exon.list.1[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.4[2]), splicing_exon=paste0(exon.list.4[1], "-", exon.list.3[2], ";", exon.list.2[1], "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[1], "-", last.end)))
        }

        for (tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.4[2]), "end"])
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.4[2], exon.list.4[1], exon.list.3[2], exon.list.2[1], exon.list.1[2], exon.list.1[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.4[2]), splicing_exon=paste0(exon.list.4[1], "-", exon.list.3[2], ";", exon.list.2[1], "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[1], "-", last.end)))
        }
      }
    }else if(as.type == "RI") {
      #event id => s1:e1-s2:e2(+) , s1:e1-s2:e2(-)
      #as id => s1:e1:s2:e2(+) , e2:s2:e1:s1(-)

      result.as.type <- "RI"

      exon.list.1 <- strsplit(event.id.list[4], "-")[[1]]

      if(strand == "+") {
        for (tmp.transcript in as.transcript) {
          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, event.id.list[3], exon.list.1[1], exon.list.1[2], event.id.list[5], sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(event.id.list[3], "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[1], "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[2], "-", event.id.list[5])))
        }

        for (tmp.transcript in diff.transcript) {
          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, event.id.list[3], exon.list.1[1], exon.list.1[2], event.id.list[5], sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(event.id.list[3], "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[1], "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[2], "-", event.id.list[5])))
        }
      }else if(strand == "-") {
        for (tmp.transcript in as.transcript) {
          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, event.id.list[5], exon.list.1[2], exon.list.1[1], event.id.list[3], sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(event.id.list[5], "-", exon.list.1[2]), splicing_exon=paste0(exon.list.1[2], "-", exon.list.1[1]), downstream_exon=paste0(exon.list.1[1], "-", event.id.list[3])))
        }

        for (tmp.transcript in diff.transcript) {
          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, event.id.list[5], exon.list.1[2], exon.list.1[1], event.id.list[3], sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(event.id.list[5], "-", exon.list.1[2]), splicing_exon=paste0(exon.list.1[2], "-", exon.list.1[1]), downstream_exon=paste0(exon.list.1[1], "-", event.id.list[3])))
        }
      }
    }else if(as.type == "AF") {
      #event id => s1:e1-s3:s2:e2-s3(+) , e1-s2:e2:e1-s3:e3(-)
      #as id => s1:e1:s2:e2:s3:e3(+) , e3:s3:e2:s2:e1:s1(-)

      result.as.type <- "AF"

      if(strand == "+") {
        exon.list.1 <- strsplit(event.id.list[4], "-")[[1]]
        exon.list.2 <- strsplit(event.id.list[6], "-")[[1]]

        for (tmp.transcript in as.transcript) {
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.2[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, event.id.list[3], exon.list.1[1], event.id.list[5], exon.list.2[1], exon.list.2[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(event.id.list[3], "-", exon.list.1[1]), splicing_exon=paste0(event.id.list[5], "-", exon.list.2[1]), downstream_exon=paste0(exon.list.2[2], "-", last.end)))
        }

        for (tmp.transcript in diff.transcript) {
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.2[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, event.id.list[3], exon.list.1[1], event.id.list[5], exon.list.2[1], exon.list.2[2], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(event.id.list[3], "-", exon.list.1[1]), splicing_exon=paste0(event.id.list[5], "-", exon.list.2[1]), downstream_exon=paste0(exon.list.2[2], "-", last.end)))
        }
      }else if(strand == "-") {
        exon.list.1 <- strsplit(event.id.list[3], "-")[[1]]
        exon.list.2 <- strsplit(event.id.list[5], "-")[[1]]

        for (tmp.transcript in as.transcript) {
          last.end <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, event.id.list[6], exon.list.2[2], event.id.list[4], exon.list.1[2], exon.list.1[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(event.id.list[6], "-", exon.list.2[2]), splicing_exon=paste0(event.id.list[4], "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[1], "-", last.end)))
        }

        for (tmp.transcript in diff.transcript) {
          last.end <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, event.id.list[6], exon.list.2[2], event.id.list[4], exon.list.1[2], exon.list.1[1], last.end, sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(event.id.list[6], "-", exon.list.2[2]), splicing_exon=paste0(event.id.list[4], "-", exon.list.1[2]), downstream_exon=paste0(exon.list.1[1], "-", last.end)))
        }
      }
    }else if(as.type == "AL") {
      #event id => e1-s2:e2:e1-s3:e3(+) , s1:e1-s3:s2:e2-s3(-)
      #as id => s1:e1:s2:e2:s3:e3(+) , e3:s3:e2:s2:e1:s1(-)

      result.as.type <- "AL"

      if(strand == "+") {
        exon.list.1 <- strsplit(event.id.list[3], "-")[[1]]
        exon.list.2 <- strsplit(event.id.list[5], "-")[[1]]

        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[1], exon.list.1[2], event.id.list[4], exon.list.2[2], event.id.list[6], sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[2], "-", event.id.list[4]), downstream_exon=paste0(exon.list.2[2], "-", event.id.list[6])))
        }

        for (tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$end == exon.list.1[1]), "start"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.1[1], exon.list.1[2], event.id.list[4], exon.list.2[2], event.id.list[6], sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.1[1]), splicing_exon=paste0(exon.list.1[2], "-", event.id.list[4]), downstream_exon=paste0(exon.list.2[2], "-", event.id.list[6])))
        }
      }else if(strand == "-") {
        exon.list.1 <- strsplit(event.id.list[4], "-")[[1]]
        exon.list.2 <- strsplit(event.id.list[6], "-")[[1]]

        for (tmp.transcript in as.transcript) {
          first.start <- as.numeric(as.transcript.gtf[which(as.transcript.gtf$seqid == chrom & as.transcript.gtf$transcript_id == tmp.transcript & as.transcript.gtf$start == exon.list.2[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.2[2], exon.list.2[1], event.id.list[5], exon.list.1[1], event.id.list[3], sep=":"), transcript_id=tmp.transcript, in_or_out="in", upstream_exon=paste0(first.start, "-", exon.list.2[2]), splicing_exon=paste0(exon.list.2[1], "-", event.id.list[5]), downstream_exon=paste0(exon.list.1[1], "-", event.id.list[3])))
        }

        for (tmp.transcript in diff.transcript) {
          first.start <- as.numeric(diff.transcript.gtf[which(diff.transcript.gtf$seqid == chrom & diff.transcript.gtf$transcript_id == tmp.transcript & diff.transcript.gtf$start == exon.list.2[2]), "end"])

          as.id.list <- rbind(as.id.list, data.frame(as_id=paste(chrom, first.start, exon.list.2[2], exon.list.2[1], event.id.list[5], exon.list.1[1], event.id.list[3], sep=":"), transcript_id=tmp.transcript, in_or_out="out", upstream_exon=paste0(first.start, "-", exon.list.2[2]), splicing_exon=paste0(exon.list.2[1], "-", event.id.list[5]), downstream_exon=paste0(exon.list.1[1], "-", event.id.list[3])))
        }
      }
    }

    as.id <- unique(as.id.list$as_id)

    for (tmp.as.id in as.id) {
      tmp.in.as.id <- as.id.list[(as.id.list$as_id == tmp.as.id & as.id.list$in_or_out == "in"),2]
      tmp.out.as.id <- as.id.list[(as.id.list$as_id == tmp.as.id & as.id.list$in_or_out == "out"),2]
      tmp.upstream.exon <- unique(as.id.list[as.id.list$as_id == tmp.as.id, "upstream_exon"])
      tmp.splicing.exon <- unique(as.id.list[as.id.list$as_id == tmp.as.id, "splicing_exon"])
      tmp.downstream.exon <- unique(as.id.list[as.id.list$as_id == tmp.as.id, "downstream_exon"])

      if (length(tmp.in.as.id) >= 1 & length(tmp.out.as.id) >= 1) {
        if (length(tmp.in.as.id) == 1 & length(tmp.out.as.id) > 1) {
          tid.set <- paste(tmp.in.as.id, paste(tmp.out.as.id, collapse="/"), sep=",")
        } else if (length(tmp.in.as.id) > 1 & length(tmp.out.as.id) == 1) {
          tid.set <- paste(paste(tmp.in.as.id, collapse="/"), tmp.out.as.id, sep=",")
        } else {
          tid.set <- paste(paste(tmp.in.as.id, collapse="/"), paste(tmp.out.as.id, collapse="/"), sep=",")
        }

        if (length(tmp.upstream.exon) > 1 | length(tmp.splicing.exon) > 1 | length(tmp.downstream.exon) > 1) {
          print(tmp.as.id)
        }

        as.id.result <- rbind(as.id.result, data.frame(as_id = tmp.as.id, tid_set = tid.set, upstream_exon=tmp.upstream.exon, splicing_exon=tmp.splicing.exon, downstream_exon=tmp.downstream.exon))
      } #else if (length(tmp.in.as.id) == 0)
      #		{
      #			tid.set <- paste(tmp.out.as.id, collapse="/")
      #			as.id.result <- rbind(as.id.result, data.frame(as_id = tmp.as.id, tid_set = tid.set))
      #		} else if (length(tmp.out.as.id) == 0)
      #		{
      #			tid.set <- paste(tmp.in.as.id, collapse="/")
      #			as.id.result <- rbind(as.id.result, data.frame(as_id = tmp.as.id, tid_set = tid.set))
      #		}
    }

    if (length(as.id.result) == 0) {
      next
    }

    #"chr", "strand", "as_type", "upstream_exon", "splicing exon", "downstream exon", "gene_id", "gene_name", "P-value", "dPSI", "as_id"
    for (as.index in 1:length(as.id.result$as_id))
    {
      tmp.as <- as.id.result[as.index,]
      #print(tmp.as)

      #as.dump.result <- rbind(as.dump.result, data.frame(gene_model=gene.model, genome_version=genome.version, chr=chrom, gene_id=gene.id, strand=strand, as_type=result.as.type, as_id=tmp.as$as_id, tid_set=tmp.as$tid_set, event_id=suppa.event.id))
      #as.dump.result <- rbind(as.dump.result, data.frame(gene_model=gene.model, chr=chrom, gene_id=gene.id, strand=strand, as_type=result.as.type, as_id=tmp.as$as_id, tid_set=tmp.as$tid_set, event_id=suppa.event.id, gene_name=gene.name, pvalue=tmp.pvalue, dPSI=tmp.dpsi))
      as.dump.result <- rbind(as.dump.result, data.frame(chr=chrom, strand=strand, as_type=result.as.type, upstream_exon=tmp.as$upstream_exon, splicing_exon=tmp.as$splicing_exon, downstream_exon=tmp.as$downstream_exon, gene_name=gene.name, pvalue=tmp.pvalue, dPSI=tmp.dpsi, as_id=tmp.as$as_id, gene_id=gene.id, tid_set=tmp.as$tid_set, suppa_id=suppa.event.id))
    }

    #print(as.id)
    #print(i)
  }

  #result.file.path <- paste("/comm/home/bond0709/project/06.ASDB.update/03.AS.ID/", gene.model, ".", genome.version, "/", gene.model, ".", genome.version, ".", result.as.type, ".asid.txt", sep="")
  #write.table(as.dump.result, result.file.name, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  #as.dump.result <- as.dump.result[, 1:9]
  return(as.dump.result)
}
