#' converting result mapping to ASDB
#' 
#' DNA and RNA annotation related to aliternative splicing from ASpedia DB(ASDB) add to asr_converter (rMATS_converter, SUPPA_converter, or spliceR_converter) result. And gene enrichment test result are provided between annotation result gene list and 
#' knowledge-based database gene list
#'
#' @param converter.result
#' asr_converter(rMATS_converter, SUPPA_converter, or spliceR_converter) result
#' @param gene.model
#' gene model of reference. One of Refseq, Ensembl, or GENCODE.
#' @param genome.version
#' genome version of reference. One of hg18, GRCh19, or GRCh38.
#' @param gsea.gene.list
#' optional. reference gene list for gene enrichment test with annotation result gene list and knowledge-based database gene list. If gene list is empty, use all genes in reference.
#' @param result.dir
#' directory where annotation result(.tsv file) is saved
#'
#' @return annotation result
#' @export
#'
#' @examples
#' annotation.result.dir <- system.file(“extdata/annotation_result”, package=“ASpediaR”)
#' annotation.result <- asr_annotation(rmats.converter.result, “Ensembl”, “GRCh38”,
#'                                       result.dir=annotation.result.dir)

asr_annotation <- function(converter.result, gene.model="Ensembl", genome.version="GRCh38", gsea.gene.list="", result.dir="") {
  if(result.dir == "") {
    print("*** ERROR MESSAGE: No such output directory.")
    return()
  }
  
  loaded.packages <- tolower((.packages()))
  
  if(("RSQLite" %in% loaded.packages) == FALSE) {
    library(RSQLite)
  }
  
  if(("stringr" %in% loaded.packages) == FALSE) {
    library(stringr)
  }
  
  if(file.exists(data.dir) == FALSE) {
    dir.create(data.dir)
  }
  
  if(!file.exists(result.dir)) {
    dir.create(result.dir)
  }
  
  data.dir <- paste0(.libPaths()[1], "/ASpediaR/data")
  
  gene.list.file.name <- paste0(data.dir, "/", gene.model, ".", genome.version, ".gene.txt")

  if(gsea.gene.list == "") {
    if(!file.exists(gene.list.file.name)) {
      url =  paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/gene_list/", gene.model, ".", genome.version, ".gene.txt")
      download.file(url, gene.list.file.name, method="auto")
    }
    
    gene.list <- read.table(gene.list.file.name, stringsAsFactors=FALSE, header=FALSE)
    whole.gene <- unique(gene.list$V1)
  } else {
    whole.gene <- gava.gene.list
  }
  
  db.file.name <- paste0(data.dir, "/", gene.model, "_", genome.version, ".sqlite")
  
  if(!file.exists(db.file.name)) {
    url =  paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/sqlite/", gene.model, "_", genome.version, ".sqlite")
    download.file(url, db.file.name, method="auto", mode="wb", header=options(timeout=6000))
  }
  
  db.connection <- dbConnect(SQLite(), dbname=db.file.name)

  as.result <- data.frame()

  a3ss.list <- converter.result[converter.result$as_type == "A3SS", ]
  a5ss.list <- converter.result[converter.result$as_type == "A5SS", ]
  se.list <- converter.result[converter.result$as_type == "SE", ]
  mxe.list <- converter.result[converter.result$as_type == "MXE", ]
  ri.list <- converter.result[converter.result$as_type == "RI", ]
  af.list <- converter.result[converter.result$as_type == "AF", ]
  al.list <- converter.result[converter.result$as_type == "AL", ]

  if(nrow(a3ss.list) > 0) {
    upstream.exon <- as.data.frame(matrix(unlist(str_split(a3ss.list$upstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(upstream.exon) <- c("exon_start", "exon_end")

    splicing.exon <- as.data.frame(matrix(unlist(str_split(a3ss.list$splicing_exon, "-")), ncol=2, byrow=TRUE))
    colnames(splicing.exon) <- c("exon_start", "exon_end")

    downstream.exon <- as.data.frame(matrix(unlist(str_split(a3ss.list$downstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(downstream.exon) <- c("exon_start", "exon_end")

    as.id.list <- paste0(a3ss.list$chr, ":", upstream.exon$exon_start, ":", upstream.exon$exon_end, ":", splicing.exon$exon_start, ":", downstream.exon$exon_start, ":", downstream.exon$exon_end)
    query.as.id <- str_replace_all(toString(as.id.list), ", ", "', '")
    query.as.id <- paste0("'", query.as.id, "'")

    as.query <- paste0("select * from as_data where as_id in (", query.as.id, ") and as_type='A3SS'")
    as.query.result <- dbGetQuery(db.connection, as.query)

    as.result <- rbind(as.result, as.query.result)
  }

  if(nrow(a5ss.list) > 0) {
    upstream.exon <- as.data.frame(matrix(unlist(str_split(a5ss.list$upstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(upstream.exon) <- c("exon_start", "exon_end")

    splicing.exon <- as.data.frame(matrix(unlist(str_split(a5ss.list$splicing_exon, "-")), ncol=2, byrow=TRUE))
    colnames(splicing.exon) <- c("exon_start", "exon_end")

    downstream.exon <- as.data.frame(matrix(unlist(str_split(a5ss.list$downstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(downstream.exon) <- c("exon_start", "exon_end")

    as.id.list <- paste0(a5ss.list$chr, ":", splicing.exon$exon_start, ":", upstream.exon$exon_end, ":", splicing.exon$exon_end, ":", downstream.exon$exon_start, ":", downstream.exon$exon_end)
    query.as.id <- str_replace_all(toString(as.id.list), ", ", "', '")
    query.as.id <- paste0("'", query.as.id, "'")

    as.query <- paste0("select * from as_data where as_id in (", query.as.id, ") and as_type='A5SS'")
    as.query.result <- dbGetQuery(db.connection, as.query)

    as.result <- rbind(as.result, as.query.result)
  }

  if(nrow(se.list) > 0) {
    upstream.exon <- as.data.frame(matrix(unlist(str_split(se.list$upstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(upstream.exon) <- c("exon_start", "exon_end")

    splicing.exon <- as.data.frame(matrix(unlist(str_split(se.list$splicing_exon, "-")), ncol=2, byrow=TRUE))
    colnames(splicing.exon) <- c("exon_start", "exon_end")

    downstream.exon <- as.data.frame(matrix(unlist(str_split(se.list$downstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(downstream.exon) <- c("exon_start", "exon_end")

    as.id.list <- paste0(se.list$chr, ":", upstream.exon$exon_start, ":", upstream.exon$exon_end, ":", splicing.exon$exon_start, ":", splicing.exon$exon_end, ":", downstream.exon$exon_start, ":", downstream.exon$exon_end)
    query.as.id <- str_replace_all(toString(as.id.list), ", ", "', '")
    query.as.id <- paste0("'", query.as.id, "'")

    as.query <- paste0("select * from as_data where as_id in (", query.as.id, ") and as_type='SE'")
    as.query.result <- dbGetQuery(db.connection, as.query)

    as.result <- rbind(as.result, as.query.result)
  }

  if(nrow(mxe.list) > 0) {
    upstream.exon <- as.data.frame(matrix(unlist(str_split(mxe.list$upstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(upstream.exon) <- c("exon_start", "exon_end")

    splicing.exon <- as.data.frame(matrix(unlist(str_split(mxe.list$splicing_exon, ";")), ncol=2, byrow=TRUE))
    colnames(splicing.exon) <- c("exon_1", "exon_2")

    splicing.exon.1 <- as.data.frame(matrix(unlist(str_split(splicing.exon$exon_1, "-")), ncol=2, byrow=TRUE))
    colnames(splicing.exon.1) <- c("exon_start", "exon_end")

    splicing.exon.2 <- as.data.frame(matrix(unlist(str_split(splicing.exon$exon_2, "-")), ncol=2, byrow=TRUE))
    colnames(splicing.exon.2) <- c("exon_start", "exon_end")

    downstream.exon <- as.data.frame(matrix(unlist(str_split(mxe.list$downstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(downstream.exon) <- c("exon_start", "exon_end")

    as.id.list <- paste0(mxe.list$chr, ":", upstream.exon$exon_start, ":", upstream.exon$exon_end, ":", splicing.exon.1$exon_start, ":", splicing.exon.1$exon_end, ":", splicing.exon.2$exon_start, ":", splicing.exon.2$exon_end, ":", downstream.exon$exon_start, ":", downstream.exon$exon_end)
    query.as.id <- str_replace_all(toString(as.id.list), ", ", "', '")
    query.as.id <- paste0("'", query.as.id, "'")

    as.query <- paste0("select * from as_data where as_id in (", query.as.id, ") and as_type='MXE'")
    as.query.result <- dbGetQuery(db.connection, as.query)

    as.result <- rbind(as.result, as.query.result)
  }

  if(nrow(ri.list) > 0) {
    upstream.exon <- as.data.frame(matrix(unlist(str_split(ri.list$upstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(upstream.exon) <- c("exon_start", "exon_end")

    downstream.exon <- as.data.frame(matrix(unlist(str_split(ri.list$downstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(downstream.exon) <- c("exon_start", "exon_end")

    as.id.list <- paste0(ri.list$chr, ":", upstream.exon$exon_start, ":", upstream.exon$exon_end, ":", downstream.exon$exon_start, ":", downstream.exon$exon_end)
    query.as.id <- str_replace_all(toString(as.id.list), ", ", "', '")
    query.as.id <- paste0("'", query.as.id, "'")

    as.query <- paste0("select * from as_data where as_id in (", query.as.id, ") and as_type='RI'")
    as.query.result <- dbGetQuery(db.connection, as.query)

    as.result <- rbind(as.result, as.query.result)
  }

  if(nrow(af.list) > 0) {
    upstream.exon <- as.data.frame(matrix(unlist(str_split(af.list$upstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(upstream.exon) <- c("exon_start", "exon_end")

    splicing.exon <- as.data.frame(matrix(unlist(str_split(af.list$splicing_exon, "-")), ncol=2, byrow=TRUE))
    colnames(splicing.exon) <- c("exon_start", "exon_end")

    downstream.exon <- as.data.frame(matrix(unlist(str_split(af.list$downstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(downstream.exon) <- c("exon_start", "exon_end")

    as.id.list <- paste0(af.list$chr, ":", upstream.exon$exon_start, ":", upstream.exon$exon_end, ":", splicing.exon$exon_start, ":", splicing.exon$exon_end, ":", downstream.exon$exon_start, ":", downstream.exon$exon_end)
    query.as.id <- str_replace_all(toString(as.id.list), ", ", "', '")
    query.as.id <- paste0("'", query.as.id, "'")

    as.query <- paste0("select * from as_data where as_id in (", query.as.id, ") and as_type='AF'")
    as.query.result <- dbGetQuery(db.connection, as.query)

    as.result <- rbind(as.result, as.query.result)
  }

  if(nrow(al.list) > 0) {
    upstream.exon <- as.data.frame(matrix(unlist(str_split(al.list$upstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(upstream.exon) <- c("exon_start", "exon_end")

    splicing.exon <- as.data.frame(matrix(unlist(str_split(al.list$splicing_exon, "-")), ncol=2, byrow=TRUE))
    colnames(splicing.exon) <- c("exon_start", "exon_end")

    downstream.exon <- as.data.frame(matrix(unlist(str_split(al.list$downstream_exon, "-")), ncol=2, byrow=TRUE))
    colnames(downstream.exon) <- c("exon_start", "exon_end")

    as.id.list <- paste0(al.list$chr, ":", upstream.exon$exon_start, ":", upstream.exon$exon_end, ":", splicing.exon$exon_start, ":", splicing.exon$exon_end, ":", downstream.exon$exon_start, ":", downstream.exon$exon_end)
    query.as.id <- str_replace_all(toString(as.id.list), ", ", "', '")
    query.as.id <- paste0("'", query.as.id, "'")

    as.query <- paste0("select * from as_data where as_id in (", query.as.id, ") and as_type='AL'")
    as.query.result <- dbGetQuery(db.connection, as.query)

    as.result <- rbind(as.result, as.query.result)
  }

  annotation.result <- data.frame()

  if(nrow(as.result) > 0) {
    for(i in 1:nrow(as.result)) {
      as.raw <- as.result[i, ]
      gene.symbol <- as.raw$gene_symbol
      chr <- as.raw$chr
      as.id <- as.raw$as_id
      as.rename <- as.raw$as_rename
      as.type <- as.raw$as_type
      strand <- as.raw$strand
      gene.name <- as.raw$gene_name
      locus.group <- as.raw$locus_group
      location <- as.raw$location
      gene.id <- as.raw$gene_id

      transcript.id <- ""
      in.transcript.id <- ""
      ex.transcript.id <- ""

      if(as.raw$tid_set != "") {
        transcript.id <- str_replace_all(as.raw$tid_set, "/", ",")
        in.transcript.id <- str_replace_all(str_split(as.raw$tid_set, ",")[[1]][1], "/", ",")
        ex.transcript.id <- str_replace_all(str_split(as.raw$tid_set, ",")[[1]][2], "/", ",")
      }

      go.bp <- as.raw$go_bp
      go.cc <- as.raw$go_cc
      go.mf <- as.raw$go_mf
      conservation <- as.raw$conservation
      variant.dbsnp <- as.raw$variant_dbsnp
      variant.cosmic <- as.raw$variant_cosmic
      variant.spidex <- as.raw$variant_spidex
      mirna <- as.raw$mirna
      repeats <- as.raw$repeats
      nmd.stop <- as.raw$nmd_stop
      nmd.cosmic <- as.raw$nmd_cosmic
      nmd.dbsnp <- as.raw$nmd_dbsnp
      protein.domain <- as.raw$protein_domain
      ptm <- as.raw$ptm
      rbp.sf <- as.raw$rbp_sf
      rbp <- as.raw$rbp_not_sf
      trans <- as.raw$trans
      isoform_PPI_a <- ""
      isoform_PPI_b <- ""

      if(as.raw$ppi != "") {
        tmp.ppi <- str_split(as.raw$ppi, ";")[[1]]
        tmp.ppi.list <- str_split(tmp.ppi, ",")
        ppi.transcript.id.list <- c()
        ppi.partner.list <- c()

        for(k in 1:length(tmp.ppi.list)) {
          ppi.transcript.id.list <- c(ppi.transcript.id.list, tmp.ppi.list[[k]][1])
          ppi.partner.list <- c(ppi.partner.list, str_replace_all(tmp.ppi.list[[k]][4], "<>", "/"))
        }

        isoform_PPI_a <- str_replace_all(toString(ppi.transcript.id.list), ", ", ";")
        isoform_PPI_b <- str_replace_all(toString(ppi.partner.list), ", ", ";")
      }

      isoform_subcellular_localization_id <- as.raw$isoform_subcellular_localization_id
      isoform_subcellular_localization <- as.raw$isoform_subcellular_localization

      tmp.result <- c(gene.symbol, chr, as.id, as.rename, as.type, strand, gene.name, locus.group, location, gene.id, transcript.id, in.transcript.id, ex.transcript.id, go.bp, go.cc, go.mf, conservation, variant.dbsnp, variant.cosmic, variant.spidex, mirna, repeats, nmd.stop, nmd.cosmic, nmd.dbsnp, protein.domain, ptm, rbp.sf, rbp, isoform_PPI_a, isoform_PPI_b, isoform_subcellular_localization_id, isoform_subcellular_localization, trans)

      annotation.result <- rbind(annotation.result, tmp.result)
    }

    colnames(annotation.result) <- c("gene_symbol", "chr", "as_id", "as_description_id", "as_type", "strand", "gene_name", "locus_group", "location", "gene_id", "transcript_id", "exon_inclusion_transcript_id", "exon_exclusion_transcript_id", "GO_BP", "GO_CC", "GO_MF", "conservation_score", "dbSNP_variant", "COSMIC_variant", "SPIDEX_variant", "miRNA_binding_site", "repeat", "NMD", "COSMIC_NMD", "dbSNP_NMD", "protein_domain", "protein_translational_modification", "RBP_splicing_factor", "RBP", "isoform_PPI_a", "isoform_PPI_b", "isoform_subcellular_localization_id", "isoform_subcellular_localization", "transmembrane")
  }
  
  #GSEA
  mining_gsea(unique(annotation.result$gene_symbol), whole.gene, result.dir)
  
  dbDisconnect(db.connection)
  return(annotation.result)
}