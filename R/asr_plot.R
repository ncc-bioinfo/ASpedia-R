#' visualization of annotation result
#'
#' @param annotation.result
#' asr_annotation function result.
#' @param gtf.file.name
#' a GTF format file of reference.
#' @param gene.model
#' gene model of reference. One of Refseq, Ensembl, or GENCODE.
#' @param genome.version
#' genome version of reference. One of hg18, GRCh19, or GRCh38.
#' @param gene.name
#' gene name to be visualization. 
#' @param as.id
#' list of AS ID to be visualization.
#' @param heights.list
#' positive integer vectors for track heights include height of ideogram and gene track.
#' @param plot.data.list 
#' list of DNA or RNA feature to be visualization. choose from “conservation”, “NMD”, “repeats”, “domain”, “PTM”, or “RBP”.
#' @param result.dir
#' directory where plots(.png files) are saved
#'
#' @return 
#' @export
#'
#' @examples
#' gtf.file.name <- system.file(“extdata”, “test_gtf.gtf”, package=“ASpediaR”)
#' plot.result.dir <- system.file(“extdata/plot_result”, package=“ASpediaR”)
#' 
#' ##using gene name
#' asr_plot(annotation.result, gtf.file.name, gene.model="Ensembl", genome.version="GRCh38",
#'           gene.name="FGFR2", result.dir=plot.result.dir)
#'
#' ##using AS ID
#' asr_plot(annotation.result, gtf.file.name, gene.model="Ensembl", genome.version="GRCh38",
#'           as.id="chr10:121520169:121519979:121518829:121518682:121517463:121517319",
#'           result.dir=plot.result.dir)
#' 
#' ##track lists and track heights are change
#' asr_plot(annotation.result, gtf.file.name, gene.model="Ensembl", genome.version="GRCh38",
#'           gene.name="FGFR2", heights.list=c(1, 2, 2, 1, 1, 1, 2),
#'           list.of.plot=c(“conservation”, "domain", "PTM", "repeats", “RBP”), result.dir=plot.result.dir)


asr_plot <- function(annotation.result, gtf.file.name, gene.model="Ensembl", genome.version="hg38", gene.name="", as.id="", heights.list="", plot.data.list="", result.dir="") {
  if(gene.name == "" || as.id == "") {
    print("*** ERROR MESSAGE: Input gene name or AS ID is empty.")
    return()
  }
  
  loaded.packages <- tolower((.packages()))
  
  if(("ggbio" %in% loaded.packages) == FALSE) {
    library(ggbio)
  }
  
  if(("RSQLite" %in% loaded.packages) == FALSE) {
    library(RSQLite)
  }
  
  if(("rtracklayer" %in% loaded.packages) == FALSE) {
    library(rtracklayer)
  }
  
  if(("stringr" %in% loaded.packages) == FALSE) {
    library(stringr)
  }
  
  if(("ggraph" %in% loaded.packages) == FALSE) {
    library(ggraph)
  }
  
  if(("igraph" %in% loaded.packages) == FALSE) {
    library(igraph)
  }
  
  if(!file.exists(gene.list.file.name)) {
    url =  paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/gene_list/", gene.model, ".", genome.version, ".gene.txt")
    download.file(url, gene.list.file.name, method="auto")
  }
  
  db.file.name <- paste0("data/", gene.model, "_", genome.version, ".sqlite")
  db.connection <- dbConnect(SQLite(), dbname=db.file.name)
  
  #annotation result filtered by gene name or as ID
  if(gene.name != "") {
    if(length(gene.name == 1)) {
      annotation.result <- annotation.result[annotation.result$gene_symbol == gene.name, ]
    } else {
      annotation.result <- annotation.result[annotation.result$gene_symbol %in% gene.name, ]
    }
  } else if(as.id != "") {
    if(length(as.id) == 1) {
      annotation.result <- annotation.result[annotation.result$as_id == as.id] 
    } else {
      annotation.result <- annotation.result[annotation.result$as_id %in% as.id]
    }
  }
  
  #check plot list
  has.conservation <- FALSE
  has.repeats <- FALSE
  has.nmd <- FALSE
  has.domain <- FALSE
  has.ptm <- FALSE
  has.rbp <- FALSE
  
  if(plot.data.list == "") {
    has.conservation <- TRUE
    has.repeats <- TRUE
    has.nmd <- TRUE
    has.domain <- TRUE
    has.ptm <- TRUE
    has.rbp <- TRUE
  } else {
    if("conservation" %in% tolower(plot.data.list)) {
      has.conservation <- TRUE
    }
    
    if("repeats" %in% tolower(plot.data.list)) {
      has.repeats <- TRUE
    }
    
    if("nmd" %in% tolower(plot.data.list)) {
      has.nmd <- TRUE
    }
    
    if("domain" %in% tolower(plot.data.list)) {
      has.domain <- TRUE
    }
    
    if("ptm" %in% tolower(plot.data.list)) {
      has.ptm <- TRUE
    }
    
    if("rbp" %in% tolower(plot.data.list)) {
      has.rbp <- TRUE
    }
  }
  
  gtf.data <- import(gtf.file.name)

  #ideogram
  ideo.plot <- Ideogram(genome=genome.version)
  
  for(anno.index in 1:nrow(annotation.result)) {
    tmp.anno <- annotation.result[anno.index, ]

    as.id <- tmp.anno$as_id
    as.id.split <- str_split(as.id, ":")[[1]]
    as.chr <- tmp.anno$chr
    as.type <- tmp.anno$as_type
    as.strand <- tmp.anno$strand
    gene.name <- tmp.anno$gene_symbol
    gene.id <- tmp.anno$gene_id
    exon.region <- ""
    as.region <- ""
  
    if(as.strand == "+") {
      exon.region <- c(as.numeric(as.id.split[2]), as.numeric(as.id.split[length(as.id.split)]))
    } else if(as.strand == "-") {
      exon.region <- c(as.numeric(as.id.split[length(as.id.split)]), as.numeric(as.id.split[2]))
    }
  
    #ideogram
    ideo.plot@subchr <- as.chr
    ideo.plot@zoom.region <- exon.region
  
    #gene plot
    transcript.id.list <- c(str_split(tmp.anno$exon_inclusion_transcript_id, ",")[[1]], str_split(tmp.anno$exon_exclusion_transcript_id, ",")[[1]])
    gene.data <- gtf.data[gtf.data$transcript_id %in% transcript.id.list & gtf.data$type == "exon", c("type", "gene_id", "gene_name", "transcript_id", "transcript_biotype", "exon_number", "exon_id", "protein_id")]
    gene.data$exon_type <- "exon"
  
    if(as.type == "A3SS") {
      if(as.strand == "+") {
        as.exon <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[3]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[6]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[5]), as.numeric(as.id.split[6]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[6]))
        as.exon.region <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[3]))
      } else if (as.strand == "-") {
        as.exon <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[2]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[4]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[5]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[4]))
        as.exon.region <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[2]))
      }
    } else if(as.type == "A5SS") {
      if(as.strand == "+") {
        as.exon <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[3]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[4]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[5]), as.numeric(as.id.split[6]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[4]))
        as.exon.region <- IRanges(as.numeric(as.id.split[5]), as.numeric(as.id.split[6]))
      } else if (as.strand == "-") {
        as.exon <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[2]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[2]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[5]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[2]))
        as.exon.region <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[5]))
      }
    } else if(as.type == "SE") {
      if(as.strand == "+") {
        as.exon <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[3]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[5]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[7]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[5]))
        as.exon.region <- IRanges(start=c(as.numeric(as.id.split[2]), as.numeric(as.id.split[6])), end=c(as.numeric(as.id.split[3]), as.numeric(as.id.split[7])))
      } else if (as.strand == "-") {
        as.exon <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[2]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[5]), as.numeric(as.id.split[4]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[7]), as.numeric(as.id.split[6]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(as.numeric(as.id.split[5]), as.numeric(as.id.split[4]))
        as.exon.region <- IRanges(start=c(as.numeric(as.id.split[7]), as.numeric(as.id.split[3])), end=c(as.numeric(as.id.split[6]), as.numeric(as.id.split[2])))
      }
    } else if(as.type == "MXE") {
      if(as.strand == "+") {
        as.exon <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[3]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[5]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[7]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[8]), as.numeric(as.id.split[9]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(start=c(as.numeric(as.id.split[4]), as.numeric(as.id.split[6])), end=c(as.numeric(as.id.split[5]), as.numeric(as.id.split[7])))
        as.exon.region <- IRanges(start=c(as.numeric(as.id.split[2]), as.numeric(as.id.split[8])), end=c(as.numeric(as.id.split[3]), as.numeric(as.id.split[9])))
      } else if (as.strand == "-") {
        as.exon <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[2]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[5]), as.numeric(as.id.split[4]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[7]), as.numeric(as.id.split[6]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[9]), as.numeric(as.id.split[8]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(start=c(as.numeric(as.id.split[7]), as.numeric(as.id.split[5])), end=c(as.numeric(as.id.split[6]), as.numeric(as.id.split[4])))
        as.exon.region <- IRanges(start=c(as.numeric(as.id.split[9]), as.numeric(as.id.split[3])), end=c(as.numeric(as.id.split[8]), as.numeric(as.id.split[2])))
      }
    } else if(as.type == "RI") {
      if(as.strand == "+") {
        as.exon <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[3]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[5]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[4]))
        as.exon.region <- IRanges(start=c(as.numeric(as.id.split[2]), as.numeric(as.id.split[4])), end=c(as.numeric(as.id.split[3]), as.numeric(as.id.split[5])))
      } else if (as.strand == "-") {
        as.exon <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[2]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[5]), as.numeric(as.id.split[4]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[3]))
        as.exon.region <- IRanges(start=c(as.numeric(as.id.split[5]), as.numeric(as.id.split[3])), end=c(as.numeric(as.id.split[4]), as.numeric(as.id.split[2])))
      }
    } else if(as.type == "AF") {
      if(as.strand == "+") {
        as.exon <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[3]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[5]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[7]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(start=c(as.numeric(as.id.split[2]), as.numeric(as.id.split[4])), end=c(as.numeric(as.id.split[3]), as.numeric(as.id.split[5])))
        as.exon.region <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[7]))
      } else if (as.strand == "-") {
        as.exon <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[2]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[5]), as.numeric(as.id.split[4]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[7]), as.numeric(as.id.split[6]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(start=c(as.numeric(as.id.split[5]), as.numeric(as.id.split[3])), end=c(as.numeric(as.id.split[4]), as.numeric(as.id.split[2])))
        as.exon.region <- IRanges(as.numeric(as.id.split[7]), as.numeric(as.id.split[6]))
      }
    } else if(as.type == "AL") {
      if(as.strand == "+") {
        as.exon <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[3]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[4]), as.numeric(as.id.split[5]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[6]), as.numeric(as.id.split[7]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(start=c(as.numeric(as.id.split[4]), as.numeric(as.id.split[6])), end=c(as.numeric(as.id.split[5]), as.numeric(as.id.split[7])))
        as.exon.region <- IRanges(as.numeric(as.id.split[2]), as.numeric(as.id.split[3]))
      } else if (as.strand == "-") {
        as.exon <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[2]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[5]), as.numeric(as.id.split[4]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.exon <- IRanges(as.numeric(as.id.split[7]), as.numeric(as.id.split[6]))
        
        if(length(gene.data[ranges(gene.data) == as.exon, ]) != 0) {
          gene.data[ranges(gene.data) == as.exon, ]$exon_type <- "as_exon"
        }
  
        as.region <- IRanges(start=c(as.numeric(as.id.split[7]), as.numeric(as.id.split[5])), end=c(as.numeric(as.id.split[6]), as.numeric(as.id.split[4])))
        as.exon.region <- IRanges(as.numeric(as.id.split[3]), as.numeric(as.id.split[2]))
      }
    }
  
    #gene plot
    gene.plot.data <- split(gene.data, gene.data$transcript_id)
    gene.text.position <- mean(exon.region)
    gene.plot <- autoplot(gene.plot.data, label.color="black", aes(fill=exon_type)) + scale_fill_manual(values=c("exon"="#619CFF", "as_exon"="#F8766D")) + ylab(paste0("Gene\n(", gene.id, ",", gene.name, ")")) + xlim(exon.region) + guides(y="none")
    gene.plot@ggplot$layers[[4]]$data$midpoint <- gene.text.position
  
    as.track <- tracks(ideo.plot, gene.plot)
  
    if(length(transcript.id.list) <= 4) {
      as.track@heights <- c(1, 1)
    }else if(length(transcript.id.list) <= 6) {
      as.track@heights <- c(1, 2)
    }else if(length(transcript.id.list) <= 8) {
      as.track@heights <- c(1, 3)
    }else {
      as.track@heights <- c(1, 4)
    }
  
  
    #conservation
    if((tmp.anno$conservation_score != "") && (has.conservation == TRUE)) {
      conservation.query <- paste0("select conservation_raw from as_data where as_id='", as.id , "' and as_type='", as.type, "'")
      conservation.result <- dbGetQuery(db.connection, conservation.query)
      tmp.conservation.data <- as.data.frame(matrix(unlist(str_split(str_split(conservation.result, ";")[[1]], ",")), ncol=5, byrow=TRUE))
      colnames(tmp.conservation.data) <- c("conservation_db", "region_type", "start", "end", "score")
  
      conservation.data <- data.frame()
  
      for(i in 1:nrow(tmp.conservation.data)) {
        tmp.conservation <- tmp.conservation.data[i, ]
        conservation.data <- rbind(conservation.data, c(tmp.conservation$conservation_db, tmp.conservation$start, tmp.conservation$score))
        conservation.data <- rbind(conservation.data, c(tmp.conservation$conservation_db, tmp.conservation$end, tmp.conservation$score))
      }
  
      colnames(conservation.data) <- c("conservation_db", "position", "score")
      conservation.data$position <- as.numeric(conservation.data$position)
      conservation.data$score <- as.numeric(conservation.data$score)
  
      conservation.db.list <- unique(conservation.data$conservation_db)
  
      for(conservation.db in conservation.db.list) {
        tmp.conservation.data <- conservation.data[conservation.data$conservation_db==conservation.db, ]
        conservation.plot <- ggplot(data=tmp.conservation.data, aes(x=position, y=score)) + geom_area(fill="#F8766D") + ylab(conservation.db) + ylim(c(0.0, 1.0))
  
        prev.track.height <- as.track@heights
        as.track <- as.track + tracks(conservation.plot)
        as.track@heights <- c(prev.track.height, 1)
      }
    }
  
    #protein domain
    if((tmp.anno$protein_domain != "") && (has.domain == TRUE)) {
      tmp.anno$protein_domain <- str_replace_all(tmp.anno$protein_domain, ", ", " ")
      tmp.domain.data.list <-str_split(str_split(tmp.anno$protein_domain, ";")[[1]], ",")
      #tmp.domain.data <- as.data.frame(matrix(unlist(str_split(str_split(tmp.anno$protein_domain, ";")[[1]], ",")), ncol=5, byrow=TRUE))
      tmp.domain.data <- data.frame()
      
      for(domain.index in 1:length(tmp.domain.data.list)){
        tmp.domain <- tmp.domain.data.list[[domain.index]]
        tmp.domain.length <- length(tmp.domain)
        
        if(tmp.domain.length > 5) {
          diff.length <- tmp.domain.length - 5
          tmp.domain.desc <- paste(tmp.domain[3:(3 + diff.length)], collapse=",")
          tmp.domain.data <- rbind(tmp.domain.data, c(tmp.domain[1:2], tmp.domain.desc, tmp.domain[(3 + diff.length + 1):tmp.domain.length]))
        } else {
          tmp.domain.data <- rbind(tmp.domain.data, tmp.domain)
        }
      }
        
      domain.position.list <- as.data.frame(matrix(unlist(str_split(tmp.domain.data[, 5], "-")), ncol=2, byrow=TRUE))
      domain.data <- cbind(tmp.domain.data[, 1:4], domain.position.list)
      colnames(domain.data) <- c("transcript_id", "pfam_id", "desc", "protein_position", "start", "end")
  
      domain.data <- domain.data[domain.data$transcript_id %in% transcript.id.list, ]
      
      if(nrow(domain.data) > 0) {
        pfam.id.list <- unique(domain.data$pfam_id)
        domain.plot.data <- data.frame()
        domain.text.data <- data.frame()
    
        for(pfam.id in pfam.id.list) {
          tmp.domain.plot.data <- domain.data[domain.data$pfam_id == pfam.id, ]
          pfam.desc <- tmp.domain.plot.data[1, "desc"]
    
          domain.region <- IRanges(start=as.numeric(tmp.domain.plot.data$start), end=as.numeric(tmp.domain.plot.data$end))
          domain.plot.region <- reduce(domain.region)
    
          for(i in 1:length(domain.plot.region)) {
            tmp.domain.region <- domain.plot.region[i]
    
            domain.plot.data <- rbind(domain.plot.data, c(start(tmp.domain.region), end(tmp.domain.region), (i - 1), (i - 1 + 0.5), "domain"))
            as.inter <- intersect(tmp.domain.region, as.region)
            exon.inter <- intersect(tmp.domain.region, as.exon.region)
    
            if(length(as.inter) > 0) {
              for(j in 1:length(as.inter)) {
                domain.plot.data <- rbind(domain.plot.data, c(start(as.inter[j]), end(as.inter[j]), (i - 1), (i - 1 + 0.5), "as"))
              }
            }
    
            if(length(exon.inter) > 0) {
              for(j in 1:length(exon.inter)) {
                domain.plot.data <- rbind(domain.plot.data, c(start(exon.inter[j]), end(exon.inter[j]), (i - 1), (i - 1 + 0.5), "exon"))
              }
            }
    
            if(start(tmp.domain.region) < exon.region[1]) {
              domain.text.data <- rbind(domain.text.data, c(exon.region[1], (i - 1 + 0.6), pfam.desc))
            } else {
              domain.text.data <- rbind(domain.text.data, c(start(tmp.domain.region), (i - 1 + 0.6), pfam.desc))
            }
          }
        }
    
        colnames(domain.plot.data) <- c("xmin", "xmax", "ymin", "ymax", "domain_type")
        colnames(domain.text.data) <- c("x", "y", "domain_desc")
    
        domain.plot.data$xmin <- as.numeric(domain.plot.data$xmin)
        domain.plot.data$xmax <- as.numeric(domain.plot.data$xmax)
        domain.plot.data$ymin <- as.numeric(domain.plot.data$ymin)
        domain.plot.data$ymax <- as.numeric(domain.plot.data$ymax)
    
        domain.text.data$x <- as.numeric(domain.text.data$x)
        domain.text.data$y <- as.numeric(domain.text.data$y)
    
        domain.plot <- ggplot() + geom_rect(data=domain.plot.data, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=domain_type)) + geom_text(data=domain.text.data, mapping=aes(x=x, y=y, label=domain_desc), hjust=0, vjust=0) + scale_fill_manual(values=c("domain"="#CAC7ED", "exon"="#9590FF", "as"=I("pink"))) + xlim(exon.region) + ylim(0, max(domain.plot.data$ymax) + 0.5) + ylab("Protein domain") + guides(y="none")
        domain.track <- tracks(domain.plot)
    
        prev.track.height <- as.track@heights
        as.track <- as.track + domain.track
        as.track@heights <- c(prev.track.height, 1)
      }
    }
  
    #PTM
    if((tmp.anno$protein_translational_modification != "") && has.ptm == TRUE){
      ptm.data <- as.data.frame(matrix(unlist(str_split(str_split(tmp.anno$protein_translational_modification, ";")[[1]], ",")), ncol=3, byrow=TRUE))
      ptm.data <- cbind(ptm.data[,c(1, 2)], as.data.frame(matrix(unlist(str_split(ptm.data$V3, ":")), ncol=2, byrow=TRUE)))
      colnames(ptm.data) <- c("ptm_type", "protein_position", "chr", "position")
  
      ptm.plot.data <- cbind(ptm.data[, c(1,3)], as.data.frame(matrix(unlist(str_split(ptm.data$position, "-")), ncol=2, byrow=TRUE)))
      colnames(ptm.plot.data) <- c("ptm_type", "chr", "start", "end")
  
      ptm.dot.plot.data <- ptm.plot.data[(as.numeric(ptm.plot.data$end) - as.numeric(ptm.plot.data$start)) <= 10, ]
      ptm.rect.plot.data <- ptm.plot.data[(as.numeric(ptm.plot.data$end) - as.numeric(ptm.plot.data$start)) > 10, ]
  
      ptm.dot.plot.count <- nrow(ptm.dot.plot.data)
      ptm.rect.plot.count <- nrow(ptm.rect.plot.data)
  
      if(ptm.dot.plot.count > 0) {
        ptm.dot.plot.data <- cbind(ptm.dot.plot.data, y=c(1:ptm.dot.plot.count))
        ptm.dot.plot.data$start <- as.numeric(ptm.dot.plot.data$start)
        ptm.dot.plot.data$end <- as.numeric(ptm.dot.plot.data$end)
        ptm.dot.plot.data$y <- as.numeric(ptm.dot.plot.data$y)
        ptm.dot.plot.data$y <- ptm.dot.plot.data$y + 0.5
      }
  
      if(ptm.rect.plot.count > 0) {
        ptm.rect.plot.data <- cbind(ptm.rect.plot.data, ymin=c((ptm.dot.plot.count + 1):(ptm.dot.plot.count + ptm.rect.plot.count)))
        ptm.rect.plot.data <- cbind(ptm.rect.plot.data, ymax=(ptm.rect.plot.data$ymin + 0.5))
        ptm.rect.plot.data$start <- as.numeric(ptm.rect.plot.data$start)
        ptm.rect.plot.data$end <- as.numeric(ptm.rect.plot.data$end)
        ptm.rect.plot.data$ymin <- as.numeric(ptm.rect.plot.data$ymin)
        ptm.rect.plot.data$ymax <- as.numeric(ptm.rect.plot.data$ymax)
  
        ptm.rect.text.data <- data.frame()
  
        for(i in 1:ptm.rect.plot.count) {
          if(ptm.rect.plot.data[i, "start"] < exon.region[1]) {
            ptm.rect.text.data <- rbind(ptm.rect.text.data, c(exon.region[1], ptm.rect.plot.data[i, "ymax"], ptm.rect.plot.data[i, "ptm_type"]))
          } else {
            ptm.rect.text.data <- rbind(ptm.rect.text.data, c(ptm.rect.plot.data[i, "start"], ptm.rect.plot.data[i, "ymax"], ptm.rect.plot.data[i, "ptm_type"]))
          }
        }
  
        colnames(ptm.rect.text.data) <- c("x", "y", "ptm_type")
        ptm.rect.text.data$x <- as.numeric(ptm.rect.text.data$x)
        ptm.rect.text.data$y <- as.numeric(ptm.rect.text.data$y)
      }
  
      ptm.plot <- NULL
  
      if(ptm.dot.plot.count == 0) {
        if(ptm.rect.plot.count > 0) {
          ptm.plot <- ggplot() + geom_rect(data=ptm.rect.plot.data, mapping=aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax), fill=I("pink")) + geom_text(data=ptm.rect.text.data, mapping=aes(x=x, y=(y + 0.1), label=ptm_type), hjust=0, vjust=0) + xlim(exon.region) + guides(y="none") + ylab("Protein Translational Modification") + ylim(1, max(ptm.rect.plot.data$ymax) + 0.5)
        }
      } else {
        if(ptm.rect.plot.count == 0) {
        ptm.plot <- ggplot() + geom_dotplot(data=ptm.dot.plot.data, mapping=aes(x=start, y=y), binaxis="y", stackdir="center", dotsize=2, fill="#9590FF") + geom_text(data=ptm.dot.plot.data, mapping=aes(x=start, y=y, label=ptm_type), hjust=0, vjust=0) + xlim(exon.region) + guides(y="none") + ylab("Protein Translational Modification") + ylim(1, max(ptm.dot.plot.data$y) + 0.5)
        } else {
          ptm.plot <- ggplot() + geom_dotplot(data=ptm.dot.plot.data, mapping=aes(x=start, y=y), binaxis="y", stackdir="center", dotsize=2, fill="#9590FF") + geom_text(data=ptm.dot.plot.data, mapping=aes(x=start, y=y, label=ptm_type), hjust=0, vjust=0) + geom_rect(data=ptm.rect.plot.data, mapping=aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax), fill=I("pink")) + geom_text(data=ptm.rect.text.data, mapping=aes(x=x, y=(y + 0.1), label=ptm_type), hjust=0, vjust=0) + xlim(exon.region) + guides(y="none") + ylab("Protein Translational Modification") + ylim(1, max(ptm.rect.plot.data$ymax) + 0.5)
        }
      }
  
      if(is.null(ptm.plot) == FALSE) {
        ptm.track <- tracks(ptm.plot)
  
        prev.track.height <- as.track@heights
        as.track <- as.track + ptm.track
        as.track@heights <- c(prev.track.height, 1)
      }
    }
    
    #NMD
    nmd.stop.plot <- ""
    nmd.cosmic.plot <- ""
    nmd.dbsnp.plot <- ""
    
    if(tmp.anno$NMD != "" && has.nmd == TRUE) {
      nmd.stop <- as.data.frame(matrix(unlist(str_split(str_split(tmp.anno$NMD, ";")[[1]], ":")), ncol=2, byrow=TRUE))
      colnames(nmd.stop) <- c("chr", "position")
      nmd.stop <- unique(nmd.stop)
      nmd.stop.y.position <- 1:nrow(nmd.stop)
      
      nmd.stop.plot.data <- data.frame(x=as.numeric(nmd.stop$position), y=as.numeric(nmd.stop.y.position))
      
      nmd.stop.plot <- ggplot() + geom_dotplot(data=nmd.stop.plot.data, mapping=aes(x=x, y=y), binaxis="y", stackdir="center", dotsize=2, fill="#9590FF") + xlim(exon.region) + guides(y="none") + ylab("NMD COSMIC") + ylim(1, max(nmd.stop.plot.data$y) + 0.5)
    }
    
    if(tmp.anno$COSMIC_NMD && has.nmd == TRUE) {
      nmd.cosmic <- as.data.frame(matrix(unlist(str_split(str_split(tmp.anno$COSMIC_NMD, ";")[[1]], ":")), ncol=2, byrow=TRUE))
      colnames(nmd.cosmic) <- c("chr", "position")
      nmd.cosmic <- unique(nmd.cosmic)
      nmd.cosmic.y.position <- 1:nrow(nmd.cosmic)
      
      nmd.cosmic.plot.data <- data.frame(x=as.numeric(nmd.cosmic$position), y=as.numeric(nmd.cosmic.y.position))
      
      nmd.cosmic.plot <- ggplot() + geom_dotplot(data=nmd.cosmic.plot.data, mapping=aes(x=x, y=y), binaxis="y", stackdir="center", dotsize=2, fill="#9590FF") + xlim(exon.region) + guides(y="none") + ylab("NMD COSMIC") + ylim(1, max(nmd.cosmic.plot.data$y) + 0.5)
    }
    
    if(tmp.anno$dbSNP_NMD && has.nmd == TRUE) {
      nmd.dbsnp <- as.data.frame(matrix(unlist(str_split(str_split(tmp.anno$dbSNP_NMD, ";")[[1]], ":")), ncol=2, byrow=TRUE))
      colnames(nmd.dbsnp) <- c("chr", "position")
      nmd.dbsnp <- unique(nmd.dbsnp)
      nmd.dbsnp.y.position <- 1:nrow(nmd.dbsnp)
      
      nmd.dbsnp.plot.data <- data.frame(x=as.numeric(nmd.dbsnp$position), y=as.numeric(nmd.dbsnp.y.position))
      
      nmd.dbsnp.plot <- ggplot() + geom_dotplot(data=nmd.dbsnp.plot.data, mapping=aes(x=x, y=y), binaxis="y", stackdir="center", dotsize=2, fill="#9590FF") + xlim(exon.region) + guides(y="none") + ylab("NMD COSMIC") + ylim(1, max(nmd.dbsnp.plot.data$y) + 0.5)
    }
    
    #repeats
    if((tmp.anno$repeats != "") && (has.repeats == TRUE)) {
      tmp.anno$repeats <- str_replace_all(tmp.anno$repeats, ", ", " ")
      tmp.repeats.data.list <-str_split(str_split(tmp.anno$repeats, ";")[[1]], ",")
      #tmp.domain.data <- as.data.frame(matrix(unlist(str_split(str_split(tmp.anno$protein_domain, ";")[[1]], ",")), ncol=5, byrow=TRUE))
      tmp.domain.data <- data.frame()
      
      for(domain.index in 1:length(tmp.domain.data.list)){
        tmp.domain <- tmp.domain.data.list[[domain.index]]
        tmp.domain.length <- length(tmp.domain)
        
        if(tmp.domain.length > 5) {
          diff.length <- tmp.domain.length - 5
          tmp.domain.desc <- paste(tmp.domain[3:(3 + diff.length)], collapse=",")
          tmp.domain.data <- rbind(tmp.domain.data, c(tmp.domain[1:2], tmp.domain.desc, tmp.domain[(3 + diff.length + 1):tmp.domain.length]))
        } else {
          tmp.domain.data <- rbind(tmp.domain.data, tmp.domain)
        }
      }
      
      domain.position.list <- as.data.frame(matrix(unlist(str_split(tmp.domain.data[, 5], "-")), ncol=2, byrow=TRUE))
      domain.data <- cbind(tmp.domain.data[, 1:4], domain.position.list)
      colnames(domain.data) <- c("transcript_id", "pfam_id", "desc", "protein_position", "start", "end")
      
      domain.data <- domain.data[domain.data$transcript_id %in% transcript.id.list, ]
      
      if(nrow(domain.data) > 0) {
        pfam.id.list <- unique(domain.data$pfam_id)
        domain.plot.data <- data.frame()
        domain.text.data <- data.frame()
        
        for(pfam.id in pfam.id.list) {
          tmp.domain.plot.data <- domain.data[domain.data$pfam_id == pfam.id, ]
          pfam.desc <- tmp.domain.plot.data[1, "desc"]
          
          domain.region <- IRanges(start=as.numeric(tmp.domain.plot.data$start), end=as.numeric(tmp.domain.plot.data$end))
          domain.plot.region <- reduce(domain.region)
          
          for(i in 1:length(domain.plot.region)) {
            tmp.domain.region <- domain.plot.region[i]
            
            domain.plot.data <- rbind(domain.plot.data, c(start(tmp.domain.region), end(tmp.domain.region), (i - 1), (i - 1 + 0.5), "domain"))
            as.inter <- intersect(tmp.domain.region, as.region)
            exon.inter <- intersect(tmp.domain.region, as.exon.region)
            
            if(length(as.inter) > 0) {
              for(j in 1:length(as.inter)) {
                domain.plot.data <- rbind(domain.plot.data, c(start(as.inter[j]), end(as.inter[j]), (i - 1), (i - 1 + 0.5), "as"))
              }
            }
            
            if(length(exon.inter) > 0) {
              for(j in 1:length(exon.inter)) {
                domain.plot.data <- rbind(domain.plot.data, c(start(exon.inter[j]), end(exon.inter[j]), (i - 1), (i - 1 + 0.5), "exon"))
              }
            }
            
            if(start(tmp.domain.region) < exon.region[1]) {
              domain.text.data <- rbind(domain.text.data, c(exon.region[1], (i - 1 + 0.6), pfam.desc))
            } else {
              domain.text.data <- rbind(domain.text.data, c(start(tmp.domain.region), (i - 1 + 0.6), pfam.desc))
            }
          }
        }
        
        colnames(domain.plot.data) <- c("xmin", "xmax", "ymin", "ymax", "domain_type")
        colnames(domain.text.data) <- c("x", "y", "domain_desc")
        
        domain.plot.data$xmin <- as.numeric(domain.plot.data$xmin)
        domain.plot.data$xmax <- as.numeric(domain.plot.data$xmax)
        domain.plot.data$ymin <- as.numeric(domain.plot.data$ymin)
        domain.plot.data$ymax <- as.numeric(domain.plot.data$ymax)
        
        domain.text.data$x <- as.numeric(domain.text.data$x)
        domain.text.data$y <- as.numeric(domain.text.data$y)
        
        domain.plot <- ggplot() + geom_rect(data=domain.plot.data, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=domain_type)) + geom_text(data=domain.text.data, mapping=aes(x=x, y=y, label=domain_desc), hjust=0, vjust=0) + scale_fill_manual(values=c("domain"="#CAC7ED", "exon"="#9590FF", "as"=I("pink"))) + xlim(exon.region) + ylim(0, max(domain.plot.data$ymax) + 0.5) + ylab("Protein domain") + guides(y="none")
        domain.track <- tracks(domain.plot)
        
        prev.track.height <- as.track@heights
        as.track <- as.track + domain.track
        as.track@heights <- c(prev.track.height, 1)
      }
    }
    
    #RBP
    #if(tmp.anno$)
    as.track <- as.track + theme_tracks_sunset(bg="white") + theme(legend.position="none", axis.title.y=element_text(angle=0, vjust=0.5))
    
    if(heights.list != "") {
      as.track@heights <- heights.list
    }
    
    print(as.track)
  
    #PPI plot
    if(tmp.anno$isoform_PPI_a != "" && tmp.anno$isoform_PPI_b != "") {
      in.transcript.id <- str_split(tmp.anno$exon_inclusion_transcript_id, ",")[[1]]
      ex.transcript.id <- str_split(tmp.anno$exon_exclusion_transcript_id, ",")[[1]]
  
      ppi.transcript.list <- str_split(tmp.anno$isoform_PPI_a, ";")[[1]]
      ppi.list <- str_split(tmp.anno$isoform_PPI_b, ";")[[1]]
  
      ppi.node <- data.frame(node=in.transcript.id, type="exon inclusion transcript")
      ppi.node <- rbind(ppi.node, data.frame(node=ex.transcript.id, type="exon exclusion transcript"))
      ppi.relation <- data.frame()
  
      for(i in 1:length(ppi.transcript.list)) {
        ppi.transcript <- ppi.transcript.list[i]
        ppi.partner.list <- str_split(ppi.list[i], "/")[[1]]
  
        for(j in 1:length(ppi.partner.list)) {
          ppi.partner <- toupper(ppi.partner.list[j])
  
          if(ppi.partner %in% ppi.node$node) {
            next
          }
  
          ppi.node <- rbind(ppi.node, data.frame(node=ppi.partner, type="gene"))
          ppi.relation <- rbind(ppi.relation, data.frame(transcript_id=ppi.transcript, gene=ppi.partner))
        }
      }
  
      ppi.net <- graph_from_data_frame(d=ppi.relation, vertices=ppi.node, directed=FALSE)
      #ggraph(ppi.net, layout="circle") + geom_edge_link(color="black") + geom_node_point(aes(color=type), size=5) + theme_void() + theme(legend.position="bottom", legend.text = element_text(size=16), legend.title=element_blank(), legend.margin=margin(b=5))
      ggraph(ppi.net, layout="circle") + geom_edge_link(color="black") + geom_node_point(aes(color=type), size=5) + geom_node_text(aes(label=name), size=6, vjust=2) + theme_void() + theme(legend.position="bottom", legend.text = element_text(size=14), legend.title=element_blank(), legend.margin=margin(b=5))
    }
  }
}
