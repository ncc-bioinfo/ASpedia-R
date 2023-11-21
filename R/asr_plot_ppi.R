#' visualization of PPI network
#' 
#' If asr_annotation function result has protein protein interaction(PPI) information, we provide PPI network plot.
#'
#' @param annotation.result
#' asr_annotation function result.
#' @param gene.name
#' gene name to be visualization. 
#' @param as.id
#' list of AS ID to be visualization.
#' @param result.dir 
#' directory where PPI plots(.png files) are saved
#' 
#' @return 
#' @export
#'
#' @examples
#' ppi.result.dir <- system.file(“extdata/ppi_result”, package=“ASpediaR”)
#' 
#' ##using gene name
#' asr_plot_ppi(annotation.result, gene.name="FGFR2", result.dir=ppi.result.dir)
#' 
#' ##using AS ID
#' asr_plot_ppi(annotation.result, 
#'             as.id="chr10:121520169:121519979:121518829:121518682:121517463:121517319",
#'             result.dir=ppi.result.dir)

asr_plot_ppi <- function(annotation.result="", gene.name="", as.id="", result.dir="") {
  if(class(annotation.result) == "character" && annotation.result == "") {
    message("*** ERROR MESSAGE: Input annotation result is empty. Please check input annotation result.")
    return()
  } else if(is.null(nrow(annotation.result)) == TRUE || nrow(annotation.result) == 0) {
    message("*** ERROR MESSAGE: Input annotation result is empty. Please check input annotation result.")
    return()
  }
  
  if(gene.name == "" && as.id == "") {
    message("*** ERROR MESSAGE: Input gene name or AS ID is required. Please check input gene name or AS ID")
    return()
  }
  
  if(result.dir == "") {
    message("*** ERROR MESSAGE: No such output directory.")
    return()
  }
  
  if(!file.exists(result.dir)) {
    dir.create(result.dir)
  }
  
  loaded.packages <- tolower((.packages()))
  
  if(("stringr" %in% loaded.packages) == FALSE) {
    library(stringr)
  }
  
  if(("ggraph" %in% loaded.packages) == FALSE) {
    library(ggraph)
  }
  
  if(("igraph" %in% loaded.packages) == FALSE) {
    library(igraph)
  }
  
  #annotation result filtered by gene name or as ID
  if(gene.name != "") {
    if(length(gene.name == 1)) {
      annotation.result <- annotation.result[annotation.result$gene_symbol == gene.name, ]
    } else {
      annotation.result <- annotation.result[annotation.result$gene_symbol %in% gene.name, ]
    }
  } else if(as.id != "") {
    if(length(as.id) == 1) {
      annotation.result <- annotation.result[annotation.result$as_id == as.id, ] 
    } else {
      annotation.result <- annotation.result[annotation.result$as_id %in% as.id, ]
    }
  }
  
  for(anno.index in 1:nrow(annotation.result)) {
    tmp.anno <- annotation.result[anno.index, ]

    as.id <- tmp.anno$as_id
    gene.name <- tmp.anno$gene_symbol
    
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
      ppi.plot <- ggraph(ppi.net, layout="circle") + geom_edge_link(color="black") + geom_node_point(aes(color=type), size=5) + geom_node_text(aes(label=name), size=6, vjust=2) + theme_void() + theme(legend.position="bottom", legend.text = element_text(size=14), legend.title=element_blank(), legend.margin=margin(b=5))
      print(ppi.plot)
      
      result.file.name <- paste0(result.dir, "/", str_replace_all(as.id, ":", "_"), "_PPI.png")
      
      png(file=result.file.name, width=1024, height=768)
      print(ppi.plot)
      dev.off()
    } else {
      print(paste0(gene.name, "\t", as.id, " has not PPI."))
    }
  }
}
