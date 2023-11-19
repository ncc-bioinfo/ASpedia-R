#' gene enrichment test between annotation gene list and knowledge-based database gene list.
#' 
#' Gene enrichment test result are provided to tsv and plot format between annotation gene list and knowledge-based database gene list.
#'
#' @param annotation.gene.list
#' gene list from asr_annotation function result
#' @param gsea.gene.list 
#' gene list from reference
#' @param result.dir
#' directory where GSEA result(.tsv and .png file) are saved
#' @param gene.model
#' gene model of reference. One of Refseq, Ensembl, or GENCODE. 
#' @param genome.version 
#' genome version of reference. One of hg18, GRCh19, or GRCh38.
#' @return
#' @export
#'
#' @examples
#' ## gene model and genome version
#' annotation.gene.list <- unique(annotation.result$gene_symbol)
#' gsea.result.dir <- paste0(system.file("extdata", package="ASpediaR"), "/gsea_result")
#' mining_gsea(annotation.gene.list, gene.model="Ensembl", genome.version="GRCh38",
#'               result.dir=gsea.result.dir)
#' 
#' ## reference gene list from user input
#' test.gene.list.file.name <- system.file(“extdata”, “test_whole_gene.txt”, package=“ASpediaR”)
#' test.gene.list <- read.table(test.gene.list.file.name, header=FALSE, stringsAsFactors=FALSE)
#' gsea.result.dir <- paste0(system.file(“extdata”, package=“ASpediaR”), “/gsea_result”)
#' mining_gsea(annotation.gene.list, gsea.gene.list=test.gene.list, result.dir=gsea.result.dir)


mining_gsea <- function(annotation.gene.list="", gsea.gene.list="", gene.model="", genome.version="", result.dir="") {
  data.dir <- paste0(.libPaths()[1], "/ASpediaR/data")
  
  if(file.exists(data.dir) == FALSE) {
    dir.create(data.dir)
  }
  
  if(class(annotation.gene.list) == "NULL") {
    print("*** ERROR MESSAGE: Input annotation gene list is empty vector. Please check input annotation gene list.")
    return()
  } else if(class(annotation.gene.list) != "character") {
    print("*** ERROR MESSAGE: Input annotation gene list is not character vector. Please check input annotation gene list.")
    return()
  } else {
    if (length(annotation.gene.list) == 1 && annotation.gene.list[1] == "") {
      print("*** ERROR MESSAGE: Input annotation gene list is empty. Please check input annotation gene list.")
      return()
    }
  }
  
  reference.gene.list <- ""
  
  if(gene.model == "" || genome.version == ""){
    if(class(gsea.gene.list) == "NULL") {
      print("*** ERROR MESSAGE: Input gene model, genome version and reference gene list is empty. Please check your input. reference gene list or gene model and genome version are required.")
      return()
    } else if(class(gsea.gene.list) != "character") {
      print("*** ERROR MESSAGE: Input reference gene list is not character vector. Please check input reference gene list.")
      return()
    } else {
      if (length(gsea.gene.list) == 1 && gsea.gene.list[1] == "") {
        print("*** ERROR MESSAGE: Input gene model, genome version and reference gene list is empty. Please check your input. reference gene list or gene model and genome version are required.")
        return()
      } else {
        reference.gene.list <- gsea.gene.list
      }
    }
  } else {
    reference.gene.list.file.name <- paste0(data.dir, "/", gene.model, ".", genome.version, ".gene.txt")
    
    if(!file.exists(reference.gene.list.file.name)) {
      url =  paste0("http://combio.hanyang.ac.kr/aspedia_v2/data/gene_list/", gene.model, ".", genome.version, ".gene.txt")
      download.file(url, reference.gene.list.file.name, method="auto")
    }
    
    reference.gene.list <- read.table(reference.gene.list.file.name, stringsAsFactors=FALSE, header=FALSE)
    reference.gene.list <- unique(reference.gene.list$V1)
  }
  
  if(result.dir == "") {
    print("*** ERROR MESSAGE: No such output directory.")
    return()
  }
  
  if(!file.exists(result.dir)) {
    dir.create(result.dir)
  }
  
  loaded.packages <- tolower((.packages()))
  
  if(("epitools" %in% loaded.packages) == FALSE) {
    library(epitools)
  }
  
  if(("stringr" %in% loaded.packages) == FALSE) {
    library(stringr)
  }
  
  if(("ggplot2" %in% loaded.packages) == FALSE) {
    library(ggplot2)
  }
  
  db.file.name <- paste0(data.dir, "/mining.sqlite")

  if(!file.exists(db.file.name)) {
    url =  "http://combio.hanyang.ac.kr/aspedia_v2/data/sqlite/mining.sqlite"
    download.file(url, db.file.name, method="auto", mode="wb")
  }

  db.connection <- dbConnect(SQLite(), dbname=db.file.name)

  contigency_table <- function(seta,setb,setall) {
    seta = intersect(seta, setall);
    setb = intersect(setb, setall);
    iset = intersect(seta, setb);
    useta = setdiff(seta, setb);
    usetb = setdiff(setb, seta);
    mtx = matrix(c(length(iset),length(useta),length(usetb),length(setall)-length(iset)-length(useta)-length(usetb)),nrow=2);
  }

  mining.result <- dbGetQuery(db.connection, "select * from mining")
  pathway.list <- unique(mining.result$origin_pathway)
  result.data <- data.frame()

  for(pathway in pathway.list) {
    mining.data <- mining.result[mining.result$origin_pathway == pathway, ]
    mining.gene <- unique(mining.data$gene)
    split_pathway <- (str_split(pathway, "_")[[1]])
    mining.pathway <- paste(split_pathway[2:length(split_pathway)], collapse="_")

    stat_CP <- chisq.test(t(contigency_table(mining.gene, annotation.gene.list, reference.gene.list)))$p.value

    result.data <- rbind(result.data, c(pathway, mining.pathway, stat_CP))
  }

  colnames(result.data) <- c("origin_pathway", "pathway", "CP")

  result.data <- cbind(result.data, adjP=p.adjust(result.data$CP, method="BY"))
  result.data$CP[is.na(result.data$CP)] <- 1
  #result.data <- cbind(result.data, "log_CP"=(-log10(as.numeric(result.data$CP))))
  result.data <- cbind(result.data, "log_CP"=(-log10(as.numeric(result.data$adjP))))
  result.data <- result.data[order(result.data$log_CP, decreasing=TRUE), ]
  
  write.table(result.data[, c("pathway", "CP", "adjP", "log_CP")], paste0(result.dir, "/GSEA_result.tsv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  
  result.data <- result.data[1:7, ]

  gsea.result.plot <- ggplot(result.data, aes(x=log_CP, y=reorder(pathway, log_CP))) + geom_bar(stat="identity") + theme_light() + theme(axis.title.x = element_text(colour="black", size=20, face="bold"), panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), axis.line = element_line(colour="black"), axis.text.x = element_text(colour="black", size=16), axis.text.y = element_text(colour="black", size=16, face="bold", vjust=0), axis.ticks = element_blank(), plot.margin.x = NULL, legend.position='none') + xlab("-log10(Pvalue)") + ylab('') + labs(fill='') + scale_x_continuous(expand=c(0, 0), limits=c(0, max(result.data$log_CP) + 1))
  # + scale_y_continuous(expand=c(0, 0), limits=c(0, y.max.value + 1))
  #, position="dodge"
  #axis.text.x = element_text(colour="black", size=16, angle=30, face="bold", hjust=1)
  #y=reorder(pathway, -log_CP)
  
  print(gsea.result.plot)
  
  png(file=paste0(result.dir, "/GSEA_result.png"), width=800, height=600)
  print(gsea.result.plot)
  dev.off()
}
