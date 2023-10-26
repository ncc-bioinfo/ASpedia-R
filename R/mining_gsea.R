library(epitools)
library(stringr)
library(ggplot2)

#' Title
#'
#' @param gene.list 
#' @param whole.gene.list 
#'
#' @return
#' @export
#'
#' @examples
mining_gsea <- function(gene.list, whole.gene.list) {
  db.file.name <- paste0(system.file(packages="ASpediaR"), "/data/mining.sqlite")

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

    stat_CP <- chisq.test(t(contigency_table(mining.gene, gene.list, whole.gene.list)))$p.value

    result.data <- rbind(result.data, c(pathway, mining.pathway, stat_CP))
  }

  colnames(result.data) <- c("origin_pathway", "pathway", "CP")

  result.data <- cbind(result.data, adjP=p.adjust(result.data$CP, method="BY"))
  result.data$CP[is.na(result.data$CP)] <- 1
  #result.data <- cbind(result.data, "log_CP"=(-log10(as.numeric(result.data$CP))))
  result.data <- cbind(result.data, "log_CP"=(-log10(as.numeric(result.data$adjP))))
  result.data <- result.data[order(result.data$log_CP, decreasing=TRUE), ]
  result.data <- result.data[1:7, ]

  gsea.result.plot <- ggplot(result.data, aes(x=log_CP, y=reorder(pathway, log_CP))) + geom_bar(stat="identity") + theme_light() + theme(axis.title.x = element_text(colour="black", size=20, face="bold"), panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), axis.line = element_line(colour="black"), axis.text.x = element_text(colour="black", size=16), axis.text.y = element_text(colour="black", size=16, face="bold", vjust=0), axis.ticks = element_blank(), plot.margin.x = NULL, legend.position='none') + xlab("-log10(Pvalue)") + ylab('') + labs(fill='') + scale_x_continuous(expand=c(0, 0), limits=c(0, max(result.data$log_CP) + 1))
  # + scale_y_continuous(expand=c(0, 0), limits=c(0, y.max.value + 1))
  #, position="dodge"
  #axis.text.x = element_text(colour="black", size=16, angle=30, face="bold", hjust=1)
  #y=reorder(pathway, -log_CP)
  
  print(gsea.result.plot)
}
