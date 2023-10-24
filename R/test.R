source("R/asr_converter.R")
source("R/asr_annotation.R")
source("R/asr_plot.R")

gene.model <- "Ensembl"
genome.version <- "hg38"
gtf.file.name <- "D:/project/ASpedia-R/00.reference/Homo_sapiens.GRCh38.99.rm.mt.add.chr.gtf"
ioe.file.name <- "D:/project/ASpedia-R/SUPPA_result/Homo_sapiens.GRCh38.99.rm.mt.add.chr.ioe"
pvalue.cutoff <- 0.05
das.cutoff <- 0.1
result.dir <- "D:/project/ASpedia-R/test_result/anno_test"

program.list <- c("rMATS", "SUPPA", "spliceR")
#program.list <- c("SUPPA", "spliceR")
#program.list <- c("spliceR")
as.type.list <- c("A3SS", "A5SS", "SE", "MXE", "RI")

for(program in program.list) {
  if(program == "rMATS") {
    total.start.time <- Sys.time()
    
    tmp.result.dir <- paste0(result.dir, "/", program)
    
    if(file.exists(tmp.result.dir) == FALSE) {
      dir.create(tmp.result.dir)
    }
    
    log.file <- file(paste0(tmp.result.dir, "/run_log.txt"))
    sink(log.file, append=TRUE, type="output")
    
    start.time <- Sys.time()
    
    converter.result <- data.frame()
    
    for(as.type in as.type.list) {
     input.file.name <- paste0("D:/project/ASpedia-R/rMATS_result/ESRP/", as.type, ".MATS.JC.txt")
     tmp.converter.result <- asr_converter(input.file.name=input.file.name, program=program, as.type=as.type, gene.model=gene.model, genome.version=genome.version, pvalue.cutoff=pvalue.cutoff, das.cutoff=das.cutoff)
     converter.result <- rbind(converter.result, tmp.converter.result)
    }
    
    result.file.name <- paste0(tmp.result.dir, "/rMATS.converter.result")
    write.table(converter.result, result.file.name, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    
    end.time <- Sys.time()
    print("rMATS convert : ")
    print(nrow(converter.result))
    print(end.time - start.time)
    
    start.time <- Sys.time()
    annotation.result <- asr_annotation(converter.result, gene.model, genome.version)
    result.file.name <- paste0(tmp.result.dir, "/rMATS.annotation.result")
    write.table(annotation.result, result.file.name, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    end.time <- Sys.time()
    print("rMATS annotation : ")
    print(end.time - start.time)
    
    #start.time <- Sys.time()
    #asr_plot(annotation.result, gtf.file.name, gene.model, genome.version)
    #end.time <- Sys.time()
    #print("rMATS plot : ")
    #print(end.time - start.time)
    
    total.end.time <- Sys.time()
    print("rMATS total :")
    print(total.end.time - total.start.time)
  }else if(program == "SUPPA") {
    total.start.time <- Sys.time()
    
    tmp.result.dir <- paste0(result.dir, "/", program)
    
    if(file.exists(tmp.result.dir) == FALSE) {
      dir.create(tmp.result.dir)
    }
    
    log.file <- file(paste0(tmp.result.dir, "/run_log.txt"))
    sink(log.file, append=TRUE, type="output")
    
    start.time <- Sys.time()
    result.file.name <- paste0(tmp.result.dir, "/SUPPA.converter.result")
    input.file.name <- "D:/project/ASpedia-R/SUPPA_result/ESRP_dPSI_empirical.dpsi"
    converter.result <- asr_converter(input.file.name=input.file.name, program=program, gene.model=gene.model, genome.version=genome.version, pvalue.cutoff=pvalue.cutoff, das.cutoff=das.cutoff, gtf.file.name=gtf.file.name, ioe.file.name=ioe.file.name)
    write.table(converter.result, result.file.name, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    end.time <- Sys.time()
    print("SUPPA convert : ")
    print(nrow(converter.result))
    print(end.time - start.time)
    
    start.time <- Sys.time()
    annotation.result <- asr_annotation(converter.result, gene.model, genome.version)
    result.file.name <- paste0(tmp.result.dir, "/SUPPA.annotation.result")
    write.table(annotation.result, result.file.name, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    end.time <- Sys.time()
    print("SUPPA annotation : ")
    print(end.time - start.time)
    
    #start.time <- Sys.time()
    #asr_plot(annotation.result, gtf.file.name, gene.model, genome.version)
    #end.time <- Sys.time()
    #print("SUPPA plot : ")
    #print(end.time - start.time)
    
    total.end.time <- Sys.time()
    print("SUPPA total :")
    print(total.end.time - total.start.time)
  }else if(program =="spliceR") {
    total.start.time <- Sys.time()
    
    tmp.result.dir <- paste0(result.dir, "/", program)

    if(file.exists(tmp.result.dir) == FALSE) {
      dir.create(tmp.result.dir)
    }
    
    log.file <- file(paste0(tmp.result.dir, "/run_log.txt"))
    sink(log.file, append=TRUE, type="output")
    
    start.time <- Sys.time()
    result.file.name <- paste0(tmp.result.dir, "/spliceR.converter.result")
    input.file.name <- "D:/project/ASpedia-R/spliceR_result/valr_test/ESRP_Ensembl_add_tss_p_id_spliceR_result.txt"
    converter.result <- asr_converter(input.file.name=input.file.name, program=program, gene.model=gene.model, genome.version=genome.version, pvalue.cutoff=pvalue.cutoff)
    write.table(converter.result, result.file.name, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    end.time <- Sys.time()
    print("SUPPA convert : ")
    print(nrow(converter.result))
    print(end.time - start.time)
    
    start.time <- Sys.time()
    annotation.result <- asr_annotation(converter.result, gene.model, genome.version)
    result.file.name <- paste0(tmp.result.dir, "/spliceR.annotation.result")
    write.table(annotation.result, result.file.name, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    end.time <- Sys.time()
    print("SUPPA annotation : ")
    print(end.time - start.time)
    
    #start.time <- Sys.time()
    #asr_plot(annotation.result, gtf.file.name, gene.model, genome.version)
    #end.time <- Sys.time()
    #print("SUPPA plot : ")
    #print(end.time - start.time)
    
    total.end.time <- Sys.time()
    print("spliceR total :")
    print(total.end.time - total.start.time)
  }
}