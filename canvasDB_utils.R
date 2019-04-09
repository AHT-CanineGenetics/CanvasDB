

importRetNetGenes <- function(tmpAnnotationDir=tmpfileDir, TALK=FALSE){
  library("rvest")
  library('biomaRt')

  cat("Downloading RetNet genes webpage...\n")
  url <- "https://sph.uth.edu/retnet/sym-dis.htm"
  page <- read_html(url)

  links <- html_nodes(page, "table")[3] %>% html_nodes("a")
  gene.names <- html_text(links)
  genes.df <- data.frame(gene.names)
  
  cat("Querying Ensembl BioMart to get CF ensembl identifiers for gene names...\n")
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- genes.df$gene.names
  
  hs.results <- getBM(filters = "hgnc_symbol", attributes = c("hgnc_symbol", "ensembl_gene_id"), values = genes, mart = mart)
  genes2 <-merge(genes.df, hs.results, by.x = "gene.names", by.y = "hgnc_symbol", all.x = TRUE)
  
  ens.ids <- hs.results$ensembl_gene_id
  cf.result <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","cfamiliaris_homolog_ensembl_gene"),values=ens.ids,mart=mart)
  
  gene.list <- merge(genes2, cf.result, by.x = "ensembl_gene_id", by.y="ensembl_gene_id", all.x = TRUE)
  retnet.genes <- gene.list$cfamiliaris_homolog_ensembl_gene
  
  a<-retnet.genes[retnet.genes != ""]
  retnet.genes <- a[!is.na(a)]
  
  tmpDatafile <- paste(tmpAnnotationDir,"/tmp_retNetGenes.txt",sep="")
  write.table(retnet.genes, file=tmpDatafile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
 # gene.list[gene.list['gene.names'] == "TEAD1"]
  
  cat("Loading ensembl ids into database...\n")
  con <- connectToInhouseDB(annotationDatabaseName)
  
  retNet.table <- annotTables[["retNetGenes"]]
  
  query <- paste("DROP TABLE IF EXISTS ",retNet.table, sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  query <- paste("CREATE TABLE ",retNet.table,
                 "(ens_id varchar(20), ",
                 "PRIMARY KEY(ens_id)",
                 ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
  
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' REPLACE INTO TABLE ",retNet.table,";",sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  dbDisconnect(con)
  
  #cleanUpTmpfiles(tmpAnnotationDir)
  
}