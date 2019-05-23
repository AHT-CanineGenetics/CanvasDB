#################################################################
##
## File: canvasDB_annotate.R
##
## Author: Adam Ameur, Uppsala Genome Center
##
## Description: Functions for annotating SNPs and small indels
##
#################################################################


## Import dbSNP data, SIFT scores and PolyPhen2 scores into local databases to enable quick annotations
setupAnnotationDBtables <- function(tmpAnnotationDir=tmpfileDir, dbSNPversion=151, TALK=FALSE){
 
  ## Import dbSNP rsIds if not already in database...
  dbSNP.table <- annotTables[["dbsnp"]]
  
  con <- connectToInhouseDB(annotationDatabaseName)
  query <- paste("SHOW TABLES LIKE '%",dbSNPversion,"_single';",sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  if(nrow(tmp) == 0){
    query <- paste("CREATE TABLE ",dbSNP.table,
                   "(SNP_id varchar(20), ",
                   "name varchar(20), ",
                   "PRIMARY KEY(SNP_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
    
    tmp <- dbGetQuery_E(con,query,TALK=TALK)
  }
  dbDisconnect(con)
  
  con <- connectToInhouseDB()
  query <- paste("SELECT * FROM ",dbSNP.table,";",sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  if(nrow(tmp) == 0){
    
    cat("Importing dbSNP ids into database...\n")
    
    dbsnpFile <- paste(tmpAnnotationDir,"cf3_snp",dbSNPversion,"_parsed.txt",sep="")
    
    if(!file.exists(dbsnpFile)){
      dbsnpRawFile <- paste(ANNOVARpathDB,"cf3_snp",dbSNPversion,".txt",sep="")
      
      if(!file.exists(dbsnpRawFile)){
        stop("dbsnp file missing! Download using ANNOVAR.\n")
      }
      
      parseCmd <- paste("utils/prepare_annot_dbsnp_tables.pl ",dbsnpRawFile," > ",dbsnpFile,sep="")
      system(parseCmd)
    }
    
    query <- paste("CREATE TABLE IF NOT EXISTS ",dbSNP.table,
                   "(SNP_id varchar(20), ",
                   "name varchar(20), ",
                   "PRIMARY KEY(SNP_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
    
    tmp <- dbGetQuery_E(con,query,TALK=TALK)
    
    query <- paste("LOAD DATA LOCAL INFILE '",dbsnpFile,"' REPLACE INTO TABLE ",dbSNP.table,";",sep="")
    
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    
  }
  dbDisconnect(con)
  
  ## Import dbSNP indel rsIds if not already in database...
  dbSNP.tableIndels <- annotTables[["dbsnpIndels"]]
  
  con <- connectToInhouseDB(annotationDatabaseName)
  query <- paste("SHOW TABLES LIKE '%",dbSNPversion,"_indels';",sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  if(nrow(tmp) == 0){
    query <- paste("CREATE TABLE ",dbSNP.tableIndels,
                   "(indel_id varchar(500), ",
                   "name varchar(20), ",
                   "PRIMARY KEY(indel_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
    
    tmp <- dbGetQuery_E(con,query,TALK=TALK)
  }
  dbDisconnect(con)
  
  con <- connectToInhouseDB()
  query <- paste("SELECT * FROM ",dbSNP.tableIndels,";",sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  if(nrow(tmp) == 0){
    
    cat("Importing dbSNP indel ids into database...\n")
    
    dbsnpFileIndels <- paste(tmpAnnotationDir,"cf3_snp",dbSNPversion,"_indels_parsed.txt",sep="")
    
    if(!file.exists(dbsnpFileIndels)){
      dbsnpRawFile <- paste(ANNOVARpathDB,"cf3_snp",dbSNPversion,".txt",sep="")
      
      if(!file.exists(dbsnpRawFile)){
        stop("dbsnp file missing! Download using ANNOVAR.\n")
      }
      
      parseCmd <- paste("utils/prepare_annot_dbsnp_tables_indels.pl ",dbsnpRawFile," > ",dbsnpFileIndels,sep="")
      system(parseCmd)
    }
    
    query <- paste("CREATE TABLE IF NOT EXISTS ",dbSNP.tableIndels,
                   "(indel_id varchar(500), ",
                   "name varchar(20), ",
                   "PRIMARY KEY(indel_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
    
    tmp <- dbGetQuery_E(con,query,TALK=TALK)
    
    query <- paste("LOAD DATA LOCAL INFILE '",dbsnpFileIndels,"' REPLACE INTO TABLE ",dbSNP.tableIndels,";",sep="")
    
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    
  }
  dbDisconnect(con)
  
  ## Import SNP scores if not already in database...
  score.table <- annotTables[["score"]]
  query <- paste("SHOW TABLES LIKE '",score.table,"';",sep="")
  con <- connectToInhouseDB()
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  if(nrow(tmp) == 0){
    
    cat("Imporing scores into database...\n")
    
#    scoreFile <- paste(tmpAnnotationDir,"hg19_ljb_all_parsed.txt",sep="")
#    
#    if(!file.exists(scoreFile)){
#      scoreRawFile <- paste(ANNOVARpathDB,"hg19_ljb_all.txt",sep="")
#      
#      if(!file.exists(scoreRawFile)){
#        stop("ljb_all file missing! Download using ANNOVAR.\n")
#      }
#      
#      parseCmd <- paste("utils/prepare_annot_score_tables.pl ",scoreRawFile," > ",scoreFile,sep="")
#      system(parseCmd)
#    }
#    
#    
#    query <- paste("CREATE TABLE ",score.table,
#                   "(SNP_id varchar(20), ",
#                   "sift double(5,2), ",
#                   "polyphen double(5,2), ",
#                   "phylop double(5,2), ",
#                   "lrt double(5,2), ",
#                   "mut_taster double(5,2), ",
#                   "gerp double(5,2), ",
#                   "PRIMARY KEY(SNP_id)",
#                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
    
    query <- paste("CREATE TABLE ",score.table,
                   "(SNP_id varchar(20), ",
                   "sift double(5,2), ",
                   "PRIMARY KEY(SNP_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
    
    tmp <- dbGetQuery_E(con,query,TALK=TALK)
    
    siftDataDir <- paste(snpEffPath, "/data/", siftDB, sep="")
    files <- list.files(path=siftDataDir, pattern="*.gz", full.names=TRUE, recursive=FALSE)
    lapply(files, function(x) {
      t <- read.table(x, header=FALSE, sep="\t")
      colnames(t)<- c("Position","Ref_allele","New_allele","Transcript_id","Gene_id","Gene_name","Region","Ref_amino_acid","New_amino_acid","Position_of_amino_acid_substitution","SIFT_score","SIFT_median_sequence_info","Num_seqs_at_position","dbSNP_id")
      t <- t[-1*grep("ref",t[,14]),]
      
      chr <- gsub(".gz", "", basename(x))
      SNPids <- paste(chr,t[,"Position"],t[,"Ref_allele"],t[,"New_allele"],sep="|")
      
      sift.scores <- as.character(t[,"SIFT_score"])
      tmp <- cbind(SNPids, sift.scores)
      
      siftScorefile <- paste(tmpfileDir,"sift_scores.txt",sep="")
      write.table(tmp,file=siftScorefile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
      
      query <- paste("LOAD DATA LOCAL INFILE '",siftScorefile,"' REPLACE INTO TABLE ",score.table,";",sep="")
      tmp <- dbGetQuery_E(con,query,TALK=TALK)
      
    })
    
#    query <- paste("LOAD DATA LOCAL INFILE '",scoreFile,"' REPLACE INTO TABLE ",score.table,";",sep="")
    
#    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    
    dbDisconnect(con)
  }
  
  
  
  addVEPAnnotation(dbSNPversion = dbSNPversion, eVersion = eVersion);

}

addSnpEffAnnotation <- function(tmpAnnotationDir=tmpfileDir, dbSNPversion=151, snpEffDb="CanFam3.1.86", TALK=FALSE){
  
  ## Annotate SNPS with VEP if not already in database...
  cat("Importing SnpEff Annotation into database...\n")
  tmpDatafile <- paste(tmpAnnotationDir,"/tmp_variantsToBeAnnotated.vcf",sep="")
  
  dbSNP.table <- annotTables[["dbsnp"]]
  con <- connectToInhouseDB()
  query <- paste("SELECT SNP_id FROM ",dbSNP.table,";",sep="")
  inputSNPs <- dbGetQuery_E(con, query, TALK=TALK)
  dbDisconnect(con)
  
  inputSNPids <- as.character(inputSNPs[,"SNP_id"])
  inputSNPidsSplit <- strsplit(inputSNPids,"\\|")
  inputChr <- sapply(inputSNPidsSplit,function(x){x[[1]]})
  inputPos <- sapply(inputSNPidsSplit,function(x){x[[2]]})
  inputRef <- sapply(inputSNPidsSplit,function(x){x[[3]]})
  inputAlt <- sapply(inputSNPidsSplit,function(x){x[[4]]})
  
  SNPdata <- unique(cbind(inputChr,inputPos,inputSNPids,inputRef,inputAlt,".",".","."))
  colnames(SNPdata) <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  
  write.table(SNPdata, file=tmpDatafile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
}

addVEPAnnotation <- function(tmpAnnotationDir=tmpfileDir, eVer=eVersion, TALK=FALSE){
  
  ## Annotate SNPS with VEP if not already in database...
  cat("Importing VEP Annotation into database...\n")
  
  dbSNP.table <- annotTables[["dbsnp"]]
  dbSNPIndel.table <- annotTables[["dbsnpIndels"]]
  vepSNP.table <- annotTables[["vepSNP"]]
  sumSNP.table <- dbTables[["SNP.summary"]]
  sumIndel.table <- dbTables[["indel.summary"]]
  
  updateDbSnps <- FALSE
  
  con <- connectToInhouseDB(annotationDatabaseName)
  query <- paste("SHOW TABLES LIKE 'annot_vep",eVer,"_single';",sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  if(nrow(tmp) == 0){
    updateDbSnps <- TRUE
    query <- paste("CREATE TABLE IF NOT EXISTS ",vepSNP.table,
                   "(SNP_id varchar(20), ",
                   "position varchar(20), ",
                   "ref_allele varchar(5), ",
                   "cf_gene varchar(20), ",
                   "cf_feat varchar(20), ",
                   "feature_type varchar(20), ",
                   "consequence varchar(50), ",
                   "cDNA_pos varchar(5), ",
                   "cds_pos varchar(5), ",
                   "prot_pos varchar(5), ",
                   "amino_acids varchar(20), ",
                   "codons varchar(20), ",
                   "dbsnp varchar(20), ",
                   "cf_symbol varchar(20), ",
                   "severity int(1), ",
                   "biotype varchar(50), ",
                   "impact varchar(50), ",
                   "sift_text varchar(50), ",
                   "sift_val decimal(5,3) default NULL, ", 
                   "annotation text, ",
                   "hs_gene varchar(20), ",
                   "PRIMARY KEY(SNP_id, cf_feat)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
    tmp <- dbGetQuery_E(con,query,TALK=TALK)
  }
  dbDisconnect(con)
  
  if (updateDbSnps){
    cat("Annotating dbSNP entries with VEP annotation...\n")
    con <- connectToInhouseDB()
    query <- paste("SELECT SNP_id FROM ",dbSNP.table,";",sep="")
    dbSNPs <- dbGetQuery_E(con, query, TALK=TALK)
    query <- paste("SELECT SNP_id FROM ",dbSNPIndel.table,";",sep="")
    dbSNPIndels <- dbGetQuery_E(con, query, TALK=TALK)
    dbDisconnect(con)
    
    runVEP(inputSNPs=dbSNPs, eVer=eVer)
    runVEP(inputSNPs=dbSNPIndels, eVer=eVer)
  }
  
  cat("Checking for un-annotated SNPs/Indels in the summary tables.....")
  
  con <- connectToInhouseDB()
  query <- paste("SELECT SNP_id FROM ",vepSNP.table,";",sep="")
  okSNPs <- dbGetQuery_E(con, query, TALK=TALK)
  
  query2 <- paste("SELECT SNP_id FROM ",sumSNP.table,";",sep="")
  summarySNPs <- dbGetQuery_E(con, query2, TALK=TALK)
  
  query3 <- paste("SELECT indel_id as SNP_id FROM ",sumIndel.table,";",sep="")
  summaryIndels <- dbGetQuery_E(con, query3, TALK=TALK)
  dbDisconnect(con)
  
  vepSNPs <- subset(summarySNPs, !(SNP_id %in% okSNPs$SNP_id))
  vepIndels <- subset(summaryIndels, !(SNP_id %in% okSNPs$SNP_id))
  
  cat(paste(": SNPS=", nrow(vepSNPs), "; Indels=", nrow(vepIndels), "\n" ,sep=""))
  if(nrow(vepSNPs) > 0){
    runVEP(inputSNPs=vepSNPs, eVer=eVer)
    runVEP(inputSNPs=vepIndels, eVer=eVer)
  }
  
}

runVEP <- function(inputSNPs, tmpAnnotationDir=tmpfileDir, eVer=eVersion, TALK=FALSE){
  
  tmpDatafile <- paste(tmpAnnotationDir,"/tmp_variantsToBeAnnotated.vcf",sep="")
  tmpVEPfile <- paste(tmpAnnotationDir,"/tmp_variantsAnnotated_vep",eVer,".txt",sep="")
  dockerDatafile <- "/data/tmp_variantsToBeAnnotated.vcf"
  dockerVEPfile <- paste("/data/tmp_variantsAnnotated_vep",eVer,".txt",sep="")
  vepSNP.table <- annotTables[["vepSNP"]]
  
  inputSNPids <- as.character(inputSNPs[,"SNP_id"])
  inputSNPidsSplit <- strsplit(inputSNPids,"\\|")
  inputChr <- sapply(inputSNPidsSplit,function(x){x[[1]]})
  inputPos <- sapply(inputSNPidsSplit,function(x){x[[2]]})
  inputRef <- sapply(inputSNPidsSplit,function(x){x[[3]]})
  inputAlt <- sapply(inputSNPidsSplit,function(x){x[[4]]})
  
  SNPdata <- unique(cbind(inputChr,inputPos,inputSNPids,inputRef,inputAlt,".",".","."))
  colnames(SNPdata) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  
  write.table(SNPdata, file=tmpDatafile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  vepCommand <- paste(VEPpath, " --cache --offline --species canis_familiaris --format vcf ",
                      "--force_overwrite --no_stats --everything --no_headers --fork 4 ",
                      "--input_file ",dockerDatafile,
                      "--output_file ", dockerVEPfile)
  print(vepCommand)
  system(vepCommand)
  
  tmpfile <- paste(tmpAnnotationDir,"/tmp_vep.txt", sep="")
  system(paste("perl ",utilsDir, "/parse_vep_output.pl ",tmpVEPfile," > ",tmpfile,sep=""))
  variantFunctions <- read.table(tmpfile, header=FALSE, as.is=TRUE, sep="\t", comment.char="#")
  
  con <- connectToInhouseDB()
  
  query <- paste("LOAD DATA LOCAL INFILE '",tmpfile,"' REPLACE INTO TABLE ",vepSNP.table," (SNP_id, position, ref_allele, cf_gene, cf_feat, feature_type, consequence, cDNA_pos, cds_pos, prot_pos, amino_acids, codons, dbsnp, cf_symbol, severity, biotype, impact, sift_text, sift_val, annotation);",sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  file.remove(tmpVEPfile)
  file.remove(tmpfile)
  
  dbDisconnect(con)

}




## Remove temporary files created by annovar
cleanUpTmpfiles <- function(tmpAnnotationDir){
  
  tmpfiles <- paste(tmpAnnotationDir,
                    c("tmp_exonicVariantsToBeAnnotated.txt",
                      "tmp_exonicVariantsToBeAnnotated.txt.log",
                      "tmp_variantsToBeAnnotated.txt",
                      "tmp_variantsToBeAnnotated.txt.exonic_variant_function",
                      "tmp_variantsToBeAnnotated.txt.invalid_input",
                      "tmp_variantsToBeAnnotated.txt.log",
                      "tmp_variantsToBeAnnotated.txt.variant_function"),sep="")
  
  for(tmpfile in tmpfiles){
    if(file.exists(tmpfile)){
      file.remove(tmpfile)
    }
  }
}





## Annotate SNPs with the following information
##
## 1. dbSNP and dbSNPCommon - Assign rsId if available
## 2. refgene - Determine effect on transcript level
## 3. Scores - Amino acid substitution effects
annotateSNPs <- function(inputSNPs, tmpAnnotationDir=".", dbSNPversion=151, TALK=FALSE){
  
  ps <- proc.time()[3]
  
  SNPseverity <- array()
  SNPseverity["stopgain"] <- 5
  SNPseverity["nonsynonymous"] <- 4
  SNPseverity["stoploss"] <- 4
  SNPseverity["splicing"] <- 4
  SNPseverity["exonic;splicing"] <- 4
  SNPseverity["ncRNA_exonic"] <- 2
  SNPseverity["ncRNA_intronic"] <- 1
  SNPseverity["intergenic"] <- 1
  SNPseverity["intronic"] <- 1
  SNPseverity["synonymous"] <- 1
  SNPseverity["UTR5"] <- 1
  SNPseverity["UTR3"] <- 1
  SNPseverity["downstream"] <- 1
  SNPseverity["upstream"] <- 1
  SNPseverity["ncRNA_UTR3"] <- 1
  SNPseverity["upstream;downstream"] <- 1
  SNPseverity["ncRNA_splicing"] <- 1
  SNPseverity["exonic"] <- 2
  SNPseverity["ncRNA_UTR5"] <- 1
  SNPseverity["UTR5;UTR3"] <- 1
  
  inputSNPids <- as.character(inputSNPs[,"SNP_id"])
  inputSNPidsSplit <- strsplit(inputSNPids,"\\|")
  inputChr <- sapply(inputSNPidsSplit,function(x){x[[1]]})
  inputPos <- sapply(inputSNPidsSplit,function(x){x[[2]]})
  inputRef <- sapply(inputSNPidsSplit,function(x){x[[3]]})
  inputAlt <- sapply(inputSNPidsSplit,function(x){x[[4]]})
  inputNrSamples <- as.numeric(inputSNPs[,"nr_samples"])
  inputSampleStr <- as.character(inputSNPs[,"sample_str"])
  
  colnames <- c("chr","pos","pos","ref","alt","SNP_id")
  
  tmpDatafile <- paste(tmpAnnotationDir,"/tmp_variantsToBeAnnotated.txt",sep="")
  SNPdata <- unique(cbind(inputChr,inputPos,inputPos,inputRef,inputAlt,inputSNPids,inputNrSamples,inputSampleStr))
  colnames(SNPdata) <- c("chr","pos","pos","ref","alt","SNP_id","nr_samples","samples")
  dataToBeAnnotated <- SNPdata[,colnames,drop=FALSE]
  dataToBeAnnotated[,"chr"] <- sub("chr","",dataToBeAnnotated[,"chr"])
  write.table(dataToBeAnnotated, file=tmpDatafile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  con <- connectToInhouseDB()
  
  tmpSNP.table <- "tmp_SNP_to_be_annotated"
  
  query <- paste("DROP TABLE IF EXISTS ",tmpSNP.table, sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  query <- paste("CREATE TABLE ",tmpSNP.table,
                 "(chr varchar(10), ",
                 "start integer(11), ",
                 "end integer(11), ",
                 "ref varchar(1), ",
                 "alt varchar(1), ",
                 "SNP_id varchar(20), ",
                 "PRIMARY KEY(SNP_id)",
                 ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
  
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",tmpSNP.table,";",sep="")
  
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  dbDisconnect(con)
  
#  dataToBeAdded <- matrix("", nrow=nrow(SNPdata), ncol=14)
#  colnames(dataToBeAdded) <- c("SNP_id","chr","pos","ref","alt","nr_samples","samples","dbSNP","class","severity","gene","details","sift","polyphen")#  dataToBeAdded[, c("SNP_id","chr","pos","ref","alt","nr_samples","samples")] <- as.matrix(SNPdata[,c("SNP_id","chr","pos","ref","alt","nr_samples","samples")])
  dataToBeAdded <- matrix("", nrow=nrow(SNPdata), ncol=8)
  colnames(dataToBeAdded) <- c("SNP_id","chr","pos","ref","alt","nr_samples","samples","dbSNP")
  dataToBeAdded[, c("SNP_id","chr","pos","ref","alt","nr_samples","samples")] <- as.matrix(SNPdata[,c("SNP_id","chr","pos","ref","alt","nr_samples","samples")])
  
  
  
  ## Step 1. Filter against dbSNP and dbSNPCommon
  cat("  - Annotating SNPs with dbSNP rs-ids...",sep="")
  con <- connectToInhouseDB()
  query <- paste("SELECT t1.SNP_id,t2.name FROM ",tmpSNP.table," as t1, ",annotTables[["dbsnp"]]," as t2 WHERE t1.SNP_id=t2.SNP_id;", sep="")
  cat(query)
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  dbDisconnect(con)
  if(nrow(tmp)>0){
    dataToBeAdded[match(tmp[,"SNP_id"],dataToBeAdded[,"SNP_id"]),"dbSNP"] <- tmp[,"name"]
  }
  
  cat(proc.time()[3] - ps,"s\n");
  
  
#   ## Step 2. Annotate exonic variants against score table
#   cat("  - Annotating SNPs with SNP scores...",sep="")
#   con <- connectToInhouseDB()
#   query <- paste("SELECT t1.SNP_id,t2.sift FROM ",tmpSNP.table," as t1, ",annotTables[["score"]]," as t2 WHERE t1.SNP_id=t2.SNP_id;", sep="")
#   cat(query)
#   tmp <- dbGetQuery_E(con, query, TALK=TALK)
#   dbDisconnect(con)
#   if(nrow(tmp)>0){
#     idsToBeAdded <- match(tmp[,"SNP_id"],dataToBeAdded[,"SNP_id"])
#     dataToBeAdded[idsToBeAdded,"sift"] <- tmp[,"sift"]
#   }
#   cat(proc.time()[3] - ps,"s\n");
#   
#   con <- connectToInhouseDB()
#   query <- paste("DROP TABLE IF EXISTS ",tmpSNP.table, sep="")
#   tmp <- dbGetQuery_E(con, query, TALK=TALK)
#   dbDisconnect(con)
#   
#   ## Step 3. Annotate against refSeq genes
# #  cat("  - Annotating SNPs against refGene...",sep="")
# #  ANNOVARcmd <- paste(ANNOVARpath," -geneanno -buildver canFam3 -dbtype refgene ",tmpDatafile," ",ANNOVARpathDB,sep="")
#   cat("  - Annotating SNPs against ensGene...",sep="")
#   ANNOVARcmd <- paste(ANNOVARpath," -geneanno -buildver canFam3 -dbtype ensgene ",tmpDatafile," ",ANNOVARpathDB,sep="")
#   cat(ANNOVARcmd)
#   system(ANNOVARcmd,ignore.stderr=FALSE)
#   
#   variantFunctionFile <- paste(tmpDatafile,".variant_function",sep="")
#   exonicVariantFunctionFile <- paste(tmpDatafile,".exonic_variant_function",sep="")
#   
#   variantFunctionFile
#   file.info(variantFunctionFile)
#   
#   if(file.info(variantFunctionFile)$size > 0){
#     variantFunctions <- read.table(variantFunctionFile, as.is=TRUE)
#     variantCategory <- as.character(variantFunctions[,1])
#     names(variantCategory) <- as.character(variantFunctions[,8])
#     variantGene <- as.character(variantFunctions[,2])
#     names(variantGene) <- as.character(variantFunctions[,8])
#     
#     dataToBeAdded[match(names(variantGene),dataToBeAdded[,"SNP_id"]),"gene"] <- variantGene
#     dataToBeAdded[match(names(variantCategory),dataToBeAdded[,"SNP_id"]),"class"] <- variantCategory
#   }
#   
#   if(file.info(exonicVariantFunctionFile)$size > 0){
#     exonicVariantFunctions <- read.table(exonicVariantFunctionFile, sep="\t")
#     exonicVariantFunctions <- exonicVariantFunctions[which(exonicVariantFunctions[,2] != "unknown"),]
#     exonicVariantCategory <- as.character(exonicVariantFunctions[,2])
#     names(exonicVariantCategory) <- as.character(exonicVariantFunctions[,ncol(exonicVariantFunctions)])
#     exonicVariantDetail <- as.character(exonicVariantFunctions[,3])
#     names(exonicVariantDetail) <- as.character(exonicVariantFunctions[,ncol(exonicVariantFunctions)])
#     
#     dataToBeAdded[match(names(exonicVariantCategory),dataToBeAdded[,"SNP_id"]),"class"] <- exonicVariantCategory
#     dataToBeAdded[match(names(exonicVariantDetail),dataToBeAdded[,"SNP_id"]),"details"] <- exonicVariantDetail
#   }
#   
#   dataToBeAdded[,"class"] <- gsub(" SNV","",dataToBeAdded[,"class"])
#   dataToBeAdded[, "severity"] <- SNPseverity[dataToBeAdded[, "class"]]
#   
#   cat(proc.time()[3] - ps,"s\n");
#   
#   con <- connectToInhouseDB()
#   query <- paste("DROP TABLE IF EXISTS ",tmpSNP.table, sep="")
#   tmp <- dbGetQuery_E(con, query, TALK=TALK)
#   dbDisconnect(con)
  
  cleanUpTmpfiles(tmpAnnotationDir)
  
  return(dataToBeAdded)
  
}



##
## Annotate indels with the following information
##
## 1. dbSNP and dbSNPCommon - Assign rsId if available
## 2. refgene - Determine effect on transcript level
annotateIndels <- function(inputIndels, tmpAnnotationDir=".", dbSNPversion=151, TALK=FALSE){
  
  ps <- proc.time()[3]
  
  indelSeverity <- array()
  indelSeverity["frameshift"] <- 5
  indelSeverity["frameshift deletion"] <- 5
  indelSeverity["frameshift insertion"] <- 5
  indelSeverity["stopgain"] <- 5
  indelSeverity["stopgain SNV"] <- 5
  indelSeverity["stoploss SNV"] <- 5
  indelSeverity["nonframeshift"] <- 3
  indelSeverity["nonframeshift deletion"] <- 3
  indelSeverity["nonframeshift insertion"] <- 3
  indelSeverity["exonic"] <- 3
  indelSeverity["exonic;splicing"] <- 3
  indelSeverity["splicing"] <- 4
  indelSeverity["UTR3"] <- 2
  indelSeverity["UTR5"] <- 2
  indelSeverity["ncRNA_exonic"] <- 2
  indelSeverity["ncRNA_splicing"] <- 2
  indelSeverity["ncRNA_UTR5"] <- 1
  indelSeverity["ncRNA_UTR3"] <- 1
  indelSeverity["ncRNA_intronic"] <- 1
  indelSeverity["intronic"] <- 1
  indelSeverity["intergenic"] <- 1
  indelSeverity["upstream"] <- 1
  indelSeverity["upstream;downstream"] <- 1
  indelSeverity["downstream"] <- 1
  indelSeverity["unknown"] <- 1
  
  inputIndelIds <- as.character(inputIndels[,"indel_id"])
  inputIndelIdsSplit <- strsplit(inputIndelIds,"\\|")
  inputChr <- sapply(inputIndelIdsSplit,function(x){x[[1]]})
  inputStart <- sapply(inputIndelIdsSplit,function(x){x[[2]]})
  inputEnd <- sapply(inputIndelIdsSplit,function(x){x[[3]]})
  inputRef <- sapply(inputIndelIdsSplit,function(x){x[[4]]})
  inputAlt <- sapply(inputIndelIdsSplit,function(x){x[[5]]})
  inputNrSamples <- as.numeric(inputIndels[,"nr_samples"])
  inputSampleStr <- as.character(inputIndels[,"sample_str"])
  inputType <- rep(NA,length(inputIndelIds))
  inputSize <- rep(NA,length(inputIndelIds))
  
  colnames <- c("chr","start","end","ref","alt","indel_id")
  
  tmpDatafile <- paste(tmpAnnotationDir,"/tmp_variantsToBeAnnotated.txt",sep="")
  indelData <- unique(cbind(inputChr,inputStart,inputEnd,inputRef,inputAlt,inputIndelIds,inputNrSamples,inputSampleStr,inputType,inputSize))
  colnames(indelData) <- c("chr","start","end","ref","alt","indel_id","nr_samples","samples","type","size")
  
  indelDataTmp <- indelData
  
  ## Create synthetic DNA strings for ref/alt of correct length. Should be ok since ANNOVAR anyway doesn't use the actual sequence..
  insertIds <- which(indelDataTmp[,"ref"] == "-")
  deletionIds <- which(indelDataTmp[,"alt"] == "-")
  
  if(length(insertIds)>0){
    indelData[insertIds,"type"] <- "ins"
    indelData[insertIds,"size"] <- sapply(insertIds, function(i){nchar(indelDataTmp[i,"alt"])})
  }
  
  if(length(deletionIds)>0){
    indelData[deletionIds,"type"] <- "del"
    indelData[deletionIds,"size"] <- sapply(deletionIds, function(i){nchar(indelDataTmp[i,"ref"])})
    indelDataTmp[deletionIds,"start"] <- as.numeric(indelDataTmp[deletionIds,"start"])+1 ## Required for ANNOVAR annotation of deletions.
  }
  
  otherIds <- which(!(indelData[,"type"] %in% c("ins","del")))
  indelData[otherIds,"type"] <- "indel"
  
  dataToBeAnnotated <- indelDataTmp[,colnames,drop=FALSE]
  dataToBeAnnotated[,"chr"] <- sub("chr","",dataToBeAnnotated[,"chr"])
  
  write.table(dataToBeAnnotated, file=tmpDatafile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  con <- connectToInhouseDB()
  
  tmpIndel.table <- "tmp_indels_to_be_annotated"
  
  query <- paste("DROP TABLE IF EXISTS ",tmpIndel.table, sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  query <- paste("CREATE TABLE ",tmpIndel.table,
                 "(chr varchar(10), ",
                 "start integer(11), ",
                 "end integer(11), ",
                 "ref varchar(1), ",
                 "alt varchar(1), ",
                 "indel_id varchar(500), ",
                 "PRIMARY KEY(indel_id)",
                 ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")
  
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",tmpIndel.table,";",sep="")
  
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  
  dbDisconnect(con)
  
  dataToBeAdded <- matrix("", nrow=nrow(dataToBeAnnotated), ncol=11)
  colnames(dataToBeAdded) <- c("indel_id","chr","start","end","ref","alt","type","size","nr_samples","samples","dbSNP")
  dataToBeAdded[, c("indel_id","chr","start","end","ref","alt")] <- as.matrix(dataToBeAnnotated[,c("indel_id","chr","start","end","ref","alt")])
  
  ## Step 1. Filter against dbSNP and dbSNPCommon
  cat("  - Annotating indels with dbSNP rs-ids...",sep="")
  con <- connectToInhouseDB()
  query <- paste("SELECT t1.indel_id,t2.name FROM ",tmpIndel.table," as t1, ",annotTables[["dbsnpIndels"]]," as t2 WHERE t1.indel_id=t2.indel_id;", sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  dbDisconnect(con)
  if(nrow(tmp)>0){
    dataToBeAdded[match(tmp[,"indel_id"],dataToBeAdded[,"indel_id"]),"dbSNP"] <- tmp[,"name"]
  }
  
#  con <- connectToInhouseDB()
#  query <- paste("SELECT t1.indel_id,t2.name FROM ",tmpIndel.table," as t1, ",annotTables[["dbsnpCommonIndels"]]," as t2 WHERE t1.indel_id=t2.indel_id;", sep="")
#  tmp <- dbGetQuery_E(con, query, TALK=TALK)
#  dbDisconnect(con)
#  if(nrow(tmp)>0){
#    dataToBeAdded[match(tmp[,"indel_id"],dataToBeAdded[,"indel_id"]),"dbSNP"] <- tmp[,"name"]
#    dataToBeAdded[match(tmp[,"indel_id"],dataToBeAdded[,"indel_id"]),"dbSNPcommon"] <- tmp[,"name"]
#  }
  
  cat(proc.time()[3] - ps,"s\n");
  
  # ## Step 2. Annotate against refSeq genes
  # 
  # cat("  - Annotating indels against refGene...",sep="")
  # ANNOVARcmd <- paste(ANNOVARpath," -geneanno -buildver canFam3 -dbtype refgene ",tmpDatafile," ",ANNOVARpathDB,sep="")
  # system(ANNOVARcmd,ignore.stderr=TRUE)
  # cat(proc.time()[3] - ps,"s\n");
  # 
  # variantFunctionFile <- paste(tmpDatafile,".variant_function",sep="")
  # exonicVariantFunctionFile <- paste(tmpDatafile,".exonic_variant_function",sep="")
  # 
  # ## #################################################
  # ##
  # ## Combine all annotations into a single output file
  # ##
  # ## #################################################
  # 
  # cat("  - Combining annotations and preparing output...",sep="")
  # 
  # dataToBeAdded[, c("indel_id","chr","start","end","ref","alt","type","size","nr_samples","samples")] <- as.matrix(indelData[,c("indel_id","chr","start","end","ref","alt","type","size","nr_samples","samples")])
  # 
  # ## Add refSeq annotations
  # if(file.info(variantFunctionFile)$size > 0){
  #   variantFunctions <- read.table(variantFunctionFile, as.is=TRUE)
  #   variantCategory <- as.character(variantFunctions[,1])
  #   names(variantCategory) <- as.character(variantFunctions[,8])
  #   variantGene <- as.character(variantFunctions[,2])
  #   names(variantGene) <- as.character(variantFunctions[,8])
  #   
  #   dataToBeAdded[match(names(variantGene),dataToBeAdded[,"indel_id"]),"gene"] <- variantGene
  #   dataToBeAdded[match(names(variantCategory),dataToBeAdded[,"indel_id"]),"class"] <- variantCategory
  # }
  # 
  # if(file.info(exonicVariantFunctionFile)$size > 0){
  #   exonicVariantFunctions <- read.table(exonicVariantFunctionFile, sep="\t")
  #   exonicVariantFunctions <- exonicVariantFunctions[which(exonicVariantFunctions[,2] != "unknown"),]
  #   exonicVariantCategory <- as.character(exonicVariantFunctions[,2])
  #   names(exonicVariantCategory) <- as.character(exonicVariantFunctions[,ncol(exonicVariantFunctions)])
  #   exonicVariantDetail <- as.character(exonicVariantFunctions[,3])
  #   names(exonicVariantDetail) <- as.character(exonicVariantFunctions[,ncol(exonicVariantFunctions)])
  #   
  #   dataToBeAdded[match(names(exonicVariantCategory),dataToBeAdded[,"indel_id"]),"class"] <- exonicVariantCategory
  #   dataToBeAdded[match(names(exonicVariantDetail),dataToBeAdded[,"indel_id"]),"details"] <- exonicVariantDetail
  # }
  # 
  # dataToBeAdded[, "severity"] <- indelSeverity[dataToBeAdded[, "class"]]

  cat(proc.time()[3] - ps,"s\n");
  
  con <- connectToInhouseDB()
  query <- paste("DROP TABLE IF EXISTS ",tmpIndel.table, sep="")
  tmp <- dbGetQuery_E(con, query, TALK=TALK)
  dbDisconnect(con)
  
  cleanUpTmpfiles(tmpAnnotationDir)
  
  return(dataToBeAdded)
  
}


rebuildSNPSummary <- function(){
  
  ps <- proc.time()[3]
  
  con <- connectToInhouseDB()
  sample.table <- dbTables[["sample"]]
  SNPsummary.table <- dbTables[["SNP.summary"]]
  chunkSize <- 500000
  
  sampleIds <- list()
  SNPidsInSummary <- NULL
  
  query <- paste("SELECT SNP_id FROM ",SNPsummary.table,";")
  currSNPids <- dbGetQuery_E(con,query,TALK=FALSE)
  if(nrow(currSNPids)>0){ SNPidsInSummary <- currSNPids[,"SNP_id"] }
  
  query <- paste("UPDATE ",SNPsummary.table," SET nr_samples=0, samples='';",sep="")
  res <- dbGetQuery_E(con,query,TALK=FALSE)
  
  query <- paste("SELECT sample_id FROM ",sample.table,";",sep="")
  res <- dbGetQuery_E(con,query,TALK=FALSE)
  sampleIds <- res[,1]
  
  for(sampleId in sampleIds){
    
    cat("Processing sample ",sampleId," of ",length(sampleIds),"...\n",sep="")
    
    SNPdata.table <- paste("snp_data_",sampleId,sep="")
    variantsForSample <- list()
    SNPsToBeInserted <- list()
    
    SNPsToBeInserted[["nrSamples"]] <- array(NA,0)
    SNPsToBeInserted[["sampleStr"]] <- array(NA,0)
    
    variantsForSample[["sample_id"]] <- sampleId
    variantsForSample[["SNPs"]] <- NULL
    variantsForSample[["indels"]] <- NULL
    
    query <- paste("SELECT SNP_id FROM ",SNPdata.table,";",sep="")
    snpIds <- dbGetQuery_E(con,query,TALK=FALSE)
    variantsForSample[["SNPs"]] <- snpIds[,1]
    
    SNPids <- variantsForSample[["SNPs"]]
    
    sampleStr <- paste(",",sampleId,sep="")
    
    seenSNPs <- (names(SNPsToBeInserted[["nrSamples"]]) %in% SNPids)
    SNPsToBeInserted[["sampleStr"]][seenSNPs] <- paste(SNPsToBeInserted[["sampleStr"]][seenSNPs],sampleStr,sep="")
    SNPsToBeInserted[["nrSamples"]][seenSNPs] <- SNPsToBeInserted[["nrSamples"]][seenSNPs]+1
    
    newSNPids <- SNPids[!(SNPids %in% names(SNPsToBeInserted[["nrSamples"]]))]
    
    sampleStrsToBeAdded <- rep(as.character(sampleId), length(newSNPids))
    names(sampleStrsToBeAdded) <- newSNPids
    nrSamplesToBeAdded <- rep(1, length(newSNPids))
    names(nrSamplesToBeAdded) <- newSNPids
    SNPsToBeInserted[["sampleStr"]] <- c(SNPsToBeInserted[["sampleStr"]], sampleStrsToBeAdded)
    SNPsToBeInserted[["nrSamples"]] <- c(SNPsToBeInserted[["nrSamples"]], nrSamplesToBeAdded)
    
    SNPobjectsize <- object.size(SNPsToBeInserted)
    cat("Object sizes, SNPs:",SNPobjectsize,"\n")
    
    SNPsToBeInsertedTable <- cbind(names(SNPsToBeInserted[["nrSamples"]]), SNPsToBeInserted[["nrSamples"]], SNPsToBeInserted[["sampleStr"]])
    colnames(SNPsToBeInsertedTable) <- c("SNP_id","nr_samples","sample_str")
    
    totalNrToInsert <- nrow(SNPsToBeInsertedTable)
    
    ## Divide large data into smaller chunks, to make sure MySQL doesn't crash
    chunkStart <- 1
    chunkEnd <- min(chunkSize,totalNrToInsert)
    moreChunksToInsert <- TRUE
    
    while(moreChunksToInsert){
      
      SNPsToBeInsertedTableChunk <- SNPsToBeInsertedTable[chunkStart:chunkEnd,]
      
      alreadyInSummary <- (SNPsToBeInsertedTableChunk[,"SNP_id"] %in% SNPidsInSummary)
      SNPsToBeUpdated <- SNPsToBeInsertedTableChunk[alreadyInSummary,,drop=FALSE]
      SNPsToBeAdded <- SNPsToBeInsertedTableChunk[!alreadyInSummary,,drop=FALSE]
      SNPidsInSummary <- c(SNPidsInSummary,SNPsToBeAdded[,"SNP_id"])
      
      cat(" chunk:",chunkStart,"-",chunkEnd," total:",nrow(SNPsToBeInsertedTable)," to update:",nrow(SNPsToBeUpdated)," to add:",nrow(SNPsToBeAdded),proc.time()[3] - ps,"s\n");
      updateSNPsummaryTable(SNPsToBeUpdated, SNPsToBeAdded)
      
      if(chunkEnd<totalNrToInsert){
        chunkStart <- chunkEnd+1
        chunkEnd <-  min(chunkEnd+chunkSize,totalNrToInsert)
      }
      else{
        moreChunksToInsert <- FALSE
      }
    }
  }
  
  query <- "update snp_summary set samples=substr(samples,2);"
  res <- dbGetQuery_E(con,query,TALK=FALSE)
}

rebuildIndelSummary <- function(){
  
  ps <- proc.time()[3]
  
  con <- connectToInhouseDB()
  sample.table <- dbTables[["sample"]]
  Indelsummary.table <- dbTables[["indel.summary"]]
  chunkSize <- 500000
  
  sampleIds <- list()
  IndelidsInSummary <- NULL
  
  query <- paste("SELECT indel_id FROM ",Indelsummary.table,";")
  currIndelIds <- dbGetQuery_E(con,query,TALK=FALSE)
  if(nrow(currIndelIds)>0){ IndelIdsInSummary <- currIndelIds[,"indel_id"] }
  
  query <- paste("UPDATE ",Indelsummary.table," SET nr_samples=0, samples='';",sep="")
  res <- dbGetQuery_E(con,query,TALK=FALSE)
  
  query <- paste("SELECT sample_id FROM ",sample.table,";",sep="")
  res <- dbGetQuery_E(con,query,TALK=FALSE)
  sampleIds <- res[,1]
  
  for(sampleId in sampleIds){
    
    cat("Processing sample ",sampleId," of ",length(sampleIds),"...\n",sep="")
    
    IndelData.table <- paste("indel_data_",sampleId,sep="")
    variantsForSample <- list()
    IndelsToBeInserted <- list()
    
    IndelsToBeInserted[["nrSamples"]] <- array(NA,0)
    IndelsToBeInserted[["sampleStr"]] <- array(NA,0)
    
    variantsForSample[["sample_id"]] <- sampleId
    variantsForSample[["Indels"]] <- NULL
    variantsForSample[["indels"]] <- NULL
    
    query <- paste("SELECT indel_id FROM ",IndelData.table,";",sep="")
    indelIds <- dbGetQuery_E(con,query,TALK=FALSE)
    variantsForSample[["indels"]] <- indelIds[,1]
    
    IndelIds <- variantsForSample[["indels"]]
    
    sampleStr <- paste(",",sampleId,sep="")
    
    seenIndels <- (names(IndelsToBeInserted[["nrSamples"]]) %in% IndelIds)
    IndelsToBeInserted[["sampleStr"]][seenIndels] <- paste(IndelsToBeInserted[["sampleStr"]][seenIndels],sampleStr,sep="")
    IndelsToBeInserted[["nrSamples"]][seenIndels] <- IndelsToBeInserted[["nrSamples"]][seenIndels]+1
    
    newIndelIds <- IndelIds[!(IndelIds %in% names(IndelsToBeInserted[["nrSamples"]]))]
    
    sampleStrsToBeAdded <- rep(as.character(sampleId), length(newIndelIds))
    names(sampleStrsToBeAdded) <- newIndelIds
    nrSamplesToBeAdded <- rep(1, length(newIndelIds))
    names(nrSamplesToBeAdded) <- newIndelIds
    IndelsToBeInserted[["sampleStr"]] <- c(IndelsToBeInserted[["sampleStr"]], sampleStrsToBeAdded)
    IndelsToBeInserted[["nrSamples"]] <- c(IndelsToBeInserted[["nrSamples"]], nrSamplesToBeAdded)
    
    IndelObjectsize <- object.size(IndelsToBeInserted)
    cat("Object sizes, Indels:",IndelObjectsize,"\n")
    
    IndelsToBeInsertedTable <- cbind(names(IndelsToBeInserted[["nrSamples"]]), IndelsToBeInserted[["nrSamples"]], IndelsToBeInserted[["sampleStr"]])
    colnames(IndelsToBeInsertedTable) <- c("indel_id","nr_samples","sample_str")
    
    totalNrToInsert <- nrow(IndelsToBeInsertedTable)
    
    ## Divide large data into smaller chunks, to make sure MySQL doesn't crash
    chunkStart <- 1
    chunkEnd <- min(chunkSize,totalNrToInsert)
    moreChunksToInsert <- TRUE
    
    while(moreChunksToInsert){
      
      IndelsToBeInsertedTableChunk <- IndelsToBeInsertedTable[chunkStart:chunkEnd,]
      
      alreadyInSummary <- (IndelsToBeInsertedTableChunk[,"indel_id"] %in% IndelIdsInSummary)
      IndelsToBeUpdated <- IndelsToBeInsertedTableChunk[alreadyInSummary,,drop=FALSE]
      IndelsToBeAdded <- IndelsToBeInsertedTableChunk[!alreadyInSummary,,drop=FALSE]
      IndelIdsInSummary <- c(IndelIdsInSummary,IndelsToBeAdded[,"indel_id"])
      
      cat(" chunk:",chunkStart,"-",chunkEnd," total:",nrow(IndelsToBeInsertedTable)," to update:",nrow(IndelsToBeUpdated)," to add:",nrow(IndelsToBeAdded),proc.time()[3] - ps,"s\n");
      updateIndelSummaryTable(IndelsToBeUpdated, IndelsToBeAdded)
      
      if(chunkEnd<totalNrToInsert){
        chunkStart <- chunkEnd+1
        chunkEnd <-  min(chunkEnd+chunkSize,totalNrToInsert)
      }
      else{
        moreChunksToInsert <- FALSE
      }
    }
  }
  
  query <- "update indel_summary set samples=substr(samples,2);"
  res <- dbGetQuery_E(con,query,TALK=FALSE)
}