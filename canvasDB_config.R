############################################################################
##
## File: canvasDB_config.R
##
## Author: Adam Ameur, Uppsala Genome Center
##
## Description: Global variables and configurations needed for the canvasDB
##
############################################################################

tmpfileDir <- paste(rootDir,"tmpfile_dir/",sep="") # Directory for temporary files.

if(!file.exists(tmpfileDir)){
    dir.create(tmpfileDir)
}

#indbDatabaseName <- "canvasdb_test"
#annotationDatabaseName <- "canvasdb_test_annot_hg19"
indbDatabaseName <- "canvasdb"
annotationDatabaseName <- "canvasdb_annot_cf3"


#ANNOVARpath <- "/Volumes/Data/tools/annovar/annotate_variation.pl" ## Add path to ANNOVAR executable
#ANNOVARpathDB <- "/Volumes/Data/tools/annovar/humandb/" ## Add path to ANNOVAR annotation database
ANNOVARpath <- "/home/eschofield/local/src/annovar/annotate_variation.pl" ## Add path to ANNOVAR executable
ANNOVARpathDB <- "/home/eschofield/local/src/annovar/dogdb/" ## Add path to ANNOVAR annotation database

VEPpath <- "/opt/vep/vep"
snpEffPath <- "/home/eschofield/local/src/snpEff"
snpEffDB <- "CanFam3.1.86"
siftDB <- "CanFam3.1.83"

VEPpath <- "/Users/Ellen/Git/ensembl-vep/vep"
VEPpath <- paste("docker run -it -v /Users/Ellen/DockerData/vep_data:/opt/vep/.vep -v ", tmpfileDir, ":/data/ ensemblorg/ensembl-vep ./vep",sep="")

dbSNPversion <- 151
eVersion <- 95

mysqlConfigFile <- ".my.cnf"  # File containing MySQL user info

mysqlEngine <- "MyISAM"

## MySQL database tables
annotTables <- list()
annotTables[["score"]] <- paste(annotationDatabaseName,".annot_score",sep="")
annotTables[["dbsnp"]] <- paste(annotationDatabaseName,".annot_snp",dbSNPversion,"_single",sep="")
#annotTables[["dbsnpCommon"]] <- paste(annotationDatabaseName,".annot_snp",dbSNPversion,"_single_common",sep="")
annotTables[["dbsnpIndels"]] <- paste(annotationDatabaseName,".annot_snp",dbSNPversion,"_indels",sep="")
#annotTables[["dbsnpCommonIndels"]] <- paste(annotationDatabaseName,".annot_snp",dbSNPversion,"_indels_common",sep="")
annotTables[["vepSNP"]] <- paste(annotationDatabaseName,".annot_vep",eVersion,"_single",sep="")
annotTables[["vepIndels"]] <- paste(annotationDatabaseName,".annot_vep",eVersion,"_indels",sep="")


dbTables <- list()
dbTables[["run"]] <- "runs"
dbTables[["sample"]] <- "samples"
dbTables[["SNP.summary"]] <- "snp_summary"
dbTables[["indel.summary"]] <- "indel_summary"

## This should stop all forms of "Scientific writing, ie 1e+05"
options(scipen = 999)
