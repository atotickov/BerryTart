library(tidyverse)
library(splitstackshape)
library(ape)

readGFF <- function(gffFile, nrows = -1) {
    cat("Reading ", gffFile, ": ", sep="")
    gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
      header=FALSE, comment.char="&", nrows = nrows,
      
      # для результата TRF
      colClasses=c("character", "character", "character", "integer", 
                   "integer", "character", "character", "character", 
                   "character"),      
      
      stringsAsFactors=TRUE)
    
    # для trf
    colnames(gff) = c("chr_id", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
    
    cat("found", nrow(gff), "rows with classes:",
        paste(sapply(gff, class), collapse=", "), "\n")
    stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
    return(gff)
 }


gff2tab <- function(file, attrsep=";", fieldsep="=") {

    gff <- readGFF(file)

    ## parse attributes into table
    att <- strsplit(gff[,"attributes"], split = attrsep, fixed = TRUE)
 
    ## get all fields (attributes with = : strsplit returns 2 fields)
    fields <- unique(c(unlist(sapply(att, function(atts) {
        a <- strsplit(atts, split = fieldsep, fixed = TRUE)
        unlist(lapply(a, function(x) if (length(x)==2) x[1]))
    }))))
    attm <- matrix(NA, nrow=nrow(gff), ncol=length(fields))
    colnames(attm) <- fields
        
    
    ## fuse back fields without fieldsep
    att2 <- lapply(att, function(x) {
        a <- strsplit(x, split = fieldsep, fixed = TRUE)
        ## fuse back fields without =
        newa <- list()
        if ( length(a)>0 )
          for ( j in 1:length(a) ) {
              len <- length(newa)
              ## append if length is ok (= present)
              if ( length(a[[j]])==2 ) newa <- append(newa, list(a[[j]]))
              if ( length(a[[j]])==1 ) # or append to previous content
                newa[[len]][2] <- paste(newa[[len]][2],a[[j]],sep=";")
              if ( !length(a[[j]])%in%1:2) cat(paste("error\n"))
          }
        newa})
    for ( i in 1:nrow(attm) ) {
        key <- unlist(lapply(att2[[i]], function(y) y[1]))
        val <- unlist(lapply(att2[[i]], function(y) y[2]))
        attm[i,key] <- val
    }

    cbind.data.frame(gff[,colnames(gff)!="attributes"],attm,
                     stringsAsFactors=FALSE)
}


data <- readGFF('/mnt/tank/scratch/skliver/common/mustelidae/atotik/human/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.2.7.7.80.10.50.500.dat.gff.gz')
data_table <- gff2tab('/mnt/tank/scratch/skliver/common/mustelidae/atotik/human/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.2.7.7.80.10.50.500.dat.gff.gz')
caroline::write.delim(data_table, '/mnt/tank/scratch/skliver/common/mustelidae/atotik/human/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.2.7.7.80.10.50.500.dat.gff.gz.tab', quote = FALSE, row.names = FALSE, sep='\t')

