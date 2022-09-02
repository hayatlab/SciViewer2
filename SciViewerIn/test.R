rm(list=ls())
setwd("~/scripts/SciViewerDev/SciViewerIn/")
gc()

#RSQLite::dbWriteTable(connSc,paste0(StudyName,"_SignatureMatrix1"),m,overwrite=T)
#q <- RSQLite::dbGetQuery(connSc, paste0("SELECT * FROM RNA_abc_SignatureMatrix1"))

##--- code to reformat .gmt file; 
if(F){
  
  connGenes <- RSQLite::dbConnect(RSQLite::SQLite(),'C:/Dhawal/SHINYAPP_DATA/hs_genes.db')
  RSQLite::dbListTables(connGenes)
  query <- paste0("SELECT * FROM gAlias")
  genes <- RSQLite::dbGetQuery(connGenes, query)
  RSQLite::dbDisconnect(connGenes)
  rm(connGenes,query)
  
  out <- c()
  for(myfile in c("c2.cp.kegg.v7.5.1.symbols.gmt","c2.cp.wikipathways.v7.5.1.symbols.gmt")){
    d <- GSA::GSA.read.gmt(myfile)
    cf <- data.frame(geneset = d$geneset.names,
                     num_genes=NA,
                     genes=NA,
                     genes_with_alias=NA)
    for(i in 1:length(d$genesets)){
      cat(i,"\n")
      cf$num_genes[i] <- length(d$genesets[[i]])
      cf$genes[i] <- paste0(d$genesets[[i]],collapse = ",")
      myhgnc <- unique(genes[genes$geneSymbol%in%d$genesets[[i]],]$HGNC)
      cf$genes_with_alias[i] <- paste0(unique(genes[genes$HGNC%in%myhgnc,]$geneSymbol),collapse = ",")
      rm(myhgnc)
    }
    out <- rbind(out,cf)
    rm(d,cf)
  }
  
  write.table(out[1:50,],file = "msigdb_signatures.txt",quote = F,sep = "\t",row.names = F)
  
  
}

##--- geneal tests
if(F){
  source("input_functions.R",local = TRUE,encoding = "UTF-8")
  source("vars.R",local = TRUE,encoding = "UTF-8")
  
  ##-- from rds
  sa <- readRDS("/home/rstudio/data/SCS/TM_BoneMarrow_SmartSeq.rds")
  so <- local({
    si_assays<- Seurat::Assays(sa)
    so <- CreateSeuratObject(counts = sa@assays[[1]]@data,
                             assay = si_assays[1],
                             meta.data = sa@meta.data)
    if(length(si_assays)>1){
      message(' ..adding additional assay information to the Seurat object\n')
      for(i in 2:length(si_assays)){
        so[[si_assays[i]]] <- CreateAssayObject(counts = sa@assays[[i]]@data)
      }
    }
    scvis <- sa@reductions
    if(length(scvis)>1) message('Multiple data reductions found. Will populate them.\n')
    for(i in 1:length(scvis)){
      message(" ..adding coordinates: ",gsub("^X_","",names(scvis)[i]),"\n")
      so[[gsub("^X_","",names(scvis)[i])]] = CreateDimReducObject(embeddings = scvis[[i]]@cell.embeddings, key = gsub("^X_","",names(scvis)[i]))
    }
    so
  })
  so@active.ident <- factor(so$free_annotation)
  
  ##-- from h5ad file
  so <- local({
    ad <- anndata::read_h5ad(filename = "C:/Dhawal/scRNASeq_data/PCL/tabula-muris-senis-facs-processed-official-annotations-Aorta.h5ad")
    meta <- ad$obs
    meta$V1 <- meta$V2 <- NULL
    mat <- ad$X
    qq <- dimnames(mat)
    mat <- as(mat,'matrix.csr')
    mat <- as(mat,'dgCMatrix')
    dimnames(mat) <- qq
    rm(qq)
    mat <- t(mat)
    cat(dim(mat)[1],"\n")
    so = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0,meta.data = meta)
    scvis <- ad$obsm
    if(length(scvis)>1) message('Multiple data reductions found. Will populate them.\n')
    for(i in 1:length(scvis)){
      message(" ..adding coordinates: ",gsub("^X_","",names(scvis)[i]),"\n")
      so[[gsub("^X_","",names(scvis)[i])]] = CreateDimReducObject(embeddings = scvis[[i]], key = gsub("^X_","",names(scvis)[i]))
    }
    so
  })
  
  
  # subsample
  sampled.cells <- sample(x = colnames(so), 
                          size = 100, replace = F)
  so <- subset(so, cells = sampled.cells)
  rm(sampled.cells)
  
  ##-- meta
  meta <- local({
    coords <- data.frame(SAMPID=rownames(so@meta.data))
    for(f in names(so@reductions)){
      qq <- as.data.frame(so@reductions[[f]]@cell.embeddings)
      qq <- qq[,1:2]
      names(qq) <- c('V1','V2')
      names(qq) <- paste0(names(qq),'_',f)
      coords <- cbind(coords, qq)
      rm(qq)
    }
    my.meta <- so@meta.data
    my.meta <- cbind(my.meta,coords)
    my.meta
  })
  
  
  study <- data.frame(Database="abc",
                      Description=NA,
                      SampleSize=NA,
                      TISSUES=NA,
                      DISEASE_Variable=NA, 
                      PMID=NA,
                      GEO=NA,
                      STATUS='Public',
                      RATING='High',
                      CONTVAR=NA,
                      CATVAR='free_annotation')
  
  seurat2sqlite(so=so,
                si_study=study,
                si_reduction='umap',
                si_compute_cellmarker = F,
                si_cell_markers=NULL,
                si_cell_col = NULL,
                si_compute_diseasemarker = F,
                si_disease_markers=NULL,
                si_disease_col=NULL,
                si_celltype='free_annotation',
                si_donor = 'mouse.id',
                signature_file = NULL,
                db_address='asd',
                OUTDIR=NULL)
  
  
  si_celltype='free_annotation'
  si_donor = 'mouse.id'
  my_assay='RNA'
  signature_file <- read.delim("data/msigdb_signatures.txt",header=F,stringsAsFactors = F)
  si_reduction='umap'
  si_study=study
  si_compute_cellmarker = F
  si_cell_markers=NULL
  si_cell_col = NULL
  si_compute_diseasemarker = F
  si_disease_markers=NULL
  si_disease_col=NULL
  db_address='asd'
  OUTDIR=NULL

  
  ## quick checks
  DefaultAssay(so) = 'RNA'
  group.by = c(si_celltype,si_donor)
  
  m <- so@assays$RNA@data
  n <- m[,colnames(m)[1:20]]
  a <- apply(n,1,mean)
  
  m <- FetchData(so,vars = c(si_celltype,si_donor))
  names(m) <- c('cell_type','donor')
  m$id <- paste0(m$cell_type,"|",m$donor)
  split(m,m$id)

  
  m <- m[which(rowSums(is.na(m)) == 0), , drop = F]
  for (i in 1:ncol(m)) {
    m[, i] <- as.factor(m[, i])
  }
  num.levels <- sapply(
    X = 1:ncol(m),
    FUN = function(i) {
      length(levels(m[, i]))
    }
  )
  
  apply(so@assays$RNA@data,1,sd)
  
  
  category.matrix <- sparse.model.matrix(object = as.formula(
    object = paste0(
      '~0+',paste0("m[,",1:length(group.by),"]",collapse = ":"
      )
    )))
  colsums <- colSums(x = category.matrix)
  category.matrix <- category.matrix[, colsums > 0]
  colsums <- colsums[colsums > 0]
  colnames(category.matrix) <- sapply(
    X = colnames(x = category.matrix),
    FUN = function(name) {
      name <- gsub(pattern = "m\\[, [1-9]*\\]", replacement = "", x = name)
      return(paste0(rev(x = unlist(x = strsplit(x = name, split = ":"))), collapse = "_"))
    })
  
  
  m <- Seurat::AverageExpression(object = so,assays = my_assay,
                            group.by = c(si_celltype,si_donor))[[1]]
  
}


##--- test output
if(F){
  study = "RNA_HS_healthyCholangiocytes_10x"
  db_address = "/home/rstudio/data/SCS/scRNA_RNA_HS_healthyCholangiocytes_10x.db"
  connSc <- RSQLite::dbConnect(RSQLite::SQLite(),paste0(db_address))
  RSQLite::dbListTables(connSc)
  
  sc_study <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_study"))
  m <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_metaFeatures"))
  n <- RSQLite::dbGetQuery(conn = connSc,paste0('SELECT * FROM ',study,"_FeatureSummary"))
  
  
  #c("author_cell_type","donor_uuid","Phase","cell_type","disease",
  #  "development_stage","sample_preservation_method","sex","author_cluster")
  Categorical_Vars <- "author_cell_type,donor_uuid,Phase,cell_type,disease,development_stage,sample_preservation_method,sex,author_cluster"
  Celltype <- "author_cell_type"
  qq <- unlist(strsplit(as.character(Categorical_Vars),","))
  qq <- ifelse(qq==Celltype,as.character("SVcell_type"),as.character(qq))
  if(sum(qq%in%names(so@meta.data))!=length(qq) ){
    #cat(paste0(qq,collapse = ","),"\n",
    #    paste0(names(so@meta.data),","),"\n")
    stop(" could not find provided categorical variables in the metadata file\n") 
  }
  if(!'SVcell_type'%in%qq){
    stop(" could not find column with name 'cell_type' in the metadata file\n") 
  }
  Categorical_Vars = qq
  message('Categorical variables in the metadata: ',Categorical_Vars,"\n\n")
  
   
}