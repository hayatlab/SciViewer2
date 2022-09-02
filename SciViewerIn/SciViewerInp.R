#!/usr/bin/env Rscript
require("optparse")
require("getopt")
source("input_functions.R",local = TRUE,encoding = "UTF-8")
source("vars.R",local = TRUE,encoding = "UTF-8")

option_list = list(
  ##------- Input data files
  make_option(c("--h5ad"), type="character", default=NULL, 
              help="If the data file is in the H5AD format, then the file path [default= NULL]", metavar="character"),
  make_option(c("--seurat"), type="character", default=NULL, 
              help="If the data file is in the Seurat format, then the file path [default= NULL]", metavar="character"),
  
  ##------- Mandatory inputs
  make_option(c("--db_address"), type="character", default=NULL, 
              help="(*) full address where the database file is to be written [default= NULL]", metavar="character"),
  make_option(c("--StudyName"), type="character", default=NULL, 
              help="(*) name using which you would like to identify the study in the application [default= NULL]", metavar="character"),
  make_option(c("--celltypeCol"), type="character", default=NULL, 
              help="(*) name of cell type annotation column in the metadata file [default= NULL]", metavar="character"),
  make_option(c("--org"), type="character", default=NULL, 
              help="(*) organism name [default= NULL]", metavar="character"),
  make_option(c("--diseaseName"), type="character", default=NULL, 
              help="(*) disease term [default= NULL]", metavar="character"),
  make_option(c("--shotdescr"), type="character", default=NULL, 
              help="(*) short one-liner description of the study [default= NULL]", metavar="character"),
  
  ##------- Optional inputs from section-2
  make_option(c("--reduction"), type="character", default=NULL, 
              help="(*) Default coordinate system to display.  [default= NULL]", metavar="character"),
  make_option(c("--donorCol"), type="character", default=NULL, 
              help="Name of donor annotation column in the metadata file [default= NULL]", metavar="character"),
  
  ##------- Optional inputs from section-3
  make_option(c("--Continuous_Vars"), type="character", default=NULL, 
              help="continuous variables identified from the metadata file, this info will be used for data display [default= %default]", metavar="character"),
  make_option(c("--Categorical_Vars"), type="character", default=NULL, 
              help="categorical variables identified from the metadata file, this info will be used for data display. column name 'cell_type' is required [default= %default]", metavar="character"),
  
  ##------- Inputs for DE/Markers calculation, section-4
  make_option(c("--compute_markers"), type="logical", default=T, 
              help="whether to compute markers for cell type and any additional provided covariate [default= %default]", metavar="logical"),
  make_option(c("--covariates"), type="character", default=NULL, 
              help="Comma-separated list of covariates from the metadata file for which markers genes to be calculated [default= NULL]", metavar="character"),
  make_option(c("--marker_file"), type="character", default=NULL, 
              help="Precomputed and formatted text file path containing marker information [default= NULL]", metavar="character"),
  make_option(c("--compute_deg"), type="logical", default=F, 
              help="whether to compute differential gene expression [default= %default]", metavar="logical"),
  make_option(c("--disease_var"), type="character", default=NULL, 
              help="Column header for the test variable [default= NULL]", metavar="character"),
  make_option(c("--deg_file"), type="character", default=NULL, 
              help="Precomputed and formatted text file path containing differentially expressed genes between conditions [default= NULL]", metavar="character"),
  
  ##------- Study details, section-5.1
  make_option(c("--StudyDescr"), type="character", default=NULL, 
              help="detailed description of the study, this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("--Tissue"), type="character", default=NULL, 
              help="tissue used in the study, this info will be reflected in the summary section", metavar="character"),
  make_option(c("--PMID"), type="character", default=NULL, 
              help="publication id if available, this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("--GEO"), type="character", default=NULL, 
              help="data repository link if available, this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("--StudyStatus"), type="character", default='Internal', 
              help="status of the study (e.g. public/internal), this info will be reflected in the summary section [default= %default]", metavar="character"),
  make_option(c("--StudyRating"), type="character", default='High', 
              help="your rating for the study (e.g. high/poor/average), this info will be reflected in the summary section [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


###--------------------- Input tests
if(is.null(opt$db_address)){
  stop('path to write database is not provided!')
}
if(is.null(opt$celltypeCol)){
  stop('column name for cell type annotations not provided!')
}
if(is.null(opt$StudyName)){
  stop('name of the study to be used for writing the tables in the database is not specified!')
}

###--------------------- H5AD to data objects
if(!is.null(opt$h5ad)){
  so <- local({
    require(anndata)
    require(dplyr)
    require(Matrix)
    h5adpath=opt$h5ad
    #h5adpath="/gpfs01/home/glxaw/data/scRNASeq_datasets/Tabula_Muris/tabula-muris-senis-facs-processed-official-annotations-Aorta.h5ad"
    #h5adpath='/gpfs01/home/glxaw/data/scRNASeq_datasets/Barbry_HealthyAirways/Barbry_HealthyAirways_Lung_fromDhawal_mahmoudConvert_seuratDisk.h5ad'
    ad <- anndata::read_h5ad(filename = h5adpath)
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
  OUTDIR = dirname(opt$h5ad)

}else if(!is.null(opt$seurat)){
   so <- local({
     #-- this part simplifies the Seurat object and reduces size, removes version incompatibility etc
     sa <- readRDS(file = opt$seurat)
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
   OUTDIR = dirname(opt$seurat)
}else{
   stop(' provided input format is not supported!')
}

##--- signature file
signature_file=NULL
if(file.exists("data/msigdb_signatures.txt")){
  signature_file <- read.delim("data/msigdb_signatures.txt",header=F,stringsAsFactors = F)
  if(sum(c("V1", "V2", "V3")%in%names(signature_file))==3){
    message(" ..gene signature file exits, will be used for computing the signatures\n")
  }else{
    message(" ..gene signature file exits, but not correctly formatted. Skipping signature calculations.\n")
  }
}else{
  message(" ..gene signature file does not exits. Skipping signature calculations.\n")
}

# precomputed marker file
if(!is.null(opt$marker_file)){
  Marker_Precomputed <- local({
    ext <- opt$marker_file
    validate(need(ext %in% c("txt","csv"), "Please upload txt/csv file!"))
    cf <- NULL
    if(ext=='txt'){
      cf = read.delim(opt$marker_file,as.is = T,colClasses = "#",sep="\t")
    }else if(ext == 'csv'){
      cf = read.csv2(opt$marker_file,as.is = T,colClasses = "#")
    }
    cf
  })
}else{
  Marker_Precomputed = NULL
}

# precomputed de file
if(!is.null(opt$deg_file)){
  DE_Precomputed <- local({
    ext <- opt$deg_file
    validate(need(ext %in% c("txt","csv"), "Please upload txt/csv file!"))
    cf <- NULL
    if(ext=='txt'){
      cf = read.delim(opt$deg_file,as.is = T,colClasses = "#",sep="\t")
    }else if(ext == 'csv'){
      cf = read.csv2(opt$deg_file,as.is = T,colClasses = "#")
    }
    cf
  })
}else{
  DE_Precomputed = NULL
}




###--------------------- generate database
write.scstudy2.sqlitedb(so = so,
                        db_address=opt$db_address,
                        StudyName=opt$StudyName,
                        Donors_VariableName = opt$donorCol,
                        Celltype = opt$celltypeCol,
                        Reduction_map = opt$reduction,
                        Organism=opt$org,
                        Disease=opt$diseaseName,
                        ShortDescr=opt$shotdescr,
                        
                        Marker_Calc = opt$compute_markers,
                        Marker_Covariates = opt$covariates,
                        Marker_Precomputed=Marker_Precomputed,
                        
                        DE_Calc = opt$compute_deg,
                        Disease_VariableName=opt$disease_var,
                        DE_Precomputed=DE_Precomputed,
                        
                        StudyDescr=opt$StudyDescr,
                        Tissue=opt$Tissue,
                        PMID=opt$PMID,
                        GEO=opt$GEO,
                        StudyStatus=opt$StudyStatus,
                        StudyRating=opt$StudyRating,
                        
                        Continuous_Vars=opt$Continuous_Vars,
                        Categorical_Vars=opt$Categorical_Vars,
                        
                        OUTDIR = OUTDIR,
                        signature_file=signature_file)






