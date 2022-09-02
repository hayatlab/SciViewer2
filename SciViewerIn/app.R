source("input_functions.R",local = TRUE,encoding = "UTF-8")
source("vars.R",local = TRUE,encoding = "UTF-8")


##----------------- extra code
# javascript code to collapse box
jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"


##----------------- application
server <- function(input, output,session) {
  VARS <- reactiveValues(
    roots = ROOTS,
    outdir=NULL
  )
  observeEvent({input$ab573dfg3},{
    if(isTruthy(input$ab573dfg3)){
      if(dir.exists(input$datadir)){
        VARS$roots = ROOTS
        inputdir=ifelse(isTRUE(grep("/$",input$datadir)),
                        as.character(input$datadir),
                        as.character(paste0(input$datadir,"/")))
        names(inputdir) = 'inputdir'
        cat("inputdir:",inputdir,"\n")
        VARS$roots = append(inputdir,VARS$roots)
        output$direxists <- NULL
        rm(inputdir)
      }else{
        output$direxists <- renderPrint({ cat('Directory does not exits!\n') })
      }
    }
  })
  observeEvent({VARS$roots},{
    if(isTruthy(VARS$roots)){
      shinyFileChoose(input,'infile', session=session,roots=VARS$roots, 
                      filetypes=c('', 'rds','h5ad'))
    }
  })
  
  variables = reactiveValues(
    so = NULL,
    redmap = NULL,
    inptfile = NULL,
    si_markerfile_path = NULL,
    si_defile_path = NULL,
    si_markerfile = NULL,
    si_defile = NULL,
    simarker = NULL,
    mat = NULL,
    meta = NULL,
    db = NULL,
    studyname = NULL,
    svorg = NULL,
    svdisease = NULL,
    shotdescr = NULL,
    descr = NULL,
    tissue = NULL,
    pmid = NULL,
    geo = NULL,
    status = NULL,
    rating = NULL,
    run = NULL,
    cell_type = NULL,
    barcode = NULL,
    clusterx = NULL,
    clustery = NULL,
    donor = NULL,
    disease = NULL,
    de=NULL,
    covariates = NULL,
    contvars = NULL,
    catvars = NULL
  )
  observeEvent({input$infile},{
    inFile <- parseFilePaths(roots=VARS$roots, input$infile)
    if(isTruthy(inFile$datapath)){
      variables$so = NULL
      variables$redmap = NULL
      variables$inptfile = NULL
      variables$si_markerfile_path = NULL
      variables$si_defile_path = NULL
      variables$si_markerfile = NULL
      variables$si_defile = NULL
      variables$simarker = NULL
      variables$mat = NULL
      variables$meta = NULL 
      variables$db = NULL
      variables$studyname = NULL
      variables$svorg = NULL
      variables$svdisease = NULL
      variables$shotdescr = NULL
      variables$descr = NULL
      variables$tissue = NULL
      variables$pmid = NULL
      variables$geo = NULL
      variables$status = NULL
      variables$rating = NULL
      variables$run = NULL
      variables$cell_type = NULL
      variables$barcode = NULL
      variables$clusterx = NULL
      variables$clustery = NULL
      variables$donor = NULL
      variables$de = NULL
      variables$disease = NULL
      variables$covariates = NULL
      variables$contvars = NULL
      variables$catvars = NULL
      
      ext <- tools::file_ext(unname(inFile$datapath))
      VARS$outdir <- dirname(unname(inFile$datapath))
      message("observed file extention: ",ext,"\n")
      validate(need(ext %in% c("h5ad","rds"), "Please upload a correct file!"))
      output$direxists <- renderPrint({ cat('outdir: ',VARS$outdir) })
      so = NULL
      if(ext=='rds'){
        scProg <- shiny::Progress$new()
        on.exit(scProg$close())
        scProg$set(message = "reading the file..", value = 0)
        so <- local({
          #-- this part simplifies the Seurat object and reduces size, removes version incompatibility etc
          sa = readRDS(unname(inFile$datapath))
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
        variables$inptfile <- paste0(' --seurat \"',unname(inFile$datapath),"\"")
        scProg$set(message = "done!", value = 1)
      }else if(ext =='h5ad'){
        scProg <- shiny::Progress$new()
        on.exit(scProg$close())
        scProg$set(message = "reading the file..", value = 0)
        so <- local({
          ad <- anndata::read_h5ad(filename = unname(inFile$datapath))
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
        variables$inptfile <- paste0(' --h5ad \"',unname(inFile$datapath),"\"")
        scProg$set(message = "done!", value = 1)
      }else{
        stop('Input data is not the correct format.\n')
      }
      if(!is.null(so)){
        variables$meta <- local({
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
        output$redmap <- renderUI({
          selectInput("redmap",label="Default reduction map to display",choices = names(so@reductions),multiple = F)
        })
        variables$so = so
      }
    }
  })
  
  observeEvent({input$infile_marker},{
    inFile <- parseFilePaths(roots=VARS$roots, input$infile_marker)
    ext <- tools::file_ext(unname(inFile$datapath))
    validate(need(ext %in% c("txt","csv"), "Please upload txt/csv file!"))
    variables$si_markerfile_path = unname(inFile$datapath)
    cf <- NULL
    if(ext=='txt'){
      cf = read.delim(unname(inFile$datapath),as.is = T,colClasses = "#",sep="\t")
    }else if(ext == 'csv'){
      cf = read.csv2(unname(inFile$datapath),as.is = T,colClasses = "#")
    }
    variables$si_markerfile = cf
  })
  
  observeEvent({input$infile_de},{
    inFile <- parseFilePaths(roots=VARS$roots, input$infile_de)
    ext <- tools::file_ext(unname(inFile$datapath))
    validate(need(ext %in% c("txt","csv"), "Please upload txt/csv file!"))
    variables$si_defile_path = unname(inFile$datapath)
    cf <- NULL
    if(ext=='txt'){
      cf = read.delim(unname(inFile$datapath),as.is = T,colClasses = "#",sep="\t")
    }else if(ext == 'csv'){
      cf = read.csv2(unname(inFile$datapath),as.is = T,colClasses = "#")
    }
    variables$si_defile = cf
  })
  
  
  output$vars <- renderUI({
    orderInput('variables', 'Available variables in the meta-data file', items = colnames(variables$meta),
               item_class = 'info',
               connect = c('cell_type','donor','disease','covariates','clusterx','clustery')) #'barcode',
  })
  output$metatable = DT::renderDataTable({
    DT::datatable(variables$meta,rownames = F,
                  options = list(searching = TRUE,pageLength = 5,
                                 scrollX = TRUE,scrollCollapse = TRUE))
  })
  output$metasummary = renderPrint({
    summary(variables$meta)
  })
  output$cell_type <- renderUI({
    orderInput('cell_type', 'Cell cluster annotation column', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables')
  })
  output$clusterx <- renderUI({
    orderInput('clusterx', 'X-coordinate of reduction method', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables')
  })
  output$clustery <- renderUI({
    orderInput('clustery', 'Y-coordinate of reduction method', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables')
  })
  
  output$vars1 <- renderUI({
    qq <- colnames(variables$meta)
    if(isTruthy(input$cell_type_order[[2]])){
      qq <- qq[qq!=input$cell_type_order[[2]]]
    }
    if(isTruthy(input$clusterx_order[[2]])){
      qq <- qq[qq!=input$clusterx_order[[2]]]
    }
    if(isTruthy(input$clustery_order[[2]])){
      qq <- qq[qq!=input$clustery_order[[2]]]
    }
    orderInput('variables1', 'Available variables in the meta-data file', items = qq,
               item_class = 'warning',
               connect = c('contvars','catvars'))
  })
  output$contvars <- renderUI({
    orderInput('contvars', 'Continuous variables in the data, that you want to use in the display', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables1')
  })
  output$catvars <- renderUI({
    orderInput('catvars', 'Categorical variables in the data, that you want to use in the display', items = input$cell_type_order[[2]], placeholder = 'Drag item/s here...', connect = 'variables1')
  })
  
  output$vars2 <- renderUI({
    orderInput('variables2', 'Available variables in the meta-data file', items = colnames(variables$meta),
               item_class = 'info',
               connect = c('disease','donor','covariates')) #'barcode',
  })
  output$donor <- renderUI({
    orderInput('donor', 'Donor annotation column', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables2')
  })
  output$disease <- renderUI({
    orderInput('disease', 'Disease/Test/Experiment annotation column', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables2')
  })
  output$covariates <- renderUI({
    orderInput('covariates', 'Variables for computing marker genes', items = NULL, placeholder = 'Drag item/s here...', connect = 'variables2')
  })
  
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
  
  
  ##--------- display command
  observeEvent({input$cmd},{
    if(isTruthy(input$cmd)){
      
      if(isTruthy(input$db)) { variables$db = input$db}
      if(isTruthy(input$studyname)) { variables$studyname = input$studyname}
      if(isTruthy(input$redmap)) { variables$redmap = input$redmap}
      if(isTruthy(input$descr)) { variables$descr = input$descr}
      if(isTruthy(input$tissue)) { variables$tissue = input$tissue}
      if(isTruthy(input$svorg)) { variables$svorg = input$svorg}
      if(isTruthy(input$svdisease)) { variables$svdisease = input$svdisease}
      if(isTruthy(input$shotdescr)) { variables$shotdescr = input$shotdescr}
      if(isTruthy(input$pmid)) { variables$pmid = input$pmid}
      if(isTruthy(input$geo)) { variables$geo = input$geo}
      if(isTruthy(input$status)) { variables$status = input$status}
      if(isTruthy(input$rating)) { variables$rating = input$rating}
      if(isTruthy(input$de)) { variables$de = input$de}
      if(isTruthy(input$simarker)) { variables$simarker = input$simarker}
      if(isTruthy(input$cell_type_order[[2]])) { variables$cell_type = input$cell_type_order[[2]][1]}
      if(isTruthy(input$donor_order[[2]])) { variables$donor = input$donor_order[[2]][1]}
      if(isTruthy(input$disease_order[[2]])) { variables$disease = input$disease_order[[2]][1]}
      if(isTruthy(input$covariates_order[[2]])) { variables$covariates = paste0(input$covariates_order[[2]],collapse = ",")}
      if(isTruthy(input$contvars_order[[2]])) { variables$contvars = paste0(input$contvars_order[[2]],collapse = ",")}
      if(isTruthy(input$catvars_order[[2]])) { variables$catvars = paste0(input$catvars_order[[2]],collapse = ",")}
      
      message("... 01: categorical variables before running functions: ",variables$catvars,"\n")
      
      output$cmd_out <- renderText({
        paste0("Rscript --vanilla SciViewerInp.R ",variables$inptfile,"
               --db_address \"",variables$db ,"\"
               --StudyName \"",variables$studyname,"\"
               --celltypeCol \"",variables$cell_type,"\" 
               --reduction \"",variables$redmap,"\"
               --donorCol \"",variables$donor,"\" 
               --Continuous_Vars \"",variables$contvars,"\" 
               --Categorical_Vars \"",variables$catvars,"\" 
               --compute_markers \"",variables$simarker,"\"
               --covariates \"",variables$covariates,"\" 
               --marker_file \"",variables$si_markerfile_path,"\" 
               --compute_deg \"",variables$de,"\" 
               --disease_var \"",variables$disease,"\" 
               --deg_file \"",variables$si_defile_path,"\"
               --StudyDescr \"",variables$descr,"\" 
               --Tissue \"",variables$tissue,"\" 
               --org \"",variables$svorg,"\"
               --diseaseName \"",variables$svdisease,"\"
               --shotdescr \"",variables$shotdescr,"\"
               --PMID \"",variables$pmid,"\" 
               --GEO \"",variables$geo,"\" 
               --StudyStatus \"",variables$status,"\" 
               --StudyRating \"",variables$rating,"\"")
      })
    }
  })
  
  ##--------- calculations
  observeEvent({input$run},{
    if(isTruthy(input$run)){
      
      if(isTruthy(input$db)) { variables$db = input$db}
      if(isTruthy(input$studyname)) { variables$studyname = input$studyname}
      if(isTruthy(input$redmap)) { variables$redmap = input$redmap}
      if(isTruthy(input$descr)) { variables$descr = input$descr}
      if(isTruthy(input$tissue)) { variables$tissue = input$tissue}
      if(isTruthy(input$svorg)) { variables$svorg = input$svorg}
      if(isTruthy(input$svdisease)) { variables$svdisease = input$svdisease}
      if(isTruthy(input$shotdescr)) { variables$shotdescr = input$shotdescr}
      if(isTruthy(input$pmid)) { variables$pmid = input$pmid}
      if(isTruthy(input$geo)) { variables$geo = input$geo}
      if(isTruthy(input$status)) { variables$status = input$status}
      if(isTruthy(input$rating)) { variables$rating = input$rating}
      if(isTruthy(input$de)) { variables$de = input$de}
      if(isTruthy(input$simarker)) { variables$simarker = input$simarker}
      if(isTruthy(input$cell_type_order[[2]])) { variables$cell_type = input$cell_type_order[[2]][1]}
      if(isTruthy(input$clusterx_order[[2]])) { variables$clusterx = input$clusterx_order[[2]][1]}
      if(isTruthy(input$clustery_order[[2]])) { variables$clustery = input$clustery_order[[2]][1]}
      if(isTruthy(input$donor_order[[2]])) { variables$donor = input$donor_order[[2]][1]}
      if(isTruthy(input$disease_order[[2]])) { variables$disease = input$disease_order[[2]][1]}
      if(isTruthy(input$covariates_order[[2]])) { variables$covariates = paste0(input$covariates_order[[2]],collapse = ",")}
      if(isTruthy(input$contvars_order[[2]])) { variables$contvars = paste0(input$contvars_order[[2]],collapse = ",")}
      if(isTruthy(input$catvars_order[[2]])) { variables$catvars = paste0(input$catvars_order[[2]],collapse = ",")}
      
      message("... 02: categorical variables before running functions: ",variables$catvars,"\n")
      
      output$finalSel <- renderText({
        paste0("<b> Output database: </b> ",variables$db," </font><br>",
               "<b> Study identifier: </b> ",variables$studyname," </font><br>",
               "<b> Organism: </b> ",variables$svorg," </font><br>",
               "<b> Disease: </b> ",variables$svdisease," </font><br>",
               "<b> One-liner description: </b> ",variables$shotdescr," </font><br>",
               "<b> Description: </b> ",variables$descr," </font><br>",
               "<b> Tissue source: </b> ",variables$tissue," </font><br>",
               "<b> PMID: </b> ",variables$pmid," </font><br>",
               "<b> GEO: </b> ",variables$geo," </font><br>",
               "<b> Availability status: </b> ",variables$status," </font><br>",
               "<b> User rating: </b> ",variables$rating," </font><br>",
               "<b> Whether to log normalize the data: </b> ",variables$lognorm," </font><br>",
               "<b> ColumnID for cell types: </b> ",variables$cell_type," </font><br>",
               "<b> ColumnID for cell barcodes: </b> ",variables$barcode," </font><br>",
               "<b> ColumnID for donor: </b> ",variables$donor," </font><br>",
               "<b> ColumnID for disease: </b> ",variables$disease," </font><br>",
               "<b> ColumnID for X/Y cluster coordinates: </b> ",variables$clusterx,",",variables$clustery," </font><br>",
               "<b> ColumnID for covariates in DE modeling: </b> ",variables$covariates," </font><br>",
               "<b> ColumnID for continuous variables: </b> ",variables$contvars," </font><br>",
               "<b> ColumnID for categorical variables: </b> ",variables$catvars," </font><br>"
        )
      })
      
      ##------------------------
      write.scstudy2.sqlitedb(so = variables$so,
                              db_address=variables$db,
                              StudyName=variables$studyname,
                              Celltype = variables$cell_type,
                              Reduction_map = variables$redmap,
                              Donors_VariableName = variables$donor,
                              DE_Calc = variables$de,
                              DE_Precomputed=variables$si_defile,
                              Disease_VariableName=variables$disease,
                              Marker_Calc = variables$simarker,
                              Marker_Covariates = variables$covariates,
                              Marker_Precomputed=variables$si_markerfile,
                              StudyDescr=variables$descr,
                              Tissue=variables$tissue,
                              Organism=variables$svorg,
                              Disease=variables$svdisease,
                              ShortDescr=variables$shotdescr,
                              PMID=variables$pmid,
                              GEO=variables$geo,
                              StudyStatus=variables$status,
                              StudyRating=variables$rating,
                              Continuous_Vars=variables$contvars,
                              Categorical_Vars=variables$catvars,
                              OUTDIR = VARS$outdir,
                              signature_file=signature_file)
      tryCatch({
        sendSweetAlert(session = session,
                       title = "DONE",
                       text = "Generated SQLite database file for provided study!",
                       type = "success")
      },error = function(e){
        sendSweetAlert(session = session,
                       title = "Failed",
                       text = "Unable to generate database file at this point!",
                       type = "error")
      })
      
      
    }
  })
}

ui <- fluidPage(
  shiny::HTML(
    "<div style = 'background-color:#ffffff;color: #6a51a3;font-size:30px;font-weight:bold; vertical-align:middle'>
     <img src = 'logo.png' align = 'left'  height = '55px' width = '200px'>
        Single-cell Interactive Viewer: Data Input 
     </div>"
  ),
  tags$style(HTML("
           .box.box-solid.box-primary>.box-header {
           color:#fff;
           background:#e5f5e0
          }
          .box.box-solid.box-primary{
           border-bottom-color:#a1d99b;
           border-left-color:#a1d99b;
           border-right-color:#a1d99b;
           border-top-color:#a1d99b;
           background:#e5f5e0
          }")
  ),
  tags$hr(),
  tags$br(),
  h2('1.Input file'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  shiny::HTML('<i>The application accepts scRNA-seq files in .h5ad or .rds formats.</i><br> 
              <b><i>.h5ad file: </b></i> should contain at least the RNA expression layer with metadata and cluster coordinate information.<br>
              <b><i>.rds file: </b></i> this is a seurat object with count data, metadata and cluster coordinates<br>'),
  tags$br(),
  
  shinyFilesButton('infile', 'File Select', 'Please select a dataset', FALSE),
  actionButton("filenotfound", "Can't locate the file",icon =icon("question-circle"),width = '300px',class = "btn-warning"),
  
  conditionalPanel(
    condition = "input.filenotfound>'0'",
    box(id="sdf",title = "", 
        solidHeader = T,collapsed = F,width = 12,status = 'primary',collapsible = T,
        shiny::HTML('<i>please provide the directory path where .h5ad/.rds files are available. this path is added to the file browser option above.<br>
              if left blank, following directories are available for browsing as default (<b>/home </b>; </b> .</b>; <b>C:/</b>  ).</i><br>'),
        fluidRow(column(width=8,textInput("datadir", "directory path", NULL)),
                 column(width=2,tags$br(),actionButton("ab573dfg3","Go!",icon = icon("play-circle")))
        ),
        verbatimTextOutput('direxists')
    )
  ),
  
  tags$br(),
  
  tags$br(),
  tags$br(),
  
  h2('2. Sneakpeak into the metadata'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  conditionalPanel(
    condition = "input.infile>'0'",
    box(id="m1",title = "metadata table", 
        solidHeader = T,collapsed = T,width = 12,status = 'primary',collapsible = T,
        DT::dataTableOutput("metatable")
    ),
    box(id="m2",title = " metadata summary", 
        solidHeader = T,collapsed = T,width = 12,status = 'primary',collapsible = T,
        verbatimTextOutput('metasummary')
    ),
    tags$br(),
    tags$br(),
  ),
  tags$br(),
  
  h2('3.Defining variables'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  uiOutput("vars"), 
  tags$br(),
  fluidRow(column(width = 3,uiOutput("cell_type"),tags$br()),
           column(width = 3,uiOutput("redmap"),tags$br()),
           column(width = 3,uiOutput("donor"),tags$br()),
           column(width = 3)
           #column(width = 3),
           #column(width = 3,  uiOutput("clustery"),tags$br())
  ),
  tags$br(),
  
  h4('3.1.Continuous/Categorical variables'),
  tags$hr(style = "border-top: 1px solid #d9d9d9;"),
  shiny::HTML('<i> continuous variables are grouped into 5 distinct groups. while, for each categorical variable arepresented by distinct colors.</i>'),
  uiOutput("vars1"),
  tags$br(),
  fluidRow(column(width = 6,uiOutput("contvars"),tags$br()),
           column(width = 6,uiOutput("catvars"),tags$br())
  ),
  tags$br(),
  
  
  h2('4.Calculating markers/ differential expression'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  shiny::HTML('<i>please note calculating the markers or differentially expressed genes can take a significant amount of time.
              </i><br><br>'),
  tags$br(),
  uiOutput("vars2"), 
  tags$br(),
  shiny::HTML("<i>Markers are computed using <u>FindMarkers</u> function from Seurat package.
                           The default non-parametric Wilcox test is used. The cells can be grouped by cell type, disease or any other covariate/s.
                           The provided function will group the cells based on a selected covariate and compute differences against the rest of the cells.</i><br>"),
  fluidRow(column(width = 4,  
                  selectInput("simarker",label="Compute Marker genes",choices = c('TRUE','FALSE'),multiple = F)
  ),
  column(width = 4,
         conditionalPanel(
           condition = "input.simarker == 'TRUE'",
           uiOutput("covariates"),
           shiny::HTML("<br><br>")
         )
  ),
  column(width = 4)
  ),
  conditionalPanel(
    condition = "input.simarker == 'FALSE'",
    shiny::HTML("<i> Do you consider uploading the user defined markers? If yes, please use the file link below to upload the marker data. 
                 You can also <b> skip</b> this step.<br><br>
                 The file should contain following mandatory headers (order does not matter):<u>logFC,AveExpr,t,P.Value,adj.P.Val,B,Test,Tag,geneSymbol</u><br>
                 Test represents the comparison label for the logFC, Tag represents the cell type. If one/more columns are not available, please add NA values for those columns. 
                 </i><br>"),
    
    
    shinyFilesButton('infile_marker', 'File Select', 'Please select pre-computed result file', TRUE),
    tags$br()
  ),
  shiny::HTML("<br><br>"),
  
  shiny::HTML("<i>Differentially expressed genes are computed using <u>FindMarkers</u> function from Seurat package.
                           The default non-parametric Wilcox test is used. The cells can be grouped by cell type, disease or any other covariate/s.
                           The provided function will group the cells based on a selected cell types and compute differences using the user-provided contrast.</i><br>"),
  fluidRow(column(width = 4,  
                  selectInput("de",label="Compute differential expression",choices = c('TRUE','FALSE'),multiple = F)
  ),
  column(width = 4,
         conditionalPanel(
           condition = "input.de == 'TRUE'",
           uiOutput("disease"),
         )
  ),
  column(width = 4)
  ),
  conditionalPanel(
    condition = "input.de == 'FALSE'",
    shiny::HTML("<i> Do you consider uploading the user defined differential genes? If yes, please use the file link below to upload the data. 
                 You can also <b> skip</b> this step.<br><br>
                 The file should contain following mandatory headers (order does not matter):<u>logFC,AveExpr,t,P.Value,adj.P.Val,B,Test,Tag,geneSymbol</u><br>
                 Test represents the comparison label for the logFC, Tag represents the cell type. If one/more columns are not available, please add NA values for those columns. 
                 </i><br>"),
    shinyFilesButton('infile_de', 'File Select', 'Please select pre-computed result file', TRUE),
    tags$br()
    #fluidRow(column(width = 4,uiOutput("donor"),tags$br()),
    #         column(width = 4,uiOutput("disease"),tags$br()),
    #         column(width = 4,uiOutput("covariates"),tags$br())
    #)
  ),
  #verbatimTextOutput('covariatesout'),
  tags$br(),
  
  
  h2('5.Study details'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  fluidRow(column(width = 6,textInput("db", "Output database name", NULL)),
           column(width = 6,textInput("studyname", "Study name", NULL))
  ),
  tags$br(),
  fluidRow(column(width = 6,textInput("svorg", "Organism", NULL)),
           column(width = 6,textInput("svdisease", "Disease/Developmental stage", NULL))
  ),
  tags$br(),
  fluidRow(column(width = 6,textInput("tissue", "Tissue used for the scRNA-seq", NULL)),
           column(width = 6,textInput("shotdescr", "One-line description of the study", NULL))
           
  ),
  tags$br(),
  
  h4('5.1 Optional study details'),
  tags$hr(style = "border-top: 1px solid #d9d9d9;"),
  shiny::HTML('<i> these details are listed next to the study in the visualizations.</i>'),
  fluidRow(column(width = 11,textInput("descr", "Brief description of the study/Abstract", NULL)),
           column(width=1)
  ),
  tags$br(),
  fluidRow(column(width = 6,textInput("pmid", "If the study is published, then the PMID #", NULL)),
           column(width = 6,textInput("geo", "If the study is published, link to data (e.g. GEO)", NULL))
  ),
  tags$br(),
  fluidRow(column(width = 6,
                  selectInput("status",label="Select status of the study",choices = c('Internal','Publicly Available'),multiple = F)),
           column(width = 6,
                  selectInput("rating",label="Select your rating for the study",choices = c('High','Medium','Poor'),multiple = F))
  ),
  tags$br(),
  tags$br(),
  
  tags$hr(style = "border-top: 1px solid #000000;"),
  fluidRow(  column(width = 4,  
                    actionButton("cmd", "show batch command",icon =icon("eye"),width = '300px',class = "btn-secondary"),
                    shiny::HTML('If you want to run multiple studies in a batch, use the above botton to retrieve a command line for a given study. You can run multiple such commands using batch submission approach.')
                    ),
             column(width = 1),  
             column(width = 7,
                    htmlOutput("cmd_out"))
  ),
  tags$hr(style = "border-top: 1px solid #000000;"),
  fluidRow(column(width = 4),
           column(width = 4,  actionButton("run", "Go!",icon =icon("play"),width = '300px',class = "btn-warning")),
           column(width = 4)
  ),
  tags$hr(style = "border-top: 1px solid #000000;"),
  
  tags$br(),
  h2('6.Run details'),
  tags$hr(style = "border-top: 1px solid #000000;"),
  htmlOutput("finalSel"),
  
  
  tags$br()
  
)

shinyApp(ui, server)