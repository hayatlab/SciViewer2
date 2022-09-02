source("vars.R",local = TRUE,encoding = "UTF-8")
source("plot_functions.R",local = TRUE,encoding = "UTF-8")
source("MMdbLib.R",local = TRUE,encoding = "UTF-8")

##----------------------------------------------
warningMessage <- paste0("<font color=\"#e6550d\">
                           <br><br><br>
                           Requested plot/table is unavailable. <br>
                           Reasons:<br> 
                           1) Gene not present in the underlying database <br>
                           2) Gene locus is unavaialbe in Ensembl GRCh37 transcript database. <br>
                           3) A technical glitch with plot rendering? <br>
                           Try refreshing the page or contact the Data science team/ developers for details.
                           </font><br><br><br>")

##----------------------------------------------
js <- "
    function(el, x, inputName){
      var id = el.getAttribute('id');
      var gd = document.getElementById(id);
      var xScrollPos = el.scrollLeft || document.documentElement.scrollLeft;
      var yScrollPos = el.scrollTop || document.documentElement.scrollTop;
      var d3 = Plotly.d3;
      Plotly.plot(id).then(attach);
        function attach() {
          var xaxis = gd._fullLayout.xaxis;
          var yaxis = gd._fullLayout.yaxis;
          var l = gd._fullLayout.margin.l;
          var t = gd._fullLayout.margin.t;
          var coordinates = [null, null]
          
          gd.addEventListener('mousemove', 
          function(evt) {
              var coordinates = [xaxis.p2c(evt.pointerX-l), 
                                 yaxis.p2c(evt.pointerY-t)];
              Shiny.setInputValue(inputName, coordinates);
          });

        };
  }
  "

##----------------------------------------------

server <- function(input, output, session) {
  
  cat("DhawalJain\n")
  
  gc()
  Prog <- shiny::Progress$new()
  on.exit(Prog$close())
  Prog$set(message = "initializing the app..", value = 0)
  observe({ cat("tab item: ",input$AppTab,"\n") })

  ROOTS=c(workdir='.',
          datadir='C:/',
          home='/home/',
          shinydata='/srv/shiny-server/data/')
  
  ###---- 0) meta
  VARS <- reactiveValues(
    connGenes = list(),
    connList = list(),
    dbb = NULL,
    roots = ROOTS,
    sc_studies = c(),
    availGenes = NULL
  )
  
  variables <- reactiveValues(
    gene = NULL,
    sc_study = NULL,
    connID_1 = NULL,
    sc_study_attribs = NULL,
    donor = NULL,
    celltype = NULL,
    catvars = NULL,
    catvars_valuelist = NULL,
    contvars = NULL,
    availcelltypes = NULL,
    sc_louvain = NULL,
    sel_catvar = NULL,
    sel_catvarval = NULL,
    expndf=NULL,
    avail_markerTests = NULL,
    marker_cellTypeDf = NULL,
    marker_geneDf = NULL,
    cellmark_upenrich = list(),
    cellmark_dnenrich = list(),
    expndf_markers = NULL,
    avail_deTests = NULL,
    de_cellTypeDf = NULL,
    de_geneDf = NULL,
    cellde_upenrich = list(),
    cellde_dnenrich = list(),
    expndf_de = NULL
  ) 
  compvariables <- reactiveValues(
    gene = NULL,
    geneExp = c(),
    geneMarker = c(),
    geneBioMarker = c(),
    geneExpSummarized = NULL
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
        cat('Directory does not exits!\n')
        output$direxists <- renderPrint({ cat('Directory does not exits!\n') })
      }
    }
  })
  observeEvent({VARS$roots},{
    if(isTruthy(VARS$roots)){
      shinyFileChoose(input,'infile', session=session,roots=VARS$roots, 
                      filetypes=c('', 'db'))
    }
  })
  
  ###---  1) genes database
  Prog$set(message = "gathering available studies..", value = 0.3)
  VARS$connGenes <- ({
    RSQLite::dbConnect(RSQLite::SQLite(),GENEFILE)
  })
  
  ###---- 2) set the sqlite databases with access 
  VARS$dbb <- ({
    tryCatch({
      read.delim(file = SCDBFILE,header = F,comment.char = "#",stringsAsFactors = F)$V1
    },error= function(e){
    })
  })
  observeEvent({input$infile},{
    VARS$dbb = NULL
    VARS$dbb = tryCatch({
      read.delim(file = SCDBFILE,header = F,comment.char = "#",stringsAsFactors = F)$V1
    },error= function(e){
    })
    cat(VARS$roots,"\n")
    inFile <- parseFilePaths(roots=VARS$roots, input$infile)
    if(isTruthy(inFile$datapath)){
      ext <- tools::file_ext(unname(inFile$datapath))
      validate(need(ext %in% c("db"), "Please select files with .db extension!"))
      VARS$dbb <- append(VARS$dbb,unname(inFile$datapath))
      VARS$dbb <- unique(VARS$dbb)
      rm(ext)
    }
  })
  
  ##--- 3) check if remote databases are allowed
  observe({
    if(USE_REMOTE_DB==TRUE){
      VARS$dbb = NULL
      cf <- GetStudyName(AppName=REPO_NAME)
      cf <- as.character(unlist(cf[,1]))
      cf <- cf[grep("^scRNA_",cf)]
      if(isTruthy(cf)){
        VARS$dbb <- unique(cf)
      }
      Prog$set(message = "Using remote databases..", value = 0.4)
      rm(cf)
    }else{
      Prog$set(message = "Using local databases..", value = 0.4)
      print("Currently working with local databases\n")
    }
  })
  
  
  ##--- 4) populate meta object
  observeEvent({VARS$dbb},{
    if(isTruthy(VARS$dbb)){
      VARS$connList <- ({
        connList <- list()
        if(length(VARS$dbb)>0){
          for(f in VARS$dbb){
            id = pddExpn::randomStringGenerator(n = 1,lenght = 8)
            if(USE_REMOTE_DB==TRUE){
              connList[[id]] <- f #gsub("^scRNA_","",f)
            }else{
              connList[[id]] <- RSQLite::dbConnect(RSQLite::SQLite(),f)
            }
            rm(id,f)
          }
        }
        connList
      })
      VARS$sc_studies <- ({
        sc_studies <- c()
        for(i in 1:length(connList)){
          if(USE_REMOTE_DB==F){
            tablist <- RSQLite::dbListTables(connList[[i]])
            tablist <- tablist[grep("_study",tablist)]
          }else{
            tablist <- paste0(connList[[i]],"_study")
            #tablist <- gsub("^scRNA_","",tablist)
          }
          for(j in tablist){
            cat(i,"\t",j,"\n")
            query <- paste0("SELECT * FROM ",j)
            sc_studies <- rbind(sc_studies,
                                cbind(
                                  queryDB(HANDLER=connList[[i]], QUERY=query,
                                          REPO_NAME=REPO_NAME,
                                          USE_REMOTE_DB=USE_REMOTE_DB),
                                  ObjID=names(connList)[i]
                                )
            )
          }
          rm(i,j)
        }
        sc_studies
      })
    }
  })
  
  
  ##--- 5) list available studies
  output$browse_studytable <- DT::renderDataTable({
    if(isTruthy(nrow(VARS$sc_studies)>0)){
      cf <- VARS$sc_studies
      cf <- cf[,c("ORGANISM","DISEASE","TISSUES","STATUS","SampleSize","SHORTDESCR")]
      names(cf) <- c('Species','Disease/Condition','Tissue',"Status",'Number of Donors','Brief')
      datatable(cf,rownames = F,selection='single',extensions = 'Buttons',
                options = list(pageLength = 20, autoWidth = TRUE))
    }
  })
  

  
  source(file = "server_df.R",local = TRUE,encoding = "UTF-8")$value
  Prog$set(message = "retrieved data for vizualisations..", value = 0.7)
  
  source(file = "server_sc.R",local = TRUE,encoding = "UTF-8")$value
  Prog$set(message = "prepared data for vizualisations..", value = 0.8)
  
  source(file = "server_compare.R",local = TRUE,encoding = "UTF-8")$value
  Prog$set(message = "retrieved data for cross study comparison..", value = 0.9)
  
  Prog$set(message = "done!", value = 1)
  
}