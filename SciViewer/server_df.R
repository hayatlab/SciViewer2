## server_df.R

##-----------------------------------------------
##--- 1) study details, available genes etc
##-----------------------------------------------
observeEvent({input$browse_studytable_rows_selected},{
  if(isTruthy(input$browse_studytable_rows_selected)){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = "resetting/initializing dataframe extractions..", value = 0)
    cat(" Initializing the data\n")
    
    variables$gene = NULL
    variables$sc_study = NULL
    variables$connID_1 = NULL
    variables$sc_study_attribs = NULL
    variables$donor = NULL
    variables$celltype = NULL
    variables$catvars = NULL
    variables$catvars_valuelist = NULL
    variables$contvars = NULL
    variables$availcelltypes = NULL
    variables$sc_louvain = NULL
    variables$sel_catvar = NULL
    variables$sel_catvarval = NULL
    variables$expndf = NULL
    variables$avail_markerTests = NULL
    variables$marker_cellTypeDf = NULL
    variables$marker_geneDf = NULL
    variables$cellmark_upenrich = list()
    variables$cellmark_dnenrich = list()
    variables$expndf_markers = NULL
    variables$avail_deTests = NULL
    variables$de_cellTypeDf = NULL
    variables$de_geneDf = NULL
    variables$cellde_upenrich = list()
    variables$cellde_dnenrich = list()
    variables$expndf_de = NULL

    compvariables$gene = input$sel_genex
    compvariables$geneExp = c()
    compvariables$geneMarker = c()
    compvariables$geneBioMarker = c()
    compvariables$geneExpSummarized = NULL
    
    #-- get the study data frame
    cf <- VARS$sc_studies[input$browse_studytable_rows_selected,]
    
    #-- description
    variables$sc_study_attribs <- ({
      paste("<font color=\"#FF0000\"><b> STUDY STATUS: </b> ",as.character(cf$STATUS),"! </font><br>",
            "<b> Species: </b>", as.character(cf$ORGANISM),"<br>",
            "<b> Number of Cells: </b>", format(as.numeric(cf$CellCount),big.mark=",")," ",
            "<b> Number of features: </b>", format(as.numeric(cf$FeatureCount),big.mark=","),"<br>",
            "<b> TISSUE: </b>", as.character(cf$TISSUES)," ",
            "<b> DISEASE/CONDITION: </b>", as.character(cf$DISEASE)," ",
            "<b> SAMPLE SIZE: </b>", as.character(cf$SampleSize),"<br>",
            "<b> PUBMED: </b>", a(as.character(cf$PMID),href=paste0("https://pubmed.ncbi.nlm.nih.gov/",as.character(cf$PMID)), target="_blank"),
            "<b> ACCSSION: </b>", a("Data",href=paste0(as.character(cf$GEO)), target="_blank"),"<br>",
            "<b> STUDY ABSTRACT: </b> <font color=\"#bdbdbd\">", as.character(cf$Description),"</font><br>"
      )
    })
    
    #--- study
    variables$sc_study = as.character(cf$Database)
    
    #--- connID
    variables$connID_1 <- ({as.character(cf$ObjID) })
    
    #--- base dataframe
    variables$sc_louvain <- ({
      query <- paste0("SELECT * FROM ",variables$sc_study,"_metaFeatures")
      louvain <- queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                         QUERY=query,REPO_NAME=REPO_NAME,
                         USE_REMOTE_DB=USE_REMOTE_DB)
      louvain$value <- NA
      louvain
    })
    scProg$set(message = "retrieved required base dataframe..", value = 0.5)
    
    #--- celltype
    x <- unlist(strsplit(as.character(cf$cell_type),","))
    if(!is.null(x)){
      variables$celltype=x
    }else{
      scProg$set(message = "can't locate user-defined cell type column in the database. Exiting!", value = 0.5)
      stop("can't locate user-defined cell type column in the database. Exiting!")
    }
    rm(x)
    
    #--- donor_var
    x <- unlist(strsplit(as.character(cf$donor),","))
    if(!is.null(x)){
      variables$donor=x
    }
    
    #--- categorical variables
    x <- unlist(strsplit(as.character(cf$CATVAR),","))
    if(!is.null(variables$donor)){x <- x[x!=variables$donor]}
    x <- x[x!=variables$celltype]
    if(length(x)>1){
      variables$catvars=x
    }
    rm(x)
    
    #-- call available unique values for catvars
    if(!is.null(variables$catvars)){
      for(f in variables$catvars){
        if(f %in% names(variables$sc_louvain)){
          variables$catvars_valuelist[[f]] <- unique(variables$sc_louvain[,f])
        }else{
          variables$catvars <- variables$catvars[variables$catvars!=f]
        }
      }
    }

    #--- continuous variables
    x <- unlist(strsplit(as.character(cf$CONTVAR),","))
    if(!is.null(x)){
      variables$contvars=x
    }
    rm(x)
    
    #-- available cell types
    variables$availcelltypes <- ({
      unique(as.character(variables$sc_louvain[,as.character(variables$celltype)]))
    })
    
    #-- genes
    cat(" getting gene lists for: ", as.character(cf$Database), " objID: ", as.character(cf$ObjID),"\n")
    query <- paste0("SELECT geneSymbol FROM ",as.character(cf$Database),"_data")
    availGenes <- queryDB(HANDLER=VARS$connList[[as.character(cf$ObjID)]], 
                          QUERY=query,REPO_NAME=REPO_NAME,
                          USE_REMOTE_DB=USE_REMOTE_DB)
    updateSelectizeInput(session, 'sel_gene', choices = availGenes[,1], 
                         server = TRUE,selected ="ITGB8")
    updateSelectizeInput(session, 'sel_genex', choices = availGenes[,1], 
                         server = TRUE,selected ="ITGB8")

    #--- gene expression dataframe initialization
    variables$expndf = data.frame(SAMPID=variables$sc_louvain[,'SAMPID'])
    rownames(variables$expndf) = variables$expndf[,1]
    variables$expndf[,1] = NULL
    variables$expndf_markers = variables$expndf
    variables$expndf_de = variables$expndf
    
    scProg$set(message = "done..", value = 1)
  }
})

##-----------------------------------------------
##--- 2) alias gene names
##-----------------------------------------------
gAlias <- reactive({
  query <- paste0("SELECT * FROM gAliasProtein")
  gGenes <- RSQLite::dbGetQuery(VARS$connGenes, query)
  gAlias <- list()
  for(f in input$sel_gene){
    hgnc <- gGenes[gGenes$geneSymbol==f,]$HGNC %>% as.character
    gg <- as.character(gGenes[gGenes$HGNC%in%hgnc,]$geneSymbol)
    gAlias[[f]] <- unique(c(f,gg))
  }
  gAlias
})
observeEvent({variables$sc_study_attribs},{
  if(isTruthy(variables$sc_study_attribs)){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = "retrieving alias for gene names, if any..", value = 0)
    #-- messages
    scProg$set(message = paste0(" study1: ", variables$sc_study), value = 0.1)
    cat(" study1: ", variables$sc_study,"\n")
    scProg$set(message = paste0(" study1 connID: ", variables$connID_1), value = 0.2)
    cat(" study1 connID: ", variables$connID_1,"\n")
    scProg$set(message = paste0(" cell type column: ",variables$celltype), value = 0.3)
    cat(" cell type column: ",variables$celltype,"\n")
    scProg$set(message = paste0(" Donor column: ",variables$donor), value = 0.4)
    cat(" Donor column: ",variables$donor,"\n")
    scProg$set(message = paste0(" Categorical variables: ",variables$catvars), value = 0.5)
    cat(" Categorical variables: ",variables$catvars,"\n")
    scProg$set(message = paste0(" Continuous variables: ",variables$contvars), value = 0.6)
    cat(" Continuous variables: ",variables$contvars,"\n")
    scProg$set(message = paste0(" Available celltypes:  ",variables$availcelltypes), value = 0.7)
    cat(" Available celltypes:  ",variables$availcelltypes,"\n")
    scProg$set(message = paste0(" Alias gene names: ", paste0(unlist(gAlias()),collapse = ";")), value = 0.8)
    cat(" Alias gene names: ", paste0(unlist(gAlias()),collapse = ";"),"\n")
    cat("\n")
    scProg$set(message = "done..", value = 1)
  }
})


##-----------------------------------------------
##--- 3) update selection choices
##-----------------------------------------------



##-----------------------------------------------
##--- 4) gene expression dataframes
##-----------------------------------------------
observeEvent({input$sel_genego},{
    if(isTruthy(variables$sc_louvain) & isTruthy(input$sel_genego) &  isTruthy(input$sel_gene)){
      if(is.null(variables$gene) | sum(variables$gene!=input$sel_gene)>0 ){
        scProg <- shiny::Progress$new()
        on.exit(scProg$close())
        scProg$set(message = "retrieving plot data for multiple genes..", value = 0)
        cat(" Generating dataframe list for multiple genes..\n")
        
        variables$gene=input$sel_gene  
        
        ##------------ generate and update expndf
        for(g in input$sel_gene){
          y <- g
          if(!y %in% names(variables$expndf)){ #names(variables$pldf)
            pl <- NULL
            cntr =1 
            while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
              cat("  gene/alias name: ", gAlias()[[g]][cntr],"\n")
              pl <-  get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                                      variables$sc_study,
                                      gAlias()[[g]][cntr],
                                      variables$sc_louvain,
                                      REPO_NAME=REPO_NAME,
                                      USE_REMOTE_DB=USE_REMOTE_DB)
              cntr = cntr +1 
            }
            if(isTRUE(nrow(pl)>0)){
              scProg$set(message = paste0('  getting plot df for gene ',g, ":",y,":",names(variables$expndf)), value = 0.2)
              cat('  getting plot df for gene ',g, ":",y,":",names(variables$expndf),"\n")
              variables$expndf <- cbind(variables$expndf,pl$value)
              names(variables$expndf)[ncol(variables$expndf)] <- g
            }
          }else{
            scProg$set(message = paste0('plot dataframe for gene ',g, " is already calculated"), value = 0.2)
            cat(' plot dataframe for gene ',g, " is already calculated\n")
          }
        }
        rm(y,cntr,g)
        z <- paste0(input$sel_gene)
        for(y in names(variables$expndf)){ 
          if(!y%in%z){
            scProg$set(message = paste0('  removing plot dataframe for gene ',y), value = 0.5)
            cat('  removing plot dataframe for gene ',y, "\n")
            if(y %in%names(variables$expndf)){
              variables$expndf[,y] <- NULL
            }
          }
        }
        scProg$set(message = "done..", value = 1)
      }
    }
  })


##-----------------------------------------------
##--- 5) marker gene dataframes
##-----------------------------------------------
observeEvent(c(variables$gene),{
  if(isTruthy(variables$gene) & isTruthy(input$sel_genego)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("querying selected gene as marker across available features"), value = 0)
    cat("querying selected gene as marker across available features\n")
    variables$marker_geneDf = NULL
    
    cf <- local({
      cf <- c()
      for(mygene in variables$gene){
        query <- paste0("SELECT * FROM ",variables$sc_study,"_Marker WHERE geneSymbol = '",mygene,"'")
        scProg$set(message = paste0("querying gene marker: ",mygene), value = 0.5)
        cat("querying gene marker: ",mygene,"\n")
        pl <- tryCatch({
          queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                  QUERY=query,REPO_NAME=REPO_NAME,
                  USE_REMOTE_DB=USE_REMOTE_DB)
        },error=function(e){
          c()
        })
        cf <- rbind(cf,pl)
        rm(pl)
      }
      cf
    })
    if(isTRUE(nrow(cf)>0)){
      scProg$set(message = paste0("marker dataframe for ", variables$gene," has ",unique(cf$Test)," unique tests"), value = 0.7)
      cat("marker dataframe for ", variables$gene," has ",unique(cf$Test)," unique tests","\n")
      cf$Test <- paste0(cf$Test,"\n(",cf$geneSymbol,")")
      variables$marker_geneDf <- cf
    }
    rm(cf)
    gc()
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$sel_markercelltype),{
  if(isTruthy(input$sel_markercelltype) & isTruthy(variables$sc_study)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("retrieving marker genes for ",input$sel_markercelltype), value = 0)
    cat("retrieving marker genes for ",input$sel_markercelltype,"\n")
    #--- resetting the dataframes
    variables$avail_markerTests = NULL
    variables$marker_cellTypeDf = NULL
    variables$expndf_markers = NULL
    variables$expndf_markers = data.frame(SAMPID=variables$sc_louvain[,'SAMPID'])
    rownames(variables$expndf_markers) = variables$expndf_markers[,1]
    variables$expndf_markers[,1] = NULL

    query <- paste0("SELECT * FROM ",variables$sc_study,"_Marker WHERE Tag = '",input$sel_markercelltype,"'")
    pl <- tryCatch({
      queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
              QUERY=query,REPO_NAME=REPO_NAME,
              USE_REMOTE_DB=USE_REMOTE_DB)
    },error=function(e){
      return(NULL)
    })
    if(isTRUE(nrow(pl)>0)){
      cat("  marker dataframe for ", input$sel_markercelltype,": ",nrow(pl),"genes; test: ",unique(pl$Test),"\n")
      variables$avail_markerTests = unique(pl$Test)
      variables$marker_cellTypeDf <- pl
    }
    rm(pl)
    gc()
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$marker_tab_rows_selected),{
  if(isTruthy(input$marker_tab_rows_selected) & isTruthy(variables$marker_cellTypeDf) &
     isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("extracting expression values for selected marker genes from the table"), value = 0)
    cat("extracting expression values for selected marker genes from the table\n")
    
    ID = variables$marker_cellTypeDf[input$marker_tab_rows_selected,]$geneSymbol %>%
      as.character()
    scProg$set(message = paste0("marker gene/s ",ID), value = 0.2)
    cat("marker gene/s ",ID,"...\n")
    ##------------ generate and update expression data frame
    for(g in ID){
      if(!g %in% names(variables$expndf_markers)){
        pl <- NULL
        pl <- get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                               variables$sc_study,
                               g,
                               variables$sc_louvain,
                               REPO_NAME=REPO_NAME,
                               USE_REMOTE_DB=USE_REMOTE_DB)
        if(isTRUE(nrow(pl)>0)){
          cat('  getting plot df for gene ',g,"\n")
          variables$expndf_markers <- cbind(variables$expndf_markers,pl$value)
          names(variables$expndf_markers)[ncol(variables$expndf_markers)] <- g
        }
      }else{
        scProg$set(message = paste0("expression entry for marker gene ",g, " is already extracted"), value = 0.5)
        cat('expression entry for marker gene ',g, " is already extracted\n")
      }
    }
    rm(y,g)
    z <- ID
    for(y in names(variables$expndf_markers)){ 
      if(!y%in%z){
        scProg$set(message = paste0("removing expression entry for selected marker gene ",y), value = 0.5)
        cat('removing expression entry for selected marker gene ',y, "\n")
        if(y %in%names(variables$expndf_markers)){
          variables$expndf_markers[,y] <- NULL
        }
      }
    }
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$sc_markeropgo),{
  if(isTruthy(input$sc_markeropgo) & isTruthy(variables$sc_study) & 
     isTruthy(input$sel_markercelltype) &  isTRUE(nrow(variables$marker_cellTypeDf)>0) & 
     isTruthy(input$sel_markerfdrslider) & isTruthy(input$sel_markercomp)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("computing enrichment of the marker genes"), value = 0)
    cat("computing enrichment of the marker genes\n")
    
    variables$cellmark_upenrich = list()
    variables$cellmark_dnenrich = list()
    pl <- variables$marker_cellTypeDf 
    pl <- pl[pl$Test==input$sel_markercomp,]
    ## Enrichment
    if(isTRUE(nrow(pl)>0)){
      variables$cellmark_upenrich <- local({
        qf <-  pl
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$sel_markerfdrslider) & qf$logFC>0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Upregulated (wrt rest cells)'
          XX = gprofiler2::gost(query = qf[1],organism = 'hsapiens',#ORGANISM[[variables$sc_study[,"ORGANISM"]]],
                           ordered_query = F,exclude_iea = T,
                           significant = F,user_threshold = 0.05,correction_method = 'g_SCS')
          XX[['info']]  <-     shiny::HTML("<u><b>Celltype: </b>",input$sel_markercelltype,"       
                 <b>Comparison: </b>",input$sel_markercomp,"   
                 <b>FDR: </b>",input$sel_markerfdrslider,"            
                 <b>Direction of change: </b>Over-expressed </u><br>")
          XX
        }else{ NULL }
      })
      scProg$set(message = paste0("calculated enrichment for upregulated genes"), value = 0.5)
      cat("calculated enrichment for upregulated genes\n")
      variables$cellmark_dnenrich <- local({
        qf <-  pl
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$sel_markerfdrslider) & qf$logFC < 0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Downregulated (wrt rest cells)'
          XX <- gprofiler2::gost(query = qf[1],organism = 'hsapiens',#ORGANISM[[variables$sc_study[,"ORGANISM"]]],
                           ordered_query = F,exclude_iea = T,
                           significant = F,user_threshold = 0.05,correction_method = 'g_SCS')
          XX[['info']]  <-     shiny::HTML("<u><b>Celltype: </b>",input$sel_markercelltype,"       
                 <b>Comparison: </b>",input$sel_markercomp,"  
                 <b>FDR: </b>",input$sel_markerfdrslider,"            
                 <b>Direction of change: </b>Under-expressed </u><br>")
          XX
        }else{ NULL }
      })
      scProg$set(message = paste0("calculated enrichment for downregulated genes"), value = 0.9)
      cat("calculated enrichment for downregulated genes\n")
      gc()
      scProg$set(message = "done..", value = 1)
    }
  }
})



##-----------------------------------------------
##--- 6) biological marker gene dataframes
##-----------------------------------------------
observeEvent(c(variables$gene),{
  if(isTruthy(variables$gene) & isTruthy(input$sel_genego)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("querying selected gene as biological marker across available comparisons"), value = 0)
    cat("querying selected gene as biological marker across available comparisons\n")
    variables$de_geneDf = NULL
    
    cf <- local({
      cf <- c()
      for(mygene in variables$gene){
        query <- paste0("SELECT * FROM ",variables$sc_study,"_DEG WHERE geneSymbol = '",mygene,"'")
        scProg$set(message = paste0("querying gene marker: ",mygene), value = 0.5)
        cat("querying gene bioloical marker: ",mygene,"\n")
        pl <- tryCatch({
          queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                  QUERY=query,REPO_NAME=REPO_NAME,
                  USE_REMOTE_DB=USE_REMOTE_DB)
        },error=function(e){
          c()
        })
        cf <- rbind(cf,pl)
        rm(pl)
      }
      cf
    })
    if(isTRUE(nrow(cf)>0)){
      scProg$set(message = paste0("biological marker dataframe for ", variables$gene," has ",unique(cf$Test)," unique tests"), value = 0.7)
      cat("biological marker dataframe for ", variables$gene," has ",unique(cf$Test)," unique tests","\n")
      cf$Test <- paste0(cf$Test,"\n(",cf$geneSymbol,")")
      variables$de_geneDf <- cf
    }
    rm(cf)
    gc()
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$sel_decelltype),{
  if(isTruthy(input$sel_decelltype) & isTruthy(variables$sc_study)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("retrieving biological marker genes for ",input$sel_decelltype), value = 0)
    cat("retrieving biological marker genes for ",input$sel_decelltype,"\n")
    #--- resetting the dataframes
    variables$avail_deTests = NULL
    variables$de_cellTypeDf = NULL
    variables$expndf_de = NULL
    variables$expndf_de = data.frame(SAMPID=variables$sc_louvain[,'SAMPID'])
    rownames(variables$expndf_de) = variables$expndf_de[,1]
    variables$expndf_de[,1] = NULL
    
    query <- paste0("SELECT * FROM ",variables$sc_study,"_DEG WHERE Tag = '",input$sel_decelltype,"'")
    pl <- tryCatch({
      queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
              QUERY=query,REPO_NAME=REPO_NAME,
              USE_REMOTE_DB=USE_REMOTE_DB)
    },error=function(e){
      return(NULL)
    })
    if(isTRUE(nrow(pl)>0)){
      cat("Biological marker dataframe for ", input$sel_decelltype,": ",nrow(pl),"genes; test: ",unique(pl$Test),"\n")
      variables$avail_deTests = unique(pl$Test)
      variables$de_cellTypeDf <- pl
    }
    rm(pl)
    gc()
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$de_tab_rows_selected),{
  if(isTruthy(input$de_tab_rows_selected) & isTruthy(variables$de_cellTypeDf) &
     isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("extracting expression values for selected biological marker genes from the table"), value = 0)
    cat("extracting expression values for selected biological marker genes from the table\n")
    
    ID = variables$de_cellTypeDf[input$de_tab_rows_selected,]$geneSymbol %>%
      as.character()
    scProg$set(message = paste0("marker gene/s ",ID), value = 0.2)
    cat("marker gene/s ",ID,"...\n")
    ##------------ generate and update expression data frame
    for(g in ID){
      if(!g %in% names(variables$expndf_de)){
        pl <- NULL
        pl <- get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                               variables$sc_study,
                               g,
                               variables$sc_louvain,
                               REPO_NAME=REPO_NAME,
                               USE_REMOTE_DB=USE_REMOTE_DB)
        if(isTRUE(nrow(pl)>0)){
          cat('  getting plot df for gene ',g,"\n")
          variables$expndf_de <- cbind(variables$expndf_de,pl$value)
          names(variables$expndf_de)[ncol(variables$expndf_de)] <- g
        }
      }else{
        scProg$set(message = paste0("expression entry for biological marker gene ",g, " is already extracted"), value = 0.5)
        cat('expression entry for biological marker gene ',g, " is already extracted\n")
      }
    }
    rm(y,g)
    z <- ID
    for(y in names(variables$expndf_de)){ 
      if(!y%in%z){
        scProg$set(message = paste0("removing expression entry for selected biological marker gene ",y), value = 0.5)
        cat('removing expression entry for selected biological marker gene ',y, "\n")
        if(y %in%names(variables$expndf_de)){
          variables$expndf_de[,y] <- NULL
        }
      }
    }
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(input$sc_deopgo),{
  if(isTruthy(input$sc_deopgo) & isTruthy(variables$sc_study) & 
     isTruthy(input$sel_decelltype) &  isTRUE(nrow(variables$de_cellTypeDf)>0) & 
     isTruthy(input$sel_defdrslider) & isTruthy(input$sel_decomp)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("computing enrichment of the biological marker genes"), value = 0)
    cat("computing enrichment of the biological marker genes\n")
    
    variables$cellde_upenrich = list()
    variables$cellde_dnenrich = list()
    pl <- variables$de_cellTypeDf 
    pl <- pl[pl$Test==input$sel_decomp,]
    ## Enrichment
    if(isTRUE(nrow(pl)>0)){
      variables$cellde_upenrich <- local({
        qf <-  pl
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$sel_defdrslider) & qf$logFC>0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Upregulated (wrt rest cells)'
          XX = gprofiler2::gost(query = qf[1],organism = 'hsapiens',#ORGANISM[[variables$sc_study[,"ORGANISM"]]],
                                ordered_query = F,exclude_iea = T,
                                significant = F,user_threshold = 0.05,correction_method = 'g_SCS')
          XX[['info']]  <-     shiny::HTML("<u><b>Celltype: </b>",input$sel_decelltype,"       
                 <b>Comparison: </b>",input$sel_decomp,"   
                 <b>FDR: </b>",input$sel_defdrslider,"            
                 <b>Direction of change: </b>Over-expressed </u><br>")
          XX
        }else{ NULL }
      })
      scProg$set(message = paste0("calculated enrichment for upregulated genes"), value = 0.5)
      cat("calculated enrichment for upregulated genes\n")
      variables$cellde_dnenrich <- local({
        qf <-  pl
        qf <- qf[as.numeric(qf$adj.P.Val) <= as.numeric(input$sel_defdrslider) & qf$logFC < 0,]$geneSymbol
        qf <- unique(qf)
        if(length(qf)>2){
          qf = list(A=qf)
          names(qf) <- 'Downregulated (wrt rest cells)'
          XX <- gprofiler2::gost(query = qf[1],organism = 'hsapiens',#ORGANISM[[variables$sc_study[,"ORGANISM"]]],
                                 ordered_query = F,exclude_iea = T,
                                 significant = F,user_threshold = 0.05,correction_method = 'g_SCS')
          XX[['info']]  <-     shiny::HTML("<u><b>Celltype: </b>",input$sel_decelltype,"       
                 <b>Comparison: </b>",input$sel_decomp,"  
                 <b>FDR: </b>",input$sel_defdrslider,"            
                 <b>Direction of change: </b>Under-expressed </u><br>")
          XX
        }else{ NULL }
      })
      scProg$set(message = paste0("calculated enrichment for downregulated genes"), value = 0.9)
      cat("calculated enrichment for downregulated genes\n")
      gc()
      scProg$set(message = "done..", value = 1)
    }
  }
})


##-----------------------------------------------
##--- 7) gene across studies
##-----------------------------------------------
observeEvent(c(input$sel_genexgo),{
  if(input$AppTab=='compare' & isTruthy(input$sel_genexgo) & 
     isTruthy(input$sel_genex) & isTruthy(VARS$sc_studies)){
    cat(" extracting data for the gene",input$sel_genex," across studies\n")
    
    compvariables$gene = input$sel_genex
    compvariables$geneExp = c()
    compvariables$geneMarker = c()
    compvariables$geneBioMarker = c()
    compvariables$geneExpSummarized = NULL
    
    gAliasX <- local({
      query <- paste0("SELECT * FROM gAliasProtein WHERE geneSymbol = '",compvariables$gene,"'")
      gGenes <- RSQLite::dbGetQuery(VARS$connGenes, query)
      hgnc <- gGenes$HGNC %>% as.character()
      gg <- as.character(gGenes[gGenes$HGNC%in%hgnc,]$geneSymbol)
      gAliasX <- list()
      gAliasX[[compvariables$gene]] <- unique(c(compvariables$gene,gg))
      gAliasX
    })
    
    for(i in 1:nrow(VARS$sc_studies)){
      pl <- NULL
      ql <- NULL
      gl <- NULL
      cntr =1 
      cf <- VARS$sc_studies[i,]
      cat(" extracting data for the gene",compvariables$gene," in ",cf$Database,"\n")
      while(cntr<=length(gAliasX[[compvariables$gene]]) && !isTRUE(nrow(pl)>0)){
        cat("  gene/alias name: ", gAliasX[[compvariables$gene]][cntr],"\n")
        ##-- 01) Expression
        pl <- local({
          query <- paste0("SELECT * FROM ",as.character(cf$Database),"_FeatureSummary WHERE feature = '",compvariables$gene,"'")
          tryCatch({
            queryDB(HANDLER=VARS$connList[[as.character(cf$ObjID)]], 
                    QUERY=query,REPO_NAME=REPO_NAME,
                    USE_REMOTE_DB=USE_REMOTE_DB)
          },error=function(e){ NULL
          })
        })
        ##-- 02) Marker
        gl <- local({
          query <- paste0("SELECT * FROM ",as.character(cf$Database),"_Marker WHERE geneSymbol = '",compvariables$gene,"'")
          tryCatch({
            queryDB(HANDLER=VARS$connList[[as.character(cf$ObjID)]], 
                    QUERY=query,REPO_NAME=REPO_NAME,
                    USE_REMOTE_DB=USE_REMOTE_DB)
          },error=function(e){ NULL
          })
        })
        ##-- 03) BioMarker
        ql <- local({
          query <- paste0("SELECT * FROM ",as.character(cf$Database),"_DEG WHERE geneSymbol = '",compvariables$gene,"'")
          tryCatch({
            queryDB(HANDLER=VARS$connList[[as.character(cf$ObjID)]], 
                    QUERY=query,REPO_NAME=REPO_NAME,
                    USE_REMOTE_DB=USE_REMOTE_DB)
          },error=function(e){ NULL
          })
        })
        ## counter
        cntr = cntr +1 
      }
      if(isTRUE(nrow(pl)>0)){
        pl$Database <- as.character(cf$Database)
        compvariables$geneExp <- rbind(compvariables$geneExp, pl)
        cat("Expression Comparison: A study added\n")
      }
      if(isTRUE(nrow(gl)>0)){
        gl$Database <- as.character(cf$Database)
        compvariables$geneMarker <- rbind(compvariables$geneMarker, gl)
        cat("Marker Expression Comparison: A study added\n")
      }
      if(isTRUE(nrow(ql)>0)){
        ql$Database <- as.character(cf$Database)
        compvariables$geneBioMarker <- rbind(compvariables$geneBioMarker, ql)
        cat("Biological Marker Expression Comparison: A study added\n")
      }
      rm(pl,ql,gl)
      gc()
    }
    
    ## summarized table
    myfun01 <- function(x){ log2(mean(2^x,na.rm=T)) }
    compvariables$geneExpSummarized <- ({
      compvariables$geneExp %>%
        dplyr::group_by(cell_type,feature,Database) %>%
        dplyr::summarize(AveExpr=myfun01(norm_avg_priorLT))
    })

    
  }
})






###----------------------------------------------------------
###--------- trash
###----------------------------------------------------------

if(F){
  
  ###----------------------------------------------------------
  ###---------------------- Study 1
  ###----------------------------------------------------------
  ## sc_louvain, assign variables based on study name
  observeEvent({input$sc_study},{
    if(isTruthy(input$sc_study)  &
       (is.null(variables$sc_study) | isTruthy(input$sc_study != variables$sc_study))
    ){
      
      sc_studies = VARS$sc_studies 
      variables$sc_study = input$sc_study
      
      
      cat("###----------------------------------\n")
      cat(" Inside bracket-1\n setting the study values to null\n\n")
      
      variables$connID_1 = NULL
      variables$feature = NULL
      variables$add_features = NULL
      variables$cellTypes = NULL
      variables$add_feature_list = list()
      variables$sc_attribute = NULL
      variables$sc_louvain = NULL
      variables$sc_df = NULL
      variables$pldf = list()
      variables$pldf_marker = list()
      variables$pldf_deg = list()
      variables$cordf = NULL
      variables$corrdf_sel=NULL
      variables$marker_dot = c()
      variables$marker_volc = NULL
      variables$gene = NULL
      variables$geneA = NULL
      variables$markenrich_up = list()
      variables$markenrich_dn = list()
      variables$degenrich_up = list()
      variables$degenrich_dn = list()
      variables$markerTests = NULL
      variables$degTests = NULL
      gc()
      
      scProg$set(message = "Calculated louvain DF..", value = 0.5)
      
      scProg$set(message = "done..", value = 1)
    }
  })
  
  ## generate data frame for genes
  observeEvent({input$geneA
    variables$sc_louvain},{
      
      if(isTruthy(input$sc_study) & isTruthy(input$geneA) & 
         isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain)){
        
        if(is.null(variables$gene) | isTruthy(variables$gene!=input$geneA[1])
        ){
          variables$sc_df = NULL
          variables$gene=input$geneA[1]
          variables$marker_dot = c()
          variables$deg_dot = c()
          
          cat(" Inside bracket-2\nvariables$gene=", variables$gene,
              "\nconnID: ", variables$connID_1,"\n")
          
          z <- paste0(input$sc_study,'_',variables$gene)
          ##-- use all gene Aliases to search the database 
          cntr = 1
          cat(" xxx: ",length(gAlias()[[variables$gene]])," ", isTRUE(nrow(variables$sc_df))," ",gAlias()[[variables$gene]][cntr],"\n")
          
          while(cntr <= length(gAlias()[[variables$gene]]) &&
                !isTRUE(nrow(variables$sc_df)>0)){
            cat("  gene/alias name: ", gAlias()[[variables$gene]][cntr],"\n")
            variables$sc_df <- variables$pldf[[z]] <- 
              get_plot_df_sc04(connSc = VARS$connList[[variables$connID_1]],
                               study = input$sc_study,
                               genename = gAlias()[[variables$gene]][cntr],
                               louvain = variables$sc_louvain,
                               REPO_NAME=REPO_NAME,
                               USE_REMOTE_DB=USE_REMOTE_DB)
            cntr = cntr+1
          }
          #if(is.null(variables$sc_df)){
          #  variables$sc_df <- NULL
          #  variables$pldf[[z]] <- NULL
          #}
          cat(" names(variables$pldf):",names(variables$pldf),"\n")
          cat("calculated first-gene DF for plot\n\n")
          #scProg$set(message = "done..", value = 1)
        }
        
        if(is.null(variables$geneA) | 
           sum(variables$geneA!=input$geneA)>0
        ){
          
          variables$geneA=input$geneA  
          variables$marker_dot = c()
          variables$deg_dot = c()
          
          scProg <- shiny::Progress$new()
          on.exit(scProg$close())
          scProg$set(message = "calculating DF for multiple genes..", value = 0)
          
          cat(" Inside bracket-3\ncalculating for multiple genes..\n")
          ##------------ generate and update pldf
          for(g in input$geneA){
            y <- paste0(input$sc_study,'_',g)
            if(!y %in% names(variables$pldf)){
              pl <- NULL
              cntr =1 
              while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
                cat("  gene/alias name: ", gAlias()[[g]][cntr],"\n")
                pl <-  get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                                        input$sc_study,
                                        gAlias()[[g]][cntr],
                                        variables$sc_louvain,
                                        REPO_NAME=REPO_NAME,
                                        USE_REMOTE_DB=USE_REMOTE_DB)
                cntr = cntr +1 
              }
              if(isTRUE(nrow(pl)>0)){
                cat('  getting plot df for gene ',g, ":",y,":",names(variables$pldf),"\n")
                variables$pldf[[y]] <- pl
              }
            }else{
              cat('  plot df for gene ',g, " is already calculated\n")
            }
          }
          rm(y,cntr,g)
          z <- paste0(input$sc_study,'_',input$geneA)
          for(y in names(variables$pldf)){
            if(!y%in%z){
              cat('  removing plot df for gene ',y, "\n")
              if(y %in%names(variables$pldf)){
                variables$pldf[[y]] <- NULL
              }
            }
          }
          
          ###------- multigene corr 
          variables$corrdf <- NULL
          variables$corrdf_sel <- NULL
          if(length(variables$pldf)>1){
            scProg$set(message = "Calculating multigene correlation matrix..", value = 0.8)
            cat('  Calculating multigene correlation matrix. Available cell group: ', variables$feature,"\n")
            cat(" ", names(variables$pldf),"\n")
            qw <- data.frame(cell_type=variables$sc_df[,'cell_type'])
            if(nrow(qw)>0){
              for(i in 1:length(variables$pldf)){
                if(nrow(variables$pldf[[i]])>0){
                  cat(" corr matrix: i=",i," gene=",names(variables$pldf)[i],"\n")
                  qw <- cbind(qw,variables$pldf[[i]][,'value'])
                }
              }
              gnames <- gsub(paste0(input$sc_study,"_"),"",names(variables$pldf))
              names(qw)[2:ncol(qw)] <- gnames
              variables$corrdf <- qw
              variables$corrdf_sel <- qw[,2:ncol(qw)]
            }
          }
          cat("\n")
          
          ###---------- get marker results
          scProg$set(message = "Retrieving Markers..", value = 0.85)
          for(g in input$geneA){
            pl <- NULL
            cntr <-1
            while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
              cat("  Marker gene/alias name: ", gAlias()[[g]][cntr],"\n")
              query <- paste0("SELECT * FROM ",input$sc_study,"_Marker WHERE geneSymbol = '",gAlias()[[g]][cntr],"'")
              pl <- tryCatch({
                queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                        QUERY=query,REPO_NAME=REPO_NAME,
                        USE_REMOTE_DB=USE_REMOTE_DB)
                
              },error= function(e){
                return(NULL)
              })
              cntr = cntr +1 
            }
            if(isTRUE(nrow(pl)>0)){
              cat( " Marker rows for dotplot: ",nrow(pl),"\n")
              variables$marker_dot <- rbind(variables$marker_dot,pl)
            }
          }
          rm(cntr,g,pl)
          
          ###---------- get deg results
          scProg$set(message = "Retrieving DEG", value = 0.9)
          for(g in input$geneA){
            pl <- NULL
            cntr <-1
            while(cntr<=length(gAlias()[[g]]) && !isTRUE(nrow(pl)>0)){
              cat("  Marker gene/alias name: ", gAlias()[[g]][cntr],"\n")
              query <- paste0("SELECT * FROM ",input$sc_study,"_DEG WHERE geneSymbol = '",gAlias()[[g]][cntr],"'")
              pl <- tryCatch({
                queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                        QUERY=query,REPO_NAME=REPO_NAME,
                        USE_REMOTE_DB=USE_REMOTE_DB)
                
              },error= function(e){
                return(NULL)
              })
              cntr = cntr +1 
            }
            if(isTRUE(nrow(pl)>0)){
              cat( " DEG rows for dotplot: ",nrow(pl),"\n")
              variables$deg_dot <- rbind(variables$deg_dot,pl)
            }
          }
          rm(cntr,g,pl)
          
          ###---------- end
          scProg$set(message = "done..", value = 1)
          
        }
      }
    })
  
  ## select cell types for correlation matrices
  observeEvent({input$sc_sel_3a},{
    if(isTruthy(input$sc_sel_3a) & isTruthy(variables$corrdf)){
      cat("Inside bracket-4\n\n")
      qw <- variables$corrdf
      qw <- qw[qw$cell_type%in%input$sc_sel_3a,]
      cat("sc_sel_3a ", input$sc_sel_3a,"  nrow:",nrow(qw),"\n")
      cat(paste0(qw[1,],collapse = ";"),"\n")
      variables$corrdf_sel <- qw[,2:ncol(qw)]
    }else if(isTruthy(variables$corrdf)){
      qw <- variables$corrdf
      cat("NULLVal: sc_sel_3a ", input$sc_sel_3a,"  nrow:",nrow(qw),"\n")
      variables$corrdf_sel <- qw[,2:ncol(qw)]
    }
  },ignoreNULL = F)
  
  
  ###----------------------------------------------------------
  ###---------------------- Enrichment analyses for Markers/DEG
  ###----------------------------------------------------------
  ## select cell types for displaying Marker volcano plots
  observeEvent({input$sc_celltype4a},{
    if(isTruthy(input$sc_celltype4a)){
      variables$pldf_marker = list()
      variables$markerTests = NULL
      query <- paste0("SELECT * FROM ",input$sc_study,"_Marker WHERE Tag = '",input$sc_celltype4a,"'")
      pl <- tryCatch({
        queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                QUERY=query,REPO_NAME=REPO_NAME,
                USE_REMOTE_DB=USE_REMOTE_DB)
      },error=function(e){
        return(NULL)
      })
      if(isTRUE(nrow(pl)>0)){
        cat("  getting the Marker data frame for volcano plot\t", input$sc_celltype4a,"\t",nrow(pl),"\n")
        variables$markerTests = unique(pl$Test)
        variables$marker_volc <- pl
      }
    }
  })
  
  ## select cell types for displaying DEG volcano plots
  observeEvent({input$sc_celltype5a},{
    if(isTruthy(input$sc_celltype5a)){
      variables$pldf_deg = list()
      variables$degTests = NULL
      query <- paste0("SELECT * FROM ",input$sc_study,"_DEG WHERE Tag = '",input$sc_celltype5a,"'")
      pl <- tryCatch({
        queryDB(HANDLER=VARS$connList[[variables$connID_1]], 
                QUERY=query,REPO_NAME=REPO_NAME,
                USE_REMOTE_DB=USE_REMOTE_DB)
      },error=function(e){
        return(NULL)
      })
      if(isTRUE(nrow(pl)>0)){
        cat("  getting the DEG data frame for volcano plot\t", input$sc_celltype5a,"\t",nrow(pl),"\n")
        variables$degTests = unique(pl$Test)
        variables$deg_volc <- pl
      }
    }
  })
  
  
  ###----------------------------------------------------------
  ###---------------------- Option selection in the app
  ###----------------------------------------------------------
  ## Update input choices across the app
  observe({
    sc_studies = VARS$sc_studies
    cat(" Inside bracket-5\n")
    cat(" SelectInputbar population_1..\n")
    cat("\n")
    updateSelectizeInput(session, 'sc_celltype1a', 
                         choices = variables$cellTypes, 
                         server = TRUE,selected =NULL)
    updateSelectizeInput(session, 'sc_celltype2a', 
                         choices = variables$cellTypes, 
                         server = TRUE,selected =NULL)
    updateSelectizeInput(session, 'sc_celltype3a', 
                         choices = variables$cellTypes, 
                         server = TRUE,selected =NULL)
    updateSelectizeInput(session, 'sc_celltype4a', 
                         choices = variables$cellTypes, 
                         server = TRUE,selected =NULL)
    updateSelectizeInput(session, 'sc_celltype5a', 
                         choices = variables$cellTypes, 
                         server = TRUE,selected =NULL)
    updateSelectizeInput(session, 'sc_celltype6a', 
                         choices = variables$cellTypes, 
                         server = TRUE,selected =NULL)
    updateSelectizeInput(session, 'sc_sel_3a', 
                         choices = variables$cellTypes, 
                         server = TRUE,selected =NULL)
    ##
    feat <- unlist(strsplit(as.character(sc_studies[sc_studies$Database==input$sc_study,]$CATVAR),","))
    updateSelectizeInput(session, 'sc_sel_2a',  
                         choices = feat, 
                         server = TRUE,selected =feat[1])
    #updateSelectizeInput(session, 'sc_sel_1n',  
    #                     choices = feat, 
    #                     server = TRUE,selected =feat[1])
    
    ## additional pca display
    if(!is.null(variables$add_feature)){
      updateSelectizeInput(session, 'sc_sel_1',  
                           choices = variables$add_feature, 
                           server = TRUE,selected =variables$add_feature[1])
    }
  })
  
  ## select 1a
  observe({
    cat(" Inside bracket-6\n")
    cat(" SelectInputbar population_2..\n")
    cat("\n")
    updateSelectizeInput(session, 'sc_sel_1a', 
                         choices = variables$add_feature_list[[input$sc_sel_1]], 
                         server = TRUE,selected =NULL)
  })
  ### expression color scale choices
  observe({
    cat("Inside bracket-7\n\n")
    updateSelectizeInput(session, 'c1t2expnscale', 
                         choices = c("GrayBlue","Clay","magma","Viridis",
                                     "GreenYellow","Purple","Reds",
                                     "Oranges"), 
                         server = TRUE,selected =NULL)
  })
  observe({
    cat("Inside bracket-7a\n\n")
    updateSelectizeInput(session, 'sc_sel_4a', 
                         choices = variables$markerTests, 
                         server = TRUE,selected =NULL)
  })
  observe({
    cat("Inside bracket-7b\n\n")
    updateSelectizeInput(session, 'sc_sel_5a', 
                         choices = variables$degTests, 
                         server = TRUE,selected =NULL)
  })
  
  ###----------------------------------------------------------
  ###---------------------- Multigene display for celltype marker genes
  ###----------------------------------------------------------
  observeEvent({input$marker_tab_rows_selected
    variables$sc_louvain},{
      
      if(isTruthy(input$marker_tab_rows_selected) & 
         isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain)){
        
        ID = variables$marker_volc[input$marker_tab_rows_selected,]$geneSymbol %>%
          as.character()
        
        cat(" Inside bracket-6Marker\ncalculating for multiple genes..\n")
        ##------------ generate and update pldf
        for(g in ID){
          y <- paste0(input$sc_study,'_',g)
          if(!y %in% names(variables$pldf_marker)){
            pl <- NULL
            pl <- get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                                   input$sc_study,
                                   g,
                                   variables$sc_louvain,
                                   REPO_NAME=REPO_NAME,
                                   USE_REMOTE_DB=USE_REMOTE_DB)
            if(isTRUE(nrow(pl)>0)){
              cat('  getting plot df for gene ',g, ":",y,":",names(variables$pldf_marker),"\n")
              variables$pldf_marker[[y]] <- pl
            }
          }else{
            cat('  plot df for gene ',g, " is already calculated\n")
          }
        }
        rm(y,g)
        z <- paste0(input$sc_study,'_',ID)
        for(y in names(variables$pldf_marker)){
          if(!y%in%z){
            cat('  removing plot df for gene ',y, "\n")
            if(y %in%names(variables$pldf_marker)){
              variables$pldf_marker[[y]] <- NULL
            }
          }
        }
        
      }
    })
  
  ###----------------------------------------------------------
  ###---------------------- Multigene display for celltype DEG genes
  ###----------------------------------------------------------
  observeEvent({input$deg_tab_rows_selected
    variables$sc_louvain},{
      
      if(isTruthy(input$deg_tab_rows_selected) & 
         isTruthy(variables$sc_study) & isTruthy(variables$sc_louvain)){
        
        ID = variables$deg_volc[input$deg_tab_rows_selected,]$geneSymbol %>%
          as.character()
        
        cat(" Inside bracket-6DEG\ncalculating for multiple genes..\n")
        ##------------ generate and update pldf
        for(g in ID){
          y <- paste0(input$sc_study,'_',g)
          if(!y %in% names(variables$pldf_deg)){
            pl <- NULL
            pl <- get_plot_df_sc04(VARS$connList[[variables$connID_1]],
                                   input$sc_study,
                                   g,
                                   variables$sc_louvain,
                                   REPO_NAME=REPO_NAME,
                                   USE_REMOTE_DB=USE_REMOTE_DB)
            if(isTRUE(nrow(pl)>0)){
              cat('  getting plot df for gene ',g, ":",y,":",names(variables$pldf_deg),"\n")
              variables$pldf_deg[[y]] <- pl
            }
          }else{
            cat('  plot df for gene ',g, " is already calculated\n")
          }
        }
        rm(y,g)
        z <- paste0(input$sc_study,'_',ID)
        for(y in names(variables$pldf_deg)){
          if(!y%in%z){
            cat('  removing plot df for gene ',y, "\n")
            if(y %in%names(variables$pldf_deg)){
              variables$pldf_deg[[y]] <- NULL
            }
          }
        }
        
      }
    })
  
  ###----------------------------------------------------------
  ###---------------------- Study 2
  ###----------------------------------------------------------
  
  observeEvent({input$comp_study
    input$t2action1},{
      if(isTruthy(input$comp_study) & isTruthy(variables$gene) &
         input$AppTab=='compare'
      ){
        
        sc_studies = VARS$sc_studies
        cat("setting the comp-study values to null\n")
        variables$comp_study = input$comp_study
        variables$comp_feature = NULL
        variables$comp_add_features = NULL
        variables$comp_cellTypes = NULL
        variables$comp_add_feature_list = list()
        variables$comp_attribute = NULL
        variables$comp_louvain = NULL
        variables$comp_df = NULL
        
        scProg <- shiny::Progress$new()
        on.exit(scProg$close())
        scProg$set(message = "initializing DF generation for comparing study..", value = 0)
        
        cat("Inside bracket-8\ncalculating df for comparison study..\n")
        
        connID2 <- reactive({
          sc_studies[sc_studies$Database==input$comp_study,]$ObjID %>%
            as.character
        })
        
        variables$comp_louvain <- local({
          query <- paste0("SELECT * FROM ",input$comp_study,"_metaFeatures")
          louvain <- queryDB(HANDLER=VARS$connList[[connID2()]], 
                             QUERY=query,REPO_NAME=REPO_NAME,
                             USE_REMOTE_DB=USE_REMOTE_DB)
          louvain$value <- 0
          louvain
        })
        scProg$set(message = "Calculated louvain DF..", value = 0.5)
        
        variables$comp_attribute <- 
          paste("<font color=\"#FF0000\"><b> STUDY STATUS: </b> ",sc_studies[sc_studies$Database==input$comp_study,]$STATUS,"! </font><br>",
                "<b> # of cells: </b>",format(nrow(variables$comp_louvain) ,big.mark = ",",scientific = F),"<br>",
                "<b> TISSUE: </b>", sc_studies[sc_studies$Database==input$comp_study,]$TISSUES," ",
                "<b> SAMPLE SIZE: </b>", sc_studies[sc_studies$Database==input$comp_study,]$SampleSize,"<br>",
                "<b> PUBMED: </b>", a(sc_studies[sc_studies$Database==input$comp_study,]$PMID,href=paste0("https://pubmed.ncbi.nlm.nih.gov/",sc_studies[sc_studies$Database==input$comp_study,]$PMID), target="_blank"),
                "<b> ACCSSION: </b>", a("COVID19",href=paste0(sc_studies[sc_studies$Database==input$comp_study,]$GEO), target="_blank"),"<br>",
                "<b> STUDY ABSTRACT: </b> <font color=\"#bdbdbd\">", sc_studies[sc_studies$Database==input$comp_study,]$Description,"</font><br>"
          )
        
        ## feature column, this part should go! update the database column name  
        x <- unlist(strsplit(as.character(sc_studies[sc_studies$Database==input$comp_study,]$CATVAR),","))
        if(!is.null(x)){
          variables$comp_feature=x[1]
        }else{
          variables$comp_feature='cell_type'
        }
        
        #####-----------------------------------------------
        ##-- use all gene Aliases to search the database 
        cntr = 1
        while(cntr <= length(gAlias()[[variables$gene]]) &&
              !isTRUE(nrow(variables$comp_df)>0)){
          cat("  compStudy gene/alias name: ", gAlias()[[variables$gene]][cntr],"\n")
          variables$comp_df <- 
            get_plot_df_sc04(connSc = VARS$connList[[connID2()]],
                             study = input$comp_study,
                             genename = gAlias()[[variables$gene]][cntr],
                             louvain = variables$comp_louvain,
                             REPO_NAME=REPO_NAME,
                             USE_REMOTE_DB=USE_REMOTE_DB)
          cntr = cntr+1
        }
        rm(cnt)
        scProg$set(message = "Calculated plot DF..", value = 0.75)
        
        ### cell type in study 1
        updateSelectizeInput(session, 'sc_celltype21', 
                             choices = variables$cellTypes, 
                             server = TRUE,selected =NULL)
        ### cell type in study 2
        cat( 'study b feature: ', variables$comp_feature,'\n')
        b_celltypes <- unique(variables$comp_df[,variables$comp_feature])
        updateSelectizeInput(session, 'sc_celltype22', 
                             choices = b_celltypes, 
                             server = TRUE,selected =NULL)
        cat("\n")
        scProg$set(message = "done..", value = 1)
        
      }
    })
  
  
}


