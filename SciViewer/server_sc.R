## server_sc.R

##--------------------------------------------------------------
##--- Biological Marker expression
##--------------------------------------------------------------
output$dedisplayselcelltype <- renderText({ 
  input$sel_celltype[1]
})
observe({
  updateSelectizeInput(session, 'sel_decelltype',choices = variables$availcelltypes,
                       server = TRUE,selected =NULL)
})
observe({
  updateSelectizeInput(session, 'sel_decomp',
                       choices = variables$avail_deTests, #variables$sel_catvar
                       selected =NULL,server=TRUE)
})
observeEvent(c(input$sel_genego),{
  if(isTruthy(variables$gene) & isTruthy(input$sel_genego) & isTRUE(nrow(variables$de_geneDf)>0)){
    cat(" preparing the marker dotplot\n")
    output$sc_dedotplot <- renderUI({ 
      output$tempj <- renderPlotly({
        sc_dge_dotplot(variables$de_geneDf,x='geneSymbol',y='Test',
                       colCol='logFC',sizeCol='adj.P.Val',testCol='Test',configlayout=T)
      })
      plotlyOutput("tempj")
    })
  }else{
    output$sc_dedotplot <- renderUI({
      output$tempj <- renderPlotly({empty_plot("biological marker information unavailable")})
      plotlyOutput("tempj")
    })
  }
})

observeEvent(c(variables$de_cellTypeDf,input$sel_defdrslider,input$sel_decomp),{ 
  if(isTRUE(nrow(variables$de_cellTypeDf)>0) & isTruthy(input$sel_defdrslider) &
     isTruthy(input$sel_decomp)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("displaying biological marker data table for ",input$sel_decelltype), value = 0)
    cat("displaying biological marker data table for ",input$sel_decelltype,"\n")
    
    ## Spreadsheet
    output$de_tab = DT::renderDataTable({
      cf <- variables$de_cellTypeDf
      cf <- cf[as.numeric(cf$adj.P.Val) <= as.numeric(input$sel_defdrslider),]
      cf <- cf[cf$Test==input$sel_decomp,]
      if(nrow(cf)>0){
        cf$logFC = round(cf$logFC,2)
        cf$AveExpr = round(cf$AveExpr,2)
        cf$t = round(cf$t,2)
        cf$B = round(cf$B,2)
        cf <- cf[,c("geneSymbol","logFC", "adj.P.Val", "Test")]
        datatable(cf,rownames = F,extensions = 'Buttons',
                  options = list(searching = TRUE,pageLength = 5,dom = 'Bfrtip', 
                                 buttons = list(list(extend = 'csv',  filename = paste0("BioMarkers_FDR_",input$sel_defdrslider,"_",variables$sc_study,"_",input$sel_decelltype,"")),
                                                list(extend ='excel', filename = paste0("BioMarkers_FDR_",input$sel_defdrslider,"_",variables$sc_study,"_",input$sel_decelltype,""))), 
                                 scrollX = TRUE,scrollCollapse = TRUE
                  )) %>% formatSignif(columns = c('adj.P.Val'),digits = 2)
      }
    })
    
    ## Volcano plot
    scProg$set(message = paste0("displaying volcano plot for ",input$sel_decelltype), value = 0.5)
    cat("displaying volcano plot ",input$sel_decelltype,"\n")
    output$sc_deVolcano <- renderPlotly({
      cf <- variables$de_cellTypeDf
      cf <- cf[cf$Test==input$sel_decomp,]
      cf <- cf[order(cf[,"adj.P.Val"],decreasing = F),]
      cf <- cf[1:min(c(nrow(cf),500)),]
      pddExpn::sc_dge_volcanoplot(pl=cf)
    })
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(variables$expndf_de),{
  if(isTRUE(ncol(variables$expndf_de)>0)){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("displaying biological marker gene expression plots"), value = 0)
    cat("displaying biological marker gene expression plots. ",ncol(variables$expndf_de),"\n")
    
    #-- multigene umap plot
    output$sc_deumap <- renderUI({ 
      myL <- local({
        myL <- list()
        for(i in 1:ncol(variables$expndf_de)){
          cat(i,": multigene umap plot: ", ncol(variables$expndf_de),"\n")
          vals <- variables$expndf_de[,i]
          cf <- data.frame(V1=variables$sc_louvain[,'V1'],
                           V2=variables$sc_louvain[,'V2'],
                           value=ifelse(!is.finite(vals),0,vals))
          myL[[names(variables$expndf_de)[i]]] <- cf
        }
        myL
      })
      output$tmpk <- renderPlotly({
        plot_multi_umap_raster(myL,genenames = names(variables$expndf_de),legendPlot = F)
      })
      plotlyOutput("tmpk")
    })
    scProg$set(message = "done..", value = 1)
  }else{
    output$tmpk <- renderPlotly({ empty_plot("") })
    plotlyOutput("tmpk")
  }
})

observeEvent(c(variables$cellde_upenrich,variables$cellde_dnenrich, 
               input$sc_deopgo,input$sel_dedirchange),{
                 if(input$AppTab=='scstudy' & isTruthy(variables$cellde_upenrich) & isTruthy(input$sel_dedirchange) &
                    isTruthy(variables$cellde_dnenrich) & isTruthy(input$sel_dedirchange)){
                   
                   scProg <- shiny::Progress$new()
                   on.exit(scProg$close())
                   scProg$set(message = paste0("displaying biological marker gene enrichment plot"), value = 0)
                   cat("displaying biological marker gene enrichment plot\n")
                   
                   rr1 <- reactive({ variables$cellde_upenrich})
                   if(input$sel_dedirchange==2){
                     rr1 <- reactive({ variables$cellde_dnenrich})
                   }
                   
                   output$enrichde_plot <- renderUI({ 
                     if(isTruthy(rr1())){
                       output$ewascjyxy <- renderPlotly({
                         gostplot(gostres = rr1(),capped = T,interactive = T)
                       })
                       plotlyOutput("ewascjyxy")
                     }else{
                       renderPlotly({empty_plot("information unavailable")})
                     }  
                   })
                   output$enrichde_tab <- DT::renderDataTable({ 
                     if(isTruthy(rr1())){
                       rr1()[[1]] %>%
                         select(c("source","term_id" ,"term_name","term_size", "intersection_size","p_value")) %>%
                         filter(intersection_size>1 & p_value<0.05) %>%
                         arrange(p_value) %>%
                         datatable(rownames = F,extensions = 'Buttons',
                                   options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                                                  buttons = list(list(extend = 'csv',  filename = paste0("Enrichments_FDR")),
                                                                 list(extend ='excel', filename = paste0("Enrichments_FDR"))), 
                                                  scrollX = TRUE,scrollCollapse = TRUE)) %>% 
                         formatSignif(columns = c('p_value'),digits = 2)
                     }else{ DT::renderDataTable({}) }  
                   })
                   output$deEnrichmentText <- renderText({ rr1()[['info']] })
                   scProg$set(message = "done..", value = 1)
                 }
               })


##--------------------------------------------------------------
##--- Marker expression
##--------------------------------------------------------------
output$markerdisplayselcelltype <- renderText({ 
  input$sel_celltype[1]
})
observe({
  updateSelectizeInput(session, 'sel_markercelltype',choices = variables$availcelltypes,
                       server = TRUE,selected =NULL)
})
observe({
  updateSelectizeInput(session, 'sel_markercomp',
                       choices = variables$avail_markerTests, #variables$sel_catvar
                       selected =NULL,server=TRUE)
})
observeEvent(c(input$sel_genego),{
    if(isTruthy(variables$gene) & isTruthy(input$sel_genego) & isTRUE(nrow(variables$marker_geneDf)>0)){
      cat(" preparing the marker dotplot\n")
      output$sc_markerdotplot <- renderUI({ 
        output$tempe <- renderPlotly({
          sc_dge_dotplot(variables$marker_geneDf,x='geneSymbol',y='Test',
                         colCol='logFC',sizeCol='adj.P.Val',testCol='Test',configlayout=T)
        })
        plotlyOutput("tempe")
      })
    }else{
      output$sc_markerdotplot <- renderUI({
        output$tempe <- renderPlotly({empty_plot("marker information unavailable")})
        plotlyOutput("tempe")
      })
    }
  })

observeEvent(c(variables$marker_cellTypeDf,input$sel_markerfdrslider,input$sel_markercomp),{ 
  if(isTRUE(nrow(variables$marker_cellTypeDf)>0) & isTruthy(input$sel_markerfdrslider) &
     isTruthy(input$sel_markercomp)){
    
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("displaying marker data table for ",input$sel_markercelltype), value = 0)
    cat("displaying marker data table for ",input$sel_markercelltype,"\n")
    
    ## Spreadsheet
    output$marker_tab = DT::renderDataTable({
      cf <- variables$marker_cellTypeDf
      cf <- cf[as.numeric(cf$adj.P.Val) <= as.numeric(input$sel_markerfdrslider),]
      cf <- cf[cf$Test==input$sel_markercomp,]
      if(nrow(cf)>0){
        cf$logFC = round(cf$logFC,2)
        cf$AveExpr = round(cf$AveExpr,2)
        cf$t = round(cf$t,2)
        cf$B = round(cf$B,2)
        cf <- cf[,c("geneSymbol","logFC", "adj.P.Val", "Test")]
        datatable(cf,rownames = F,extensions = 'Buttons',
                  options = list(searching = TRUE,pageLength = 5,dom = 'Bfrtip', 
                                 buttons = list(list(extend = 'csv',  filename = paste0("Markers_FDR_",input$sel_markerfdrslider,"_",variables$sc_study,"_",input$sel_markercelltype,"")),
                                                list(extend ='excel', filename = paste0("Markers_FDR_",input$sel_markerfdrslider,"_",variables$sc_study,"_",input$sel_markercelltype,""))), 
                                 scrollX = TRUE,scrollCollapse = TRUE
                  )) %>% formatSignif(columns = c('adj.P.Val'),digits = 2)
      }
    })
    
    ## Volcano plot
    scProg$set(message = paste0("displaying volcano plot for ",input$sel_markercelltype), value = 0.5)
    cat("displaying volcano plot ",input$sel_markercelltype,"\n")
    output$sc_markerVolcano <- renderPlotly({
      cf <- variables$marker_cellTypeDf
      cf <- cf[cf$Test==input$sel_markercomp,]
      cf <- cf[order(cf[,"adj.P.Val"],decreasing = F),]
      cf <- cf[1:min(c(nrow(cf),500)),]
      pddExpn::sc_dge_volcanoplot(pl=cf)
    })
    scProg$set(message = "done..", value = 1)
  }
})

observeEvent(c(variables$expndf_markers),{
  if(isTRUE(ncol(variables$expndf_markers)>0)){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("displaying marker gene expression plots"), value = 0)
    cat("displaying marker gene expression plots. ",ncol(variables$expndf_markers),"\n")

    #-- multigene umap plot
    output$sc_markerumap <- renderUI({ 
      myL <- local({
        myL <- list()
        for(i in 1:ncol(variables$expndf_markers)){
          cat(i,": multigene umap plot: ", ncol(variables$expndf_markers),"\n")
          vals <- variables$expndf_markers[,i]
          cf <- data.frame(V1=variables$sc_louvain[,'V1'],
                           V2=variables$sc_louvain[,'V2'],
                           value=ifelse(!is.finite(vals),0,vals))
          myL[[names(variables$expndf_markers)[i]]] <- cf
        }
        myL
      })
      output$tmpf <- renderPlotly({
        plot_multi_umap_raster(myL,genenames = names(variables$expndf_markers),legendPlot = F)
      })
      plotlyOutput("tmpf")
    })
    scProg$set(message = "done..", value = 1)
  }else{
    output$tmpf <- renderPlotly({ empty_plot("") })
    plotlyOutput("tmpf")
  }
})

observeEvent(c(variables$cellmark_upenrich,variables$cellmark_dnenrich, 
               input$sc_markeropgo,input$sel_markerdirchange),{
    if(input$AppTab=='scstudy' & isTruthy(variables$cellmark_upenrich) & isTruthy(input$sel_markerdirchange) &
       isTruthy(variables$cellmark_dnenrich) & isTruthy(input$sel_markerdirchange)){
      
      scProg <- shiny::Progress$new()
      on.exit(scProg$close())
      scProg$set(message = paste0("displaying marker gene enrichment plot"), value = 0)
      cat("displaying marker gene enrichment plot\n")
      
      rr1 <- reactive({ variables$cellmark_upenrich})
      if(input$sel_markerdirchange==2){
        rr1 <- reactive({ variables$cellmark_dnenrich})
      }
      
      output$enrichmarkers_plot <- renderUI({ 
        if(isTruthy(rr1())){
          output$ewascjyx <- renderPlotly({
            gostplot(gostres = rr1(),capped = T,interactive = T)
          })
          plotlyOutput("ewascjyx")
        }else{
          renderPlotly({empty_plot("information unavailable")})
        }  
      })
      output$enrichmarker_tab <- DT::renderDataTable({ 
        if(isTruthy(rr1())){
          rr1()[[1]] %>%
            select(c("source","term_id" ,"term_name","term_size", "intersection_size","p_value")) %>%
            filter(intersection_size>1 & p_value<0.05) %>%
            arrange(p_value) %>%
            datatable(rownames = F,extensions = 'Buttons',
                      options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                                     buttons = list(list(extend = 'csv',  filename = paste0("Enrichments_FDR")),
                                                    list(extend ='excel', filename = paste0("Enrichments_FDR"))), 
                                     scrollX = TRUE,scrollCollapse = TRUE)) %>% 
            formatSignif(columns = c('p_value'),digits = 2)
        }else{ DT::renderDataTable({}) }  
      })
      output$markerEnrichmentText <- renderText({ rr1()[['info']] })
      scProg$set(message = "done..", value = 1)
    }
  })


##--------------------------------------------------------------
##--- Multi gene expression
##--------------------------------------------------------------
observe({
  updateSelectizeInput(session, 'sc_catvar_expndot',
                       choices = c(variables$celltype,variables$donor,variables$catvars), #variables$sel_catvar
                       selected =variables$celltype,server=TRUE)
})
observe({
  updateSelectizeInput(session, 'sc_corrcelltype',choices = variables$availcelltypes,
                       server = TRUE,selected =NULL)
})

observeEvent(c(variables$gene,variables$expndf),{
  if(isTruthy(variables$gene) & isTruthy(variables$expndf)){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("displaying multi gene expression plot"), value = 0)
    cat("displaying multi gene expression plot\n")
    
    output$sc_multigeneumap <- renderUI({ 
      myL <- local({
        myL <- list()
        for(i in 1:ncol(variables$expndf)){
          vals <- variables$expndf[,i]
          cf <- data.frame(V1=variables$sc_louvain[,'V1'],
                           V2=variables$sc_louvain[,'V2'],
                           value=ifelse(!is.finite(vals),0,vals))
          myL[[names(variables$expndf)[i]]] <- cf
        }
        myL
      })
      
      output$tmpa <- renderPlotly({
        plot_multi_umap_raster(myL,genenames = names(variables$expndf),legendPlot = F)
      })
      plotlyOutput("tmpa")
    })
    scProg$set(message = "done..", value = 1)
    
  }
})

observeEvent((input$sc_catvar_expndotgo),{
  if(isTruthy(input$sc_catvar_expndotgo) & isTruthy(input$sc_catvar_expndot) & isTruthy(variables$gene) & isTruthy(variables$expndf)){
    scProg <- shiny::Progress$new()
    on.exit(scProg$close())
    scProg$set(message = paste0("displaying gene expression dot plot"), value = 0)
    cat("displaying gene expression dot plot\n")
    sc_expndotplot <- reactive({
      gnames <- gsub(paste0(variables$sc_study,"_"),"",names(variables$expndf))
      myL <- local({
        myL <- list()
        for(i in 1:ncol(variables$expndf)){
          vals <- variables$expndf[,i]
          cf <- data.frame(V1=variables$sc_louvain[,'V1'],
                           V2=variables$sc_louvain[,'V2'],
                           value=vals)
          cf <- cbind(variables$sc_louvain[,c(variables$catvars,variables$celltype,variables$donor)],cf)
          myL[[names(variables$expndf)[i]]] <- cf
        }
        myL
      })
      plot_multigene_grouped_heatmap01(myL,
                                     genenames = gnames,
                                     feature = input$sc_catvar_expndot,
                                     x = 'V1',y = 'V2',ycol = 'value',
                                     log.transform = F)
    })
    output$sc_expndotplot <- renderUI({ 
      output$tempc <- renderPlotly(sc_expndotplot()[[2]])
      plotlyOutput("tempc",width = '600px',height = '450px')
    })
    scProg$set(message = "done..", value = 1)
    
  }
})


observeEvent(c(variables$expndf,input$sc_corrcelltypego),{
    if(input$AppTab=='scstudy' & isTruthy(input$sc_corrcelltypego) & 
       isTruthy(variables$expndf) & isTruthy(variables$sc_louvain) & isTRUE(ncol(variables$expndf)>1)){
      
      scProg <- shiny::Progress$new()
      on.exit(scProg$close())
      scProg$set(message = paste0("displaying correlation plot"), value = 0)
      cat("displaying correlation plot\n")
      
      z <- cbind(celltype=variables$sc_louvain[,variables$celltype],variables$expndf)
      z <-  z[rowSums(is.na(z[,-1]))!=(ncol(z)-1),]
      if(isTruthy(input$sc_corrcelltype)){
        z <- z[z[,1]%in%input$sc_corrcelltype,]
      }
      output$sc_corrplot <- renderUI({ 
        hc <- hchart_cor(cor(z[,-1],use = 'pairwise.complete.obs'),sourceref = SOURCEREF)
        rows <- max(c(300,ncol(z[,-1]) *100))
        rows <- paste0(rows,"px")
        output$tempd <- renderHighchart(hc)
        highchartOutput("tempd")
      })
      scProg$set(message = "done..", value = 1)
      
    }
})



##--------------------------------------------------------------
##--- Single gene expression
##--------------------------------------------------------------
observe({
  updateSelectizeInput(session, 'sc_gcolscale', 
                       choices = c("GrayBlue","Clay","magma","Viridis",
                                   "GreenYellow","Purple","Reds","Oranges"), 
                       server = TRUE,selected =NULL)
})
observe({
  updateSelectizeInput(session, 'sc_expncatvar',
                       choices = c(variables$celltype,variables$donor,variables$catvars), #variables$sel_catvar
                       selected =variables$celltype,server=TRUE)
})
observe({
  updateSelectizeInput(session, 'sc_catvar_expncount',
                       choices = c(variables$celltype,variables$donor,variables$catvars), #variables$sel_catvar
                       selected =variables$celltype,server=TRUE)
})

observeEvent(c(variables$expndf, input$sc_gcolscale),{
    if(input$AppTab=='scstudy' & isTruthy(input$sel_genego) & 
       isTruthy(variables$expndf) & isTruthy(variables$gene) & isTruthy(input$sc_gcolscale)){
      
      scProg <- shiny::Progress$new()
      on.exit(scProg$close())
      scProg$set(message = paste0("displaying gene expression for ", variables$gene), value = 0)
      cat("displaying gene expression for ", variables$gene,"\n")
      
      myL <- local({
        myL <- list()
        for(i in 1){
          vals <- variables$expndf[,i]
          cf <- data.frame(V1=variables$sc_louvain[,'V1'],
                           V2=variables$sc_louvain[,'V2'],
                           value=ifelse(!is.finite(vals),0,vals))
          myL[[names(variables$expndf)[i]]] <- cf
        }
        myL
      })
      if(names(variables$expndf)[1]!=variables$gene[1]){
        stop("the gene order on expndf and variable$gene does not match\n")
      }
      output$expnumapplot <- renderPlotly({ 
        ql <- select_colScale(input$sc_gcolscale)
        suppressWarnings({
          plot_multi_umap_raster(myL,
                                 genenames = variables$gene[1],
                                 legendPlot = F,
                                 colScale=ql[[1]],
                                 minValCol=ql[[2]]) %>%
            layout(height = 300, width = 350)
        })
      })
      scProg$set(message = "done..", value = 1)
      
    }
  })

observeEvent(c(input$sel_sgexpngo,input$sc_expncatvar,input$rad_expngraphchoice,input$rad_expnptgrp,input$rad_expnchoice),{
    if(input$AppTab=='scstudy' & isTruthy(input$sc_expncatvar) & isTruthy(input$sel_sgexpngo) &
       isTruthy(input$rad_expngraphchoice) & isTruthy(input$rad_expnptgrp) & isTruthy(input$rad_expnchoice)){
      
      cat("generating gene expression summary using barplots\n")
      
      output$sc_gebarplot <- renderHighchart({
        
        myL <- local({
          myL <- list()
          for(i in 1){
            vals <- variables$expndf[,i]
            cf <- data.frame(V1=variables$sc_louvain[,'V1'],
                             V2=variables$sc_louvain[,'V2'],
                             value=vals)
            cf <- cbind(variables$sc_louvain[,c(variables$catvars,variables$celltype,variables$donor)],cf)
            myL[[names(variables$expndf)[i]]] <- cf
          }
          myL
        })
        pl = myL[[1]]
        if(input$rad_expnchoice==2){
          pl <- pl[!is.na(pl$value),]
        }else{
          pl$value <- ifelse(is.na(pl$value),0,pl$value)
        }
        features = input$sc_expncatvar
        min.cell.number = 10
        ycol = 'value'
        ylab = "Normalized expression"
        sampleType = 'cells'
        if(input$rad_expnptgrp==2 & isTruthy(variables$donor)){
          features = c(input$sc_expncatvar[1],variables$donor)
          ycol = 'value'
          pl <- pl %>% dplyr::group_by_at(features) %>% dplyr::select(dplyr::all_of(ycol)) %>% 
            dplyr::summarise_all(c("avg_expn"=mean))
          features = input$sc_expncatvar[1]
          min.cell.number = 0
          ycol = 'avg_expn'
          ylab = "Average expression per donor"
          sampleType = 'donors'
        }
        if(input$rad_expngraphchoice==1){
          dhc_columnPlot(pl=pl,
                         features = features,
                         plotType = 'expn',min.cell.number = min.cell.number,
                         ycol = ycol,
                         log.transform = F,main = paste0(variables$gene[1]," expression"),
                         xlab = "",ylab = ylab,sampleType = sampleType,
                         sourceref=SOURCEREF)
        }else if(input$rad_expngraphchoice==2){
          dhc_boxPlot(pl=pl,
                      features = features,
                      ycol = ycol,
                      log.transform = F,
                      main = paste0(variables$gene[1]," expression"),
                      xlab = "",ylab = ylab,
                      sourceref=SOURCEREF)
        }
      })
    }
  })
observeEvent({input$sc_catvar_expncount},{
  if(isTruthy(variables$sc_louvain) & isTruthy(input$sc_catvar_expncount)){
    output$sc_expncountbar <- renderHighchart({
      
      myL <- local({
        myL <- list()
        for(i in 1){
          vals <- variables$expndf[,i]
          cf <- data.frame(V1=variables$sc_louvain[,'V1'],
                           V2=variables$sc_louvain[,'V2'],
                           value=vals)
          cf <- cbind(variables$sc_louvain[,c(variables$catvars,variables$celltype,variables$donor)],cf)
          myL[[names(variables$expndf)[i]]] <- cf
        }
        myL
      })
      pl = myL[[1]]
      pl <- pl[!is.na(pl$value),]
      dhc_columnPlot(pl=pl,
                     features = input$sc_catvar_expncount,
                     ycol = 'value',plotType = 'count',
                     log.transform = F,main = paste0("% of 'expressing' cells <br>",variables$gene[1],""),
                     xlab = "",ylab = "% of 'total'expressing' cells",
                     sourceref=SOURCEREF)
    })
  }
})





##--------------------------------------------------------------
##--- Base UMAP and selections 
##--------------------------------------------------------------
 ## issues A) HOVER does not work!
 ##        B) Categorical variable value selection does not work!
observeEvent({variables$sc_study_attribs},{
  if(isTruthy(variables$sc_study_attribs)){
    output$sc_study_attribs <- renderText({ variables$sc_study_attribs })
    output$sc_study_attribs01 <- renderText({ variables$sc_study_attribs })
    
  }
})
observe({
  updateSelectizeInput(session, 'sel_celltype',choices = variables$availcelltypes,
                       server = TRUE,selected =NULL)
  updateSelectizeInput(session, 'sel_donor',choices =unique(variables$sc_louvain[,variables$donor]),
                       server = TRUE,selected =NULL)
  if(!is.null(variables$catvars)){
    updateSelectizeInput(session, 'sel_catvars',choices = variables$catvars,
                         server = TRUE,selected =NULL)
  }
})
observe({
  updateSelectizeInput(session, 'sel_catvarsval',
                       choices = variables$catvars_valuelist[[input$sel_catvars]], #variables$sel_catvar
                       selected =NULL,server=TRUE)
})
observe({
  cat("01.variables$sel_catvarval: ",input$sel_catvarsval,"\t(",input$sel_catvars,")\n")
})

observeEvent(c(input$sel_metadatago),{
  if(isTruthy(variables$sc_louvain) & isTruthy(input$sel_metadatago) & isTruthy(input$sel_catvars)){
    cat("02.variables$sel_catvarval: ",input$sel_catvarsval,"\t(",input$sel_catvars,")\n")
    umapplot <- reactive({
      plot_umap_raster(pl=variables$sc_louvain,
                       feature = input$sel_catvars,
                       genename = paste0(input$sel_catvars," clusters"),
                       legendPlot = F,
                       selectedFeature = NULL,#input$sel_catvarsval,
                       source.val='umap3',js=js)
    })
    output$umapplot_hover <- renderText({
      getHoverText(variables$sc_louvain,colids=c('V1','V2',input$sel_catvars),
                   hovercoord = input$umap3 )
    })
    output$umapplot <- renderPlotly({ 
      suppressWarnings({umapplot() %>% layout(height = 450, width = 450)})
    })
  }
})

observeEvent(c(variables$sc_louvain,input$sel_celltypego),{
  if(isTruthy(variables$sc_louvain)){
    cat("rendering UMAP cluster plot\n")
    umapplot <- reactive({
      plot_umap_raster(pl=variables$sc_louvain,
                       feature = variables$celltype,
                       genename = "Cell type clusters",
                       legendPlot = F,
                       selectedFeature = input$sel_celltype,
                       source.val='umap1',js=js)
    })
    output$umapplot_hover <- renderText({
      getHoverText(variables$sc_louvain,colids=c('V1','V2',variables$celltype),
                   hovercoord = input$umap1 )
    })
    cat("hover::",input$umap1,"\n")
    output$umapplot <- renderPlotly({ 
      suppressWarnings({
        umapplot() %>% layout(height = 450, width = 450)
      })
    })
  }
})
observeEvent(c(variables$sc_louvain,input$sel_donorgo),{
  if(isTruthy(variables$sc_louvain) & isTruthy(input$sel_donorgo)){
    cat("rendering UMAP cluster plot for donors\n")
    umapplot <- reactive({
      plot_umap_raster(pl=variables$sc_louvain,
                       feature = variables$donor,
                       genename = "Donor clusters",
                       legendPlot = F,
                       selectedFeature = input$sel_donor,
                       source.val='umap2',js=js)
    })
    output$umapplot_hover <- renderText({
      getHoverText(variables$sc_louvain,colids=c('V1','V2',variables$donor),
                   hovercoord = input$umap2 )
    })
    output$umapplot <- renderPlotly({ 
      suppressWarnings({
        umapplot() %>% layout(height = 450, width = 450)
      })
    })
  }
})

##--------------------------------------------------------------
##--- QC plots 
##--------------------------------------------------------------
observe({
  updateSelectizeInput(session, 'qc_attrib01',
                       choices = c(variables$celltype,variables$donor,variables$catvars), #variables$sel_catvar
                       selected =variables$celltype,server=TRUE)
})
observe({
  updateSelectizeInput(session, 'qc_attrib02',
                       choices = c(variables$contvars), #variables$sel_catvar
                       selected =variables$contvars[1],server=TRUE)
})
observeEvent({input$qc_attrib01},{
  if(isTruthy(variables$sc_louvain) & isTruthy(input$qc_attrib01)){
    output$sc_qcbarplot <- renderHighchart({
      dhc_columnPlot(pl=variables$sc_louvain,
                     features = input$qc_attrib01,
                     ycol = 'value',plotType = 'count',
                     log.transform = F,main = paste0("Cell type proportions<br>(",input$qc_attrib01,")"),
                     xlab = "",ylab = "% of total",
                     sourceref=SOURCEREF)
    })
  }
})
observeEvent({input$qc_attrib02},{
  if(isTruthy(variables$sc_louvain) & isTruthy(input$qc_attrib02)){
    output$sc_qcscatter <- renderUI({
      tmpL <- list()
      for(i in 1:length(input$qc_attrib02)){
        tmpL[[i]] <- variables$sc_louvain
        tmpL[[i]]$value <- tmpL[[i]][,input$qc_attrib02[i]]
        tmpL[[i]]$value <- as.numeric(as.character( tmpL[[i]]$value))
      }
      output$tmpa <- renderPlotly({
        plot_multi_umap_raster(tmpL,
                               genenames = input$qc_attrib02,
                               legendPlot = F
                               #xrange =range2a()[[1]],
                               #yrange=range2a()[[2]],
                               #colScale=ql[[1]],
                               #minValCol=ql[[2]]
        )
      })
      plotlyOutput("tmpa")
    })
  }
})






