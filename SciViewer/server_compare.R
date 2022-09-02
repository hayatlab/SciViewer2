## server_compare.R

##-----------------------------------------------------------------
##    gene expression across available studies
##-----------------------------------------------------------------

observeEvent(c(compvariables$geneExpSummarized),{
  if(isTRUE(nrow(compvariables$geneExpSummarized)>0)){
    ## data table
    cat(" displaying summary expression table across all available studies for gene ",compvariables$gene,'\n')
    output$comparativeExpnTable <- DT::renderDataTable({
      pf <- compvariables$geneExpSummarized
      pf$AveExpr <- round(pf$AveExpr,2)
      pf <- pf[order(pf$AveExpr,decreasing = T),
               ]
      datatable(pf,rownames = F,extensions = 'Buttons',
                options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                               buttons = list(list(extend = 'csv',  filename = paste0("ExpressionComparison_",compvariables$gene)),
                                              list(extend ='excel', filename = paste0("ExpressionComparison_",compvariables$gene))), 
                               scrollX = TRUE,scrollCollapse = TRUE)) %>% 
        formatStyle('AveExpr',
                    background = styleColorBar(range(pf$AveExpr,na.rm=T), 'lightblue'),
                    backgroundRepeat = 'no-repeat',backgroundPosition = 'left')
    })
  }
})


observeEvent(c(input$comparativeExpnTable_rows_selected),{
  if(isTruthy(input$comparativeExpnTable_rows_selected) & isTRUE(nrow(compvariables$geneExpSummarized)>0)){
   
    output$comparativeExpnPlot <- renderUI({
      output$tmpl <- renderPlotly({
        db <- compvariables$geneExpSummarized[input$comparativeExpnTable_rows_selected,]$Database %>%
          as.character()
        pl <- compvariables$geneExp
        pl <- pl[pl$Database%in%db,]
        if(nrow(pl)>0){
          P <- pl %>%
            ggplot(aes(x=cell_type,y=norm_avg,fill=cell_type))+
            geom_boxplot(alpha=0.2) + facet_wrap(~Database,ncol=2,scales = 'free')+
            theme_bw()+
            theme(legend.position = 'none',strip.background = element_blank())+
            xlab("")+ylab("")+
            coord_flip()
          P %>% ggplotly()
        }else{
          empty_plot("select row from above table")
        }
      })
      plotlyOutput("tmpl")
    })
  }
})


##-----------------------------------------------------------------
##    gene summarized as marker across studies
##-----------------------------------------------------------------

observeEvent(c(compvariables$geneMarker),{
  if(isTRUE(nrow(compvariables$geneMarker)>0)){
    ## data table
    cat(" displaying summary of gene as marker across all available studies for gene ",compvariables$gene,'\n')
    output$comparativeMarkerTable <- DT::renderDataTable({
      cf <- compvariables$geneMarker
      cf$logFC = round(cf$logFC,2)
      cf$AveExpr = round(cf$AveExpr,2)
      cf$t = round(cf$t,2)
      cf$B = round(cf$B,2)
      cf <- cf[,c("Database","geneSymbol", "Test","logFC", "adj.P.Val")]
      datatable(cf,rownames = F,extensions = 'Buttons',
                options = list(searching = TRUE,pageLength = 5,dom = 'Bfrtip', 
                               buttons = list(list(extend = 'csv',  filename = paste0("Markers_FDR_CrossStudies")),
                                              list(extend ='excel', filename = paste0("Markers_FDR_CrossStudies"))), 
                               scrollX = TRUE,scrollCollapse = TRUE
                )) %>% 
        formatSignif(columns = c('adj.P.Val'),digits = 2) %>%
        formatStyle('logFC',
                    background = styleColorBar(range(cf$logFC,na.rm=T), 'lightblue'),
                    backgroundRepeat = 'no-repeat', backgroundPosition = 'right')
    })
  }else{
    output$comparativeMarkerTable <- DT::renderDataTable({c()})
  }
})







##-----------------------------------------------------------------
##    gene summarized as biological marker across studies
##-----------------------------------------------------------------


observeEvent(c(compvariables$geneBioMarker),{
  if(isTRUE(nrow(compvariables$geneBioMarker)>0)){
    ## data table
    cat(" displaying summary of gene as biomarker across all available studies for gene ",compvariables$gene,'\n')
    output$comparativeBioMarkerTable <- DT::renderDataTable({
      cf <- compvariables$geneBioMarker
      cf$logFC = round(cf$logFC,2)
      cf$AveExpr = round(cf$AveExpr,2)
      cf$t = round(cf$t,2)
      cf$B = round(cf$B,2)
      cf <- cf[,c("Database","geneSymbol", "Test","logFC", "adj.P.Val")]
      datatable(cf,rownames = F,extensions = 'Buttons',
                options = list(searching = TRUE,pageLength = 5,dom = 'Bfrtip', 
                               buttons = list(list(extend = 'csv',  filename = paste0("BioMarkers_FDR_CrossStudies")),
                                              list(extend ='excel', filename = paste0("BioMarkers_FDR_CrossStudies"))), 
                               scrollX = TRUE,scrollCollapse = TRUE
                )) %>% 
        formatSignif(columns = c('adj.P.Val'),digits = 2) %>%
        formatStyle('logFC',
                    background = styleColorBar(range(cf$logFC,na.rm=T), 'lightblue'),
                    backgroundRepeat = 'no-repeat', backgroundPosition = 'right')
    })
  }else{
    output$comparativeBioMarkerTable <- DT::renderDataTable({NULL})
  }
})

