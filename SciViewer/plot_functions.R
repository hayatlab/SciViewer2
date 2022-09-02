## following functions are highly customized, they should become part of pddExpn eventually
queryDB <- function(HANDLER=NULL, QUERY=NULL,
                    REPO_NAME=NULL,USE_REMOTE_DB=FALSE){
  
  if(is.null(HANDLER) | is.null(QUERY)){
    return(NULL)
  }
  if(USE_REMOTE_DB==FALSE){
    return(RSQLite::dbGetQuery(HANDLER, QUERY)) 
  }
  else if(USE_REMOTE_DB==TRUE){
    #QUERY = paste0("SELECT * FROM    dfg WHRE x = 'from'")
    #HANDLER= 'qqq'
    #gsub("FROM\\s+",paste0("FROM ",HANDLER,"."),QUERY)
    QUERY = gsub("FROM\\s+",paste0("FROM ",HANDLER,"."),QUERY)
    return(GetData(AppName=REPO_NAME, sSQL=QUERY))
  }
  else{
    return(NULL)
  }
}


empty_plot <- function(title = NULL){
  p <- plotly_empty(type = "scatter", mode = "markers") %>%
    layout(
      title = list(text = title,y = 0.2),
      xaxis = list(showline = TRUE,mirror = "ticks",linecolor = toRGB("black"),linewidth = 1,range = c(-3,3)),
      yaxis = list(showline = TRUE,mirror = "ticks",linecolor = toRGB("black"),linewidth = 1,range = c(-3,3))
    )
  return(p)
} 

get_marker_df <- function(connSc,study){
  
  markers <- matrix(NA,nrow = 0,ncol=11)
  markers <- as.data.frame(markers)
  names(markers) <- c("cluster", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", 
                      "gene", "cell_type", "test", "id", "top10")
  
  query <- paste0('SELECT * FROM ',study,"_markers")
  cf <- tryCatch({
    queryDB(HANDLER=connSc, QUERY=query,REPO_NAME=REPO_NAME,USE_REMOTE_DB=USE_REMOTE_DB)
    #RSQLite::dbGetQuery(connSc,query)
  },error = function(e){
    cf <- markers
    return(cf)
  })
  return(cf)
  
}

hchart_cor <- function(mcor,sourceref='SciViewer') {
  
  mcor <- round(mcor,2)
  mcor[is.na(mcor)] = 0
  
  #fntltp <- JS("function(){
  #                return this.series.xAxis.categories[this.point.x] + ' ~ ' +
  #                       this.series.yAxis.categories[this.point.y] + ': <b>' +
  #                       Highcharts.numberFormat(this.point.value, 2)+'</b>';
  #             ; }")
  cor_colr <- list( list(0, '#FF5733'),
                    list(0.5, '#F8F5F5'),
                    list(1, '#2E86C1')
  )
  
  hs <- hchart(mcor) %>%
    hc_plotOptions(
      series = list(
        boderWidth = 0,
        dataLabels = list(enabled = TRUE)
      ))   %>%
    hc_legend(align = "right", layout = "vertical") %>%
    hc_colorAxis(  stops= cor_colr,min=-1,max=1)
  hs <- hs %>%
    hc_credits(
      enabled = TRUE,
      text = paste0("Source:",sourceref),
      href = "",
      style = list(fontSize = "11px")
    )%>%
    hc_exporting(enabled=T
    ) %>%
    hc_scrollbar(
      barBackgroundColor = "gray",
      barBorderRadius = 7,
      barBorderWidth = 0,
      buttonBackgroundColor = "gray",
      buttonBorderWidth = 0,
      buttonArrowColor = "yellow",
      buttonBorderRadius = 7,
      rifleColor = "yellow",
      trackBackgroundColor = "white",
      trackBorderWidth = 1,
      trackBorderColor = "silver",
      trackBorderRadius = 7
    ) %>%
    hc_chart(zoomType = "xy",
             borderColor = "#EBBA95",
             borderRadius = 10,
             borderWidth = 3,
             plotShadow = F,
             allowForce = T,
             turboThreshold =1,
             allowForce=T,
             animation=T,
             boostThreshold = 1,
             usePreallocated = T,
             useGPUTranslations =T,
             seriesThreshold = 2
    )
  
  return(hs)
}


get_plot_df_sc04 <- function(connSc,study="Madisson_LungTissue",
                             genename="WNT4",louvain,REPO_NAME=NULL,USE_REMOTE_DB=FALSE){
  #louvain must have a value column
  if(!'value'%in%names(louvain)){
    louvain$value <- NA
  }
  
  query <- paste0("SELECT * FROM ",study,"_data WHERE geneSymbol = '", genename,"'")
  gdf <- queryDB(HANDLER=connSc, 
                 QUERY=query,REPO_NAME=REPO_NAME,
                 USE_REMOTE_DB=USE_REMOTE_DB)
  
  if(nrow(gdf)>0){
    cnta <- matrix(NA,nrow = nrow(louvain),ncol=nrow(gdf)) %>% as.data.frame
    for(i in 1:nrow(gdf)){
      gdfa <- data.frame(col_index= as.numeric(unlist(strsplit(gdf[i,]$col_index,","))),
                         value = as.numeric(unlist(strsplit(gdf[i,]$value,","))))
      gdfa <- gdfa[order(gdfa$col_index,decreasing = F),]
      cnta[gdfa$col_index,i] <- gdfa$value
      rm(gdfa,i)
    }
    if(ncol(cnta)>1){
      cnta <- rowMeans(cnta)
    }else{
      cnta <- cnta$V1
    }
    louvain$value <- cnta
  }else{
    louvain <- louvain[0,]
  }
  return(louvain)
}

plot_multigene_grouped_heatmap01 <- function (pldf, genenames = NULL, feature, x = "V1", y = "V2", 
          ycol = "value", xrange = NULL, yrange = NULL, log.transform = F, 
          colScale = c("#fcfbfd", "#9e9ac8", "#3f007d"), minValCol = "#f0f0f0") 
{
  if (is.null(genenames) | length(pldf) != length(genenames)) {
    stop("provided number of gene names and plot data tables do not match\n")
  }
  if (!is.null(xrange)) {
    for (i in 1:length(pldf)) {
      g <- pldf[[i]]
      g <- subset(g, g[, x] >= xrange[1] & g[, x] <= xrange[2])
      pldf[[i]] <- g
      rm(g)
    }
  }
  if (!is.null(yrange)) {
    for (i in 1:length(pldf)) {
      g <- pldf[[i]]
      g <- subset(g, g[, y] >= yrange[1] & g[, y] <= yrange[2])
      pldf[[i]] <- g
      rm(g)
    }
  }
  if (log.transform == T) {
    for (i in 1:length(pldf)) {
      pldf[[i]][, ycol] <- as.numeric(pldf[[i]][, ycol])
      pldf[[i]][, ycol] <- log2(pldf[[i]][, ycol] + 1)
    }
  }
  pl <- c()
  for (i in 1:length(pldf)) {
    funx <- function(x) {
      sum(!is.na(x))
    }
    funy <- function(x) {
      o <- mean(x,na.rm=T)
      ifelse(!is.finite(o),0,o)
    }
    z <- pldf[[i]] %>% dplyr::group_by_at(feature) %>% dplyr::select(dplyr::all_of(ycol)) %>% 
      dplyr::summarise_all(c(funy, length, funx))
    z <- as.data.frame(z)
    names(z) <- c("group", "average", "count", "expnperc")
    z$name <- genenames[i]
    pl <- rbind(pl, z)
    rm(z)
  }
  pl$expnperc <- round(pl$expnperc * 100/pl$count, 3)
  pal <- (grDevices::colorRampPalette(colScale))(length(unique(pl$average)))
  pal[1] <- minValCol
  require(heatmaply)
  z <- reshape2::dcast(pl, group ~ name, value.var = "average")
  rownames(z) <- z$group
  z$group <- NULL
  if (ncol(z) > 1) {
    z <- apply(z, 1, function(x) {
      x <- x * 100/sum(x)
      x <- ifelse(!is.finite(x), 0, x)
      return(x)
    }) %>% t %>% as.data.frame
  }
  heat <- heatmaply::heatmaply(z, colors = pal, Colv = NULL, 
                               show_dendrogram = c(T, F))
  o <- heat[["x"]][["data"]][[1]][["text"]][, 1]
  o <- gsub("<br>.*$", "", gsub("row: ", "", o))
  pl$group <- factor(pl$group, levels = o)
  rows <- min(c(800, ceiling(nrow(pl)/50) * 400))
  columns <- min(c(800, ceiling(ncol(pl)/2) * 200))
  P <- ggplot(pl, aes(x = name, y = group, col = average, size = expnperc, 
                      text = paste0("Gene: ", name, 
                                    "<br>Cell type: ",group, 
                                    "<br>Expression: ", formatC(average, digits = 2,format = "e"), 
                                    "<br>Cell count: ", count, 
                                    "<br>#Cells Expressing: ", round((count*expnperc/100)), 
                                    "<br>%Cell Expressing: ", expnperc, "%<br>"))) + 
    geom_point() + theme(axis.text.x = element_text(color = "black", 
                                                    angle = 90, hjust = 0.5)) + 
    scale_color_gradient2(limits = range(pl$average,na.rm = T), low = "#ffffe5", mid = "#78c679", high = "#004529",space = "Lab") +
    theme_bw() + xlab("") + ylab("")
  P <- P + labs(color = "normalized \nexpression")
  P <- P %>% ggplotly(tooltip = "text")
  P <- P %>% 
    layout(width = columns, height = rows,
           xaxis = list(tickangle=30, tickfont = list(color='crimson')),
           yaxis = list(tickfont = list(color='crimson')))
  plotlist <- list(heat, P)
  return(plotlist)
}

