# SciViewer

Single Cell Interactive Viewer is an R shiny application that allows users to interactively visualize single cell datasets. The application contains two modules - 1) SciViewerIn: The input application that processes single cell data from .h5ad or .seurat format into either a local or remote database. 2) SciViewer: The main application that allows interactive visualization. 

**Installation**
- The application is developed in R 4.0 and uses several dependancies listed below. 
- R Shiny libraries
    * library(shiny)
    * library(shinydashboard)
    * library(shinyFiles)
    * library(shinyWidgets)
    * library(shinyjs)
    * library(dashboardthemes)
    * library(shinycssloaders)
    * library(shinyBS)
- R data handling libraries
    * library(Seurat)
    * library(reticulate)
    * library(data.table)
    * library(Matrix)
    * library(SparseM)
    * library(dplyr)
    * library(preprocessCore)
    * library(DT)
    * library(optparse)
    * library(gprofiler2)
- R data visualization libraries
    * library(highcharter)
    * library(plotly)
    * library(rasterly)
    * library(png)
    * library(heatmaply)
- R database libraries
    * library(RSQLite)
    * library(RMySQL)
    * library(DBI)
- Custom R library
    * library(pddExpn) : Available in this repo

- Set up R locally or on your R-server. 
- Clone the repo and install all above-listed dependencies.
- Run the **SciViewerIn** applications as following
```
runApp(appDir = "/path/to/SciViewerIn/",launch.browser = T)
```
- Run the **SciViewer** applications as following
```
runApp(appDir = "/path/to/SciViewer/",launch.browser = T)
```



