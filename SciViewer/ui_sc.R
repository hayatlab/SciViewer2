tabItem(tabName = "scstudy",
        tags$hr(),
        fluidRow(column(width=3,
                        div(style='width:300px;overflow-x: auto;height:400px;overflow-y: auto;',
                            box(
                              selectInput("sel_celltype",label="select user-provided celltype",choices = NULL,multiple = T),
                              div(style="display: inline-block;vertical-align:top; horizontal-aling:center;",
                                  actionButton("sel_celltypego"," ",icon = icon("arrows-rotate"))
                              ),
                              solidHeader = T, collapsible = F,collapsed = F, width = 12
                            ),
                            box(
                              selectInput("sel_donor",label="select donors",choices = NULL,multiple = T),
                              div(style="display: inline-block;vertical-align:top; horizontal-aling:center;",
                                  actionButton("sel_donorgo"," ",icon = icon("arrows-rotate"))
                              ),
                              solidHeader = T, collapsible = F,collapsed = F, width = 12
                            ),
                            box(
                              selectInput("sel_catvars",label="feature",choices = NULL,multiple = F),
                              selectInput("sel_catvarsval",label="sub-feature",choices = NULL,multiple = T),
                              div(style="display: inline-block;vertical-align:top; horizontal-aling:center;",
                                  actionButton("sel_metadatago"," ",icon = icon("arrows-rotate"))
                              ),
                              solidHeader = T, collapsible = F,collapsed = F, width = 12
                            )
                        ),
        ),
        column(width=6,
               plotlyOutput('umapplot'),
               #-- get plot 
               tags$br()
        ),
        column(width=3,
               fluidRow(column(width=12,
                               htmlOutput("umapplot_hover")),
                        style = "height:120px; background-color: white;"),
               tags$hr(),
               box(
                 selectizeInput("sel_gene",label="Select gene",choices = NULL,multiple = T),
                 actionButton("sel_genego","Visualize!",icon = icon("arrows-rotate")),
                 solidHeader = T, collapsible = F,collapsed = F, width = 12
               )
        )
        ),
        tags$hr(),
        tabBox(title = "",width = 14,side='left',
               tabPanel(title=" "),
               tabPanel(title="Data QC",
                        verticalTabsetPanel(
                          color = "gray",
                          contentWidth = 11,
                          verticalTabPanel(title=h6("Study description", style = 'font-size:13px;color:black;'),
                                           icon=icon("info"),
                                           box_height = "80px",
                                           box_width = "80px",
                                           htmlOutput("sc_study_attribs01")
                          ),
                          verticalTabPanel(title=h6("Cell counts", style = 'font-size:13px;color:black;'),
                                           icon=icon("signal", lib = "glyphicon"),
                                           box_height = "100px",
                                           box_width = "80px",
                                           fluidRow(column(width = 9,
                                                           highchartOutput("sc_qcbarplot",width = "600px",height = "300px") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                           ),
                                           column(width = 3,
                                                  box(
                                                    selectizeInput("qc_attrib01",label="Select attribute",choices = NULL,multiple = F),
                                                    solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                  ),
                                                  tags$br()
                                           )
                                           )
                          ),
                          verticalTabPanel(title=h6("QC plots", style = 'font-size:13px;color:black;'),
                                           box_height = "80px",
                                           box_width = "80px",
                                           icon=icon("eye-open", lib = "glyphicon"),
                                           tags$style(HTML(".select-dropdown { font-size: 11px; }")),
                                           fluidRow(column(width = 9,
                                                           div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;',
                                                               uiOutput("sc_qcscatter"))%>% shinycssloaders::withSpinner(color="#0dc5c1")
                                           ),
                                           column(width = 3,
                                                  box(
                                                    selectizeInput("qc_attrib02",label="Select QC attribute",choices = NULL,multiple = T),
                                                    solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                  )
                                           )
                                           ),
                                           tags$br()
                          )
                        ),
                        tags$br()
               ),
               tabPanel(title="Single-Gene",
                        verticalTabsetPanel(
                          color = "gray",
                          contentWidth = 11,
                          verticalTabPanel(title=h6("Scatter", style = 'font-size:13px;color:black;'),
                                           icon=icon("signal", lib = "glyphicon"),
                                           box_height = "80px",
                                           box_width = "80px",
                                           fluidRow(column(width=9, 
                                                           #-- plot
                                                           plotlyOutput('expnumapplot')
                                           ),
                                           column(width=3,
                                                  box(
                                                    selectInput("sc_gcolscale",label="select color-scale",choices = NULL,multiple = F),
                                                    solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                  )
                                           )
                                           ),
                                           tags$br()
                          ),
                          verticalTabPanel(title=h6("Summary", style = 'font-size:13px;color:black;'),
                                           icon=icon("signal", lib = "glyphicon"),
                                           box_height = "80px",
                                           box_width = "80px",
                                           fluidRow(column(width = 9,
                                                           highchartOutput("sc_gebarplot") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                                           tags$br(),
                                                           shiny::HTML("When expression is summarized per donor, only first attribute from the selected <b> Metadata feature</b> is used for summarizing the data.<br>")
                                           ),
                                           column(width = 3,
                                                  box(
                                                    selectInput("sc_expncatvar",label="Metadata feature",choices = NULL,multiple = T),
                                                    radioButtons("rad_expnchoice", label = shiny::HTML("<b>Cells to use</b>"),
                                                                 choices = list("All" = 1, "Only expressing" = 2), 
                                                                 selected = 1,inline = F),
                                                    radioButtons("rad_expngraphchoice", label = shiny::HTML("<b>Graph type</b>"),
                                                                 choices = list("Barplot" = 1, "Boxplot (outlier omitted)" = 2), 
                                                                 selected = 1,inline = F),
                                                    radioButtons("rad_expnptgrp", label = shiny::HTML("<b>Group by</b>"),
                                                                 choices = list("Cell" = 1, "Donor" = 2), 
                                                                 selected = 1,inline = F),
                                                    actionButton("sel_sgexpngo","",icon = icon("arrows-rotate")),
                                                    solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                  )
                                           )
                                           ),
                                           tags$br()
                          ),
                          verticalTabPanel(title=h6("Counts", style = 'font-size:13px;color:black;'),
                                           icon=icon("signal", lib = "glyphicon"),
                                           box_height = "80px",
                                           box_width = "80px",
                                           fluidRow(column(width = 9,
                                                           highchartOutput("sc_expncountbar") %>% shinycssloaders::withSpinner(color="#0dc5c1") ,
                                                           tags$br(),
                                                           shiny::HTML("The plot summarizes % of total cells expressing the genes along selected <b>Metadata feature</b><br>")
                                           ),
                                           column(width = 3,
                                                  box(
                                                    selectInput("sc_catvar_expncount",label="Metadata feature",choices = NULL,multiple = T),
                                                    solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                  )
                                           )
                                           ),
                                           tags$br()
                          )
                        )
               ),
               tabPanel(title="Multi-Gene",
                        #shiny::HTML("Select multiple genes in the <b><i>Select gene</b></i> section above."),
                        verticalTabsetPanel(
                          color = "gray",
                          contentWidth = 11,
                          verticalTabPanel(title=h6("Features", style = 'font-size:13px;color:black;'),
                                           icon=icon("signal", lib = "glyphicon"),
                                           box_height = "80px",
                                           box_width = "80px",
                                           div(style='width:700px;overflow-x: auto;height:400px;overflow-y: auto;',
                                               uiOutput("sc_multigeneumap")),
                                           shiny::HTML("<b>NOTE:</b> The expression scale is arbritrary and the visualization should be used for qualitative comparison<br>"),
                                           tags$br()
                          ),
                          verticalTabPanel(title=h6("Dotplot", style = 'font-size:13px;color:black;'),
                                           icon=icon("signal", lib = "glyphicon"),
                                           box_height = "80px",
                                           box_width = "80px",
                                           fluidRow(column(width = 9,
                                                           div(style='width:700px;overflow-x: auto;height:400px;overflow-y: auto;',
                                                               uiOutput("sc_expndotplot")%>% shinycssloaders::withSpinner(color="#0dc5c1")),
                                                           shiny::HTML("The dotplot summarize average gene expression for selected genes.
                                                                                 Size of the point represent % of cells that express the gene. While color represents expression scale.<br>"), 
                                                           tags$br()
                                           ),
                                           column(width = 3,
                                                  box(
                                                    selectInput("sc_catvar_expndot",label="Metadata feature",choices = NULL,multiple = F),
                                                    actionButton("sc_catvar_expndotgo","",icon = icon("arrows-rotate")),
                                                    solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                  )
                                           )
                                           )
                          ),
                          verticalTabPanel(title=h6("Correlations", style = 'font-size:13px;color:black;'),
                                           icon=icon("signal", lib = "glyphicon"),
                                           box_height = "80px",
                                           box_width = "80px",
                                           fluidRow(column(width = 9,
                                                           div(style='width:500px;overflow-x: auto;height:400px;overflow-y: auto;',
                                                               uiOutput("sc_corrplot") %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                                                           shiny::HTML("Pearson correlation coefficients are calculated between the genes along the cells expressing the genes.
                                                                     Based on the selection of cell types, the correlation coefficients are updated.<br>
                                                                     <b>Note:</b> Minimum 2 genes are required for calculating correlation coefficients.<br><br><br>")
                                           ),
                                           column(width = 3,
                                                  box(
                                                    selectInput("sc_corrcelltype",label="select cell types",choices = NULL,multiple = T),
                                                    actionButton("sc_corrcelltypego","",icon = icon("arrows-rotate")),
                                                    solidHeader = T, collapsible = F,collapsed = F, width = 12
                                                  )
                                           )
                                           )
                          )
                        )
               ),
               tabPanel(title="Cell Markers",
                        fluidRow(column(width=10,
                                        verticalTabsetPanel(
                                          color = "gray",
                                          contentWidth = 11,
                                          verticalTabPanel(title=h6("Dotplot", style = 'font-size:13px;color:black;'),
                                                           icon=icon("signal", lib = "glyphicon"),
                                                           box_height = "80px",
                                                           box_width = "80px",
                                                           shiny::HTML("<u>Select genes from the above section and <b> hit the refresh button.</b></u><br><br>"),
                                                           div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;',
                                                               uiOutput("sc_markerdotplot")),
                                                           shiny::HTML("<i>blank plot suggests that selected genes do not qualify as marker at FDR 20%</i>")
                                          ),
                                          verticalTabPanel(title=h6("Tables", style = 'font-size:13px;color:black;'),
                                                           icon=icon("signal", lib = "glyphicon"),
                                                           box_height = "80px",
                                                           box_width = "80px",
                                                           shiny::HTML("<u>Make selection using dropdown menu on right.</u><br><br>"),
                                                           fluidRow(column(width=6,
                                                                           shiny::HTML("<b>Marker table</b><br>Select one or more rows from the table below to generate gene expression plots.<br>"),
                                                                           DT::dataTableOutput("marker_tab") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                           ),
                                                           column(width=6,
                                                                  shiny::HTML("<b>Volcano plot</b><br>The volcano plot summarizes top 500 genes sorted by their FDR values<br>"),
                                                                  plotlyOutput("sc_markerVolcano") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                           )#,
                                                           #style = "height:250px; background-color: white;"
                                                           ),
                                                           div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;',
                                                               uiOutput("sc_markerumap"))
                                          ),
                                          verticalTabPanel(title=h6("Enrichments", style = 'font-size:13px;color:black;'),
                                                           icon=icon("signal", lib = "glyphicon"),
                                                           box_height = "80px",
                                                           tags$br(),
                                                           box_width = "80px",
                                                           uiOutput('markerEnrichmentText'),
                                                           tags$br(),
                                                           uiOutput("enrichmarkers_plot")%>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                                           DT::dataTableOutput("enrichmarker_tab"),
                                                           tags$br(),
                                                           shiny::HTML("The genes at selected FDR cutoff (on right) are used for performing enrichment analyses. The app uses
                                                     <a href=\"https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html\",target=\"_blank\">
                                                     gProfiler2 </a>for calculating the enrichments and the displays.</i><br><br>")
                                          ),
                                          tags$br()
                                        )
                        ),
                        column(width = 2,
                               box(
                                 tags$br(),
                                 selectInput("sel_markercelltype",label="select celltype",choices = NULL,multiple = F),
                                 selectInput("sel_markercomp",label="select comparison",choices = NULL,multiple = F),
                                 tags$br(),
                                 sliderInput("sel_markerfdrslider", label = shiny::HTML("<b>FDR cutoff</b>"), min=0,max=0.2, value = 0.05,step = 0.05),
                                 tags$br(),
                                 radioButtons("sel_markerdirchange", label = shiny::HTML("<b>Direction of change</b>"),
                                              choices = list("Over-expressed" = 1, "Under-expressed" = 2), 
                                              selected = 1,inline = F),
                                 tags$hr(),
                                 actionButton("sc_markeropgo",shiny::HTML("calculate<br>enrichments!"),icon = icon("arrows-rotate")),
                                 solidHeader = T, collapsible = F,collapsed = F, width = 12
                               ),
                        )
                        )
               ),
               tabPanel(title="Biological Markers",
                        fluidRow(column(width=10,
                                        verticalTabsetPanel(
                                          color = "gray",
                                          contentWidth = 11,
                                          verticalTabPanel(title=h6("Dotplot", style = 'font-size:13px;color:black;'),
                                                           icon=icon("signal", lib = "glyphicon"),
                                                           box_height = "80px",
                                                           box_width = "80px",
                                                           shiny::HTML("<u>Select genes from the above section and <b> hit the refresh button.</b></u><br><br>"),
                                                           div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;',
                                                               uiOutput("sc_dedotplot")),
                                                           shiny::HTML("<i>blank plot suggests that selected genes do not qualify as biological markers at FDR 20%</i>")
                                          ),
                                          verticalTabPanel(title=h6("Tables", style = 'font-size:13px;color:black;'),
                                                           icon=icon("signal", lib = "glyphicon"),
                                                           box_height = "80px",
                                                           box_width = "80px",
                                                           shiny::HTML("<u>Make selection using dropdown menu on right.</u><br><br>"),
                                                           fluidRow(column(width=6,
                                                                           shiny::HTML("<b>Biological Marker table</b><br>Select one or more rows from the table below to generate gene expression plots.<br>"),
                                                                           DT::dataTableOutput("de_tab") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                           ),
                                                           column(width=6,
                                                                  shiny::HTML("<b>Volcano plot</b><br>The volcano plot summarizes top 500 genes sorted by their FDR values<br>"),
                                                                  plotlyOutput("sc_deVolcano") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                                                           )#,
                                                           #style = "height:250px; background-color: white;"
                                                           ),
                                                           div(style='width:700px;overflow-x: auto;height:450px;overflow-y: auto;',
                                                               uiOutput("sc_deumap"))
                                          ),
                                          verticalTabPanel(title=h6("Enrichments", style = 'font-size:13px;color:black;'),
                                                           icon=icon("signal", lib = "glyphicon"),
                                                           box_height = "80px",
                                                           tags$br(),
                                                           box_width = "80px",
                                                           uiOutput('deEnrichmentText'),
                                                           tags$br(),
                                                           uiOutput("enrichde_plot")%>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                                           DT::dataTableOutput("enrichde_tab"),
                                                           tags$br(),
                                                           shiny::HTML("The genes at selected FDR cutoff (on right) are used for performing enrichment analyses. The app uses
                                                     <a href=\"https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html\",target=\"_blank\">
                                                     gProfiler2 </a>for calculating the enrichments and the displays.</i><br><br>")
                                          ),
                                          tags$br()
                                        )
                        ),
                        column(width = 2,
                               box(
                                 tags$br(),
                                 selectInput("sel_decelltype",label="select celltype",choices = NULL,multiple = F),
                                 selectInput("sel_decomp",label="select comparison",choices = NULL,multiple = F),
                                 tags$br(),
                                 sliderInput("sel_defdrslider", label = shiny::HTML("<b>FDR cutoff</b>"), min=0,max=0.2, value = 0.05,step = 0.05),
                                 tags$br(),
                                 radioButtons("sel_dedirchange", label = shiny::HTML("<b>Direction of change</b>"),
                                              choices = list("Over-expressed" = 1, "Under-expressed" = 2), 
                                              selected = 1,inline = F),
                                 tags$hr(),
                                 actionButton("sc_deopgo",shiny::HTML("calculate<br>enrichments!"),icon = icon("arrows-rotate")),
                                 solidHeader = T, collapsible = F,collapsed = F, width = 12
                               ),
                        )
                        )
               )## add here
        )
)
