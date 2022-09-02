tabItem(tabName = "compare",
        tags$hr(),
        fluidRow(
          column(width=3,
                 div(style="display: inline-block;vertical-align:center; horizontal-aling:right;",
                     selectizeInput("sel_genex",label="Select gene",choices = NULL,multiple = F)
                 )
          ),
          column(width=3,
                 tags$br(),
                 div(style="display: inline-block;vertical-align:center; horizontal-aling:left;",
                     actionButton("sel_genexgo","Explore across studies!",icon = icon("arrows-rotate"))
                 )
          ),
          column(width=6)
        ),
        tabBox(title = "",width = 14,side='left',
               tabPanel(title="Expression comparison",
                        fluidRow(
                          tags$hr(),
                          tags$br(),
                          div(style='width:800px;overflow-x: auto;height:400px;overflow-y: auto;',
                              DT::dataTableOutput("comparativeExpnTable") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                          )
                        ),
                        fluidRow(
                          div(style='width:800px;overflow-x: auto;height:600px;overflow-y: auto;',
                              uiOutput("comparativeExpnPlot")
                          ),
                          shiny::HTML("X-axis represents average normalized expression per individual. Please note that the presented expression values are log-transformed after averaging!")
                        )
               ),
               tabPanel(title="Gene as Marker",
                        fluidRow(
                          tags$hr(),
                          tags$br(),
                          div(style='width:800px;overflow-x: auto;height:400px;overflow-y: auto;',
                              DT::dataTableOutput("comparativeMarkerTable") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                          )
                        )
               ),
               tabPanel(title="Gene as Biological Marker",
                        fluidRow(
                          tags$hr(),
                          tags$br(),
                          div(style='width:800px;overflow-x: auto;height:400px;overflow-y: auto;',
                              DT::dataTableOutput("comparativeBioMarkerTable") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                          )
                        )
               )
        ),
        tags$br()  
) ## End of compare
