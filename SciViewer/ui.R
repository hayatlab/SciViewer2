source("vars.R",local = TRUE,encoding = "UTF-8")
#source("formats.R",local = TRUE,encoding = "UTF-8")
options(shiny.maxRequestSize = 100*1024^2)

# javascript code to collapse box
jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"

## New theme with shiny dashboard
shinyUI(
  tagList(
    tags$style(HTML(".selectize-input { font-size: 13px; }")),
    tags$style(HTML(".select-dropdown { font-size: 13px; }")),
    tags$style(HTML(".main-sidebar { font-size: 13px; }")),
    tags$head(HTML('
                 <!-- Global site tag (gtag.js) - Google Analytics -->
                 <script async src="https://www.googletagmanager.com/gtag/js?id=UA-18253988-3"></script>
                 <script>
                 window.dataLayer = window.dataLayer || [];
                 function gtag(){dataLayer.push(arguments);}
                 gtag("js", new Date());
                 
                 gtag("config", "UA-18253988-3");
                 </script>')),
    
    dashboardPage(
      
      dashboardHeader(
        titleWidth = "100%",
        title = shiny::HTML(
          "<div style = 'background-color:#ffffff;color: #6a51a3;font-size:30px;font-weight:bold; vertical-align:middle'>
           <img src = 'logo.png' align = 'left'  height = '55px' width = '200px'>
            Single cell Interactive Viewer 
           </div>")
      ),

      dashboardSidebar(collapsed = F,
                       width=100,
        sidebarMenu(id = 'AppTab',
                    tags$br(),
                    menuItem(text = shiny::HTML("<br>Browse<br>"),tabName = "scsel"),
                    menuItem(text = shiny::HTML("<br>Visualize<br>"),tabName = "scstudy"),
                    menuItem(text = shiny::HTML("<br>Compare<br>"),tabName = "compare"),
                    menuItem(text = shiny::HTML("<br>Help<br>"),tabName = "help"),
                    tags$br(),
                    tags$br(),
                    tags$br(),
                    tags$br(),
                    tags$br()
        )
      ),
      dashboardBody(
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
          }")),
        
        shinyDashboardThemes(theme = "poor_mans_flatly"),#"purple_gradient" poor_mans_flatly
        tabItems(          
          
          source(file = 'ui_browse.R',local = T,encoding = "UTF-8")$value,
          
          source(file = 'ui_sc.R',local = T,encoding = "UTF-8")$value,
          
          source(file = 'ui_compare.R',local = T,encoding = "UTF-8")$value,
          
          source(file = 'ui_help.R',local = T,encoding = "UTF-8")$value
          
        ) # end of tabItems
      ) # end of dashboard body
    ),# end of dashboard
    
    tags$footer(align="center",
                tags$hr(),
                shiny::HTML(
                  "<div vertical-align:middle'>
           <img src = 'PDD_logo.png' align = 'left'  height = '100px' width = '130px'>
           </div>"),
                tags$p("Copyright (c) 2020",color="black"), 
                tags$a(" Bayer US LLC.; ", href = "http://www.bayer.us/",target="_blank"), 
                tags$a(" PDD labs, Strategic collaboration - Speciality lung diseases", href = "https://www.brighamhealth.org/home", target="_blank"), 
                tags$p("Version 1.0"),
                #tags$p("App logo created using www.designevo.com"),
                style = "bottom:0;width:100%;color: #6a51a3;padding: 10px;
                background-color: #ffffff;z-index: 1000;"
    )
  )  #End of taglist
) #End




