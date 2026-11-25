install.packages("shiny", quietly = TRUE)
install.packages("DBI", quietly = TRUE)
install.packages("RMariaDB", quietly = TRUE)
install.packages("ggplot2", quietly = TRUE)
install.packages("tidyverse", quietly = TRUE)

library(shiny)
library(DBI)
library(RMariaDB)
library(ggplot2)
library(tidyverse)


# There are 3 visualisations that they require.
#(1) They want the ability to select a particular knockout mouse and visualise
# the statistical scores of all phenotypes tested. The plot should intuitively
# show the phenotypes that are significantly affected by the gene knockout.
# (2) The collaborator also wishes to visualise the
#statistical scores of all knockout mice for a selected phenotype. 
#(3) Lastly, they would like
#to visualise clusters of genes with similar phenotype scores. 



# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("Placeholder"),

    # For (0) - they want to select a particular knockout mouse 
    selectInput( "knockout_mouse",
      "Select particular knockout mouse:",
      choices = NULL,
      multiple = FALSE
      ),
    plotOutput("plot_1"),
    verbatimTextOutput("debug"),
     selectInput( "phenotype",
      "Select particular phenotype:",
      choices = NULL,
      multiple = FALSE
      ),
    plotOutput("plot_2") 
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

output$debug <- renderPrint({
  cat("Selected:", input$knockout_mouse, "\n")
  d <- df_phynotype_statistical_scores_plot_1()
  cat("Rows:", nrow(d), "\n")
  str(d)
})

  df_phynotype_statistical_scores <- data.frame(
  phenotype = rep(c("Phenotype1", "Phenotype2", "Phenotype3"), each = 4),
  knockout  = rep(c("Mouse1", "Mouse2", "Mouse3", "Mouse4"), times = 3),
  p_value   = c(0.01, 0.30, 0.04, 0.80,
                0.05, 0.001, 0.20, 0.60,
                0.001, 0.02, 0.50, 0.03)
)
  
  updateSelectInput(session, "knockout_mouse", choices = unique( df_phynotype_statistical_scores$knockout ) )
  updateSelectInput(session, "phenotype", choices = unique( df_phynotype_statistical_scores$phenotype ) )
  
  df_phynotype_statistical_scores_plot_1 <- reactive (df_phynotype_statistical_scores %>% filter(. , knockout == input$knockout_mouse) )
  # df_phynotype_statistical_scores <- Here we should query the database to get a dataframe of p value + phenotype + knockout
    output$plot_1 <- renderPlot({
      # p-value - numerical continious. visualising over the phenotypes and gene knockout. Violin plot
      ggplot() + geom_violin(data = df_phynotype_statistical_scores_plot_1(), aes(
        x = phenotype ,
        y = p_value
      ))
    })
    df_plot_2 <- reactive({
    req(input$phenotype)
    df_phynotype_statistical_scores %>%
      filter(phenotype == input$phenotype) %>%
      mutate(neg_log_p = -log10(p_value),
             significant = p_value < 0.05)
  })

  output$plot_2 <- renderPlot({
    ggplot(df_plot_2(), aes(x = knockout, y = neg_log_p, fill = significant)) +
      geom_col() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
      labs(x = "Knockout mouse", y = "-log10(p-value)", fill = "p < 0.05",
           title = paste("Phenotype:", input$phenotype)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
}

# Run the application 
runApp(shinyApp(ui = ui, server = server), launch.browser = TRUE)
