library(shiny)
library(DT)
library(shinyjs)
library(shinyalert)

# UI
######################################################################################
######################################################################################
landscapeTab <- {
  shiny::tabPanel(
    "Landscape",
    shiny::br(),
    shiny::fluidRow(
      column(
        width = 6,
        shiny::selectInput(
          inputId = "landscape",
          label = "Landscape structure (field boundaries)",
          choices = list(
            "Landscape 1" = 1,
            "Landscape 2" = 2,
            "Landscape 3" = 3,
            "Landscape 4" = 4,
            "Landscape 5" = 5
          ),
          selected = 1,
        )
      ),
      shinyBS::bsTooltip("landscape",title="Landscape Shapefile", placement = "top", trigger="hover"),
      column(
        width = 6,
        shiny::selectInput(
          inputId = "aggregLevel",
          label = "Spatial aggregation of croptypes",
          choices = list(
            "Highly fragmented" = "low",
            "Balanced" = "medium",
            "Highly aggregated" = "high"
          ),
          selected = "low",
        )
      )
    ),
    shinyBS::bsTooltip("aggregLevel",title="Level of spatial aggregation of the landscape", placement = "top", trigger="hover"),
    hr(),
    shiny::fluidRow(
      tags$div(lang="en",
                 h3("Croptypes"),
                 editableDTUI(id = "croptypes")
      ),
      hr(),
      column(
        width = 4,
        IntegerInput(
          inputId = "rotationPeriod",
          label = "Rotation period (years)",
          value = 0,
          max = VALUEMAX
        )
      ),
      shinyBS::bsTooltip("rotationPeriod",title=ROTATION_PERIOD, placement = "bottom", trigger="hover"),
      column(
        width = 8,
        align = "left",
        htmlOutput("rotationText")
        )
    ),
    hr(),
    shiny::fluidRow(
      column(
        width = 4,
        IntegerInput(
          inputId = "nYear",
          label = "Simulation duration (years)",
          value = 10,
          max = 50
        )
      ),
      shinyBS::bsTooltip("nYear",title="Number of cropping seasons (e.g. years)", placement = "bottom", trigger="hover"),
      column(
        width = 4,
        IntegerInput(
          inputId = "nTSpY",
          label = "Time steps per year (days)",
          value = 120,
          max = 365
        )
      ),
      shinyBS::bsTooltip("nTSpY",title="Number of time steps per cropping season (e.g. days)", placement = "bottom", trigger="hover"),
      column(
        width = 4,
        IntegerInput(
          inputId = "seed",
          label = "Seed (RNG)",
          value = 1,
          max = 99999
        )
      ),
      shinyBS::bsTooltip("seed",title="Seed value for Random Number Generator", placement = "bottom", trigger="hover"),
    )
  )
}
######################################################################################
cultivarTab <- {
  shiny::tabPanel(
    "Cultivars and Genes",
    h3("Genes"),
    editableDTUI(id = "genes"),
    h3("Cultivars"),
    editableDTUI(id = "cultivars"),
    h3("Cultivars and Genes"),
    editableDTUI(id = "cultivarsgenes"),
  )
}

######################################################################################
pathogenTab <- {
  shiny::tabPanel(
    "Pathogen",
    shiny::div(
      shiny::selectInput(
        inputId = "defaultPathogen",
        label = "Default Pathogen",
        choices = list(
          "Rust" = "Rust"
        ),
        width = "25%"
      ),
      align = "center"
    ),
    shiny::fluidRow(
      shiny::numericInput(
        inputId = "inoculum",
        label = "Initial probability of infection for the first susceptible host",
        value = 0.0001,
        min = 0.0,
        max = 1.0,
        step = 0.0001
      ),
      shinyBS::bsTooltip("inoculum",title=INOCULUM, placement = "bottom", trigger="hover"),
    ),
    shiny::fluidRow(
      column(
        width = 4,
        shiny::numericInput(
          inputId = "patho_survival_prob",
          label = "Off-season survival probability",
          value = 0.0001,
          min = 0.0001,
          max = 1.0,
          step = 0.0001
        ),
        shinyBS::bsTooltip("patho_survival_prob",title=SURVIVAL_PROB, placement = "right", trigger="hover"),
        # shiny::numericInput(
        #   inputId = "patho_repro_sex_prob",
        #   label = "Prob. for an infectious host to reporduce via sex rather than clonal",
        #   value = 0,
        #   min = 0.0,
        #   max = 1.0,
        #   step = 0.10
        # ),
        shiny::numericInput(
          inputId = "patho_infection_rate",
          label = "Infection rate",
          value = 0.4,
          min = 0.0,
          max = 2.0,
          step = 0.1
        ),
        shinyBS::bsTooltip("patho_infection_rate",title=INFECTION_RATE, placement = "right", trigger="hover"),
        shiny::numericInput(
          inputId = "patho_propagule_prod_rate",
          label = "Reproduction rate (clonal)",
          value = 3.125,
          min = 0.0,
          step = 0.4
        ),
        shinyBS::bsTooltip("patho_propagule_prod_rate",title=PROPAGULE_PROD_RATE, placement = "right", trigger="hover")
      ),
      column(
        width = 4,
        shiny::numericInput(
          inputId = "patho_latent_period_exp",
          label = "Latent period duration",
          value = 10,
          min = 0,
          max = 100,
          step = 1
        ),
        shinyBS::bsTooltip("patho_latent_period_exp",title=LATENT_PERIOD_EXP, placement = "left", trigger="hover"),
        shiny::numericInput(
          inputId = "patho_latent_period_var",
          label = "Variance of the latent period duration",
          value = 9,
          min = 0.0,
          max = 100,
          step = 1
        ),
        shinyBS::bsTooltip("patho_latent_period_var",title=LATENT_PERIOD_VAR, placement = "left", trigger="hover"),
        shiny::numericInput(
          inputId = "patho_infectious_period_exp",
          label = "Infectious period duration",
          value = 24,
          min = 0,
          max = 365,
          step = 1
        ),
        shinyBS::bsTooltip("patho_infectious_period_exp",title=INFECTIOUS_PERIOD_EXP, placement = "left", trigger="hover"),
        shiny::numericInput(
          inputId = "patho_infectious_period_var",
          label = "Variance of the infectious period duration",
          value = 105,
          min = 0,
          step = 1
        ),
        shinyBS::bsTooltip("patho_infectious_period_var",title=INFECTIOUS_PERIOD_VAR, placement = "left", trigger="hover")
      ),
      column(
        width = 4,
        shiny::numericInput(
          inputId = "patho_sigmoid_kappa",
          label = "Contamination function: Kappa",
          value = 5.333,
          min = 0.0001,
          max = 10,
          step = 0.01
        ),
        shinyBS::bsTooltip("patho_sigmoid_kappa",title=SIGMOID_KAPPA, placement = "left", trigger="hover"),
        shiny::numericInput(
          inputId = "patho_sigmoid_sigma",
          label = "Contamination function: Sigma",
          value = 3,
          min = 0.0,
          max = 100,
          step = 1
        ),
        shinyBS::bsTooltip("patho_sigmoid_sigma",title=SIGMOID_SIGMA, placement = "left", trigger="hover")
        #,
        # shiny::numericInput(
        #   inputId = "patho_sigmoid_plateau",
        #   label = "Plateau parameter of the sigmoid contamination function",
        #   value = 1,
        #   min = 0,
        #   max = 100,
        #   step = 1
        # )
      )
    )
  )
}


######################################################################################
inputUi <- {
  shiny::sidebarPanel(
    id = "inputpanel",
    shiny::h3("Input", align = "center"),
    shiny::div(
      shiny::selectInput(
        inputId = "demo",
        label = "Default Strategies",
        choices = list(
          "Mosaic" = "MO",
          "Mixture" = "MI",
          "Rotation" = "RO",
          "Pyramiding" = "PY"
        ),
        width = "25%"
      ),
      align = "center"
    ),
    shinyBS::bsTooltip("demo",title="Load existing parameters or click Advanced Mode", placement = "top", trigger="hover"),
    shiny::tabsetPanel(id = "inputtabpanel", landscapeTab, cultivarTab, pathogenTab),
    width = 12,
    align = "center"
  )
}
######################################################################################
outputUi <- {
  shiny::mainPanel(
    shiny::h3("Output"),
    shiny::plotOutput(outputId = "landscapeimg", dblclick = "plot_landscapeimg"),
    shiny::uiOutput(outputId = "video"),
    width = 12,
    align = "center",
    id = "outputpanel"
  )
}
######################################################################################
ui <- {
  shiny::fluidPage(
    tags$html(lang="en"),
    title = "Landsepi Demo",
    shinyjs::useShinyjs(),
    useShinyalert(),
    tags$head(
      tags$link(href = "style.css", rel = "stylesheet"),
      tags$script(src = "script.js"),
      tags$meta(name = "description", content = "A stochastic, spatially-explicit, demo-genetic model 
                simulating the spread and evolution of a plant pathogen in a heterogeneous landscape 
                to assess resistance deployment strategies. It is based on a spatial geometry for describing 
                the landscape and allocation of different cultivars, a dispersal kernel for the 
                dissemination of the pathogen, and a SEIR ('Susceptible-Exposed-Infectious-Removed’) 
                structure with a discrete time step. It provides a useful tool to assess the performance 
                of a wide range of deployment options with respect to their epidemiological, 
                evolutionary and economic outcomes."),
      tags$meta(name = "title", content = "Landsepi: Landscape Epidemiology and Evolution"),
      tags$meta(name = "author", content = "Jean-François Rey"),
      tags$meta(name = "keywords", content = "R, Shiny, INRAE, landscape, epidemiology, strategies, deployment, hosts, pathogens, stochastic, spatio-temporal, demonstration")
    ),
    fluidRow(
      titlePanel(div(img(src="landsepi-logo.png", width="60" ), "Landsepi : Landscape Epidemiology and Evolution")),
      actionButton("About", "About"),
      actionButton("Mode", "Advanced Mode On/Off", icon = icon("exchange-alt")),
      shinyBS::bsTooltip("Mode",title="Edit all input parameters", placement = "top", trigger="hover"),
      align="center"),
    shiny::br(),
    fluidRow(
      actionButton("showInputside", label = "", icon = icon("wpforms")),
      shinyBS::bsTooltip("showInputside",title="Only show input parameters", placement = "bottom", trigger="hover"),
      actionButton("showBothside", label = "", icon = icon("columns")),
      shinyBS::bsTooltip("showBothside",title="Show both,  input parameters and output", placement = "bottom", trigger="hover"),
      actionButton("showOutputside", label = "", icon = icon("chart-line") ),
      shinyBS::bsTooltip("showOutputside",title="Only show output", placement = "bottom", trigger="hover"),
      align="center"
    ),
    shiny::br(),
    shiny::fluidRow(
      column(6, id = "inputside", inputUi),
      column(6, id = "outputside", outputUi)
    ),
    #shiny::sidebarLayout(inputUi, outputUi),
    shiny::fluidRow(shiny::div(
      shiny::actionButton(inputId = "generateLandscape", label = "Generate the landscape"),
      align = "center"
    )),
    shinyBS::bsTooltip("generateLandscape",title=GENERATE_LANDSCAPE, placement = "top", trigger="hover"),
    shiny::br(),
    shiny::fluidRow(shiny::div(
      shiny::actionButton(inputId = "runSimulation", label = "Run simulation"),
      shiny::actionButton(inputId = "stopSimulation", label = "Stop simulation"),
      align = "center"
    )),
    shinyBS::bsTooltip("runSimulation",title=RUN_SIMULATION, placement = "top", trigger="hover"),
    shinyBS::bsTooltip("stopSimulation",title=STOP_SIMULATION, placement = "top", trigger="hover"),
    shiny::br(),
    shiny::fluidRow(shiny::div(
      shiny::downloadButton(outputId = "export", label = "Export simulation"),
      align = "center"
    )),
    shinyBS::bsTooltip("export",title=EXPORT_SIMULATION, placement = "top", trigger="hover"),
    shiny::br()
  )
}