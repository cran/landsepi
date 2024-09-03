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
      shinyBS::bsTooltip("landscape", title = "Landscape Shapefile", placement = "top", trigger = "hover"),
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
    shinyBS::bsTooltip("aggregLevel", title = "Level of spatial aggregation of the landscape", placement = "top", trigger = "hover"),
    hr(),
    shiny::fluidRow(
         shiny::checkboxInput("rotationcheck", "Active Rotation strategy", value = FALSE),
         shinyBS::bsTooltip("rotationcheck", title = "Active Rotation strategy", placement = "top", trigger = "hover"),
        column(
            width = 2,
            IntegerInput(
                inputId = "rotationPeriod",
                label = "Rotation period (years)",
                value = 0,
                max = VALUEMAX
            )
        ),
        shinyBS::bsTooltip("rotationPeriod", title = ROTATION_PERIOD, placement = "bottom", trigger = "hover"),
        column(
            width = 10,
            align = "left",
            htmlOutput("rotationText"),
        )
    ),
    #hr(),
    shiny::fluidRow(
      shiny::div(
        #style="margin:auto;display: inline-block;vertical-align:middle;text-align:left",
        lang = "en",
        h3("Croptypes"),
        editableDTUI(id = "croptypes")
      )),
    shiny::fluidRow(
      shiny::div(
        style="margin:auto;display: inline-block;vertical-align:middle;text-align:left",
        shiny::actionButton(inputId = "addcrop",label = "Add Croptype"),
        shinyBS::bsTooltip("addcrop", title = "Add a croptype type", placement = "right", trigger = "hover")
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
      shinyBS::bsTooltip("nYear", title = "Number of cropping seasons (e.g. years)", placement = "bottom", trigger = "hover"),
      column(
        width = 4,
        IntegerInput(
          inputId = "nTSpY",
          label = "Time steps per year (days)",
          value = 120,
          max = 365
        )
      ),
      shinyBS::bsTooltip("nTSpY", title = "Number of time steps per cropping season (e.g. days)", placement = "bottom", trigger = "hover"),
      column(
        width = 4,
        IntegerInput(
          inputId = "seed",
          label = "Seed (RNG)",
          value = 1,
          max = 99999
        )
      ),
      shinyBS::bsTooltip("seed", title = "Seed value for Random Number Generator", placement = "bottom", trigger = "hover"),
    )
  )
}
######################################################################################
cultivarTab <- {
  shiny::tabPanel(
    "Cultivars and Genes",
    h3("Genes"),
    editableDTUI(id = "genes"),
    shiny::fluidRow(
      column(12,
    shiny::div(
      style="margin:auto;display: inline-block;vertical-align:middle;text-align:left",
      shiny::selectInput(
        inputId = "listgenes",
        label = "",
        choices = listGenes,
        selected = 1,
        width = "120px",
      ),
      shinyBS::bsTooltip("listgenes", title = "Select a gene type", placement = "top", trigger = "hover"),
    ),
    shiny::div(
      style="margin:auto;display: inline-block;vertical-align:middle;text-align:left",
      shiny::actionButton(inputId = "addgene",label = "Add Gene"),
      shinyBS::bsTooltip("addgene", title = "Add selected Gene type", placement = "right", trigger = "hover")
    )
    )),
    # h3("Cultivars"),
    editableDTUI(id = "cultivars"),
    h3("Cultivars and Genes"),
    editableDTUI(id = "cultivarsgenes"),
    shiny::fluidRow(
      column(12,
             shiny::div(
               style="margin:auto;display: inline-block;vertical-align:middle;text-align:left",
               shiny::selectInput(
                 inputId = "listcultivarstype",
                 label = "",
                 choices = listCultivarsType,
                 selected = 1,
                 width = "150px",
               ),
               shinyBS::bsTooltip("listcultivarstype", title = "Select a cultivars type", placement = "top", trigger = "hover"),
             ),
             shiny::div(
               style="margin:auto;display: inline-block;vertical-align:middle;text-align:left",
               shiny::actionButton(inputId = "addcultivar",label = "Add a cultivar"),
               shinyBS::bsTooltip("addcultivar", title = "Add selected cultivar type", placement = "right", trigger = "hover")
             )
      ))
  )
}

######################################################################################
pathogenTab <- {
  shiny::tabPanel(
    "Pathogen",
    # shiny::div(
    #   shiny::selectInput(
    #     inputId = "defaultPathogen",
    #     label = "Default Pathogen",
    #     choices = list(
    #       "Rust of wheat" = "Rust",
    #       "Grapevine downy mildew" = "Mildew",
    #       "Banana black sigatoka" = "Sigatoka",
    #       "No Pathogen" = "No Pathogen"
    #     ),
    #     width = "25%"
    #   ),
    #   align = "center",
    #   shinyBS::bsTooltip("defaultPathogen", title = "Select a Pathogen <br><b>!!! That will reload all simulation parameters in normal Mode !!!</b>", placement = "top", trigger = "hover"),
    # ),
    shiny::fluidRow(
      column(6,
             shiny::numericInput(
               inputId = "inoculum",
               label = "Initial probability of infection for susceptible cultivar",
               value = 0.0005,
               min = 0.0,
               max = 1.0,
               step = 0.0001
             ),
             shinyBS::bsTooltip("inoculum", title = INOCULUM, placement = "bottom", trigger = "hover"),
             shiny::h5(id="inoculum_message", "(Local initial probability of infection for susceptible cultivar should be upper than 0.01)")
      ),
      shiny::column(
        width = 6,
        shiny::radioButtons("inoculumstrategy", "Inoculum localisation", c("Global"=1,"Local"=0)),
        shinyBS::bsTooltip("inoculumstrategy", title = "Global : all susceptibles cultivars fields are infected | Local : a random susceptible field is infected", placement = "top", trigger = "hover"),
      ),
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
        shinyBS::bsTooltip("patho_survival_prob", title = SURVIVAL_PROB, placement = "right", trigger = "hover"),
        # shiny::numericInput(
        #   inputId = "patho_repro_sex_prob",
        #   label = "Prob. for an infectious host to reproduce via sex rather than clonal",
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
        shinyBS::bsTooltip("patho_infection_rate", title = INFECTION_RATE, placement = "right", trigger = "hover"),
        shiny::numericInput(
          inputId = "patho_propagule_prod_rate",
          label = "Reproduction rate (clonal)",
          value = 3.125,
          min = 0.0,
          step = 0.4
        ),
        shinyBS::bsTooltip("patho_propagule_prod_rate", title = PROPAGULE_PROD_RATE, placement = "right", trigger = "hover")
      ),
      column(
        width = 4,
        shiny::numericInput(
          inputId = "patho_latent_period_mean",
          label = "Latent period duration",
          value = 10,
          min = 0,
          max = 100,
          step = 1
        ),
        shinyBS::bsTooltip("patho_latent_period_mean", title =  LATENT_PERIOD_MEAN, placement = "left", trigger = "hover"),
        # shiny::numericInput(
        #   inputId = "patho_latent_period_var",
        #   label = "Variance of the latent period duration",
        #   value = 9,
        #   min = 0.0,
        #   max = 100,
        #   step = 1
        # ),
        # shinyBS::bsTooltip("patho_latent_period_var", title = LATENT_PERIOD_VAR, placement = "left", trigger = "hover"),
        shiny::numericInput(
          inputId = "patho_infectious_period_mean",
          label = "Infectious period duration",
          value = 24,
          min = 0,
          max = 365,
          step = 1
        ),
        shinyBS::bsTooltip("patho_infectious_period_mean", title = INFECTIOUS_PERIOD_MEAN, placement = "left", trigger = "hover"),
        # shiny::numericInput(
        #   inputId = "patho_infectious_period_var",
        #   label = "Variance of the infectious period duration",
        #   value = 105,
        #   min = 0,
        #   step = 1
        # ),
        # shinyBS::bsTooltip("patho_infectious_period_var", title = INFECTIOUS_PERIOD_VAR, placement = "left", trigger = "hover")
      ),
      # column(
      #   width = 4,
      #   shiny::numericInput(
      #     inputId = "patho_sigmoid_kappa",
      #     label = "Contamination function: Kappa",
      #     value = 5.333,
      #     min = 0.0001,
      #     max = 10,
      #     step = 0.01
      #   ),
      #   shinyBS::bsTooltip("patho_sigmoid_kappa", title = SIGMOID_KAPPA, placement = "left", trigger = "hover"),
      #   shiny::numericInput(
      #     inputId = "patho_sigmoid_sigma",
      #     label = "Contamination function: Sigma",
      #     value = 3,
      #     min = 0.0,
      #     max = 100,
      #     step = 1
      #   ),
      #   shinyBS::bsTooltip("patho_sigmoid_sigma", title = SIGMOID_SIGMA, placement = "left", trigger = "hover"),
        # ,
        # shiny::numericInput(
        #   inputId = "patho_sigmoid_plateau",
        #   label = "Plateau parameter of the sigmoid contamination function",
        #   value = 1,
        #   min = 0,
        #   max = 100,
        #   step = 1
        # )
      # )
    ),
    # shiny::fluidRow(
    #   shiny::column(
    #     width = 4,
    #     shiny::checkboxInput("patho_repro_sex_active", "Active Sexual Reproduction", value = FALSE),
    #     shinyBS::bsTooltip("patho_repro_sex_active", title = "Active sexual reproduction at the end of the season", placement = "top", trigger = "hover"),
    #   ),
    #   column(
    #     width = 4,
    #     shiny::numericInput(
    #       inputId = "patho_sex_propagule_release_mean",
    #       label = "Average number of seasons before release of sexual propagules",
    #       value = 1.0,
    #       min = 1.0,
    #       step = 1.0
    #     ),
    #     shinyBS::bsTooltip("patho_sex_propagule_release_mean", title = SEX_PROPAGULE_RELEASE_MEAN, placement = "left", trigger = "hover")
    #   ),
    #   column(
    #     width = 4,
    #     shiny::numericInput(
    #       inputId = "patho_sex_propagule_viability_limit",
    #       label = "Sexual propagule viability limit: nb of seasons (years)",
    #       value = 5,
    #       min = 1,
    #       step = 1
    #     ),
    #     shinyBS::bsTooltip("patho_sex_propagule_viability_limit", title = SEX_PROPAGULE_VIABILITY_LIMIT, placement = "left", trigger = "hover"),
    #   )
    # )
  )
}

######################################################################################
treatmentTab <- {
  shiny::tabPanel(
    "Treatment",
    shiny::br(),
    shiny::fluidRow(
      column(12,
        align = "center",
        shiny::checkboxInput("treatment_active", "Active Cultivars Treatment", value = FALSE),
        shinyBS::bsTooltip("treatment_active", title = "Active cultivars treatment", placement = "top", trigger = "hover"),
      )
    ),
    shiny::fluidRow(
      column(
        width = 4,
        align = "center",
      #   shiny::numericInput(
      #     inputId = "treatment_day_start",
      #     label = "First Day of treatment",
      #     value = 1,
      #     min = 1,
      #     step = 1
      #   ),
      #   shinyBS::bsTooltip("treatment_day_start", title = "First day of treatment", placement = "top", trigger = "hover"),
        shiny::numericInput(
          inputId = "treatment_days_interval",
          label = "Days (step) between treatment",
          value = 14,
          min = 1,
          step = 1
        ),
        shinyBS::bsTooltip("treatment_days_interval", title = "Step between treatment", placement = "top", trigger = "hover"),
      ),
      column(
        width = 4,
        align = "center",
      selectInput("treatment_cultivars_select", "Cultivars to Treat", c(NULL),
                  multiple = TRUE,
                  selectize = TRUE
      ),
      shinyBS::bsTooltip("treatment_cultivars_select", title = "Select cultivars to treat", placement = "top", trigger = "hover")#,
      # shiny::numericInput(
      #     inputId = "treatment_cost",
      #     label = "Treatment cost",
      #     value = 0.0,
      #     min = 0.0,
      #     max = 1000,
      #     step = 0.1
      # ),
      # shinyBS::bsTooltip("treatment_cost", title = "cost of a single treatments (monetary units/ha)", placement = "top", trigger = "hover"),
      ),
      column(
        width = 4,
        align = "center",
      #   shiny::numericInput(
      #     inputId = "treatment_degradation_rate",
      #     label = "Treatment degradation rate",
      #     value = 0.1,
      #     min = 0.01,
      #     step = 0.1
      #   ),
        shinyBS::bsTooltip("treatment_degradation_rate", title = "Degradation per time step of treatment concentration. 0.10 for protectant fungicides, 0.07 for locally systemic fungicides, and 0.06 to 0.05 for systemic fungicides.", placement = "top", trigger = "hover"),
        shiny::numericInput(
          inputId = "treatment_efficiency",
          label = "Treatment efficacy",
          value = 0.0,
          min = 0.0,
          max = 1.0,
          step = 0.1
        ),
        shinyBS::bsTooltip("treatment_efficiency", title = "Maximal efficiency of chemical treatments (i.e. fractional reduction of pathogen infection at the application date) => (0, null efficacy, 1, total efficacy)", placement = "top", trigger = "hover"),
      )
      )
  )
}

######################################################################################
inputUi <- {
  shiny::sidebarPanel(
    id = "inputpanel",
    shiny::h3("Input", align = "center"),
    shiny::fluidRow(
      column(12,
      shiny::div(
        style="margin:auto;display: inline-block;vertical-align:middle;text-align:center",
        shiny::selectInput(
          inputId = "defaultPathogen",
          label = "Pathosystem",
          choices = list(
            "Rust of wheat" = "rust",
            "Grapevine downy mildew" = "mildew",
            "Banana black sigatoka" = "sigatoka"
            #"No Pathogen" = "no Pathogen"
          ),
          width = "200px"
        ),
        #align = "center",
        shinyBS::bsTooltip("defaultPathogen", title = "Select a Pathosystem <br>", placement = "top", trigger = "hover"),
      ),
      shiny::div(
          style="margin:auto;display: inline-block;vertical-align:middle;text-align:center",
          shiny::selectInput(
              inputId = "demo",
              label = "Strategy",
              choices = list(
                  "Mosaic" = "MO",
                  "Mixture" = "MI",
                  "Rotation" = "RO",
                  "Pyramiding" = "PY"
              ),
              width = "130px"
          ),
          shinyBS::bsTooltip("demo", title = "Select a strategy", placement = "top", trigger = "hover"),
      ),
      shiny::div(
        style="margin:auto;display: inline-block;vertical-align:middle;text-align:center",
        actionButton("load", "Load", icon = icon("rotate-right", lib = "font-awesome")),
        shinyBS::bsTooltip("load", title = "Load the strategy and the pathosystem  | Erase all modifications", placement = "top", trigger = "hover"),
      ),
      align="center"
      )
    ),
    shiny::tabsetPanel(id = "inputtabpanel", cultivarTab, landscapeTab, treatmentTab, pathogenTab ),
    width = 12,
    align = "center"
  )
}
######################################################################################
## OUTPUT UI

landscapeTab <- shiny::tabPanel(
  title = "Landscape",
  value = "landscapeTab",
  shiny::br(),
  shiny::plotOutput(outputId = "landscapeimg", dblclick = "plot_landscapeimg")
)

## Video Tab

videoTab <- shiny::tabPanel(
  title = "Video",
  value = "videoTab",
  shiny::br(),
  shiny::uiOutput(outputId = "video")
  )

## Table and graphics TAB

tableTab <- shiny::tabPanel(
  title = "Tables and Graphics",
  value = "tableTab",
  shiny::br(),
  fluidRow(
    column(
      width=6,
      shiny::h4("Epidemiological control (disease severity)"),
      DT::DTOutput("audpctable")
    ),
    column(
      width=6,
      shiny::h4("Evolutionary control (resistance durability)"),
      DT::DTOutput("durabilitytable")
    )
  ),
  fluidRow(
      column(
          width=10,
          br(),
          shiny::h4("Frequency of pathogen genotypes"),
          shiny::imageOutput(outputId = "freqpathogenotypes", inline = FALSE, width = "50%")
  ))
)

## OUTPUT
outputUi <- {
  shiny::mainPanel(
    shiny::h3("Output"),
    # shiny::plotOutput(outputId = "landscapeimg", dblclick = "plot_landscapeimg"),
    # shiny::uiOutput(outputId = "video"),
    shiny::tabsetPanel(id = "outputtabpanel", landscapeTab, videoTab, tableTab ),
    width = 12,
    align = "center",
    id = "outputpanel"
  )
}
######################################################################################
ui <- {
  shiny::fluidPage(
    tags$html(lang = "en"),
    title = "Landsepi Demo",
    shinyjs::useShinyjs(),
    useShinyalert(force=TRUE),
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
      titlePanel(div(img(src = "landsepi-logo.png", width = "60"), "Landsepi : Landscape Epidemiology and Evolution")),
      actionButton("About", "About"),
      shinyBS::bsTooltip("About", title = "About Landsepi", placement = "top", trigger = "hover"),
      #actionButton("Mode", "Advanced Mode On/Off", icon = icon("arrow-right-arrow-left", lib = "font-awesome")),
      #shinyBS::bsTooltip("Mode", title = "Edit all input parameters", placement = "top", trigger = "hover"),
      align = "center"
    ),
    shiny::br(),
    fluidRow(
      actionButton("showInputside", label = "", icon = icon("wpforms")),
      shinyBS::bsTooltip("showInputside", title = "Only show input parameters", placement = "bottom", trigger = "hover"),
      actionButton("showBothside", label = "", icon = icon("table-columns", lib = "font-awesome")),
      shinyBS::bsTooltip("showBothside", title = "Show both,  input parameters and output", placement = "bottom", trigger = "hover"),
      actionButton("showOutputside", label = "", icon = icon("chart-line")),
      shinyBS::bsTooltip("showOutputside", title = "Only show output", placement = "bottom", trigger = "hover"),
      align = "center"
    ),
    shiny::br(),
    shiny::fluidRow(
      column(6, id = "inputside", inputUi),
      column(6, id = "outputside", outputUi)
    ),
    # shiny::sidebarLayout(inputUi, outputUi),
    shiny::fluidRow(shiny::div(
      shiny::actionButton(inputId = "generateLandscape", label = "Generate the landscape"),
      align = "center"
    )),
    shinyBS::bsTooltip("generateLandscape", title = GENERATE_LANDSCAPE, placement = "top", trigger = "hover"),
    shiny::br(),
    shiny::fluidRow(shiny::div(
      shiny::actionButton(inputId = "runSimulation", label = "Run simulation"),
      shiny::actionButton(inputId = "stopSimulation", label = "Stop simulation"),
      align = "center"
    )),
    shinyBS::bsTooltip("runSimulation", title = RUN_SIMULATION, placement = "top", trigger = "hover"),
    shinyBS::bsTooltip("stopSimulation", title = STOP_SIMULATION, placement = "top", trigger = "hover"),
    shiny::br()#,
    # shiny::fluidRow(shiny::div(
    #   shiny::downloadButton(outputId = "export", label = "Download simulation parameters"),
    #   align = "center"
    # )),
    # shinyBS::bsTooltip("export", title = EXPORT_SIMULATION, placement = "top", trigger = "hover"),
    # shiny::br()
  )
}
