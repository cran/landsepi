simul_params <- createSimulParams(outputDir = paste0(ROOT_PATH, "/www/tmp/"))

## Pathogen parameters
simul_params <- setPathogen(simul_params, loadPathogen("rust"))
## Initial conditions
simul_params <- setInoculum(simul_params, 5e-4)

## Outputs
simul_params <- setOutputs(simul_params, list(
  epid_outputs = "audpc_rel", evol_outputs = "",
  thres_breakdown = 50000,
  GLAnoDis = 1.48315,
  audpc100S = 0.76
))


# Server
server <- function(input, output, session) {

  ##########################
  # Buttons enable/disable
  ##########################
  # reactive values that indicate if values is ok
  can_gen_landscape <- shiny::reactiveValues(
    proportions = TRUE,
    croptypeID = TRUE,
    rotation = TRUE,
    seed = TRUE
  )

  can_run_simul <- shiny::reactiveValues(
    landscape = FALSE,
    seed = TRUE,
    nYear = TRUE,
    nTSpY = TRUE,
    croptypes = TRUE,
    cultivars = TRUE,
    cultivarsgenes = TRUE,
    genes = TRUE,
    patho_infectious_rate = TRUE,
    patho_survival_prob = TRUE,
    patho_repro_sex_prob = TRUE,
    patho_propagule_prod_rate = TRUE,
    patho_latent_period_exp = TRUE,
    patho_latent_period_var = TRUE,
    patho_infectious_period_exp = TRUE,
    patho_infectious_period_var = TRUE,
    patho_sigmoid_kappa = TRUE,
    patho_sigmoid_sigma = TRUE,
    patho_sigmoid_plateau = TRUE,
    inoculum = TRUE
  )

  ## Observe reactiveValues

  # Can run landscape generation
  # Enable / disable buttons
  observe({
    if (can_gen_landscape$proportions &&
      can_gen_landscape$croptypeID &&
      can_gen_landscape$rotation &&
      can_gen_landscape$seed
    ) {
      shinyjs::enable(id = "generateLandscape")
    }
    else {
      shinyjs::disable(id = "generateLandscape")
      shinyjs::disable(id = "runSimulation")
      shinyjs::disable(id = "stopSimulation")
      shinyjs::disable(id = "export")
    }
  })

  # Can run simulation
  # Enable / disable buttons
  observe({
    if (can_run_simul$landscape &&
      can_run_simul$seed &&
      can_run_simul$nYear &&
      can_run_simul$nTSpY &&
      can_run_simul$croptypes &&
      can_run_simul$cultivars &&
      can_run_simul$cultivarsgenes &&
      can_run_simul$genes &&
      can_run_simul$patho_infectious_rate &&
      can_run_simul$patho_survival_prob &&
      can_run_simul$patho_repro_sex_prob &&
      can_run_simul$patho_propagule_prod_rate &&
      can_run_simul$patho_latent_period_exp &&
      can_run_simul$patho_latent_period_var &&
      can_run_simul$patho_infectious_period_exp &&
      can_run_simul$patho_infectious_period_var &&
      can_run_simul$patho_sigmoid_kappa &&
      can_run_simul$patho_sigmoid_sigma &&
      can_run_simul$patho_sigmoid_plateau &&
      can_run_simul$inoculum
    ) {
      shinyjs::enable(id = "runSimulation")
      shinyjs::disable(id = "stopSimulation")
      shinyjs::enable(id = "export")
    }
    else {
      shinyjs::disable(id = "runSimulation")
      shinyjs::disable(id = "stopSimulation")
      shinyjs::disable(id = "export")
    }
  })

  # Test if the croptypes proportion sum is 1
  ## TODO remove input$ and move to global.R
  ProportionValidation <- function() {
    if ((advanced_mode() == FALSE && input$demo == "RO") || (advanced_mode() && !is.na(input$rotationPeriod) && input$rotationPeriod > 0)) {
      sum_prop <-
        ((croptypes_proportions()[1] + croptypes_proportions()[2]) + (croptypes_proportions()[1] + croptypes_proportions()[3])) / 2
    }
    else {
      sum_prop <- sum(as.numeric(croptypes_proportions()))
    }

    shiny::removeUI(selector = "#propError")
    if (!isTRUE(all.equal(sum_prop, 1)) ||
      is.na(sum_prop)) {
      showErrorMessage(
        id = "propError", selectorafter = "#generateLandscape",
        message = "The sum of the proportions of all croptypes must be equal to 1 (100%)"
      )
      return(invisible(FALSE))
    }
    return(invisible(TRUE))
  }
  
  ## Print Rotation labels
  setRotationText <- function(list_name=NULL) {
    text <- paste0("<u>1st configuration</u> : croptypes 0 (<b>",list_name[1],"</b>) and 1 (<b>",list_name[2],"</b>)")
    text <- paste0(text, "<br/><u>2nd configuration</u> : croptypes 0 (<b>",list_name[1],"</b>) and 2 (<b>",list_name[3],"</b>)")
    text  
  }

  #############################################################################
  # Observe EVENT
  #############################################################################

  # About
  # Modal Dialog
  observeEvent(input$About, {
    showModal(modalDialog(
      title = paste0("About : Landsepi V", packageVersion("landsepi")),
      easyClose = TRUE,
      size = "l",
      footer = NULL,
      div(HTML("<h1>Landsepi: Landscape Epidemiology and Evolution</h1><img src='landsepi-logo.png'  align='right' alt='' width='120'/>
                          <h3> A stochastic, spatially-explicit, demo-genetic model
                         simulating the spread and evolution of a plant pathogen in a heterogeneous landscape
                         to assess resistance deployment strategies. It is based on a spatial geometry for describing
                         the landscape and allocation of different cultivars, a dispersal kernel for the
                         dissemination of the pathogen, and a SEIR ('Susceptible-Exposed-Infectious-Removed’)
                         structure with a discrete time step. It provides a useful tool to assess the performance
                         of a wide range of deployment options with respect to their epidemiological,
                         evolutionary and economic outcomes.</h3>
                          <h3> Authors:</h3> J-L Gaussen, J. Papaïx, J-F Rey, L. Rimbaud
                          <h3>Package project:</h3><a href='https://CRAN.R-project.org/package=landsepi' target='_blank'> CRAN package</a><br/><a href='https://gitlab.paca.inra.fr/CSIRO-INRA/landsepi' target='_blank'> Source code</a>
                          <br/> License GPL-3
                          <h3> How to cite the package:</h3> <b>Rimbaud L, Papaïx J, Rey J-F (2019).</b> landsepi: Landscape Epidemiology and Evolution. R package version 0.1.0, &lt;URL: https://cran.r-project.org/package=landsepi&gt;.
                          <h3> Full model description:</h3> <b>Rimbaud L, Papaïx J, Rey J-F, Barrett LG and Thrall PH. 2018.</b> Assessing the durability and efficiency of landscape-based strategies to deploy plant resistance to pathogens. PLoS Computational Biology 14(4): e1006067. <a href='https://doi.org/10.1371/journal.pcbi.1006067' target='_blank'>https://doi.org/10.1371/journal.pcbi.1006067</a>
                          <div>
                          <img src='Republique_Francaise_RVB.jpg' alt='RF' style='width:50px; margin-left: 10px;' />
                          <img src='LogoINRAE_Fr_rvb_web.png' alt='INRAE' style='width:50px; margin-left: 10px;' />
                          <img src='logoBIOSP.jpeg' alt='BioSP' style='width:50px; margin-left: 10px;'/>
                          <img src='PATHO_inra_logo.png' alt='Pathologie végétale' style='width:50px; margin-left: 10px;'/>
                          <img src='CSIRO_Logo.png' alt='CSIRO' style='width:40px; margin-left: 10px;'/>
                          </div>
                "))
    ))
  })

  ## User Mode button switch between mode
  observeEvent(input$Mode, {
    advanced_mode(!advanced_mode())
    if (advanced_mode()) {
      printVerbose("enable mode edition", level=3)
      removeCssClass("Mode","btn-default")
      updateTabsetPanel(session,inputId = "inputtabpanel", selected = "Cultivars and Genes")
      shinyjs::disable(id = "demo")
      shinyjs::disable(id = "rotationPeriod")
      shiny::updateNumericInput(session, "rotationPeriod", value = 0)
      shinyjs::enable(id = "patho_infection_rate")
      shinyjs::enable(id = "patho_propagule_prod_rate")
      shinyjs::enable(id = "patho_latent_period_exp")
      shinyjs::enable(id = "patho_latent_period_var")
      shinyjs::enable(id = "patho_infectious_period_exp")
      shinyjs::enable(id = "patho_infectious_period_var")
      shinyjs::enable(id = "patho_survival_prob")
      shinyjs::enable(id = "patho_repro_sex_prob")
      shinyjs::enable(id = "patho_sigmoid_kappa")
      shinyjs::enable(id = "patho_sigmoid_sigma")
      shinyjs::enable(id = "patho_sigmoid_plateau")
    }
    else {
      printVerbose("disable mode edition", level=3)
      addCssClass("Mode","btn-default")
      shinyjs::disable(id = "rotationPeriod")
      shiny::updateNumericInput(session, "rotationPeriod", value = 0)
      shinyjs::enable(id = "demo")
      shinyjs::disable(id = "patho_infection_rate")
      shinyjs::disable(id = "patho_propagule_prod_rate")
      shinyjs::disable(id = "patho_latent_period_exp")
      shinyjs::disable(id = "patho_latent_period_var")
      shinyjs::disable(id = "patho_infectious_period_exp")
      shinyjs::disable(id = "patho_infectious_period_var")
      shinyjs::disable(id = "patho_survival_prob")
      shinyjs::disable(id = "patho_repro_sex_prob")
      shinyjs::disable(id = "patho_sigmoid_kappa")
      shinyjs::disable(id = "patho_sigmoid_sigma")
      shinyjs::disable(id = "patho_sigmoid_plateau")
    }
  })

  # Inputs / Outputs
  ######################################################################################
  # Landscape
  shiny::observeEvent(input$landscape, {
    can_gen_landscape$proportions <<- ProportionValidation()
    can_run_simul$landscape <<- FALSE

    shinyjs::show(id = "landscapeimg")
    output$landscapeimg <- renderPlot({
      plot(loadLandscape(input$landscape))
    })
  })
  shiny::observeEvent(input$aggregLevel, {
    can_gen_landscape$proportions <<- ProportionValidation()
    can_run_simul$landscape <<- FALSE
  })

  ######################################################################################
  # Rotation period validation
  shiny::observeEvent(input$rotationPeriod, {
    can_gen_landscape$rotation <<- TRUE
    can_run_simul$landscape <<- FALSE
    shiny::removeUI(selector = "#rotationPeriodError")
    if (input$demo == "RO" && advanced_mode() == FALSE) {
      if (input$rotationPeriod < 1 ||
        input$rotationPeriod > input$nYear ||
        is.na(input$rotationPeriod)) {
        showErrorMessage(
          id = "rotationPeriodError", selectorafter = "#generateLandscape",
          message = paste0(
            "The rotation period should be between 1 and ", input$nYear, " (the simulation duration)"
          )
        )
        can_gen_landscape$rotation <<- FALSE
      }
    } else {
      if (input$rotationPeriod < 0 ||
          input$rotationPeriod > input$nYear ||
          is.na(input$rotationPeriod)) {
        showErrorMessage(
          id = "rotationPeriodError", selectorafter = "#generateLandscape",
          message = paste0(
            "The rotation period should be between 1 and ", input$nYear, " (the simulation duration) or 0 for none"
          )
        )
        can_gen_landscape$rotation <<- FALSE
      }
    }
    can_gen_landscape$proportions <<- ProportionValidation()
    can_run_simul$landscape <<- FALSE
  })
  ######################################################################################
  # nYear validation
  shiny::observeEvent(input$nYear, {
    can_run_simul$landscape <<- FALSE

    shiny::removeUI(selector = "#nYearError")
    if (input$nYear < 1 || input$nYear > 100 || is.na(input$nYear)) {
      showErrorMessage(
        id = "nYearError", selectorafter = "#generateLandscape",
        message = "The simulation duration should be between 1 and 100"
      )
      can_run_simul$nYear <<- FALSE
    }
    else {
      simul_params <<- setTime(simul_params, Nyears = input$nYear, nTSpY = input$nTSpY)
      can_run_simul$nYear <<- TRUE
    }
  })
  ######################################################################################
  # nTSpY validation
  shiny::observeEvent(input$nTSpY, {
    can_run_simul$landscape <<- FALSE

    shiny::removeUI(selector = "#nTSpYError")
    if (input$nTSpY < 1 || input$nTSpY > 365 || is.na(input$nTSpY)) {
      showErrorMessage(
        id = "nTSpYError", selectorafter = "#generateLandscape",
        message = "The time step should be between 1 and 365"
      )
      can_run_simul$nTSpY <<- FALSE
    }
    else {
      simul_params <<- setTime(simul_params, Nyears = input$nYear, nTSpY = input$nTSpY)
      can_run_simul$nTSpY <<- TRUE
    }
  })
  ######################################################################################
  # seed validation
  shiny::observeEvent(input$seed, {
    can_run_simul$seed <<- TRUE
    can_run_simul$landscape <<- FALSE


    shiny::removeUI(selector = "#seedError")
    if (input$seed < 0 || input$seed > 99999 || is.na(input$seed)) {
      showErrorMessage(
        id = "seedError", selectorafter = "#generateLandscape",
        message = "The seed value should be between 0 and 99999"
      )
      can_gen_landscape$seed <<- FALSE
      can_run_simul$seed <<- FALSE
    }
    else {
      simul_params <<- setSeed(simul_params, input$seed)
      can_run_simul$seed <<- TRUE
      can_gen_landscape$seed <<- TRUE
    }
  })

  ######################################################################################
  #
  # Patho Tab Observe
  #
  ######################################################################################
  # Select Pathogen
  shiny::observeEvent(input$defaultPathogen, {
    simul_params <<- setPathogen(simul_params, loadPathogen(disease = tolower(input$defaultPathogen)))
    updateNumericInput(session = session, inputId = "patho_survival_prob", value = simul_params@Pathogen$survival_prob)
    updateNumericInput(session = session, inputId = "patho_repro_sex_prob", value = simul_params@Pathogen$repro_sex_prob)
     updateNumericInput(session = session, inputId = "patho_infection_rate", value = simul_params@Pathogen$infection_rate)
    updateNumericInput(session = session, inputId = "patho_propagule_prod_rate", value = simul_params@Pathogen$propagule_prod_rate)
    updateNumericInput(session = session, inputId = "patho_latent_period_exp", value = simul_params@Pathogen$latent_period_exp)
    updateNumericInput(session = session, inputId = "patho_latent_period_var", value = simul_params@Pathogen$latent_period_var)
    updateNumericInput(session = session, inputId = "patho_infectious_period_exp", value = simul_params@Pathogen$infectious_period_exp)
    updateNumericInput(session = session, inputId = "patho_infectious_period_var", value = simul_params@Pathogen$infectious_period_var)
    updateNumericInput(session = session, inputId = "patho_sigmoid_kappa", value = simul_params@Pathogen$sigmoid_kappa)
    updateNumericInput(session = session, inputId = "patho_sigmoid_sigma", value = simul_params@Pathogen$sigmoid_sigma)
    updateNumericInput(session = session, inputId = "patho_sigmoid_plateau", value = simul_params@Pathogen$sigmoid_plateau)
  })

  # inoculum
  shiny::observeEvent(input$inoculum, {
    shiny::removeUI(selector = "#pathoInoculumError")
    if (input$inoculum > 1 || input$inoculum < 0 || is.na(input$inoculum)) {
      showErrorMessage(
        id = "pathoInoculumError", selectorafter = "#generateLandscape",
        message = "The probability of initial infection should be between 0 and 1"
      )
      can_run_simul$inoculum <<- FALSE
    }
    else {
      simul_params <<- setInoculum(simul_params, input$inoculum)
      can_run_simul$inoculum <<- TRUE
    }
  })

  # survival prob
  shiny::observeEvent(input$patho_survival_prob, {
    shiny::removeUI(selector = "#pathoSurProbError")
    if (input$patho_survival_prob > 1 || input$patho_survival_prob < 0 || is.na(input$patho_survival_prob)) {
      showErrorMessage(
        id = "pathoSurProbError", selectorafter = "#generateLandscape",
        message = "The probability for a propagule to survive the off-season should be between 0 and 1"
      )
      can_run_simul$patho_survival_prob <<- FALSE
    }
    else {
      simul_params@Pathogen$survival_prob <<- input$patho_survival_prob
      can_run_simul$patho_survival_prob <<- TRUE
    }
  })

  # repro_sex_prob
  # shiny::observeEvent(input$patho_repro_sex_prob, {
  #   shiny::removeUI(selector = "#pathoReproSexProbError")
  #   if (input$patho_repro_sex_prob > 1 || input$patho_repro_sex_prob < 0 || is.na(input$patho_repro_sex_prob)) {
  #     showErrorMessage(
  #       id = "pathoReproSexProbError", selectorafter = "#generateLandscape",
  #       message = "The probability for an infectious host to reproduce via sex rather than via cloning should be between 0 and 1"
  #     )
  #     can_run_simul$patho_repro_sex_prob <<- FALSE
  #   }
  #   else {
  #     simul_params@Pathogen$repro_sex_prob <<- input$patho_repro_sex_prob
  #     can_run_simul$patho_repro_sex_prob <<- TRUE
  #   }
  # })

  # propagule_prod_rate
  shiny::observeEvent(input$patho_propagule_prod_rate, {
    shiny::removeUI(selector = "#pathoProdRateError")
    if (input$patho_propagule_prod_rate > VALUEMAX || input$patho_propagule_prod_rate < 0 || is.na(input$patho_propagule_prod_rate)) {
      showErrorMessage(
        id = "pathoProdRateError", selectorafter = "#generateLandscape",
        message = paste0("The maximal expected effective propagule production rate of an infectious host per time step should be between 0 and ",VALUEMAX)
      )
      can_run_simul$patho_propagule_prod_rate <<- FALSE
    }
    else {
      simul_params@Pathogen$propagule_prod_rate <<- input$patho_propagule_prod_rate
      can_run_simul$patho_propagule_prod_rate <<- TRUE
    }
  })

  # latent_period_exp
  shiny::observeEvent(input$patho_latent_period_exp, {
    shiny::removeUI(selector = "#pathoLatPerExpError")
    if (input$patho_latent_period_exp > VALUEMAX || input$patho_latent_period_exp < 0 || is.na(input$patho_latent_period_exp)) {
      showErrorMessage(
        id = "pathoLatPerExpError", selectorafter = "#generateLandscape",
        message = paste0("The minimal expected duration of the latent period should be between 0 and ",VALUEMAX)
      )
      can_run_simul$patho_latent_period_exp <<- FALSE
    }
    else {
      simul_params@Pathogen$latent_period_exp <<- input$patho_latent_period_exp
      can_run_simul$patho_latent_period_exp <<- TRUE
    }
  })

  # latent_period_var
  shiny::observeEvent(input$patho_latent_period_var, {
    shiny::removeUI(selector = "#pathoLatPerVarError")
    if (input$patho_latent_period_var > VALUEMAX || input$patho_latent_period_var < 0 || is.na(input$patho_latent_period_var)) {
      showErrorMessage(
        id = "pathoLatPerVarError", selectorafter = "#generateLandscape",
        message = paste0("The variance of the infectious period duration should be between 0 and ",VALUEMAX)
      )
      can_run_simul$patho_latent_period_var <<- FALSE
    }
    else {
      simul_params@Pathogen$latent_period_var <<- input$patho_latent_period_var
      can_run_simul$patho_latent_period_var <<- TRUE
    }
  })

  # infectious_period_exp
  shiny::observeEvent(input$patho_infectious_period_exp, {
    shiny::removeUI(selector = "#pathoInfPerExpError")
    if (input$patho_infectious_period_exp > VALUEMAX || input$patho_infectious_period_exp < 0 || is.na(input$patho_infectious_period_exp)) {
      showErrorMessage(
        id = "pathoInfPerExpError", selectorafter = "#generateLandscape",
        message = paste0("The maximal expected duration of the infectious period should be between 0 and ",VALUEMAX)
      )
      can_run_simul$patho_infectious_period_exp <<- FALSE
    }
    else {
      simul_params@Pathogen$infectious_period_exp <<- input$patho_infectious_period_exp
      can_run_simul$patho_infectious_period_exp <<- TRUE
    }
  })

  # infectious_period_var
  shiny::observeEvent(input$patho_infectious_period_var, {
    shiny::removeUI(selector = "#pathoInfPerVarError")
    if (input$patho_infectious_period_var > VALUEMAX || input$patho_infectious_period_var < 0 || is.na(input$patho_infectious_period_var)) {
      showErrorMessage(
        id = "pathoInfPerVarError", selectorafter = "#generateLandscape",
        message = paste0("The variance of the infectious period duration should be between 0 and ",VALUEMAX)
      )
      can_run_simul$patho_infectious_period_var <<- FALSE
    }
    else {
      simul_params@Pathogen$infectious_period_var <<- input$patho_infectious_period_var
      can_run_simul$patho_infectious_period_var <<- TRUE
    }
  })

  # sigmoid_kappa
  shiny::observeEvent(input$patho_sigmoid_kappa, {
    shiny::removeUI(selector = "#pathoSigKapError")
    if (input$patho_sigmoid_kappa > VALUEMAX || input$patho_sigmoid_kappa < 0 || is.na(input$patho_sigmoid_kappa)) {
      showErrorMessage(
        id = "pathoSigKapError", selectorafter = "#generateLandscape",
        message = paste0("The kappa parameter of the sigmoid contamination function should be between 0 and ",VALUEMAX)
      )
      can_run_simul$patho_sigmoid_kappa <<- FALSE
    }
    else {
      simul_params@Pathogen$sigmoid_kappa <<- input$patho_sigmoid_kappa
      can_run_simul$patho_sigmoid_kappa <<- TRUE
    }
  })

  # sigmoid_sigma
  shiny::observeEvent(input$patho_sigmoid_sigma, {
    shiny::removeUI(selector = "#pathoSigSigError")
    if (input$patho_sigmoid_sigma > VALUEMAX || input$patho_sigmoid_sigma < 0 || is.na(input$patho_sigmoid_sigma)) {
      showErrorMessage(
        id = "pathoSigSigError", selectorafter = "#generateLandscape",
        message = paste0("The sigma parameter of the sigmoid contamination function should be between 0 and ",VALUEMAX)
      )
      can_run_simul$patho_sigmoid_sigma <<- FALSE
    }
    else {
      simul_params@Pathogen$sigmoid_sigma <<- input$patho_sigmoid_sigma
      can_run_simul$patho_sigmoid_sigma <<- TRUE
    }
  })

  # sigmoid_plateau
  # shiny::observeEvent(input$patho_sigmoid_plateau, {
  #   shiny::removeUI(selector = "#pathoSigPlaError")
  #   if (input$patho_sigmoid_plateau > 10 || input$patho_sigmoid_plateau < 0 || is.na(input$patho_sigmoid_plateau)) {
  #     showErrorMessage(
  #       id = "pathoSigPlaError", selectorafter = "#generateLandscape",
  #       message = "The plateau parameter of the sigmoid contamination function should be between 0 and ?"
  #     )
  #     can_run_simul$patho_sigmoid_plateau <<- FALSE
  #   }
  #   else {
  #     simul_params@Pathogen$sigmoid_plateau <<- input$patho_sigmoid_plateau
  #     can_run_simul$patho_sigmoid_plateau <<- TRUE
  #   }
  # })

  # infection rate validation
  shiny::observeEvent(input$patho_infection_rate, {
    shiny::removeUI(selector = "#pathoInfRateError")
    if (input$patho_infection_rate > 1 || input$patho_infection_rate < 0 || is.na(input$patho_infection_rate)) {
      showErrorMessage(
        id = "pathoInfRateError", selectorafter = "#generateLandscape",
        message = "The maximal expected infection rate of a propagule on a healthy host should be between 0 and 1"
      )
      can_run_simul$path_infection_rate <<- FALSE
    }
    else {
      simul_params@Pathogen$infection_rate <<- input$patho_infection_rate
      can_run_simul$patho_infection_rate <<- TRUE
    }
  })


  ######################################################################
  # Handle the download gpkg button
  ######################################################################
  output$export <-
    shiny::downloadHandler(
      filename = "landsepi_landscape.gpkg",
      content <- function(file) {
        simul_params <<- saveDeploymentStrategy(simul_params)
        file.copy(file.path(simul_params@OutputDir, simul_params@OutputGPKG), file)
      },
      contentType = "application/x-sqlite3"
    )
  ######################################################################################
  # Handle the "Generate the landscape" button
  shiny::observeEvent(input$generateLandscape, {
    withProgress(message = "Generating Landscape, Please wait...", value = 0, {
      shinyjs::disable(id = "generateLandscape")
      shinyjs::disable(id = "export")
      output$video <- NULL
      # if(!dir.exists(paste0(ROOT_PATH,"/www/tmp/"))) dir.create(paste0(ROOT_PATH,"/www/tmp/"))
      # setwd(paste0(ROOT_PATH,"/www/tmp/"))

      # Remove old files
      cleanDir(simul_params@OutputDir)

      #print(simul_params@Croptypes)
      # print(simul_params@Cultivars)
      # print(simul_params@CultivarsGenes)
      # print(simul_params@Genes)

      # Croptypes Rotation
      if ( (advanced_mode() == FALSE && input$demo == "RO") || (advanced_mode() && input$rotationPeriod > 0)) {
        rotation_period <- input$rotationPeriod
        prop <- list(
          c(croptypes_proportions()[1], croptypes_proportions()[2]),
          c(croptypes_proportions()[1], croptypes_proportions()[3])
        )
        # aggregLevel = strtoi(input$aggregLevel)
        rotation_sequence <- list(
          c(simul_params@Croptypes$croptypeID[1], simul_params@Croptypes$croptypeID[2]),
          c(simul_params@Croptypes$croptypeID[1], simul_params@Croptypes$croptypeID[3])
        )
      }
      else {
        rotation_period <- 0
        rotation_sequence <- list(c(simul_params@Croptypes$croptypeID))
        if (input$demo == "PY") {
          prop <- list(croptypes_proportions()[1:2])
        } else {
          prop <- list(croptypes_proportions())
        }
      }

      simul_params <<- setSeed(simul_params, input$seed)

      incProgress(0.4)
      # Run the landscape generation
      simul_params <<- setLandscape(simul_params, loadLandscape(input$landscape))
      ## Dispersal parameters
      simul_params <<- setDispersalPathogen(simul_params, loadDispersalPathogen(input$landscape))
      ## Dispersal parameters
      disp_host <- loadDispersalHost(simul_params, type = "no")
      simul_params <<- setDispersalHost(simul_params, disp_host)

      ## Define the value of aggreg from aggregLevel

      switch(input$aggregLevel,
        "low" = {
          aggreg <- 0.07
          algo <- "periodic"
        },
        "medium" = {
          aggreg <- 0.25
          algo <- "exp"
        },
        "high" = {
          aggreg <- 10
          algo <- "periodic"
        },
        {
          aggreg <- 0.25
          algo <- "exp"
        }
      )

      simul_params <<- allocateLandscapeCroptypes(simul_params,
        rotation_period = rotation_period,
        rotation_sequence = rotation_sequence,
        rotation_realloc = FALSE,
        prop = prop,
        aggreg = aggreg,
        algo = algo,
        graphic = TRUE
      )

      setwd(ROOT_PATH)

      incProgress(0.5)
      # Print the image of the landscape
      # TODO Loop images for rotation demo
      # output$landscape <- shiny::renderImage({
      #  list(
      #    src = file.path("www/tmp", "landscape_year1.png"),
      #    contentType = 'image/png',
      #    width = "70%",
      #    height = "auto",
      #   alt = "Landscape"
      #  )
      # }, deleteFile = FALSE)

      shinyjs::show(id = "landscapeimg")
      output$landscapeimg <- renderPlot({
        imgs <- normalizePath(list.files(simul_params@OutputDir, pattern = ".png", full.names = TRUE))
        pngs <- lapply(imgs, readPNG)
        asGrobs <- lapply(pngs, rasterGrob)
        p <- grid.arrange(grobs = asGrobs, nrow = 1)
      })

      # Using slick : trouble with images size...
      # output$landscape <- renderSlickR({
      #  imgs <- list.files("www",pattern=".png",full.names = TRUE)
      #  slickR(imgs, slickOpts=list(adaptiveHeight=FALSE, respondTo="min"), height = "auto")
      # })

      shinyjs::enable(id = "generateLandscape")
      shinyjs::enable(id = "export")
      can_run_simul$landscape <<- TRUE
      shinyjs::click("showBothside")
    })
  })

  #############################################################################
  # Stop simulation : Kill future promise
  future_process <- NULL
  observeEvent(input$stopSimulation, {
    cat(file = stderr(), "STOP button -> stop process ", future_process$job$pid, "\n")
    tools::pskill(future_process$job$pid, signal = tools::SIGTERM)
    tools::pskill(future_process$job$pid, signal = tools::SIGKILL)
  })

  ######################################################################################
  # Handle the "Run simulation" button

  shiny::observeEvent(input$runSimulation, {
    
    printVerbose(simul_params, level=2)
    withProgress(message = "Running Simulation, please wait...", value = 0, {
      progressBar <- Progress$new()
      progressBar$set(value = NULL, message = "Running Simulation, please wait...")
      # setwd(paste0(ROOT_PATH,"/www/tmp/"))

      shinyjs::disable(id = "generateLandscape")
      shinyjs::disable(id = "runSimulation")
      shinyjs::disable(id = "export")
      shinyjs::enable(id = "stopSimulation")
      
      shinyjs::disable("showInputside")
      shinyjs::disable("showBothside")
      #shinyjs::click("showOutputside") # seems not working -> force it
      shinyjs::showElement(id = "outputside")
      shinyjs::hideElement(id = "inputside")
      removeCssClass("inputside", "col-sm-12")
      removeCssClass("inputside", "col-sm-6")
      addCssClass("inputside", "col-sm-0")
      removeCssClass("outputside", "col-sm-0")
      removeCssClass("outputside", "col-sm-6")
      addCssClass("outputside", "col-sm-12")
      shinyBS::removeTooltip(session, "runSimulation") ## avoid tooltip to stay active

      progressBar$set(value = 0.4)

      plan(list(multicore, multisession))

      future_process <<- future({
        res <- landsepi::runSimul(simul_params,
          graphic = FALSE, videoMP4 = TRUE
        )
      })

      then(future_process,
        onFulfilled = function(value) {
          progressBar$set(value = 0.8, message = "Simulation ended : making video...")

          shinyjs::enable(id = "generateLandscape")
          # shinyjs::enable(id = "runSimulation")
          shinyjs::enable(id = "export")
          shinyjs::enable(id = "runSimulation")
          shinyjs::disable(id = "stopSimulation")

          output$landscapeimg <- NULL
          hide(id = "landscapeimg")

          output$video <-
            shiny::renderUI(
              tags$video(
                id = "video",
                type = "video/mp4",
                src = paste0("tmp/", basename(simul_params@OutputDir), "/video.mp4?rand=", as.integer(Sys.time())),
                controls = "controls",
                width = "100%",
                height = "auto"
              )
            )
          shinyjs::enable("showInputside")
          shinyjs::enable("showBothside")
          shinyjs::click("showOutputside")
          shinyBS::addTooltip(session, "runSimulation",title=RUN_SIMULATION, placement = "top", trigger="hover")
          # setwd(dirname(getwd()))
        },
        onRejected = function(err) {
          setwd(ROOT_PATH)
          cat(file = stderr(), "\n ### KILL'EM ALL ### -> Kill simulation \n")
          cleanDir(simul_params@OutputDir)
          shinyjs::enable("showInputside")
          shinyjs::enable("showBothside")
          shinyjs::click("showBothside")
          shinyjs::enable(id = "generateLandscape")
          shinyBS::addTooltip(session, "runSimulation",title=RUN_SIMULATION, placement = "top", trigger="hover")
          can_run_simul$landscape <<- FALSE
          future_process <- NULL
          
          shinyalert::shinyalert(
            "Oups! Something went wrong !",
            "Please check inputs",
            type = "error",
            size = "m",
            closeOnEsc = TRUE,
            showCancelButton = TRUE, showConfirmButton = FALSE
          )
          0
        }
      ) %>%
        finally(~ progressBar$close())

      setwd(ROOT_PATH)
    })
  })
  ######################################################################################
  # Handle the demo list
  shiny::observeEvent(input$demo, {
    # Cultivar tab
    switch(input$demo,
      MO = {
        simul_params <<- loadDemoMO(simul_params)
      },
      MI = {
        simul_params <<- loadDemoMI(simul_params)
      },
      RO = {
        simul_params <<- loadDemoRO(simul_params)
      },
      PY = {
        simul_params <<- loadDemoPY(simul_params)
      },
      {
        # Default case
        print("input$demo : Unknown input$demo")
      }
    )

    default_gene <<- simul_params@Genes[1,]
    default_cultivar <<- simul_params@Cultivars[1,]
    default_croptype <<- simul_params@Croptypes[1,c(1,2)]
    
    simul_params_croptypes(simul_params@Croptypes)
    simul_params_cultivars(simul_params@Cultivars)
    simul_params_cultivarsgenes(simul_params@CultivarsGenes)
    simul_params_genes(simul_params@Genes)
    checkAllTables()

    can_gen_landscape$proportions <<- TRUE
    can_gen_landscape$croptypeID <<- TRUE
    can_gen_landscape$rotation <<- TRUE
    can_gen_landscape$seed <<- TRUE
    can_run_simul$landscape <<- FALSE
    can_run_simul$seed <<- TRUE
    can_run_simul$nYear <<- TRUE
    can_run_simul$nTSpY <<- TRUE
    can_run_simul$croptypes <<- TRUE
    can_run_simul$cultivars <<- TRUE
    can_run_simul$cultivarsgenes <<- TRUE
    can_run_simul$genes <<- TRUE

    # Landscape tab
    shiny::updateSelectInput(session, "landscape", selected = 1)
    simul_params <<- setLandscape(simul_params, loadLandscape(1))
    shiny::updateSelectInput(session, "aggregLevel", selected = "low")

    # Enable all the conditionnal inputs by default, we disable it later if needed
    shinyjs::disable(id = "rotationPeriod")
    shiny::updateNumericInput(session, "rotationPeriod", value = 0)

    if (input$demo == "MO") {
      croptypes_proportions(c(0.33, 0.33, 0.34))
      shiny::updateSelectInput(session, "aggregLevel", selected = "high")
    }
    else if (input$demo == "MI" || input$demo == "PY") {
      croptypes_proportions(c(0.5, 0.5))
      shiny::updateSelectInput(session, "aggregLevel", selected = "high")
    }
    else if (input$demo == "RO") {
      croptypes_proportions(c(0.5, 0.5, 0.5))
      shinyjs::enable(id = "rotationPeriod")
      shiny::updateNumericInput(session, "rotationPeriod", value = 2)
      shiny::updateSelectInput(session, "aggregLevel", selected = "medium")
    }
    else if (input$demo == "PY") {
      shiny::updateSelectInput(session, "aggregLevel", selected = "low")
    }
    
    output$rotationText <- renderUI({HTML(setRotationText(simul_params@Croptypes[,2]))})
  })

  ###################
  ### TABS TABLES ###
  ###################

  #### croptypes table ####
  simul_params_croptypes(simul_params@Croptypes)
  croptypesTable <- editableDTServer(
    id = "croptypes",
    DTdata = reactive({
      return(cbind(simul_params_croptypes(), data.frame(Proportions = croptypes_proportions())))
    }),
    disableCol = shiny::reactive({
      if (isTRUE(advanced_mode())) {
        c("croptypeID")
      } else {
        names(simul_params_croptypes())
        # print(names(simul_params_croptypes()))
      }
    }),
    canRm = advanced_mode,
    rownames = FALSE,
    tooltips = c("Croptype index (starts at 0)","Croptype name"),
    row.default = default_croptype,
    row.cols = 1:2,
    row.inc = c(1,2)
  )

  ##### croptypes table modification #####
  shiny::observeEvent(croptypesTable$data,
    {
      message("Croptypes update")

      if (sum(is.na(croptypesTable$value))) {
        return(1)
      }

      # message("data ", croptypesTable$data)
      # message("value ", croptypesTable$value)
      # message("i ", croptypesTable$row)
      # message("j ", croptypesTable$col)

      croptypes_proportions(croptypesTable$data[, "Proportions"])
      can_gen_landscape$proportions <<- ProportionValidation()
      
      if (can_gen_landscape$proportions == FALSE) can_run_simul$landscape <<- FALSE

      if (isTRUE(advanced_mode())) {
        shiny::isolate(simul_params_croptypes(croptypesTable$data[, 1:(ncol(croptypesTable$data) - 2)]))
        
        if (nrow(croptypesTable$data) == 0 
            || checkCroptypesTable(croptypesTable$data[, -which(colnames(croptypesTable$data) %in% c("Proportions", "delete"))]) == FALSE) {
          can_run_simul$croptypes <<- FALSE
          can_gen_landscape$croptypeID <<- FALSE
        }
        else {
          croptypesTable$data[,"croptypeID"] <- seq(1:nrow(croptypesTable$data))-1
          output$rotationText <- renderUI({HTML(setRotationText(croptypesTable$data[, 2]))})
          simul_params <<- setCroptypes(simul_params, croptypesTable$data[, 1:(ncol(croptypesTable$data) - 2)])
          can_run_simul$croptypes <<- TRUE
          can_gen_landscape$croptypeID <<- TRUE
        }
      }
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
  )

  #### cultivars table ####
  simul_params_cultivars(simul_params@Cultivars)
  cultivarsTable <- editableDTServer(
    id = "cultivars",
    DTdata = shiny::reactive(simul_params_cultivars()),
    disableCol = shiny::reactive({
      if (isTRUE(advanced_mode())) {
        c("reproduction_rate","death_rate")
      } else {
        names(simul_params_cultivars())
      }
    }),
    canRm = advanced_mode,
    rownames = FALSE,
    tooltips = CULTIVARS_TOOLTIP,
    row.default = default_cultivar,
    row.inc = c(1),
    col.hidden = which( names(simul_params_cultivars()) %in% c("reproduction_rate","death_rate")) -1
  )

  ##### cultivars table modification #####
  shiny::observeEvent(cultivarsTable$data,
    {
      message("Cultivars update")

      if (sum(is.na(cultivarsTable$value))) {
        return(1)
      }
      
      #message("data ", cultivarsTable$data)
      #message("value ", cultivarsTable$value)
      #message("i ", cultivarsTable$row)
      #message("j", cultivarsTable$col)

      if (isTRUE(advanced_mode())) {
        if (#nrow(cultivarsTable$data) == 0 ||
             checkCultivarsTable(cultivarsTable$data[, -which(colnames(cultivarsTable$data) %in% c("delete"))]) == FALSE) {
          can_run_simul$cultivars <<- FALSE
        }
        else {
          # Change croptypes and genes cultivars names
          # here croptypes can be invalid we use simul_params_croptypes
          # cultivarsGenes is always a reference
          # rm cultivars
          if (cultivarsTable$col == 0 && nrow(simul_params_cultivars()) > nrow(cultivarsTable$data)) {
            shiny::isolate(simul_params_croptypes(simul_params_croptypes()[, -which(colnames(simul_params_croptypes()) == cultivarsTable$value[, 1])]))
            simul_params@CultivarsGenes <<- simul_params@CultivarsGenes[-c(cultivarsTable$row), , drop = FALSE]
          }
          # add cultivars
          if (nrow(simul_params_cultivars()) < nrow(cultivarsTable$data)) {
            shiny::isolate(simul_params_croptypes(cbind(simul_params_croptypes(), rep(0, nrow(simul_params_croptypes())))))
            simul_params@CultivarsGenes <<- rbind(simul_params@CultivarsGenes, rep(0, ncol(simul_params@CultivarsGenes)))
          }
          # rename a cultivars in croptypes
          crop <- simul_params_croptypes()
          colnames(crop) <- c(colnames(simul_params_croptypes())[1:2], cultivarsTable$data[, 1])
          if(nrow(cultivarsTable$data) != 0) simul_params <<- setCroptypes(simul_params,crop)
          simul_params_croptypes(crop)
          
          colnames(simul_params@CultivarsGenes) <<- genesTable$data[, 1]
          # if (nrow(simul_params@CultivarsGenes) != 0)
            rownames(simul_params@CultivarsGenes) <<- c(cultivarsTable$data[, "cultivarName"])
            printVerbose(paste0("update CultivarsGenes",simul_params@CultivarsGenes))
          simul_params_cultivarsgenes(simul_params@CultivarsGenes)

          # update cultivars
          simul_params <<- setCultivars(simul_params, cultivarsTable$data[, -ncol(cultivarsTable$data)])
          shiny::isolate(simul_params_cultivars(simul_params@Cultivars))

          can_run_simul$cultivars <<- TRUE
        }
      }
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
  )

  #### cultivars genes table ####
  simul_params_cultivarsgenes(simul_params@CultivarsGenes)
  cultivars_genesTable <- editableDTServer(
    id = "cultivarsgenes",
    DTdata = shiny::reactive(simul_params_cultivarsgenes()),
    disableCol = shiny::reactive(c()),
    canRm = shiny::reactive({
      FALSE
    }),
    rownames = TRUE
  )

  ##### cultivars genes table modification #####
  shiny::observeEvent(cultivars_genesTable$data,
    {
      message("Cultivars genes update")

      if (sum(is.na(cultivars_genesTable$value))) {
        return(1)
      }

      # message("data ", cultivars_genesTable$data)
      # message("value ", cultivars_genesTable$value)
      # message("i ", cultivars_genesTable$row)
      # message("j", cultivars_genesTable$col)

      #if (isTRUE(advanced_mode())) {
        if (checkCultivarsGenesTable(cultivars_genesTable$data) == FALSE) {
          can_run_simul$cultivarsgenes <<- FALSE
        }
        else {
          simul_params@CultivarsGenes <<- cultivars_genesTable$data
          #print(simul_params@CultivarsGenes)
          simul_params_cultivarsgenes(simul_params@CultivarsGenes)
          can_run_simul$cultivarsgenes <<- TRUE
        }
      #}
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )

  #### genes table ####
  simul_params_genes(simul_params@Genes)
  genesTable <- editableDTServer(
    id = "genes",
    DTdata = shiny::reactive(simul_params_genes()),
    disableCol = shiny::reactive({
      if (isTRUE(advanced_mode())) {
        c()
      } else {
        names(simul_params_genes())
      }
    }),
    canRm = advanced_mode,
    rownames = FALSE,
    tooltips = GENES_TOOLTIP,
    row.default = default_gene,
    row.inc = c(1)
  )

  ##### Genes table modification #####
  shiny::observeEvent(genesTable$data,
    {
      message("Genes update")

      if (sum(is.na(genesTable$value))) {
        return(1)
      }

      # message("data ", genesTable$data)
      # message("value ", genesTable$value)
      # message("i ", genesTable$row)
      # message("j", genesTable$col)

      if (isTRUE(advanced_mode())) {
        if (#nrow(genesTable$data) == 0  ||
            checkGenesTable(genesTable$data[, -which(colnames(genesTable$data) %in% c("delete"))]) == FALSE) {
          can_run_simul$genes <<- FALSE
        }
        else {
          # rename genes in cultivars genes table
          # remove line -> remove genes in cultivars genes
          if (genesTable$col == 0 && nrow(simul_params@Genes) > nrow(genesTable$data)) {
            #print("remove here")
            simul_params@CultivarsGenes <<- simul_params@CultivarsGenes[, -c(genesTable$row), drop = FALSE]
            printVerbose(paste0("set Cultivars Genes ",simul_params@CultivarsGenes))
          }
          # add line -> add a genes in cultivars genes
          if (nrow(simul_params@Genes) < nrow(genesTable$data)) {
            #print("add here")
            simul_params@CultivarsGenes <<- cbind(simul_params@CultivarsGenes, rep(0, nrow(simul_params@CultivarsGenes)))
          }
          colnames(simul_params@CultivarsGenes) <<- genesTable$data[, 1]
          simul_params_cultivarsgenes(simul_params@CultivarsGenes)

          simul_params <<- setGenes(simul_params, genesTable$data[, -ncol(genesTable$data)])
          simul_params_genes(simul_params@Genes)
          can_run_simul$genes <<- TRUE
        }
      }
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )

  # More parameters tab
  shiny::updateNumericInput(session, "nYear", value = 10)
  shiny::updateNumericInput(session, "nTSpY", value = 120)
  simul_params <<- setTime(simul_params, Nyears = 10, nTSpY = 120)
  shiny::updateNumericInput(session, "seed", value = 1)
  simul_params <<- setSeed(simul_params, 1)

  ## Patho tabs default
  shinyjs::disable(id = "patho_infection_rate")
  shinyjs::disable(id = "patho_propagule_prod_rate")
  shinyjs::disable(id = "patho_latent_period_exp")
  shinyjs::disable(id = "patho_latent_period_var")
  shinyjs::disable(id = "patho_infectious_period_exp")
  shinyjs::disable(id = "patho_infectious_period_var")
  shinyjs::disable(id = "patho_survival_prob")
  shinyjs::disable(id = "patho_repro_sex_prob")
  shinyjs::disable(id = "patho_sigmoid_kappa")
  shinyjs::disable(id = "patho_sigmoid_sigma")
  shinyjs::disable(id = "patho_sigmoid_plateau")

  # Remove image
  output$landscapeimg <- renderPlot({
    plot(loadLandscape(input$landscape))
  })
  # output$landscapeimg <- NULL


  ############################################################
  # Screen split
  # Layout buttons split screen (left, middle, right)
  observeEvent(input$showOutputside, {
    shinyjs::showElement(id = "outputside")
    shinyjs::hideElement(id = "inputside")
    removeCssClass("inputside", "col-sm-12")
    removeCssClass("inputside", "col-sm-6")
    addCssClass("inputside", "col-sm-0")
    removeCssClass("outputside", "col-sm-0")
    removeCssClass("outputside", "col-sm-6")
    addCssClass("outputside", "col-sm-12")
  })
  observeEvent(input$showBothside, {
    removeCssClass("inputside", "col-sm-12")
    removeCssClass("inputside", "col-sm-0")
    removeCssClass("outputside", "col-sm-12")
    removeCssClass("outputside", "col-sm-0")
    addCssClass("outputside", "col-sm-6")
    addCssClass("inputside", "col-sm-6")
    shinyjs::showElement(id = "outputside")
    shinyjs::showElement(id = "inputside")
  })
  observeEvent(input$showInputside, {
    removeCssClass("inputside", "col-sm-6")
    removeCssClass("inputside", "col-sm-0")
    removeCssClass("outputside", "col-sm-12")
    removeCssClass("outputside", "col-sm-6")
    addCssClass("outputside", "col-sm-0")
    addCssClass("inputside", "col-sm-12")
    shinyjs::showElement(id = "inputside")
    shinyjs::hideElement(id = "outputside")
  })

  observeEvent(input[["inputtabpanel"]], {
    if (input[["inputtabpanel"]] == "Cultivars and Genes") {
      removeCssClass("inputside", "col-sm-6")
      removeCssClass("inputside", "col-sm-0")
      removeCssClass("outputside", "col-sm-12")
      removeCssClass("outputside", "col-sm-6")
      addCssClass("outputside", "col-sm-0")
      addCssClass("inputside", "col-sm-12")
      shinyjs::showElement(id = "inputside")
      shinyjs::hideElement(id = "outputside")
    }
  })
  # })
}
