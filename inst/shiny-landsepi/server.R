simul_params <- createSimulParams(outputDir = paste0(ROOT_PATH, "/www/tmp/"))

## fake time
simul_params <- setTime(simul_params, Nyears = 10, nTSpY = 120)
## Pathogen parameters
simul_params <- setPathogen(simul_params, loadPathogen("rust"))
simul_pathogen("rust")
## Initial conditions
# simul_params <- setInoculum(simul_params, 5e-4)
## Outputs
outputs <- loadOutputs(epid_outputs = "audpc_rel"
                       , evol_outputs = c("durability", "evol_aggr")
                       , disease="rust")
# simul_params <- setOutputs(simul_params, list(
#   epid_outputs = "audpc_rel", evol_outputs = c("durability", "evol_aggr"),
#   thres_breakdown = 50000,
#   GLAnoDis = 1.48315,
#   audpc100S = 0.76
# ))
simul_params <- setOutputs(simul_params, outputs)

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
    patho_latent_period_mean = TRUE,
    patho_latent_period_var = TRUE,
    patho_infectious_period_mean = TRUE,
    patho_infectious_period_var = TRUE,
    patho_sigmoid_kappa = TRUE,
    patho_sigmoid_sigma = TRUE,
    patho_sigmoid_plateau = TRUE,
    patho_sex_propagule_viability_limit = TRUE,
    patho_sex_propagule_release_mean = TRUE,
    inoculum = TRUE,
    treatment = TRUE
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
      can_run_simul$patho_latent_period_mean &&
      can_run_simul$patho_latent_period_var &&
      can_run_simul$patho_infectious_period_mean &&
      can_run_simul$patho_infectious_period_var &&
      can_run_simul$patho_sigmoid_kappa &&
      can_run_simul$patho_sigmoid_sigma &&
      can_run_simul$patho_sigmoid_plateau &&
      can_run_simul$patho_sex_propagule_viability_limit &&
      can_run_simul$patho_sex_propagule_release_mean &&
      can_run_simul$inoculum &&
      can_run_simul$treatment
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

    if (input$rotationcheck && !is.na(input$rotationPeriod) && input$rotationPeriod > 0) {
      sum_prop <-
        ((croptypes_proportions()[1] + croptypes_proportions()[2]) + (croptypes_proportions()[1] + croptypes_proportions()[3])) / 2
      message = paste0("The sum of the proportions of rotation croptypes must be equal to 1 (100%) <br/> example: ",
                       croptypesTable$data[1, 2],"=0.5 , ",croptypesTable$data[2, 2], "=0.5, ",croptypesTable$data[3, 2],"=0.5")
    }
    else {
      sum_prop <- sum(as.numeric(croptypes_proportions()))
      message = "The sum of the proportions of all croptypes must be equal to 1 (100%)"
    }

    shiny::removeUI(selector = "#propError")
    if (!isTRUE(all.equal(sum_prop, 1)) ||
      is.na(sum_prop)) {
      showErrorMessage(
        id = "propError", selectorafter = "#generateLandscape",
        message = message
      )
      return(invisible(FALSE))
    }
    return(invisible(TRUE))
  }

  ## Print Rotation labels
  setRotationText <- function(list_name = NULL) {
    text <- paste0("<u>1st configuration</u> : croptypes 0 (<b>", list_name[1], "</b>) and 1 (<b>", list_name[2], "</b>) , Proportion sum = ",(croptypes_proportions()[1] + croptypes_proportions()[2]), " (have to be 1)")
    text <- paste0(text, "<br/><u>2nd configuration</u> : croptypes 0 (<b>", list_name[1], "</b>) and 2 (<b>", list_name[3], "</b>) , Proportion sum = ",(croptypes_proportions()[1] + croptypes_proportions()[3]), " (have to be 1)")
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
      div(HTML(aboutText))
    ))
  })

  # ## User Mode button switch between mode
  # observeEvent(input$Mode, {
  #   advanced_mode(!advanced_mode())
  #   if (advanced_mode()) {
  #     printVerbose("enable mode edition", level = 3)
  #     removeCssClass("Mode", "btn-default")
  #     updateTabsetPanel(session, inputId = "inputtabpanel", selected = "Cultivars and Genes")
  #     shinyjs::disable(id = "demo")
  #     shinyjs::disable(id = "rotationPeriod")
  #     shiny::updateNumericInput(session, "rotationPeriod", value = 0)
  #     shinyjs::enable(id = "patho_infection_rate")
  #     shinyjs::enable(id = "patho_propagule_prod_rate")
  #     shinyjs::enable(id = "patho_latent_period_mean")
  #     shinyjs::enable(id = "patho_latent_period_var")
  #     shinyjs::enable(id = "patho_infectious_period_mean")
  #     shinyjs::enable(id = "patho_infectious_period_var")
  #     shinyjs::enable(id = "patho_survival_prob")
  #     shinyjs::enable(id = "patho_repro_sex_prob")
  #     shinyjs::enable(id = "patho_sigmoid_kappa")
  #     shinyjs::enable(id = "patho_sigmoid_sigma")
  #     shinyjs::enable(id = "patho_sigmoid_plateau")
  #     shinyjs::enable(id = "patho_sex_propagule_viability_limit")
  #     shinyjs::enable(id = "patho_sex_propagule_release_mean")
  #   }
  #   else {
  #     printVerbose("disable mode edition", level = 3)
  #     addCssClass("Mode", "btn-default")
  #     shinyjs::disable(id = "rotationPeriod")
  #     shiny::updateNumericInput(session, "rotationPeriod", value = 0)
  #     shinyjs::enable(id = "demo")
  #     shinyjs::disable(id = "patho_infection_rate")
  #     shinyjs::disable(id = "patho_propagule_prod_rate")
  #     shinyjs::disable(id = "patho_latent_period_mean")
  #     shinyjs::disable(id = "patho_latent_period_var")
  #     shinyjs::disable(id = "patho_infectious_period_mean")
  #     shinyjs::disable(id = "patho_infectious_period_var")
  #     shinyjs::disable(id = "patho_survival_prob")
  #     shinyjs::disable(id = "patho_repro_sex_prob")
  #     shinyjs::disable(id = "patho_sigmoid_kappa")
  #     shinyjs::disable(id = "patho_sigmoid_sigma")
  #     shinyjs::disable(id = "patho_sigmoid_plateau")
  #     shinyjs::disable(id = "patho_sex_propagule_viability_limit")
  #     shinyjs::disable(id = "patho_sex_propagule_release_mean")
  #   }
  # })

  # Inputs / Outputs
  ######################################################################################
  # Landscape
  shiny::observeEvent(input$landscape, {
    can_gen_landscape$proportions <<- ProportionValidation()
    can_run_simul$landscape <<- FALSE

    #shinyjs::show(id = "landscapeimg")
    output$landscapeimg <- renderPlot({
      plot(loadLandscape(input$landscape))
    })
  })
  shiny::observeEvent(input$aggregLevel, {
    can_gen_landscape$proportions <<- ProportionValidation()
    can_run_simul$landscape <<- FALSE
  })

  ######################################################################################
  # Rotation checkbox
  # shiny::observeEvent(input$rotationcheck, {
  #   if(input$rotationcheck == 1) {
  #     shinyjs::enable("rotationPeriod")
  #     shinyjs::enable("rotationText")
  #     can_gen_landscape$proportions <<- ProportionValidation()
  #   } else {
  #     shinyjs::disable("rotationPeriod")
  #     shinyjs::disable("rotationText")
  #     can_gen_landscape$proportions <<- ProportionValidation()
  #   }
  # 
  # })
  
  ######################################################################################
  # Rotation period validation
  shiny::observeEvent(input$rotationPeriod, {
    can_gen_landscape$rotation <<- TRUE
    can_run_simul$landscape <<- FALSE
    shiny::removeUI(selector = "#rotationPeriodError")
    if (simul_demo() == "RO" && advanced_mode() == FALSE) {
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
      simul_params <<- setTime(simul_params, Nyears = as.numeric(input$nYear), nTSpY = as.numeric(input$nTSpY))
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
      simul_params <<- setTime(simul_params, Nyears = as.numeric(input$nYear), nTSpY = as.numeric(input$nTSpY))
      can_run_simul$nTSpY <<- TRUE
      # if (input$patho_repro_sex_active == TRUE) {
      #   simul_params <<- updateReproSexProb(simul_params, c(rep(0, simul_params@TimeParam$nTSpY), 1))
      # } else {
        simul_params <<- updateReproSexProb(simul_params, rep(0, simul_params@TimeParam$nTSpY + 1))
      # }
    }
  }, ignoreNULL = FALSE)
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
  shiny::observeEvent(input$defaultPathogen, {simul_pathogen(tolower(input$defaultPathogen))})
  
  # Update pathogen and pathosystem parameters
  update_pathogen <- function(){
    simul_params <<- setPathogen(simul_params, loadPathogen(disease = simul_pathogen()))
    simul_params <<- setOutputs(simul_params, loadOutputs(epid_outputs = "audpc_rel"
                                                          , evol_outputs = c("durability", "evol_aggr")
                                                          , disease = simul_pathogen()))
    #simul_pathogen(tolower(input$defaultPathogen))
    #printVerbose(tolower(input$defaultPathogen))
    # if (advanced_mode() == FALSE) {
    #   update_demo()
    # }
    updateNumericInput(session = session, inputId = "patho_survival_prob", value = simul_params@Pathogen$survival_prob)
    updateNumericInput(session = session, inputId = "patho_repro_sex_prob", value = simul_params@Pathogen$repro_sex_prob)
    updateNumericInput(session = session, inputId = "patho_infection_rate", value = simul_params@Pathogen$infection_rate)
    updateNumericInput(session = session, inputId = "patho_propagule_prod_rate", value = simul_params@Pathogen$propagule_prod_rate)
    updateNumericInput(session = session, inputId = "patho_latent_period_mean", value = simul_params@Pathogen$latent_period_mean)
    #updateNumericInput(session = session, inputId = "patho_latent_period_var", value = simul_params@Pathogen$latent_period_var)
    updateNumericInput(session = session, inputId = "patho_infectious_period_mean", value = simul_params@Pathogen$infectious_period_mean)
    #updateNumericInput(session = session, inputId = "patho_infectious_period_var", value = simul_params@Pathogen$infectious_period_var)
    #updateNumericInput(session = session, inputId = "patho_sigmoid_kappa", value = simul_params@Pathogen$sigmoid_kappa)
    #updateNumericInput(session = session, inputId = "patho_sigmoid_sigma", value = simul_params@Pathogen$sigmoid_sigma)
    #updateNumericInput(session = session, inputId = "patho_sigmoid_plateau", value = simul_params@Pathogen$sigmoid_plateau)
    #updateCheckboxInput(session = session, inputId = "patho_repro_sex_active", value = FALSE)
    #updateNumericInput(session = session, inputId = "patho_sex_propagule_viability_limit", value = simul_params@Pathogen$sex_propagule_viability_limit)
    #updateNumericInput(session = session, inputId = "patho_sex_propagule_release_mean", value = simul_params@Pathogen$sex_propagule_release_mean)
    hostName = disease2CultivarType(simul_pathogen())
    cvalues <- c("crop", "nonCrop")
    names(cvalues) <- c(hostName, "nonCrop (city, forest, etc.)")
    updateSelectInput(session, inputId = "listcultivarstype", choices=cvalues)
    
    options <- pathosystemParams(input$defaultPathogen)
    shiny::updateNumericInput(session, inputId = "nTSpY", value = options$nTSpY)
  }

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
      #simul_params <<- setInoculum(simul_params, input$inoculum)
      # Have to be set when all landscape, cultivars and pathogen have been set
      can_run_simul$inoculum <<- TRUE
    }
  })

  shiny::observeEvent(input$inoculumstrategy, {
    if (input$inoculumstrategy == 0) {
      message("PI0 Local")
      shinyjs::show("inoculum_message")
      # showErrorMessage(
      #     id = "inoculumstrategyWarn", selectorafter = "#generateLandscape",
      #     message = "Pathogen : Initial probability of infection for susceptible cultivar should be upper than 0.01"
      # )
    } else {
      message("PI0 Global default behaviour")
      #shiny::removeUI(selector = "#inoculumstrategyWarn")
      shinyjs::hide("inoculum_message")
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
        message = paste0("The maximal expected effective propagule production rate of an infectious host per time step should be between 0 and ", VALUEMAX)
      )
      can_run_simul$patho_propagule_prod_rate <<- FALSE
    }
    else {
      simul_params@Pathogen$propagule_prod_rate <<- input$patho_propagule_prod_rate
      can_run_simul$patho_propagule_prod_rate <<- TRUE
    }
  })

  # latent_period_mean
  shiny::observeEvent(input$patho_latent_period_mean, {
    shiny::removeUI(selector = "#pathoLatPerExpError")
    if (input$patho_latent_period_mean > VALUEMAX || input$patho_latent_period_mean < 0 || is.na(input$patho_latent_period_mean)) {
      showErrorMessage(
        id = "pathoLatPerExpError", selectorafter = "#generateLandscape",
        message = paste0("The minimal expected duration of the latent period should be between 0 and ", VALUEMAX)
      )
      can_run_simul$patho_latent_period_mean <<- FALSE
    }
    else {
      simul_params@Pathogen$latent_period_mean <<- input$patho_latent_period_mean
      can_run_simul$patho_latent_period_mean <<- TRUE
    }
  })

  # latent_period_var
  # shiny::observeEvent(input$patho_latent_period_var, {
  #   shiny::removeUI(selector = "#pathoLatPerVarError")
  #   if (input$patho_latent_period_var > VALUEMAX || input$patho_latent_period_var < 0 || is.na(input$patho_latent_period_var)) {
  #     showErrorMessage(
  #       id = "pathoLatPerVarError", selectorafter = "#generateLandscape",
  #       message = paste0("The variance of the infectious period duration should be between 0 and ", VALUEMAX)
  #     )
  #     can_run_simul$patho_latent_period_var <<- FALSE
  #   }
  #   else {
  #     simul_params@Pathogen$latent_period_var <<- input$patho_latent_period_var
  #     can_run_simul$patho_latent_period_var <<- TRUE
  #   }
  # })

  # infectious_period_mean
  shiny::observeEvent(input$patho_infectious_period_mean, {
    shiny::removeUI(selector = "#pathoInfPerExpError")
    if (input$patho_infectious_period_mean > VALUEMAX || input$patho_infectious_period_mean < 0 || is.na(input$patho_infectious_period_mean)) {
      showErrorMessage(
        id = "pathoInfPerExpError", selectorafter = "#generateLandscape",
        message = paste0("The maximal expected duration of the infectious period should be between 0 and ", VALUEMAX)
      )
      can_run_simul$patho_infectious_period_mean <<- FALSE
    }
    else {
      simul_params@Pathogen$infectious_period_mean <<- input$patho_infectious_period_mean
      can_run_simul$patho_infectious_period_mean <<- TRUE
    }
  })

  # infectious_period_var
  # shiny::observeEvent(input$patho_infectious_period_var, {
  #   shiny::removeUI(selector = "#pathoInfPerVarError")
  #   if (input$patho_infectious_period_var > VALUEMAX || input$patho_infectious_period_var < 0 || is.na(input$patho_infectious_period_var)) {
  #     showErrorMessage(
  #       id = "pathoInfPerVarError", selectorafter = "#generateLandscape",
  #       message = paste0("The variance of the infectious period duration should be between 0 and ", VALUEMAX)
  #     )
  #     can_run_simul$patho_infectious_period_var <<- FALSE
  #   }
  #   else {
  #     simul_params@Pathogen$infectious_period_var <<- input$patho_infectious_period_var
  #     can_run_simul$patho_infectious_period_var <<- TRUE
  #   }
  # })

  # sigmoid_kappa
  # shiny::observeEvent(input$patho_sigmoid_kappa, {
  #   shiny::removeUI(selector = "#pathoSigKapError")
  #   if (input$patho_sigmoid_kappa > VALUEMAX || input$patho_sigmoid_kappa < 0 || is.na(input$patho_sigmoid_kappa)) {
  #     showErrorMessage(
  #       id = "pathoSigKapError", selectorafter = "#generateLandscape",
  #       message = paste0("The kappa parameter of the sigmoid contamination function should be between 0 and ", VALUEMAX)
  #     )
  #     can_run_simul$patho_sigmoid_kappa <<- FALSE
  #   }
  #   else {
  #     simul_params@Pathogen$sigmoid_kappa <<- input$patho_sigmoid_kappa
  #     can_run_simul$patho_sigmoid_kappa <<- TRUE
  #   }
  # })
  # 
  # # sigmoid_sigma
  # shiny::observeEvent(input$patho_sigmoid_sigma, {
  #   shiny::removeUI(selector = "#pathoSigSigError")
  #   if (input$patho_sigmoid_sigma > VALUEMAX || input$patho_sigmoid_sigma < 0 || is.na(input$patho_sigmoid_sigma)) {
  #     showErrorMessage(
  #       id = "pathoSigSigError", selectorafter = "#generateLandscape",
  #       message = paste0("The sigma parameter of the sigmoid contamination function should be between 0 and ", VALUEMAX)
  #     )
  #     can_run_simul$patho_sigmoid_sigma <<- FALSE
  #   }
  #   else {
  #     simul_params@Pathogen$sigmoid_sigma <<- input$patho_sigmoid_sigma
  #     can_run_simul$patho_sigmoid_sigma <<- TRUE
  #   }
  # })

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

  # # sexual propagules viability limit
  # shiny::observeEvent(input$patho_repro_sex_active, {
  #   if (input$patho_repro_sex_active == TRUE) {
  #     # force check
  #     updateNumericInput(session = session, inputId = "patho_sex_propagule_viability_limit", value = simul_params@Pathogen$sex_propagule_viability_limit)
  #     updateNumericInput(session = session, inputId = "patho_sex_propagule_release_mean", value = simul_params@Pathogen$sex_propagule_release_mean)
  #     simul_params <<- updateReproSexProb(simul_params, c(rep(0, simul_params@TimeParam$nTSpY), 1))
  #   }
  #   else {
  #     shiny::removeUI(selector = "#pathoDorLimError")
  #     shiny::removeUI(selector = "#pathoMuExpError")
  #     simul_params <<- updateReproSexProb(simul_params, rep(0, simul_params@TimeParam$nTSpY + 1))
  #   }
  # })
  # 
  # # sexual propagules viability limit
  # shiny::observeEvent(input$patho_sex_propagule_viability_limit, {
  #   shiny::removeUI(selector = "#pathoDorLimError")
  #   if (input$patho_repro_sex_active == TRUE) {
  #     if (input$patho_sex_propagule_viability_limit < 0 || input$patho_sex_propagule_viability_limit > simul_params@TimeParam$nTSpY) {
  #       showErrorMessage(
  #         id = "pathoDorLimError", selectorafter = "#generateLandscape",
  #         message = "Wrong value for pathogen sexual propagules viability limit"
  #       )
  #       can_run_simul$path_sex_propagule_viability_limit <<- FALSE
  #     }
  #     else {
  #       simul_params@Pathogen$sex_propagule_viability_limit <<- input$patho_sex_propagule_viability_limit
  #       can_run_simul$patho_sex_propagule_viability_limit <<- TRUE
  #     }
  #   }
  # })

  # sexual propagules mu exp average release
  # shiny::observeEvent(input$patho_sex_propagule_release_mean, {
  #   shiny::removeUI(selector = "#pathoMuExpError")
  #   if (input$patho_repro_sex_active == TRUE) {
  #     if (input$patho_sex_propagule_release_mean < 1) {
  #       showErrorMessage(
  #         id = "pathoMuExpError", selectorafter = "#generateLandscape",
  #         message = "Pathogen sexual propagules average number of cropping seasons value have to be > 0"
  #       )
  #       can_run_simul$patho_sex_propagule_release_mean <<- FALSE
  #     }
  #     else {
  #       simul_params@Pathogen$patho_sex_propagule_release_mean <<- input$patho_sex_propagule_release_mean
  #       can_run_simul$patho_sex_propagule_release_mean <<- TRUE
  #     }
  #   }
  # })

  ######################################################################################
  #
  # Treatment Tab Observe
  #
  ######################################################################################
  updateTreatment <- function() {
    shiny::removeUI(selector = "#treatmentError")
    
    if (treatment_is_active()) {
      # if (length(input$treatment_cultivars_select) < 1 || 
      #   input$treatment_days_interval > input$nTSpY ||
      #    input$treatment_efficiency < 0 || input$treatment_efficiency > 1
      #   ) {
      #   showErrorMessage(
      #     id = "treatmentError", selectorafter = "#generateLandscape",
      #     message = "Trouble in treatment values: at least one cultivar must be selected, interval between two treatments must be < 100, and efficiency must be between 0 and 1"
      #   )
      # }
      if (length(input$treatment_cultivars_select) < 1) {
        showErrorMessage(id = "treatmentError", selectorafter = "#generateLandscape",
                         message = "Trouble in treatment values: at least one cultivar must be selected")
        can_run_simul$treatment <<- FALSE
      }else if (input$treatment_days_interval > input$nTSpY) {
        showErrorMessage(id = "treatmentError", selectorafter = "#generateLandscape",
                         message = "Trouble in treatment values: interval between two treatments must be < 100 days")
        can_run_simul$treatment <<- FALSE
      } else if (input$treatment_efficiency < 0 || input$treatment_efficiency > 1) {
        showErrorMessage(id = "treatmentError", selectorafter = "#generateLandscape",
                         message = "Trouble in treatment values: efficiency must be between 0 and 1")
        can_run_simul$treatment <<- FALSE
      }
      else {
        days_list <- seq(1, as.numeric(input$nTSpY), input$treatment_days_interval)
        cults <- which(simul_params_cultivars()[, 1] %in% input$treatment_cultivars_select)
        cults <- cults -1 #id cultivars start at 0
        treatment <- loadTreatment(tolower(input$defaultPathogen))
        treatment$treatment_efficiency = input$treatment_efficiency
        treatment$treatment_timesteps = days_list
        treatment$treatment_cultivars = cults
        treatment$treatment_application_threshold = rep(0.0,length(cults))
        simul_params <<- setTreatment(simul_params, treatment)
        can_run_simul$treatment <<- TRUE
      }
    }
    else {
      simul_params <<- setTreatment(simul_params, loadTreatment())
      can_run_simul$treatment <<- TRUE
    }
  }

  # Active Treatment
  shiny::observeEvent(input$treatment_active, {
    treatment_is_active(input$treatment_active)
    if (treatment_is_active()) {
      shinyjs::enable(id = "treatment_days_interval")
      # shinyjs::enable(id = "treatment_day_start")
      shinyjs::enable(id = "treatment_cultivars_select")
      shinyjs::enable(id = "treatment_efficiency")
      # shinyjs::enable(id = "treatment_degradation_rate")
      #shinyjs::enable(id = "treatment_cost")
      updateSelectInput(session, "treatment_cultivars_select", choices = c(Choose='',simul_params_cultivars()[, 1]), selected=NULL)
    }
    else {
      shinyjs::disable(id = "treatment_days_interval")
      # shinyjs::disable(id = "treatment_day_start")
      shinyjs::disable(id = "treatment_cultivars_select")
      shinyjs::disable(id = "treatment_efficiency")
      # shinyjs::disable(id = "treatment_degradation_rate")
      #shinyjs::disable(id = "treatment_cost")
    }
    updateTreatment()
  })

  # # day start Treatment
  # shiny::observeEvent(input$treatment_day_start, {
  #   updateTreatment()
  # })
  # days between Treatment
  shiny::observeEvent(input$treatment_days_interval, {
    updateTreatment()
  })

  # Cultivars Treatment
  shiny::observeEvent(input$treatment_cultivars_select, {
    #printVerbose(input$treatment_cultivars_select)
    updateTreatment()
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  # # beta Treatment
  # shiny::observeEvent(input$treatment_degradation_rate, {
  #   updateTreatment()
  # })

  # trait red Treatment
  shiny::observeEvent(input$treatment_efficiency, {
    updateTreatment()
  })
  
  # shiny::observeEvent(input$treatment_cost, {
  #     updateTreatment()
  # })

  ######################################################################
  # Handle the download gpkg button
  ######################################################################
  # output$export <-
  #   # shiny::downloadHandler(
  #   #   filename = "landsepi_landscape.gpkg",
  #   #   content <- function(file) {
  #   #     simul_params <<- saveDeploymentStrategy(simul_params)
  #   #     file.copy(file.path(simul_params@OutputDir, simul_params@OutputGPKG), file)
  #   #   },
  #   #   contentType = "application/x-sqlite3"
  #   # )
  #   shiny::downloadHandler(
  #     filename = "landsepi_landscape.zip",
  #     content <- function(file) {
  #       simul_params <<- saveDeploymentStrategy(simul_params)
  #       filels <- file.path(simul_params@OutputDir, simul_params@OutputGPKG)
  #       filetxt <- list.files(simul_params@OutputDir, pattern = "*.txt", full.names = TRUE)
  #       filels <- c(filels, filetxt)
  #       zip(zipfile = file, files = filels, extras = "-j")
  #     },
  #     contentType = "application/zip"
  #   )
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

      # print(simul_params@Croptypes)
      # print(simul_params@Cultivars)
      # print(simul_params@CultivarsGenes)
      # print(simul_params@Genes)

      # Croptypes Rotation
      if (input$rotationcheck && input$rotationPeriod > 0) {
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
        prop <- list(croptypes_proportions())
      }
      
      ## TODO check PY
      # if (simul_demo() == "PY") {
      #     prop <- list(croptypes_proportions()[1:2])
      #   } else {
      #    prop <- list(croptypes_proportions())
      # }

      simul_params <<- setSeed(simul_params, input$seed)

      incProgress(0.4)
      # Run the landscape generation
      simul_params <<- setLandscape(simul_params, loadLandscape(input$landscape))
      ## Dispersal parameters
      simul_params <<- setDispersalPathogen(simul_params, loadDispersalPathogen(input$landscape)[[1]])
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

      #shinyjs::show(id = "landscapeimg")
      output$landscapeimg <- renderPlot({
        imgs <- normalizePath(list.files(simul_params@OutputDir, pattern = ".png", full.names = TRUE))
        pngs <- lapply(imgs, readPNG)
        asGrobs <- lapply(pngs, rasterGrob)
        p <- grid.arrange(grobs = asGrobs, nrow = 1)
      })
      shiny::updateTabsetPanel(session, "outputtabpanel", selected = "landscapeTab")

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
    printVerbose(simul_params, level = 2)
    
    withProgress(message = "Running, please wait...", value = 0, {

      printVerbose("set inoculm",level = 1)
      if(input$inoculumstrategy == 0){
        pI0 <- inoculumRandomLocal(simul_params,input$inoculum,1)
        #print(inoculumToMatrix(pI0))
        printVerbose(pI0, level=3) 
      }
      else { 
        pI0 <- input$inoculum
        printVerbose(pI0, level=3) 
      }
      simul_params <<- setInoculum(simul_params, pI0)
      

      
      progressBar <- Progress$new()
      progressBar$set(value = NULL, message = "Running Simulation, please wait...")
      # setwd(paste0(ROOT_PATH,"/www/tmp/"))

      shinyjs::disable(id = "generateLandscape")
      shinyjs::disable(id = "runSimulation")
      shinyjs::disable(id = "export")
      shinyjs::enable(id = "stopSimulation")

      shinyjs::disable("showInputside")
      shinyjs::disable("showBothside")
      # shinyjs::click("showOutputside") # seems not working -> force it
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

      plan(list(multisession, multicore))

      future_process <<- future({
        res <- landsepi::runSimul(simul_params,
          graphic = TRUE, videoMP4 = TRUE
        )
      }, seed=input$seed)

      then(future_process,
        onFulfilled = function(value) {
          progressBar$set(value = 0.8, message = "Simulation ended : making video...")

          shinyjs::enable(id = "generateLandscape")
          # shinyjs::enable(id = "runSimulation")
          shinyjs::enable(id = "export")
          shinyjs::enable(id = "runSimulation")
          shinyjs::disable(id = "stopSimulation")

          #output$landscapeimg <- NULL
          #hide(id = "landscapeimg")

          output$video <-
            shiny::renderUI(
              tags$video(
                id = "video",
                type = "video/webm",
                src = paste0("video/", basename(simul_params@OutputDir), "/video.webm?rand=", as.integer(Sys.time())),
                controls = "controls",
                width = "100%",
                height = "auto"
              )
            )
          shiny::updateTabsetPanel(session, "outputtabpanel", selected = "videoTab")
          shinyjs::enable("showInputside")
          shinyjs::enable("showBothside")
          shinyjs::click("showOutputside")
          shinyBS::addTooltip(session, "runSimulation", title = RUN_SIMULATION, placement = "top", trigger = "hover")
          # setwd(dirname(getwd()))

          ## Table and graphics tab
          audpc <- read.csv2(paste0(simul_params@OutputDir,"/audpc_rel.txt"), sep = ",")
          audpcTable <- sapply(names(audpc), function(name) {
            # round(sum(as.numeric(audpc[,name]))/nrow(audpc), 2)
            round(mean(as.numeric(audpc[,name]), na.rm=TRUE), 2)
          })
          output$audpctable <- DT::renderDT(data.frame("AUDPC" = audpcTable),
                                            options= list(dom = 't'))
          
          durability <- read.csv2(paste0(simul_params@OutputDir,"/durability.txt"), sep = ",")
          durabilityTable <- data.frame(Durability = round(as.numeric(durability / input$nTSpY), 2), Erosion = NA)
          rownames(durabilityTable) <- simul_params@Genes$geneName[ which(apply(simul_params@CultivarsGenes, MARGIN = 2, sum) != 0)]
          
          index_notBroken <- which(durabilityTable$Durability > input$nYear)
          index_broken <- which(durabilityTable$Durability <= input$nYear)
          # durabilityTable$Durability[is.na(durabilityTable$Durability)] <- ""  ## (not necessary: NA are not displayed)
          durabilityTable$Durability[index_notBroken] <- "Not broken"
          durabilityTable$Durability[index_broken] <- paste(as.character(durabilityTable$Durability[index_broken]), "years")
          
          ## Compute level of erosion
          for (g in which(apply(simul_params@CultivarsGenes, MARGIN = 2, sum) != 0)){ #1:nrow(simul_params@Genes)){ #use affected genes
            nLevels <- simul_params@Genes$Nlevels_aggressiveness[g] - 1 ## -1 because first level does not count
            if (nLevels>0){ 
              evol_aggr_g <- read.csv2(paste0(simul_params@OutputDir,"/evol_aggr_",simul_params@Genes$geneName[g],".txt"), sep = ",",header = TRUE)
              lastLevel <- rev(which(evol_aggr_g$R_invasion < input$nTSpY*input$nYear))[1] - 1  ## last level reached (-1 because first level doesn't count)
              if (is.na(lastLevel)) {
                durabilityTable$Erosion[g] <- "No erosion"
              }else{
                durabilityTable$Erosion[g] <- paste(as.character(round(lastLevel / nLevels * 100, 2)), "%")
              }
            # }else{ ## else nLevels==0
              # durabilityTable$Erosion[g] <- ""  ## (not necessary: NA are not displayed)
            } ## if nlevels >0
          } ## for g
          output$durabilitytable <- DT::renderDT(data.frame(durabilityTable),
                                                 options= list(dom = 't'))
          
          output$freqpathogenotypes <- renderPlot({
            freqPathoGenotypes_file <- normalizePath(paste0(simul_params@OutputDir,"/freqPathoGenotypes.tiff"))
            #img <- tiff::readTIFF(freqPathoGenotypes_file)
            tiffs <- lapply(freqPathoGenotypes_file, tiff::readTIFF)
            asGrobs <- lapply(tiffs, rasterGrob)
            p <- grid.arrange(grobs = asGrobs, nrow = 1)
          })
          
        },
        onRejected = function(err) {
          setwd(ROOT_PATH)
          cat(file = stderr(), "\n ### KILL'EM ALL ### -> Kill simulation \n")
          cleanDir(simul_params@OutputDir)
          shinyjs::enable("showInputside")
          shinyjs::enable("showBothside")
          shinyjs::click("showBothside")
          shinyjs::enable(id = "generateLandscape")
          shinyBS::addTooltip(session, "runSimulation", title = RUN_SIMULATION, placement = "top", trigger = "hover")
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
  simul_demo <- reactiveVal("MO")
  shiny::observeEvent(input$demo, {})
  
  shiny::observeEvent(input$load, {
      update_pathogen()
      simul_demo(input$demo)
      update_demo()
  },ignoreInit = FALSE)
  

  ## Load a demo strategy
  update_demo <- function() {
    # Cultivar tab
    switch(simul_demo(),
      MO = {
        simul_params <<- loadDemoMO(simul_params, disease = simul_pathogen())
      },
      MI = {
        simul_params <<- loadDemoMI(simul_params, disease = simul_pathogen())
      },
      RO = {
        simul_params <<- loadDemoRO(simul_params, disease = simul_pathogen())
      },
      PY = {
        simul_params <<- loadDemoPY(simul_params, disease = simul_pathogen())
      },
      {
        # Default case
        print(paste("simul_demo() : Unknown ", simul_demo()))
      }
    )

    default_gene <<- simul_params@Genes[1, ]
    default_cultivar <<- simul_params@Cultivars[1, ]
    default_croptype <<- simul_params@Croptypes[1, c(1, 2)]

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
    shinyjs::hide(id = "rotationPeriod")
    shinyjs::hide(id = "rotationText")
    shinyjs::disable(id = "rotationPeriod")
    shinyjs::disable(id = "rotationText")
    shinyjs::enable("addcrop")
    can_rm_croptypes(TRUE)
    
    shiny::updateCheckboxInput(session, "rotationcheck",value = 0)
    shiny::updateNumericInput(session, "rotationPeriod", value = 0)

    if (simul_demo() == "MO") {
      croptypes_proportions(c(0.33, 0.33, 0.34))
      shiny::updateSelectInput(session, "aggregLevel", selected = "high")
    }
    else if (simul_demo() == "MI" || simul_demo() == "PY") {
      croptypes_proportions(c(0.5, 0.5))
      shiny::updateSelectInput(session, "aggregLevel", selected = "high")
    }
    else if (simul_demo() == "RO") {
      croptypes_proportions(c(0.5, 0.5, 0.5))
      shiny::updateCheckboxInput(session, "rotationcheck",value = 1)
      shinyjs::show(id = "rotationPeriod")
      shinyjs::show(id = "rotationText")
      shinyjs::enable(id = "rotationPeriod")
      shinyjs::enable(id = "rotationText")
      shiny::updateNumericInput(session, "rotationPeriod", value = 2)
      shiny::updateSelectInput(session, "aggregLevel", selected = "medium")
      shinyjs::disable("addcrop")
      can_rm_croptypes(FALSE)
    }
    else if (simul_demo() == "PY") {
      shiny::updateSelectInput(session, "aggregLevel", selected = "low")
    }

    output$rotationText <- renderUI({
      HTML(setRotationText(simul_params@Croptypes[, 2]))
    })
  }

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
    disableCol = shiny::reactive({c("croptypeID")}),
    canRm = can_rm_croptypes,
    rownames = FALSE,
    tooltips = c("Croptype index (starts at 0)", "Croptype name"),
    row.default = default_croptype,
    row.cols = 1:2,
    row.inc = c(1, 2)
  )

  ##### croptypes table modification #####
  shiny::observeEvent(croptypesTable$data,
    {
      printVerbose("Croptypes update", level=2)

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
        if(isTRUE(can_rm_croptypes())) shiny::isolate(simul_params_croptypes(croptypesTable$data[, 1:(ncol(croptypesTable$data) - 2)]))
        else shiny::isolate(simul_params_croptypes(croptypesTable$data[, 1:(ncol(croptypesTable$data) - 1)]))
        
        if (nrow(croptypesTable$data) == 0 ||
          checkCroptypesTable(croptypesTable$data[, -which(colnames(croptypesTable$data) %in% c("Proportions", "delete"))]) == FALSE) {
          can_run_simul$croptypes <<- FALSE
          can_gen_landscape$croptypeID <<- FALSE
        }
        else {
          croptypesTable$data[, "croptypeID"] <- seq(1:nrow(croptypesTable$data)) - 1
          output$rotationText <- renderUI({
            HTML(setRotationText(croptypesTable$data[, 2]))
          })
          #simul_params <<- setCroptypes(simul_params, croptypesTable$data[, 1:(ncol(croptypesTable$data) - 2)])
          if(isTRUE(can_rm_croptypes())) simul_params <<- setCroptypes(simul_params, croptypesTable$data[, 1:(ncol(croptypesTable$data) - 2)])
          else simul_params <<- setCroptypes(simul_params, croptypesTable$data[, 1:(ncol(croptypesTable$data) - 1)])
          can_run_simul$croptypes <<- TRUE
          can_gen_landscape$croptypeID <<- TRUE
        }
      }
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
  )

  #### Croptypes Add ####
  shiny::observeEvent(input$addcrop,{
    
    printVerbose("Add Croptypes", level=2)
    index <- 1
    repeat{
      newCroptypeId <- sum(length(grep("Susceptible_crop", croptypesTable$data[,"croptypeName"])))+index
      susceptibleName <- paste0("Susceptible_crop","_",newCroptypeId)
      if( length((grep(susceptibleName, croptypesTable$data[,"croptypeName"]))) == 0) break;
      index <- index + 1
    }
    newSusceptible <- max(croptypesTable$data[,"croptypeID"])+1
    if(is.infinite(newSusceptible)) newSusceptible <- 1
    
    simul_params@Croptypes[nrow(croptypesTable$data)+1,1] <<- newSusceptible
    names(simul_params@Croptypes[nrow(croptypesTable$data)+1,]) <<- newSusceptible
    simul_params@Croptypes[nrow(croptypesTable$data)+1,2] <<- susceptibleName
    simul_params@Croptypes[nrow(croptypesTable$data)+1,-c(1,2)] <<- rep(0,ncol(croptypesTable$data)-4)

    shiny::isolate(simul_params_croptypes(simul_params@Croptypes))
    can_run_simul$croptypes <<- FALSE

    croptypes_proportions(c(croptypesTable$data[, "Proportions"],0))
    can_gen_landscape$proportions <<- ProportionValidation()
  })
  
  #### cultivars table ####
  simul_params_cultivars(simul_params@Cultivars)
  cultivarsTable <- editableDTServer(
    id = "cultivars",
    DTdata = shiny::reactive(simul_params_cultivars()),
    disableCol = shiny::reactive({
      if (isTRUE(advanced_mode())) {
        c("reproduction_rate")
      } else {
        names(simul_params_cultivars())
      }
    }),
    canRm = advanced_mode,
    rownames = FALSE,
    tooltips = CULTIVARS_TOOLTIP,
    row.default = default_cultivar,
    row.inc = c(1),
    #col.hidden = which(names(simul_params_cultivars()) %in% c("reproduction_rate")) - 1
    col.hidden = c(1,2,3,4,5,6,7,8,9,10)
  )

  ##### cultivars table modification #####
  shiny::observeEvent(cultivarsTable$data,
    {
      message("Cultivars update")

      if (sum(is.na(cultivarsTable$value))) {
        return(1)
      }

      # message("data ", cultivarsTable$data)
      # message("value ", cultivarsTable$value)
      # message("i ", cultivarsTable$row)
      # message("j", cultivarsTable$col)

      if (isTRUE(advanced_mode())) {
        if ( # nrow(cultivarsTable$data) == 0 ||
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
          if (nrow(cultivarsTable$data) != 0) simul_params <<- setCroptypes(simul_params, crop)
          simul_params_croptypes(crop)
          # rename cultivars in genes
          colnames(simul_params@CultivarsGenes) <<- genesTable$data[, 1]
          # if (nrow(simul_params@CultivarsGenes) != 0)
          rownames(simul_params@CultivarsGenes) <<- c(cultivarsTable$data[, "cultivarName"])
          printVerbose(paste0("update CultivarsGenes", simul_params@CultivarsGenes))
          simul_params_cultivarsgenes(simul_params@CultivarsGenes)
          # update cultivars
          simul_params <<- setCultivars(simul_params, cultivarsTable$data[, -ncol(cultivarsTable$data),drop=FALSE])
          shiny::isolate(simul_params_cultivars(simul_params@Cultivars))
          can_run_simul$cultivars <<- TRUE

          # update treatment cultivars list
          updateSelectInput(session, "treatment_cultivars_select", choices = c(Choose='',simul_params_cultivars()[, 1]))
        }
      }
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
  )
  
  shiny::observeEvent(input$addcultivar,{
    
    message("\n Add Cultivar")
      
    # print(simul_params@Cultivars)
    # print(cultivarsTable$data)
    # print(simul_params@CultivarsGenes)
    # print(cultivars_genesTable$data)
    # update cultivars
    cultivar <- disease2CultivarType(tolower(input$defaultPathogen), input$listcultivarstype)
    #cultivar <- input$listcultivarstype
    
    index <- 1
    repeat{
      newCultivarId <- sum(length(grep(cultivar, cultivarsTable$data[,1])))+index
      cultivarName <- paste0(cultivar,"_",newCultivarId)
      if( length((grep(cultivarName, cultivarsTable$data[,1]))) == 0) break;
      index <- index + 1
    }
    newCultivars <- loadCultivar(cultivarName, type=cultivar)
    simul_params <<- setCultivars(simul_params, rbind(simul_params@Cultivars, newCultivars))
    shiny::isolate(simul_params_cultivars(simul_params@Cultivars))
    can_run_simul$cultivars <<- TRUE
    
    # update croptypes and cultivarsGenes
    shiny::isolate(simul_params_croptypes(cbind(simul_params_croptypes(), rep(0, nrow(simul_params_croptypes())))))
    simul_params@CultivarsGenes <<- rbind(simul_params@CultivarsGenes, rep(0, ncol(simul_params@CultivarsGenes)))
    # rename a cultivars in croptypes
    crop <- simul_params_croptypes()
    #colnames(crop) <- c(colnames(simul_params_croptypes())[1:2], cultivarName)
    colnames(crop) <- c(colnames(simul_params_croptypes())[1:2], simul_params@Cultivars$cultivarName)
    if (nrow(cultivarsTable$data) != 0) simul_params <<- setCroptypes(simul_params, crop)
    simul_params_croptypes(crop)
    
    # rename cultivars in genes
    colnames(simul_params@CultivarsGenes) <<- genesTable$data[, 1]
    rownames(simul_params@CultivarsGenes) <<- c(cultivarsTable$data[, "cultivarName"],cultivarName)
    printVerbose(paste0("update CultivarsGenes ", simul_params@CultivarsGenes))
    simul_params_cultivarsgenes(simul_params@CultivarsGenes)
    # update treatment cultivars list
    updateSelectInput(session, "treatment_cultivars_select", choices = c(Choose='',simul_params_cultivars()[, 1]))
  })

  #### cultivars genes table ####
  simul_params_cultivarsgenes(simul_params@CultivarsGenes)
  cultivars_genesTable <- editableDTServer(
    id = "cultivarsgenes",
    DTdata = shiny::reactive(simul_params_cultivarsgenes()),
    #disableCol = shiny::reactive(c()),
    canRm = shiny::reactive({
      TRUE
    }),
    rownames = TRUE
  )

  ##### cultivars genes table modification #####
  shiny::observeEvent(cultivars_genesTable$data,
    {
      message("Cultivars Genes update")

      if (sum(is.na(cultivars_genesTable$value))) {
        return(1)
      }

      # message("data ", cultivars_genesTable$data)
      # message("data rownames ", rownames(cultivars_genesTable$data))
      # message("value ", cultivars_genesTable$value)
      # message("i ", cultivars_genesTable$row)
      # message("j", cultivars_genesTable$col)
      
      # rm line
      if( nrow(cultivars_genesTable$data) < nrow(simul_params@CultivarsGenes) ){
          shinyjs::runjs(paste0('Shiny.setInputValue("cultivars-deletePressed", "button_', cultivars_genesTable$row,'", {priority: "event"})'))
          #cultivarsTable$data <- cultivarsTable$data[-cultivars_genesTable$row,]
          return()
      }
      
      # rownames (cultivars)
      if( !is.na(cultivars_genesTable$col) && cultivars_genesTable$col == 0
          && nrow(cultivars_genesTable$data) == nrow(simul_params@CultivarsGenes)
          && is.character(cultivars_genesTable$value)) {
          #update cultivars table and exit
          cultivarsTable$data$cultivarName[cultivars_genesTable$row] <- cultivars_genesTable$value
          return()
      }
      
      # if (isTRUE(advanced_mode())) {
      if (checkCultivarsGenesTable(cultivars_genesTable$data) == FALSE) {
        can_run_simul$cultivarsgenes <<- FALSE
        
      }
      else {
        simul_params@CultivarsGenes <<- cultivars_genesTable$data[,-ncol(cultivars_genesTable$data),drop=FALSE]
        # print(simul_params@CultivarsGenes)
        simul_params_cultivarsgenes(simul_params@CultivarsGenes)
        can_run_simul$cultivarsgenes <<- TRUE
      }
      # }
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )

  #### genes table ####
  genesHiddenCols <- c(2,3,5,7,8,9,10)
  simul_params_genes(simul_params@Genes)
  genesTable <- editableDTServer(
    id = "genes",
    DTdata = shiny::reactive(simul_params_genes()),
    disableCol = shiny::reactive({c()}),
    canRm = advanced_mode,
    rownames = FALSE,
    tooltips = GENES_TOOLTIP,
    row.default = default_gene,
    row.inc = c(1),
    col.hidden = genesHiddenCols
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
        if ( # nrow(genesTable$data) == 0  ||
          checkGenesTable(genesTable$data[, -which(colnames(genesTable$data) %in% c("delete"))]) == FALSE) {
          can_run_simul$genes <<- FALSE
        }
        else {
          # rename genes in cultivars genes table
          # remove line -> remove genes in cultivars genes
          if (genesTable$col == 0 && nrow(simul_params@Genes) > nrow(genesTable$data)) {
            # print("remove here")
            simul_params@CultivarsGenes <<- simul_params@CultivarsGenes[, -c(genesTable$row), drop = FALSE]
            printVerbose(paste0("set Cultivars Genes ", simul_params@CultivarsGenes))
          }
          # add line -> add a genes in cultivars genes
          if (nrow(simul_params@Genes) < nrow(genesTable$data)) {
            simul_params@CultivarsGenes <<- cbind(simul_params@CultivarsGenes, rep(0, nrow(simul_params@CultivarsGenes)))
          }
          colnames(simul_params@CultivarsGenes) <<- genesTable$data[, 1]
          simul_params_cultivarsgenes(simul_params@CultivarsGenes)

          simul_params <<- setGenes(simul_params, genesTable$data[, -ncol(genesTable$data), drop=FALSE])
          simul_params_genes(simul_params@Genes)
          can_run_simul$genes <<- TRUE
          
          if( nrow(simul_params@Genes) == 0) shinyjs::disable("addcultivar")
          else shinyjs::enable("addcultivar")
          
        }
      }
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )
  shiny::observeEvent(input$addgene,{
  
    index <- 1
    repeat{
      newGeneId <- sum(length(grep(input$listgenes, genesTable$data[,1])))+index
      geneName <- paste0(input$listgenes,"_",newGeneId)
      if( length((grep(geneName, genesTable$data[,1]))) == 0) break;
      index <- index + 1
    }
    newGene <- loadGene(geneName, type=input$listgenes)
    simul_params@CultivarsGenes <<- cbind(simul_params@CultivarsGenes, rep(0, nrow(simul_params@CultivarsGenes)))

    colnames(simul_params@CultivarsGenes) <<- newGene[, 1]
    simul_params_cultivarsgenes(simul_params@CultivarsGenes)
    
    simul_params <<- setGenes(simul_params, rbind(simul_params@Genes,newGene, stringsAsFactors = FALSE))
    simul_params_genes(simul_params@Genes)
    can_run_simul$genes <<- TRUE
    #simul_params_genes(simul_params@Genes)
  })

  # More parameters tab
  shiny::updateNumericInput(session, "nYear", value = 10)
  #shiny::updateNumericInput(session, "nTSpY", value = 120)
  shinyjs::hide("nTSpY")
  #simul_params <<- setTime(simul_params, Nyears = 10, nTSpY = 120)
  shiny::updateNumericInput(session, "seed", value = 1)
  simul_params <<- setSeed(simul_params, 1)

  ## Patho tabs default
  # shinyjs::disable(id = "patho_infection_rate")
  # shinyjs::disable(id = "patho_propagule_prod_rate")
  # shinyjs::disable(id = "patho_latent_period_mean")
  # shinyjs::disable(id = "patho_latent_period_var")
  # shinyjs::disable(id = "patho_infectious_period_mean")
  # shinyjs::disable(id = "patho_infectious_period_var")
  # shinyjs::disable(id = "patho_survival_prob")
  # shinyjs::disable(id = "patho_repro_sex_prob")
  # shinyjs::disable(id = "patho_sigmoid_kappa")
  # shinyjs::disable(id = "patho_sigmoid_sigma")
  # shinyjs::disable(id = "patho_sigmoid_plateau")
  # shinyjs::disable(id = "patho_sex_propagule_viability_limit")
  # shinyjs::disable(id = "patho_sex_propagule_release_mean")

  shinyjs::disable(id = "treatment_days_interval")
  #shinyjs::disable(id = "treatment_day_start")
  shinyjs::disable(id = "treatment_cultivars_select")
  shinyjs::disable(id = "treatment_efficiency")
  #shinyjs::disable(id = "treatment_degradation_rate")
  #shinyjs::disable(id = "treatment_cost")
  
  shinyjs::hide(id = "rotationcheck")
  
  shinyjs::hide(id = "cultivars-tableEDT")
  outputOptions(output, 'cultivars-tableEDT', suspendWhenHidden=FALSE)

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

#  observeEvent(input[["inputtabpanel"]], {
#    if (input[["inputtabpanel"]] == "Cultivars and Genes") {
#      removeCssClass("inputside", "col-sm-6")
#      removeCssClass("inputside", "col-sm-0")
#      removeCssClass("outputside", "col-sm-12")
#      removeCssClass("outputside", "col-sm-6")
#      addCssClass("outputside", "col-sm-0")
#      addCssClass("inputside", "col-sm-12")
#      shinyjs::showElement(id = "inputside")
#      shinyjs::hideElement(id = "outputside")
#    }
#  })

  
  ## Load the strategy after loading page
  isolate(update_pathogen())
  isolate(update_demo())
  shinyjs::click("About")
  
}
