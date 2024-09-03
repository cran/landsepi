# Module DT editable #
# Permet d'afficher un tableau
# et d'editer les valeurs
#
# ui.R : editableDTUI("montableau")
# server.R :
# mon_tableau_modifie <- callModule(editableDT, id = "montableau", DTdata = mon_tableau, disableCol = c("colonnes","non","modifiable"))
# observeEvent(mon_tableau_modifie$data, { message("mon tableau a ete modifie") })
#
####
## !!! ATTENTION !!!
## DT input$ retourne les indices de colonnes en commençant par 0
## Cela est du au faite qu'on n'affiche pas le nom de la ligne rownames = FALSE
## donc l'indice 2 dans R sera à 1 pour DT
####

# Ajout un bouton de suppression a chaque ligne du DT
# Param df le data.frame
# Param id un id pour le boutton
# Param mns le ns du module avec l'event id "module_clickbutton"
deleteButton <- function(df, id, mns, ...) {
  f <- function(i) {
    as.character(
      shiny::actionButton(
        paste(id, i, sep = "_"),
        label = NULL,
        icon = shiny::icon("trash"),
        onclick = paste0('Shiny.setInputValue(\"', mns, '\", this.id, {priority: \"event\"})')
      )
    )
  }

  deleteCol <- character(0)
  if (nrow(df) > 0) {
    deleteCol <- unlist(lapply(seq_len(nrow(df)), f))
  }

  return(data.frame(delete = deleteCol))
}



# UI Part
editableDTUI <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(outputId = ns("tableEDT"))
}

# Server Part
# Param DTdata a reactive data.frame
# Param disableCol reactiveVal for colnames not editable
# Param canRM reactiveVal TRUE add delete button otherwise not
# Param rownames TRUE if show rownames FALSE otherwise
# Param tooltips header tooltips message
# Param row.default a default row (when adding a new line) depending of row.cols if not all cols existing
# Param row.colsid list of columns that apply row.default
# Param row.inc list of columns to incremente value at new line
# Param col.hidden index of columns to hide.
editableDTServer <- function(id, DTdata, disableCol = shiny::reactiveVal(c()), canRm = shiny::reactiveVal(TRUE), rownames = FALSE, tooltips = NULL, row.default = NULL, row.colsid = NULL, row.inc = NULL, col.hidden = c()) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
    
      # la donnee devient reactive pour pouvoir la retourner et l'invalider
      # en gros on peut utiliser rv en sortie de editableDT
      # et surveiller rv$data coté serveur pour voir les modifs
      # rv <- shiny::reactiveValues(data = isolate({DTdata()}), value = "", row = NA, col = NA)
      rv <- shiny::reactiveValues(data = NULL, value = NULL, row = NA, col = NA)
      # render DT pour UI
      output$tableEDT <- DT::renderDT(
        {
          #rv$data <- DTdata()
          if (canRm()) {
            rv$data <- cbind(DTdata(), deleteButton(DTdata(), "button", ns("deletePressed")))
          } else {
            rv$data <- DTdata()
          }
        },
        rownames = rownames,
        ## if show rownames cols indice start at 0
        ## otherwise it start at 1
        editable = { if( rownames ){
          if(canRm()) temp = list(target = "cell", disable = list(columns = c((match(disableCol(), names(DTdata()))), ncol(DTdata())+1)))
          else temp = list(target = "cell", disable = list(columns = c(match(disableCol(), names(DTdata())))))
        }
        else {
          if(canRm()) temp = list(target = "cell", disable = list(columns = c((match(disableCol(), names(DTdata()))-1), ncol(DTdata()))))
          else temp = list(target = "cell", disable = list(columns = c(match(disableCol(), names(DTdata()))-1)))
        }
          #print(temp)
          temp
        },
        #editable = list(target = "cell"),
        selection = "none",
        escape = FALSE,
        server = TRUE,
        #extensions = list('FixedColumns'=NULL, 'Buttons'=NULL),
        extensions = list( 'FixedColumns'=NULL, 'Buttons'=NULL),
        options = list(
          #dom = if (canRm()) { 'tB'} else "t",
          dom = "t",
          scrollX = TRUE,
          fixedColumns = TRUE,
          columnDefs = list(list(visible=FALSE, targets=col.hidden)),
          buttons = list(list( extend = "collection", text="Add line", icon = shiny::icon("plus"),
                                action = DT::JS(paste0("function ( e, dt, node, config ) {
                                      Shiny.setInputValue('",id,"-addLine","', true,{priority: 'event'});
                                   }"))))
        ),
        callback = if( !is.null(tooltips)) {
          JS(paste0("var tips = ['",paste0(tooltips,sep='',collapse='\',\''),"'],
                    header = table.columns().header();
                    for (var i = 0; i < tips.length; i++) {
                    $(header[i]).attr('title', tips[i]);
                    }",
                    "table.on('keydown', function(e){",
                    "  var keys = [9,13,37,38,39,40];",
                    "  if(e.target.localName == 'input' && keys.indexOf(e.keyCode) > -1){",
                    "    $(e.target).trigger('blur');",
                    "  }",
                    "});"
                    
                    ))
                  }
        else {JS("table.on('keydown', function(e){",
                 "  var keys = [9,13,37,38,39,40];",
                 "  if(e.target.localName == 'input' && keys.indexOf(e.keyCode) > -1){",
                 "    $(e.target).trigger('blur');",
                 "  }",
                 "});")}
      )

    
      # Le proxy pour les MAJ
      proxy <- DT::dataTableProxy(outputId = "tableEDT", session = session)
    
      # Une cellule est édité
      shiny::observeEvent(input$tableEDT_cell_edit, {
        message(input$tableEDT_cell_edit)
    
        thecell <- input$tableEDT_cell_edit
        thecell$value <- gsub("[^a-zA-Z0-9\\.\\-]", "_", thecell$value,)

        print(thecell$value)
        
        isolate({
          i <- thecell$row
          if(rownames){
            j <- thecell$col
            if(j == 0) rownames(rv$data)[i] <- thecell$value
            else rv$data[i, j] <-  DT::coerceValue(thecell$value, rv$data[i, j])
            rv$value <- thecell$value
            rv$row <- i
            rv$col <- j
          }
          else {
            j <- thecell$col + 1 ### ATTENTION DT retourne les indices de columns en JS ca commence à 0
            rv$data[i, j] <- DT::coerceValue(thecell$value, rv$data[i, j])
            rv$value <- rv$data[i, j]
            rv$row <- i
            rv$col <- j
          }
        })
        
        # on met à jour rv$data pour être à jour coté serveur
        # shiny::isolate({
        #   rv$data[i, j] <- DT::coerceValue(thecell$value, rv$data[i, j])
        #   rv$value <- rv$data[i, j]
        #   rv$row <- i
        #   rv$col <- j
        # })
        ## force reactive value return
        #shiny::observe({ rv$value <- rv$data[i, j]})
        # on met a jour le tableau cote client (normalement il n'y a pas de changement mais ca sera à jour)
        #DT::replaceData(proxy = proxy, data = rv$data, resetPaging = TRUE, rownames = rownames)
      })
    
      # On assume que du moment qu'on peut ajouter une ligne on peut en supprimer une via la 
      # colonne delete
      shiny::observeEvent(input$addLine, {
        print("addLine")

        nbcol <-  ncol( rv$data)
        namecol <- colnames(rv$data)
        if(! 'delete' %in% names(rv$data) )
        {
          nbcol <-  nbcol+1
          namecol <- c(colnames(rv$data),"delete")
        }
        
        newline <- matrix(rep(c("0"), nbcol), byrow = TRUE, ncol = nbcol)
        # if default value exists
        if( !is.null(row.default) ){
          if( is.null(row.colsid)) newline <- row.default
          else {
            j <- 1
            lapply(row.colsid, FUN = function(x){newline[x] <<- row.default[j];j<<-j+1;})
          }
        }
        #print(newline)
        
        newline[nbcol] <- as.character(
          shiny::actionButton(
            paste("button", nrow(rv$data)+1, sep = "_"),
            label = NULL,
            icon = shiny::icon("trash"),
            onclick = paste0('Shiny.setInputValue(\"', ns("deletePressed"), '\", this.id, {priority: \"event\"})')
          )
        )
        newline <- as.data.frame(newline)
        colnames(newline) <- namecol
        
        # incremente cols
        if( !is.null(row.inc) ) {
          lapply(row.inc, function(x){
            if( is.na(as.numeric(newline[x])) ) {
              j <- nrow(rv$data)+1
              while( paste0(newline[x],j) %in% rv$data[,x] ) j <- j+1
              newline[x] <<- paste0(newline[x],j)
            } else {
              val <- max(sort(as.numeric(rv$data[,x]))) +1
              newline[x] <<- val
            }
          })
        }
        
        
        sapply(1:nbcol, function(i){class(newline[,i]) <<- class(rv$data[,i]); mode(newline[,i]) <<- mode(rv$data[,i]); })
        
        #print(sapply(newline,mode))
        #DT::addRow(proxy,newline) # addRow bug, du coup on met a jour tout le tableau "server=FALSE"
        shiny::isolate(rv$data <- rbind(rv$data, newline))
        #print(sapply(rv$data,mode))
        # if (canRm() == TRUE) {
        #   proxy %>%
        #     DT::replaceData(data = cbind(rv$data, deleteButton(rv$data, "button", ns("deletePressed"))), resetPaging = FALSE, rownames = FALSE)
        # } else {
          proxy %>%
            DT::replaceData(data = rv$data, resetPaging = TRUE, rownames=rownames)
        # }
      })
    
      shiny::observeEvent(input$deletePressed, {
        id <- as.integer(sub(".*_([0-9]+)", "\\1", input$deletePressed))
        # shinyalert::shinyalert(
        #   paste0("-> Remove line ",id," !"),
        #   #paste0("Values : ", paste(rv$data[id,-ncol(rv$data)],collapse = " | "), "\n\nAre you sure ?"),
        #   paste0("Values : <b>", paste(rownames(rv$data[id,]),collapse = " | "), "</b>\n\nAre you sure ?"),
        #   type = "warning",
        #   closeOnEsc = FALSE,
        #   showCancelButton = TRUE,
        #   html = TRUE,
        #   size = "m",
        #   callbackR = function(x) {
        #     if (x == TRUE) {
              rv$row <- id
              rv$value <- rv$data[rv$row, ]
              rv$data <- rv$data[-rv$row, ,drop = FALSE]
              rv$col <- 0
              
              shiny::isolate(rv$data <- cbind(rv$data[,-ncol(rv$data), drop=FALSE], deleteButton(rv$data, "button", ns("deletePressed"))))
              proxy %>%
                DT::replaceData(data = rv$data, resetPaging = TRUE, rownames=rownames)
        #     }
        #   }
        # )
      })
      
      return(rv)
      #return(reactive({c(rv$data, rv$row, rv$col, rv$value)}))
    }
  )
}
