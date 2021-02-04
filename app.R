# .............................................................................
# Projeto de alocação de voluntários em turmas de inglês,
# conforme preferências e restricões apresentadas na proposta.
# 
# Desenvolvido por Vinicius Aguiar.
# .............................................................................

# Bibliotecas
library(memoise)
library(gsubfn)
library(GA)
library(shinythemes)
library(shiny)
library(DT)

# .............................................................................
# Algoritmo Genético ..........................................................
# .............................................................................

geneticAlg <- function(vol, vol.eq, class, class.eq, 
                       class.prog, class.level, pop, iter){
  
  # Pré-processamento .........................................................
  
  # Renomeia as entradas
  colnames(class.eq) <- c("id", "id2")
  colnames(class.prog) <- c("current", "future")
  colnames(class.level) <- c("current", "future", "age")
  colnames(class) <- c("id", "age", "level", "period")
  colnames(vol) <- c("id", "name", "gender", "job", "date", "owner", 
                     "current", "commitment", "level", "period", "class","age")
  
  # Seleciona EAs, owners e alocáveis
  ea <- grep("EA", vol$job)
  own <- which(vol$owner == 1)
  alloc <- setdiff(1:nrow(vol), union(ea, own))
  
  # Remove elementos repetidos de vol.eq
  vol.eq <- t(apply(vol.eq, 1, sort))
  vol.eq <- data.frame(unique(vol.eq))
  colnames(vol.eq) <- c("id", "id2")
  
  # Remove EAs de vol.eq
  vol.eq <- vol.eq[-which(vol.eq$id %in% ea), ]
  vol.eq <- vol.eq[-which(vol.eq$id2 %in% ea), ]
  
  # Transforma nível de inglês em letras em vol e class
  new <- data.frame("Teens 4.2", "D")
  colnames(new) <- c("id", "id2")
  class.eq <- rbind(class.eq, new)
  
  for(i in 1:nrow(class.eq)){
    j <- grep(class.eq$id[i], vol$level)
    if(length(j) > 0){
      vol$level[j] <- class.eq$id2[i]
    }
    
    j <- grep(class.eq$id[i], class$level)
    if(length(j) > 0){
      class$level[j] <- class.eq$id2[i]
    }
  }
  
  # Armazena constantes
  mean.class <- length(alloc) / nrow(class)
  mean.comm <- mean(vol$commitment[union(own, alloc)])
  
  new.vol = max(vol$date)
  vet <- intersect(union(own, alloc), which(vol$date != new.vol))
  fresh <- intersect(union(own, alloc), which(vol$date == new.vol))
  ratio.vet <- length(vet)/length(fresh)
  
  wom <- intersect(union(own, alloc), grep("Fem", vol$gender))
  men <- intersect(union(own, alloc), grep("Masc", vol$gender))
  ratio.men <- length(men)/length(wom)
  
  # Cria tabela de preferências
  pref <- matrix(0, length(alloc), nrow(class))
  
  # Contabiliza as preferencias ..
  for(i in 1:length(alloc)){
    j <- alloc[i]
    
    # .. por periodo
    if(vol$period[j] != "Indiferente"){
      k <- grep(vol$period[j], class$period)
      if(length(k) > 0){
        pref[i, -k] <- pref[i, -k] - 3
        pref[i, k] <- pref[i, k] + 3
      }
    }
    
    # .. por turma
    if(substr(vol$class[j], 1, 10) == "Acompanhar"){
      k <- grep(vol$current[j], class.prog$current)
      if(length(k) > 0){
        k <- class.prog$future[k]
        k <- grep(k, class$id)
        pref[i, k] <- 1
      }
    } else if(substr(vol$class[j], 1, 10) == "Permanecer"){
      k <- grep(vol$current[j], class$id)
      if(length(k) > 0){
        k <- class$level[k]
        k <- grep(k, class$level)
        pref[i, k] <- 1
      }
    }
    
    # .. por faixa etaria
    if(vol$age[j] != "Indiferente"){
      k <- grep(vol$age[j], class$age)
      if(length(k) > 0){
        pref[i, -k] <- pref[i, -k] - 1
        pref[i, k] <- pref[i, k] + 1
      }
    }
    
    # .. por nivel de ingles
    k <- which((class$level) <= (vol$level[j]))
    pref[i, -k] <- pref[i, -k] - 1
    pref[i, k] <- pref[i, k] + 1
  }
  
  # Algoritmo Genético ........................................................
  
  # Avalia a adaptabilidade do cromossomo
  fitness <- function(x){
    len <- length(x)
    size.class <- table(x)
    
    # Maximiza a preferencias
    f <- 0
    for(i in 1:len){
      f <- f + pref[i, x[i]]
    }
    
    # Considera positiva quando as 
    # restricoes de mesma escala sao satisfeitas
    g <- 0
    for(k in 1:nrow(vol.eq)){
      i <- vol.eq$id[k]
      j <- vol.eq$id2[k]
      
      if(i %in% alloc){
        i <- which(alloc == i)
        i <- (class$period[x[i]])
      } else {
        i <- grep(vol$current[i], class$id)
        i <- (class$period[i])
      }
      
      if(j %in% alloc){
        j <- which(alloc == j)
        j <- (class$period[x[j]])
      } else {
        j <- grep(vol$current[j], class$id)
        j <- (class$period[j])
      }
      
      if(i == j){
        g <- g + 1
      }
    }
    
    # Penaliza turmas com tamanho distante da media
    c1 <- sum(abs(size.class - mean.class))
    c1 <- c1 +(nrow(class) - length(size.class)) * (mean.class)
    
    # Penaliza turmas com ..
    c2 <- 0
    c3 <- 0
    c4 <- 0
    for(i in 1:nrow(class)){
      j <- own[grep(class$id[i], vol$current[own])]
      j <- union(j, alloc[which(x == i)])
      
      # .. comprometimento distante da media
      c2 <- c2 + abs(mean.comm - mean(vol$commitment[j]))
      
      # .. proporcao vet:cal distante da media
      ratio <- length(intersect(vet, j))/length(intersect(fresh, j))
      c3 <- c3 + abs(ratio.vet - ratio)
      
      # .. proporcao de hom:mul distante da media
      ratio <- length(intersect(men, j))/length(intersect(wom, j))
      c4 <- c4 + abs(ratio.men - ratio)
    }
    
    return(f + g - c1 - c2 - c3 - c4)
  }
  
  mfitness <- memoise(fitness)
  
  # Gera a populacao inicial
  population <- function(obj, ...){
    pop <- matrix(sample(1:nrow(class), obj@popSize*length(obj@lower), replace = T),
                  nrow = obj@popSize)
    
    return(pop)
  }
  
  # Controla o processo de recombinacao
  crossover <- function(obj, prnts, ...){
    len <- length(obj@lower)
    i <- sample(1:len, 1)
    
    father <- obj@population[prnts[1], i:len]
    mother <- obj@population[prnts[2], i:len]
    obj@population[prnts[1], i:len] <- mother
    obj@population[prnts[2], i:len] <- father
    
    children <- rbind(obj@population[prnts[1], ],
                      obj@population[prnts[2], ])
    
    fit <- c(NA, NA)
    if(all(father == mother))
      fit <- c(obj@fitness[prnts[1]], obj@fitness[prnts[2]])
    
    return(list("children" = children, "fitness" = fit))
  }
  
  # Transforma um alelo do cromossomo
  mutation <- function(obj, prnt, ...){
    i <- sample(1:length(obj@lower), 1)
    j <- sample(1:nrow(class), 1)
    obj@population[prnt, i] <- j
    
    return(obj@population[prnt, ])
  }
  
  # Executa o algoritmo genetico
  GA <- ga(type = "real-value",
           fitness = mfitness,
           lower = rep(1, length(alloc)),
           upper = rep(nrow(class), length(alloc)),
           population = population,
           crossover = crossover,
           mutation = mutation,
           popSize = pop, maxiter = iter, monitor = T)
  
  # Pos-processamento .........................................................
  
  sol.ga <- GA@solution
  sol.vol <- data.frame("Turma" = vol$current,
                        "Voluntario" = vol$name,
                        "Funcao" = vol$job,
                        "Owner" = vol$owner,
                        "Nivel" = rep("", nrow(vol)),
                        "Pref.Periodo" = rep("", nrow(vol)),
                        "Pref.Turma" = rep("", nrow(vol)),
                        "Pref.Faixa.Etaria" = rep("", nrow(vol)),
                        stringsAsFactors = F)
  sol.class <- data.frame("Turma" = class$id,
                          "Tamanho" = rep(0, nrow(class)),
                          "Comprometimento" = rep(0, nrow(class)),
                          "Proporcao Vet:Cal" = rep(0, nrow(class)),
                          "Proporcao Hom:Mul" = rep(0, nrow(class)))
  
  # Popula os voluntarios alocados com ..
  for(i in 1:length(sol.ga)){
    j <- alloc[i]   # vol
    k <- sol.ga[i]  # class
    
    # .. nome da turma
    sol.vol[j, 1] <- (class$id[k])
    
    # .. se esta no nivel adequado
    if((class$level[k]) <= (vol$level[j])){
      sol.vol[j, 5] <- "Ok"
    } else {
      sol.vol[j, 5] <- "--"
    }
    
    # .. pref de periodo
    if(vol$period[j] != "Indiferente"){
      if((vol$period[j]) == (class$period[k])){
        sol.vol[j, 6] <- "Ok"
      } else {
        sol.vol[j, 6] <- "--"
      }
    }
    
    # .. pref de turma
    if(substr(vol$class[j], 1, 10) == "Acompanhar"){
      p <- grep(vol$current[j], class.prog$current)
      if((class$id[k]) == (class.prog$future[p])){
        sol.vol[j, 7] <- "Ok"
      } else {
        sol.vol[j, 7] <- "--"
      }
    } else if(substr(vol$class[j], 1, 10) == "Permanecer"){
      p <- grep(vol$current[j], class$id)
      if(class$level[k] == class$level[p]){
        sol.vol[j, 7] <- "Ok"
      } else {
        sol.vol[j, 7] <- "--"
      }
    }
    
    # .. pref de faixa etaria
    if(vol$age[j] != "Indiferente"){
      if((vol$age[j]) == (class$age[k])){
        sol.vol[j, 8] <- "Ok"
      } else {
        sol.vol[j, 8] <- "--"
      }
    }
  }
  
  # Para cada turma gera estatisticas de ..
  for(i in 1:nrow(class)){
    j <- grep(class$id[i], sol.vol[, 1])
    
    # .. tamanho
    sol.class[i, 2] <- length(j)
    
    # .. comprometimento
    sol.class[i, 3] <- round(mean(vol$commitment[j]), 1)
    
    # .. proporcao vet:cal
    sol.class[i, 4] <- round(
      length(intersect(vet, j))/length(intersect(fresh, j)), 1)
    
    # .. proporcao hom:mul
    sol.class[i, 5] <- round(
      length(intersect(men, j))/length(intersect(wom, j)), 1)
  }
  
  # Ordena a tabela de voluntarios por turma
  sol.vol <- sol.vol[order(sol.vol[, 1]), ]
  
  return(list(sol.class, sol.vol))
}

# .............................................................................
# Servidor ....................................................................
# .............................................................................

server <- function(input, output, session) {
  sol.vol <- NULL
  
  # Aciona o botao de otimizacao
  observeEvent(input$go, {
    
    if(is.null(input$vol) || 
       is.null(input$class) || 
       is.null(input$vol.eq) || 
       is.null(input$class.eq) || 
       is.null(input$class.prog) || 
       is.null(input$class.level)){
      vol <- read.csv("input/tb_vol_cadastro-1.csv", 
                      stringsAsFactors = F, fileEncoding = "latin1")
      vol.eq <- read.csv("input/tb_vol_mesmaescala.csv", 
                         stringsAsFactors = F, fileEncoding = "latin1")
      class <- read.csv("input/tb_tur_cadastro.csv", 
                        stringsAsFactors = F, fileEncoding = "latin1")
      class.eq <- read.csv("input/tb_dm_equivalencia.csv", 
                           stringsAsFactors = F, fileEncoding = "latin1")
      class.prog <- read.csv("input/tb_tur_progressao.csv", 
                             stringsAsFactors = F, fileEncoding = "latin1")
      class.level <- read.csv("input/tbl_dm_tipodeturmas-1.csv", 
                              stringsAsFactors = F, fileEncoding = "latin1")
    } else {
      vol <- read.csv(input$vol$datapath, 
                      stringsAsFactors = F, fileEncoding = "latin1")
      class <- read.csv(input$class$datapath, 
                        stringsAsFactors = F, fileEncoding = "latin1")
      vol.eq <- read.csv(input$vol.eq$datapath, 
                         stringsAsFactors = F, fileEncoding = "latin1")
      class.eq <- read.csv(input$class.eq$datapath, 
                           stringsAsFactors = F, fileEncoding = "latin1")
      class.prog <- read.csv(input$class.prog$datapath, 
                             stringsAsFactors = F, fileEncoding = "latin1")
      class.level <- read.csv(input$class.level$datapath, 
                              stringsAsFactors = F, fileEncoding = "latin1")
    }
    
    list[sol.class, sol.vol] <- 
      geneticAlg(vol, vol.eq, class, class.eq, 
                 class.prog, class.level, input$pop, input$iter)
    
    sol.vol <<- sol.vol
    
    output$sol.class <-
      renderDT(sol.class,
               server = T,
               editable = T,
               rownames = F,
               selection = "none",
               options = list(info = F, 
                              paging = F, 
                              searching = F,
      columnDefs = list(list(className = 'dt-center', targets = 0:4))))
    
    defs <- list(list(className = 'dt-center', targets = 0:7))
    output$sol.vol <-
      renderDT(sol.vol,
               server = T,
               editable = T,
               rownames = F,
               selection = "none",
               options = list(info = F, 
                              paging = F, 
                              searching = F,
      columnDefs = list(list(className = 'dt-center', targets = 0:7))))
  })
  
  # Proxy
  proxy <- dataTableProxy("sol.vol")
  
  # Atualiza sol.vol
  observeEvent(input$sol.vol_cell_edit, {
    info <- input$sol.vol_cell_edit
    i <- info$row
    j <- info$col + 1
    v <- info$value
    sol.vol[i, j] <<- DT::coerceValue(v, sol.vol[i, j])
    replaceData(proxy, sol.vol, resetPaging = F, rownames = F)
  })
  
  # Aciona o botao de download csv
  output$csv <- downloadHandler(
    filename = function(){
      paste0("Alocacao.", format(Sys.Date(), format="%d.%m.%Y"), ".csv")
    },
    content = function(file){
      sol.vol <- sol.vol[order(sol.vol[, 1]), ]
      write.csv(sol.vol[, 1:2], file, row.names = F)
    }
  )
  
  # Fecha o programa ao fechar sua aba no navegador
  session$onSessionEnded(function(){
    stopApp()
    #q("no")
  })
}

# .............................................................................
# Interface ...................................................................
# .............................................................................

# Define o tipo de arquivo de entrada
csv <- c(".csv", "text/csv",
         "text/comma-separated-values,text/plain")

# Seta a interface
ui <- fluidPage(
  title = "Alocação de Voluntários",
  tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
  
  # Muda o estilo da pagina
  tags$br(),
  theme = shinytheme("lumen"),
  tags$style(".form-group {margin-bottom: -10px}"),
  tags$style(".shiny-file-input-progress {display: none}"),
  
  # Sidebar
  sidebarPanel(
    
    # Entrada de CSVs
    fileInput("vol", "", accept = csv,
              buttonLabel = "Voluntários", 
              placeholder = "tb_vol_cadastro-1.csv"),
    fileInput("vol.eq", "", accept = csv,
              buttonLabel = "Mesma escala", 
              placeholder = "tb_vol_mesmaescala.csv"),
    fileInput("class", "", accept = csv,
              buttonLabel = "Turma", 
              placeholder = "tb_tur_cadastro.csv"),
    fileInput("class.level", "", accept = csv,
              buttonLabel = "Tipo de turma", 
              placeholder = "tbl_dm_tipodeturmas-1.csv"),
    fileInput("class.prog", "", accept = csv,
              buttonLabel = "Progressão de turma", 
              placeholder = "tb_tur_progressao.csv"),
    fileInput("class.eq", "", accept = csv,
              buttonLabel = "Equivalência de turma", 
              placeholder = "tb_dm_equivalencia.csv"),
    
    # Linha horizontal
    tags$hr(),
    
    tags$h4("Algoritmo Genético"),
    
    # Sliders
    sliderInput("pop", "População:",
                min = 16, max = 64, value = 16),
    tags$br(),
    sliderInput("iter", "Gerações:",
                min = 128, max = 8192, value = 128),
    
    # Aviso
    tags$br(),
    tags$em("Aperte em otimizar para executar um exemplo."),
    tags$br(),
    tags$br(),
    
    # Inicia a otimizacao
    actionButton("go", "Otimizar"),
    downloadButton("csv", "Gerar CSV")
  ),
  
  # Painel direito
  mainPanel(
    tabsetPanel(
      tabPanel("Turmas", DTOutput("sol.class")),
      tabPanel("Voluntários", DTOutput("sol.vol"))
    )
  )
)

# .............................................................................
# Execucao ....................................................................
# .............................................................................

shinyApp(ui, server)
