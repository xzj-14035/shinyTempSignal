#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import ggtree
#' @import ggplot2
#' @import stringr
#' @import ggprism
#' @import Cairo
#' @import ggpmisc
#
#' @importFrom shinyjs toggle
#' @rawNamespace import(ggpubr, except = rotate)
#' @rawNamespace import(ape, except = rotate)
#' @importFrom ggtree rotate
#' @importFrom treeio read.beast
#' @importFrom treeio read.codeml
#' @importFrom treeio merge_tree
#' @importFrom treeio rescale_tree
#' @importFrom ape extract.clade
#' @importFrom treeio read.phylip.tree
#' @importFrom forecast forecast
#' @importFrom stats cor
#' @importFrom stats lm
#' @importFrom stats qqnorm
#' @importFrom stats rstudent
#' @importFrom stats pt
#' @importFrom stats shapiro.test
#' @importFrom stats ts
#' @importFrom stats acf
#' @importFrom DescTools RunsTest
#' @importFrom stats na.omit
#' @importFrom stats as.formula
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @noRd

app_server <- function(input, output, session)  {
  options(shiny.maxRequestSize = 4000*1024^2)
  category <- ..eq.label.. <- ..rr.label.. <- NULL
  click_node <- function(x, y, tr) {
    sq_dx <- (x - tr$x)^2
    sq_dy <- (y - tr$y)^2
    i <- which.min(sq_dx + sq_dy)
    node <- tr$node[i]
    return(node)
  }

  observeEvent(input$plotClick, {
    p <- reactive(ggtree(tree()))
    x <- reactive(as.numeric(input$plotClick$x))
    y <- reactive(as.numeric(input$plotClick$y))
    node <- reactive(click_node(x(), y(), p()$data) |> as.character())
    updateTextInput(session,"node",value = node())
  })

  observeEvent(input$node, {
    req(input$node != "")
    tree2plot <- reactive(extract.clade(tree(), node = as.numeric(input$node)))
    output$plot1 <- renderPlot({
    if (input$tip) {
      p <- reactive(ggtree(tree2plot(),color=input$color3, size=input$size) + geom_tiplab()+geom_nodelab(aes(label=node),hjust=-.3))
      if (!is.null(all_data())) {
        p()%<+%all_data() + geom_tippoint(aes(shape=category,color=category))+
          geom_tiplab(aes(color=category))+
          scale_color_manual(values = c("up" = "red", "down" = "blue"))
      } else {
        p()
      }
    } else {
      ggtree(tree2plot(),color=input$color3, size=input$size)+geom_nodelab(aes(label=node),hjust=-.3)
    }
    }, height = reactive(input$height))})

  #1.读入树文件
  tree <- reactive({
    req(input$treefile)
    if (input$filetype=="Newick") {
      tree <- read.tree(input$treefile$datapath) %>% as.phylo()
    } else if (input$filetype=="beast") {
      tree <- read.beast(input$treefile$datapath) %>% as.phylo()
    } else if (input$filetype=="Nexus") {
      tree <- read.nexus(input$treefile$datapath)
    } else if (input$filetype=="phylip") {
      tree <- read.phylip.tree(input$treefile$datapath) %>% as.phylo()
    }
  })
#      if (as.numeric(input$node)<length(as.phylo(tree)$tip.label)) {
#      stop("it is a tip label")
#    }
  # initilize tree data
  date <- reactive({
    req(tree())
    dateType3(tree = tree(), pattern = input$regression) |>
      dateNumeric(format = input$format)
  })
  divergence <- reactive({
    req(tree())
    getdivergence(tree = tree())})
  label <- reactive({as.phylo(tree())$tip.label})
  df <- reactive(data.frame(label=label(),
    date=date(),divergence=divergence()
  ))
  # prepare up.table
  pd <- reactive({
    req(df())
    estimate(df(),p=input$pvalue)})
  up.table <- reactive(
    {req(df()); temp <- df()[pd()$up,,drop=F]; temp$category <- rep("up", times = length(pd()$up)); temp})
#  up.table()$category <- reactive({req(up.table()); rep("up",nrow(up.table()))})
  # prepare down.table
  down.table <- reactive({
    req(pd() & pd())
    temp <- df()[pd()$down,,drop=F]
    temp$category <- rep("down", times = length(pd()$down))
    temp
  })
#  down.table()$category <- reactive({req(down.table()); rep("down",nrow(down.table()))})
  d <- reactive({req(up.table() & down.table()); rbind(up.table(), down.table())})
  # prepare keep.table
  keep <- reactive({req(pd()); pd()$up==pd()$down})
  keep.table <- reactive({req(df() & keep()); df()[keep(),,drop=F]})
  # prepare exclude.table
  keep <- reactive({req(pd()); pd()$up==pd()$down})
  exclude.table <- reactive({req(df() & keep()); df()[!keep(),,drop=F]})
  # prepare all
  all_data <- reactive({req(d() & df()); merge(d(), df(), by = "label", all = T)})

  regression <- reactive({
    req(input$regression_btn)
    lm(as.formula(paste(input$y_var, "~", input$x_var)), all_data())
  })
  a <- reactive({req(tree()); length(tree()$tip.label) + 1})
  b <- reactive({req(tree()); length(tree()$tip.label) + tree()$Nnode})
  table2 <- reactive({
    req(a() & b())
    generate_dataframe(a = a(), b = b(),
      tree = tree(), regression = regression(), format = input$format
  )})

  data2 <- reactive({
    req(input$outdata)
    read.csv(input$outdata$datapath)
  })
  
  merged_df <- reactive({
    req(data2())
    # 将上传的表格与现有的表格合并
    merge(df(), data2()) %>% unique() %>% na.omit()
  })

  estimate2 <- function(lm,p){
    rst <- rstudent(lm)
    down <- 0.5-abs(0.5-pt(rst,lm$df.residual))< p / 2&rst<0
    up <- 0.5-abs(0.5-pt(rst,lm$df.residual))< p/ 2&rst>0
    return(list(down=down,up=up))}

  need.df <- reactive({
    req(input$x_var & input$y_var)
    merged_df()[,c("label",input$x_var,input$y_var)]}
  )
  regression2 <- reactive({
    req(input$regression_btn)
    lm(as.formula(paste(input$y_var, "~", input$x_var)), need.df())
  })
  need.pd <- reactive({
    req(regression2() & input$pvalue)
    estimate2(lm=regression2(),p=input$pvalue)}
  )
  need.up.table <- reactive({
    req(need.df() & pd())
    temp <- need.df()[need.pd()$up,,drop=F]
    temp$category <- rep("up", times = length(need.pd()$up))
    temp
  })
  need.down.table <- reactive({
    req(need.pd() & need.df())
    temp <- need.df()[need.pd()$down, , drop = F]
    temp$category <- rep("down", times = length(need.pd()$down))
    temp
  })
  need.keep <- reactive({
    req(need.pd())
    need.pd()$up == need.pd()$down
  })
  need.keep.table <- reactive({
    req(need.df() & need.keep())
    need.df()[need.keep(), , drop = F]
  })
  need.exclude.table <- reactive({
    req(need.df() & need.keep())
    need.df()[!need.keep(), , drop = F]}
  )
  table3 <- reactive({
    req(a() & b() & regression2() & tree())
    generate_dataframe(tree = tree(), regression = regression2(),
      a = a(), b = b(), format = input$format)
})





  #全部在外面取出来不就好了
  estimate <- function(df,p) {
    lm=lm(df$divergence~ df$date)
    rst <- rstudent(lm)
    down <- 0.5-abs(0.5-pt(rst,lm$df.residual))< p / 2&rst<0
    up <- 0.5-abs(0.5-pt(rst,lm$df.residual))< p/ 2&rst>0
    return(list(down=down,up=up))
  }

  output$plot1 <- renderPlot({
    all <- reactive(merge(d(), df(),by='label',all = T))
    if (input$tip) {
      p <- reactive(ggtree(tree(),color=input$color3, size=input$size) + geom_tiplab()+geom_nodelab(aes(label=node),hjust=-.3))
      p()%<+%all_data()+ geom_tippoint(aes(shape=category,color=category))+
        geom_tiplab(aes(color=category))+
        scale_color_manual(values = c("up" = "red", "down" = "blue"))
    }else{
      ggtree(tree(),color=input$color3, size=input$size)+geom_nodelab(aes(label=node),hjust=-.3)
    }
    }, height = reactive(input$height))
  #2.取出日期
  output$datetable <- renderDataTable({
    data.frame(label=label(),
        date=date(),divergence=divergence()
    )})



  output$download1 <- downloadHandler(
    filename = function(){
      "Sample-Dates.csv"
    },
    content = function(file){
      write.csv(df(),file,row.names = FALSE)
    }
  )
  
  #3.取出divergence，回归分析
  output$plot2 <- renderPlot({
    ggplot(df(), aes(x = date, y = divergence)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, formula = y ~ x,colour=input$color2) +
      geom_point(data = down.table(), aes(x = date, y = divergence), color = 'blue') +
      geom_point(data = up.table(),aes(x = date, y = divergence), color ="red")+
      # geom_text(data = d, aes(x = date, y = divergence, label = label)) +
      theme_bw() +
      stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
  })
  output$outliers <- renderDataTable({
    rbind(up.table(), down.table())
  })
  observeEvent(input$delete,{
    output$plot2 <- renderPlot({
      ggplot(keep.table(), aes(x = date, y = divergence)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, formula = y ~ x,colour=input$color2) +
        geom_point(data = exclude.table(), aes(x = date, y = divergence), color = 'gray') +
        # geom_text(data = d, aes(x = date, y = divergence, label = label)) +
        theme_bw() +
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
    })
    output$plot3 <- renderPlot({
      ggplot(keep.table(), aes_string(x = input$x_var, y = input$y_var)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, formula = y ~ x,colour=input$color2) +
        geom_point(data = exclude.table(), aes_string(x = input$x_var, y = input$y_var), color = 'gray') +
        # geom_text(data = d, aes(x = date, y = divergence, label = label)) +
        theme_bw() +
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
    })
    print("here do delete")
    output$plot1 <- renderPlot({
      bowser()
      to_drop <- reactive(c(down.table()$label,up.table()$label))
      tip_reduced <- reactive(drop.tip(tree(), to_drop()))
      if (input$tip) {
        ggtree(tip_reduced()) + geom_tiplab()+geom_nodelab(aes(label=node()),hjust = -.3)
      }else{
        ggtree(tip_reduced()) + geom_nodelab(aes(label=node()),hjust = -.3)
      }
    }, height = reactive(input$height))
  })
  observeEvent(input$reset,{
    output$plot1 <- renderPlot({
      if (input$tip) {
        p <- ggtree(tree(),color=input$color3, size=input$size) + geom_tiplab()+geom_nodelab(aes(label=node),hjust=-.3)
        p%<+%all_data()+ geom_tippoint(aes(shape=category,color=category))+
          geom_tiplab(aes(color=category))+
          scale_color_manual(values = c("up" = "red", "down" = "blue"))
      }else{
        ggtree(tree(),color=input$color3, size=input$size)+geom_nodelab(aes(label=node),hjust=-.3)
      }
    }, height = reactive(input$height))
    
    output$plot2 <- renderPlot({
      ggplot(df(), aes(x = date, y = divergence)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, formula = y ~ x,colour=input$color2) +
        geom_point(data = down.table(), aes(x = date, y = divergence), color = 'blue') +
        geom_point(data = up.table(),aes(x = date, y = divergence), color ="red")+
        # geom_text(data = d, aes(x = date, y = divergence, label = label)) +
        theme_bw() +
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
    })
    output$plot3 <- renderPlot({
      if(is.null(merged_df())) {
        ggplot(df(), aes_string(x = input$x_var, y =input$y_var)) +
          geom_point() +
          geom_smooth(method = "lm", se = FALSE, formula = y ~ x,colour=input$color2) +
          geom_point(data = need.down.table(), aes_string(x = input$x_var, y =input$y_var), color = 'blue') +
          geom_point(data = need.up.table(),aes_string(x = input$x_var, y =input$y_var), color ="red")+
          # geom_text(data = d, aes(x = date, y = divergence, label = label)) +
          theme_bw() +
          stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
      } else {
        ggplot(merged_df(), aes_string(x = input$x_var, y =input$y_var)) +
          geom_point() +
          geom_smooth(method = "lm", se = FALSE, formula = y ~ x,colour=input$color2) +
          geom_point(data = need.down.table(), aes_string(x = input$x_var, y =input$y_var), color = 'blue') +
          geom_point(data = need.up.table(),aes_string(x = input$x_var, y =input$y_var), color ="red")+
        # geom_text(data = d, aes(x = date, y = divergence, label = label)) +
        theme_bw() +
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
      }
    })
    print("here do reset")
  })
  
  generate_dataframe <- function(a, b, tree, regression, format) {
    sub.tree <- vector(mode = "list", length = length(a - b + 1))
    date <- vector(mode = "list", length = length(a - b + 1))
    divergence <- vector(mode = "list", length = length(a - b + 1))
    df_td <- vector(mode = "list", length = length(a - b + 1))
    model.results <- vector(mode = "list", length = length(a - b + 1))
    for (i in a:b) {
      f=i-a+1
      sub.tree[[f]] <- extract.clade(tree,node = i)
      date[[f]] <- dateType3(sub.tree[[f]],pattern = regression) %>% dateNumeric(format = format) 
      divergence[[f]] <- getdivergence(tree = sub.tree[[f]])
      df_td[[f]] <- as.data.frame(cbind(date=date[[f]],divergence=divergence[[f]]))
      if (unique(is.na(divergence[[f]]))) {
        model.results[[f]] <- NA
      } else{
        m <- lm(divergence ~ date, data = df_td[[f]])
        rst <- rstudent(m)
        upval <- c((0.5 - abs(pt(
          rst, m$df.residual
        ) - 0.5))
        < input$pvalue/ 2 & rst > 0)
        downval <- c((0.5 - abs(pt(
          rst, m$df.residual
        ) - 0.5))
        < input$pvalue / 2 & rst < 0)
        modele <-
          summary(lm(divergence ~ date, data = df_td[[f]]))
        model.results[[f]] <- data.frame(
          node = i,
          tip.number=Ntip(sub.tree[[f]]),
          r.squared = modele$r.squared,
          adj.r.squared = modele$adj.r.squared,
          pvalue = modele$coefficients[nrow(modele$coefficients), ncol(modele$coefficients)],
          slope = modele$coefficients[nrow(modele$coefficients), 1],
          intercept = modele$coefficients[1, 1],
          up = length(which(upval == T)),
          down = length(which(downval == T)),
          total_abnormal = length(which(upval == T)) + length(which(downval == T))
          
        )
      }}
    df <- na.omit(do.call(rbind, model.results))
    df[order(df$total_abnormal, decreasing = T), ]
  }

  output$dataframe <- renderDataTable({
    table2()
  })

  output$download2.table <- downloadHandler(
    filename = "Subtree_regression_intergration.csv",
    content = function(file){
      write.csv(table2(), file, row.names = FALSE)
    }
  )
  

  observeEvent(merged_df(), {
    updateSelectInput(session, "x_var", choices = colnames(merged_df()))
    updateSelectInput(session, "y_var", choices = colnames(merged_df()))
  })
  
  output$data_table <- renderDataTable({
    req(merged_df())
    merged_df()
  })

  
  
  output$plot3 <- renderPlot({
    req(merged_df())
    ggplot(merged_df(), aes_string(x = input$x_var, y =input$y_var)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, formula = y ~ x,colour=input$color2) +
      geom_point(data = need.down.table(), aes_string(x = input$x_var, y =input$y_var), color = 'blue') +
      geom_point(data = need.up.table(),aes_string(x = input$x_var, y =input$y_var), color ="red")+
      # geom_text(data = d, aes(x = date, y = divergence, label = label)) +
      theme_bw() +
      stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
  })
  
  output$outliers2 <- renderDataTable({
    rbind(need.up.table(),need.down.table())
  })
  
  output$out_dataframe <- renderDataTable({
    table3()
  })

  output$download3.table <- downloadHandler(
    filename = function(){
      "out_data__regression_intergration.csv"
    },
    content = function(file){
      write.csv(table3(),file,row.names = FALSE)
    }
  )
  
}

