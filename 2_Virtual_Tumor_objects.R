library(R6)
library(tidyr)
library(dplyr)
library(purrr)
library(rlang)
library(here)
library(cowplot)

# Author: Thibaud Derippe
# File containing the Virtual Object coded as an R6 object.

source(here('0_Lindner_model_PaSM_config.r'))
# Some first initial funcions are outside the VT object, for instance related to 

data_VT <- read.table(here('data/data_cell_viability.csv'), sep = ";", header = T) %>%
  as_tibble



# Goal: analyse the curve to sample VC
curve_sampling <- function(data = data_VT, filter_data , uniqueIDperOutcome = T, celltheque, drug = 1:2){
  #
  # allconc1 <- c(0,0.08,0.16,0.32,0.64,1.3,2.60,5,10,20)
  # allconc2 <- c(0,5,10,15)


  filter_data <- enexpr(filter_data)




  # Just to know if we want to keep all id or just the minmal
  # Unless doing thing manual, users already works with minimal celltheque anyway
  # if(uniqueIDperOutcome == T){
  #   uniqueID <- celltheque_spread(returnOnIDperGroup = T, celltheque = celltheque, drug = drug)$cellid
  # }else{

  uniqueID <- unique(celltheque$cellid)
  # }



  # Ok let's start

  # Take the observation data
  data %>%
    # apply the filter (cell line, drug,....)
    filter(!!filter_data) %>%
    # reduce and nest
    group_by(ID) %>%
    nest() %>%
    # pull(data) -> x; x <- x[[1]]
    mutate(data = map(data, function(x){


      # we need to guess the main drug (one with on different conc per line)
      x %>%
        select(starts_with("conc")) %>%
        gather("key","value") %>%
        distinct() %>%
        group_by(key) %>%
        tally %>%
        filter(n == nrow(x)) %>%
        pull(key) -> maindrug



      coltemp <- parse_expr(maindrug)

      # get all the concentrations, pus (for surviving)
      allconc1 <-  x[[maindrug]] %>% unique()
      allconc1 <- c(allconc1, Inf)


      # Ok so x is one curve profile
      x %>%
        # first replace value above 100 by 100
        mutate(Value = if_else(Value > 100, 100, Value)) -> x

      # Then prevent value to reincrease by comparing a value and the previous one
      # If the value is higher, replace by the previous one
      for(a in 2:nrow(x)){

        if(x$Value[[a]] >  x$Value[[a - 1]]) x$Value[[a]] <-  x$Value[[a - 1]]

      }

      # take the data
      x %>%
        # add a new row to handle From last conc to Inf
        bind_rows({

          x %>%
            slice(1) %>%
            mutate(Value = 0) -> tempp

          tempp[maindrug] <- Inf


          tempp

        }) %>%
        # make a lag to compute the difference (providing n, pct of cell that died during the step)
        mutate(Value2 = lag(Value)) %>%
        mutate(Value2 = if_else(is.na(Value2), 100, Value2)) %>%
        select(Value, Value2, everything()) %>%
        mutate(n = Value2 - Value) %>%
        # slice(-1) %>%
        mutate(n = round(n)) -> temp



      # now constitute the bags, hardest part to code!

      temp %>%
        rowid_to_column() %>%
        group_by(rowid) %>%
        nest() %>%
        # pull(data) -> y; y <- y[[8]]
        mutate(sample = map(data, function(y){


          # take concentration A12
          A12 <- y$conc4

          # for each row
          if(maindrug == "conc1"){

            colnam <- paste0("Veneto_", A12)
          }else{

            colnam <- paste0("A11_", A12)
          }

          # print(colnam)
          # print(celltheque %>% tail)
          temp<-  celltheque[celltheque[[colnam]] ==y[[maindrug]],  ] %>%
            # filter(!!parse_expr(colnam) == y[[maindrug]]) %>%
            pull(cellid)
          # print("la")




        }))



    })) %>%
    unnest() %>%
    unnest(data) %>%
    select(n , sample, everything())



}

cell_harvest <- function( sample_df = curve_sampling(), ncol = n, samplecol = sample, returnPlot = F, xlog = T,  plot_by_con2 = F,wrap = T){


  # ncol <- expr(n); samplecol <- expr(sample)
  ncol<- enexpr(ncol)
  samplecol<- enexpr(samplecol)


  # sample_df %>%
  #   select(n, sample) %>%
  #   unnest()

  sample_df %>%
    # slice(6:7) %>%
    mutate(sample = map2(!!ncol, !!samplecol, function(nn, samplee){
      # print(nn)
      if(class(samplee) == "try-error" ) return(NA)
      if(length(samplee) == 0 ) return(NA)
      sample(c(samplee, samplee), size = max(nn,0), replace = T)

    } )) %>%
    unnest(sample) %>%
    filter(!is.na(sample)) %>%
    pull(sample) -> harvest
  # %>%
  # pull(sample)
  if(length(harvest) < 100){

    harvest <- c(harvest, sample(harvest, size = 100 - length(harvest), replace = T))

  }

  if(returnPlot == F) return(harvest)
  # Verification by reconstituting the plot



  return(list(bag = harvest, plot = reconstitution_plot(harvest, xlog = xlog,  plot_by_con2 = plot_by_con2, wrap = wrap)))
}


# Define VT ---------------------------------------------------------------
# Main VT object, containing the harvest (current 100 VCs), the celltheque (once cell per bag of interest),
# and many function for plotting, getting the OF, and optimizing the VT (ie, replacing VCs among the harvest to improve OF)
#' @export
VT2 <- R6Class("VT2", list(

  drugs = NULL,
  filter_data = NULL,
  curve_sampling = NULL,
  harvest = NULL,
  celltheque = NULL,
  cellthequespread = NULL,
  name = NULL,
  pen_nbag = NULL,
  pen_max = NULL,
  cell_by_cell_last = NULL,

  initialize = function(filter_data, name, pen_nbag = 0, pen_max = 0, prev_VT = NULL){
# print("what")

    if(!is.null(prev_VT)){

      self$drugs <- prev_VT2$drugs
      self$filter_data <- prev_VT2$filter_data
      self$curve_sampling <- prev_VT2$curve_sampling
      self$harvest <- prev_VT2$harvest
      self$celltheque <- prev_VT2$celltheque
      self$cellthequespread <- prev_VT2$cellthequespread
      self$name <- prev_VT2$name
      self$pen_nbag <- prev_VT2$pen_nbag
      self$pen_max <- prev_VT2$pen_max
      self$cell_by_cell_last <- prev_VT2$cell_by_cell_last


    }else{

    filter_data <- enexpr(filter_data)

    # first guess drugs
    conc_name <- names(data_VT)[grepl("^conc", names(data_VT))]

    data_VT %>%
      filter(!!filter_data) %>%
      select(starts_with("conc")) %>%
      imap(function(x,y){

        x[x!=0] <- y
        x[x == "0"] <- ""
        x
      }) %>%
      as.tibble() %>%
      mutate(all = paste0(!!!parse_exprs(conc_name))) %>%
      distinct(all) %>%
      filter(all !="") %>%
      mutate(drugs = map(all, ~as.double( str_split(.x, "conc")[[1]][-1]))) %>%
      mutate(ndrug = map_dbl(drugs, ~length(.x))) %>%
      arrange(desc(ndrug)) -> drugs


    drug <- drugs$drugs[which(drugs$ndrug == max(drugs$ndrug))]

    if(length(drug) == 1) drug <- drug[[1]]
    # analyses drugs

    print(paste0("We found the following drugs: ", paste0(drug, collapse = ", "), " (according celltheque will be loaded)"))

    # self$name <- name
    self$drugs <- drug


    self$filter_data <- deparse(filter_data)

    drugvec <- unique(reduce(drug, c))
    drugvec <- drugvec[order(drugvec)]

    if(length(drugvec) == 2){

      if(drugvec[[1]] == 1 & drugvec[[2]] == 4){

        # celltheques <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/new_celltheque_QSPVP/one_per_bag_combined/Veneto.RDS")
        celltheques <- readRDS(here('Virtual_cells', 'one_per_bag_combined', 'Veneto.RDS'))
      }else if( drugvec[[1]] == 2 & drugvec[[2]] == 4){

        # celltheques <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/new_celltheque_QSPVP/one_per_bag_combined/A11.RDS")
        celltheques <- readRDS(here('Virtual_cells', 'one_per_bag_combined', 'A11.RDS'))
        
      }else {

        stop("Two drug but not the combination accepted")

      }
    }else if (length(drugvec) == 3){

      # celltheques <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/new_celltheque_QSPVP/one_per_bag_combined/both.RDS")

      celltheques <- readRDS(here('Virtual_cells', 'one_per_bag_combined', 'both.RDS'))
      
    }
    celltheques <-  celltheques %>%
                      rowid_to_column("cellid")#cellthequeLoad(drug = drugvec, return = 3)


    self$curve_sampling <- curve_sampling( celltheque =celltheques, filter = !!filter_data,  drug = drugvec)

    harvest <-  cell_harvest(sample_df =    self$curve_sampling)

    if(length(  harvest) > 100)   harvest <- sample(harvest, size = 100)

    self$harvest <- harvest


    #
    self$celltheque <- celltheques
    # self$cellthequespread <-  NULL#celltheques[[2]] %>% ungroup
    self$name <- name
    self$pen_nbag <- pen_nbag
    self$pen_max <- pen_max

    }
  }

))






# Reconstitute ----------------------------------------------------------------

VT2$set("public", "reconstitute", function(drug = NULL,  xlog = T,  wrap = T,
                                            comparisonData = T , aera = T, accesRes = F , harvestAlt = NULL){




  if(is.null(drug)) drug <- self$drugs
  if(is.list(drug)){

    return(
      self$drugs %>%
        map(~ self$reconstitute( drug = .x,
                         comparisonData = comparisonData , aera = aera, harvestAlt = harvestAlt) %>%
              mutate(drugs = paste0(.x, collapse = "_"))) %>%
        bind_rows()

    )

  }


  conc1 <- parse_expr(paste0("conc", drug[[1]]))
  proj_conc1 <- unique(data_VT[deparse(conc1)])
  if(is.null(harvestAlt)){

    harvest <- self$harvest
  }else{

    harvest <- harvestAlt
  }
celltheque <- self$celltheque

  if(length(drug) == 1){

    conc2 <- expr(conc2)
    proj_conc2 <-  0
    celltheque <- celltheque %>%
      mutate(conc2 = 0)

  }else{

    conc2 <- parse_expr(paste0("conc", drug[[2]]))
    proj_conc2 <- unique(data_VT[deparse(conc2)])
  }





 # Create all combination of drug we want !
  temp  <-   tibble(cellid = harvest) %>%
    mutate(conc1 = list(proj_conc1)) %>%
    unnest(conc1) %>%
    mutate(concaaaa = list(proj_conc2)) %>%
    unnest(concaaaa)

  # need to take conc==0 for all other concs

  # allconccelltheque <- names(celltheque)[grepl("conc", names(celltheque))]
  #
  # for(a in allconccelltheque[!allconccelltheque %in% paste0("conc", drug)]){
  #
  #   temp <- temp %>%
  #     mutate(!!a := 0)
  #
  # }

  if(drug[[1]] == 1){

tomerge <-     celltheque %>%
      ungroup() %>%
      select(cellid, Veneto_0, Veneto_5, Veneto_10, Veneto_15) %>%
      gather(-cellid, key = "conc4", value = "FirstDeath") %>%
      mutate(conc4 = as.double(gsub("Veneto_","", conc4)))

temp <- temp %>%
    left_join(tomerge, by =  c("cellid", "conc4")) %>%
  mutate(res = if_else(conc1 >=FirstDeath, T, F)) %>%
  group_by(conc1, conc4  ) %>%
  summarise(res = mean (!res) *  100, .groups = 'drop')
  # group_by(res) %>% tally

allconccelltheque <- exprs(conc1, conc4)
  }else{



    tomerge <-     celltheque %>%
      ungroup() %>%
      select(cellid, A11_0, A11_5, A11_10, A11_15) %>%
      gather(-cellid, key = "conc4", value = "FirstDeath") %>%
      mutate(conc4 = as.double(gsub("A11_","", conc4)))

    temp <-   temp %>%
      left_join(tomerge, by = c("cellid", "conc4")) %>%
      mutate(res = if_else(conc2 >=FirstDeath, T, F)) %>%
      group_by(conc2, conc4  ) %>%
      summarise(res = mean (!res) *  100, .groups = 'drop')
    # group_by(res) %>% tally
    allconccelltheque <- exprs(conc2, conc4)
  }

  # temp
  #
  #
  #
  #
  # temp <- temp %>%
  #   left_join(celltheque, by = c("cellid",allconccelltheque)) %>%
  #   group_by(!!!parse_exprs(allconccelltheque)) %>%
  #   summarise(res = mean (!res) *  100, .groups = 'drop')  %>%
  #   mutate(group = paste0(!!conc2)) #%>%
  # mutate(Drug = if_else(Drug == 2, "A-1155463", "Venetoclax"))


    data_VT %>%
      filter(!!conc1 > 0) %>%
      pull(ID) %>%
      unique() -> IDpossible



    filter <- parse_expr(paste0("(", self$filter_data, ") & ID %in% IDpossible"))

    dataobs <- data_VT %>%  #%>% filter(conc2 %in% unique(harvest$conc2))
      filter(!!filter) #%>%
    # mutate(Drug = "Drug1")#○%>%
    # mutate(Value = if_else(Value >100, 100, Value))


    forarea <- temp %>%
      ungroup() %>%
      left_join(dataobs %>% select(!!!allconccelltheque, Value), by = c(deparse(conc1), "conc4")) %>%
      gather(res, Value, key = "Reconst", value = "res")

    return(forarea)


})
# Plotting ----------------------------------------------------------------

VT2$set("public", "plot", function(drug = NULL,  xlog = T,  wrap = T,
                                  comparisonData = T , aera = T, accesRes = F , harvestAlt = NULL){

  # need to analyse the curve to see the subplots, and use the function recursively !

  if(is.null(drug)) drug <- self$drugs

  if(is.list(drug)){

    return(
      self$drugs %>%
        map(~ self$plot( drug = .x, xlog = xlog,  wrap = wrap,
                         comparisonData = comparisonData , aera = aera, harvestAlt = harvestAlt)+
              guides(col = F)+
              ggtitle(paste0("Drug: ",paste0(.x, collapse = " & ")))) %>%
        invoke(.fn = plot_grid)

    )

  }

  # case for one group

  celltheque <- self$celltheque
  if(!is.null(harvestAlt)){
    harvest  <- harvestAlt
  }else{
    harvest <- self$harvest
  }

  if(is.null(drug)) drug <- self$drugs



data <- self$reconstitute()


conc1 <- parse_expr(paste0("conc", drug[[1]]))
proj_conc1 <- unique(data_VT[deparse(conc1)])

if(length(drug) == 1){

  conc2 <- expr(conc2)
  proj_conc2 <-  0
  celltheque <- celltheque %>%
    mutate(conc2 = 0)

}else{

  conc2 <- parse_expr(paste0("conc", drug[[2]]))
  proj_conc2 <- unique(data_VT[deparse(conc2)])
}
# data %>%
  # ggplot()+
  temp <- data %>%
    mutate(Reconst = if_else(Reconst == "res", " Recon-\nstructed", "Observed")) %>%
    mutate(Drug = if_else(drug[[1]] == 2, "A-1155463", "Venetoclax")) %>%
    mutate(group = "A-12")

    #%>%
    # mutate(group = paste0(group, Reconst))

  colexpr <- expr(Reconst)

  plotoutput <- temp %>%
  distinct() %>%
  ggplot()+
  geom_line(aes(!!conc1, res, group = Reconst, col = !!colexpr))+
  geom_point(aes(!!conc1, res, group = Reconst, col = !!colexpr))+
  theme_bw()

# title <- Cellname
if(aera == T & comparisonData == T){


  ribbondata <- data %>%
    spread(Reconst, res)

  OF <- self$OF(detail = T,   drug = drug)



  plotoutput <- plotoutput +
    geom_text(data = OF %>% filter(!is.na(!!conc2)), aes(x = 1, y = 100, label = paste0("Area = ",round(area)))) +
    geom_ribbon(data = ribbondata, aes(x = !!conc1, ymin = Value, ymax = res), alpha = 0.3)
}


plotoutput <- plotoutput + ggtitle(self$name)

if(wrap == T) plotoutput <-  plotoutput+
  facet_wrap(~conc4)

if(xlog == T) plotoutput <- plotoutput + scale_x_log10()
return(plotoutput + labs(col = "", x = "Concentration (µM)", y = "Cell viability (%)")+ #, title = title
         theme(plot.title  = element_text(hjust = 0.5)))




})


# OF Function -------------------------------------------------------

drug = NULL; filtre = NULL; detail = F; perconc1 = F; harvestAlt = NULL; pen_nbag = NULL; pen_delta = NULL
VT2$set("public", "OF", function(drug = NULL, filtre = NULL, detail = F, perconc1 = F, harvestAlt = NULL, pen_nbag = NULL, pen_delta = NULL){

  if(is.null(pen_nbag)) pen_nbag <-self$pen_nbag
  if(is.null(pen_delta)) pen_delta <-self$pen_max


  if(is.null(drug)) drug <- self$drugs

  if(is.list(drug)){

    # Do it for each combination of drugs
    drug %>%
      map(~ self$OF(drug = .x,  detail = T,
                    harvestAlt = harvestAlt , pen_nbag = pen_nbag, pen_delta = pen_delta)) -> tempmultidrug

    # get only the highest first penalty terms (max Delta)
    penaltyterms <- map_dbl(tempmultidrug, ~ .x %>% slice((nrow(.x) - 1)) %>% pull(area))

    maxpenalty <- max(penaltyterms)

    return(

      maxpenalty +  map_dbl(tempmultidrug , ~ .x %>%
                slice(-( nrow(.x) - 1)) %>%
                summarise(sum = sum(area)) %>%
                pull(sum)) %>%
        sum

    )

  }


if(is.null(harvestAlt)){

  harvest <- self$harvest
}else{

  harvest <- harvestAlt
}

  Cellline =self$cell_line


  conc1expr <- parse_expr(paste0("conc", drug[[1]]))
  # proj_conc1 <- unique(data_VT[deparse(conc1)])

  if(length(drug) == 1){

    conc2expr <- expr(conc2)
    celltheque <- celltheque %>%
      mutate(conc2 = 0)
    # proj_conc2 <-  0

  }else{

    conc2expr <- parse_expr(paste0("conc", drug[[2]]))
    # conc2expr <- unique(data_VT[deparse(conc2)])
  }



  #
  #   if(perconc1 == F){
  #
  #     conc1expr <- expr(conc1)
  #     conc2expr <- expr(conc2)
  #   }else{
  #
  #     conc1expr <- expr(conc2)
  #     conc2expr <- expr(conc1)
  #   }

  filtre <- enexpr(filtre)

  reconst <- self$reconstitute(drug = drug, harvestAlt = harvest)


  reconst %>%
    distinct() %>%
    spread(Reconst, res) %>%
    ungroup() %>%
    # filter(conc2==0) %>%
    # group_by(conc2) %>%
    mutate(conc = log(!!conc1expr)) %>%
    mutate(conc = if_else(conc == - Inf, log(0.04), conc)) %>%
    group_by(!!conc2expr) %>%
    nest() %>%
    mutate(area = map_dbl(data, function(x){

      x %>%
        mutate( conclag = lag(conc), reslag = lag(res), Valuelag = lag(Value)) %>%
        slice(-1) %>%
        mutate(area = pmap_dbl(list(Value,Valuelag, res, reslag, conc, conclag), function(Value,Valuelag, res, reslag, conc, conclag){

          # Value = 17.25; Valuelag =  46.69;  res = 18; reslag = 47; conc1 = 2.9957323; conc1lag = 2.302585

          h <- conc - conclag
          a <- sqrt(h ^2 + abs(reslag-Valuelag)^ 2)
          b <- sqrt(h ^2 + abs(res -Value)^ 2)
          area <- (a + b) * h /2

          if((Value > res) != (Valuelag > reslag)) area <- area / 2
          area
        })) %>%
        pull(area) %>%
        sum

    })) %>%
    select(-data) -> temp

  # adding penalties terms
  maxdelta <-
    reconst %>%
      spread(Reconst, res) %>%
      mutate(diff = abs(res - Value)) %>%
      arrange(desc(diff)) %>%
      slice(1) %>%
      pull(diff)



  temp <- temp %>%
    ungroup %>%
    add_row(!!conc2expr := NA, area =   maxdelta * pen_nbag) %>%
    add_row(!!conc2expr := NA, area =   max(temp$area)  * pen_delta)




  if(!is.null(filtre)) temp <- temp %>% filter(!!filtre)

  if(detail == T) return(temp)

  temp %>% pull(area) %>% sum
})



# Distrib -----------------------------------------------------------------

VT2$set("public", "distrib", function() {

  harvest <- self$harvest
  celltheque <- self$celltheque

  celltheque %>%
    select(-starts_with("cellid"), -rowid, -source,-from, -starts_with("conc"), -res) %>%
    names -> characcell

  tibble(cellid = harvest) %>%
    left_join(celltheque %>% distinct(cellid, !!!parse_exprs(characcell))) -> temp

  temp %>%
    gather(-cellid, key = "key", value = "value") %>%
    distinct(key, value) %>%
    group_by(key) %>%
    tally %>%
    filter(n == 1) %>%
    pull(key) ->torem

  torem <- c(torem, "cellid")

  psych::pairs.panels(temp[!names(temp) %in% torem],
                      method = "pearson", # correlation method
                      hist.col = "#00AFBB",
                      density = TRUE,  # show density plots
                      ellipses = FALSE # show correlation ellipses
  )

})




# ModuleDiffnn ------------------------------------------------------------


VT2$set("public", "moduleDiffn", function(){

  harvest <- self$harvest
  curvDecompo <- self$curve_sampling

  curvDecompo %>%
    mutate(n = n * length(harvest)/100) %>%
    # filter(conc2 >= 0 )  %>%
    mutate(neff = map_dbl(sample, ~ sum(harvest %in% .x))) %>%
    mutate(ndif = neff-n) %>%
    select(n, sample, neff, ndif,ID)-> temp

  # see doc ?sample
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  # those in addition
  temp %>%
    filter(ndif >= 0) %>%
    arrange(desc(ndif)) %>%
    ungroup() %>%
    slice(1:3) %>% # take the three lines with the most excess
    # sample_n(3) %>%
    # slice(1) %>% pull(sample) -> x; x <- x[[1]]
    mutate(sampletorem = map2(sample, ndif,  function(x,y){

      temp2 <- harvest[harvest %in% x]
      resample(temp2, size = max(y, length(temp2)))

    })) %>%
    pull(sampletorem) %>% reduce(c) -> to_remove #9


  # Those missing
  temp %>%
    filter(ndif < 0) %>%
    arrange(ndif) %>%
    select(neff, ndif, everything()) %>%
    ungroup() %>%
    slice(1:3) %>% # take the three lines with the most excess
    # sample_n(3) %>%
    mutate(sampletoadd = map2(sample, ndif, function(x,z){

      temp <- x
      # print(z)
      resample(x, size = -z, replace = T)
    })) %>%
    pull(sampletoadd) %>% reduce(c) -> to_add # those to add


  # to_remove <- unique(to_remove) #otherwise too many

  to_remove <- resample(to_remove)
  to_add <- resample(to_add)

  to_modify <-   min(c(length(to_remove), length(to_add)))
  to_modify <- min(to_modify, 20)
  to_modify <- sample(to_modify, 1)

  to_remove <- to_remove[1:to_modify]
  to_add <- to_add[1:to_modify]

## difficulty: remove only one if they are several to remov
  candidate <- harvest


for(a in to_remove) if(a %in% candidate) candidate <- candidate[- (which(candidate == a)[1])]

  to_add <- to_add[1:(100-length(candidate))]
  # candidate <- harvest[- which(harvest %in% to_remove)[1:to_modify] ]

  candidate <- c(candidate, to_add)


  # candidate <- harvest[- sample(which(harvest %in% to_remove),sample(length(to_remove), 1))]

  # candidate <- c(candidate, sample(to_add, sample(length(to_add), 1), replace = T))

  candidate

})

#
# test <- c(3,3,3,4)
#
# torem <- c(3, 3)
#
# for(a in )
#
# test[ - which(test %in% torem)]

# Optim -------------------------------------------------------------------


VT2$set("public", "optim", function( nmax = Inf, pen_nbag = NULL, pen_delta = NULL, drug = NULL){




  harvest <- self$harvest
  celltheque <- self$celltheque
  if(is.null(drug)) drug <- self$drugs
  # drug = list(c(1,4),c(2,4))#
  Cellline = self$cell_line
  if(is.null(pen_nbag)) pen_nbag <- self$pen_nbag
  if(is.null(pen_delta)) pen_delta <- self$pen_max


  OF <- self$OF(detail = F,  pen_nbag = pen_nbag, pen_delta = pen_delta, drug = drug)
  print(paste0("Initial OF = ", OF))

  n <- 0
  # print("here")

  while(n < nmax){ #
    # print(n)

    try({

      harvest_candidate <- self$moduleDiffn()
      OF_candidate <-  self$OF( harvestAlt = harvest_candidate,detail = F,
                                pen_nbag = pen_nbag, pen_delta = pen_delta , drug = drug)

      if(OF_candidate < OF){

        print("############################## NEW CHAMPION ##############################")
        print(paste0("New OF: ", round(OF_candidate,1), " - Diff = ",  round(OF - OF_candidate, 1 ), " - ncells = ", length(harvest_candidate)))

        OF <- OF_candidate
        harvest <- harvest_candidate
        self$harvest <- harvest_candidate

        # eval(expr((!!self_name)$harvest <<- harvest_candidate))
      }
    })

    n <- n +1

  }
})


# optim2 ------------------------------------------------------------------



VT2$set("public", "optim2", function(   pen_nbag = NULL, pen_delta = NULL , drug = NULL, saveEveryXmin = NULL, pathSave = NULL){

  library(progress)


  time0 <- Sys.time()

  harvest_to_beat <- self$harvest
  celltheque <- self$celltheque
  if(is.null(drug)) drug = self$drugs
  Cellline = self$cell_line
  if(is.null(pen_nbag)) pen_nbag <- self$pen_nbag
  if(is.null(pen_delta)) pen_delta <- self$pen_max
  cellids <- unique(celltheque$cellid)

  todo <- cellids
  if(length(self$cell_by_cell_last) > 0 ) todo <- todo[which(todo == self$cell_by_cell_last):length(todo)]




  harvest_candidate <- harvest_to_beat
  OF <-   self$OF(  detail = F, harvestAlt = harvest_candidate,    pen_nbag = pen_nbag, pen_delta = pen_delta, drug = drug)

  pb <- progress_bar$new(total = length(cellids) * length(harvest_to_beat))
  pb$update(ratio = (length(cellids) - length(todo))/length(cellids)) #

  for(a in todo){

    for(b in harvest_to_beat){

      harvest_candidate <- harvest_to_beat

      harvest_candidate[harvest_candidate == b][[1]] <- a
      OF_candidate <-   self$OF( detail = F, harvestAlt = harvest_candidate,   pen_nbag = pen_nbag, pen_delta = pen_delta, drug = drug)
      # print(OF)
      # print(OF_candidate)
      if(OF_candidate < OF){

        print("############################## NEW CHAMPION ##############################")
        print(paste0("New OF: ", round(OF_candidate,1), " - Diff = ",  round(OF - OF_candidate, 1 ), " - ncells = ", length(harvest_candidate)))

        OF <- OF_candidate
        harvest_to_beat <-  harvest_candidate

        self$harvest <- harvest_candidate

      }

      pb$tick()
      self$cell_by_cell_last <- a

    }


if(!is.null(saveEveryXmin) & !is.null(pathSave)){

  deltaTime <- difftime( Sys.time() ,time0, units = "min") %>% as.double()

if(deltaTime >saveEveryXmin){

  saveRDS(self, pathSave)
  time0 <- Sys.time()
}
} # end saving system

    }
})



# Save --------------------------------------------------------------------



VT2$set("public", "save", function(){


  saveRDS(self, file.path(active_VT_project,  "3_virtual_tumors", paste0(self$name, ".RDS")))

})



# extract_bags ------------------------------------------------------------


VT2$set("public", "extract_bags", function(){


  drug <- self$drugs
  harvest <- self$harvest
  harvestspread <- self$cellthequespread %>%
    filter(cellid %in% harvest)



  stpreadnew <- NULL

  root <- file.path(active_VT_project, "2_celltheques")

  pathspread <- paste0(root, "/celltheque_one_per_cell_spread_drug_", paste0(drug, collapse = "_"))

  listspread <- list.files(pathspread)

  for(a in listspread){


    # need to create and save the full spread

    print(a)



    temp_spread <- load_spread(path_celltheque = a, drug = drug,returnOnIDperGroup = F, update = F ) %>%
      mutate(from = a)

    if(is.null(stpreadnew )){
      stpreadnew <- temp_spread
    }else{
      stpreadnew <- bind_rows(stpreadnew, temp_spread )
    }
  }


  # Here cellid are based on source, but we need to reattribute new cell id
  stpreadnew <- stpreadnew %>%
    rename(cellidfromsource= cellid) %>%
    rowid_to_column(var = "cellidnew")

  stpreadnew %>%
    left_join(

      # we want to remove all profile not present in our virtual tumor
      # with a left_join and removing the NA
      harvestspread %>%
        select(cellid, starts_with("conc")) %>%
        rename(cellid_in_harvest = cellid)
    ) %>%
    filter(!is.na(cellid_in_harvest)) -> newtemp


  #` some cells are here twice, both with harvestspread and newtemp...`
  # harvestpreadnew <- bind_rows(harvestspread ,
  #                              newtemp %>%
  #                                mutate(cellid =  cellidnew)) %>%
  #   select(-n, -starts_with("cellid0"))
  # %>   mutate(cellid =  cellidnew+ max(harvestspread$cellid))

  names_temp <-   names(newtemp)
  names_temp <-    names_temp[grepl("conc", names(newtemp))]

  names_temp <- paste0("`", names_temp, "`")


  bags <- newtemp %>%
    group_by( !!!parse_exprs(names_temp)) %>%
    nest() %>%
    select(data, everything()) %>%
    rowid_to_column(var = "bag") %>%
    unnest() %>%
    ungroup()


  # bags %>% select(cellidfromsource, from)
  # bags <- map(bags, ~ .x %>%    distinct(cellid,  Bak ,  Bax,  Bcl2, Bclxl , Mcl1  , eta, from ))

  return(bags)
  # saveRDS(self, paste0("D:/these/Second_project/QSP/modeling_work/Virtual_Tumor/", self$name, ".RDS"))

})






# tibble ------------------------------------------------------------------



VT2$set("public", "tibble", function(addconc = F, drug = NULL){

  temp <- tibble(cellid = VT2$harvest) %>%
    left_join(VT2$cellthequespread) %>%
    select(cellid, Bak, Bax, Bcl2, Bclxl, Mcl1, contains("BIM"), contains("PUMA"), contains("NOXA"), contains("eta"))

  if(addconc == T){

    temp %>%
      mutate(conc1 = list(proj_conc1)) %>%
      unnest() %>%
      mutate(conc2 = list(proj_conc2)) %>%
      unnest() -> temp

  }

  if(!is.null(drug)){

    temp %>%
      mutate(drug = list(drug)) %>%
      unnest -> temp

  }

  return(temp)
})




# pdf ---------------------------------------------------------------------




VT2$set("public", "pdf",  function(){

  celltheque <- self$celltheque

  harvest <- self$harvest

  bags <- self$extract_bags()

  path <- file.path(active_VT_project,  "3_virtual_tumors", paste0(self$name, ".pdf"))



  pdf(path, width = 12, height = 10)
  print(self$plot())
  print(self$plot(comparisonData = F, wrap  = F ))



  tibble(a = harvest) %>%
    group_by(a) %>%
    tally %>%
    arrange(desc(n)) %>%
    # slice(1:2) %>% # just for trying the formula
    pull(a) %>%
    map(~ print(self$plots_cell(cellID = .x,bags = bags)))
  dev.off()
  shell.exec(path)
})





# plots_cell --------------------------------------------------------------



VT2$set("public", "plots_cell", function( cellID , bags = NULL){



  celltheque <- self$celltheque

  harvest <- self$harvest

  drug <- self$drugs
  # cellID <- harvest[[20]]

  conc1 <- parse_expr(paste0("conc", drug[[1]]))

  if(length(drug) == 1){

    conc2 <- expr(conc2)


  }else{

    conc2 <- parse_expr(paste0("conc", drug[[2]]))

  }


  plot1 <-  matrix_cell(celltheque = celltheque, cellID = cellID, drug = drug)

  ## plot 2: add points to reconstitution plots

  celltheque %>%
    filter(cellid == cellID) %>%
    group_by(!!conc2) %>%
    filter(res == T) %>%
    slice(1) %>%
    left_join(
      ## add the pct of viability at according concentrations of reproduction
      self$plot(comparisonData = F, wrap  = F, accesRes = T) %>%
        select(!!conc1, !!conc2, res) %>% rename(y = res)

    )  -> temp_for_points



  plot2 <- self$plot( comparisonData = F, wrap  = F ) +
    geom_point(data = temp_for_points, aes(x = !!conc1, y = y), fill = "red", shape = 21, size = 3)

  # get the bag
  if(is.null(bags)) bags <-  self$extract_bags()

  bags %>%
    filter(cellid_in_harvest == cellID) %>% pull(bag) %>% unique() -> bagn

  bag <- bags %>%
    filter(bag == bagn & !is.na(from))


  # plot 3: distributions inside the VT


  celltheque %>%
    select(-starts_with("cellid"), -rowid, -source,-from, -starts_with("conc"), -res) %>%
    names -> characcell

  tibble(cellid = harvest) %>%
    left_join(celltheque %>% distinct(cellid, !!!parse_exprs(characcell))) -> temp

  temp %>%
    gather(-cellid, key = "key", value = "value") %>%
    distinct(key, value) %>%
    group_by(key) %>%
    tally %>%
    filter(n == 1) %>%
    pull(key) ->torem

  temp <- temp[ , ! names(temp) %in% torem]

  temp <- temp %>%
    select(-starts_with("Drug"), - starts_with("group"))

  bagspread <- bag %>%
    rename(cellid = cellid_in_harvest) %>%
    select(!!!parse_exprs(names(temp))) %>%
    # slice(1) %>%
    gather("key", "value", -cellid)

  plot3 <- temp %>%
    # select(-!!!parse_expr(torem))
    gather("key", "value", -cellid) %>%
    ggplot()+
    geom_vline( data = celltheque %>%
                  filter(cellid == cellID) %>%
                  select(!!!parse_exprs(names(temp))) %>%
                  slice(1) %>%
                  gather("key", "value", -cellid),
                aes(xintercept = value), size = 2, col = "red", alpha = 0.5)+

    geom_histogram(aes(value))+
    facet_wrap(~key, scales = "free")+
    theme_bw()+
    ggtitle("Cell place inside the VT")+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_rect(data = bagspread %>%
                group_by(key) %>%
                summarise(min = min(value), max = max(value)),
              aes(xmin = min, xmax = max, ymin = 0, ymax = Inf), fill = "red", alpha = 0.1)

  # Fig 4

  plot4 <-  bagspread %>%
    ggplot() +
    geom_histogram(aes(value))+
    facet_wrap(~key, scales = "free")+
    theme_bw()+
    geom_vline( data = celltheque %>%
                  filter(cellid == cellID) %>%
                  select(!!!parse_exprs(names(temp))) %>%
                  slice(1) %>%
                  gather("key", "value", -cellid),
                aes(xintercept = value), size = 2, col = "red", alpha = 0.5)+
    ggtitle(paste0(nrow(bag), " cells with same profile"))+
    theme(plot.title = element_text(hjust = 0.5))

  plot_grid(plot1, plot2, plot3, plot4)

})



# description -------------------------------------------------------------




VT2$set("public", "description",  function(returnAna = F, drug = NULL){

  celltheque <- self$celltheque

  harvest <- self$harvest

  if(is.null(drug)) drug <- self$drugs



  if(is.list(drug)){


   return( drug %>%
      map(~   self$description(drug = .x, returnAna = returnAna))
   )


  }

  if("group" %in% names(celltheque)){

    celltheque <- celltheque %>%
      filter(group == paste0("Drug", paste0(drug, collapse = "_")))

  }
  # cellID <- harvest[[20]]

  conc1 <- parse_expr(paste0("conc", drug[[1]]))

  if(length(drug) == 1){

    conc2 <- expr(conc2)


  }else{

    conc2 <- parse_expr(paste0("conc", drug[[2]]))

  }

  ## now let's compute the analysis


  tibble(id = unique(harvest))  %>%
    mutate(analyse = map(id, function(id){
      # print(id)
      #empty tibble to fill for futur unnest
      results <- tibble(a = NA)

      # look at first drug
      matrix_drug1 <- matrix_cell(celltheque = celltheque, id, drug = drug, plot = F)
      drug1 <- sum(matrix_drug1$`0`)
      nconc1 <- nrow(matrix_drug1)

      results$drug1 <- case_when(drug1 == nconc1 ~ "Native_dead",
                                 drug1 >= nconc1 / 2 ~ "High_sensitiv",
                                 drug1 > 0 ~ "Sensitiv",
                                 drug1 == 0 ~ "Resistant"
      )
      results$drug1Level <- nconc1 - drug1 + 1
      if(results$drug1Level == nconc1 + 1) results$drug1Level <- Inf
      # look at second drug
      if(length(drug) == 2){

        matrix_drug2 <- matrix_cell(celltheque = celltheque, id, drug = c(drug[[2]], drug[[1]]), plot = F)
        drug2 <- sum(matrix_drug2$`0`)
        nconc2 <- nrow(matrix_drug2)

        results$drug2 <- case_when(drug2 == nconc2 ~ "Native_dead",
                                   drug2 >= nconc2 / 2 ~ "High_sensitiv",
                                   drug2 > 0 ~ "Sensitiv",
                                   drug2 == 0 ~ "Resistant"
        )
        results$drug2Level <- nconc2 - drug2 + 1
        if(results$drug2Level == nconc2 + 1) results$drug2Level <- Inf
      }
      # Finally, the total

      results$total <- map_dbl(matrix_drug1[,-1 ], ~sum(.x )) %>% sum

      results %>% select(-a)


    })) %>%
    unnest() -> analyse



  dimm <- dim(matrix_cell(celltheque = celltheque, harvest[[1]], drug = drug, plot = F))
  nconc <- dimm[1] * (dimm[2]-1)

  tibble(id = harvest) %>%
    left_join(analyse) %>%
    mutate(profile = case_when(drug1 %in% c("High_sensitiv", "Sensitiv") & drug2 == "Resistant" ~ "Drug1_dependant",
                               drug2 %in% c("High_sensitiv", "Sensitiv") & drug1 == "Resistant" ~ "Drug2_dependant",
                               drug1  %in% c("High_sensitiv", "Sensitiv") & drug2 %in% c("High_sensitiv", "Sensitiv") ~ "Both",
                               # drug1 == "Resistant" & drug2 == "Resistant" & total < nconc/4 ~
                               total == 0 ~ "Unkilable",
                               total == nconc ~ "Die anyway",
                               drug1 == "Resistant" &drug2 == "Resistant"~ "Synergism",
                               T ~ "hm?"
    )) -> comput
  if(returnAna == T) return(comput)

  comput %>%
    group_by(profile) %>%
    tally

})




# prot_expr ---------------------------------------------------------------



VT2$set("public", "prot_expr", function(return_median_i = F){

  celltheque <- self$celltheque

  harvest <- self$harvest

  drug <- self$drugs

  if("group" %in% names(celltheque)){

    celltheque <- celltheque %>%
      filter(group == paste0("Drug", paste0(drug, collapse = "_")))

  }
  # cellID <- harvest[[20]]

  conc1 <- parse_expr(paste0("conc", drug[[1]]))

  if(length(drug) == 1){

    conc2 <- expr(conc2)


  }else{

    conc2 <- parse_expr(paste0("conc", drug[[2]]))

  }

  bags <-  self$extract_bags()

  # bags_median

  bags %>%
    select(-starts_with("conc"), -"cellidnew", -"cellidfromsource", -"from") %>%
    distinct() %>%
    gather("prot", "value", -cellid_in_harvest, -bag) %>%
    group_by(bag, prot,cellid_in_harvest) %>%
    summarise(median = median(value)) %>%
    spread(key = prot, value = median) -> median_i

  if(return_median_i == T) return(median_i)

  tibble(cellid_in_harvest = harvest) %>%
    # left_join(bags %>% filter(is.na(from)) %>% select(cellid, bag)) %>%
    left_join(median_i) %>%
    gather("prot", "value", -cellid_in_harvest, -bag) %>%
    group_by(prot) %>%
    summarise(median = median(value)) %>%
    spread(key = prot, value = median)

})





# Load Virtual Tumor ------------------------------------------------------

# Define VT ---------------------------------------------------------------
#'
#' @export
load_VT <- function(){

  path <- file.path(active_VT_project, "3_virtual_tumors")

  files <- list.files(path)

  files <- files[grepl(".rds$", tolower(files))]

  files <- data.frame(files = files)

  print(files)

  choice <- readline(prompt = "Give the number of the virtual tumor you want: ")

  pathRDS <- files %>%
    slice(as.double(choice)) %>%
    pull(files)



  return(readRDS(file.path(path, pathRDS)))
}

# Test --------------------------------------------------------------------



#
# VT2 <- VT2$new(filter_data = Drug %in% c(1,2)  & Cell_line_bin == 1, name = "SUDHL4-three_drugs")
# VT2$OF()
# VT2$plot()
# VT2$moduleDiffn()
# VT2$optim()
# VT2$optim2()
