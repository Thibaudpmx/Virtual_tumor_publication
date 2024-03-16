library(R6)
library(tidyr)
library(dplyr)
library(purrr)
library(rlang)
library(here)
library(cowplot)
library(RxODE)
library(progress)

# Author: Thibaud Derippe
# source("D:/these/Second_project/QSP/VirtualTumor/R/Virtual_Tumor2.R")
source(here('0_Lindner_model_PaSM_config.r'))

# Load the cell line 1 virtual tumor (calibrated with 3 or 2 drug)
VT3d <- readRDS(here('calibrated_VT/VT_both_cell_line_1.RDS'))
VT2d <- readRDS(here('calibrated_VT/VT_Venetoclax_cell_line_1.RDS'))

# Check they are the one fitted
VT3d$plot()
VT2d$plot()

# In VT celltheque, there is only one cell per bag... but in reality we 
# have computed several cells that can belong in each bag..
# so lets retake all the possible cells needed:
allbags3d <- readRDS(here('calibrated_VT/sudh_all_bags.RDS'))
allbags2d <- readRDS(here('calibrated_VT/sudh2d_all_bags.RDS'))

# ABM functions Functions ---------------------------------------------------------------

add_events_veneto <- tibble( cmt = c("Veneto_gut"), time = 700 + seq(0,21*24,24), amt = c(50),evid = 1)

add_events <- add_events_veneto

simulations <- function( VT, add_events, modifParam = character(), random = NULL){

  if(!is.null(random)){

  to_use <-  tibble(cellid = VT$harvest) %>%
      mutate(new = map(cellid, function(x){

        random %>%
          filter(cellid == x) %>%
          sample_n(1) %>%
          select(Bcl20, Bclxl0, Mcl10,  BIM0, PUMA0, NOXA0,  BAK0, BAXc0)


      })) %>%
      unnest() %>%
    select(-cellid) %>%
  rowid_to_column("cellid")
  }else{

    to_use <- tibble(cellid = VT$harvest) %>%
      left_join(VT$celltheque)

    }

  simulations <-
    to_use %>%
    distinct() %>%
    # filter(cellid == 42471) %>% slice(1) %>%
    select(cellid, BAK0, BAXc0, Bcl20, Bclxl0,  Mcl10, BIM0, PUMA0, NOXA0 )%>%
    rowid_to_column("id")%>%
    mutate(!!!parameters_default_values) %>%
    mutate(Mcl1_I0 = 0, Bclxl_I0 = 0, Bcl2_I0=0 ) %>%
    mutate(ratioTumor = 1)

  # modification values
  if(length(modifParam) > 0){
    for(a in 1:length(modifParam)){

      a <- modifParam[a]
      simulations[[names(a)]] <- a
    }
  }

  events <- crossing(id = unique(simulations$id),add_events %>%
                       bind_rows(tibble(cmt = "Veneto_gut", time = 700:800, amt = 0, evid = 0)) %>%
                       arrange(time))

  res <- model_RxODE$solve(simulations, events) %>% as_tibble() %>%
    left_join(simulations %>% distinct(id, cellid))


  # res %>% ggplot()+geom_line(aes(time, Pore, group = id)) +scale_y_log10() + geom_hline(yintercept = 10, col = "red")


  timeperid <- res %>%
    filter(Pore >= 10) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup() %>%
    select(cellid, time)


  tibble(cellid = to_use$cellid) %>%
    left_join(timeperid) %>%
    mutate(time = if_else(is.na(time), Inf, time)) %>%
    rename(timedeath = time) %>%
    mutate(timedeath = (timedeath - 700) / 24)

}

# simple ABM (without senescence)
ABM_simple <- function(simulation, nsimul = 10, timeend = 30){




  pb <- progress::progress_bar$new(clear = FALSE,total = nsimul ,
                                   format = "[:bar] :current/:total (:percent) in :elapsed")

  res_final  <- tibble()

  for(a in 1:nsimul){


    t0 <- simulation %>% select(cellid, timedeath) %>%
      mutate(timedeath0 = timedeath) %>%
      mutate(timedeath = timedeath  + 1)  # +1 because adm after day 1

    tim <- 0
    step <- 0.1
    lambda0 <- 0.131
    probmult <- lambda0 * step

    res <- tibble(time = 0, TV = 100)

    while(tim < timeend){
      tim  <- tim + step


      t0 %>%
        mutate(testMult =  sample(x = c(T, F), prob = c(probmult, 1 - probmult), size = nrow(t0), replace = T)) %>%
        mutate(testDeath =  timedeath < tim) %>%
        filter(testDeath == F) -> temp

      death <- nrow(t0) - nrow(temp)
      # if(death > 0) print(paste0(death," death"))
      # print
      # -> temp

      t0 <- bind_rows(temp, temp %>% filter(testMult == T) %>% mutate(timedeath =  timedeath0 + tim)) # - 1 because we added 1 before
      res <- res %>% bind_rows(tibble(time = tim, TV = nrow(t0)))
    }

    res_final <- bind_rows(res_final, res %>% mutate(simul = a))

    pb$tick()

  }

  res_final
}

#  ABM with senescence
ABM_non_prolif <- function(simulation, nsimul = 10, timeend = 21, pctNP = 0.5, nres = 0){




  pb <- progress::progress_bar$new(clear = FALSE,total = nsimul ,
                                   format = "[:bar] :current/:total (:percent) in :elapsed")

  res_final  <- tibble()

  for(a in 1:nsimul){

    #
    # simulation <-   simulations_veneto
    # pctNP <- 0.3
    t0 <- simulation %>% select(cellid, timedeath) %>%
      mutate(timedeath0 = timedeath) %>%
      mutate(timedeath = timedeath  + 1) %>%   # +1 because adm after day 1
      mutate(NP = sample(c(T,F),prob = c(pctNP, 1 - pctNP),
                         nrow(simulation), replace = T))

    # nres <- 2
    if(nres > 0){

      nInf0 <-  t0 %>% filter(timedeath == Inf) %>% nrow()

      ntodo <- nres -  nInf0

      if(ntodo >0){

        tochange <- which(t0$timedeath != Inf)[1:ntodo]

        t0$timedeath[tochange] <- Inf
        t0$timedeath0[tochange] <- Inf
        t0$NP[tochange] <- F

      }
    }

    tim <- 0
    step <- 0.1
    lambda0 <- 0.131
    probmult <- lambda0 * step

    res <- tibble(time = 0, TV = 100)

    while(tim < 30){
      tim  <- tim + step


      t0 %>%
        mutate(testMult =  sample(x = c(T, F), prob = c(probmult, 1 - probmult), size = nrow(t0), replace = T)) %>%
        mutate(testDeath =  timedeath < tim) %>%
        filter(testDeath == F | NP == T) -> temp

      death <- nrow(t0) - nrow(temp)
      # if(death > 0) print(paste0(death," death"))
      # print
      # -> temp

      t0 <- bind_rows(temp, temp %>% filter(testMult == T & (NP == F | tim < timedeath)) %>% mutate(timedeath =  timedeath0 + tim)) # - 1 because we added 1 before
      res <- res %>% bind_rows(tibble(time = tim, TV = nrow(t0)))
    }
    res_final <- bind_rows(res_final, res %>% mutate(simul = a))

    pb$tick()

  }

  res_final
}

plot_ABM <- function(..., name = character(), obs = 1:4){

  # datas <- exprs(reswithA15, rescontrol); name <- c("A15", "Control")
  datas <- enexprs(...)

  allDatas <- eval(expr(bind_rows(!!!datas, .id = "Group")))

  if(length(name) >0){

    for(a in 1:length(name)){

      newname <- name[[a]]
      allDatas$Group[allDatas$Group == a] <- newname
    }

  }

  # compar <- read.table("D:/these/Second_project/QSP/modeling_work/In_vivo/SU_DHL4_full.csv", header = T, sep = ";", na.strings = ".")
  compar <- read.table(here('data', 'mice_SU_DHL4_full.csv'), header = T, sep = ";", na.strings = ".")
 
  OBS <-
    compar %>% filter(ID %in% c(1:4)) %>%
    mutate(OBS = case_when(ID == 1 ~ "Control",
                           ID == 2 ~ "Venetoclax",
                           ID == 3 ~ "A-15",
                           T ~ "Combo")) %>%
    filter(ID %in% obs)

  allDatas %>%
    group_by(time, Group) %>%
    summarise(q5 = quantile(TV, 0), q95 = quantile(TV, 1) ) %>%
    mutate(q5 = q5 * 2.57, q95 = q95 * 2.57) %>%
    ggplot()+
    geom_ribbon(aes(x = time, ymin = q5, ymax = q95, fill = Group), alpha = 0.5) +
    # geom_line(aes(time, TV* 2.57, col = Dose))+
    scale_y_log10()+
    geom_line(data = OBS,aes(TIME, DV, group = factor(Dose), lty = fct_reorder(OBS, ID)))+
    theme_bw()+
    labs(x = "Time (days)", y = "Tumor Volume", fill = "Simulations", lty = "Observations") +
    scale_linetype_manual(values = c(1,2,3,4)[obs])+
    geom_vline(xintercept = 1, lty = 3)


}

add_events <- add_events_veneto <- tibble( cmt = c("Veneto_gut"), time = 700 + seq(0,21*24,24), amt = c(50),evid = 1)

add_events_A15 <- tibble( cmt = c("A15"), time =  700 +  c(0,48,96, 168, 216, 264 , 336, 384, 432),
                          amt = c(1.5), evid = 1)
add_events_combo <-  bind_rows( add_events_A15, add_events_veneto) %>% arrange(time)

# Control -----------------------------------------------------------------




# control
nsimul <- 10

rescontrol <- ABM_simple(tibble(cellid = 1:100, timedeath = Inf), nsimul)




# simple ABM with Three drugs -------------------------------------------------------------
VT3d$plot()

# Venetoclax
map(1:10, function(x){

  sim <- simulations(VT3d, add_events_veneto, modifParam = c( ratioTumor = 0.3), random = allbags3d)
  ABM_simple(sim %>% mutate(timedeath = timedeath), nsimul ) %>% mutate(nrandom = x)
}) %>%
  bind_rows() -> reswithVeneto3d

reswithVeneto3d %>%
  filter(nrandom <=8) %>%
  group_split(nrandom) %>%

  map(function(x){

    temp <<- x
    print(temp)
    plot_ABM(temp, rescontrol, name = c("Veneto", "Control"), obs = 1:2)

    }) %>%
  invoke(.fn = plot_grid)

simulations_veneto3d <- simulations(VT3d, add_events_veneto, modifParam = c( ratioTumor = 0.3))
reswithVeneto3d <- ABM_simple(simulations_veneto3d %>% mutate(timedeath = timedeath), nsimul )
plotA3d <- plot_ABM(reswithVeneto3d, rescontrol, name = c("Veneto", "Control"), obs = 1:2);plotA3d

# A15
map(1:10, function(x){

  sim <- simulations(VT3d, add_events_A15, modifParam = c( ke_Mcl1_I = 0,   hillA15 = 5,
                                                              EC50A15 = 1E-6), random = allbags3d)
  ABM_simple(sim %>% mutate(timedeath = timedeath), 1 ) %>% mutate(simul = x)
}) %>%
  bind_rows() -> reswithA153d

simulations_A153d <- simulations(VT3d, add_events_A15, modifParam = c( ke_Mcl1_I = 0,   hillA15 = 5,
                                                                       EC50A15 = 1E-6))
reswithA153d  <- ABM_simple(simulations_A153d , nsimul )
plotB3d <- plot_ABM(rescontrol,reswithA153d, name = c(" Control", "A-15"), obs = c(1,3));plotB3d



# Combo
simulations_combo3d <- simulations(VT3d, add_events_combo, modifParam = c(  ratioTumor = 0.3, ke_Mcl1_I = 0,   hillA15 = 5,
                                                                            EC50A15 = 1E-6))
reswithCombo3d <- ABM_simple(simulations_combo3d, nsimul )
plot_ABM(rescontrol,reswithCombo3d, name = c(" Control", "Combo"), obs = c(1,4))
reswithCombo3d_0.4 <- ABM_non_prolif(simulations_combo3d, nsimul, pctNP = 0.5 )
plot_ABM(rescontrol,reswithCombo3d_0.4, name = c(" Control", "Combo"), obs = c(1,4))
# All

plot_ABM(rescontrol,reswithA153d, reswithVeneto3d, reswithCombo3d, name = c(" Control", "A-15", "Veneto", "Combo"))

# compar <- read.table("D:/these/Second_project/QSP/modeling_work/In_vivo/SU_DHL4_full.csv", header = T, sep = ";", na.strings = ".")

compar <- read.table(here('data', 'mice_SU_DHL4_full.csv'), header = T, sep = ";", na.strings = ".")
OBS <-
  compar %>% filter(ID %in% c(1:4)) %>%
  filter(ID != 1) %>%
  mutate(group = case_when(ID == 1 ~ "Control",
                           ID == 2 ~ "Venetoclax",
                           ID == 3 ~ "A-1592668",
                           T ~ "Combo"))
OBS <- OBS %>%
  bind_rows(
    OBS %>%
      filter(ID == 4) %>%
      mutate(group = "Combo 50% senescence")

  )

# bind_rows(
#
#   reswithVeneto3d %>% filter(time <21) %>% mutate(group = "Venetoclax", drug = T),
#   reswithA153d %>% filter(time <21) %>% mutate(group = "A-1592668", drug = T),
#   reswithCombo3d %>% filter(time <29) %>% mutate(group = "Combo", drug = T),
#   reswithCombo3d_0.4  %>% filter(time <29)  %>% mutate(group = "Combo 50% senescence", drug = T),
#   crossing(rescontrol %>% filter(time <15), group = c("Venetoclax","A-1592668", "Combo", "Combo 50% senescence"), drug = F )
# ) %>%
#   filter(!(drug %in% c("A-1592668", "Venetoclax") & time >15)) %>%
#
#   group_by(time, group, drug) %>%
#   summarise(q5 = quantile(TV, 0.05, na.rm = T), q95 = quantile(TV, 0.95, na.rm = T) ) %>%
#   mutate(q5 = q5 * 2.57, q95 = q95 * 2.57) %>%
#   mutate(drug=if_else(drug == T, "Treatment", "Control")) %>%
#
#   ggplot()+
#   geom_ribbon(aes(x = time, ymin = q5, ymax = q95, fill = drug), alpha = 0.5) +
#   # geom_line(aes(time, TV* 2.57, col = Dose))+
#   facet_wrap(~fct_relevel(group, "Venetoclax"), scales = "free")+
#   scale_y_log10()+
#   geom_line(data = OBS,aes(TIME, DV))+
#   geom_line(data = compar %>% filter(ID %in% c(1)) ,aes(TIME, DV))+
#   theme_bw()+
#   labs(x = "Time (days)", y = "Tumor Volume (mm3)", fill = "Group", lty = "Observations") +
#   # scale_linetype_manual(values = c(1,2,3,4)[obs])+
#   geom_vline(xintercept = 1, lty = 3)



bind_rows(

  reswithVeneto3d %>% filter(time <21) %>% mutate(group = "Venetoclax", drug = T),
  reswithA153d %>% filter(time <21) %>% mutate(group = "A-1592668", drug = T),
  reswithCombo3d %>% filter(time <29) %>% mutate(group = "Combo", drug = T),
  reswithCombo3d_0.4  %>% filter(time <29)  %>% mutate(group = "Combo 50% senescence", drug = T),
  crossing(rescontrol %>% filter(time <15), group = c("Venetoclax","A-1592668", "Combo", "Combo 50% senescence"), drug = F )
) %>%
  write.table(file = "ABT3D.csv")


read.table("ABT3D.csv", header = T) %>%
  filter(!(drug %in% c("A-1592668", "Venetoclax") & time >15)) %>%

  group_by(time, group, drug) %>%
  summarise(q5 = quantile(TV, 0, na.rm = T), q95 = quantile(TV, 1, na.rm = T) ) %>%
  mutate(q5 = q5 * 2.57, q95 = q95 * 2.57) %>%
  mutate(drug=if_else(drug == T, "Treatment", "Control")) %>%

  ggplot()+
  geom_ribbon(aes(x = time, ymin = q5, ymax = q95, fill = drug), alpha = 0.5) +
  # geom_line(aes(time, TV* 2.57, col = Dose))+
  facet_wrap(~fct_relevel(group, "Venetoclax"), scales = "free", nrow = 1)+
  # scale_y_log10()+
  geom_line(data = OBS,aes(TIME, DV, lty = "Observations"))+
  geom_line(data = compar %>% filter(ID %in% c(1)) ,aes(TIME, DV))+
  theme_bw()+

  labs(x = "Time (days)", y = "Tumor Volume (mm3)", fill = "ABM simulations", lty = "") +
  # scale_linetype_manual(values = c(1,2,3,4)[obs])+
  geom_vline(xintercept = 1, lty = 3) -> whole3dplot


# simple ABM with Two drugs -------------------------------------------------------------
VT2d <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/new_VT_qspvp/VT_Venetoclax_cell_line_1.RDS")
VT2d$plot()

# Venetoclax
# Venetoclax
map(1:10, function(x){

  sim <- simulations(VT2d, add_events_veneto, modifParam = c( ratioTumor = 0.3), random = allbags2d)
  ABM_simple(sim %>% mutate(timedeath = timedeath), nsimul ) %>% mutate(nrandom = x)
}) %>%
  bind_rows() -> reswithVeneto2d

reswithVeneto2d %>%
  filter(nrandom <=8) %>%
  group_split(nrandom) %>%

  map(function(x){

    temp <<- x
    print(temp)
    plot_ABM(temp, rescontrol, name = c("Veneto", "Control"), obs = 1:2)

  }) %>%
  invoke(.fn = plot_grid)


simulations_veneto2d <- simulations(VT2d, add_events_veneto, modifParam = c( ratioTumor = 0.3))
reswithVeneto2d <- ABM_simple(simulations_veneto2d %>% mutate(timedeath = timedeath), nsimul )
plotA2d <- plot_ABM(reswithVeneto2d, rescontrol, name = c("Veneto", "Control"), obs = 1:2);plotA2d

# reswithVeneto2d05 <- ABM_non_prolif(simulations_veneto2d %>% mutate(timedeath = timedeath), nsimul , pctNP = 0.3)
# plot_ABM(reswithVeneto2d05, rescontrol, name = c("Veneto", "Control"), obs = 1:2)

# A15
simulations_A152d <- simulations(VT2d, add_events_A15, modifParam = c( ke_Mcl1_I = 0,   hillA15 = 5,
                                                                       EC50A15 = 1E-6))
reswithA152d  <- ABM_simple(simulations_A152d , nsimul )
plotB2d <- plot_ABM(rescontrol,reswithA152d, name = c(" Control", "A-15"), obs = c(1,3));plotB2d



# Combo
simulations_combo2d <- simulations(VT2d, add_events_combo, modifParam = c(  ratioTumor = 0.3, ke_Mcl1_I = 0,   hillA15 = 5,
                                                                            EC50A15 = 1E-6))
reswithCombo2d <- ABM_simple(simulations_combo2d, nsimul )
plotC2d <-  plot_ABM(rescontrol,reswithCombo2d, name = c(" Control", "Combo"), obs = c(1,4))

# reswithCombo2d <- ABM_simple(simulations_combo2d, nsimul )

reswithCombo2d_0.4 <- ABM_non_prolif(simulations_combo2d, nsimul, pctNP = 0.5 )

plotD2d <-  plot_ABM(rescontrol,reswithCombo2d_0.4, name = c(" Control", "Combo"), obs = c(1,4));plotD2d
# All

# simple  ABM - Final plot --------------------------------------------------------------


plot_ABM(rescontrol,reswithA152d, reswithVeneto2d, reswithCombo2d, name = c(" Control", "A-15", "Veneto", "Combo"))

plot_grid(plotA2d, plotB2d, plotC2d, plotD2d)

compar <- read.table(here('data', 'mice_SU_DHL4_full.csv'), header = T, sep = ";", na.strings = ".")
  # read.table("D:/these/Second_project/QSP/modeling_work/In_vivo/SU_DHL4_full.csv", header = T, sep = ";", na.strings = ".")

OBS <-
  compar %>% filter(ID %in% c(1:4)) %>%
  filter(ID != 1) %>%
  mutate(group = case_when(ID == 1 ~ "Control",
                         ID == 2 ~ "Venetoclax",
                         ID == 3 ~ "A-1592668",
                         T ~ "Combo"))
OBS <- OBS %>%
  bind_rows(
    OBS %>%
      filter(ID == 4) %>%
      mutate(group = "Combo 50% senescence")

  )

# bind_rows(
#
#   reswithVeneto2d %>% filter(time <21) %>% mutate(group = "Venetoclax", drug = T),
#   reswithA152d %>% filter(time <21) %>% mutate(group = "A-1592668", drug = T),
#   reswithCombo2d %>% filter(time <29) %>% mutate(group = "Combo", drug = T),
#   reswithCombo2d_0.4  %>% filter(time <29)  %>% mutate(group = "Combo 50% senescence", drug = T),
#   crossing(rescontrol %>% filter(time <15), group = c("Venetoclax","A-1592668", "Combo", "Combo 50% senescence"), drug = F )
# ) %>%
#   filter(!(drug %in% c("A-1592668", "Venetoclax") & time >15)) %>%
#
#
#   group_by(time, group, drug) %>%
#   summarise(q5 = quantile(TV, 0.05, na.rm = T), q95 = quantile(TV, 0.95, na.rm = T) ) %>%
#   mutate(q5 = q5 * 2.57, q95 = q95 * 2.57) %>%
#   mutate(drug=if_else(drug == T, "Treatment", "Control")) %>%
#
#   ggplot()+
#   geom_ribbon(aes(x = time, ymin = q5, ymax = q95, fill = drug), alpha = 0.5) +
#   # geom_line(aes(time, TV* 2.57, col = Dose))+
#   facet_wrap(~fct_relevel(group, "Venetoclax"), scales = "free")+
#   # scale_y_log10()+
#   geom_line(data = OBS,aes(TIME, DV))+
#   geom_line(data = compar %>% filter(ID %in% c(1)) ,aes(TIME, DV))+
#   theme_bw()+
#   labs(x = "Time (days)", y = "Tumor Volume (mm3)", fill = "Group", lty = "Observations") +
#   # scale_linetype_manual(values = c(1,2,3,4)[obs])+
#   geom_vline(xintercept = 1, lty = 3)

bind_rows(

  reswithVeneto2d %>% filter(time <21) %>% mutate(group = "Venetoclax", drug = T),
  reswithA152d %>% filter(time <21) %>% mutate(group = "A-1592668", drug = T),
  reswithCombo2d %>% filter(time <29) %>% mutate(group = "Combo", drug = T),
  reswithCombo2d_0.4  %>% filter(time <29)  %>% mutate(group = "Combo 50% senescence", drug = T),
  crossing(rescontrol %>% filter(time <15), group = c("Venetoclax","A-1592668", "Combo", "Combo 50% senescence"), drug = F )
) %>%
  filter(!(drug %in% c("A-1592668", "Venetoclax") & time >15)) %>%
  write.table(file = "ABT2D.csv")


read.table("ABT2D.csv", header = T) %>%
  group_by(time, group, drug) %>%
  summarise(q5 = quantile(TV, 0, na.rm = T), q95 = quantile(TV, 1, na.rm = T) ) %>%
  mutate(q5 = q5 * 2.57, q95 = q95 * 2.57) %>%
  mutate(drug=if_else(drug == T, "Treatment", "Control")) %>%

  ggplot()+
  geom_ribbon(aes(x = time, ymin = q5, ymax = q95, fill = drug), alpha = 0.5) +
  # geom_line(aes(time, TV* 2.57, col = Dose))+
  facet_wrap(~fct_relevel(group, "Venetoclax"), scales = "free", nrow = 1)+
  # scale_y_log10()+
  geom_line(data = OBS,aes(TIME, DV, lty ="Observations"))+
  geom_line(data = compar %>% filter(ID %in% c(1)) ,aes(TIME, DV))+
  theme_bw()+
  labs(x = "Time (days)", y = "Tumor Volume (mm3)", fill = "ABM simulations", lty = "") +
  # scale_linetype_manual(values = c(1,2,3,4)[obs])+
  geom_vline(xintercept = 1, lty = 3)->  whole2dplot


# linear scale
plot_grid(whole2dplot+ggtitle("Virtual Tumor calibrated with two drugs")+
            theme(plot.title = element_text(hjust = 0.5)), whole3dplot+ggtitle("Virtual Tumor calibrated with three drugs")+
            theme(plot.title = element_text(hjust = 0.5)), ncol = 1, labels = c("A", "B"))

# log scale
plot_grid(whole2dplot+ggtitle("Virtual Tumor calibrated with two drugs")+
            theme(plot.title = element_text(hjust = 0.5))+   scale_y_log10(), whole3dplot+ggtitle("Virtual Tumor calibrated with three drugs")+
            theme(plot.title = element_text(hjust = 0.5))+   scale_y_log10(), ncol = 1, labels = c("A", "B"))


# both
plot_grid(whole2dplot   +ggtitle("Virtual Tumor calibrated with two drugs") +  theme(plot.title = element_text(hjust = 0.5)),
          whole2dplot+   scale_y_log10() , whole3dplot+ggtitle("Virtual Tumor calibrated with three drugs")+
            theme(plot.title = element_text(hjust = 0.5)),
          whole3dplot + scale_y_log10(),ncol = 1, labels = c("A", "", "B", ""))
# reswithA15_NR <- ABM_non_prolif(simulations_A15, nsimul,pctNP = 0.4 )
# Simulations venetoclax  -------------------------------------------------


## protocol

## Compute time death

# simulation <- simulations_veneto
## Perform ABM





reswithVeneto_NP0.5 <- ABM_non_prolif(simulations_veneto, nsimul,pctNP = 0.4)
plot_ABM(reswithVeneto_NP0.5, rescontrol, name = c("Veneto", "Control"), obs = 1:2)

# Find KPD A-1592 -------------------------------------------------

# Add events


# Simulations
simulations_A15 <- simulations(VT, add_events_A15, modifParam = c( ke_Mcl1_I = 0,   hillA15 = 5,
                                                                   EC50A15 = 1E-6))

reswithA15 <- ABM_simple(simulations_A15, nsimul )

reswithA15_NR <- ABM_non_prolif(simulations_A15, nsimul,pctNP = 0.4 )

plot_ABM(rescontrol,reswithA15, name = c(" Control", "A-15"), obs = c(1,3))
plot_ABM(rescontrol,reswithA15_NR, name = c(" Control", "A-15"), obs = c(1,3))


# Combo -------------------------------------------------------------------


# Add events



# add_events_Veneto <- tibble(Proto = c("1"), cmt = c("Veneto_gut"), time = c(0), amt = c(50), method = c("add"), ADM = c("1"), evid = c(1))
#
# add_events_Veneto <- add_events_Veneto %>%
#   slice(1) %>%
#   select(-time) %>%
#   crossing(time = seq(0,21*24,24))
#
# add_events_combo <-  bind_rows( add_events_mcl1, add_events_Veneto) %>% arrange(time) %>%
#   mutate(time = time + 700)

# Simulations
simulations_combo <- simulations(VT, add_events_combo, modifParam = c(  ratioTumor = 0.3))

reswithCombo <- ABM_simple(simulations_combo, nsimul )

plot_ABM(rescontrol,reswithCombo, name = c(" Control", "Combo"), obs = c(1,4))


reswithCombo_NP0.5 <- ABM_non_prolif(simulations_combo, nsimul,pctNP = 0.4, nres=3)
plot_ABM(reswithCombo_NP0.5, rescontrol, name = c("Combo", "Control"), obs  = c(1,4))


reswithCombo_NP <- ABM_simple(simulations_combo, nsimul )

# allplot reunited --------------------------------------------------------

plot_ABM(rescontrol,reswithA15, reswithVeneto_NP0.5, reswithCombo_NP0.5, name = c(" Control", "A-15", "Veneto", "Combo"))



# delay time death impact -------------------------------------------------
# compar <- read.table("D:/these/Second_project/QSP/modeling_work/In_vivo/SU_DHL4_full.csv", header = T, sep = ";", na.strings = ".")
compar <- read.table(here('data', 'mice_SU_DHL4_full.csv'), header = T, sep = ";", na.strings = ".")

death_impact <- tibble(delay = c(0,5:10)) %>%
  mutate(simulation = map(delay,function(x){

    print(x)
    simulations_combo %>% mutate(timedeath = timedeath +  x) %>%
      ABM_simple(nsimul = 1, timeend = 30)

  } )) %>%
  unnest()


death_impact %>%
  filter(TV * 2.57 < 2000) %>%
  ggplot()+
  geom_line(aes(time, TV  * 2.57, col = factor(delay)))+
  scale_y_log10()+
  geom_line(data = compar %>% filter(ID %in% c(1,4)), aes(TIME, DV, group = ID))+
  geom_vline(xintercept = 1, lty = 3)

# At the end it's parallel coz log(2)/0.131 = 5.29, max(simul) =  2.125

# same with delay between 0:x

death_impact <- tibble(delay = c(0,1,2,3,4,5,10)) %>%
  mutate(simulation = map(delay,function(x){

    print(x)
    simulations_combo %>% mutate(timedeath = timedeath + sample(0:x, replace = T, size = 100)) %>%
      ABM_simple(nsimul = 1, timeend = 30)

  } )) %>%
  unnest()



death_impact %>%
  filter(TV * 2.57 < 2000) %>%
  ggplot()+
  geom_line(aes(time, TV  * 2.57, col = factor(delay)))+
  scale_y_log10()+
  geom_line(data = compar %>% filter(ID %in% c(1,4)), aes(TIME, DV, group = ID))+
  geom_vline(xintercept = 1, lty = 3)


# With non multiplicative cells -------------------------------------------

crossing(pctNP = seq(0,0.5,0.1), nres = 0:10) %>%
  mutate(data = map2(pctNP, nres, function(pctNP, nres){

    set.seed(135643)


    simulation <-   simulations_combo

    t0 <- simulation %>% select(cellid, timedeath) %>%
      mutate(timedeath0 = timedeath) %>%
      mutate(timedeath = timedeath  + 1) %>%   # +1 because adm after day 1
      mutate(NP = sample(c(T,F),prob = c(pctNP, 1 - pctNP),
                         nrow(simulations_combo), replace = T))

    # nres <- 2
    t0$timedeath[1:nres] <- Inf
    t0$timedeath0[1:nres] <- Inf
    t0$NP[1:nres] <- F

    tim <- 0
    step <- 0.1
    lambda0 <- 0.131
    probmult <- lambda0 * step

    res <- tibble(time = 0, TV = 100)

    while(tim < 30){
      tim  <- tim + step


      t0 %>%
        mutate(testMult =  sample(x = c(T, F), prob = c(probmult, 1 - probmult), size = nrow(t0), replace = T)) %>%
        mutate(testDeath =  timedeath < tim) %>%
        filter(testDeath == F | NP == T) -> temp

      death <- nrow(t0) - nrow(temp)
      # if(death > 0) print(paste0(death," death"))
      # print
      # -> temp

      t0 <- bind_rows(temp, temp %>% filter(testMult == T & NP == F) %>% mutate(timedeath =  timedeath0 + tim)) # - 1 because we added 1 before
      res <- res %>% bind_rows(tibble(time = tim, TV = nrow(t0)))
    }

    res
  })) -> temp


temp %>%
  filter(nres !=0) %>%
  unnest() %>%
  ggplot()+
  geom_line(aes(time, TV  * 2.57), col = "red")+
  scale_y_log10()+
  geom_line(data = compar %>% filter(ID %in% c(1,4)), aes(TIME, DV, group = ID))+
  geom_vline(xintercept = 1, lty = 3)+
  facet_grid(pctNP ~ nres)



pctNP <- 0.3

simulation <-   simulations_combo

t0 <- simulation %>% select(cellid, timedeath) %>%
  mutate(timedeath0 = timedeath) %>%
  mutate(timedeath = timedeath  + 1) %>%   # +1 because adm after day 1
  mutate(NP = sample(c(T,F),prob = c(pctNP, 1 - pctNP),
                     nrow(simulations_combo), replace = T))

nres <- 2
t0$timedeath[1:nres] <- Inf
t0$timedeath0[1:nres] <- Inf
t0$NP[1:nres] <- F

tim <- 0
step <- 0.1
lambda0 <- 0.131
probmult <- lambda0 * step

res <- tibble(time = 0, TV = 100)

while(tim < 30){
  tim  <- tim + step


  t0 %>%
    mutate(testMult =  sample(x = c(T, F), prob = c(probmult, 1 - probmult), size = nrow(t0), replace = T)) %>%
    mutate(testDeath =  timedeath < tim) %>%
    filter(testDeath == F | NP == T) -> temp

  death <- nrow(t0) - nrow(temp)
  # if(death > 0) print(paste0(death," death"))
  # print
  # -> temp

  t0 <- bind_rows(temp, temp %>% filter(testMult == T & NP == F) %>% mutate(timedeath =  timedeath0 + tim)) # - 1 because we added 1 before
  res <- res %>% bind_rows(tibble(time = tim, TV = nrow(t0)))
}



res %>%
  filter(TV * 2.57 < 2000) %>%
  ggplot()+
  geom_line(aes(time, TV  * 2.57), col = "red")+
  scale_y_log10()+
  geom_line(data = compar %>% filter(ID %in% c(1,4)), aes(TIME, DV, group = ID))+
  geom_vline(xintercept = 1, lty = 3)



# Same for venetoclax -----------------------------------------------------


simulation <-   simulations_veneto
pctNP <- 0.3
t0 <- simulation %>% select(cellid, timedeath) %>%
  mutate(timedeath0 = timedeath) %>%
  mutate(timedeath = timedeath  + 1) %>%   # +1 because adm after day 1
  mutate(NP = sample(c(T,F),prob = c(pctNP, 1 - pctNP),
                     nrow(simulations_combo), replace = T))

# nres <- 2
t0$timedeath[1:nres] <- Inf
t0$timedeath0[1:nres] <- Inf
t0$NP[1:nres] <- F

tim <- 0
step <- 0.1
lambda0 <- 0.131
probmult <- lambda0 * step

res <- tibble(time = 0, TV = 100)

while(tim < 30){
  tim  <- tim + step


  t0 %>%
    mutate(testMult =  sample(x = c(T, F), prob = c(probmult, 1 - probmult), size = nrow(t0), replace = T)) %>%
    mutate(testDeath =  timedeath < tim) %>%
    filter(testDeath == F | NP == T) -> temp

  death <- nrow(t0) - nrow(temp)
  # if(death > 0) print(paste0(death," death"))
  # print
  # -> temp

  t0 <- bind_rows(temp, temp %>% filter(testMult == T & NP == F) %>% mutate(timedeath =  timedeath0 + tim)) # - 1 because we added 1 before
  res <- res %>% bind_rows(tibble(time = tim, TV = nrow(t0)))
}

# compar <- read.table("D:/these/Second_project/QSP/modeling_work/In_vivo/SU_DHL4_full.csv", header = T, sep = ";", na.strings = ".")
compar <- read.table(here('data', 'mice_SU_DHL4_full.csv'), header = T, sep = ";", na.strings = ".")
OBS <-
  compar %>% filter(ID %in% c(1:4)) %>%
  mutate(OBS = case_when(ID == 1 ~ "Control",
                         ID == 2 ~ "Venetoclax",
                         ID == 3 ~ "A-15",
                         T ~ "Combo")) %>%
  filter(ID %in% 1:2)


res %>%
  ggplot()+
  geom_line(aes(x = time, y = TV * 2.57), col = "red") +
  # geom_line(aes(time, TV* 2.57, col = Dose))+
  scale_y_log10()+
  geom_line(data = OBS,aes(TIME, DV, group = factor(Dose), lty = fct_reorder(OBS, ID)))+
  theme_bw()+
  labs(x = "Time (days)", y = "Tumor Volume", fill = "Simulations", lty = "Observations") +
  scale_linetype_manual(values = c(1,2))+
  geom_vline(xintercept = 1, lty = 3)
