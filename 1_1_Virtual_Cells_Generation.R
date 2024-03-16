library(tidyr)
library(dplyr)
library(purrr)
library(rlang)
library(PaSM)
library(here)
# Author: Thibaud Derippe
# Objectif: generating virtual Cells using PaSM algorithm.
# Two methods: 1) equidistant parameter sampling (meshing) and 2) random values
# After this script, use the file 1_2_Rearrange_VC_files to have a better organisations of the VCs

# Method 1 Equidistant parameter sampling (meshing) -----------------------------------------------

# 9 values per monotonic parameter (10 for Mcl1), 9^5 * 10 = 590,490 possible combination
ToaddBrut <- crossing( Bcl20 = seq(20,1020,112), Bclxl0 = seq(50,1500, 170),
                       Mcl10 = seq(2,150, 16), BIM0 = seq(0,200,25), PUMA0 = seq(0,200,25), NOXA0 =  seq(0,200,25))


# 8 combination values per non-monotonic parameter (Bax and Bak), 8 * 590,490 = 4,723,920 total cells!
BakBaxcombo <- crossing(BAK0 = c(0,500,1000), BAXc0 = c(0,500,1000)) %>%
  filter(!(BAK0 == 0 & BAXc0 == 0));BakBaxcombo


# Path of the PaSM config file (includes model ODE system and monotonic parameter sets definition)
sourcefile <- here("0_Lindner_model_PaSM_config.R")

# Folder containing all the generated virtual cells QSP results
# path_output <- "D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/celltheque_verif_article" # C:/Users/thibaudd/Desktop/new_celltheque_QSPVP
path_output <- here("Virtual_Cells")
dir.create(path_output)

## Lets go !

for(a in 1:nrow(BakBaxcombo)){ # For each of the 8 possible Bax/Bak combination


# Define the targets corresponding to "apoptosis triggering"
# Note the only protocol we will use is "Dose0", because dose are
# administered as initial value of compartment (monotonic parameters)
  targets <- tribble(~protocol, ~cmt, ~time, ~min, ~max,
                     "dose0", "TimeAbove",  748 , 1E-4, Inf)

# Add the two non-monontonic parameter for the loop turn to the df of unique VC
  WholeDF <-   ToaddBrut %>% crossing(BakBaxcombo %>% slice(a))

# Extract the BaxC0 and Bak value to create the name of the file to save
  BAXc0_a <- unique(WholeDF$BAXc0)
  BAK0 <- unique(WholeDF$BAK0)
  nameBaxBAK <- paste0("Bax",BAXc0_a,
                       "_BaK", BAK0)

# Step 1: extract and remove the VP with spontaneous death without any drug

  # name of the file storing all spontaneous apoptosis VC
  file_alldead <- file.path(path_output, paste0(nameBaxBAK, "_spont_death.RDS"))

  if(!file.exists(file_alldead)){


    # Create a new PaSM object
    searchSpont <- PaSM::VP_proj_creator$new(sourcefile)

    # Add the targets
    searchSpont$set_targets(manual = targets )

    # Use PaSM algorithm  without any drug
    searchSpont$add_VP(WholeDF %>% mutate(Bcl2_I0 = 0 ,Bclxl_I0 = 0,Mcl1_I0 = 0), npersalve = 1000, use_green_filter = T)


    # Extract only the cells without spontaneous death ("Accepted" VPs) and save them.
    spontDeath <- searchSpont$poolVP %>% select(Bcl20, Bclxl0, Mcl10,  BIM0, PUMA0, NOXA0,  BAK0, BAXc0)


    spontDeath %>%
      saveRDS(file_alldead)




  }else{

    spontDeath <- readRDS(file_alldead)

  }


 # Remove the previous VCs that spontaneously trigger apoptosis
  WholeDF %>%
    left_join(spontDeath %>% mutate(toRem = T)) %>%
    filter(is.na(toRem)) %>%
    select(-toRem) %>%
    rowid_to_column() %>%
    mutate(group = floor(rowid/4E4) + 1 ) %>% ####  And create subgroups for the analysis
    select(-rowid) -> cancerDF


  for(b in unique(cancerDF$group)){ # For each of this subgroup

    groupdf <-  cancerDF %>% filter(group == b) # Select only the group to analyse

    concx <-  c(0,0.08,0.16,0.32,0.64,1.3,2.6,5,10,20) # Add concentrations of main drugs (Veneto or A-11)
    concy <-  c(0 , 5, 10, 15) # Same for A-12




    # Determine fates of  A11 +/- A12 drugs concentrations


    file_groupb_A11 <- file.path(path_output, paste0(nameBaxBAK, "_group", b, "_A11.RDS")) # name of the subgroup files


    if(!file.exists(file_groupb_A11)){


      VP_df <- crossing(groupdf, Bclxl_I0 = concx, Mcl1_I0 = concy, Bcl2_I0 = 0 ) # Add the drugs concentrations to try (as temporar "distinct VPs")
      self2 <- QSPVP::VP_proj_creator$new(sourcefile) # Create a new PaSM object


      self2$set_targets(manual = targets ) # Add targets

      # Launch PaSM main algorithm
      self2$add_VP(VP_df, npersalve = 1000, use_green_filter = T, keep = c("Pore", "TimeAbove"),
                   timeSave = c(700,120,748))


      # Take back all initial df
      VP_df %>%
        left_join(self2$poolVP %>% select(!!!self2$param,rowid )) %>% # left_join PaSM results
        mutate(status = if_else(is.na(rowid), F, T)) %>% # If VP is poolVP, then apotosis, otherwise (isnarowid) survival
        select(-rowid) -> restemp



      # Simplify the notation by only taking the first A-11 drug concentration of apoptosis for each A-12 values
      restemp %>%
        bind_rows(
          crossing(groupdf,Bclxl_I0 = Inf,  Mcl1_I0 = concy, status =  T, Bcl2_I0 = 0)
        ) %>%
        arrange(Mcl1_I0, Bclxl_I0) %>%
        filter(status) %>%
        group_by(!!!parse_exprs(names(WholeDF)), Mcl1_I0) %>%
        slice(1) %>%
        select(-status) %>%
        mutate(Mcl1_I0 = paste0("A11_", Mcl1_I0)) %>%
        spread(key = Mcl1_I0, value = Bclxl_I0) %>%
        select(A11_0, A11_5, A11_10, A11_15, everything() )  -> restemp2


       # Save the file
      restemp2 %>% select(-group) %>%
        saveRDS(file_groupb_A11)

    } # end  group A11



    # Compute Venetoclax + A-12 (same code as for A-11 + A12)...

    file_groupb_Veneto <- file.path(path_output, paste0(nameBaxBAK, "_group", b, "_veneto.RDS"))


    if(!file.exists(file_groupb_Veneto)){

      # ... except now  Bclxl_I0 is set to 0 and  Bcl2_I0 = concx
      VP_df <- crossing(groupdf, Bclxl_I0 = 0, Mcl1_I0 = concy, Bcl2_I0 = concx )
      self2 <- QSPVP::VP_proj_creator$new(sourcefile)


      self2$set_targets(manual = targets )

      self2$add_VP(VP_df, npersalve = 1000, use_green_filter = T, keep = c("Pore", "TimeAbove"),
                   timeSave = c(700,120,748))


      VP_df %>%
        left_join(self2$poolVP %>% select(!!!self2$param,rowid )) %>%
        mutate(status = if_else(is.na(rowid), F, T)) %>%
        select(-rowid) -> restemp



      # Simplify the notation by only taking the first Venetoclax drug concentration of apoptosis for each A-12 values

      restemp %>%
        bind_rows(
          crossing(groupdf,Bclxl_I0 = 0,  Mcl1_I0 = concy, status =  T, Bcl2_I0 = Inf)
        ) %>%
        arrange(Mcl1_I0, Bcl2_I0) %>%
        filter(status) %>%
        group_by(!!!parse_exprs(names(WholeDF)), Mcl1_I0) %>%
        slice(1) %>%
        select(-status) %>%
        mutate(Mcl1_I0 = paste0("Veneto_", Mcl1_I0)) %>%
        spread(key = Mcl1_I0, value = Bcl2_I0) %>%
        select(Veneto_0, Veneto_5, Veneto_10, Veneto_15, everything() )  -> restemp2 #

      restemp2 %>%
        ungroup() %>%
        group_by( Veneto_0, Veneto_5, Veneto_10, Veneto_15 ) %>% tally()


      restemp2 %>% select(-group) %>%
        saveRDS(file_groupb_Veneto)

    } # end  group Veneto



  } # end for each group

} # end for each Bax/Bak combo




# Method 2 Random value -----------------------------------------------





while(T){ # Unlimited loop, will generate cells until we stop the algorithm

  name <- Sys.time()  %>% gsub(pattern = " |-|:", replacement = "_")

  size_per_it <-  5E4


  # Sampling in the uniform low
  # Note: all the rest is similater to the code above

  WholeDF <-  tibble(

    Bcl20 = runif(size_per_it, 0, 2000),
    Bclxl0 = runif(size_per_it, 0, 2000),
    Mcl10 = runif(size_per_it, 0, 500),
    BIM0 = runif(size_per_it, 0, 1000),
    PUMA0 = runif(size_per_it, 0, 1000),
    NOXA0 =  runif(size_per_it, 0, 1000),
    BAK0 =   runif(1, 0, 2000),
    BAXc0 =   runif(1, 0, 2000)


  ) %>%
    map_df(~ round(.x, 2))




  targets <- tribble(~protocol, ~cmt, ~time, ~min, ~max,
                     "dose0", "TimeAbove",  748 , 1E-4, Inf)

  file_alldead <- file.path(path_output, paste0(name, "_spont_death.RDS"))

  searchSpont <- QSPVP::VP_proj_creator$new(sourcefile)




  searchSpont$set_targets(manual = targets )


  searchSpont$add_VP(WholeDF %>% mutate(Bcl2_I0 = 0 ,Bclxl_I0 = 0,Mcl1_I0 = 0), npersalve = 1000, use_green_filter = T)



  spontDeath <- searchSpont$poolVP %>% select(Bcl20, Bclxl0, Mcl10,  BIM0, PUMA0, NOXA0,  BAK0, BAXc0)


  spontDeath %>%
    saveRDS(file_alldead)




  WholeDF %>%
    left_join(spontDeath %>% mutate(toRem = T)) %>%
    filter(is.na(toRem)) %>%
    select(-toRem) %>%
    rowid_to_column() %>%
    mutate(group = floor(rowid/4E4) + 1 ) %>% #### here the number per group !!
    select(-rowid) -> cancerDF








  for(b in unique(cancerDF$group)){

    concx <-  c(0,0.08,0.16,0.32,0.64,1.3,2.6,5,10,20)
    concy <-  c(0 , 5, 10, 15)


    groupdf <-  cancerDF %>% filter(group == b)

    # Compute A11

    file_groupb_A11 <- file.path(path_output, paste0(name, "_group", b, "_A11.RDS"))


    if(!file.exists(file_groupb_A11)){


      VP_df <- crossing(groupdf, Bclxl_I0 = concx, Mcl1_I0 = concy, Bcl2_I0 = 0 )
      self2 <- QSPVP::VP_proj_creator$new(sourcefile)


      self2$set_targets(manual = targets )

      self2$add_VP(VP_df, npersalve = 1000, use_green_filter = T, keep = c("Pore", "TimeAbove"),
                   timeSave = c(700,120,748))


      VP_df %>%
        left_join(self2$poolVP %>% select(!!!self2$param,rowid )) %>%
        mutate(status = if_else(is.na(rowid), F, T)) %>%
        select(-rowid) -> restemp



      # try so simplify its notation
      restemp %>%
        # group_by(!!!parse_exprs(names(VP_df2))) %>% tally %>% filter(n != 40)
        bind_rows(
          crossing(groupdf,Bclxl_I0 = Inf,  Mcl1_I0 = concy, status =  T, Bcl2_I0 = 0)
        ) %>%
        # group_by(!!!parse_exprs(names(VP_df2))) %>% tally %>% filter(n != 44)
        arrange(Mcl1_I0, Bclxl_I0) %>%
        filter(status) %>%
        group_by(!!!parse_exprs(names(WholeDF)), Mcl1_I0) %>%
        # summarise(sum = sum(status)) %>% filter(sum == 1)
        slice(1) %>%
        select(-status) %>%
        mutate(Mcl1_I0 = paste0("A11_", Mcl1_I0)) %>%
        spread(key = Mcl1_I0, value = Bclxl_I0) %>%
        select(A11_0, A11_5, A11_10, A11_15, everything() )  -> restemp2 #

      # restemp2 %>%
      # ungroup() %>%
      # group_by(  A11_0  ,A11_5, A11_10, A11_15 ) %>% tally()


      restemp2 %>% select(-group) %>%
        saveRDS(file_groupb_A11)

    } # end  group A11



    # Compute Venetoclax

    file_groupb_Veneto <- file.path(path_output, paste0(name, "_group", b, "_veneto.RDS"))


    if(!file.exists(file_groupb_Veneto)){


      VP_df <- crossing(groupdf, Bclxl_I0 = 0, Mcl1_I0 = concy, Bcl2_I0 = concx )
      self2 <- QSPVP::VP_proj_creator$new(sourcefile)


      self2$set_targets(manual = targets )

      self2$add_VP(VP_df, npersalve = 1000, use_green_filter = T, keep = c("Pore", "TimeAbove"),
                   timeSave = c(700,120,748))


      VP_df %>%
        left_join(self2$poolVP %>% select(!!!self2$param,rowid )) %>%
        mutate(status = if_else(is.na(rowid), F, T)) %>%
        select(-rowid) -> restemp



      # try so simplify its notation
      restemp %>%
        # group_by(!!!parse_exprs(names(VP_df2))) %>% tally %>% filter(n != 40)
        bind_rows(
          crossing(groupdf,Bclxl_I0 = 0,  Mcl1_I0 = concy, status =  T, Bcl2_I0 = Inf)
        ) %>%
        # group_by(!!!parse_exprs(names(VP_df2))) %>% tally %>% filter(n != 44)
        arrange(Mcl1_I0, Bcl2_I0) %>%
        filter(status) %>%
        group_by(!!!parse_exprs(names(WholeDF)), Mcl1_I0) %>%
        # summarise(sum = sum(status)) %>% filter(sum == 1)
        slice(1) %>%
        select(-status) %>%
        mutate(Mcl1_I0 = paste0("Veneto_", Mcl1_I0)) %>%
        spread(key = Mcl1_I0, value = Bcl2_I0) %>%
        select(Veneto_0, Veneto_5, Veneto_10, Veneto_15, everything() )  -> restemp2 #

      restemp2 %>%
        ungroup() %>%
        group_by( Veneto_0, Veneto_5, Veneto_10, Veneto_15 ) %>% tally()


      restemp2 %>% select(-group) %>%
        saveRDS(file_groupb_Veneto)

    } # end  group Veneto



  } # end for each group


} # end_while_loop



