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
# Source the file 2_Virtual_Tumor_objects.R
source(here("2_Virtual_Tumor_objects.R"))

# Same process for each VT (for each cell line, optimixing vene, A11 or both togheter)
# Venetoclax - SUDHL4 -----------------------------------------------

# Initiate the VT. 
#Drug 1= Venetoclax, Drug 2 = A-1155463
# Cell_line_bin 1 == SU-DHL-4 , cell_line_bin 2 == KARPAS-422 (look at data_VT object if doubt)
self <- VT2$new(filter_data = Drug %in% 1 & Cell_line_bin == 1 , pen_nbag = 10,  name = "test")

# Just here as a demo, the VT object self-contains the celltheque (with number of bags depending on how
# the drugs remaining in the dataset to calibrate) and the current 100 VCs composing the VT
self$harvest # The cellid of the 100 current VCs
self$celltheque # The celltheque corresponds to one cell per possible bag (with cellid to match current harvest)

# The function reconstitue provide the  results with current 100 VCs (self$harvest)
self$reconstitute()

# You can plot the current fitting and get current OF
self$plot()
self$OF()
# Now you can use first method to optimize. Note it is a endless loop so you'll have to choose
# when to stop it (eg, when no longer improvement of OF printed in console)
self$optim()

#Then you can reupdate the previous plot and OF to see the improvment
self$plot()
self$OF()

# Now you can use the second method to optimize. Take times, so provide a path where to save the 
# the current Virtual tumor object, and frequency of saving
self$cell_by_cell_last <- NULL # useless the first time, allow to reset the algorithm if needed
self$optim2(saveEveryXmin = 10,pathSave = here("calibrated_VT/VT_Venetoclax_cell_line_1.RDS" ))
# Finally, save the final Virtual tumor
saveRDS(self, here("calibrated_VT/VT_Venetoclax_cell_line_1.RDS"))

# The same process was applied for all VT below
# Venetoclax - KARPAS422 -----------------------------------------------

self <- VT2$new(filter_data = Drug %in% 1 & Cell_line_bin == 2 , pen_nbag = 10,  name = "test")
self$plot()
self$OF()
self$optim()
self$cell_by_cell_last <- NULL
self$optim2(saveEveryXmin = 10,pathSave = "calibrated_VT/VT_Veneto_cell_line_2.RDS" )


saveRDS(self, here( "calibrated_VT/VT_Venetoclax_cell_line_1.RDS"))

# A11 - SUDH4 -----------------------------------------------
self <- VT2$new(filter_data = Drug %in% 2 & Cell_line_bin == 1 , pen_nbag = 10,  name = "test")
self$plot()
self$OF()
self$optim()
self$cell_by_cell_last <- NULL
self$optim2(saveEveryXmin = 10,pathSave = "calibrated_VT/VT_A11_cell_line_1.RDS" )
saveRDS(self, here(self, "calibrated_VT/VT_A11_cell_line_1.RDS"))

# A11 - KARPAS -----------------------------------------------
self <- VT2$new(filter_data = Drug %in% 2 & Cell_line_bin == 2 , pen_nbag = 10,  name = "test")
self$plot()
self$OF()
self$optim()
self$cell_by_cell_last <- NULL
self$optim2(saveEveryXmin = 10,pathSave = "calibrated_VT/VT_A11_cell_line_2.RDS" )
saveRDS(self, here(self, "calibrated_VT/VT_A11_cell_line_2.RDS"))
# 3 drugs - SUDHL4 --------------------------------------------------------



self <- VT2$new(filter_data = Drug %in% c(1,2) & Cell_line_bin == 1 , pen_nbag = 10,  name = "test")

# prev_harvest <- readRDS("calibrated_VT/cell1prevharvest.RDS")
# self$harvest <- prev_harvest
# self <- readRDS("calibrated_VT/VT_both_cell_line_1.RDS" )

self$plot()
self$OF(detail = T)
self$optim()
self$optim2(saveEveryXmin = 10,pathSave = "calibrated_VT/VT_both_cell_line_1.RDS" )
saveRDS(self, here( "calibrated_VT/VT_both_cell_line_1.RDS"))


# 3 drugs KARPAS ----------------------------------------------------------

self <- VT2$new(filter_data = Drug %in% c(1,2) & Cell_line_bin == 2, pen_nbag = 10,  name = "test")
self$plot()
self$OF(detail = T)
self$optim()
self$optim2(saveEveryXmin = 10, pathSave = "calibrated_VT/VT_both_cell_line_2.RDS")
saveRDS(self, here( "calibrated_VT/VT_both_cell_line_2.RDS"))



# Comparison drug --------------------------------------------

VT <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/3_virtual_tumors/SU-DHL-4-Venetoclax.RDS")

VT2 <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/3_virtual_tumors/KARPAS-Venetoclax.RDS")


VT_description(VT) %>%
  rename(SUDH =n) %>%
  full_join(VT_description(VT2) %>%
              rename(Karpas = n))

VT_prot_expr(VT) %>% mutate(name = VT@name) %>%
  bind_rows(

    VT_prot_expr(VT2) %>% mutate(name = VT2@name)
  )


VT_prot_expr(VT,return_median_i = T) %>% mutate(name = VT@name) %>%
  bind_rows(

    VT_prot_expr(VT2,return_median_i = T) %>% mutate(name = VT2@name)
  ) -> temp

temp %>%
  mutate(name = gsub("-Ve.*", "", name)) %>%
  select(-BAK0, -BAXc0,- cellid_in_harvest, - Drug) %>%
  gather("key", "value", -bag, - name) %>%
  ggplot()+
  geom_density(aes(value, fill = name), alpha = 0.5)+
  facet_wrap(~key, scales = "free")




# Recreating bags ------------------------------------------------------------
# In VT celltheque, there is only one cell per bag... but in reality we 
# have computed several cells that can belong in each bag..
# so lets retake all the possible cells needed:
# The file generated are further used in the '4_agent_based_model.R' file

sudh2d <-  readRDS(here('calibrated_VT','VT_Venetoclax_cell_line_1.RDS'))  #readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/new_VT_qspvp/")
allbags <- function(self){
  
  allprofiles <- tibble(cellid = self$harvest) %>% 
    left_join(self$celltheque) %>% 
    distinct() %>% 
    select(-Bcl20, - Bclxl0, - Mcl10  , -BIM0, - PUMA0, - NOXA0 , -  BAK0, - BAXc0 )
  
  
  
  path <- here('Virtual_Cells')#"D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/new_celltheque_QSPVP"
  setwd(path)
  
  allfiles <- list.files(); allfiles <- allfiles[grepl(".RDS",allfiles)]
  allfiles <- allfiles[!grepl('spont_death',allfiles)]
  a <- allfiles[[1]]
  
  allcells <- tibble()
  
  for(a in allfiles){
    
    print(paste0(which(allfiles == a), "/", length(allfiles)))
    
    tempDF <- readRDS(a)
    
    tempDF <- tempDF  %>% 
      left_join(allprofiles) %>% 
      filter(!is.na(cellid))
    
    allcells <- bind_rows(allcells, tempDF)
    print(nrow(allcells))
  }
  
  allcells
}

karpas_bags <- allbags(karpas)
saveRDS(karpas_bags, here('calibrated_VT', 'karpas_all_bags.RDS'))
sudh_bags <- allbags(sudh)
saveRDS(sudh_bags,  here('calibrated_VT', 'sudh_all_bags.RDS'))

sudh2d_bags <- allbags(sudh2d)
saveRDS(sudh2d_bags, here('calibrated_VT', 'sudh2d_all_bags.RDS'))

