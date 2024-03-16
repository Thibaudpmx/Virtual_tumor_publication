# Author: Thibaud Derippe 
# The goal of this file is to merge and reduce the results
# And creates the "Bags" for each drug combination


# path_VC = path of the folder used in the file 1_1_Virtual_cells_Generation
path_VC <- here("Virtual_Cells")
setwd(path_VC)



# Merging A11 and Veneto File -------------------------------------------------------------


# Goal: merge the A11 and Veneto file for getting full res profile for each VC
tomerge <- list.files(path_VC)[grepl("_A11.RDS",list.files(path_VC))] # All File of A-11


a <- tomerge[[1]]


for(a in tomerge){ # For each _A11 file

  venetofile <- gsub("A11", "veneto", a) # Get the equivalent Veneto file


  A11 <- readRDS(file.path(path_VC,a)) %>% select(-Bcl2_I0) # Read A11 res

  veneto <- readRDS(file.path(path_VC,venetofile)) %>% # Read Veneto Res
    select(-Bclxl_I0)


  A11 %>%
    left_join(veneto) %>% # Join by all parameters values to make sure it is the same VC
    saveRDS(file.path(path_VC, gsub("_A11", "_both", a))) # save the file finishing by "Both"

  # So now both as the followng format
  # A11_0  A11_5 A11_10 A11_15 Bcl20 Bclxl0 Mcl10  BIM0 PUMA0 NOXA0  BAK0 BAXc0 Veneto_0 Veneto_5 Veneto_10 Veneto_15

  unlink(file.path(path_VC,a)) # Delete the file of A11 and Veneto,  because "_both"
  unlink(file.path(path_VC,venetofile)) # already contains all the information !
}


# Merging all files from same group ---------------------------------------


# Goal: merge all results from the same same groups

tomerge <- list.files(path_VC)[grepl("_group",list.files(path_VC))]


# a <- tomerge[[1]]

while(sum(grepl("_group",list.files(path_VC))) >0){ # Until there is _group remaining

  a <- list.files(path_VC)[grepl("_group",list.files(path_VC))][[1]] # Take the first "group" remaining

  key <- gsub("_group.+", "", a) # get the key to gather (ex: key = "Bax0_BaK1000")


  allfiles <- list.files(path_VC)[grepl(key,list.files(path_VC))] # Get all file with this key
  allfiles <- allfiles[grepl("_both.RDS", allfiles)] #make sure to keep only the "_both" files


 merged <-  map(allfiles, ~readRDS( file.path(path_VC,.x))) %>% bind_rows() # and combine them!


 merged %>%
   saveRDS(file.path(path_VC, paste0(key, "_both.RDS"))) # save the merged dataset


  map(allfiles, ~ unlink(file.path(path_VC,.x))) # and remove individual group files
  # (untile all files with _group are removed)
}



# Minimal bag per file ------------------------------------------------------------



allfiles <- list.files(path_VC)[grepl("_both.RDS",list.files(path_VC))]

a <- allfiles[[1]]

for(a in allfiles){ # For each "_both" file

df <- readRDS(a) # read the file

# reduce per A11
if(!dir.exists('A11_one_per_bag')) dir.create('A11_one_per_bag')
if(!dir.exists('venetoclax_one_per_bag')) dir.create('venetoclax_one_per_bag')
if(!dir.exists('both_one_per_bag')) dir.create('both_one_per_bag')
df %>%
  group_by(A11_0, A11_5, A11_10, A11_15) %>%
  slice(1) %>%
  select(-starts_with("Veneto")) %>%
  saveRDS(file.path("A11_one_per_bag", a))

# reduce per Veneto

df %>%
  group_by(Veneto_0, Veneto_5, Veneto_10, Veneto_15) %>%
  slice(1) %>%
  select(-starts_with("A11")) %>%
  saveRDS(file.path("venetoclax_one_per_bag", a))

# reduce per both

df %>%
  group_by(Veneto_0, Veneto_5, Veneto_10, Veneto_15, A11_0, A11_5, A11_10, A11_15) %>%
  slice(1) %>%
  saveRDS(file.path("both_one_per_bag", a))
}





# Now creates the final bags: one cell for each concentrations sen --------


dir.create("one_per_bag_combined")

# spont death (just add one line)

spontdeath <- readRDS("Bax500_Bak0_spont_death.RDS") %>%
  slice(1) %>%
  mutate(A11_0 = 0, A11_5 = 0, A11_10 = 0, A11_15 = 0,
         Veneto_0 = 0, Veneto_5 = 0, Veneto_10 = 0, Veneto_15 = 0);spontdeath


# A11: read all files, merged them, reduce


A11all <- list.files("A11_one_per_bag") #


map(A11all, ~readRDS(file.path("A11_one_per_bag", .x))) %>%
  bind_rows() %>%
  group_by(A11_0, A11_5, A11_10, A11_15) %>%
  slice(1) %>%
  bind_rows(spontdeath %>% select(-Veneto_0, -Veneto_10, - Veneto_15, - Veneto_5)) %>%
  saveRDS(file.path("one_per_bag_combined", "A11.RDS"))


# Veneto: read all files, merged them, reduce

Venetoall <- list.files("venetoclax_one_per_bag")

map(Venetoall, ~readRDS(file.path("venetoclax_one_per_bag", .x))) %>%
  bind_rows() %>%
  group_by(Veneto_0, Veneto_5, Veneto_10, Veneto_15) %>%
  slice(1) %>%
  bind_rows(spontdeath %>% select(-A11_0, -A11_5, - A11_10, - A11_15)) %>%
  saveRDS(file.path("one_per_bag_combined", "Veneto.RDS"))

# Combination: read all files, merged them, reduce

bothAll <- list.files("both_one_per_bag")

map(bothAll, ~readRDS(file.path("both_one_per_bag", .x))) %>%
  bind_rows() %>%
  group_by(Veneto_0, Veneto_5, Veneto_10, Veneto_15, A11_0, A11_5, A11_10, A11_15) %>%
  slice(1) %>%
  bind_rows(spontdeath) %>%
  # is.na() %>% sum
  saveRDS(file.path("one_per_bag_combined", "both.RDS"))

# And voilà !
