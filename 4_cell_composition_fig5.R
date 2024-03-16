library(progress)


# Author: Thibaud Derippe 
# Goal: analyse the cell composition/heterogeneity in the final calibrated VTs
# Final drugs -------------------------------------------------------------

karpas <- readRDS(here('calibrated_VT', 'VT_both_cell_line_2.RDS'))
sudh <- readRDS(here('calibrated_VT', 'VT_both_cell_line_1.RDS'))
self <- karpas
# 
# self$plot()


compositionplot <- function(self, labels = c("A", "B", "C")){
  fullcells <- tibble(cellid = self$harvest) %>%
    left_join(self$celltheque)
  
  fullcells2 <- fullcells %>%
    mutate(A11unk = if_else(A11_0 == Inf, T, F)) %>%
    mutate(Venetounk = if_else(Veneto_0 == Inf, T, F)) %>%
    mutate(A12unk = if_else(A11_15 > 0, T, F)) %>%
    mutate(A11_A12_unk = if_else(A11_15 == Inf, T, F))%>%
    mutate(Veneto_A12_unk = if_else(Veneto_15 == Inf, T, F)) %>%
    mutate(Veneto_only = if_else(A11unk == T & A12unk == T & Venetounk == F, T, F)) %>%
    mutate(A11_only = if_else(A11unk == F & A12unk == T & Venetounk == T, T, F)) %>%
    mutate(A12_only = if_else(A11unk == T & A12unk == F & Venetounk == T, T, F)) %>%
    mutate(Veneto_or_A11 = if_else(A11unk == F & Venetounk == F & A12unk == T , T, F)) %>%
    mutate(Veneto_or_A12 = if_else(A11unk == T & Venetounk == F & A12unk == F , T, F)) %>%
    mutate(A11_or_A12 = if_else(A11unk == F & Venetounk == T & A12unk == F , T, F)) %>%
    mutate(any = if_else(A11unk == F & Venetounk == F & A12unk == F , T, F)) %>%
    mutate(unkill_solo = if_else(A11unk == T & Venetounk == T & A12unk == T , T, F)) %>%
    mutate(A11_A12_needed = if_else(A11_A12_unk == F & A12unk == T  & A11unk == T & Venetounk == T, T, F)) %>% # Combo effective while not solo?
    mutate(Veneto_A12_needed = if_else(Veneto_A12_unk == F & A12unk == T & A11unk == T & Venetounk == T, T, F)) %>% # Combo effective while not solo?
    mutate(Veneto_A12_only = if_else(A11_A12_needed == F & Veneto_A12_needed == T, T, F)) %>%
    mutate(A11_A12_only = if_else(A11_A12_needed == T & Veneto_A12_needed == F, T, F)) %>%
    mutate(any_combo = if_else(A11_A12_needed == T & Veneto_A12_needed == T, T, F)) %>%
    mutate(A11_A12_needed_VI = if_else(A11_A12_unk == F & A12unk == T  & A11unk == T, T, F)) %>% # Veneto ind
    mutate(Veneto_A12_needed_A11I = if_else(Veneto_A12_unk == F & A12unk == T  & Venetounk == T, T, F)) %>%
    mutate(resistAllCombo = if_else(A11_15 == Inf & Veneto_15 == Inf, T, F)) %>% # Veneto ind
    mutate(resistA11_A12 = if_else(A11_15 == Inf, T, F)) %>% # Veneto ind
    mutate(resistVeneto_A12 = if_else(Veneto_15 == Inf, T, F))  # Vene
  
  compo <- fullcells2 %>%
    summarise(A11unk = mean(A11unk), Venetounk = mean(Venetounk), A12unk=  mean(A12unk),
              A11_A12_unk = mean(A11_A12_unk), Veneto_A12_unk = mean(Veneto_A12_unk),
              Veneto_only = mean(Veneto_only), A11_only = mean(A11_only), A12_only = mean(A12_only),
              Veneto_or_A11 = mean(Veneto_or_A11), Veneto_or_A12 = mean(Veneto_or_A12), A11_or_A12 = mean(A11_or_A12), any = mean(any),
              unkill_solo = mean(unkill_solo), A11_A12_needed = mean(A11_A12_needed), Veneto_A12_needed = mean(Veneto_A12_needed),
              Veneto_A12_only = mean(Veneto_A12_only), A11_A12_only = mean(A11_A12_only), any_combo = mean(any_combo), Veneto_A12_needed_A11I = mean(Veneto_A12_needed_A11I),
              A11_A12_needed_VI = mean(A11_A12_needed_VI), resistAllCombo = mean(resistAllCombo),
              resistA11_A12= mean(resistA11_A12) , resistVeneto_A12 = mean(resistVeneto_A12)
    ); compo
  
  # # Overall single treatment
  # compo %>%
  #   select(A11unk, Venetounk, A12unk ) %>%
  #   gather(key = "Pheno", value = "Pct") %>%
  #   mutate(Pct = (1-Pct) * 100) %>%
  #   arrange(desc(Pct))#
  # 
  # # Decomposition
  # compo %>%
  #   select(Veneto_only, A11_only, A12_only, Veneto_or_A11, Veneto_or_A12, A11_or_A12,   any, unkill_solo) %>%
  #   gather(key = "Pheno", value = "Pct") %>%
  #   mutate(Pct = Pct * 100) %>%
  #   arrange(desc(Pct))#%>%pull(Pct) %>% sum
  # 
  # # Resit
  # compo %>%
  #   select(Veneto_A12_only, A11_A12_only, any_combo) %>%
  #   gather(key = "Pheno", value = "Pct") %>%
  #   mutate(Pct = Pct * 100) %>%
  #   arrange(desc(Pct))#
  # 
  # 
  # compo %>%
  #   select(A11_A12_needed_VI, Veneto_A12_needed_A11I)
  # 
  
  library(eulerr)
  
  fullcells %>%
    mutate(A11 = if_else(A11_0 == Inf, F, T),
           Veneto = if_else(Veneto_0 == Inf, F, T),
           A12 = if_else(A11_15 == 0, T, F)) %>%
    select(A11, Veneto, A12) %>%
    mutate(resist = T)-> temp
  
  
  # if_else(!A11 & !Veneto & !A12, T, F)
  plotA <- plot(euler(combinations = temp %>%
                        rename(`A-1155463` = A11, `A-1210477` = A12, Venetoclax = Veneto, Resistant =resist)),  quantities = TRUE,
                fills = c("darkgreen", "blue", "red", "black"),alpha = 0.6,edges  = 1, main = "Single therapies" );plotA
  
  # self$reconstitute(drug = c(1,4)) %>%
  #   arrange(res)
  # library(VennDiagram)
  # library(RColorBrewer)
  # myCol <- brewer.pal(3, "Pastel2")
  #
  # # single treatment plot
  # overrideTriple=T
  #
  # any <-  compo$any
  # g = draw.triple.venn(area1 = (1-compo$Venetounk) * 100,
  #                      area2 = (1 - compo$A12unk)  * 100,
  #                      area3 = (1 - compo$A11unk)  * 100,
  #                      n12 =  (compo$Veneto_or_A12 + any)  * 100,
  #                      n23 = (any + compo$A11_or_A12) * 100,
  #                      n13 = (any + compo$Veneto_or_A11) * 100,
  #                      n123 = any * 100,
  #                      category = c(paste0("Venetoclax\n(", (1-compo$Venetounk) * 100, "%)"),
  #                                   paste0("A-1210477\n(", (1-compo$A12unk) * 100, "%)"),
  #                                   paste0("A-1155463 (" , (1-compo$A11unk) * 100, "%)")), fill = c("blue", "darkgreen", "red"), col = c("blue", "darkgreen", "red"),
  #                      , print.mode = c("raw"), sigdig=2, ind = T, scaled = T, euler.d = T, cat.col  = c("blue", "darkgreen", "red"))
  #
  # plotA <- plot_grid(grid.arrange(gTree(children=g)));plotA #, top="Single treatment"
  
  
  # pie_chart Veneto + A=12
  
  # Combo from where?
  fullcells2 %>%
    filter(Veneto_A12_needed_A11I) %>%
    select(cellid, Veneto_A12_unk, A11unk, A12unk, Venetounk) %>%
    # filter(!A11unk) %>%
    group_by(A11unk) %>%
    tally -> tempfora
  
  nametemp <- paste0( "Combo only\n(", tempfora$n[tempfora$A11unk], " resist mono;\n",
                      sum(tempfora$n) - tempfora$n[tempfora$A11unk], " sensitive A-11 )")
  
  data <- tribble(~group, ~ value, ~ group2,
                  "Veneto only", (compo$Veneto_only + compo$Veneto_or_A11 ) * 100, "Single",
                  "A-1210477 only",  (compo$A12_only + compo$A11_or_A12 ) * 100, "Single",
                  "Both single",  (compo$Veneto_or_A12 + compo$any ) * 100,"Single",
                  nametemp,  (compo$Veneto_A12_needed_A11I ) * 100, "Combo",
                  "Resistant", compo$resistVeneto_A12 * 100, "Combo") %>%
    arrange(desc(group))
  
  data <- data %>%
    mutate(prop = value / sum(data$value) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
    mutate(label = paste0(group,"\n(", value, "%)"))
  
  # plot()
  # plotB <-webr::PieDonut(fullcells2 %>% mutate(group =
  #                                                           case_when(Veneto_only | Veneto_or_A11 ~ "Veneto only",
  #                                                                     A12_only | A11_or_A12 ~ "A-1210477 only",
  #                                                                     Veneto_or_A12 | any ~ "Both single",
  #                                                                     Veneto_A12_needed_A11I ~ "Combo",
  #                                                                     T ~ "Resistant")
  # 
  # ),
  # aes(x = group,y = A11unk));plotB
  
  
  # plotB <- ggplot(data %>% mutate(group3 = if_else(grepl("Combo only", group), 1, 0)), aes(x="", y=value, fill=group)) +
  #   geom_bar(stat="identity", width=1, alpha= 0.6, col = "black") +
  #   
  #   # geom_segment(aes(x = , xend = 1, y = 0, yend = data$value[[1]]))+
  #   # geom_segment(aes(x = 1, xend = 1, y =data$value[1:3] %>% sum(), yend = 100))+
  #   # scale_color_manual(values = c("black","red"))
  #   coord_polar("y", start=0)+
  #   theme_void() +
  #   labs(fill = "Killed by", title = "Venetoclax + A-1210477")+
  #   theme(plot.title = element_text(hjust = 0.5))+
  #   geom_text(aes(y = ypos, label = value), color = "black", size=4)+
  #   # ggtitle("Venetoclax & A-1210477") +
  #   scale_fill_manual(values = c("red", "deeppink3", "darkgrey", "black", "blue"))+
  #   theme(plot.title = element_text(hjust = 0.5)); plotB
  
  
  plotB <- ggplot(data %>% mutate(group3 = if_else(grepl("Combo only", group), 1, 0))) +
    geom_point(aes(x = 0.5, y = 0.5, shape = "Covered by\nmonotherapies" ), size = 8, alpha= 1)+
    
    geom_bar( aes(x=1, y=value, fill=group, alpha = group2), stat="identity", width=1,  col = "black") +
      
    geom_segment(data=tibble(x =seq(0.5,1.5,0.1)[-6]), aes(x = x, xend = x, y = 0, yend = data$value[[1]], lty = "Covered by\nMonotherapy"))+
    geom_segment(data=tibble(x =seq(0.5,1.5,0.1)[-6]), aes(x = x, xend = x, y =data$value[1:3] %>% sum(), yend = 100))+
       # geom_segment(data=tibble(y =c(seq(0,50,5))), aes(x = 0.5, xend = 1.5, y = y, yend = y))+
    # scale_color_manual(values = c("black","red"))
     scale_alpha_manual(values = c(1,0.6))+
    coord_polar("y", start=0)+
    scale_shape_manual(values = 1)+
    theme_void() +
    labs(fill = "Killed by", title = "Venetoclax + A-1210477", lty = "", shape = "")+
    guides(lty = F, alpha = F)+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_text(aes(x = 1, y = ypos, label = value, fill=group), color = "black", size=4)+
    # ggtitle("Venetoclax & A-1210477") +
    scale_fill_manual(values = c("red", "deeppink3", "darkgrey", "grey26", "blue"))+
    theme(plot.title = element_text(hjust = 0.5)); plotB
  
  
  # pie_chart A11 + A=12
  
  # fullcells2 %>%
  #   filter(A11_A12_needed_VI) %>%
  #   select(cellid, Veneto_A12_unk, A11unk, A12unk, Venetounk) %>%
  #   # filter(!A11unk) %>%
  #   group_by(Venetounk) %>%
  #   tally -> tempfora
  
  nametemp <- paste0( "Combo only\n(", tempfora$n[tempfora$Venetounk], " resist mono;\n",
                      sum(tempfora$n) - tempfora$n[tempfora$Venetounk], " sensitive A-11 )")
  
  data2 <- tribble(~group, ~ value, ~ group2,
                   "A-1155463 only", (compo$A11_only + compo$Veneto_or_A11 ) * 100, "Single",
                   "A-1210477 only",  (compo$A12_only + compo$Veneto_or_A12 ) * 100, "Single",
                   "Both single",  (compo$A11_or_A12 + compo$any ) * 100,"Single",
                   nametemp,  (compo$A11_A12_needed_VI ) * 100, "Combo",
                   "Resistant", compo$resistA11_A12 * 100, "Combo") %>%
    arrange(desc(group))
  
  data2 <- data2 %>%
    mutate(prop = value / sum(data$value) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
    mutate(label = paste0(group,"\n(", value, "%)"))
  
  
  plotC <- ggplot(data2) +
    geom_point(aes(x = 0.5, y = 0.5, shape = "Covered by\nmonotherapies" ), size = 8, alpha= 1)+
    
    
    geom_bar(aes(x=1, y=value, fill=group, alpha = group2), stat="identity", width=1, col ="black") +
    geom_segment(data=tibble(x =seq(0.5,1.5,0.1)[-6]), aes(x = x, xend = x, y = data2$value[1:2] %>% sum, yend = 100, lty = "Covered by\nMonotherapy"))+
    # geom_segment(data=tibble(x =seq(0.5,1.5,0.1)[-6]), aes(x = x, xend = x, y =data2$value[1:3] %>% sum(), yend = 100))+
      scale_alpha_manual(values = c(1,0.6))+
      coord_polar("y", start=0)+
      scale_shape_manual(values = 1)+
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Killed by", title = "A-1155463 + A-1210477")+
    theme(plot.title = element_text(hjust = 0.5))+
    # ggtitle("A-1155463 & A-1210477")+
    labs(fill = "Killed by", lty = "", shape="")+
    guides(lty = F, alpha = F)+
    geom_text(aes(x = 1, y = ypos, label = value), color = "black", size=4) +
    scale_fill_manual(values = c("chartreuse4", "red","darkgoldenrod4", "darkgrey", "grey26"));plotC
  

  # plot_grid(plotA, plot_grid(plotB,
  # plotC, ncol =  1, labels = c("B", "C")), labels = c("A", "", ""))
  
  plot_grid(plotA, plotB, plotC, labels = labels, nrow = 1)
}

# plot_grid(compositionplot(sudh), 
#           compositionplot(karpas, labels = c("D", "E", "F")), ncol = 1)
# 
# 
# # plotC, ncol =  1, labels = c("B", "C")), labels = c("A", "", ""))
# # combo
# # combo treatment drug
# # g = draw.triple.venn(area1 = 98, area2 = 100-17, area3 = 98, n12 = 0, n23 = 45, n13 = 14, n123 = 0,
# #                      category = c("Venetoclax\n78%", "A12\n29%", "A11 26%"), fill = c("blue", "darkgreen", "red"), col = c("blue", "darkgreen", "red"),
# #                      , print.mode = c("raw"), sigdig=2, ind = T, scaled = T, euler.d = T, cat.col  = c("blue", "darkgreen", "red"))
# # grid.arrange(gTree(children=g))
# 
# plot_grid(grid.arrange(gTree(children=g)),grid.arrange(gTree(children=g)))
# 

plot_grid(compositionplot(karpas), compositionplot(sudh, labels = c("D", "E", "F")), ncol = 1)

uniquevalu <- unique(c(fullcells$A11_0, fullcells$A11_10))

crossing( fullcells %>% select(cellid) %>% rowid_to_column(), uniquevalu) %>% 
  left_join(fullcells %>%  select(cellid, A11_0)) %>% 
  mutate(res = if_else(uniquevalu >= A11_0,T, F)) %>% 
  mutate(uniquevalu = paste0("A11", uniquevalu)) %>% 
  unique() %>% 
  spread(key = uniquevalu, value = res) -> temp

crossing( fullcells %>% select(cellid) %>% rowid_to_column(), uniquevalu) %>% 
  left_join(fullcells %>%  select(cellid, Veneto_0)) %>% 
  mutate(res = if_else(uniquevalu >= Veneto_0,T, F)) %>% 
  mutate(uniquevalu = paste0("Veneto", uniquevalu)) %>% 
  unique() %>% 
  spread(key = uniquevalu, value = res) -> temp2

fullcells %>%  select(cellid, A11_0, A11_5, A11_10, A11_15) %>% gather("key", "value", -cellid ) %>% 
  filter(value == 0) %>% 
  mutate(A12_death = gsub("A11_", "", key) %>% as.double()) %>% 
  distinct(cellid, A12_death) %>% 
  arrange(A12_death) %>% 
  group_by(cellid) %>% 
  slice(1) -> A12deathlimit

crossing( fullcells %>% select(cellid) %>% rowid_to_column(), uniquevalu = c(0,5,10,15)) %>% 
  left_join(A12deathlimit) %>%
  mutate(A12_death=if_else(is.na(A12_death), Inf, A12_death)) %>% 
  mutate(res = if_else(uniquevalu >= A12_death,T, F)) %>% 
  mutate(uniquevalu = paste0("A12", uniquevalu)) %>% 
  select(-A12_death) %>% 
  unique() %>% 
  spread(key = uniquevalu, value = res) -> temp3

  
final <- temp %>% left_join(temp2)%>% left_join(temp3) %>% 
  select(-rowid, -cellid, - A11_0, -Veneto_0) %>% 
  select(-A110.16, - A110.64, -A1110, - A110, -Veneto0, -Veneto0.16, -Veneto0.64, -Veneto10) %>% 
  select(A110.08,A111.3, A1120, A115, A11Inf, Veneto0.08,Veneto1.3,  Veneto5, Veneto20);final


plot(euler(final , shape = "ellipse"),  quantities = TRUE)

tibble(A11 = Veneto)
fullcells

fullcells %>%
  mutate(A11 = if_else(A11_0 == Inf, F, T),
         Veneto = if_else(Veneto_0 == Inf, F, T),
         A12 = if_else(A11_15 == 0, T, F)) %>%
  select(A11, Veneto, A12) %>%
  mutate(resist = if_else(!A11 & !Veneto & !A12, T, F), Allcells = T) -> temp

# plotA <- 
  
  plot(euler(shape = "ellipse", combinations = temp %>%
                      rename(`A-1155463` = A11, `A-1210477` = A12, Venetoclax = Veneto, Resistant =resist)),  quantities = TRUE,
              fills = c("red", "blue", "darkgreen", "black"), alpha = 0.6,edges  = 1)

