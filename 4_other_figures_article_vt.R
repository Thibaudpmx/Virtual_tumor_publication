source(here('2_Virtual_Tumor_objects.R'))
theme_set(theme_bw())
# create_VT_Project("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg")

# Author: Thibaud Derippe
# From CPT PSP guideline: Accepted figure files include PDF, TIFF,
# and EPS in PC format, preferably from Photoshop or Illustrator software
# Resolution 125 dpi (dots per inch)
# let's do pdf!


# Some functions

plotFinal2D <- function(cellline, title = "", xlabb = NA){

  main_drug <- cellline$drugs[[1]]
  drug_name = if_else(main_drug ==2, 'A-1155463', 'Venetoclax' )

  OFS <-   cellline$OF(drug = cellline$drugs, detail = T) %>%
    slice(1:4) %>%
    conc4cor()

if(is.na(xlabb)) xlabb <-  paste0("Concentration ", drug_name," (µM)")

  cellline$reconstitute() %>%
    mutate(conc = !!parse_expr(paste0('conc', main_drug))) %>%
    mutate(drugs = drug_name) %>%
    mutate(Reconst= if_else(Reconst == "res", " Recon-\nstructed", "Observed")) %>%
    conc4cor() %>%
    {temp <<- .} %>%

    ggplot()+
    geom_line(aes(conc, res, col = Reconst))+
    geom_point(aes(conc, res, col = Reconst))+
    geom_ribbon(data = temp  %>% spread(key = Reconst, value = res), aes(x = conc, ymin = Observed, ymax =` Recon-\nstructed`), alpha = 0.3)+
    facet_wrap(~conc4)+
    geom_text(data = OFS, aes(1, 95, label = paste0("Area = ", round(area))))+
    labs(x = xlabb, y = "Cell Viability (%)", title =title, col = "")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_log10()
}



# Fig1: used data Datas -----------------------------------------------------------

data_VT %>%
  filter(Drug %in% c(1,2)) %>%
  mutate(conc = if_else(Drug == 1, conc1, conc2)) %>%
  mutate(conc42 = factor(paste(conc4, " µM"))) %>%
  mutate(Drug1 = if_else(Drug1 == "Venetoclax",  "Venetoclax\n(Bcl2 inhibitor)", "A-1155463\n(Bclxl inhibitor)")) %>%
  ggplot()+
  geom_line(aes(conc, Value, col = fct_reorder(conc42, conc4)))+
  geom_point(aes(conc, Value, col = fct_reorder(conc42, conc4)))+
  facet_grid(Cell_line  ~ Drug1)+
  scale_x_log10()+
  theme_bw()+
  labs(x = "Concentration (µM)", col = "  A-1210477\n(Mcl1 inhibitor)",y= "Cell Viability (%)") -> pk_data

read.table(here('data/mice_SU_DHL4_full.csv'), header = T, sep = ";") %>%
  as.tibble %>%
  mutate(group = case_when(Drug == 0 ~ "Vehicle",
                           Drug == 1 ~ "Venetoclax\n(Bcl2 inhibitor)\n",
                           Drug == 2 ~ "A-1592668\n(Mcl1 indirect\ninhibitor)\n",
                           T ~ "Combination")) %>%
  ggplot()+
  geom_line(aes(TIME, DV, col = fct_reorder(group, Drug)))+
  geom_point(aes(TIME, DV, col = fct_reorder(group, Drug)))+
  labs(x = "Time (days)", y = "Tumor volume (mm3)", col = "", lty ="")+
  scale_linetype_manual(values = 2)+
  theme_bw()+
  geom_vline(data = tibble(x = c(1,22)), aes(xintercept = x, lty = "QD\nTreatment\nWindow"))+
  facet_wrap(~"SU-DHL-4")+
  # geom_rect(aes(xmin = 1, xmax = 22, ymin = 0, ymax = Inf), alpha = 0.005)+
  scale_y_log10() -> pd_data


plot_grid(pk_data, pd_data,labels  = c("A", "B"))



# Fig 2: Sampling system ---------------------------------------------------------


# Plot A1
ploaa <- crossing(Bcl20 = seq(20,1020,112), Bclxl0 = seq(50,1500, 170)) %>%
  ggplot()+
  geom_point(aes(Bcl20, Bclxl0))+
  geom_hline(aes(yintercept = Bclxl0))+
  geom_vline(aes(xintercept = Bcl20))

ploaa2 <- crossing(Mcl10 = seq(2,150, 16), BIM0 = seq(0,200,25)) %>%
  ggplot()+
  geom_point(aes(Mcl10, BIM0))+
  geom_hline(aes(yintercept = BIM0))+
  geom_vline(aes(xintercept = Mcl10))

ploaa3 <- crossing( PUMA0 = seq(0,200,25), NOXA0 =  seq(0,200,25)) %>%
  ggplot()+
  geom_point(aes(PUMA0, NOXA0))+
  geom_hline(aes(yintercept = NOXA0))+
  geom_vline(aes(xintercept = PUMA0))

ploaB <-  crossing(BAK0 = c(0,500,1000), BAXc0 = c(0,500,1000)) %>%
  filter(!(BAK0 == 0 & BAXc0 == 0)) %>%
  ggplot()+
  geom_point(aes(BAK0, BAXc0))+
  geom_hline(aes(yintercept = BAXc0))+
  geom_vline(aes(xintercept = BAK0))

plotA1 <- plot_grid(ploaa, ploaa2, ploaa3, ploaB, nrow = 2);plotA1

# plot A2
tibble(param = c("Bcl20", "Bclxl0", "Mcl10", "BIM0","PUMA0", "NOXA0"), min = 0, max = c(2000,2000,500,1000,1000,1000)) %>%
  mutate(sample=  map2(min, max,  ~ runif(n = 2000, min = .x, max = .y) )) %>%
  unnest() %>%
  ggplot()+
  # geom_density(aes(sample))+
  geom_histogram(aes(sample))+
  facet_wrap(~param, scales = "free")+
  scale_x_continuous(breaks = c(0,500,1000,2000))+

  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Sampled parameter values")+
  labs(title = "10,000 to 50,000 sampling")  ->plotsamplA; plotsamplA


tibble(param = c("BAK0", "BAX0"), sample = c(543,1323)) %>%
  # mutate(sample=  map2(min, max,  ~ runif(n = 1, min = .x, max = .y) )) %>%
  # # unnest() %>%
  ggplot()+
  # stat_function(fun = function(x) 0.001)+
  # geom_density(aes(sample),alpha = 0.1)+
  geom_histogram(aes(sample), size = 200)+
  geom_rect(aes(xmin = sample -100, xmax = sample +100, ymin = 0, ymax = 1))+
  # geom_vline(data = tibble(param = c("BAK0", "BAX0"), value = c(200,600)), aes(xintercept = value), col = "red", size = 10)+
  facet_wrap(~param, ncol = 1, scales = 'free')+
  scale_x_continuous(breaks = c(0,543, 1323, 2000))+
  coord_cartesian(xlim = c(0,2000))+
  scale_y_continuous(breaks = c(0,1))+
  theme(axis.title.y=element_blank())+
  labs(title = "Unique", x = "Sampled") ->plotsamplB; plotsamplB

plotgridA <- plot_grid(plotsamplA, plotsamplB, rel_widths = c(3,1))




plotgridC <- plot_grid(plotA1, plotgridA, ncol = 1, labels = c("A.1", "A.2"));plotgridC


# Colomn C


values <- c(8,4,2,1,10,9,9,8)
makeprofile <- function(values= c(8,4,2,1,10,9,9,8)){
  temp <- crossing(Conc1 = c(unique(data_VT$conc1)), drug = c("A-1155463", "Venetoclax"), Conc2 = c(0,5,10,15)) %>%
    mutate(Fate = "Death")

  if(values[[1]] > 0) temp$Fate[temp$drug == "A-1155463" & temp$Conc2 == 0][1:values[[1]]] <- "Survival"
  if(values[[2]] > 0)  temp$Fate[temp$drug == "A-1155463" & temp$Conc2 == 5][1:values[[2]]] <- "Survival"
  if(values[[3]] > 0)  temp$Fate[temp$drug == "A-1155463" & temp$Conc2 == 10][1:values[[3]]] <- "Survival"
  if(values[[4]] > 0) temp$Fate[temp$drug == "A-1155463" & temp$Conc2 == 15][1:values[[4]]] <- "Survival"

  if(values[[5]] > 0) temp$Fate[temp$drug != "A-1155463" & temp$Conc2 == 0][1:values[[5]]] <- "Survival"
  if(values[[6]] > 0) temp$Fate[temp$drug != "A-1155463" & temp$Conc2 == 5][1:values[[6]]] <- "Survival"
  if(values[[7]] > 0) temp$Fate[temp$drug != "A-1155463" & temp$Conc2 == 10][1:values[[7]]] <- "Survival"
  if(values[[8]] > 0) temp$Fate[temp$drug != "A-1155463" & temp$Conc2 == 15][1:values[[8]]] <- "Survival"
  temp
}


cell1 <- makeprofile()
cell2  <- makeprofile( c(8,4,4,0,4,2,1,0))

test <- bind_rows(cell1, cell2, .id = "name") %>%
  mutate(name = if_else(name == 1, "Cell 1 (Bcl2 = 143,...)","Cell 2 (Bcl2 = 341,...)" ))

test %>%
  # mutate(testcol = sample(c("a", "b","c"), 1:nrow(test),replace = T)) %>%
  ggplot()+
  # facet_wrap(~profile, scales = "free", ncol = 1)+
  geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black")+
  geom_rect(data =test %>%
              filter(drug != "Venetoclax") %>% group_by(name) %>% slice(1), aes(xmin = 0.5, xmax = 10.5, ymin = 0.5, ymax = 1.5), col = "blue", alpha = 0.2, size = 2)+

  # geom_rect(data =test %>%
  # filter(drug == "Venetoclax" & name == "Cell 2 (Bcl2 = 341, Mcl1 = 68 ...)"), aes(xmin = 3.5, xmax = 5.5, ymin = 0.5, ymax = 1.5), col = "black", alpha = 0, size = 2)+
  geom_text(data =test %>%
              filter(drug != "Venetoclax") %>% group_by(name) %>% slice(1), aes(5, 1, label = "Same bag\n\"A-1155463 alone\""), col = "blue", size = 3)+
  scale_y_discrete()+ #•labels = temp$delta
  facet_grid(drug ~ name)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  # scale_x_discrete(labels =)+
  labs(x = "Concentration main drug (µM)", y = "Concentration A-1210477 (µM)" ) -> plot_simulations;plot_simulations


# Colomn C

data_VT %>%
  filter(Drug %in% c(1,2)) %>%
  filter(Cell_line_bin == 1) %>%
  mutate(conc = if_else(Drug == 1, conc1, conc2)) %>%
  mutate(conc42 = factor(paste(conc4, " µM"))) %>%
  mutate(Drug1 = if_else(Drug1 == "Venetoclax",  "Venetoclax\n(Bcl2 inhibitor)", "A-1155463\n(Bclxl inhibitor)")) -> tempplot

tempplot %>%
ggplot()+
  geom_line(aes(conc, Value, col = fct_reorder(conc42, conc4)))+

  facet_wrap(  ~ Drug1, ncol =  1)+
  scale_x_log10(breaks = unique(data_VT$conc1))+
  theme_bw()+
  geom_line(data = tempplot %>%
              filter(!(conc4 == 0 & ! conc2 %in% c(5,10) )) %>%
              filter(!(conc4 == 5 & ! conc2 %in% c(0.32,0.64) ))%>%
              filter(!(conc4 == 10 & ! conc2 %in% c(0.08,0.16) ))%>%
              filter(!(conc4 == 15 & ! conc2 %in% c(0,0.08) )) %>%
              filter(Drug == 2),
            aes(conc, Value, col = fct_reorder(conc42, conc4)), size = 3 )+
  geom_line(data = tempplot %>%
              filter(!(conc4 == 0 & ! conc1 %in% c(Inf) )) %>%
              filter(!(conc4 == 5 & ! conc1 %in% c(10,20) ))%>%
              filter(!(conc4 == 10 & ! conc1 %in% c(10,20) ))%>%
              filter(!(conc4 == 15 & ! conc1 %in% c(5,10) )) %>%
              filter(Drug == 1),
            aes(conc, Value, col = fct_reorder(conc42, conc4)), size = 3 )+
  # geom_rect(data =tempplot %>% filter(Drug == 1) %>% slice(1),
  #           aes(xmin = 0.3,xmax = 0.66, ymin = 85, ymax = 95), alpha =0, col = "black", size = 2)+
  geom_line(data = tempplot %>%
              filter(!(conc4 == 0 & ! conc1 %in% c(0.32,0.64) )) %>%
              filter(!(conc4 == 5 & ! conc1 %in% c(0.08,0.16) ))%>%
              filter(!(conc4 == 10 & ! conc1 %in% c(0,0.08) ))%>%
              filter(!(conc4 == 15 & ! conc1 %in% c(0) )) %>%
              filter(Drug == 1),
            aes(conc, Value, group = conc4, alpha = "" ), col ="black", size = 1, lty = 1 )+
  # geom_line(data = tempplot %>% filter(Drug == 1 & conc4 == 0 & conc1 == 20) %>%
  #             crossing(d = 1:2) %>%
  #             mutate(Value = if_else(d == 2, 0, Value)) %>%
  #             mutate(conc = if_else(d == 2, Inf, conc)),
  #            aes(conc, Value, col = conc42), lty = 2) +
  scale_alpha_manual(values = 1)+
  # geom_segment(data = tibble(Drug1 =  "A-1155463\n(Bclxl inhibitor)"),col = "blue", alpha = 0.2,
               # aes(x =8 , y = 90, xend = 10,yend = 100 ),arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))+
  # geom_text(data = tibble(Drug1 =  "Venetoclax\n(Bcl2 inhibitor)"), aes(x = 25,y= 10, label = "?"), col = "red", alpha = 0.8)+
  geom_line(data = tempplot %>%
              filter(!(conc4 == 0 & ! conc2 %in% c(5,10) )) %>%
              filter(!(conc4 == 5 & ! conc2 %in% c(0.32,0.64) ))%>%
              filter(!(conc4 == 10 & ! conc2 %in% c(0.32,0.64) ))%>%
              filter(!(conc4 == 15 & ! conc2 %in% c(0) )) %>%
              filter(Drug == 2),
            aes(conc, Value, group = conc4 ), col ="black", size = 1 )+

  geom_point(aes(conc, Value, col = fct_reorder(conc42, conc4)))+
  geom_point(data = tempplot %>%
               filter(conc4 == 15 & conc1 == 0 & conc2 == 0),
             aes(conc, Value), col ="black", size = 4 )+
  # geom_rect(aes(xmin = 0.5, xmax = 1, ymin = 10, ymax = 15))+
  # geom_segment(aes(x = 5, xend = 10, y = 90, yend = 80 ), col = "red")+
  labs(x = "Concentration main drug (µM)", alpha = "Cell2\ndeath\nzones", col = "Cell1 death\nzones for \nA-1210477:",y= "Cell Viability (%)")  -> plotpkonly;plotpkonly




#
# plot_grid(plot_simulations,plotpkonly, labels = LETTERS )

#final plooooooooooooot

pdf('figures/fig2.pdf', width = 16, height = 7)
plot_grid(plotgridC, plot_simulations,plotpkonly, labels = c("", "B", "C"), nrow = 1 )
dev.off()


# Some bags -------------------------------------------------------------------


map(c(1:10), ~makeprofile( c(.x,4,4,0,4,2,1,0))) %>%
      bind_rows( .id = "name") %>%
  filter(name %in% c(1,2,3,8,9,10)) %>%
  # mutate(name = as.double(name)) %>%
  # mutate(name = as.factor(name)) %>%
  filter(Conc2 == 0, drug != "Venetoclax") %>%
  ggplot()+
  # facet_wrap(~profile, scales = "free", ncol = 1)+
  geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black")+
  scale_y_discrete()+ #•labels = temp$delta
  facet_wrap(~ fct_reorder(name, as.double(name)))+
  labs(x = "Concentration main drug (µM)", y = "A-1210477 (µM)" ) -> plotaa;plotaa



map(c(6:8), ~makeprofile( c(8,.x,4,0,4,2,1,0))) %>%
  bind_rows( .id = "name") %>%
  bind_rows(  map_dfr(c(8), ~makeprofile( c(8,.x,4,1,4,2,1,0))) %>% mutate(name = "4")) %>%
  # filter(name %in% c(1,2,3,8,9,10)) %>%
  # mutate(name = as.double(name)) %>%
  # mutate(name = as.factor(name)) %>%
  filter(drug != "Venetoclax") %>%
  ggplot()+
  # facet_wrap(~profile, scales = "free", ncol = 1)+
  geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black")+
  scale_y_discrete()+ #•labels = temp$delta
  facet_wrap(~ fct_reorder(name, as.double(name)))+
  labs(x = "Concentration main drug (µM)", y = "A-1210477 (µM)" ) -> plotab;plotab

makeprofile( c(8,8,4,0,1,0,0,0)) %>% mutate(name = "1")  %>%
  bind_rows(makeprofile( c(8,8,4,0,10,8,6,5)) %>% mutate(name = "2") ) %>%
  ggplot()+
  # facet_wrap(~profile, scales = "free", ncol = 1)+
  geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black")+
  scale_y_discrete()+ #•labels = temp$delta
  facet_grid(drug ~ name)+

  labs(x = "Concentration main drug (µM)", y = "A-1210477 (µM)" ) -> plotac;plotac

plot_grid(plotaa, plotab, plotac, ncol = 1)


plot_grid(plotpkonly , plot_simulations, plot_grid(plotaa, plotab, plotac, ncol = 1), nrow = 1 )

# Fig3 Curve decomposition -----------------------------------------------------


# First try with arrow middle segment
data_VT %>%
  filter(Drug == 1, Cell_line_bin == 2) %>%
  filter(conc4 == 0 ) %>%
  select(conc1, Value) %>%
  mutate(Value = if_else(Value > 100,100, Value)) %>%
  mutate(conc1lag = lag(conc1), Valuelag = lag(Value)) %>%
  mutate(midConc = exp((log(conc1) + log(conc1lag)) / 2)) %>%
  mutate(midConc = if_else(midConc == 0, ((log(0.04) + log(0.08))/2) %>% exp(), midConc)) -> decompo

plotA <-decompo %>%
  mutate(forcol = case_when(conc1 == 0.08 ~ "blue",
                            conc1 == 0.16 ~ "chocolate",
                            conc1 == 20 ~ "darkgreen",
                            T ~ "black")) %>%
 ggplot()+
  geom_line(aes(conc1, Value))+
  geom_point(aes(conc1, Value))+
  scale_x_log10(breaks = unique(data_VT$conc1))+
  geom_rect(aes(xmin = conc1, xmax = conc1lag, ymin = Valuelag, ymax = Value ,fill = forcol), alpha = 0.2)+
  coord_cartesian(ylim = c(0,100)) +
  geom_segment(aes(x = midConc, xend = midConc, y = Value, yend = Valuelag, col = forcol),
               arrow=arrow(ends='both', length = unit(0.10,"cm"), type = "closed"))+
  geom_text(aes(x=midConc, y = (Value + Valuelag) /2, label = (Valuelag - Value ) %>% round %>%
                  paste0("%"), col = forcol),
            nudge_x = 0.2, nudge_y = 1 )+
  geom_segment(aes(x = 0, xend = 0, y = 95.9, yend = 100),
               arrow=arrow(ends='both', length = unit(0.10,"cm"), type = "closed"), col = "red")+
  geom_rect(data = mtcars %>% slice(1), aes(xmin = 20, xmax = Inf, ymin = 18, ymax = 0), fill = "red",alpha = 0.2)+
  geom_segment(aes(x = 30, xend = 30, y =  18.0, yend = 0),
               arrow=arrow(ends='both', length = unit(0.10,"cm"), type = "closed"), col = "red")+
  geom_text(data = mtcars %>% slice(1), aes(x=30, y =9, label =18 %>%
                  paste0("%")),
            nudge_x = 0.15, col ="red")+
  scale_color_manual(values = c("black", "blue", "chocolate", "darkgreen"))+
  scale_fill_manual(values = c("black", "blue", "chocolate", "darkgreen"))+
  labs(x = "Drug Concentration (µM)", y = "Cell Viability (%)")+
  guides(col = F, fill = F); plotA

# decompo %>%
# ggplot()+
#   geom_line(aes(conc1, Value))+
#   geom_point(aes(conc1, Value))+
#   scale_x_log10(breaks = unique(data_VT$conc1))


# Other plots
notouched <- VT2$new(filter_data = Drug %in% 1 & Cell_line_bin ==2 , pen_nbag = 10,  name = "test")

# plotD <- notouched$plot() + labs( x = "Concentration (µM)",title = NULL);plotD

plotD <- plotFinal2D(notouched, xlabb = 'Concentration (µM)');plotD

touched <- VT2$new(filter_data = Drug %in% 1 & Cell_line_bin == 2 , pen_nbag = 10,  name = "test")

touched$optim()
plotF <- plotFinal2D(touched, xlabb = 'Concentration (µM)');plotF



notouched2 <- VT2$new(filter_data = Drug %in% 1 & Cell_line_bin ==2 , pen_nbag = 10,  name = "test")


notouched2$curve_sampling %>%
  filter(conc4 == 0) %>%
  mutate(a = map2(n, sample, function(x, y) sample(y,x, replace = T))) %>%
  select(a) %>%
  unnest() %>% pull() -> harvesttemp

 notouched2$harvest <- harvesttemp
 plotC <- notouched2$reconstitute() %>%
  mutate(Reconst = if_else(Reconst == "res", " Recon-\nstructed", "Observed")) %>%
   conc4cor() %>%
  ggplot()+
  geom_line(aes(conc1, res, col = Reconst))+
  geom_point(aes(conc1, res, col = Reconst))+

  facet_wrap(~conc4)+
  scale_x_log10()+
    labs(x = "Concentration (µM)", y = "Cell viability (%)", col = "");plotC



 pdf('figures/fig3.pdf', width = 16, height = 7)
 plot_grid(plotA, ggplot,plotC,
           plotD,ggplot,plotF, labels = LETTERS, nrow = 2 )
 dev.off()


# Final Plot (before powerpoint)

# notouched2$plot() + labs( x = "Concentration (µM)",title = NULL)
#
# temp <- decompo %>%
#   mutate(delta = (Valuelag - Value) %>% round) %>%
#   mutate(delta = if_else(is.na(delta), (100 - Value) %>% round, delta)) %>%
#   rename(conc = conc1) %>%
#   select(delta) %>%
#   rowid_to_column("profile") %>%
#   mutate(profile = profile - 1) %>%
#   add_row(profile = 11, delta = 18)
#
# map(0:10, ~tibble(conc = c(unique(data_VT$conc1), Inf), death = c(rep(F, .x), rep(T, 11 - .x)),
#                   profile = .x)) %>%
#   bind_rows() %>%
#   left_join(temp) %>%
#   ggplot()+
#   # facet_wrap(~profile, scales = "free", ncol = 1)+
#   geom_tile(aes(factor(conc), factor(profile), fill = death), col  = "black")+
#   scale_y_discrete(labels = temp$delta)
#
#
# data_VT %>%
#   filter(Drug == 1, Cell_line_bin == 2) %>%
#   filter(conc4 == 0 ) %>%
#   select(conc1, Value) %>%
#   mutate(Value = if_else(Value > 100,100, Value)) %>%
#   mutate(conc1lag = lag(conc1), Valuelag = lag(Value)) %>%
#   mutate(midConc = exp((log(conc1) + log(conc1lag)) / 2)) %>%
#   mutate(midConc = if_else(midConc == 0, ((log(0.04) + log(0.08))/2) %>% exp(), midConc)) %>%
#   ggplot()+
#   geom_line(aes(conc1, Value))+
#   geom_point(aes(conc1, Value))+
#   scale_x_log10()+
#   coord_cartesian(ylim = c(0,100)) +
#   # geom_segment(aes(x = conc1lag, xend = conc1 , y = Valuelag +  2 , yend = Value + 2 ),
#                # arrow=arrow( length = unit(0.10,"cm"), type = "closed"), col = "red") +
#   geom_text(aes(x=midConc, y = (Value + Valuelag) /2, label = (Valuelag - Value ) %>% round),
            # nudge_x = 0.1 )



# Fig4: Final Virtual tumors 3D----------------------------------------------------

cell_line2 <- readRDS(here( "calibrated_VT/VT_both_cell_line_2.RDS" ))
cell_line1 <- readRDS(here( "calibrated_VT/VT_both_cell_line_1.RDS" ))

cell_line2$plot()

conc4cor <- function(df){

  df %>%
    mutate(conc4bis = conc4) %>%
    mutate(conc4 = paste0('[A-12] ', conc4, ' µM')) %>%
    mutate(conc4 = fct_reorder(factor(conc4), conc4bis)) %>%
    select(-conc4bis)

}

plotFinal <- function(cellline, title = ""){

OFS <-   cellline$OF(drug = c(1,4), detail = T) %>%
    slice(1:4) %>%
    mutate(drugs = "Venetoclax") %>%
    bind_rows(  cellline$OF(drug = c(2,4), detail = T) %>%
                  slice(1:4) %>%
                  mutate(drugs = "A-1155463")  ) %>%
  conc4cor()

cellline$reconstitute() %>%
  mutate(conc = if_else(is.na(conc2), conc1, conc2)) %>%
  mutate(drugs = if_else(drugs == "1_4", "Venetoclax", "A-1155463")) %>%
  mutate(Reconst= if_else(Reconst == "res", " Recon-\nstructed", "Observed")) %>%
  conc4cor() %>%
  {temp <<- .} %>%

  ggplot()+
  geom_line(aes(conc, res, col = Reconst))+
  geom_point(aes(conc, res, col = Reconst))+
  geom_ribbon(data = temp %>% select(-conc1, -conc2) %>% spread(key = Reconst, value = res), aes(x = conc, ymin = Observed, ymax =` Recon-\nstructed`), alpha = 0.3)+
  facet_grid(drugs~conc4)+
  geom_text(data = OFS, aes(1, 95, label = paste0("Area = ", round(area))))+
  labs(x = "Concentration A-1155463 or Venetoclax (µM)", y = "Cell Viability (%)", title =title, col = "")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_log10()
}

# Figure 4
pdf('figures/fig4.pdf',height = 7, width = 8)
plot_grid(
plotFinal(cell_line1, 'KARPAS-422'),
plotFinal(cell_line2, 'SU-DHL-4'), labels = c('A','B'), ncol = 1
)
dev.off()
shell.exec('figures/fig4.pdf')

cell_line1$reconstitute() %>%
  mutate(cell = 1) %>%
  bind_rows(cell_line2$reconstitute() %>% mutate(cell = 2 )) %>%
  mutate(conc = if_else(is.na(conc2), conc1, conc2)) %>%
  select(-conc1, -conc2) %>% spread(key = Reconst, value = res) %>%
  mutate(dif = abs(res - Value)) %>%
  pull(dif) %>% summary()





# cell_line1$reconstitute() %>%
#   mutate(cell = 1) %>%
#   bind_rows(cell_line2$reconstitute() %>% mutate(cell = 2 )) %>%
#   mutate(conc = if_else(is.na(conc2), conc1, conc2)) %>%
#   mutate(drugs = if_else(drugs == "1_4", "Venetoclax", "A-1155463")) %>%
#   # mutate(Reconst= if_else(Reconst == "res", " Recon-\nstructed", "Observed")) %>%
#   select(-conc1, -conc2) %>% spread(key = Reconst, value = res) %>%
#   ggplot()+
#   geom_point(aes(Value,res ))+
#   geom_abline()+
#   facet_grid(drugs~conc4)

# Fig S8 Final Virtual tumors 2D----------------------------------------------------

files <- c( here( "calibrated_VT/VT_A11_cell_line_1.RDS" ),
            here( "calibrated_VT/VT_A11_cell_line_2.RDS" ),
            here( "calibrated_VT/VT_Venetoclax_cell_line_1.RDS" ),
            here( "calibrated_VT/VT_Veneto_cell_line_2.RDS" )

)




allplots <- map2(files, c('A-1155463 | SU-DHL-4',
                          "A-1155463 | KARPAS-422",
                          "Venetoclax | SU-DHL-4",
                          "Venetoclax | KARPAS-422"
                         ), ~plotFinal2D(readRDS(.x), .y))

pdf('figures/S8.pdf', width = 11, height = 8)
invoke(plot_grid, allplots, labels = LETTERS )
dev.off()




# How many bags -----------------------------------------------------------
possiblevalues <- 0:10
twoD <- crossing(A = possiblevalues, B =  possiblevalues, C =  possiblevalues, D =  possiblevalues) %>%
  filter(A >= B, B >= C, C>= D)
twoD2 <- twoD
names(twoD2) <- paste0(names(twoD), "2")

crossing(twoD,twoD2) %>%
  filter(!(A > 0 & A2 ==0)) %>%
  filter(!(B > 0 & B2 ==0)) %>%
  filter(!(C > 0 & C2 ==0)) %>%
  filter(!(D > 0 & D2 ==0)) %>%
  filter(!(A == 0 & A2 >0)) %>%
  filter(!(B == 0 & B2 >0)) %>%
  filter(!(C == 0 & C2 >0)) %>%
  filter(!(D == 0 & D2 >0))



# How many cells in total -------------------------------------------------

setwd(here("Virtual_Cells"))

all <- list.files()[grepl("\\.RDS$", list.files())]

map_dbl(all, ~readRDS(.x) %>% nrow()) %>% sum

setwd(here("Virtual_Cells",'one_per_bag_combined'))

tibble(file = list.files()) %>%
  mutate(nrow = map_dbl(file, ~ readRDS(.x) %>% nrow))



# Fig S5 Illustration penalty terms ----------------------------------------------
# Note: no need to made this example reproducible in Git
# Just try by yourself calibrate VT with different penalty terms

previous <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/3_virtual_tumors/SUDHL4-Venetoclax_pendelta_10.RDS")
previous2 <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/3_virtual_tumors/SUDHL4-Venetoclax.RDS")

previous$plot() # 6 + 26 + 5 + 34
test <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/new_VT_qspvp/VT_A11_cell_line_1.RDS")
test$celltheque <- previous$celltheque
test$harvest <- previous$harvest
plotFinal2D(test, xlabb = '')

previous2$plot() # 19 + 19 + 6 + 27 = 71




# Illustration algo -------------------------------------------------------

OFf <- function(df){
    df %>%
     group_split(group) %>%
    map(~.x %>%
          # x %>%
      rename(Value = obs) %>%
      mutate( conclag = lag(conc), reslag = lag(res), Valuelag = lag(Value)) %>%
      slice(-1) %>%
      mutate(area = pmap_dbl(list(Value,Valuelag, res, reslag, conc, conclag), function(Value,Valuelag, res, reslag, conc, conclag){

        # Value = 17.25; Valuelag =  46.69;  res = 18; reslag = 47; conc1 = 2.9957323; conc1lag = 2.302585

        h <- conc - conclag
        a <- sqrt(h ^2 + abs(reslag-Valuelag)^ 2)
        b <- sqrt(h ^2 + abs(res -Value)^ 2)
        area <- (a + b) * h /2

        if((Value > res) != (Valuelag > reslag)) area <- area / 2
        if(reslag == Valuelag  & res == Value ) area <- 0
        area
      }))
    ) %>%
    bind_rows()
}


tibble(obs = c(100, 70,50,30), res = c(100, 78,49,24), conc = c(0, 1,2,3)) %>%
  mutate(dif = obs-res) -> temp



tibble(obs = c(100, 70,50,30,0), res = c(100, 78,49,24,0), conc = c(0, 1,2,3,4), group = 1) %>%
  bind_rows(tibble(obs = c(100, 30,20,5,0), res = c(100, 15,10,9,0), conc = c(0, 1,2,3,4), group = 2)  ) %>%
  mutate(dif = obs-res) -> bothgroup

temp <- group1 <- bothgroup %>%
  filter(group == 1)

plotf <- function(temp ){

  temp %>%
    group_split(group) %>%
    map(~ .x %>%
    mutate(lagobs = lag(obs), lagres=lag(res)) %>%
    mutate(difobs = lagobs - obs, difres = lagres-res) %>%
    mutate(obsy = (lagobs + obs)/2, resy = (lagres + res)/2) %>%
    mutate(yboth = map2_dbl(obsy, resy, ~min(c(.x, .y)) - 5))
    ) %>%
    bind_rows() -> slopes

  temp %>%
  gather(a, value, obs, res) %>%
    mutate(group = paste0(group, a)) %>%
    # filter(a = "obs") %>%
  ggplot()+
    geom_line(aes(conc, value, col = a, group = group))+
    # facet_wrap(~a)
  geom_text(data = slopes, aes(x = conc-0.5, y = obsy, label =  difobs), col = "red")+
    geom_text(aes(x = 2, y = 95, label = paste0("OF = ", round(sum(OFf(temp)$area),1))))+
    geom_text(data = slopes, aes(x = conc-0.5, y = resy, label =difres), col = "blue")+
    geom_text(data = slopes, aes(x = conc-0.5, y = yboth , label = difobs  - difres), col = "black")+
    geom_segment(data = slopes, aes(x = conc-0.4, xend = conc -0.6, y = yboth + 2, yend = yboth +2), col = "black")+
    # geom_text(data = OFf(temp), aes(x = conc-0.5, y = reslag, label = round(area,1)))+
# +
  geom_point(aes(conc, value, col = a))+
  geom_segment(data = temp , aes(x = conc, xend = conc, y = obs, yend = res), lty = 1)+ # ,
  # arrow = arrow(ends = "both", length = unit(0.1, "inches"))
  theme_bw()+
    scale_x_continuous(breaks = 0:4, labels = c(0:3, "Inf"))+
    labs(x = "log(conc)", y = "Cell viability (%)", col = "")+
  geom_ribbon(data = temp, aes(x = conc, ymin = obs, ymax = res, group = group), alpha = 0.2 ) #+
  # geom_text(data = temp %>% slice(nrows), aes(x = conc - 0.2, y = (obs + res) /2, label = paste0(dif, " VCs\n]",conc-1, ";", conc,"]" )))
};plot1 <- plotf(bothgroup);plot1

plotf <- function(temp ){

  temp %>%
    group_split(group) %>%
    map(~ .x %>%
          mutate(lagobs = lag(obs), lagres=lag(res)) %>%
          mutate(difobs = lagobs - obs, difres = lagres-res) %>%
          mutate(obsy = (lagobs + obs)/2, resy = (lagres + res)/2) %>%
          mutate(yboth = map2_dbl(obsy, resy, ~min(c(.x, .y)) - 10))
    ) %>%
    bind_rows() -> slopes

  temp %>%
    gather(a, value, obs, res) %>%
    mutate(group = paste0(group, a)) %>%
    # filter(a = "obs") %>%
    ggplot()+
    geom_line(aes(conc, value, col = a, group = group))+
    # facet_wrap(~a)
    # geom_text(data = slopes, aes(x = conc-0.8, y = obsy, label =  difobs), col = "red", size = 5)+
    geom_label(aes(x = 2, y = 95, label = paste0("OF = ", round(sum(OFf(temp)$area),1))), size = 5)+
    # geom_text(data = slopes, aes(x = conc-0.2, y = obsy, label =difres), col = "blue", size = 5)+
    # geom_text(data = slopes, aes(x = conc-0.5, y = obsy, label = "-"), col = "black", size = 5)+
    # geom_text(data = slopes, aes(x = conc+0.3, y = obsy , label = paste0(" = ",difobs  - difres)), col = "black", size = 5)+
    # geom_segment(data = slopes, aes(x = conc-0.4, xend = conc -0.6, y = yboth + 2, yend = yboth +2), col = "black")+
    # geom_text(data = OFf(temp), aes(x = conc-0.5, y = reslag, label = round(area,1)))+
    # +
    geom_point(aes(conc, value, col = a))+
    geom_segment(data = temp , aes(x = conc, xend = conc, y = obs, yend = res), lty = 1)+ # ,
    # arrow = arrow(ends = "both", length = unit(0.1, "inches"))
    theme_bw()+
    scale_x_continuous(breaks = 0:4, labels = c(0:3, "Inf"))+
    labs(x = "log(conc)", y = "Cell viability (%)", col = "")+
    geom_ribbon(data = temp, aes(x = conc, ymin = obs, ymax = res, group = group), alpha = 0.2 ) #+
  # geom_text(data = temp %>% slice(nrows), aes(x = conc - 0.2, y = (obs + res) /2, label = paste0(dif, " VCs\n]",conc-1, ";", conc,"]" )))
};plot1 <- plotf(bothgroup);plot1

modif <- function(df, max = 0){


  df %>%
    group_split(group) %>%
    map(~ .x %>%
          mutate(lagobs = lag(obs), lagres=lag(res)) %>%
          mutate(difobs = lagobs - obs, difres = lagres-res)
    ) %>%
    bind_rows() %>%
    mutate(dif =difobs  -difres) %>%
    select(obs, res, conc, group, dif)-> df

  temp2 <-  df %>%
    mutate(dif2 = dif) %>%
    mutate(difprev= lag(res) - res)

  maxt = max(temp2$dif, na.rm = T)
  mint = min(temp2$dif, na.rm = T)
  torpelace <- max(c(max,min(maxt, abs(mint))))

  temp2$dif2[temp2$dif ==maxt & !is.na(temp2$dif )][[1]] <-  temp2$dif2[temp2$dif ==maxt & !is.na(temp2$dif )][[1]] - torpelace
  temp2$dif2[temp2$dif == mint& !is.na(temp2$dif )][[1]] <- temp2$dif2[temp2$dif == mint& !is.na(temp2$dif )][[1]]  + torpelace

  for(a in 2:4){

    temp2$res[[a]] <-
      temp2$res[[a-1]] - temp2$difprev[[a]]  - temp2$dif[[a]] + temp2$dif2[[a]]

  }

  temp2 %>% select(obs, res, conc, dif2, group) %>% rename(dif = dif2)
}

plot1 <- plotf(temp %>% filter(group == 1), c(2,4));plot1

plotf(temp %>% filter(group == 1), c(2,4))

# Plots one group 1

plotA <-  plotf(group1)
plotA <-  plotf(group1);plotA

first_rows <- cowplot::plot_grid(

  plotf(group1),
    # geom_rect(aes(xmin = 0.3, xmax = 0.7, ymin = 75, ymax = 95), alpha = 0, col = "black")+
    # geom_rect(aes(xmin = 1.3, xmax = 1.7, ymin = 50, ymax = 70), alpha = 0, col = "black"),
  plotf(modif(group1)),
    # geom_rect(aes(xmin = 2.3, xmax = 2.7, ymin = 25, ymax = 45), alpha = 0, col = "black")+
    # geom_rect(aes(xmin = 3.3, xmax = 3.7, ymin = 0, ymax = 20), alpha = 0, col = "black"),
  plotf(modif(group1) %>% modif) ,
    # geom_rect(aes(xmin = 1.3, xmax = 1.7, ymin = 50, ymax = 65), alpha = 0, col = "black")+
    # geom_rect(aes(xmin = 3.3, xmax = 3.7, ymin = 6, ymax = 20), alpha = 0, col = "black"),

  plotf(modif(group1) %>% modif %>% modif), nrow = 1, labels = c("A", "", "", "")
)

first_rows <- cowplot::plot_grid(

  plotf(group1) ,
    # geom_rect(aes(xmin = 0.3, xmax = 0.7, ymin = 75, ymax = 95), alpha = 0, col = "black")+
    # geom_rect(aes(xmin = 1.3, xmax = 1.7, ymin = 50, ymax = 70), alpha = 0, col = "black"),
  plotf(modif(group1)),
    # geom_rect(aes(xmin = 2.3, xmax = 2.7, ymin = 25, ymax = 45), alpha = 0, col = "black")+
    # geom_rect(aes(xmin = 3.3, xmax = 3.7, ymin = 0, ymax = 20), alpha = 0, col = "black"),
 nrow = 1, labels = c("A", "", "", "")
)

# Plots one group 2
group2 <- bothgroup %>% filter(group == 2)
cowplot::plot_grid(

  plotf(group2),
  plotf(modif(group2)),
  plotf(modif(group2) %>% modif),

  plotf(modif(group2) %>% modif %>% modif), nrow = 2
)



# Plots both groups 2
second_rows <- cowplot::plot_grid(
  plotf(bothgroup)+
    # geom_rect(aes(xmin = 0.2, xmax = 0.7, ymin = 75, ymax = 95), alpha = 0, col = "black")+
    # geom_text(aes(x = 0.3, y = 90, label = "C"))+
    # geom_rect(aes(xmin = 0.2, xmax = 0.7, ymin = 45, ymax = 70), alpha = 0, col = "black")+
    # geom_rect(aes(xmin = 1.3, xmax = 1.7, ymin = 50, ymax = 70), alpha = 0, col = "black")+
    # geom_rect(aes(xmin = 2.2, xmax = 2.7, ymin = 0, ymax = 15), alpha = 0, col = "black")+
    geom_text(aes(x = 0.3, y = 65, label = "D")),


map(0:1, ~ tibble(conc1 = 0:3, death = c(rep(F, .x), rep(T, 4- .x)), profile = 1, id = .x  )) %>%
  map(~bind_rows(.x, tibble(conc1 = 0:3, death = c(F, T, T, T), profile = 0, id = unique(.x$id)) )) %>%
  bind_rows() %>%
  mutate(fate = if_else(death, "Death", "Survival")) %>%
  # filter( id %in% 1:4) %>%
  rename(bag = id) %>%
  ggplot()+
  geom_tile(aes(factor(conc1), factor(profile), fill = fate), col = "black")+
  facet_wrap(~bag, labeller = label_both)+
  theme_bw()+
  labs(x = "log(conc)", y = "Conc drug2"),

  map(0:4, ~ tibble(conc1 = 0:3, death = c(rep(F, .x), rep(T, 4- .x)), profile = 0, id = .x  )) %>%
    map(~bind_rows(.x, tibble(conc1 = 0:3, death = c(F, T, T, T), profile = 1, id = unique(.x$id)) )) %>%
    bind_rows() %>%
    mutate(fate = if_else(death, "Death", "Survival")) %>%
    filter( id %in% 1:4) %>%
    rename(bag = id) %>%
    ggplot()+
    geom_tile(aes(factor(conc1), factor(profile), fill = fate), col = "black")+
    facet_wrap(~bag, labeller = label_both)+
  theme_bw()+
  labs(x = "log(conc)", y = "Conc drug2"), nrow = 1, labels = c("B", "C", "D")


)

cowplot::plot_grid(first_rows, second_rows, ncol = 1)

tibble(obs = c(100, 70,50,30), res = c(100, 70,49,49 - 25), conc = c(0, 1,2,3)) %>%
  mutate(dif = obs-res) -> temp2


temp2 <- modif(modif(temp))

plot2 <- plotf(temp2, c(3:4));plot2


temp3 <- tibble(obs = c(100, 70,50,30), res = c(100, 72,43 + 7 ,37), conc = c(0, 1,2,3)) %>%
  mutate(dif = obs-res)

plot3 <- plotf(temp3, c(3:4));plot3

cowplot::plot_grid(plot1, plot2, plot3)





# mutate(testcol = sample(c("a", "b","c"), 1:nrow(test),replace = T)) %>%

forppt <-makeprofile(values = c(10,7,6,2,8,7,6,2))

forppt %>%
# filter(drug == "Venetoclax" ) %>%
ggplot()+
  geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black", fill = "White")+
  facet_wrap(~drug)+
   scale_y_discrete()+ #•labels = temp$delta
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Concentration main drug (µM)", y = "Concentration A-1210477 (µM)" )

library(scales)
show_col(hue_pal()(1))


forppt %>%
  mutate(Fate = if_else(Conc1 == 1.3 & Conc2 == 10 & drug == "Venetoclax", "Death", "NA")) %>%
  mutate(Fate = if_else(Conc1 == 0.64 & Conc2 == 10 & drug == "Venetoclax", "Survival", Fate)) %>%
  filter(drug == "Venetoclax" ) %>%
  ggplot()+
  geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black")+
  facet_wrap(~drug)+
  scale_y_discrete()+ #•labels = temp$delta
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c(hue_pal()(1), "White",hue_pal()(2)[[2]]))+
  labs(x = "Concentration main drug (µM)", y = "Concentration A-1210477 (µM)" )


forppt %>%
  mutate(Fate = if_else(Conc1 == 1.3 & Conc2 == 10 & drug == "Venetoclax", "Death", "NA")) %>%
  # mutate(Fate = if_else(Conc1 == 0.64 & Conc2 == 10 & drug == "Venetoclax", "Survival", Fate)) %>%
  filter(drug == "Venetoclax" ) %>%
  ggplot()+
  geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black")+
  facet_wrap(~drug)+
  scale_y_discrete()+ #•labels = temp$delta
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c(hue_pal()(1), "White",hue_pal()(2)[[2]]))+
  labs(x = "Concentration main drug (µM)", y = "Concentration A-1210477 (µM)" )

  guides(fill = F)


  forppt %>%
    mutate(Fate = if_else(Conc1 >= 1.3 & Conc2 >= 10 & drug == "Venetoclax", "Death", "NA")) %>%
    mutate(Fate = if_else(Conc1 == 0.64 & Conc2 == 10 & drug == "Venetoclax", "Survival", Fate)) %>%
    filter(drug == "Venetoclax" ) %>%
    ggplot()+
    geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black")+
    facet_wrap(~drug)+
    scale_y_discrete()+ #•labels = temp$delta
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values = c(hue_pal()(1), "White",hue_pal()(2)[[2]]))+
    labs(x = "Concentration main drug (µM)", y = "Concentration A-1210477 (µM)" )


  forppt %>%
    mutate(Fate = if_else(Conc1 >= 1.3 & Conc2 >= 10 & drug == "Venetoclax", "Death", "NA")) %>%
    # mutate(Fate = if_else(Conc1 <= 0.64 & Conc2 <= 10 & drug == "Venetoclax", "Survival", Fate)) %>%
    filter(drug == "Venetoclax" ) %>%
    ggplot()+
    geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black")+
    facet_wrap(~drug)+
    scale_y_discrete()+ #•labels = temp$delta
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values = c(hue_pal()(1), "White",hue_pal()(2)[[2]]))+
    labs(x = "Concentration main drug (µM)", y = "Concentration A-1210477 (µM)" )



  forppt %>%
    # mutate(Fate = if_else(Conc1 >= 1.3 & Conc2 >= 10 & drug == "Venetoclax", "Death", "NA")) %>%
    # mutate(Fate = if_else(Conc1 <= 0.64 & Conc2 <= 10 & drug == "Venetoclax", "Survival", Fate)) %>%
    filter(drug == "Venetoclax" ) %>%
    ggplot()+
    geom_tile(aes(factor(Conc1), factor(Conc2), fill = Fate), col = "black")+
    facet_wrap(~drug)+
    scale_y_discrete()+ #•labels = temp$delta
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values = c(hue_pal()(1),hue_pal()(2)[[2]]))+
    labs(x = "Concentration main drug (µM)", y = "Concentration A-1210477 (µM)" )




  # S7 Gof ---------------------------------------------------------------------

  karpas <- readRDS(here('calibrated_VT','VT_both_cell_line_2.RDS'))
  sudh <-  readRDS(here('calibrated_VT','VT_both_cell_line_1.RDS'))  
  
  self <- karpas
  self2 <- sudh

  cell_line1$reconstitute() %>%
    mutate(cell = 1) %>%
    bind_rows(cell_line2$reconstitute() %>% mutate(cell = 2 )) %>%
    mutate(cell_line = "SU-DHL-4") %>%
    bind_rows(

      cell_line2$reconstitute() %>%
        mutate(cell = 1) %>%
        bind_rows(cell_line2$reconstitute() %>% mutate(cell = 2 )) %>%
        mutate(cell_line = "KARPAS-422")
    ) %>%
    mutate(conc = if_else(is.na(conc2), conc1, conc2)) %>%
    mutate(drugs = if_else(drugs == "1_4", "Venetoclax", "A-1155463")) %>%
    select(-conc1, -conc2) %>% spread(key = Reconst, value = res) %>%
    ggplot()+facet_wrap(~cell_line)+
    geom_point(aes(Value,res, col = drugs, shape = factor(conc4 )))+
    geom_abline()+
    labs(x = "Observed cell viability (%)", y = "Reconstituted cell viability (%)", col = "Drug", shape = "Concentration\nA-1210477 (µM)" )

# S10 ---------------------------------------------------------------------



  cell_line2$reconstitute() %>%
    mutate(cell = 1) %>%
    bind_rows(cell_line2$reconstitute() %>% mutate(cell = 2 )) %>%
    mutate(conc = if_else(is.na(conc2), conc1, conc2)) %>%
    mutate(drugs = if_else(drugs == "1_4", "Venetoclax", "A-1155463")) %>%
    select(-conc1, -conc2) %>% spread(key = Reconst, value = res) %>%
    ggplot()+
    geom_point(aes(Value,res, col = drugs, shape = factor(conc4 )))+
    geom_abline()+
    labs(x = "Observed cell viability (%)", y = "Reconstituted cell viability (%)", col = "Drug", shape = "Concentration\nA-1210477 (µM)" )

  
  prot_distrib <- function(self, self2, name = "cell1", name2 = "cell2"){
    tibble(cellid = self$harvest) %>% 
      left_join(self$celltheque) %>%
      mutate(name = name) %>% 
      bind_rows(
        
        tibble(cellid = self2$harvest) %>% 
          left_join(self2$celltheque) %>%
          mutate(name = name2)
        
      ) %>% 
      ungroup %>% 
      gather(Bcl20, Bclxl0, Mcl10,  BIM0, PUMA0, NOXA0,   BAK0, BAXc0, key = "key",value = "value" ) %>% 
      mutate(value = if_else(value == 0, 1, value)) %>% 
      ggplot()+
      # geom_histogram(aes(value, y = ..density.., fill = name), alpha = 0.6)+
      geom_density(aes(value, fill = name), alpha = 0.6)+
      facet_wrap(~key, scales = "free")+
      theme_bw()+
      scale_x_log10()+
      labs(x = "Protein values (nM)", fill = "Cell line")
  };prot_distrib(self, self2, "KARPAS-422", "SU-DHL4")
  