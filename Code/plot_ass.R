#Run Rarefaction/Extrapolation with iNEXT3D


####packages####
if(!require("neonDivData"))
  install.packages('neonDivData', repos = c(
    daijiang = 'https://daijiang.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))
library(neonDivData)
library(tidyverse)
library(iNEXT.3D)
install.packages("iNEXT")
library(iNEXT)
####Beetles####
#plot-scale rarefaction

#filter for plots with at least one sample bout
good_beetle_plots<-data_beetle %>%
                    group_by(plotID) %>%
                    summarise(n_observation = n_distinct(observation_datetime))%>%
                    filter(!((n_observation<=3)))%>%#remove plots with only 1 sampling bout (can't rarefy)
                    select(plotID)
bad_plots<-c("WREF_001","BARR_036","BARR_049","BARR_083","BARR_084","RMNP_012","TEAK_010","TEAK_004")

#make a species by sampling bout incidence (1/0) matrix
beetle_df2 <- data_beetle %>% 
            mutate(plot_date = paste(plotID, observation_datetime))%>%
            select(plot_date, plotID, taxon_name) %>% 
            mutate(present = 1) %>% 
            group_by(plot_date,plotID, taxon_name) %>% 
            summarise(present = sum(present)/sum(present)) %>% 
            pivot_wider(names_from = taxon_name, values_from = present, values_fill = 0)%>%  
            ungroup()%>%
            mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
            filter(Total>1) %>%
            select (!Total)%>%
            group_by(plotID) %>% 
            add_tally() %>%
            filter(n>5) %>% 
            select (!n)%>% 
            column_to_rownames(var = "plot_date")%>% 
            filter(!plotID%in%bad_plots)#take only "good" plots


            group_by(plotID) %>% 
            select(!plot_date)%>%   
            mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
            filter(Total==2) %>%
            select (!Total)%>%
  
  add_tally() %>%
  filter(n>5) %>% 
  select (!n)
            group_by(plotID) %>% 
            add_tally() %>%
            filter(n>50) %>% 
            select (!n)

##abundance based
beetle_df2 <- data_beetle %>% 
              filter(plotID%in%good_beetle_plots$plotID)%>%#take only "good" plots
              #mutate(plot_date = paste(plotID, observation_datetime))%>%
              select( plotID, taxon_name) %>% 
              mutate(present = 1) %>% 
              group_by(plotID, taxon_name) %>% 
              summarise(present = sum(present)) %>% 
              pivot_wider(names_from = taxon_name, values_from = present, values_fill = 0)%>%  
              ungroup()

2

#format as list of lists to run in iNEXT3D
lmz<-split(beetle_df2,beetle_df2$plotID)%>%
    lapply(., function(x)x[,-1 ])%>%
    lapply(.,t) 

as.integer(lmz)
lmz[15]
zx<-Fish_incidence_data[1]


#issue
#c(WREF_001,BARR_036,BARR_049,BARR_083,BARR_084,RMNP_012,TEAK_010)

#run iNEXT3D
ass<-ObsAsy3D(zz, diversity = 'TD', q = c(0),
              datatype = "incidence_raw")



ass<-iNEXT3D(lmz, diversity = 'TD', q = c(0),datatype ="incidence_raw")


#ass<-iNEXT3D(lmz, diversity = 'TD', q = c(0),datatype ="incidence_raw")

ass<-iNEXT(lmz, q = c(0),datatype ="incidence_freq")

 colSums(lmz[259]$RMNP_012)             

CC<-ass$TDAsyEst

ggiNEXT3D(ass, type = 1, facet.var = "Assemblage")

#see if observed species richness estimate falls within 95% CIs of asymptotic estimator (True/False)
beetle_plot_rar<-ass%>%
            rename(plotID=Assemblage)%>%
            arrange(plotID)%>%
            group_by(plotID)%>%
            mutate(fit= qTD>=qTD.LCL[row_number()-1] & qTD<=qTD.UCL[row_number()-1])%>%#compare the observed to the lcl and ucl of the extrapolation
            mutate(percent= qTD/qTD[row_number()-1])%>%
            filter(Method=="Observed")%>%
            select(plotID,qTD,fit,percent)

#PLOT
ggObsAsy3D(ass$Assemblage, profile = "q", group)



write_rds(ass, "Data/beetle_plot_rar.csv")


####Birds####

good_bird<-read.csv("Data/good_bird_plots_4.csv")

#make a species by sampling bout incidence (1/0) matrix
bird_df <-  data_bird %>% 
              filter(plotID%in%good_bird$plotID)%>%#take only "good" plots
              mutate(date= as.Date(observation_datetime))%>%
              mutate(plot_date = paste(plotID, date))%>%
              select(plot_date, plotID, taxon_name) %>% 
              mutate(present = 1) %>% 
              group_by(plot_date,plotID, taxon_name) %>% 
              summarise(present = sum(present)/sum(present)) %>% 
              pivot_wider(names_from = taxon_name, values_from = present, values_fill = 0)%>%  
              ungroup()%>%
              select(!plot_date)







#format as list of lists to run in iNEXT3D
lm<-split(bird_df,bird_df$plotID)%>%
    lapply(., function(x)x[,-1 ])%>%
    lapply(.,t)  


#run iNEXT3D
ass<-ObsAsy3D(lm, diversity = 'TD', q = c(0),
              datatype = "incidence_raw")



##see if observed species richness estimate falls within 95% CIs of asymptotic estimator (True/False)
bird_plot_rar<-ass%>%
            rename(plotID=Assemblage)%>%
            arrange(plotID)%>%
            group_by(plotID)%>%
            mutate(fit= qTD>=qTD.LCL[row_number()-1] & qTD<=qTD.UCL[row_number()-1])%>%#compare the observed to the lcl and ucl of the extrapolation
            mutate(percent= qTD/qTD[row_number()-1])%>%
            filter(Method=="Observed")%>%
            select(plotID,qTD,fit,percent)

#PLOT
ggObsAsy3D(ass, profile = "q")



write.csv(bird_plot_rar, "Data/bird_plot_rar.csv")
####Plants####

good_plant<-read.csv("Data/good_plant_plots_4.csv")

#make a species by sampling bout incidence (1/0) matrix
plant_df <-  data_plant %>% 
            filter(plotID%in%good_plant$plotID)%>%#take only "good" plots
            mutate(date= as.Date(observation_datetime))%>%
            mutate(plot_date = paste(plotID, date))%>%
            select(plot_date, plotID, taxon_name) %>% 
            mutate(present = 1) %>% 
            group_by(plot_date,plotID, taxon_name) %>% 
            summarise(present = sum(present)/sum(present)) %>% 
            pivot_wider(names_from = taxon_name, values_from = present, values_fill = 0)%>%  
            ungroup()%>%
            select(!plot_date)




#format as list of lists to run in iNEXT3D
lm<-split(plant_df,plant_df$plotID)%>%
    lapply(., function(x)x[,-1 ])%>%
    lapply(.,t)  


#run iNEXT3D
ass<-ObsAsy3D(lm, diversity = 'TD', q = c(0),
              datatype = "incidence_raw")


##see if observed species richness estimate falls within 95% CIs of asymptotic estimator (True/False)
plant_plot_rar<-ass%>%
              rename(plotID=Assemblage)%>%
              arrange(plotID)%>%
              group_by(plotID)%>%
              mutate(fit= qTD>=qTD.LCL[row_number()-1] & qTD<=qTD.UCL[row_number()-1])%>%#compare the observed to the lcl and ucl of the extrapolation
              mutate(percent= qTD/qTD[row_number()-1])%>%
              filter(Method=="Observed")%>%
              select(plotID,qTD,fit,percent)



#PLOT
ggObsAsy3D(ass, profile = "q")



write.csv(plant_plot_rar, "Data/plant_plot_rar.csv")


####mammals####

good_mammal<-read.csv("Data/good_mammal_plots_12.csv")

#make a species by sampling bout incidence (1/0) matrix
mammal_df <-  data_small_mammal %>% 
              filter(plotID%in%good_mammal$plotID)%>%#take only "good" plots
              mutate(plot_date = paste(plotID, observation_datetime))%>%
              select(plot_date, plotID, taxon_name) %>% 
              mutate(present = 1) %>% 
              group_by(plot_date,plotID, taxon_name) %>% 
              summarise(present = sum(present)/sum(present)) %>% 
              pivot_wider(names_from = taxon_name, values_from = present, values_fill = 0)%>%  
              ungroup()%>%
              select(!plot_date)



#format as list of lists to run in iNEXT3D
lm<-split(mammal_df,mammal_df$plotID)%>%
    lapply(., function(x)x[,-1 ])%>%
    lapply(.,t)  

#run iNEXT3D

ass<-ObsAsy3D(lm, diversity = 'TD', q = c(0),
              datatype = "incidence_raw")



##see if observed species richness estimate falls within 95% CIs of asymptotic estimator (True/False)

mammal_plot_rar<-ass%>%
                rename(plotID=Assemblage)%>%
                arrange(plotID)%>%
                group_by(plotID)%>%
                mutate(fit= qTD>=qTD.LCL[row_number()-1] & qTD<=qTD.UCL[row_number()-1])%>%#compare the observed to the lcl and ucl of the extrapolation
                mutate(percent= qTD/qTD[row_number()-1])%>%
                filter(Method=="Observed")%>%
                select(plotID,qTD,fit,percent)

#PLOT
ggObsAsy3D(ass, profile = "q")



write.csv(mammal_plot_rar, "Data/mammal_plot_rar.csv")
