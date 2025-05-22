#Run Rarefaction/Extrapolation with iNEXT3D


####packages####
if(!require("neonDivData"))
  install.packages('neonDivData', repos = c(
    daijiang = 'https://daijiang.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

library(neonDivData)
library(tidyverse)
library(iNEXT.3D)
library(iNEXT)

####Beetles####
#make a sampling bout per plot by species matrix with count data
beetle_df <- data_beetle %>%
              mutate(plot_date = paste(plotID, observation_datetime)) %>%  # Create a unique plot_date identifie
              select(plot_date, plotID, taxon_name) %>%
              group_by(plot_date, plotID, taxon_name) %>%# 
              summarise(count = n(), .groups = 'drop') %>%# Summarize to get count data for each species per sampling bout per plot (i.e. plot_date)
              pivot_wider(names_from = taxon_name, values_from = count, values_fill = list(count = 0)) %>% # Pivot wider to create a matrix of taxa counts fill in zeros
              column_to_rownames(var = "plot_date")#make plot_date rownames

#make a  plot sampling bout per plot by species incidence (1/0) matrix
beetle_df2 <- data_beetle %>% 
              mutate(plot_date = paste(plotID, observation_datetime))%>%
              select(plot_date, plotID, taxon_name) %>% 
              mutate(present = 1) %>% 
              group_by(plot_date,plotID, taxon_name) %>% 
              summarise(present = sum(present)/sum(present)) %>% 
              pivot_wider(names_from = taxon_name, values_from = present, values_fill = 0)%>%  
              #ungroup()%>%
              #mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% # Calculate total number of taxa present in each row
              #filter(Total>1) %>%# Filter rows with more than one taxon present
              #select (!Total)%>%
              #group_by(plotID) %>% 
              #add_tally() %>%
              #filter(n>10) %>% #Filter plotIDs with more than 5 sampling bouts (i.e., rows)
              #select (!n)%>% 
              column_to_rownames(var = "plot_date")


##make a plot by species count matrix
beetle_df3 <- data_beetle %>% 
              select(plotID, taxon_name) %>% 
              mutate(present = 1) %>% 
              group_by(plotID, taxon_name) %>% 
              summarise(present = sum(present)) %>% 
              pivot_wider(names_from = taxon_name, values_from = present, values_fill = 0)%>%  
              ungroup()

#make a plot by species incidence (1/0) matrix
beetle_df4 <- data_beetle %>%
              select( plotID, taxon_name) %>%
              mutate(present = 1) %>% 
              group_by(plotID, taxon_name) %>% 
              summarise(present = sum(present)/sum(present)) %>% 
              pivot_wider(names_from = taxon_name, values_from = present, values_fill = 0)%>%  
              

####iNEXT3D for rarefaction
  
#make a  plot sampling bout per plot by species incidence (1/0) matrix with filtering out plots with few species 
beetle_df_rar <- data_beetle %>% 
                mutate(plot_date = paste(plotID, observation_datetime))%>%
                select(plot_date, plotID, taxon_name) %>% 
                mutate(present = 1) %>% 
                group_by(plot_date,plotID, taxon_name) %>% 
                summarise(present = sum(present)/sum(present)) %>% 
                pivot_wider(names_from = taxon_name, values_from = present, values_fill = 0)%>%  
                ungroup()%>%
                mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% # Calculate total number of taxa present in each row
                filter(Total>1) %>%# Filter rows with more than one taxon present
                select (!Total)%>%
                group_by(plotID) %>% 
                add_tally() %>%
                filter(n>10) %>% #Filter plotIDs with more than 10 sampling bouts (i.e., rows)
                select (!n)%>% 
                column_to_rownames(var = "plot_date")
  
  
  
#format as list of lists to run in iNEXT3D
lmz<-split(beetle_df_rar,beetle_df2$plotID)%>%
          lapply(., function(x)x[,-1 ])%>%
          lapply(.,t)




#run iNEXT3D Takes a while for whole data set
#this throws an error when running the entire data set because of plots with small numbers.
ass<-iNEXT3D(lmz, diversity = 'TD', q = c(0),datatype ="incidence_raw")


#run on a subset
zz<-lmz[c(1:10)]

ass<-iNEXT3D(zz, diversity = 'TD', q = c(0),datatype ="incidence_raw")



#plot curves. can;t see much with so many plots        

ggiNEXT3D(ass, type = 1, facet.var = "Assemblage")

#see if observed species richness estimate falls within 95% CIs of asymptotic estimator fit= (True/False) and "percent" of estimated
beetle_plot_rar<-ass$TDAsyEst%>%
                rename(plotID=Assemblage)%>%
                arrange(plotID)%>%
                group_by(plotID)%>%
                mutate(fit=TD_obs>=qTD.LCL & TD_obs<=qTD.UCL)%>%#compare the observed to the lcl and ucl of the extrapolation
                mutate(percent= TD_obs/TD_asy)%>%
                #filter(Method=="Observed")%>%
                select(plotID,qTD,fit,percent)

