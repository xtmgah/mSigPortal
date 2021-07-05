set_wd()
libztw()

# Module for Signature Association module 
source('Sigvisualfunc.R')


# Public dataset sample Size ----------------------------------------------
load('../Database/Exposure/PCAWG_WGS_exposure_refdata.RData')
data1 <- exposure_refdata %>% mutate(Key=paste0(Study,"_",Dataset)) %>% select(Key,Cancer_Type,Sample) %>% unique() %>% count(Key,Cancer_Type) 
load('../Database/Exposure/Breast560_WGS_exposure_refdata.RData')
data2 <- exposure_refdata %>% mutate(Key=paste0(Study,"_",Dataset)) %>% select(Key,Cancer_Type,Sample) %>% unique() %>% count(Key,Cancer_Type)
load('../Database/Exposure/Breast80_WGS_exposure_refdata.RData')
data3 <- exposure_refdata %>% mutate(Key=paste0(Study,"_",Dataset)) %>% select(Key,Cancer_Type,Sample) %>% unique() %>% count(Key,Cancer_Type)
data3.1<- tibble(Key="Sherlock-Lung-WGS",Cancer_Type="LCINS",n=232) 
load('../Database/Exposure/non-PCAWG_WGS_exposure_refdata.RData')
data3.2 <- exposure_refdata %>% mutate(Key=paste0(Study,"_",Dataset)) %>% select(Key,Cancer_Type,Sample) %>% unique() %>% count(Key,Cancer_Type) 
datall1 <- bind_rows(data1,data2,data3,data3.1,data3.2)
tmp <- datall1 %>% group_by(Cancer_Type) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% pull(Cancer_Type)
pall1 <- datall1%>% 
  mutate(Cancer_Type=factor(Cancer_Type,levels = tmp)) %>% 
  ggplot(aes(Cancer_Type,n,fill=Key))+geom_col()+theme_ipsum_rc(axis = "XY",axis_title_size = 16,axis_title_just = 'm',grid="Y",strip_text_face = 'bold',ticks = TRUE)+labs(x="",y="Number of samples\n",fill="WGS")+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+scale_y_continuous(expand = c(0,0),limits = c(0,500))+scale_x_discrete(expand = c(0,0))+geom_text(aes(label = n), vjust = -0.5)+scale_fill_d3()



load('../Database/Exposure/TCGA_WES_exposure_refdata.RData')
data4 <- exposure_refdata %>% mutate(Key=paste0(Study,"_",Dataset)) %>% select(Key,Cancer_Type,Sample) %>% unique() %>% count(Key,Cancer_Type) 
load('../Database/Exposure/non-PCAWG_WES_exposure_refdata.RData')
data5 <- exposure_refdata %>% mutate(Key=paste0(Study,"_",Dataset)) %>% select(Key,Cancer_Type,Sample) %>% unique() %>% count(Key,Cancer_Type) 

datall2 <- bind_rows(data4,data5)
tmp <- datall2 %>% group_by(Cancer_Type) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% pull(Cancer_Type)
pall2 <- datall2%>% 
  mutate(Cancer_Type=factor(Cancer_Type,levels = tmp)) %>% 
  ggplot(aes(Cancer_Type,n,fill=Key))+geom_col()+theme_ipsum_rc(axis = "XY",axis_title_size = 16,axis_title_just = 'm',grid="Y",strip_text_face = 'bold',ticks = TRUE)+labs(x="",y="Number of samples\n",fill="WES")+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+scale_y_continuous(expand = c(0,0),limits = c(0,1200))+scale_x_discrete(expand = c(0,0))+geom_text(aes(label = n), vjust = -0.5)+scale_fill_d3()



ggsave(filename = 'pall1.svg',plot = pall1,width = 20,height = 8,device = cairo_pdf)
ggsave(filename = 'pall2.svg',plot = pall2,width = 24,height = 8,device = cairo_pdf)





# dataset -----------------------------------------------------------------
load('../TMP/Association_Prepare/PCAWG_vardata.RData')







# Exposure Data -----------------------------------------------------------
if(Data_Source == "Public_Data"){
  # parameters for all the tables
  study_input <- "Breast560"
  dataset_input <- "WGS"
  signature_set_name_input <- "Organ-specific Cancer Signatures (SBS)" 
  
  # load exposure data files
  exposure_data_file <- paste0('../Database/Exposure/',study_input,"_",dataset_input,'_exposure_refdata.RData')
  
  if(!file.exists(exposure_data_file)){
    stop("ERROR: Exposure data is not avaiable for selected study. please select other study")
  }
  
  load(exposure_data_file)
  
  exposure_refdata_selected <- exposure_refdata %>% filter(Signature_set_name==signature_set_name_input)
  
  genome <- case_when(
    study_input == "PCAWG" ~ "GRCh37",
    study_input == "TCGA" ~ "GRCh37",
    TRUE ~ "GRCh37"
  )
  genomesize <-  genome2size(genome)
  
  cancer_types <- exposure_refdata_selected %>% pull(Cancer_Type) %>% unique()
  
}


#  
# exposure data as Number of mutations



# exposure data as Proportion

# exposure data as Presence


library(inspectdf)
  exposure_refdata %>% inspect_num() %>% show_plot()
exposure_refdata %>% ggplot(aes(log2(Exposure+1),col=Signature_name))+geom_density()+facet_wrap(~Cancer_Type )+scale_color_d3()



# Variable Data -----------------------------------------------------------



