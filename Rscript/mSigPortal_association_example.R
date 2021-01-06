set_wd()
libztw()

# Module for Signature Association module 
source('Sigvisualfunc.R')





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



