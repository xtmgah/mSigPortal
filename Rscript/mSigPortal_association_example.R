set_wd()
libztw()

# Module for Signature Association module 
source('Sigvisualfunc.R')

if(Data_Source == "Public_Data"){
  # parameters for all the tables
  study_input <- "PCAWG"
  dataset_input <- "WGS"
  signature_set_name_input <- "COSMIC v3 Signatures (SBS)" 
  cancer_type_input <- "Lung-AdenoCA"
  # load exposure data files
  exposure_data_file <- paste0('../Database/Exposure/',study_input,"_",dataset_input,'_exposure_refdata.RData')
  association_data_file <- paste0('../Database/Association/',study_input,'_vardata.RData')
  
  if(!file.exists(exposure_data_file) | !file.exists(association_data_file)){
    stop("ERROR: Exposure or assoicaiton variable data are not avaiable for selected study. please check input or select other study")
  }
  load(exposure_data_file)
  load(association_data_file)
  exposure_refdata_selected <- exposure_refdata %>% filter(Signature_set_name==signature_set_name_input, Cancer_Type == cancer_type_input)
  
  exposure_refdata_selected <- exposure_refdata_selected %>% 
    select(Sample,Signature_name,Signature_exposure=Exposure) %>% 
    group_by(Sample) %>% 
    mutate(Signature_exposure_ratio=Signature_exposure/sum(Signature_exposure),Signature_exposure_cat=if_else(Signature_exposure>100, "Observed", "Not_observed")) %>% 
    ungroup() %>%
    mutate(Signature_exposure_cat=factor(Signature_exposure_cat,level=c("Not_observed","Observed")))
  
  # output jason for the variable2 
  # vardata_refdata %>%
  #   select(data_source,data_type,variable_name,variable_value_type) %>%
  #   unique() %>%
  #   toJSON(pretty = TRUE, auto_unbox = TRUE) %>%
  #   write(paste0('../Database/Association/',study_input,'_vardata.json'))
  # 
  
  vardata_refdata_selected <- vardata_refdata %>% filter(Cancer_Type == cancer_type_input)
  
  # overlapped samples
  osamples <- intersect(unique(vardata_refdata_selected$Sample),unique(exposure_refdata_selected$Sample))
  
  vardata_refdata_selected <- vardata_refdata_selected %>% filter(Sample %in% osamples)
  exposure_refdata_selected <- exposure_refdata_selected %>% filter(Sample %in% osamples)
  
  # expsorue variant list
  Exposure_varlist <- colnames(exposure_refdata_selected)[-c(1:2)]
  #Association_varlist, # load corresponding json file as the paramters for variant  
  Exposure_varinput <- "Signature_exposure_ratio"
  Association_varinput_source <- 'genomic data'
  Association_varinput_type <- 'evolution_and_heterogeneity'
  Association_varinput_name <- 'purity'
  
  exposure_refdata_selected <- exposure_refdata_selected %>% select(Sample,Exposure_varinput)
  vardata_refdata_selected <- vardata_refdata_selected %>%
    filter(data_source==Association_varinput_source, data_type==Association_varinput_type,variable_name==Association_varinput_name)  
    
  if(unique(vardata_refdata_selected$variable_value_type) == "numeric") { vardata_refdata_selected$variable_value <- as.numeric(vardata_refdata_selected$variable_value )} 
    
  vardata_refdata_selected <- vardata_refdata_selected %>% 
    pivot_wider(id_cols = Sample,names_from = variable_name,values_from = variable_value)
  
  data_input <- left_join(vardata_refdata_selected,exposure_refdata_selected)
  
  mSigPortal_associaiton(data=data_input,Var1 = Association_varinput_name, Var2=Exposure_varinput,type = "parametric",xlab=Association_varinput_name, ylab=Exposure_varinput,filter_zero1=FALSE, filter_zero2=FALSE,log1=FALSE,log2=TRUE, type="parametric", collapse_var1=NULL, collapse_var2=NULL, file = "association_result.svg")
  
  ## asssociation_data.txt will output as download text file. 
  data_input %>% write_delim(file = 'asssociation_data.txt',delim = '\t',col_names = T,na = '')
}

