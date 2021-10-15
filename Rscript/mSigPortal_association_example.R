set_wd()
libztw()

# Module for Signature Association module 
source('Sigvisualfunc.R')

if(Data_Source == "Public_Data"){
  # parameters for all the tables
  study_input <- "PCAWG"
  dataset_input <- "WGS"
  signature_set_name_input <- "COSMIC_v3_Signatures_GRCh37_SBS96" 
  cancer_type_input <- "Panc-AdenoCA"
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
  #Signature_name_input <- 'SBS4'
  Exposure_varinput <- "Signature_exposure_ratio"
  exposure_refdata_selected <- exposure_refdata_selected%>% select(Sample,Signature_name,Exposure_varinput)
  
  
  ## extract the variable information
  clist <- vardata_refdata_selected %>% select(data_source,data_type,variable_name,variable_value_type) %>% unique()
  vardata_refdata_selected %>% write_delim('vardata_refdata_selected.txt',delim = '\t',col_names = T)
  # clist will be used for the Assocaition Variable Data and Select Variables. 
  
  if(regression == TRUE) {
    ## the threshold and collapse by default change to NA because the NULL is not working here
    # Var1 <- c('quality control', 'sequencing metrics', 'FWHM_Normal',NA, FALSE, NA)
    # Var2 <- c('genomic data', 'sv count', 'DEL',NA, TRUE, NA)
    # Var3 <- c('germline data', 'ancestry', 'EUR',NA, FALSE, NA)
    Var1 <- list(source = 'quality control', type = 'sequencing metrics', name = 'FWHM_Normal', filter = NULL, log2 = FALSE, collapse = NULL)
    Var2 <- list(source = 'genomic data', type = 'sv count', name = 'DEL', filter = NULL, log2 = TRUE, collapse = NULL)
    #Var3 <- list(source = 'germline data', type = 'ancestry', name = 'EUR', filter = NULL, log2 = FALSE, collapse = NULL)
    Var3 <- list(source = 'quality control', type = 'sequencing metrics', name = '%_of_paired_reads_mapping_to_different_chromosomes_Normal', filter = NULL, log2 = FALSE, collapse = NULL)
    
    ### add more parameters according to user's input
    listpars <- list(Var1, Var2, Var3)
    vardata_refdata_selected <- multivariable_inputs(vardata_refdata_selected, listpars)
    data_input <- left_join(exposure_refdata_selected,vardata_refdata_selected) %>% select(-Sample)
    
    
    ## change variable name if detected special chacters ## 
    colnames(data_input)[-c(1:2)] <-  str_replace_all(str_replace(str_replace_all(colnames(data_input)[-c(1:2)],"[^[:alnum:]_ ]*", ""),"^[^[:alpha:]]*",""),"  *","_")
    
    rformula = paste0(Exposure_varinput, " ~ ",paste0(colnames(data_input)[-c(1:2)],collapse = ' + '))
    ## regressionby group of signature name
    result <- mSigPortal_associaiton_group(data=data_input,Group_Var = "Signature_name",type = "glm",regression = TRUE,formula = rformula )
    ## put result as a short table above the figure
    
    signature_name_list <- unique(result[[1]]) ## drop-down list for the signature name
    signature_name_input <- signature_name_list[1]  ## by default, select the first signature name
    data_input <- data_input %>% filter(Signature_name == signature_name_input) %>% select(-Signature_name)

    mSigPortal_associaiton(data=data_input,type = "glm",regression = TRUE,formula = rformula, output_plot = "association_result.svg")
    data_input %>% write_delim(file = 'asssociation_data.txt',delim = '\t',col_names = T,na = '')
    
    ### will add codes for the regression analysis / multivariate analysis 
    
    # input_formula <- 'lm(Signature_exposure ~ '
    # 
    # 
    # exposure_refdata_selected
    # 
    # vardata_refdata_selected
    # 
    # Group_Var = "Signature_name",
    # Var1 = args$associationVar$name, Var2 = args$exposureVar$name, type = args$associationVar$type,
    # filter1 = args$associationVar$filter, filter2 = args$exposureVar$filter,
    # log1 = args$associationVar$log2, log2 = args$exposureVar$log2,
    # collapse_var1 = args$associationVar$collapse, collapse_var2 = NULL
    
    
  }else
  {
    Association_varinput_source <- 'genomic data'
    Association_varinput_type <- 'consensus cnv features'
    Association_varinput_name <- 'wgd_status'
    # Association_varinput_source <- 'genomic data'
    # Association_varinput_type <- 'evolution_and_heterogeneity'
    # Association_varinput_name <- 'purity'
    vardata_refdata_selected <- vardata_refdata_selected %>%
      filter(data_source==Association_varinput_source, data_type==Association_varinput_type,variable_name==Association_varinput_name)  
    
    if(unique(vardata_refdata_selected$variable_value_type) == "numeric") { vardata_refdata_selected$variable_value <- as.numeric(vardata_refdata_selected$variable_value )} 
    
    vardata_refdata_selected <- vardata_refdata_selected %>% 
      pivot_wider(id_cols = Sample,names_from = variable_name,values_from = variable_value)
    
    
    ## check data integration
    vardata_refdata_selected <- exposure_refdata_selected %>% select(Sample) %>% unique() %>% left_join(vardata_refdata_selected)
    ## including NA
    if(length(unique(vardata_refdata_selected[[2]])) == 1 ) {
      stop(paste0("mSigPortal Association failed: the selected variable name ",Association_varinput_name," have only unique value: ", unique(vardata_refdata_selected[[2]]),'.'))
    }
    tmpdata <- vardata_refdata_selected
    colnames(tmpdata)[2] <- 'Variable'
    tmpvalue <- tmpdata %>% count(Variable) %>% filter(n<2) %>% dim() %>% .[[1]]
    
    if(tmpvalue !=0){
      stop(paste0("mSigPortal Association failed: the selected variable name ",Association_varinput_name," have not enough obsevations for both levels."))
    }
    
    ### combined dataset
    data_input <- left_join(vardata_refdata_selected,exposure_refdata_selected) %>% select(-Sample)
    
    ## dropdown list for collapse_var1 and collapse_var2 
    collapse_var1_list <- levels(data_input[[Association_varinput_name]])
    collapse_var2_list <- levels(data_input[[Exposure_varinput]])
    
    ## association test by group of signature name
    result <- mSigPortal_associaiton_group(data=data_input,Group_Var = "Signature_name",Var1 = Association_varinput_name, Var2=Exposure_varinput,type = "parametric",filter1=NULL, filter2=NULL,log1=FALSE,log2=FALSE,  collapse_var1=NULL, collapse_var2=NULL)
    ## put result as a short table above the figure
    
    signature_name_list <- unique(result[[1]]) ## dropdown list for the signature name
    signature_name_input <- signature_name_list[1]  ## by default, select the first signature name
    
    data_input <- data_input %>% filter(Signature_name == signature_name_input) %>% select(-Signature_name)
    
    mSigPortal_associaiton(data=data_input,Var1 = Association_varinput_name, Var2=Exposure_varinput,type = "parametric",xlab=Association_varinput_name, ylab=Exposure_varinput,filter1=NULL, filter2=NULL,log1=FALSE,log2=FALSE, collapse_var1=NULL, collapse_var2=NULL, output_plot = "association_result.svg")
    
    ## asssociation_data.txt will output as download text file. 
    data_input %>% write_delim(file = 'asssociation_data.txt',delim = '\t',col_names = T,na = '')
  }
}

