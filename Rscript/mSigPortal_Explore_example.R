set_wd()
library(tidyverse)
library(cowplot)
library(hrbrthemes)
library(ggsci)
library(ggrepel)
library(ggforce)
library(ggtext)

# source msigportal function ----------------------------------------------
source('Sigvisualfunc.R')

# load reference signatures files -----------------------------------------
load('../Database/Signature/signature_refsets.RData')


# Signature Explore -------------------------------------------------------

# section 1: Current reference signatures in mSigPortal -------------------
nsig_data <-  signature_refsets %>% 
  group_by(Profile,Signature_set_name) %>% 
  summarise(N=n_distinct(Signature_name)) %>% 
  ungroup() %>% 
  mutate(Profile=factor(Profile,levels = c("SBS96","SBS192","SBS1536","ID83","DBS78","RS32")))

nsig_data <- nsig_data %>% left_join(
  nsig_data %>% group_by(Profile) %>% summarise(Total=sum(N)) 
) %>% 
  mutate(freq=N/Total) %>% 
  mutate(N2=if_else(freq>0.02,as.character(N),""))

# put the follow pie-chart on website 
signature_piechart(nsig_data,sigsetcolor,output_plot = 'tmp.svg')



# section 2: Mutational signature profile  --------------------------------------------------------------
# # add signature SVG plots for all the signatures and put on line
# sigPlt.plotDBS('COSMIC v3 Signatures (DBS).txt','COSMIC_v3_Signatures_DBS','mSigPortal','78',percentage=True)
# sigPlt.plotSBS('COSMIC v3 Signatures (SBS).txt','COSMIC_v3_Signatures_SBS','mSigPortal','96',percentage=True)
# sigPlt.plotSBS('COSMIC v2 Signatures (SBS).txt','COSMIC_v2_Signatures_SBS','mSigPortal','96',percentage=True)
# sigPlt.plotSBS('Organ-specific Cancer Signatures (SBS).txt','Organ-specific_Cancer_Signatures_SBS','mSigPortal','96',percentage=True)
# sigPlt.plotID('SignatureAnalyzer PCAWG WGS Signatures (ID).txt','SignatureAnalyzer_PCAWG_WGS_Signatures_ID','mSigPortal','83',percentage=True)
# sigPlt.plotSBS('Environmental Mutagen Signatures (SBS).txt','Environmental_Mutagen_Signatures_SBS','mSigPortal','96',percentage=True)
# #sigPlt.plotRS('Cancer Reference Signatures (RS).txt','Cancer_Reference_Signatures_RS','mSigPortal','32',percentage=True)
# sigPlt.plotSBS('Cancer Reference Signatures (SBS).txt','Cancer_Reference_Signatures_SBS','mSigPortal','96',percentage=True)
# sigPlt.plotSBS('SigProfiler PCAWG Strand Signatures (SBS).txt','SigProfiler_PCAWG_Strand_Signatures_SBS','mSigPortal','192',percentage=True)
# #sigPlt.plotRS('Organ-specific Cancer Signatures (RS).txt','Organ-specific Cancer Signatures (RS).txt','mSigPortal','32',percentage=True)
# sigPlt.plotDBS('SignatureAnalyzer PCAWG WGS Signatures (DBS).txt','SignatureAnalyzer_PCAWG_WGS_Signatures_DBS','mSigPortal','78',percentage=True)
# sigPlt.plotID('COSMIC v3 Signatures (ID).txt','COSMIC_v3_Signatures_ID','mSigPortal','83',percentage=True)
# sigPlt.plotSBS('SignatureAnalyzer PCAWG WGS Signatures (SBS).txt','SignatureAnalyzer_PCAWG_WGS_Signatures_SBS','mSigPortal','96',percentage=True)
# sigPlt.plotSBS('SigProfiler PCAWG WXS Signatures (SBS).txt','SigProfiler_PCAWG_WXS_Signatures_SBS','mSigPortal','96',percentage=True)
# sigPlt.plotSBS('SignatureAnalyzer PCAWG WGS 1536 Signatures (SBS).txt','SignatureAnalyzer_PCAWG_WGS_1536_Signatures_SBS','mSigPortal','1536',percentage=True)
# sigPlt.plotSBS('Other published signatures (SBS).txt','Other_published_signatures_SBS','mSigPortal','96',percentage=True)
# sigPlt.plotID('Other published signatures (ID).txt','Other_published_signatures_ID','mSigPortal','83',percentage=True)

path_profile <- '../Database/Signature/Reference_Signature_Profiles_SVG/'
signature_profile_files <- signature_refsets %>% select(Source,Profile,Signature_set_name,Dataset,Signature_name) %>% unique() %>% mutate(Path=str_replace_all(Signature_set_name," ","_"),Path=str_remove_all(Path,"[()]"),Path=paste0(path_profile,'Profiles/',Path,"/",Signature_name,".svg"))

## multiple selected profiles
signature_source_input <- "Reference_signatures"
profile_name_input <- "SBS96"
signatureset_name_input <- "COSMIC v3 Signatures (SBS)"
dataset_input <- "WGS"
signature_name_input <- "SBS1"
svgfile_selected <- signature_profile_files %>% filter(Source==signature_source_input,Profile==profile_name_input,Signature_set_name==signatureset_name_input,Dataset==dataset_input,Signature_name==signature_name_input) %>% pull(Path)

# link to signatures
# signature_refsets %>% select(Signature_set_name,Signature_name) %>% unique() %>%
#   mutate(Link=case_when(
#     Signature_set_name == "COSMIC v3 Signatures (SBS)" ~  paste0('https://cancer.sanger.ac.uk/cosmic/signatures/SBS/',Signature_name,'.tt'),
#   )
#   ) %>% write_csv('../Database/Signature/signature_link.csv',col_names = T,na = '')



# section3: Cosine similarities among mutational signatures -------------------------
# The parameters will be “Matrix Size”, “Reference Signature Set1” and “Reference Signature Set2”. 
profile_name_input <- "SBS96" # profile type
# the availabe options for signaturesetname1 and signaturesetname2 will be:
signature_refsets %>% filter(Profile==profile_name_input)%>% pull(Signature_set_name) %>% unique()

signatureset_name1 <- "Environmental Mutagen Signatures (SBS)" 
signatureset_name2 <- "COSMIC v3 Signatures (SBS)"


sigrefset1_data <- signature_refsets %>%
  filter(Profile==profile_name_input,Signature_set_name==signatureset_name1) %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)

sigrefset2_data <- signature_refsets %>%
  filter(Profile==profile_name_input,Signature_set_name==signatureset_name2) %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)

cos_sim_res=cos_sim_df(sigrefset1_data,sigrefset2_data) 

# put this heatmap on the web
plot_cosine_heatmap_df(cos_sim_res,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
# add a link bellow the heatmap
cos_sim_res %>% write_delim('signature_cos_sim_res.txt',delim = '\t',col_names = T)



# section4: Mutational signatures comparisons -----------------------------
## A comparison of two reference signatures
# There will be five parameters: “Profile Type”,  “Reference Signature Set1”, “Signature Name1”, “Reference Signature Set2”, “Signature Name2”; 
profile_name_input <- "SBS96" # profile type
signatureset_name1 <- "COSMIC v3 Signatures (SBS)"
signatureset_name2 <- "COSMIC v3 Signatures (SBS)"
signature_name1 <-  "SBS1"
signature_name2 <- "SBS5"

# similarly, the available options for signatureset_name1 and signatureset_name2 will be:
signature_refsets %>% filter(Profile==profile_name_input)%>% pull(Signature_set_name) %>% unique()
# the available options for the signature names will be:
signature_refsets %>% filter(Profile==profile_name_input,Signature_set_name==signatureset_name1)%>% pull(Signature_name) %>% unique()
signature_refsets %>% filter(Profile==profile_name_input,Signature_set_name==signatureset_name2)%>% pull(Signature_name) %>% unique()


profile1 <- signature_refsets %>%
  filter(Profile==profile_name_input,Signature_set_name==signatureset_name1) %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution) %>% 
  select(MutationType,one_of(signature_name1))

profile2 <- signature_refsets %>%
  filter(Profile==profile_name_input,Signature_set_name==signatureset_name2) %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution) %>% 
  select(MutationType,one_of(signature_name2))

# put this plot on the web:
plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,output_plot = 'tmp.svg')







# Exposure Exploring ------------------------------------------------------


if(Data_Source == "Public_Data"){
  
  # load exposure data files
  load('../Database/Exposure/exposure_refdata.RData')
  load('../Database/Signature/signature_refsets.RData')
  load('../Database/Seqmatrix/seqmatrix_refdata.RData')
  
  # parameters for all the tables
  study_input <- "PCAWG"
  dataset_input <- "WGS"
  signature_set_name_input <- "COSMIC v3 Signatures (SBS)" 
  exposure_refdata_selected <- exposure_refdata %>% filter(Study==study_input,Dataset==dataset_input,Signature_set_name==signature_set_name_input)
  
  genomesize <- 3101.976562
  
  ## reduce the data for signature and profile
  #signature_refsets %>% select(Profile,Signature_set_name) %>% unique()
  signature_refsets_selected <- signature_refsets %>% 
    filter(Signature_set_name==signature_set_name_input)
  seqmatrix_refdata_selected <- seqmatrix_refdata %>% filter(Study==study_input,Dataset==dataset_input,Profile == signature_refsets_selected$Profile[1])
  
}else{
  ## require for user input
  ## four parameters need: 
  ## Exposure File, Matrix File, Signature File, Genome Size, 
  exposure_refdata_selected <- read_delim('../Example_data/Sherlock_SBS96_exposure.txt',delim = '\t',col_names = T)
  signature_refsets_selected <- read_delim('../Example_data/Sherlock_SBS96_siganture.txt',delim = '\t',col_names = T)
  seqmatrix_refdata_selected <- read_delim('../Example_data/Sherlock_SBS96_matrix.txt',delim = '\t',col_names = T)
  genomesize <- 3101.976562
  cancer_type_user <- "Input" 
  
  
  ## format the input as suggest by the public dataset ##
  colnames(exposure_refdata_selected)[1] <- "Sample"
  colnames(seqmatrix_refdata_selected)[1] <- "MutationType"
  colnames(signature_refsets_selected)[1] <- "MutationType"
  seqmatrix_refdata_selected <- seqmatrix_refdata_selected %>% select(1,any_of(exposure_refdata_selected$Sample)) %>% profile_format_df()
  signature_refsets_selected <- signature_refsets_selected %>% select(MutationType,any_of(colnames(exposure_refdata_selected))) %>% profile_format_df()  
  
  exposure_refdata_selected <- exposure_refdata_selected %>% pivot_longer(cols = -Sample,names_to="Signature_name",values_to="Exposure") %>% mutate(Cancer_Type=cancer_type_user)
  signature_refsets_selected <- signature_refsets_selected %>% select(-Type,-SubType) %>% pivot_longer(cols = -MutationType,names_to="Signature_name",values_to="Contribution")
  seqmatrix_refdata_selected <- seqmatrix_refdata_selected%>% select(-Type,-SubType) %>% pivot_longer(cols = -MutationType,names_to="Sample",values_to="Mutations") %>% mutate(Cancer_Type=cancer_type_user)
  
}

## Tumor Overall Mutational Burden
data_input <- exposure_refdata_selected %>% 
  group_by(Cancer_Type,Sample) %>% 
  summarise(Burden=log10(sum(Exposure)/genomesize)) %>% 
  ungroup() 
# put this barplot on the web
TMBplot(data_input,output_plot = 'tmp.svg')


# Mutational Signature Activity
signature_name_input <- 'SBS4'

data_input <- exposure_refdata_selected %>% 
  filter(Signature_name==signature_name_input) %>% 
  group_by(Cancer_Type,Sample) %>% 
  summarise(Burden=log10(sum(Exposure)/genomesize)) %>% 
  ungroup() 
# put this barplot on the web
TMBplot(data_input,output_plot = 'tmp.svg',addnote = signature_name_input)

# Mutational Signature Assocaition
signature_name_input1 <- 'SBS2'
signature_name_input2 <- 'SBS13'
cancer_type_input <- NULL  ## toggle to select specific cancer type or combine all cancer type data (default)
signature_both <- FALSE ## toggle to choose samples with both signature detected


data_input <- left_join(
  exposure_refdata_selected %>% 
    filter(Signature_name==signature_name_input1) %>% 
    rename(Exposure1=Exposure) %>% 
    select(-Signature_name),
  
  exposure_refdata_selected %>% 
    filter(Signature_name==signature_name_input2) %>% 
    rename(Exposure2=Exposure) %>% 
    select(-Signature_name)
)

signature_association(data = data_input,cancer_type_input = cancer_type_input,signature_both = signature_both,output_plot = 'tmp.svg')


# Evaluating the Performance of Mutational Signature Decomposition --------
# load seqmatrix

exposure_refdata_input <- exposure_refdata_selected %>% mutate(Sample=paste0(Cancer_Type,"@",Sample)) %>% 
  select(Sample,Signature_name,Exposure) %>% 
  pivot_wider(id_cols = Sample,names_from=Signature_name,values_from=Exposure)

signature_refsets_input <- signature_refsets_selected %>% 
  select(MutationType,Signature_name,Contribution) %>% 
  pivot_wider(id_cols = MutationType,names_from=Signature_name,values_from=Contribution) %>% 
  arrange(MutationType) # have to sort the mutationtype

seqmatrix_refdata_input<- seqmatrix_refdata_selected %>% mutate(Sample=paste0(Cancer_Type,"@",Sample)) %>% 
  select(MutationType,Sample,Mutations) %>% 
  pivot_wider(id_cols = MutationType,names_from=Sample,values_from=Mutations) %>% 
  arrange(MutationType)  ## have to sort the mutation type


decompsite_input <- calculate_similarities(orignal_genomes = seqmatrix_refdata_input,signature = signature_refsets_input,signature_activaties = exposure_refdata_input)
decompsite_input <- decompsite_input %>% separate(col = Sample_Names,into = c('Cancer_Type','Sample'),sep = '@')
decompsite_input %>% write_delim('tmp.txt',delim = '\t',col_names = T)  ## put the link to download this table

decompsite_distribution(decompsite = decompsite_input,output_plot = 'tmp.svg') # put the distribution plot online.


# Landscape of Mutational Signature Activity
if(Data_Source == "Public_Data"){
cancer_type_input <- 'Skin-Melanoma'
}else {
  ## for input data, it will alwyas be "cancer_type_user" 
  cancer_type_input <- cancer_type_user
}

exposure_refdata_input <- exposure_refdata_selected %>% filter(Cancer_Type == cancer_type_input)%>% 
  select(Sample,Signature_name,Exposure) %>% 
  pivot_wider(id_cols = Sample,names_from=Signature_name,values_from=Exposure)

signature_refsets_input <- signature_refsets_selected %>% 
  select(MutationType,Signature_name,Contribution) %>% 
  pivot_wider(id_cols = MutationType,names_from=Signature_name,values_from=Contribution) %>% 
  arrange(MutationType) # have to sort the mutationtype

seqmatrix_refdata_input<- seqmatrix_refdata_selected %>% filter(Cancer_Type == cancer_type_input)%>% 
  select(MutationType,Sample,Mutations) %>% 
  pivot_wider(id_cols = MutationType,names_from=Sample,values_from=Mutations) %>% 
  arrange(MutationType)  ## have to sort the mutationtype
decompsite_input <- calculate_similarities(orignal_genomes = seqmatrix_refdata_input,signature = signature_refsets_input,signature_activaties = exposure_refdata_input)

cosinedata <- decompsite_input %>% select(Samples=Sample_Names,Similarity=Cosine_similarity)

data_input <- exposure_refdata_selected %>%
  filter(Cancer_Type==cancer_type_input) %>%
  select(Sample,Signature_name,Exposure) %>% 
  pivot_wider(id_cols = Sample,names_from=Signature_name,values_from=Exposure) %>% 
  rename(Samples=Sample)

data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)

sigdata <- data_input
## two parameters to add the two bars: vardata1, vardata1_cat, vardata2, vardata2_cat 
# studydata <- data_input %>% select(Samples) %>% mutate(Study=if_else((seq_along(Samples) %% 2 ==0), "A","B"))
# puritydata <-  data_input %>% select(Samples) %>% mutate(Purity=0)
# puritydata$Purity <- runif(n=length(puritydata$Purity), min=1e-12, max=.9999999999)
# highlight <-  c('SP124389','SP124273')
# Exposure_Clustering(sigdata = sigdata,studydata = studydata,studyname = "VAR1",puritydata = puritydata,purityname = 'VAR2',cosinedata = cosinedata,clustern=5,output_plot = 'tmp.svg' )
# Exposure_Clustering(sigdata = sigdata,cosinedata = cosinedata,clustern=5,output_plot = 'tmp.svg' )

# parameter: Cancer Type, Vardata_input_file
if(vardata_input_file){
  vardata_input <-  read_delim(vardata_input_file,delim = '\t',col_names = T) 
  
  vardata1_input <- vardata_input %>% select(1:2)
  colnames(vardata1_input) <- c('Samples','Study')
  vardata1_cat_input <- if_else(is.character(vardata1_input$Study),TRUE,FALSE)
  
  if(dim(vardata_input)[2]>2){
    vardata2_input <- vardata_input %>% select(1,3)
    colnames(vardata2_input) <- c('Samples','Purity')
    vardata2_cat_input <- if_else(is.character(vardata2_input$Purity),TRUE,FALSE)
  }else{
    vardata2_input <- NULL
    vardata2_cat_input <- FALSE
  }
  
}else{
  vardata_input <-  NULL
  vardata1_input <- NULL
  vardata1_cat_input <- NULL
  vardata2_input <- NULL
  vardata2_cat_input <- NULL
}

Exposure_Clustering(sigdata = sigdata,studydata = vardata1_input,studydata_cat = vardata1_cat_input,puritydata = vardata2_input,puritydata_cat = vardata2_cat_input, cosinedata = cosinedata,clustern=5,output_plot = 'tmp.svg' )

### prevalence plot
# parameter nmutation
nmutation_input <- 100
prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = 'tmp.svg')


