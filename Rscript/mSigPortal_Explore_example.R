set_wd()
library(tidyverse)
library(cowplot)
library(hrbrthemes)
library(ggsci)
library(ggrepel)
library(ggforce)

# source msigportal function ----------------------------------------------
source('Sigvisualfunc.R')

# load reference signatures files -----------------------------------------
load('signature_refsets.RData')




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

path_profile <- '~/NIH-Work/MutationSignature/mSigPortal/CBIIT/Reference_Signatures_dataset/Profiles/Reference_Signature_Profiles_SVG/'
signature_profile_files <- signature_refsets %>% select(Source,Profile,Signature_set_name,Dataset,Signature_name) %>% unique() %>% mutate(Path=str_replace_all(Signature_set_name," ","_"),Path=str_remove_all(Path,"[()]"),Path=paste0(path_profile,'Profiles/',Path,"/",Signature_name,".svg"))

## multiple selected profiles
signature_source_input <- "Reference_signatures"
profile_name_input <- "SBS96"
signatureset_name_input <- "COSMIC v3 Signatures (SBS)"
dataset_input <- "WGS"
signature_name_input <- "SBS1"
svgfile_selected <- signature_profile_files %>% filter(Source==signature_source_input,Profile==profile_name_input,Signature_set_name==signatureset_name_input,Dataset==dataset_input,Signature_name==signature_name_input) %>% pull(Path)



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


# section4: Mutational signatures comparisons
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









# load exposure data files 


# exposure association with genome data. 



