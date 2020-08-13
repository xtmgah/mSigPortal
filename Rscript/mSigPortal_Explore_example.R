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

###  Signature Explore ###
# summary of number of signatures or signature sets -----------------------
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


# put the follow piechart on website 
signature_piechart(nsig_data,sigsetcolor,output_plot = 'tmp.svg')


# add signature SVG plots for all the signatures and put on line
# three parameters to include: “Profile Type”, “Reference Signature Set”, “Signature Name” 
profiletype <- "SBS" 
signaturesetname <- "COSMIC v3 Signatures (SBS)"
signaturename <- "SBS5"
## according to these three parameters to pull out the svg files



# cosine similarity heatmap between two signatures sets -------------------
# The parameters will be “Matrix Size”, “Reference Signature Set1” and “Reference Signature Set2”. 
matrixsize <- "SBS96" # profile type
signaturesetname1 <- "Environmental Mutagen Signatures (SBS)" 
signaturesetname2 <- "COSMIC v3 Signatures (SBS)"

# the availabe options for signaturesetname1 and signaturesetname2 will be:
signature_refsets %>% filter(Profile==matrixsize)%>% pull(Signature_set_name) %>% unique()


sigrefset1 <- signature_refsets %>%
  filter(Profile==matrixsize,Signature_set_name==signaturesetname1) %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)

sigrefset2 <- signature_refsets %>%
  filter(Profile==matrixsize,Signature_set_name==signaturesetname2) %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)

cos_sim_res=cos_sim_df(sigrefset1,sigrefset2) 

# put this heatmap on the web
plot_cosine_heatmap_df(cos_sim_res,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
# add a link bellow the heatmap
cos_sim_res %>% write_delim('signature_cos_sim_res.txt',delim = '\t',col_names = T)





## A comparison of two reference signatures
# There will be five parameters: “Profile Type”,  “Reference Signature Set1”, “Signature Name1”, “Reference Signature Set2”, “Signature Name2”; 
matrixsize <- "SBS96" # profile type
signaturesetname1 <- "COSMIC v3 Signatures (SBS)"
signaturesetname2 <- "COSMIC v3 Signatures (SBS)"
signaturename1 <-  "SBS1"
signaturename2 <- "SBS5"

# similarly, the availabe options for signaturesetname1 and signaturesetname2 will be:
signature_refsets %>% filter(Profile==matrixsize)%>% pull(Signature_set_name) %>% unique()
# the available options for the signature names will be:
signature_refsets %>% filter(Profile==matrixsize,Signature_set_name==signaturesetname1)%>% pull(Signature_name) %>% unique()
signature_refsets %>% filter(Profile==matrixsize,Signature_set_name==signaturesetname2)%>% pull(Signature_name) %>% unique()


profile1 <- signature_refsets %>%
  filter(Profile==matrixsize,Signature_set_name==signaturesetname1) %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution) %>% 
  select(MutationType,one_of(signaturename1))

profile2 <- signature_refsets %>%
  filter(Profile==matrixsize,Signature_set_name==signaturesetname2) %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution) %>% 
  select(MutationType,one_of(signaturename2))

# put this plot on the web:
plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,output_plot = 'tmp.svg')





# load exposure data files 


# exposure association with genome data. 



