library(tidyverse)
library(cowplot)
library(hrbrthemes)
library(ggsci)
library(ggrepel)
library(factoextra)


# source msigportal function ----------------------------------------------
source('Sigvisualfunc.R')

# load reference signatures files -----------------------------------------
load('signature_refsets.RData')


# on the web, we need the dropdown list for signature Profile Type and reference set name --------
profilename <- "SBS"
profilename2 <- if_else(profilename == "SBS","SBS96",if_else(profilename=="DBS","DBS78",if_else(profilename=="ID","ID84",NA_character_)))
signaturesetname <- "COSMIC v3 Signatures (SBS)"

signature_refsets_input <- signature_refsets %>% filter(Profile==profilename2,Signature_set_name==signaturesetname)


# load SBS-96/DBS-78/ID-83 catalog file generated from previous result according to the profilename -----------
if(profilename == "SBS"){
  data_input <- read_delim('example_results/output/SBS/22520169-fe27-4549-ad3e-e78d1812c45e.SBS96.all',delim = '\t')
}
if(profilename == "DBS"){
  data_input <- read_delim('example_results/output/DBS/22520169-fe27-4549-ad3e-e78d1812c45e.DBS78.all',delim = '\t')
}
if(profilename == "ID"){
  data_input <- read_delim('example_results/output/ID/22520169-fe27-4549-ad3e-e78d1812c45e.ID83.all',delim = '\t')
}

## remove the false data from the result. Don't know why the web generated two additonal columns: False and seen, i think you trying to find the bug and remove. 
data_input <- data_input %>% select(-False,-seen)

# Heatmap of cosine similarity to reference set ---------------------------

sigref <- signature_refsets_input %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)

data_input
cos_sim_res=cos_sim_df(data_input,sigref) 

# put this heatmap on the web
plot_cosine_heatmap_df(cos_sim_res,cluster_rows = TRUE,plot_values = FALSE)


# on the web, dropdown list for the input sample and (specific signature or a formula), put the barplot on the web
sample_input <-  "SC420396"

# option1: specific signature
siganture_name_input <- "SBS5"
profile1 <- data_input %>% select(MutationType,one_of(sample_input))
profile2 <- sigref %>% select(MutationType,one_of(siganture_name_input))
plot_compare_profiles_diff(profile1,profile2,condensed = FALSE)

# option2: formula input
formula_input <- "0.8*SBS5;0.1*SBS1;0.1*SBS3"
profile1 <- data_input %>% select(MutationType,one_of(sample_input))
profile2 <- signature_sum_operation(sigdatabase = signature_refsets_input,sigsetname = signaturesetname,formulax = formula_input,outsigname = "Reconstructed")
plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,profile_names = c("Original","Reconstructed"))


# PCA plot ----------------------------------------------------------------

data_input
mdata_input <- as.matrix(data_input[,-1])
rownames(mdata_input) <- data_input$MutationType

res.pca <- prcomp(t(mdata_input), scale = FALSE,center = FALSE)
fviz_eig(res.pca,ncp = 10)

if((dim(data_input)[2]-1) > 20){
  fviz_pca_ind(res.pca,axes = c(1,2),geom="point",col.ind="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
}else {
  fviz_pca_ind(res.pca,axes = c(1,2),col.ind="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
}


fviz_pca_var(res.pca,axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
)

## output PCAs
res.pca$x 
res.pca$rotation

## heatmap between PCs and signatures
sigpca <-  res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "MutationType")
cos_sim_res=cos_sim_df(sigpca,sigref) 

# put this heatmap on the web
plot_cosine_heatmap_df(cos_sim_res,cluster_rows = TRUE,plot_values = FALSE)



# t-SNE plots, coming soon ###







