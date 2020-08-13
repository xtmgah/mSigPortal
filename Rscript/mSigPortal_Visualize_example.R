set_wd()
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


### Cosine Similarity tab ###
# section 1: Cosine similarity within samples # 
# two parameters need: Profile Type and Matrix Size # 
profiletype <- "SBS"
matrixsize <- "SBS-96"
matrixsize <- str_remove(matrixsize,"-")
# data input can be filed use profiletype and matrixsize from the summary file jian generated. 
data_input <- read_delim(paste0('example_results/output/',profiletype,'/22520169-fe27-4549-ad3e-e78d1812c45e.',matrixsize,'.all'),delim = '\t')
## remove the false data from the result. Don't know why the web generated two additional columns: False and seen, i think you trying to find the bug and remove. 
data_input <- data_input %>% select(-False,-seen)
# Heatmap of cosine similarity within samples  and put on the web---------------------------
cos_sim_res1=cos_sim_df(data_input,data_input) 
plot_cosine_heatmap_df(cos_sim_res1,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
# a link to download the cosine similarity bellow the plot 
# you could rename the file name if you need
cos_sim_res1 %>% write_delim('cos_sim_res1.txt',delim = '\t',col_names = T)


# section 2: Cosine similarity  to reference signatures # 
# Two parameters need: Profile Type, Reference Signature Set
# Profile Type only support SBS, DBS, ID
profiletype <- "SBS"
profilename <- if_else(profiletype == "SBS","SBS96",if_else(profiletype=="DBS","DBS78",if_else(profiletype=="ID","ID83",NA_character_)))
#all the option of Reference Signature Set
signature_refsets %>% filter(Profile==profilename) %>% pull(Signature_set_name) %>% unique() ## put on the dropdown list 
# if selected "COSMIC v3 Signatures (SBS)"
signaturesetname <- "COSMIC v3 Signatures (SBS)"
signature_refsets_input <- signature_refsets %>% filter(Profile==profilename,Signature_set_name==signaturesetname)
sigref <- signature_refsets_input %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)

matrixsize <- profilename
# load SBS-96/DBS-78/ID-83 catalog files generated from previous result according to the profilename -----------
data_input <- read_delim(paste0('example_results/output/',profiletype,'/22520169-fe27-4549-ad3e-e78d1812c45e.',matrixsize,'.all'),delim = '\t')
## remove the false data from the result. Don't know why the web generated two additional columns: False and seen, i think you trying to find the bug and remove. 
data_input <- data_input %>% select(-False,-seen)

# Heatmap of cosine similarity to reference set signature and put on the web---------------------------
cos_sim_res2=cos_sim_df(data_input,sigref) 
plot_cosine_heatmap_df(cos_sim_res2,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
# a link to download the cosine similarity bellow the plot 
# you could rename the file name if you need
cos_sim_res2 %>% write_delim('cos_sim_res2.txt',delim = '\t',col_names = T)









### Profile Comparison tab ### 
# section 1: Cosine similarity within samples # 
# three parameters need: Profile Type, Sample Name1 and Sample Name2 # 
profiletype <- "SBS"
sample_name1 <-  "SB749362"
sample_name2 <-  "SC420396"

matrixsize <- if_else(profiletype == "SBS","SBS96",if_else(profiletype=="DBS","DBS78",if_else(profiletype=="ID","ID83",NA_character_)))
data_input <- read_delim(paste0('example_results/output/',profiletype,'/22520169-fe27-4549-ad3e-e78d1812c45e.',matrixsize,'.all'),delim = '\t')
## remove the false data from the result. Don't know why the web generated two additional columns: False and seen, i think you trying to find the bug and remove. 
data_input <- data_input %>% select(-False,-seen)

profile1 <- data_input %>% select(MutationType,contains(sample_name1))
profile2 <- data_input %>% select(MutationType,contains(sample_name2))
# put the plot on the web
plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,output_plot = 'tmp.svg')


# section 2: Comparison to reference signatures # 
# four parameters need: “Profile Type”, “Sample Name”, “Reference Signature Set” and “Compare Single Signature or Combined Signatures” # 
sample_name <-  "SC420396"
profiletype <- "ID"
profilename <- if_else(profiletype == "SBS","SBS96",if_else(profiletype=="DBS","DBS78",if_else(profiletype=="ID","ID83",NA_character_)))
matrixsize <- profilename

data_input <- read_delim(paste0('example_results/output/',profiletype,'/22520169-fe27-4549-ad3e-e78d1812c45e.',matrixsize,'.all'),delim = '\t')
profile1 <- data_input %>% select(MutationType,one_of(sample_name))

#all the option of Reference Signature Set
signature_refsets %>% filter(Profile==profilename) %>% pull(Signature_set_name) %>% unique() ## put on the dropdown list 
# if selected "COSMIC v3 Signatures (SBS)"
signaturesetname <- "COSMIC v3 Signatures (ID)"
signature_refsets_input <- signature_refsets %>% filter(Profile==profilename,Signature_set_name==signaturesetname)
sigref <- signature_refsets_input %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)

user_input <- "ID1" #user_input <- "0.8*SBS5;0.1*SBS1;0.1*SBS3"

if(str_detect(user_input,";")){
  profile2 <- signature_sum_operation(sigdatabase = signature_refsets_input,sigsetname = signaturesetname,formulax = user_input,outsigname = "Reconstructed")
  profile_names = c("Original","Reconstructed")
}else {
  profile2 <- sigref %>% select(MutationType,one_of(user_input))
  profile_names = c(colnames(profile1)[2],colnames(profile2)[2])
}

plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,profile_names = profile_names,output_plot = 'tmp.svg')










### “Principal Component Analysis” ###
# Two parameters need: Profile Type, Reference Signature Set
# Profile Type only support SBS, DBS, ID
profiletype <- "SBS"
profilename <- if_else(profiletype == "SBS","SBS96",if_else(profiletype=="DBS","DBS78",if_else(profiletype=="ID","ID83",NA_character_)))
matrixsize <- profilename
#all the option of Reference Signature Set
signature_refsets %>% filter(Profile==profilename) %>% pull(Signature_set_name) %>% unique() ## put on the dropdown list 
# if selected "COSMIC v3 Signatures (SBS)"
signaturesetname <- "COSMIC v3 Signatures (SBS)"
signature_refsets_input <- signature_refsets %>% filter(Profile==profilename,Signature_set_name==signaturesetname)
sigref <- signature_refsets_input %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)


# load SBS-96/DBS-78/ID-83 catalog files generated from previous result according to the profilename -----------
data_input <- read_delim(paste0('example_results/output/',profiletype,'/22520169-fe27-4549-ad3e-e78d1812c45e.',matrixsize,'.all'),delim = '\t')
## remove the false data from the result. Don't know why the web generated two additional columns: False and seen, i think you trying to find the bug and remove. 
data_input <- data_input %>% select(-False,-seen)


# PCA plot ----------------------------------------------------------------

mdata_input <- as.matrix(data_input[,-1])
rownames(mdata_input) <- data_input$MutationType

res.pca <- prcomp(t(mdata_input), scale = FALSE,center = FALSE)

xleng <- dim(res.pca$x)[2]*0.5+1
yleng <- 4
pcap1 <- fviz_eig(res.pca,ncp = 10,main = "")
ggsave(filename = 'tmp.svg',plot = pcap1,width = xleng,height = yleng)

if((dim(data_input)[2]-1) > 20){
  pcap2 <- fviz_pca_ind(res.pca,axes = c(1,2),geom="point",col.ind="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
}else {
  pcap2 <- fviz_pca_ind(res.pca,axes = c(1,2),col.ind="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
}
ggsave(filename = 'tmp.svg',plot = pcap2,width = 10,height = 7)


## a link to download pca data1
res.pca$x %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim('PCA1.txt',delim = '\t',col_names = T)


pcap3 <- fviz_pca_var(res.pca,axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
)
ggsave(filename = 'tmp.svg',plot = pcap3,width = 10,height = 7)


# a link to download pca data2
res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim('PCA2.txt',delim = '\t',col_names = T)



## heatmap between PCs and signatures
sigpca <-  res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "MutationType")
cos_sim_res3=cos_sim_df(sigpca,sigref) 
# put this heatmap on the web
plot_cosine_heatmap_df(cos_sim_res3,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
# a link to download the cosine similarity bellow the plot 
# you could rename the file name if you need
cos_sim_res2 %>% write_delim('cos_sim_res3.txt',delim = '\t',col_names = T)










# t-SNE plots, coming soon ###







