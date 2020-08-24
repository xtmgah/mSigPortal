set_wd()
library(tidyverse)
library(cowplot)
library(hrbrthemes)
library(ggsci)
library(ggrepel)
library(factoextra)


# Source msigportal function ----------------------------------------------
source('Sigvisualfunc.R')



# define Data_Source  -----------------------------------------------------
if(Data_Source != "Public_Data"){
  ## read matrix summary list 
  foldername <- "example_results/" 
  matrixfiles <- read_csv(paste0(foldername,'matrix_files_list.txt'),col_names = T) 
  ## list of profile type and matrix size
  matrixfiles %>% select(Profile_Type,Matrix_Size)
  svgfiles <- read_csv(paste0(foldername,'/svg_files_list.txt'),col_names = T) 
  ### svgfiles and matrixfiles used as input data
  
}else{
  #parameters for the public data 
  # list for the public data
  load('Data/seqmatrix_refdata_info.RData')
  seqmatrix_refdata_info2
  foldername <- "~/NIH-Work/MutationSignature/mSigPortal/CBIIT/Reference_paper_data/The repertoire of mutational signatures in human cancer/SampleMatrix/Matrix/" 
  ## read svg summary list from public data
  svgfiles <-seqmatrix_refdata_info %>% mutate(Path=paste0(foldername,Path))
  
  ## parameters on the left panel
  study <- "PCAWG"
  cancer_type <-  'Uterus-AdenoCA'
  experimental_strategy <- 'WGS'
  svgfiles_public <- svgfiles %>% filter(Study == study, Cancer_Type==cancer_type,Dataset == experimental_strategy)
  load('Data/seqmatrix_refdata.RData')
  seqmatrix_refdata_public <- seqmatrix_refdata %>% filter(Study == study, Cancer_Type==cancer_type,Dataset == experimental_strategy)
  ### svgfiles_public and seqmatrix_refdata_public used as input data 
}



# Mutational Profiles -----------------------------------------------------
if(Data_Source != "Public_Data"){
  sample_name <-  'SC349747' 
  profile_type <- 'SBS'
  matrix_size <- '96'
  tag_filter <- 'PASS' ## not working for the public data; 
  if(tag_filter == "NA"){
    svg_selected <- svgfiles %>% filter(Sample_Name==sample_name,Profile_Type==profile_type,Matrix_Size==matrix_size,is.na(Filter)) %>% pull(Path)
  } else {
    svg_selected <- svgfiles %>% filter(Sample_Name==sample_name,Profile_Type==profile_type,Matrix_Size==matrix_size,Filter==tag_filter) %>% pull(Path)
  }
  # shown svg_selected on the web
}else {
  ## paramters on the right panel
  sample_name <- 'SP95646'
  profile_type <- 'SBS'
  matrix_size <- '96'
  tag_filter <- 'NA'  # for public data, the filter always NA
  profile_name <- paste0(profile_type,matrix_size)
  
  svg_selected <- svgfiles_public %>% filter(Profile==profile_name,Sample==sample_name) %>% pull(Path)
  # shown svg_selected on the web
}


# Load R object for the rest tables -------------------------------------------
# load public sigantures
load('Data/signature_refsets.RData')



# Cosine Similarity tab ---------------------------------------------------

# section 1: Cosine similarity within samples  ----------------------------
# two parameters need: Profile Type and Matrix Size # 
profile_type_input <- "SBS"
matrix_size_input <- "96"
profile_name <- paste0(profile_type_input,matrix_size_input)

if(Data_Source != "Public_Data"){
  # data input can be filed use profiletype and matrixsize from the summary file jian generated. 
  matrixfile_selected <- matrixfiles %>% filter(Profile_Type==profile_type_input,Matrix_Size==matrix_size_input) %>% pull(Path)
  data_input <- read_delim(matrixfile_selected,delim = '\t')
  ## removd the blank profile
  #data_input %>% select(-MutationType) %>% summarise(across(where(is.numeric),sum)) %>% select_if(~ sum(.)==0) %>% colnames()
  data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)
  
  # Heatmap of cosine similarity within samples  and put on the web---------------------------
  cos_sim_res1=cos_sim_df(data_input,data_input) 
  plot_cosine_heatmap_df(cos_sim_res1,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
  # a link to download the cosine similarity bellow the plot 
  # you could rename the file name if you need
  cos_sim_res1 %>% write_delim('cos_sim_res1.txt',delim = '\t',col_names = T)
  
} else{
  data_input <- seqmatrix_refdata_public %>% filter(Profile==profile_name) %>% 
    select(MutationType,Sample,Mutations) %>% 
    pivot_wider(id_cols = MutationType,names_from=Sample,values_from=Mutations)
  
  data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)
  cos_sim_res1=cos_sim_df(data_input,data_input) 
  plot_cosine_heatmap_df(cos_sim_res1,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
  cos_sim_res1 %>% write_delim('cos_sim_res1.txt',delim = '\t',col_names = T)
  
}



# section 2: Cosine similarity  to reference signatures -------------------
# Two parameters need: Profile Type, Reference Signature Set
# Profile Type only support SBS, DBS, ID
profile_type_input <- "SBS"
profile_name <- if_else(profile_type_input == "SBS","SBS96",if_else(profile_type_input=="DBS","DBS78",if_else(profile_type_input=="ID","ID83",NA_character_)))
#all the option of Reference Signature Set
signature_refsets %>% filter(Profile==profile_name) %>% pull(Signature_set_name) %>% unique() ## put on the dropdown list 
# if selected "COSMIC v3 Signatures (SBS)"
signatureset_name_input <- "COSMIC v3 Signatures (SBS)"
signature_refsets_data <- signature_refsets %>% filter(Profile==profile_name,Signature_set_name==signatureset_name_input)
sigref_data <- signature_refsets_data %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)

matrix_size <- str_remove(str_remove(str_remove(profile_name,"SBS"),"ID"),"DBS")

if(Data_Source != "Public_Data"){
  matrixfile_selected <- matrixfiles %>% filter(Profile_Type==profile_type_input,Matrix_Size==matrix_size) %>% pull(Path)
  data_input <- read_delim(matrixfile_selected,delim = '\t')
  data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)
  
  # Heatmap of cosine similarity to reference set signature and put on the web---------------------------
  cos_sim_res2=cos_sim_df(data_input,sigref_data) 
  plot_cosine_heatmap_df(cos_sim_res2,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
  # a link to download the cosine similarity bellow the plot 
  # you could rename the file name if you need
  cos_sim_res2 %>% write_delim('cos_sim_res2.txt',delim = '\t',col_names = T)
  
} else {
  data_input <- seqmatrix_refdata_public %>% filter(Profile==profile_name) %>% 
    select(MutationType,Sample,Mutations) %>% 
    pivot_wider(id_cols = MutationType,names_from=Sample,values_from=Mutations)
  data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)
  
  cos_sim_res2=cos_sim_df(data_input,sigref_data) 
  plot_cosine_heatmap_df(cos_sim_res2,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
  cos_sim_res2 %>% write_delim('cos_sim_res2.txt',delim = '\t',col_names = T)
}




# section 3: Cosine similarity to Public data -----------------------------
# find the common profile between data and seqmatrix

if(Data_Source != "Public_Data"){
  matrixfiles %>% mutate(Profile_Name=paste0(Profile_Type,Matrix_Size)) %>% filter(Profile_Name %in% seqmatrix_refdata_info$Profile) %>% pull(Profile_Name)
  profile_name_input <- "SBS1536"
  study_input <- "PCAWG"
  cancer_type_input <- 'Uterus-AdenoCA'
  
  ## input data
  matrixfile_selected <- matrixfiles %>% mutate(Profile_Name=paste0(Profile_Type,Matrix_Size)) %>% filter(Profile_Name == profile_name_input) %>% pull(Path)
  data_input <- read_delim(matrixfile_selected,delim = '\t')
  data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)
  
  ## seqmatrix data from public data
  sigmatrix_data <- seqmatrix_refdata %>% 
    filter(Profile==profile_name_input,Study==study_input,Cancer_Type==cancer_type_input) %>% 
    select(Sample,MutationType,Mutations) %>% 
    pivot_wider(names_from = Sample,values_from=Mutations)
  
  # Heatmap of cosine similarity to public seqmatrix data and put on the web---------------------------
  cos_sim_res3=cos_sim_df(data_input,sigmatrix_data) 
  plot_cosine_heatmap_df(cos_sim_res3,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
  # a link to download the cosine similarity bellow the plot 
  # you could rename the file name if you need
  cos_sim_res3 %>% write_delim('cos_sim_res3.txt',delim = '\t',col_names = T)
}




# Profile Comparison tab --------------------------------------------------

# section 1: Profile Comparison within samples ----------------------------
# three parameters need: Profile Type, Sample Name1 and Sample Name2 # 
profile_type_input <- "SBS"
matrix_size <- if_else(profile_type_input == "SBS","96",if_else(profile_type_input=="DBS","78",if_else(profile_type_input=="ID","83",NA_character_)))
profile_name <- paste0(profile_type_input,matrix_size_input)

if(Data_Source != "Public_Data"){
  sample_name_input1 <-  "SB749362"
  sample_name_input2 <-  "SC420396"
  
  ## input data
  matrixfile_selected <- matrixfiles %>% filter(Profile_Type==profile_type_input,Matrix_Size==matrix_size) %>% pull(Path)
  data_input <- read_delim(matrixfile_selected,delim = '\t')
  data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)
  
  profile1 <- data_input %>% select(MutationType,one_of(sample_name_input1))
  profile2 <- data_input %>% select(MutationType,one_of(sample_name_input2))
  # put the plot on the web
  plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,output_plot = 'tmp.svg')
}else {
  sample_name_input1 <-  "SP90989"
  sample_name_input2 <-  "SP90725"
  data_input <- seqmatrix_refdata_public %>% filter(Profile==profile_name) %>% 
    select(MutationType,Sample,Mutations) %>% 
    filter(Sample %in% c(sample_name_input1,sample_name_input2)) %>% 
    pivot_wider(id_cols = MutationType,names_from=Sample,values_from=Mutations)
  
  profile1 <- data_input %>% select(MutationType,one_of(sample_name_input1))
  profile2 <- data_input %>% select(MutationType,one_of(sample_name_input2))
  plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,output_plot = 'tmp.svg')
}





# section 2: Comparison to reference signatures ---------------------------

# four parameters need: “Profile Type”, “Sample Name”, “Reference Signature Set” and “Compare Single Signature or Combined Signatures” # 
profile_type_input <- "ID"
profile_name <- if_else(profile_type_input == "SBS","SBS96",if_else(profile_type_input=="DBS","DBS78",if_else(profile_type_input=="ID","ID83",NA_character_)))
matrix_size <- if_else(profile_type_input == "SBS","96",if_else(profile_type_input=="DBS","78",if_else(profile_type_input=="ID","83",NA_character_)))

#all the option of Reference Signature Set
signature_refsets %>% filter(Profile==profile_name) %>% pull(Signature_set_name) %>% unique() ## put on the dropdown list 
# if selected "COSMIC v3 Signatures (SBS)"
signatureset_name_input <- "COSMIC v3 Signatures (ID)"
signature_refsets_input <- signature_refsets %>% filter(Profile==profile_name,Signature_set_name==signatureset_name_input)
sigref_data <- signature_refsets_input %>% 
  select(Signature_name,MutationType,Contribution) %>% 
  pivot_wider(names_from = Signature_name,values_from=Contribution)

user_input <- "ID4" #user_input <- "0.8*SBS5;0.1*SBS1;0.1*SBS3"

if(str_detect(user_input,";")){
  profile2 <- signature_sum_operation(sigdatabase = signature_refsets_input,sigsetname = signatureset_name_input,formulax = user_input,outsigname = "Reconstructed")
  profile_names = c("Original","Reconstructed")
}else {
  profile2 <- sigref_data %>% select(MutationType,one_of(user_input))
  profile_names = c(colnames(profile1)[2],colnames(profile2)[2])
}

if(Data_Source != "Public_Data"){
  sample_name_input <-  "SC420396"
  ## input data
  matrixfile_selected <- matrixfiles %>% filter(Profile_Type==profile_type_input,Matrix_Size==matrix_size) %>% pull(Path)
  data_input <- read_delim(matrixfile_selected,delim = '\t')
  data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)
  profile1 <- data_input %>% select(MutationType,one_of(sample_name_input))
  
  plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,profile_names = profile_names,output_plot = 'tmp.svg')
}else {
  sample_name_input <-  "SP90989"
  profile1 <- seqmatrix_refdata_public %>% filter(Profile==profile_name) %>% 
    select(MutationType,Sample,Mutations) %>% 
    filter(Sample %in% c(sample_name_input)) %>% 
    pivot_wider(id_cols = MutationType,names_from=Sample,values_from=Mutations)
  
  plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,profile_names = profile_names,output_plot = 'tmp.svg')
}









# section 3: Profile Comparison to Public data ----------------------------
if(Data_Source != "Public_Data"){
  # find the common profile between data and seqmatrix
  matrixfiles %>% mutate(Profile_Name=paste0(Profile_Type,Matrix_Size)) %>% filter(Profile_Name %in% seqmatrix_refdata_info$Profile) %>% pull(Profile_Name)
  profile_name_input <- "SBS96"
  sample_name_input1 <-  "SC420396"
  study_input <- "PCAWG"
  cancer_type_input <- 'Uterus-AdenoCA'
  sample_name_input2 <- 'SP95454'
  
  ## input data
  matrixfile_selected <- matrixfiles %>% mutate(Profile_Name=paste0(Profile_Type,Matrix_Size)) %>% filter(Profile_Name == profile_name_input) %>% pull(Path)
  data_input <- read_delim(matrixfile_selected,delim = '\t')
  data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)
  profile1 <- data_input %>% select(MutationType,one_of(sample_name_input1))
  
  ## seqmatrix data from public data
  profile2 <- seqmatrix_refdata %>% 
    filter(Profile==profile_name_input,Study==study_input,Cancer_Type==cancer_type_input) %>% 
    select(Sample,MutationType,Mutations) %>% 
    filter(Sample==sample_name_input2) %>% 
    pivot_wider(names_from = Sample,values_from=Mutations) 
  
  # put the plot on the web
  plot_compare_profiles_diff(profile1,profile2,condensed = FALSE,output_plot = 'tmp.svg')
}



# “Principal Component Analysis” ------------------------------------------

# section 1 PCA within sample  --------------------------------------------


# parameters need: Profile Type, Matrix_Size,
# Profile Type can be any profile but compared to signature only support common types
profile_type_input <- "DBS"
matrix_size_input <- "78"
profile_name <- paste0(profile_type_input,matrix_size_input)

if(Data_Source != "Public_Data"){
  ## input data
  matrixfile_selected <- matrixfiles %>% filter(Profile_Type==profile_type_input,Matrix_Size==matrix_size_input) %>% pull(Path)
  data_input <- read_delim(matrixfile_selected,delim = '\t')
  
}else{
  data_input <- seqmatrix_refdata_public %>% filter(Profile==profile_name) %>% 
    select(MutationType,Sample,Mutations) %>% 
    pivot_wider(id_cols = MutationType,names_from=Sample,values_from=Mutations)
}

data_input <- data_input %>% select_if(~ !is.numeric(.)|| sum(.)>0)



# PCA plot ----------------------------------------------------------------
mdata_input <- as.matrix(data_input[,-1])
rownames(mdata_input) <- data_input$MutationType

res.pca <- prcomp(t(mdata_input), scale = FALSE,center = FALSE)

xleng <- dim(res.pca$x)[2]*0.25+1
xleng <- if_else(xleng>4,4,xleng)
yleng <- 4
pcap1 <- fviz_eig(res.pca,ncp = 10,main = "",addlabels = T)
ggsave(filename = 'tmp.svg',plot = pcap1,width = xleng,height = yleng)

if((dim(data_input)[2]-1) > 35){
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


# parameters:  Reference Signature Set
#all the option of Reference Signature Set
## only work if the profile included in singatures dataset
if(profile_name %in% unique(signature_refsets$Profile))
{
  signature_refsets %>% filter(Profile==profile_name) %>% pull(Signature_set_name) %>% unique() ## put on the dropdown list 
  # if selected "COSMIC v3 Signatures (DBS)"
  signatureset_name <- "COSMIC v3 Signatures (DBS)"
  signature_refsets_input <- signature_refsets %>% filter(Profile==profile_name,Signature_set_name==signatureset_name)
  sigref_data <- signature_refsets_input %>% 
    select(Signature_name,MutationType,Contribution) %>% 
    pivot_wider(names_from = Signature_name,values_from=Contribution)
  
  ## heatmap between PCs and signatures
  sigpca_data <-  res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "MutationType")
  cos_sim_res4=cos_sim_df(sigpca_data,sigref_data) 
  # put this heatmap on the web
  plot_cosine_heatmap_df(cos_sim_res4,cluster_rows = TRUE,plot_values = FALSE,output_plot = 'tmp.svg')
  # a link to download the cosine similarity bellow the plot 
  # you could rename the file name if you need
  cos_sim_res2 %>% write_delim('cos_sim_res3.txt',delim = '\t',col_names = T)
  
}




#  section 2 PCA together with public data --------------------------------
if(Data_Source != "Public_Data"){
  # parameters need: Profile Type, Matrix_Size,
  # Profile Type can be any profile but compared to signature only support common types
  # only work with overlap profile 
  matrixfiles %>% mutate(Profile_Name=paste0(Profile_Type,Matrix_Size)) %>% filter(Profile_Name %in% seqmatrix_refdata_info$Profile) %>% pull(Profile_Name)
  profile_name_input <- "SBS1536"
  study_input <- "PCAWG"
  cancer_type_input <- 'Uterus-AdenoCA'
  
  ## input data
  matrixfile_selected <- matrixfiles %>% mutate(Profile_Name=paste0(Profile_Type,Matrix_Size))%>% filter(Profile_Name==profile_name_input) %>% pull(Path)
  data_input1 <- read_delim(matrixfile_selected,delim = '\t')
  data_input1 <- data_input1 %>% select_if(~ !is.numeric(.)|| sum(.)>0)
  
  ## seqmatrix data from public data
  data_input2 <- seqmatrix_refdata %>% 
    filter(Profile==profile_name_input,Study==study_input,Cancer_Type==cancer_type_input) %>% 
    select(Sample,MutationType,Mutations) %>% 
    pivot_wider(names_from = Sample,values_from=Mutations) 
  
  data_input <- data_input1 %>% left_join(data_input2)
  
  indcolors <- c(rep("Input",length(colnames(data_input1)[-1])),rep(cancer_type_input,length(colnames(data_input2)[-1])))
  names(indcolors) <- c(colnames(data_input1)[-1],colnames(data_input2)[-1])
  
  
  # PCA plot ----------------------------------------------------------------
  mdata_input <- as.matrix(data_input[,-1])
  rownames(mdata_input) <- data_input$MutationType
  
  res.pca <- prcomp(t(mdata_input), scale = FALSE,center = FALSE)
  xleng <- dim(res.pca$x)[2]*0.25+1
  xleng <- if_else(xleng>4,4,xleng)
  yleng <- 4
  pcap1 <- fviz_eig(res.pca,ncp = 10,main = "",addlabels = T)
  ggsave(filename = 'tmp.svg',plot = pcap1,width = xleng,height = yleng)
  
  
  if((dim(data_input)[2]-1) > 35){
    pcap2 <- fviz_pca_ind(res.pca,axes = c(1,2),geom="point",col.ind=indcolors,repel = TRUE)+ scale_color_brewer(palette="Set1") 
  }else {
    pcap2 <- fviz_pca_ind(res.pca,axes = c(1,2),col.ind=indcolors,repel = TRUE)+scale_color_brewer(palette="Set1") 
  }
  ggsave(filename = 'tmp.svg',plot = pcap2,width = 10,height = 7)
  
  ## a link to download pca data1
  res.pca$x %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim('PCA3.txt',delim = '\t',col_names = T)
  
  pcap3 <- fviz_pca_var(res.pca,axes = c(1, 2),
                        col.var = "contrib", # Color by contributions to the PC
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = FALSE     # Avoid text overlapping
  )
  ggsave(filename = 'tmp.svg',plot = pcap3,width = 10,height = 7)
  
  # a link to download pca data2
  res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim('PCA4.txt',delim = '\t',col_names = T)
  
}





# Download ----------------------------------------------------------------
## need to specific download for public data in future. 

