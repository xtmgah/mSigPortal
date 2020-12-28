require(tidyverse)
require(ggtext)
require(ggforce)




# Conflicted package ------------------------------------------------------
resolve_conflicts <- function(){
  require(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
  #filter <- dplyr::filter
  #select <- dplyr::select
}


### reformat the signature name ###
signames <- 
  tibble(
    SBSname=c(
      "SBS1",
      "SBS2",
      "SBS3",
      "SBS4",
      "SBS5",
      "SBS6",
      "SBS7a",
      "SBS7b",
      "SBS7c",
      "SBS7d",
      "SBS8",
      "SBS9",
      "SBS10a",
      "SBS10b",
      "SBS11",
      "SBS12",
      "SBS13",
      "SBS14",
      "SBS15",
      "SBS16",
      "SBS17a",
      "SBS17b",
      "SBS18",
      "SBS19",
      "SBS20",
      "SBS21",
      "SBS22",
      "SBS23",
      "SBS24",
      "SBS25",
      "SBS26",
      "SBS27",
      "SBS28",
      "SBS29",
      "SBS30",
      "SBS31",
      "SBS32",
      "SBS33",
      "SBS34",
      "SBS35",
      "SBS36",
      "SBS37",
      "SBS38",
      "SBS39",
      "SBS40",
      "SBS41",
      "SBS42",
      "SBS43",
      "SBS44",
      "SBS45",
      "SBS46",
      "SBS47",
      "SBS48",
      "SBS49",
      "SBS50",
      "SBS51",
      "SBS52",
      "SBS53",
      "SBS54",
      "SBS55",
      "SBS56",
      "SBS57",
      "SBS58",
      "SBS59",
      "SBS60",
      "SBS-others"
    ),
    Subsname=c(
      "Signature Subs-01",
      "Signature Subs-02",
      "Signature Subs-03",
      "Signature Subs-04",
      "Signature Subs-05",
      "Signature Subs-06",
      "Signature Subs-07a",
      "Signature Subs-07b",
      "Signature Subs-07c",
      "Signature Subs-07d",
      "Signature Subs-08",
      "Signature Subs-09",
      "Signature Subs-10a",
      "Signature Subs-10b",
      "Signature Subs-11",
      "Signature Subs-12",
      "Signature Subs-13",
      "Signature Subs-14",
      "Signature Subs-15",
      "Signature Subs-16",
      "Signature Subs-17a",
      "Signature Subs-17b",
      "Signature Subs-18",
      "Signature Subs-19",
      "Signature Subs-20",
      "Signature Subs-21",
      "Signature Subs-22",
      "Signature Subs-23",
      "Signature Subs-24",
      "Signature Subs-25",
      "Signature Subs-26",
      "Signature Subs-27",
      "Signature Subs-28",
      "Signature Subs-29",
      "Signature Subs-30",
      "Signature Subs-31",
      "Signature Subs-32",
      "Signature Subs-33",
      "Signature Subs-34",
      "Signature Subs-35",
      "Signature Subs-36",
      "Signature Subs-37",
      "Signature Subs-38",
      "Signature Subs-39",
      "Signature Subs-40",
      "Signature Subs-41",
      "Signature Subs-42",
      "Signature Subs-43",
      "Signature Subs-44",
      "Signature Subs-45",
      "Signature Subs-46",
      "Signature Subs-47",
      "Signature Subs-48",
      "Signature Subs-49",
      "Signature Subs-50",
      "Signature Subs-51",
      "Signature Subs-52",
      "Signature Subs-53",
      "Signature Subs-54",
      "Signature Subs-55",
      "Signature Subs-56",
      "Signature Subs-57",
      "Signature Subs-58",
      "Signature Subs-59",
      "Signature Subs-60",
      "Signature Subs-others"
      
    )
  )



SBS2Subs <- signames$Subsname
names(SBS2Subs) <- signames$SBSname
Subs2SBS <- signames$SBSname
names(Subs2SBS) <- signames$Subsname


### define the signature funciton ####
Subscolor <- c(
  'Signature Subs-01'='#4a9855',
  'Signature Subs-02'='#e2a8ab',
  'Signature Subs-03'='#40004b',
  'Signature Subs-04'='#5aa1ca',
  'Signature Subs-05'='#305d39',
  'Signature Subs-06'='#785940',
  "Signature Subs-07a"='#6e70b7',
  "Signature Subs-07b"='#ff7f00',
  "Signature Subs-07c"='#fec44f',
  "Signature Subs-07d"='#846a2a',
  "Signature Subs-08"='#cab2d6',
  "Signature Subs-09"='#f4a582',
  "Signature Subs-10a"='#8dd3c7',
  "Signature Subs-10b"='#5e4fa2',
  'Signature Subs-12'='#ffed6f',
  'Signature Subs-13'='#e41a1c',
  'Signature Subs-14'='#ffffbf',
  'Signature Subs-15'='#4d4d4d',
  'Signature Subs-16'='#513276',
  'Signature Subs-17a'='#df4c7d',
  'Signature Subs-17b'='#08519c',
  'Signature Subs-18'='#b3de69',
  'Signature Subs-19'='#4d4d4d',
  'Signature Subs-20'='#b2182b',
  'Signature Subs-22'='#d9ef8b',
  'Signature Subs-21'='#e6f5d0',
  'Signature Subs-24'='#1c9099',
  'Signature Subs-25'='#35978f',
  'Signature Subs-26'='#5e4fa2',
  'Signature Subs-28'='#de77ae',
  'Signature Subs-30'='#5e4fa2',
  'Signature Subs-31'='#f781bf',
  'Signature Subs-32'='#dd1c77',
  'Signature Subs-33'='#b25d7e',
  'Signature Subs-35'='#fc8d59',
  'Signature Subs-36'='yellow',
  'Signature Subs-39'='#636363',
  'Signature Subs-40'='#b15928',
  'Signature Subs-41'='#fccde5',
  'Signature Subs-44'='#b3de69',
  'Signature Subs-46'='#e6f598',
  'Signature Subs-47'='#bababa',
  'Signature Subs-42'='#ae017e',
  'Signature Subs-54'='#fcc5c0',
  'Signature Subs-56'='#8c510a',
  'Signature Subs-others'='#cececa'
)


tmpcolor <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2','#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b')
names(tmpcolor) <- c('Signature Subs-11','Signature Subs-23','Signature Subs-27','Signature Subs-29','Signature Subs-34','Signature Subs-37','Signature Subs-38','Signature Subs-43','Signature Subs-45','Signature Subs-48','Signature Subs-49','Signature Subs-50','Signature Subs-51','Signature Subs-52','Signature Subs-53','Signature Subs-55','Signature Subs-57','Signature Subs-58','Signature Subs-59','Signature Subs-60')

Subscolor <- c(Subscolor,tmpcolor)

SBScolor <- Subscolor
names(SBScolor) <- Subs2SBS[names(SBScolor)]


### color 12 for the clustering
color12 <- rev(c(
  "#a6cee3",
  "#1f78b4",
  "#b2df8a",
  "#33a02c",
  "#fb9a99",
  "#e31a1c",
  "#fdbf6f",
  "#ff7f00",
  "#cab2d6",
  "#6a3d9a",
  "#b15928"
))


### SBS cluster funciton ####
# p4 proportion bar

SBS96_Clustering <- function(sigdata,studydata=NULL,cosinedata,cosine_cutoff=0,sigcolor=NULL,studycolor=NULL,puritydata=NULL,puritydata_cat=FALSE, puritycol=NULL,purity_cutoff=0,clustern,highlight=NULL,legendnrow=3,sampletext_size=6, filename='SBS96_clustering.pdf',height=10,width=20,hc_func='hclust',hc_metric = 'euclidean',hc_method = 'ward.D2',stand=TRUE){
  
  require(tidyverse)
  require(scales)
  require(janitor)
  require(ggsci)
  require(ggpubr)
  require(factoextra)
  require(cowplot)
  
  blankcol <- NA
  names(blankcol) <- "        "
  
  if(is.null(sigcolor)){
    sigcolorindex <- as.character(Subscolor[colnames(sigdata)[-1]])
    names(sigcolorindex) <- colnames(sigdata)[-1]
  }else{
    sigcolorindex <- sigcolor
  }
  
  # colorall <- c(sigcolorindex,"white",studycolor)
  # names(colorall) <- c(names(sigcolorindex),"        ",names(studycolor))
  # 
  colorall <- c(sigcolorindex,blankcol,studycolor)
  
  tmp=sigdata %>% adorn_percentages('row') 
  mdata=as.matrix(tmp[,-1])
  rownames(mdata) <- tmp$Samples
  
  #fviz_nbclust(mdata, kmeans, method = "gap_stat")
  kcolors <- pal_d3("category20")(clustern)
  
  res <- hcut(mdata,k = clustern,hc_func = hc_func,hc_metric = hc_metric,hc_method = hc_method,stand=stand)
  p1 <- fviz_dend(res, rect = TRUE, cex = 0.5,k_colors = kcolors,lwd = 0.5,show_labels = F)+scale_x_discrete(expand = c(0,0))+ theme(plot.margin=margin(b=-0.7,unit="cm"),title = element_blank())
  #plot.margin=margin(b=-1,unit="cm")
  
  if(!is.null(studydata)){
    studydata <- studydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order]))
    p2 <- studydata %>% ggplot(aes(Samples,1,fill=factor(Study,levels = names(studycolor))))+geom_tile(col="black")+scale_fill_manual(values =studycolor,drop=FALSE)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.4,unit="cm"),title = element_blank())+ylim(c(0,2))
  }
  p2.1 <- cosinedata %>%  mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% mutate(Similarity=if_else(Similarity<cosine_cutoff,NA_real_,Similarity)) %>% ggplot(aes(Samples,1,fill=Similarity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'C',limits = c(0.6, 1), oob = scales::squish)+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(3, "cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(title = element_blank())+ylim(c(0,2))
  p2.2 <- as_ggplot(get_legend(p2.1))
  p2.2 <- p2.2+theme(plot.margin = margin(b = 30))
  p2.1 <- cosinedata %>%  mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% mutate(Similarity=if_else(Similarity<cosine_cutoff,NA_real_,Similarity)) %>% ggplot(aes(Samples,1,fill=Similarity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'C',limits = c(0.6, 1), oob = scales::squish)+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t = 0.5,unit="cm"),title = element_blank())+ylim(c(0,2))
  
  p2.5 <- p2.2
  
  if(!is.null(puritydata)){
    puritydata <- puritydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) 
    
    if(!puritydata_cat) {
      puritydata <- puritydata %>% mutate(Purity=if_else(Purity<purity_cutoff,NA_real_,Purity))
      p2.3 <- puritydata  %>% ggplot(aes(Samples,1,fill=Purity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(3, "cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(title = element_blank())+ylim(c(0,2))
      p2.4 <- as_ggplot(get_legend(p2.3))
      p2.4 <- p2.4+theme(plot.margin = margin(b = 30))
      p2.3 <- puritydata %>% ggplot(aes(Samples,1,fill=Purity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t = 0.5,unit="cm"),title = element_blank())+ylim(c(0,2))
      #
      p2.5 <- plot_grid(p2.2, p2.4, align = "h", axis = "b", rel_widths = c(1, 1))
      
    }else { 
      names(blankcol) <- "    "
      colorall <- c(colorall,blankcol,puritycol)
      p2.3 <- puritydata %>% ggplot(aes(Samples,1,fill=factor(Purity,levels = names(puritycol))))+geom_tile(col="black")+scale_fill_manual(values =puritycol)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"),title = element_blank())+ylim(c(0,2))
      
    }
  }
  
  
  p3 <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,(Weight),fill=factor(Signature,levels = names(colorall))))+geom_bar(stat="identity",col="gray95",width=1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = "xy")+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),panel.grid.major.x=element_blank(),legend.position = "none",legend.box.spacing = unit(0,"cm"),plot.margin=margin(b=-1,t = 1,unit="cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks(),labels = comma)+scale_fill_manual(values = colorall,drop=FALSE)+guides(fill=guide_legend(nrow=2,byrow=TRUE))+xlab("")+ylab("Number of mutations \n")
  #+theme(plot.margin=margin(b=4,unit="pt"))
  
  p4 <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,Weight,fill=factor(Signature,levels = names(colorall))))+geom_bar(stat="identity",position="fill",col="gray95",width = 1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12)+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_text(size = sampletext_size, angle = 90, vjust = 0.5, hjust = 1),panel.grid.major=element_line(),legend.position = "bottom",legend.box.background = element_blank(),legend.box.spacing = unit(-0.5,"cm"),legend.key = element_rect(size = 0),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.25, "cm"),legend.key.width =unit(1, "cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks())+xlab("")+scale_fill_manual(values = colorall,drop=FALSE)+guides(fill=guide_legend(nrow=legendnrow,byrow=TRUE,label.position = "bottom"))+ylab("Signature contribution\n")
  #+theme(plot.margin=margin(t=4,unit="pt"))
  
  
  if(!is.null(highlight)){
    samhigh <- sigdata %>% mutate(Samples_high=if_else(Samples %in% highlight,paste0("*",Samples),Samples))
    p4 <-  p4+scale_x_discrete(breaks=samhigh$Samples,labels=samhigh$Samples_high)
  }
  
  #p3 <- flush_ticks(p3,flush = "Y",plot = FALSE)
  #p4 <- flush_ticks(p4,flush = "Y",plot = FALSE)
  
  
  if(!is.null(studydata)){
    if(!is.null(puritydata)){
      pall <- ggarrange(p1,p2,p2.3,p3,p2.1,p4,p2.5,ncol=1,nrow = 7,align = 'v',heights = c(2,0.1,0.1,4,0.1,7,0.1))
    }else {
      pall <- ggarrange(p1,p2,p3,p2.1,p4,p2.5,ncol=1,nrow = 6,align = 'v',heights = c(2,0.1,4,0.1,7,0.1))
    }
  }else{
    if(!is.null(puritydata)){
      pall <- ggarrange(p1,p2.3,p3,p2.1,p4,p2.5,ncol=1,nrow = 6,align = 'v',heights = c(2,0.1,4,0.1,7,0.1))
    }
  }
  
  ggsave(filename = filename,plot = pall,height = height,width = width,device = cairo_pdf)
}




### calculate_similarities ####

# cos_sim <- function(a, b) {
#   if (sum(a) == 0 | sum(b) == 0) {
#     value = 0
#   }
#   dot_product = a %*% b
#   norm_a <- norm(as.matrix(a), "2")
#   norm_b <- norm(as.matrix(b), "2")
#   value = dot_product / (norm_a * norm_b)
#   return(value)
# }

## from MutationalPattern
cos_sim <- function (x, y) 
{
  res = x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
  res = as.numeric(res)
  return(res)
}

# add the Correlation
calculate_similarities <- function(orignal_genomes, signature, signature_activaties) {
  require(entropy)
  data2 <- as.data.frame(orignal_genomes)
  #data2 <- read.delim(orignal_genomes_file,header = T,check.names = F,stringsAsFactors = F)
  #data2 <- data.frame(data2)
  data3 <-  as.data.frame(signature)
  data4 <-  as.data.frame(signature_activaties)
  #colnames(data2[,2:ncol(data2)])
  
  if(dim(data2)[1]==0 | dim(data3)[1]==0 | dim(data4)[1]==0){
    return(NA)
  }else {
    data2 <- data2[,!is.na(colSums(data2 != 0)) & colSums(data2 != 0) > 0]
    ## filter out the samples and signatures
    sample_name_total <- data4$Sample
    sample_name_tmp<- colnames(data2)[-1]
    sample_name <- sample_name_total[sample_name_total %in% sample_name_tmp]
    
    signature_name_total <- colnames(data4)[-1]
    signature_name_tmp <- colnames(data3)[-1]
    signature_name <- signature_name_total[signature_name_total %in% signature_name_tmp]
    
    data2 <- data2 %>% select(MutationType,one_of(sample_name))
    data3 <- data3 %>% select(MutationType,one_of(signature_name))
    data4 <- data4 %>% select(Sample,one_of(signature_name)) %>% filter(Sample %in% sample_name)
    
    genomes <- data2[, 2:length(data2)]
    
    #data3 <- read.delim(signature_file,header = T,check.names = F,stringsAsFactors = F)
    
    data3 <- data3[,2:ncol(data3)]
    #data4 <- read.delim(signature_activaties_file,header = T,check.names = F,stringsAsFactors = F)
    
    data4 <- data4[,2:ncol(data4)]
    est_genomes <- as.data.frame(as.matrix(data3) %*% as.matrix(t(data4)))
    
    cosine_sim_list = c()
    correlation_list = c()
    kl_divergence_list = c()
    l1_norm_list = c()
    l2_norm_list = c()
    total_mutations = c()
    relative_l1_list = c()
    relative_l2_list = c()
    for (i in 1:ncol(genomes)) {
      p_i <- as.numeric(genomes[, i])
      q_i = (est_genomes[, i])
      if(sum(q_i)>0){
        cosine_sim_list = append(cosine_sim_list, round(cos_sim(p_i, q_i), digits=3))
        correlation_list = append(correlation_list, round(cor(p_i, q_i), digits=3))
        kl_divergence_list = append(kl_divergence_list, round(KL.empirical(p_i, q_i), digits=4))
        l1_norm_list = append(l1_norm_list, round(norm(as.matrix(p_i-q_i), "1"), digits=3))
        relative_l1_list = append(relative_l1_list, round((dplyr::last(l1_norm_list)/norm(as.matrix(p_i), "1"))*100, digits=3))
        l2_norm_list = append(l2_norm_list, round(norm(as.matrix(p_i-q_i), "2"), digits=3))
        relative_l2_list = append(relative_l2_list, round((dplyr::last(l2_norm_list)/norm(as.matrix(p_i), "2"))*100, digits=3))
        total_mutations = append(total_mutations, sum(p_i))
      }else{
        cosine_sim_list = append(cosine_sim_list, NA_real_)
        correlation_list = append(correlation_list, NA_real_)
        kl_divergence_list = append(kl_divergence_list,NA_real_)
        l1_norm_list = append(l1_norm_list, NA_real_)
        relative_l1_list = append(relative_l1_list, NA_real_)
        l2_norm_list = append(l2_norm_list, NA_real_)
        relative_l2_list = append(relative_l2_list, NA_real_)
        total_mutations = append(total_mutations, sum(p_i))
      }
      
      
    }
    kl_divergence_list[!is.na(kl_divergence_list) & !is.finite(kl_divergence_list)] <- 1000
    similarities_dataframe = data.frame("Sample_Names"=sample_name,
                                        "Total_Mutations"=total_mutations,
                                        "Cosine_similarity"=cosine_sim_list,
                                        "L1_Norm"=l1_norm_list,
                                        `L1_Norm_%`=relative_l1_list,
                                        `100-L1_Norm_%`=100-relative_l1_list,
                                        "L2_Norm"=l2_norm_list,
                                        `L2_Norm_%`=relative_l2_list,
                                        `100-L2_Norm_%`=100-relative_l2_list,
                                        "KL_Divergence"= kl_divergence_list,
                                        "Correlation" = correlation_list,
                                        check.names = F)
    #write.csv(similarities_dataframe, file="MyData.csv")
    #write.table(similarities_dataframe, file="MyData.txt", sep="\t", row.names = FALSE)
    
    similarities_dataframe <- tibble("Sample_Names"=sample_name_total) %>% left_join(similarities_dataframe)
    return(similarities_dataframe)
  }
}

#calculate_similarities("Test/orignal_genomes.txt", "Test/Wsignature_sigs.txt", "Test/Wsignature_activaties.txt")
L1_normal_displot <- function(similarities_dataframe,...){
  require(ggplot2)
  require(hrbrthemes)
  similarities_dataframe %>% ggplot(aes(100-`L1_Norm_%`))+geom_histogram(color="white",binwidth = 2)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 16,axis = "xy")+scale_x_continuous(name ="\n% mutations explained by extracted signatures (100-L1_Norm %)",expand = c(0,0),limits = c(0,100))+scale_y_continuous(name = "Number of sample",expand = c(0,0))
}




#Generation of probabilities for each processes given to A mutation type 
#probabilities('Decomposed_Solution_Signatures.txt','Decomposed_Solution_Activities.txt')
probabilities <-  function(W, H){ 
  require(dplyr)
  #W=read.delim(Wfile,header = T,check.names = F,stringsAsFactors = F)
  #H=read.delim(Hfile,header = T,check.names = F,stringsAsFactors = F)
  
  MutationType <- W$MutationType
  allcolnames <- H$Samples
  
  W <- as.matrix(W[,2:ncol(W)])
  H <- t(as.matrix(H[,2:ncol(H)]))
  genomes <- W %*% H
  
  probs_all <- NULL
  
  for(i in 1:dim(H)[2]){
    probs <-(W*(H[,i])[col(W)])/(genomes[,i])[row(W)]
    rownames(probs) <- MutationType
    probs <- as.data.frame(probs) %>% rownames_to_column(var = "MutationType") %>% mutate(Sample_Names=allcolnames[i]) %>% select(Sample_Names,MutationType,everything())
    probs_all <- bind_rows(probs_all,probs)
  }
  return(probs_all)
  
}



# SBS96_plot --------------------------------------------------------------

plot_sbs_96_profile <- function(data,samplename=NULL,totalmut=NULL,samplename_plot=TRUE, totalmut_plot=TRUE,filename=NULL,percentage=TRUE,ytitle="Percentage of Single Base Subsitutions"){
  
  require(tidyverse)
  require(hrbrthemes)
  require(scales)
  require(ggpubr)  
  
  #### data is a frame with the following columns: Type, SubType, MutationType, Value ###
  
  if(dim(data)[2]==2){
    data <- data %>% mutate(Type=str_sub(MutationType,3,5),SubType=paste0(str_sub(MutationType,1,1),str_sub(MutationType,3,3),str_sub(MutationType,7,7))) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  }
  
  if(is.null(samplename)){
    samplename <- colnames(data)[4]
  }
  
  colnames(data)[4] <- "Value"
  
  data <- data %>% mutate(Seq=seq_along(MutationType)) %>% mutate(Type=fct_inorder(Type),SubType=fct_inorder(SubType),MutationType=fct_inorder(MutationType))
  
  stype <- unique(data$Type)
  
  if(!is.null(totalmut)){
    totalmut <- comma_format()(totalmut)
  }
  
  COLORS6 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE")
  names(COLORS6) <- stype
  
  data0 <- bind_cols(
    data %>% group_by(Type) %>% dplyr::slice(1) %>% ungroup() %>% mutate(Value1=Seq+0.25) %>% select(Value1,Type),
    data %>% group_by(Type) %>% dplyr::slice(16) %>% ungroup() %>% mutate(Value2=Seq+0.25) %>% select(Value2)
  )
  
  if(percentage){
    
    if(sum(data$Value)>10){
      totalmut <- sum(data$Value)
      data$Value <- data$Value/totalmut
    }
    ymax <- max(data$Value)/0.9
    ymax <- 0.04*ceiling(ymax/0.04)
    labelx <- percent
    #ylabp <- percent(pretty_breaks(n = 5)(ymax))
    ylabp <- percent(seq(0,ymax,length.out = 5))
    ytext <- ytitle
    
  }else{
    ymax <- max(data$Value)/0.9
    ymax <- ceiling(ymax)
    labelx <- comma
    #ylabp <- comma(pretty_breaks(n = 5)(ymax))
    ylabp <- comma(seq(0,ymax,length.out = 5))
    ytext <- "Number of Single Base Subsitutions"
  }
  
  
  
  #ylabp <- percent(pretty_breaks(n = 5)(data$Value))
  ylabp <- ylabp[length(ylabp)]
  
  data00 <- data %>%  mutate(B1=str_sub(SubType,1,1),B2=str_sub(SubType,2,2),B3=str_sub(SubType,3,3))
  
  p1 <- data0 %>% ggplot(aes(xmin=Value1,xmax=Value2,ymin=0,ymax=1,fill=Type))+
    geom_rect()+
    scale_fill_manual(values=COLORS6)+
    scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits=c(1,96.5))+
    scale_y_continuous(expand = c(0,0),labels =ylabp,breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color = "white"),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          #axis.line = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(b=-1.2,l = 0.2,unit="cm")
    )+
    annotate("text", x = 8.75+16*0, y = 1.5, label = stype[1],size=5,fontface =2)+
    annotate("text", x = 8.75+16*1, y = 1.5, label = stype[2],size=5,fontface =2)+
    annotate("text", x = 8.75+16*2, y = 1.5, label = stype[3],size=5,fontface =2)+
    annotate("text", x = 8.75+16*3, y = 1.5, label = stype[4],size=5,fontface =2)+
    annotate("text", x = 8.75+16*4, y = 1.5, label = stype[5],size=5,fontface =2)+
    annotate("text", x = 8.75+16*5, y = 1.5, label = stype[6],size=5,fontface =2)
  
  
  
  p2 <- data %>% 
    ggplot(aes(Seq,Value,fill=Type))+
    geom_col(color="white",width = 0.5,size=0)+
    scale_fill_manual(values=COLORS6)+
    scale_x_continuous(expand = c(0,0),labels = data$SubType,breaks = data$Seq,sec.axis = dup_axis(labels = NULL,name = NULL))+
    scale_y_continuous(expand = c(0,0),labels = labelx,breaks = seq(0,ymax,length.out = 5),limits = c(0,ymax),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=ytext)+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'gray80',size=0.6),
          axis.line.y = element_line(colour = 'gray80',size=0.6),
          plot.margin=margin(t=-1.2,b=-1,l=0.2,unit="cm")
          
    )
  
  if(samplename_plot & !is.null(samplename)){
    p2 <- p2+annotate("text", x = 1.5, y = ymax*0.9, label = samplename,size=7,fontface =2,hjust = 0)
  }
  
  if(totalmut_plot & !is.null(totalmut)){
    p2 <- p2+annotate("text", x = 78, y = ymax*0.88, label = paste0("Total Mutations: ",totalmut),size=6.5,fontface =2,hjust = 0)
  }
  
  invisible(capture.output(p2 = flush_ticks(p2)+theme(axis.text.x = element_blank())))
  
  p3 <- data00 %>% 
    ggplot(aes(Seq))+
    geom_text(aes(y=0.4,label=B1),col="gray60",angle=90,hjust = 0,vjust = 0.5,size=3.5)+
    geom_text(aes(y=1,label=B2,col=Type),angle=90,hjust = 0,vjust = 0.5,size=3.5,family = "Arial")+
    geom_text(aes(y=1.6,label=B3),col="gray60",angle=90,hjust = 0,vjust = 0.5,size=3.5)+
    scale_color_manual(values=COLORS6)+
    scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits = c(0.5,96.5))+
    scale_y_continuous(expand = c(0,0),labels =ylabp,breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,colour = "white"),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(t=-1,l=0.2,unit="cm")
    )
  #,text = element_text(family = "Arial")
  
  
  #cairo_pdf(file = 'tmp.pdf',width = 16,height = 4)
  pcomb <- ggarrange(p1,p2,p3,nrow = 3,align = "h",heights = c(1.8,10,1.5))
  #pcomb <- plot_grid(p1,p2,p3,nrow = 3,align = "v",axis = c("tb"),rel_heights = c(1.8,10,1.5))
  
  if(is.null(filename)){
    return(pcomb)
  } else {
    ggsave(filename,width = 18,height = 4,device = cairo_pdf)
  }
  
}




# PLOT_DBS_78 -------------------------------------------------------------

plot_dbs_78_profile <- function(data,samplename=NULL,samplename_plot=TRUE, filename=NULL){
  require(tidyverse)
  require(hrbrthemes)
  require(scales)
  require(cowplot)  
  
  data <- data %>% mutate(Type=paste0(str_sub(MutationType,1,3),"NN"),SubType=str_sub(MutationType,4,5)) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  
  if(is.null(samplename)){
    samplename <- colnames(data)[4]
  }
  
  colnames(data)[4] <- "Value"
  data <- data %>% mutate(Seq=seq_along(MutationType)) %>% mutate(Type=fct_inorder(Type),SubType=fct_inorder(SubType),MutationType=fct_inorder(MutationType))
  stype <- unique(data$Type)
  
  totalmut <- sum(data$Value)
  totalmut <- comma_format()(totalmut)
  
  COLORS10 = c("#03BCEE", "#0366CB", "#A1CE63", "#016601", "#FE9898", "#E32926", "#FEB166", "#FE8001", "#CB98FE", "#4C0198")
  names(COLORS10) <- stype
  
  ymax <- max(data$Value)/0.9
  ymax <- 4*ceiling(ymax/4)
  
  #ylabp <- pretty_breaks(n = 5)(data$Value)
  ylabp <- comma(seq(0,ymax,length.out = 5))
  ylabp <- ylabp[length(ylabp)]
  
  data0 <- bind_cols(
    data %>% group_by(Type) %>% filter(row_number() == 1) %>% ungroup() %>% mutate(Value1=Seq+0.25) %>% select(Value1,Type),
    data %>% group_by(Type) %>% "["(.,c(which(data$Type != lag(data$Type))-1, 78),) %>% ungroup() %>% mutate(Value2=Seq+0.25) %>% select(Value2)
  )
  
  #ylabp <- pretty_breaks(n = 5)(data$Value)
  #ylabp <- ylabp[length(ylabp)]
  
  anno.x=(data0$Value1+data0$Value2)/2
  anno.lab=data0$Type
  
  p1 <- 
    data0 %>% ggplot(aes(xmin=Value1,xmax=Value2,ymin=0,ymax=0.9,fill=Type))+
    geom_rect()+
    scale_fill_manual(values=COLORS10)+
    scale_x_continuous(expand = c(0,0),limits=c(1,78.5))+
    scale_y_continuous(expand = c(0,0),breaks = 1,labels = ylabp,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(y=" ",x=NULL)+annotate("text",x=anno.x,y=1.5,label=anno.lab,size=5,fontface =2,hjust=0.5)+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color = "white"),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          #axis.line.x = element_line(colour = 'white',size=0.6),
          #axis.line.y = element_line(colour = 'white',size=0.6),
          axis.line = element_blank(),
          plot.margin=margin(b=0,l = 0.2,unit="cm"),
    )
  #+theme(plot.background = element_rect(fill = "darkblue"))
  
  p2 <- data %>% 
    ggplot(aes(Seq,Value,fill=Type))+
    geom_col(color="white",width = 0.4,size=0)+
    scale_fill_manual(values=COLORS10)+
    scale_x_continuous(expand = c(0,0),labels = data$SubType,breaks = data$Seq,sec.axis = dup_axis(labels = NULL,name = NULL))+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,ymax,length.out = 5),limits = c(0,ymax),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y="Number of Double Base Subsitutions")+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=11,angle = 90,hjust = 1,vjust = 0.5,colour = COLORS10[data$Type]),
          #axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'gray80',size=0.6),
          axis.line.y = element_line(colour = 'gray80',size=0.6),
          plot.margin=margin(t=0,l=0.2,unit="cm")
    )
  
  if(samplename_plot & !is.null(samplename)){
    p2 <- p2+annotate("text", x = 1, y = ymax*0.9, label = paste0(samplename,": ",totalmut," double subs"),size=7,fontface =2,hjust = 0)
  }
  
  #invisible(capture.output(p2 <- flush_ticks(p2,flush = "X")))
  p2 <- p2+ theme(axis.text.y=element_text(vjust=c(0, rep(0.5, 3), 1)))
  
  
  #cairo_pdf(file = 'tmp.pdf',width = 16,height = 4)
  #pcomb <- ggarrange(p1,p2,nrow = 2,align = "h",heights = c(3,15))
  #pcomb <- grid.arrange(p1,p2,nrow=2)
  pcomb <- plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(1/6,5/6))
  
  if(is.null(filename)){
    return(pcomb)
  } else {
    ggsave(filename,pcomb,width = 18,height = 4,device = cairo_pdf)
  }
  
}
#plot_dbs_78_profile(data = data_dbs_exm1,filename = "tmp.pdf")


#save(data_sbs_exm1,data_dbs_exm1,file='../Example.RData')







# Heatmap of profiles -----------------------------------------------------

profile_heatmap_plot <- function(data,output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  # data format: Sample  Profile Mutations
  
  profile_tmp <- data %>% group_by(Profile) %>% summarise(mean=mean(Mutations)) %>% arrange(desc(mean))
  samleve_tmp <- data %>% filter(Profile==profile_tmp$Profile[1]) %>% arrange(Mutations) %>% pull(Sample)
  
  # typedata_tmp <- seqmatrix_refdata_public %>% group_by(Profile,MutationType) %>% summarise(mean=mean(Mutations,na.rm = TRUE)) %>% ungroup() %>% arrange(desc(mean)) %>% group_by(Profile) %>% slice(1:5) %>% mutate(Seq2=6-seq_along(MutationType)) %>% ungroup() %>% left_join(profile_tmp %>% select(-mean)) %>% mutate(Seq3=Seq+0.1666667*Seq2)
  # 
  # data2 <- seqmatrix_refdata_public %>% 
  #   left_join(typedata_tmp %>% select(Profile,MutationType,Seq,Seq2)) %>% 
  #   filter(!is.na(Seq)) %>% 
  #   left_join(
  #     seqmatrix_refdata_public %>% group_by(Sample,Profile) %>% summarise(total=sum(Mutations))
  #   ) %>% 
  #   mutate(Ratio=Mutations/total)
  # 
  # max_tmp <- data2 %>% group_by(Profile) %>% summarise(max=max(Ratio))
  # 
  # data2 <- data2 %>% left_join(max_tmp) %>% mutate(value=Ratio/max+Seq)
  name_max <- max(str_length(unique(data$Sample)))
  name_len <- length(unique(data$Sample))
  angle_max <- if_else(name_max < 15, 90, 30)
  vjust_max <- if_else(name_max < 15, 0.5, 1)
  
  plot_final <- data %>% left_join(profile_tmp) %>% 
    mutate(Sample=factor(Sample,levels = samleve_tmp),Profile=factor(Profile,levels = profile_tmp$Profile)) %>% 
    ggplot(aes(Sample,Profile,fill=(Mutations)))+
    geom_tile(col="white")+
    scale_fill_viridis_c(trans = "log10",label=comma_format(),na.value = 'black')+
    labs(x="",y="",fill="Number of mutations\n")+
    theme_ipsum_rc(grid = FALSE,ticks = FALSE,axis = FALSE)+
    theme(legend.key.width =unit(2, "cm"),legend.position = "top")
  
  if(name_len<=80){
    plot_final <- plot_final + theme(axis.text.x = element_text(angle = angle_max,hjust = 1,vjust = vjust_max,size = 10))
  } else{
    plot_final <- plot_final + theme(axis.text.x = element_blank())
  }
  
  ## define the length of x and y
  leng0 <- 2.5
  leng_ratio <-  0.2
  
  xleng <- leng_ratio*(name_len)+leng0+2.5
  
  xleng <- if_else(xleng>15,15,xleng)
  yleng <- leng_ratio*name_max+2.5+leng0
  
  if(is.null(output_plot)){
    return(plot_final)
  }else{
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot_final,width = plot_width,height = plot_height)
  }
  
}



# Mutational Patten funcitons ---------------------------------------------

profile_format_df <- function(data,factortype=FALSE,indel_short=FALSE){
  # format, sorting and factor the signature dataframe for SBS96,DBS78 and ID83
  
  if(dim(data)[1]==96){
    data <- data %>% mutate(Type=str_sub(MutationType,3,5),SubType=paste0(str_sub(MutationType,1,1),str_sub(MutationType,3,3),str_sub(MutationType,7,7))) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  }
  
  if(dim(data)[1]==78){
    data <- data %>% mutate(Type=paste0(str_sub(MutationType,1,3),"NN"),SubType=str_sub(MutationType,4,5)) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  }
  
  if(dim(data)[1]==83){
    idtypeorder <- c("1:Del:C","1:Del:T","1:Ins:C","1:Ins:T","2:Del:R","3:Del:R","4:Del:R","5:Del:R","2:Ins:R","3:Ins:R","4:Ins:R","5:Ins:R","2:Del:M","3:Del:M","4:Del:M","5:Del:M")
    data <- data %>% mutate(Type=str_sub(MutationType,1,7),SubType=str_sub(MutationType,9,9)) %>% select(Type,SubType,MutationType,everything()) %>% mutate(Type=factor(Type,levels = idtypeorder)) %>% arrange(Type,SubType) %>% mutate(Type=as.character(Type))
  }
  
  if(factortype){
    
    if(dim(data)[1]==83){
      tmplev <- idtypeorder
      tmplab <- str_replace(tmplev,"5","5+")
      
      if(indel_short){
        tmplab = if_else(tmplab %in% c('2:Del:M','3:Del:M','4:Del:M'),str_remove(tmplab,":.*"),tmplab)
      }
      
      data <- data %>% 
        mutate(Type=factor(Type,levels = tmplev,labels = tmplab)) %>% 
        mutate(SubType = if_else(str_detect(Type,'Del') & !str_detect(Type,"M"),as.character(as.integer(SubType)+1L),SubType)) %>% 
        mutate(SubType = if_else(str_detect(Type,'Del') & !str_detect(Type,"M"), str_replace(SubType,"6","6+"),str_replace(SubType,"5","5+")))
    }else {
      tmplev <- unique(data$Type)
      data <- data %>% mutate(Type=factor(Type,levels = tmplev))
    }
    
    
    #tmplev <- unique(data$SubType)
    #data <- data %>% mutate(SubType=factor(SubType,levels = tmplev))
    tmplev <- unique(data$MutationType)
    data <- data %>% mutate(MutationType=factor(MutationType,levels = tmplev))
  }
  
  return(data)
  
}


## need to define profile_format_df2 ## 






# Calculate cosine similarity between two signature in dataframe format --------
cos_sim_df_old <- function (mut_df1, mut_df2, output_matrix=FALSE) 
{
  colnames(mut_df1)[1] <- "MutationType"
  colnames(mut_df2)[1] <- "MutationType"
  
  mut_matrix1 <- mut_df1 %>% 
    arrange(MutationType) %>% 
    select(-MutationType) %>% 
    as.matrix()
  
  mut_matrix2 <- mut_df2 %>% 
    arrange(MutationType) %>% 
    select(-MutationType) %>% 
    as.matrix()
  
  n_samples1 = ncol(mut_matrix1)
  n_samples2 = ncol(mut_matrix2)
  res_matrix = matrix(nrow = n_samples1, ncol = n_samples2)
  for (s in 1:n_samples1) {
    signal1 = mut_matrix1[, s]
    cos_sim_vector = c()
    for (i in 1:n_samples2) {
      signal2 = mut_matrix2[, i]
      cos_sim_vector[i] = cos_sim(signal1, signal2)
    }
    res_matrix[s, ] = cos_sim_vector
  }
  rownames(res_matrix) = colnames(mut_matrix1)
  colnames(res_matrix) = colnames(mut_matrix2)
  
  if(output_matrix){
    return(res_matrix)
  }else {
    res_matrix %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
  }
}


# Calculate cosine similarity between two signature in dataframe format using coop package --------
cos_sim_df <- function (mut_df1, mut_df2, output_matrix=FALSE) 
{
  require(coop)
  
  colnames(mut_df1)[1] <- "MutationType"
  colnames(mut_df2)[1] <- "MutationType"
  
  if(identical(mut_df1,mut_df2)){
    # two identical matrix/df
    mut_matrix1 <- mut_df1 %>% 
      arrange(MutationType) %>% 
      select(-MutationType) %>% 
      as.matrix()
    res_matrix <- coop::cosine(mut_matrix1) 
    
  }else{
    # two differnt matrix/df
    df_name1 <- colnames(mut_df1)[-1]
    df_name2 <- colnames(mut_df2)[-1]
    new_df_name1 <- paste0("N1_",seq_along(df_name1))
    new_df_name2 <- paste0("N2_",seq_along(df_name2))
    colnames(mut_df1)[-1] <- new_df_name1
    colnames(mut_df2)[-1] <- new_df_name2  
    mut_df <- left_join(mut_df1,mut_df2)
    mut_matrix <- mut_df %>% 
      arrange(MutationType) %>% 
      select(-MutationType) %>% 
      as.matrix()
    
    res_matrix <- coop::cosine(mut_matrix) 
    res_matrix <- res_matrix[new_df_name1,new_df_name2]
    rownames(res_matrix) <- df_name1
    colnames(res_matrix) <- df_name2
  }
  
  if(output_matrix){
    return(res_matrix)
  }else {
    res_matrix %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
  }
  
}



# Heatmap of cosine similairty  -------------------------------------------
plot_cosine_heatmap_df <- function (cos_sim_df, col_order, cluster_rows = TRUE, method = "complete", nmax = 200L, plot_values = FALSE,output_plot = NULL,plot_width=NULL, plot_height=NULL) 
{
  colnames(cos_sim_df)[1] <- "Sample"
  # if (class(cos_sim_matrix) != "matrix") {
  #   stop("cos_sim_matrix must be a matrix")
  # }
  if (length(colnames(cos_sim_df)) == 0) {
    stop("cos_sim_df is missing colnames")
  }
  if (length(rownames(cos_sim_df)) == 0) {
    stop("cos_sim_df is missing rownames")
  }
  
  # limited to max row/col 
  
  if(dim(cos_sim_df)[1] > nmax | dim(cos_sim_df)[2] > (nmax+1)){
    nrow <- if_else(dim(cos_sim_df)[1]>nmax, nmax, dim(cos_sim_df)[1])
    ncol <- if_else(dim(cos_sim_df)[2]>nmax+1L, nmax+1L, dim(cos_sim_df)[2])
    cos_sim_df <- cos_sim_df[1:nrow,1:ncol]
  }
  
  
  # covert to matrix
  cos_sim_matrix <- as.matrix(cos_sim_df[,-1])
  rownames(cos_sim_matrix) <- cos_sim_df[[1]]
  
  if (missing(col_order)) {
    #col_order = colnames(cos_sim_df)[-1]
    hc.sample = hclust(dist(t(cos_sim_matrix)), method = method)
    col_order = rownames(t(cos_sim_matrix))[hc.sample$order]
  }
  
  if (class(col_order) != "character") {
    stop("col_order must be a character vector")
  }
  if (length(col_order) != ncol(cos_sim_df)-1) {
    stop("col_order must have the same length as the number of signatures in the explained df")
  }
  
  if (cluster_rows == TRUE) {
    hc.sample = hclust(dist(cos_sim_matrix), method = method)
    sample_order = rownames(cos_sim_matrix)[hc.sample$order]
  }
  else {
    sample_order = rownames(cos_sim_matrix)
  }
  Cosine.sim = NULL
  Signature = NULL
  Sample = NULL
  x = NULL
  y = NULL
  xend = NULL
  yend = NULL
  #cos_sim_matrix.m = melt(cos_sim_matrix)
  #colnames(cos_sim_matrix.m) = c("Sample", "Signature", "Cosine.sim")
  cos_sim_matrix.m <- cos_sim_df %>% pivot_longer(-1, names_to="Signature",values_to="Cosine.sim") %>% select(Sample,Signature,Cosine.sim)
  cos_sim_matrix.m$Signature = factor(cos_sim_matrix.m$Signature,  levels = col_order)
  cos_sim_matrix.m$Sample = factor(cos_sim_matrix.m$Sample,  levels = sample_order)
  
  cos_sim_matrix.m <- cos_sim_matrix.m %>% mutate(coslab=round(Cosine.sim, 2))
  cos_sim_matrix.m$Cosine.sim <- round(cos_sim_matrix.m$Cosine.sim,digits = 2)
  mincosine <- min(cos_sim_matrix.m$Cosine.sim)
  maxcosine <- max(cos_sim_matrix.m$Cosine.sim)
  
  if(mincosine< -0.1){
    plot_values <- TRUE
    cos_sim_matrix.m$Cosine.sim <- abs(cos_sim_matrix.m$Cosine.sim)
    mincosine <- min(cos_sim_matrix.m$Cosine.sim)
    cos_sim_matrix.m$coslab <- if_else(cos_sim_matrix.m$coslab>0,"","-")
  }
  mincosine <- if_else(mincosine<0.1,0,mincosine)
  maxcosine <- if_else(maxcosine>0.9,1,maxcosine)
  
  
  ## define the length of x and y
  name_max_x <- max(str_length(unique(cos_sim_matrix.m$Signature)))
  name_max_y <- max(str_length(unique(cos_sim_matrix.m$Sample)))
  name_len_x <- length(unique(cos_sim_matrix.m$Signature))
  name_len_y <- length(unique(cos_sim_matrix.m$Sample))
  
  angle_max_x <- if_else(name_max_x < 15, 90, 30)
  vjust_max_x <- if_else(name_max_x < 15, 0.5, 1)
  
  leng0 <- 2.5
  leng_ratio <-  0.2
  lenx <- name_max_x*0.4
  leny <- name_max_y*0.4
  xleng <- leng_ratio*name_len_x+leng0+3+leny
  yleng <- leng_ratio*name_len_y+leng0+1+lenx
  
  xleng <- if_else(xleng>30,30,xleng)
  yleng <- if_else(yleng>25,25,yleng)
  
  heatmap = ggplot(cos_sim_matrix.m, aes(x = Signature, y = Sample,  fill = Cosine.sim, order = Sample)) + 
    geom_tile(color = "white") + 
    scale_fill_viridis_c(name = "Cosine similarity\n", limits = c(mincosine, maxcosine),breaks=scales::pretty_breaks())+
    # scale_fill_distiller(palette = "YlGnBu", direction = 1,name = "Cosine similarity", limits = c(0, 1)) + 
    theme_ipsum_rc(grid = FALSE,ticks = T)+
    #theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +
    labs(x = NULL, y = NULL)+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    theme(legend.position = "top",legend.key.width = unit(2,"cm"))
  #legend.key.width = unit(1.5,"cm"),legend.key.height = unit(0.3,"cm")
  
  
  if(name_len_x<=80){
    heatmap <- heatmap + theme(axis.text.x = element_text(angle = angle_max_x,hjust = 1,vjust = vjust_max_x,size = 10))
  } else{
    heatmap <- heatmap + theme(axis.text.x = element_blank())
  }
  
  if(name_len_y>80){
    heatmap <- heatmap + theme(axis.text.y = element_blank())
  }
  
  
  if (plot_values) {
    heatmap = heatmap + geom_text(aes(label = coslab), size = 3,col="red")
  }
  if (cluster_rows == TRUE) {
    dhc = as.dendrogram(hc.sample)
    ddata = ggdendro::dendro_data(dhc, type = "rectangle")
    dendrogram = ggplot(ggdendro::segment(ddata)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + 
      scale_y_reverse(expand = c(0.1, 0)) + scale_x_continuous(expand = expansion(add=0.5)) + ggdendro::theme_dendro()
    plot_final = plot_grid(dendrogram+theme(plot.margin = margin(r = -0.2,unit = "cm")), heatmap+theme(plot.margin = margin(l = -0.2,r=0.5,unit = "cm")), align = "h",axis = "tb",rel_widths = c(0.15, 1))
  } else {
    plot_final = heatmap + ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  }
  
  
  if(is.null(output_plot)){
    return(plot_final)
  }else{
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot_final,width = plot_width,height = plot_height)
  }
  
}




# Plot two profile difference for SBS96, ID83 and DBS78 ---------------------------------------------
plot_compare_profiles_diff <- function (profile1, profile2, profile_names = NULL, profile_ymax = NULL, diff_ylim = NULL, colors = NULL, condensed = FALSE,output_plot = NULL,plot_width=NULL, plot_height=NULL) 
{
  ## profile 1 and profile 2 will be the dataframe with two columns: MutationType and value
  
  COLORS6 = c("#03BCEE", "#010101", "#E32926", "#CAC9C9", "#A1CE63", "#EBC6C4")
  COLORS10 = c("#03BCEE", "#0366CB", "#A1CE63", "#016601", "#FE9898", "#E32926", "#FEB166", "#FE8001", "#CB98FE", "#4C0198")
  COLORS16=c("#FCBD6F", "#FE8002", "#AFDC8A", "#36A02E", "#FCC9B4", "#FB896A", "#F04432", "#BB191A", "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA", "#E1E1EE", "#B5B5D7", "#8582BC", "#62409A")
  
  typelength = dim(profile1)[1]
  if (is.null(colors)) {
    if(typelength == 96){colors = COLORS6; }
    if(typelength == 78){colors = COLORS10;}
    if(typelength == 83){colors = COLORS16;}
  }
  
  ## make sure the order profile 1 and profile 2 are the same order
  
  indel_short <-  FALSE
  if(typelength == 83) {indel_short = TRUE}
  profile1 <- profile_format_df(profile1,factortype = TRUE,indel_short = indel_short) 
  profile2 <- profile_format_df(profile2,factortype = TRUE,indel_short = indel_short)
  
  
  names(colors) <- levels(profile1$Type)
  
  #print(colors)
  
  profile1[,4] <- profile1[,4]/sum(profile1[,4])  
  profile2[,4] <- profile2[,4]/sum(profile2[,4]) 
  diff = profile1[[4]] - profile2[[4]]
  RSS = sum(diff^2)
  RSS = format(RSS, scientific = TRUE, digits = 3)
  cosine_sim = cos_sim(profile1[[4]], profile2[[4]])
  cosine_sim = round(cosine_sim, 3)
  df <-  profile1 %>% left_join(profile2) %>% mutate(Difference=diff)
  if(is.null(profile_names)){
    profile_names <- colnames(df)[4:5]
  }
  colnames(df)[4:5] <- profile_names
  df <- df %>% pivot_longer(cols = -c(1,2,3)) %>% mutate(name=factor(name,levels = c(profile_names,"Difference"))) %>% mutate(Type=factor(Type,levels = names(colors)))
  
  if(is.null(profile_ymax)){
    profile_ymax <- max(c(profile1[[4]],profile2[[4]]))*1.1
  }
  
  if(is.null(diff_ylim)){
    diff_ylim <- range(diff)*1.1
  }
  dftmp = tibble(Type = rep(levels(profile1$Type)[1], 4), SubType = rep((profile1$SubType)[1], 4), name = c(profile_names, "Difference", "Difference"), value = c(profile_ymax, profile_ymax, diff_ylim[1], diff_ylim[2])) %>% mutate(name=factor(name,levels = c(profile_names,"Difference"))) %>% mutate(Type=factor(Type,levels = names(colors)))
  
  
  if (condensed) {
    plot = ggplot(data = df, aes(x = SubType, y = value,  fill = Type, width = 1)) + 
      geom_bar(stat = "identity", position = "identity", colour = "black", size = 0.2) + 
      geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
      scale_fill_manual(values = colors) + 
      facet_grid(name ~ Type, scales = "free",space = "free_x") +
      ylab("Relative contribution") +
      guides(fill = FALSE) + 
      labs(x="")+
      #scale_y_continuous(expand = c(0,0))+
      theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
      ggtitle(paste("RSS = ", RSS, "; Cosine similarity = ", cosine_sim, sep = ""))+
      theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5), strip.text.x = element_text(size = 12,hjust = 0.5,colour = "white", margin = margin()), strip.text.y = element_text(size = 12,hjust = 0.5, margin = margin()),strip.background = element_rect(fill = "#f0f0f0"), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(hjust = 0.5),axis.line.y = element_line(colour = 'black',size = 0.25))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,colour = 'black',size = 0.5)+geom_vline(xintercept = Inf,colour = 'black',size = 0.5)#axis.line.x = element_line(colour = 'black',size = 0.25),
    #panel_border(color = gray(0.5),size = 0.3)
  }
  else {
    plot = ggplot(data = df, aes(x = SubType, y = value,  fill = Type, width = 0.7)) + 
      geom_bar(stat = "identity", position = "identity", colour = "black", size = 0) + 
      geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
      scale_fill_manual(values = colors) + 
      facet_grid(name ~ Type, scales = "free",space = "free_x") +
      ylab("Relative contribution") +
      guides(fill = FALSE) + 
      labs(x="")+
      #scale_y_continuous(expand = c(0,0))+
      theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
      ggtitle(paste("RSS = ", RSS, "; Cosine similarity = ", cosine_sim, sep = ""))+
      theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5), strip.text.x = element_text(size = 12,hjust = 0.5,colour = "white", margin = margin()), strip.text.y = element_text(size = 12,hjust = 0.5, margin = margin()),strip.background = element_rect(fill = "#f0f0f0",), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(hjust = 0.5),axis.line.y = element_line(colour = 'black',size = 0.25))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,colour = 'black',size = 0.5) +geom_vline(xintercept = Inf,colour = 'black',size = 0.5) #,axis.line.x = element_line(colour = 'black',size = 0.25)
    #     panel_border(color = gray(0.5),size = 0.3)
  }
  
  
  ## add background color for strip
  require(grid)
  g <- ggplot_gtable(ggplot_build(plot))
  strip_top <- which(grepl('strip-t', g$layout$name))
  fills <- colors
  k <- 1
  for (i in strip_top) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot <- as_ggplot(g)
  
  
  
  if(is.null(output_plot)){
    return(plot)
  }else{
    xleng <- 14
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot,width = plot_width,height = plot_height)
  }
  
}


plot_compare_profiles_diff_strand <- function (profile1, profile2, profile_names = NULL, profile_ymax = NULL, diff_ylim = NULL, colors = NULL, condensed = FALSE) 
{
  
  ## check format: MutationType Strand Type  SubType value
  ## profile 1 and profile 2 will be the dataframe with two columns: MutationType and value
  
  ## make sure the order profile 1 and profile 2 are the same order
  profile1 <- profile1 %>% arrange(MutationType)
  profile2 <- profile2 %>% arrange(MutationType)
  
  profile1[,5] <- profile1[,5]/sum(profile1[,5])  
  profile2[,5] <- profile2[,5]/sum(profile2[,5]) 
  diff = profile1[[5]] - profile2[[5]]
  RSS = sum(diff^2)
  RSS = format(RSS, scientific = TRUE, digits = 3)
  cosine_sim = cos_sim(profile1[[5]], profile2[[5]])
  cosine_sim = round(cosine_sim, 3)
  df <-  profile1 %>% left_join(profile2) %>% mutate(Difference=diff)
  if(is.null(profile_names)){
    profile_names <- colnames(df)[5:6]
  }
  colnames(df)[5:6] <- profile_names
  df <- df %>% pivot_longer(cols = -c(1,2,3,4)) %>% mutate(name=factor(name,levels = c(profile_names,"Difference")))
  
  if(is.null(profile_ymax)){
    profile_ymax <- max(c(profile1[[5]],profile2[[5]]))*1.1
  }
  
  if(is.null(diff_ylim)){
    diff_ylim <- range(diff)*1.1
  }
  dftmp = tibble(Type = rep(levels(as.factor(profile1$Type))[1], 4), SubType = rep(levels(as.factor(profile1$SubType))[1], 4), name = c(profile_names, "Difference", "Difference"), value = c(profile_ymax, profile_ymax, diff_ylim[1], diff_ylim[2])) %>% mutate(Strand="T")%>% mutate(name=factor(name,levels = c(profile_names,"Difference")))
  
  if (condensed) {
    plot = ggplot(data = df, aes(x = SubType, y = value,fill=Strand,group=Strand, width = 1)) + 
      geom_col(position = "dodge2")+
      #geom_bar(stat = "identity", position = "identity", colour = "black", size = 0.2) + 
      geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
      scale_fill_manual(values = rev(pal_nejm()(2)))+
      facet_grid(name ~ Type, scales = "free") +
      ylab("Relative contribution") +
      guides(fill = FALSE) + 
      labs(x="")+
      #scale_y_continuous(expand = c(0,0))+
      theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
      ggtitle(paste("RSS = ", RSS, "; Cosine similarity = ", cosine_sim, sep = ""))+
      theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5), strip.text.x = element_text(size = 12,hjust = 0.5), strip.text.y = element_text(size = 14,hjust = 0.5),strip.background = element_rect(fill = "#f0f0f0"), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(hjust = 0.5),axis.line.x = element_line(colour = 'black',size = 0.25),axis.line.y = element_line(colour = 'black',size = 0.25))+geom_vline(xintercept = Inf,colour = 'black',size = 0.5)
    #      panel_border(color = gray(0.5),size = 0.3)
  }
  else {
    plot = ggplot(data = df, aes(x = SubType, y = value, fill=Strand,group=Strand, width = 0.7)) + 
      #geom_bar(stat = "identity", position = "identity", colour = "black", size = 0) + 
      geom_col(position = "dodge2")+
      geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
      scale_fill_manual(values = rev(pal_nejm()(2)))+
      facet_grid(name ~ Type, scales = "free") +
      ylab("Relative contribution") +
      guides(fill = FALSE) + 
      labs(x="")+
      #scale_y_continuous(expand = c(0,0))+
      theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
      ggtitle(paste("RSS = ", RSS, "; Cosine similarity = ", cosine_sim, sep = ""))+
      theme(axis.title.y = element_text(size = 12, vjust = 1), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5), strip.text.x = element_text(size = 14,hjust = 0.5), strip.text.y = element_text(size = 14,hjust = 0.5),strip.background = element_rect(fill = "#f0f0f0"), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(hjust = 0.5),axis.line.x = element_line(colour = 'black',size = 0.25),axis.line.y = element_line(colour = 'black',size = 0.25))+geom_vline(xintercept = Inf,colour = 'black',size = 0.5)
    # panel_border(color = gray(0.5),size = 0.3)
  }
  return(plot)
}




# signature_sum_operation  --------------------------------------------------------------
## use the * and ;

signature_sum_operation <- function(sigdatabase,sigsetname,formulax,outsigname="Contribution") {
  #sigsetname <- 'COSMIC v3 Signatures (SBS)'
  #formulax <- "0.4*SBS1;0.4*SBS2;0.2*SBS13"
  formulax <- str_trim(formulax)
  inforx <- (str_split(unlist(str_split(formulax,"\\;")),"\\*"))
  #sigalltmp <- tibble(MutationType=character(),Contribution=numeric())
  for(i in 1:length(inforx)){
    tmp <- inforx[[i]]
    ratio <- as.numeric(tmp[[1]])
    signame <- tmp[[2]]
    sigtmp <- sigdatabase %>% filter(Signature_set_name==sigsetname,Signature_name==signame)
    sigtmp <- sigtmp %>% select(MutationType,tmpC=Contribution) %>% mutate(tmpC=ratio*tmpC)
    if(i==1){
      colnames(sigtmp)[2] <- 'Contribution'
      sigalltmp <- sigtmp
    }else{
      sigalltmp <- left_join(sigalltmp,sigtmp) %>% mutate(Contribution=Contribution+tmpC) %>% select(MutationType,Contribution)
    }
  }
  consum <- sum(sigalltmp$Contribution)
  sigalltmp$Contribution <- sigalltmp$Contribution/consum
  colnames(sigalltmp)[2] <- outsigname
  return(sigalltmp)
}




# Define Signature Set Colors  -------------------------------------------------------
sigsetcolor <- c(
  "Cancer Reference Signatures (RS)" = "#35978f",
  "Cancer Reference Signatures (SBS)"  = "#01665e",
  "COSMIC v2 Signatures (SBS)" = "#c994c7",
  "COSMIC v3 Signatures (DBS)" = "#df65b0",
  "COSMIC v3 Signatures (ID)" = "#e7298a",
  "COSMIC v3 Signatures (SBS)" = "#d73027",
  "Environmental Mutagen Signatures (SBS)" = "#ff7f00",
  "Organ-specific Cancer Signatures (RS)" = "#74a9cf",
  "Organ-specific Cancer Signatures (SBS)" = "#0570b0",
  "Other published signatures (ID)" = "#969696",
  "Other published signatures (SBS)" = "#525252",
  "SignatureAnalyzer PCAWG WGS 1536 Signatures (SBS)" = "#d9ef8b",
  "SignatureAnalyzer PCAWG WGS Signatures (DBS)" = "#a6d96a",
  "SignatureAnalyzer PCAWG WGS Signatures (ID)" = "#66bd63",
  "SignatureAnalyzer PCAWG WGS Signatures (SBS)" = "#1a9850",
  "SigProfiler PCAWG Strand Signatures (SBS)" = "#8073ac",
  "SigProfiler PCAWG WXS Signatures (SBS)" = "#542788"
)



# Signature pie chart -----------------------------------------------------

signature_piechart <- function(data,colset, output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  # convert nsig_data to another format for the piechart
  nsig_data_pie <- data %>%
    group_by(Profile) %>% 
    arrange(N) %>%
    mutate(
      end_angle = 2*pi*cumsum(N)/sum(N),   # ending angle for each pie slice
      start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
      mid_angle = 0.5*(start_angle + end_angle),   # middle of each pie slice, for the text label
      # horizontal and vertical justifications depend on whether we're to the left/right
      # or top/bottom of the pie
      hjust = ifelse(mid_angle > pi, 1, 0),
      vjust = ifelse(mid_angle < pi/2 | mid_angle > 3*pi/2, 0, 1)
    ) %>% 
    ungroup()
  
  # radius of the pie and radius for outside and inside labels
  rpie <- 1
  rlabel_out <- 1.05 * rpie
  rlabel_in <- 0.6 * rpie
  
  
  
  p <- ggplot(nsig_data_pie) +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie, start = start_angle, end = end_angle, fill = Signature_set_name)) +
    geom_text(aes(x = rlabel_in * sin(mid_angle), y = rlabel_in * cos(mid_angle), label = N2 ), size = 14/.pt)+
    facet_wrap(~Profile,nrow = 2)+
    coord_fixed()+
    labs(x="",y="")+
    scale_fill_manual("Signature Set Name",values = colset)+
    theme_ipsum_rc(axis = FALSE, grid = FALSE)+
    theme(axis.text.y = element_blank(),axis.text.x = element_blank(),strip.text = element_text(size = 14,hjust = 0.5,face = "bold"))
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 16
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
} 



## TMB plot ###
TMBplot <- function(data,output_plot = NULL,plot_width=NULL, plot_height=NULL,addnote= NULL){
  # Cancer_Type,     Sample,   Burden
  # group Cancer Type
  data <- data %>% mutate(Burden=if_else(is.infinite(Burden),NA_real_,Burden))
  mdata <- suppressMessages(data %>% group_by(Cancer_Type) %>% summarise(median=median(Burden,na.rm = TRUE),total=n(),nsample=sum(!is.na(Burden),na.rm=TRUE)) %>% arrange(median) %>% mutate(Seq=seq_along(Cancer_Type)*10-10))
  ## remove na cancer type
  natype <- mdata %>% filter(is.na(median)) %>% pull(Cancer_Type)
  mdata <- mdata %>% filter(!is.na(median))
  data <- data %>% filter(!(Cancer_Type %in% natype))
  
  data <- suppressMessages(data %>% mutate(Cancer_Type=factor(Cancer_Type,levels = mdata$Cancer_Type)) %>% arrange(Cancer_Type,Burden) %>% left_join(mdata))
  
  data <- data %>% filter(!is.na(Burden))%>% group_by(Cancer_Type) %>% mutate(Subseq=seq_along(Sample)) %>% ungroup() %>% mutate(score=Seq+1+8*Subseq/nsample)
  mdata <- mdata %>% mutate(Seq2=Seq+10) %>% mutate(score=(Seq+Seq2)/2) %>% 
    mutate(Burden=median,Type=if_else(seq_along(Cancer_Type) %% 2 ==0,"A","B"),Label=paste0(nsample,'<br /> - <br />',total)) %>% 
    mutate(Label_col = glue::glue("<span style='color:#377eb8'>{nsample}</span><br />\u2015<br /><span style='color:#4daf4a'>{total}</span>"))
  # "<span style='color:#377eb8'>{nsample}</span><br />\u2501<br /><span style='color:#4daf4a'>{total}</span>"
  f <- function(y) seq(floor(min(y,na.rm = TRUE)), ceiling(max(y,na.rm = TRUE)))
  
  ymint <- floor(min(data$Burden,na.rm = TRUE))
  ymint <- if_else(is.infinite(ymint),-4,ymint)
  ymaxt <- ceiling(max(data$Burden,na.rm = TRUE))
  p <-  data %>% filter(!is.na(Burden)) %>% 
    ggplot(aes(score,(Burden),group=Cancer_Type))+
    geom_rect(data=mdata,aes(xmin=Seq,xmax=Seq2,ymin=-Inf,ymax=Inf,fill=Type),alpha=0.5)+
    geom_point(pch=21,size=1,stroke=0.5)+
    geom_segment(data=mdata,aes(x=Seq+1,y=(median),xend=Seq2-1,yend=(median)),col="#e41a1c")+
    scale_fill_manual(values=c('#F2F2F2','#D4D4D4'))+
    scale_x_continuous(limits = c(min(mdata$Seq),max(mdata$Seq2)),expand = c(0,0),breaks = mdata$score,labels = mdata$Label_col,sec.axis = dup_axis(labels = mdata$Cancer_Type))+
    scale_y_continuous(limits = c(ymint,ymaxt),expand = c(0,0),breaks = f)+
    labs(x="",y="Number of Mutations per Megabase\n(log10)")+
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,axis_text_size = 12,grid = "Y",ticks = FALSE)+
    theme(legend.position = "None",axis.text.x.top = element_text(size = 11,angle = -45,hjust = 1,vjust = 0),axis.text.x.bottom = element_markdown(size = 12),panel.border = element_rect(colour = "black",fill = NA))+
    coord_cartesian(clip = 'off')
  
  if(!is.null(addnote)){
    p <- p+annotate("text",x=-Inf, y = Inf, label = addnote, vjust=2, hjust=-0.2,size=7)
  }
  
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 2.5+length(unique(mdata$Cancer_Type))*0.5
    xleng <- if_else(xleng>18,18,xleng)
    yleng <- 6
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
}





decompsite_distribution <- function(decompsite,output_plot = NULL,plot_width=NULL, plot_height=NULL){
  #Cancer_Type   Sample Total_Mutations Cosine_similarity .... 
  fealist <- c('Cancer_Type','Sample','Cosine_similarity','100-L1_Norm_%','100-L2_Norm_%','KL_Divergence','Correlation')
  decompsite2 <- decompsite %>% select(one_of(fealist)) %>% pivot_longer(cols = -c(Cancer_Type,Sample))
  mtmp <- decompsite %>% group_by(Cancer_Type) %>% summarise(m=median(Cosine_similarity,na.rm = TRUE)) 
  mtmp2 <- decompsite2 %>% group_by(Cancer_Type,name) %>% summarise(m=median(value,na.rm = TRUE)) %>% ungroup() %>% mutate(m=if_else(name=="KL_Divergence",1-m,m))%>% group_by(name) %>% arrange(desc(m)) %>% mutate(Seq=seq_along(Cancer_Type)) %>% ungroup()
  ntmp <- decompsite %>% count(Cancer_Type) %>% mutate(Caner_type2=paste0(Cancer_Type," (",n,")")) %>% left_join(mtmp) %>% arrange(m)
  
  library(ggridges)
  library(scales)
  p <- decompsite2 %>% 
    left_join(mtmp2) %>% 
    left_join(ntmp) %>% 
    mutate(Cancer_Type=factor(Cancer_Type,levels = ntmp$Cancer_Type,labels =ntmp$Caner_type2 )) %>% 
    mutate(name=factor(name,levels = c('Cosine_similarity','100-L1_Norm_%','100-L2_Norm_%','KL_Divergence','Correlation'))) %>% 
    mutate(value=if_else(value<0,0,value)) %>% 
    ggplot(aes(value,y=Cancer_Type,fill=Seq,height = ..ndensity..))+
    geom_density_ridges(color="black")+
    facet_wrap(~name,scales = 'free_x',nrow = 1)+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = "X",ticks = TRUE,base_size = 12,grid = FALSE,strip_text_size = 12)+
    scale_x_continuous(name="",expand = c(0,0),breaks = pretty_breaks())+
    scale_y_discrete(name = "",expand = c(0,0))+
    scale_fill_viridis_c(direction = -1)+
    theme(strip.text.x = element_text(hjust = 0.5),legend.position = "none",panel.grid.major.x = element_line(colour = gray(0.85),size = 0.5,linetype = 2),panel.spacing = unit(0.8, "lines"),strip.background = element_rect(colour = '#cccccc',fill = c('#a1d99b')))+
    #panel_border(color = "black",size = 0.8)+
    guides(fill=guide_legend(title = "",nrow = 4))
  
  if(is.null(output_plot)){
    return(p)
  }else{
    yleng <- 2.5+length(unique(ntmp$Cancer_Type))*0.5
    yleng <- if_else(yleng>12,12,yleng)
    xleng <- 16
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
  
}




### Landscape of mutational signature cluster function ####
Exposure_Clustering <- function(sigdata,sigcolor=NULL,studydata=NULL,studydata_cat = TRUE, studycolor=NULL,study_cutoff = -Inf,studyname = NULL,puritydata=NULL,puritydata_cat=FALSE,puritycol=NULL,purity_cutoff=-Inf,purityname = NULL,cosinedata,cosine_cutoff=0,highlight=NULL,legendnrow=NULL,sampletext_size=6,output_plot=NULL,plot_height=NULL,plot_width=NULL,clustern,hc_func='hclust',hc_metric = 'euclidean',hc_method = 'ward.D2',stand=TRUE){
  require(tidyverse)
  require(scales)
  require(janitor)
  require(ggsci)
  require(ggpubr)
  require(factoextra)
  require(cowplot)
  
  # define color for categoly variable
  colset <- c(pal_npg()(10),pal_jama()(7),pal_igv()(51))
  #colset <- pal_igv()(51)
  
  colnames(sigdata)[1] <- 'Samples'
  
  ## define color for signature
  if(is.null(sigcolor)){
    uvalues <- colnames(sigdata)[-1] #sort
    if(unique(uvalues %in% names(Subscolor)) == TRUE && length(unique(uvalues %in% names(Subscolor))) == 1){
      sigcolorindex <- as.character(Subscolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else if(unique(uvalues %in% names(SBScolor)) == TRUE && length(unique(uvalues %in% names(SBScolor))) == 1){
      sigcolorindex <- as.character(SBScolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else {
      sigcolorindex <- as.character(SBScolor[1:length(uvalues)])
      names(sigcolorindex) <- uvalues
    }
  }else{
    sigcolorindex <- sigcolor
  }
  
  # cluster data
  tmp=sigdata %>% adorn_percentages('row') 
  #%>% mutate_if(is.numeric, funs(replace_na(., 0)))
  mdata=as.matrix(tmp[,-1])
  rownames(mdata) <- tmp$Samples
  
  #fviz_nbclust(mdata, kmeans, method = "gap_stat")
  
  clustern <- if_else(dim(mdata)[1]<clustern+5,2L,as.integer(clustern))
  
  kcolors <- pal_d3("category20")(clustern)
  res <- hcut(mdata,k = clustern,hc_func = hc_func,hc_metric = hc_metric,hc_method = hc_method,stand=stand)
  
  # p_cluster
  p_cluster <- fviz_dend(res, rect = TRUE, cex = 0.5,k_colors = kcolors,lwd = 0.5,show_labels = F)+scale_x_discrete(expand = c(0,0))+ theme(plot.margin=margin(b=-0.7,unit="cm"),title = element_blank())
  #plot.margin=margin(b=-1,unit="cm")
  
  ## p_mutation
  p_mutation <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,(Weight),fill=factor(Signature,levels = names(sigcolorindex))))+geom_bar(stat="identity",col="gray95",width=1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = "xy")+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),panel.grid.major.x=element_blank(),legend.position = "none",legend.box.spacing = unit(0,"cm"),plot.margin=margin(b=-1,t = 1,unit="cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks(),labels = comma)+scale_fill_manual(values = sigcolorindex,drop=FALSE)+guides(fill=guide_legend(nrow=2,byrow=TRUE))+xlab("")+ylab("Number of mutations \n")
  #+theme(plot.margin=margin(b=4,unit="pt"))
  
  # p_proportion plot
  if(is.null(legendnrow)){
    nsize <- dim(sigdata)[2]-1
    if(nsize > 15){
      legendnrow <- 2
    }else {
      legendnrow <- 1
    }
  }
  
  p_proportion <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,Weight,fill=factor(Signature,levels = names(sigcolorindex))))+geom_bar(stat="identity",position="fill",col="gray95",width = 1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,grid = FALSE)+theme(panel.background = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_text(size = sampletext_size, angle = 90, vjust = 0.5, hjust = 1),panel.grid.major=element_line(),legend.position = "bottom",legend.box.background = element_blank(),legend.box.spacing = unit(-0.5,"cm"),legend.key = element_rect(size = 0),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.25, "cm"),legend.key.width =unit(1.5, "cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks())+xlab("")+scale_fill_manual("Mutational sigantures\n",values = sigcolorindex,drop=FALSE)+guides(fill=guide_legend(nrow=legendnrow,byrow=TRUE,label.position = "bottom"))+ylab("Signature contribution\n")
  #+theme(plot.margin=margin(t=4,unit="pt"))
  #legend.title = element_blank(),
  
  p_proportion_legend <- as_ggplot(get_legend(p_proportion))
  #p_proportion_legend <- p_proportion_legend+theme(plot.margin = margin(b = 0))
  p_proportion <- p_proportion + theme(legend.position = "none")
  
  if(!is.null(highlight)){
    samhigh <- sigdata %>% mutate(Samples_high=if_else(Samples %in% highlight,paste0("*",Samples),Samples))
    p_proportion <-  p_proportion+scale_x_discrete(breaks=samhigh$Samples,labels=samhigh$Samples_high)
  }
  
  
  #p3 <- flush_ticks(p3,flush = "Y",plot = FALSE)
  #p4 <- flush_ticks(p4,flush = "Y",plot = FALSE)
  
  # cosine similarity bar
  p_cosine <- cosinedata %>%  mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% mutate(Similarity=if_else(Similarity<cosine_cutoff,NA_real_,Similarity)) %>% ggplot(aes(Samples,1,fill=Similarity))+geom_tile(col="black")+scale_fill_viridis_c("Cosine similarity\n",na.value = "#cccccc",option = 'C',limits = c(0.6, 1), oob = scales::squish)+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(2,"cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))
  #+theme(title = element_blank())
  ## p_cosine_legend = cosine similarity bar (legend)
  p_cosine_legend <- as_ggplot(get_legend(p_cosine))
  #p_cosine_legend <- p_cosine_legend+theme(plot.margin = margin(b = 0))
  p_cosine <- cosinedata %>%  mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% mutate(Similarity=if_else(Similarity<cosine_cutoff,NA_real_,Similarity)) %>% ggplot(aes(Samples,1,fill=Similarity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'C',limits = c(0.6, 1), oob = scales::squish)+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1.2,t = 0.5,unit="cm"))+ylim(c(0,2))
  #,title = element_blank()
  
  legend_size <- theme(legend.position = "bottom",legend.direction = "horizontal", legend.box.background = element_blank(),legend.key = element_rect(size = 0),legend.key.size = unit(0.25, "cm"),legend.key.width =unit(0.6, "cm"))
  #legend.box.spacing = unit(-0.5,"cm"),
  
  # study color bar
  if(!is.null(studydata)){
    colnames(studydata) <- c('Samples','Study')
    studydata <- studydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) 
    studyname <- if_else(is.null(studyname),"",studyname)
    if(!studydata_cat) {
      studyname <- paste0(studyname,"\n")
      studydata <- studydata %>% mutate(Study=if_else(Study<study_cutoff,NA_real_,Study))
      p_study <- studydata  %>% ggplot(aes(Samples,1,fill=Study))+geom_tile(col="black")+scale_fill_viridis_c(studyname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(2, "cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))
      p_study_legend <- as_ggplot(get_legend(p_study+legend_size))
      #p_study_legend <- p_study_legend+theme(plot.margin = margin(b = 0))
      p_study <- studydata %>% ggplot(aes(Samples,1,fill=Study))+geom_tile(col="black")+scale_fill_viridis_c(studyname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t = 0.5,unit="cm"))+ylim(c(0,2))
      
    }else { 
      if(is.null(studycolor)){
        studycolor <-  colset[1:length(unique(studydata$Study))]
        names(studycolor) <- unique(studydata$Study)
      }
      p_study <- studydata %>% ggplot(aes(Samples,1,fill=factor(Study,levels = names(studycolor))))+geom_tile(col="black")+scale_fill_manual(studyname,values =studycolor)+theme_minimal()+theme(legend.position = "bottom",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"))+ylim(c(0,2))+guides(fill = guide_legend(nrow = 1))
      p_study_legend <- as_ggplot(get_legend(p_study+legend_size))
      #p_study_legend <- p_study_legend+theme(plot.margin = margin(b = 0))
      p_study <- studydata %>% ggplot(aes(Samples,1,fill=factor(Study,levels = names(studycolor))))+geom_tile(col="black")+scale_fill_manual(studyname,values =studycolor)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"))+ylim(c(0,2))
    }
  }else {
    p_study <- NULL
    p_study_legend <- NULL
  }
  
  # purity color bar
  # if(is.null(purityname)){
  #   purityname <-  "              "
  # }
  
  if(!is.null(puritydata)){
    colnames(puritydata) <- c('Samples','Purity')
    puritydata <- puritydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) 
    purityname <- if_else(is.null(purityname),"",purityname)
    if(!puritydata_cat) {
      purityname <-  paste0(purityname,"\n")
      puritydata <- puritydata %>% mutate(Purity=if_else(Purity<purity_cutoff,NA_real_,Purity))
      p_purity <- puritydata  %>% ggplot(aes(Samples,1,fill=Purity))+geom_tile(col="black")+scale_fill_viridis_c(purityname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(2, "cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))
      p_purity_legend <- as_ggplot(get_legend(p_purity+legend_size))
      #p_purity_legend <- p_purity_legend+theme(plot.margin = margin(b = 0))
      p_purity <- puritydata %>% ggplot(aes(Samples,1,fill=Purity))+geom_tile(col="black")+scale_fill_viridis_c(purityname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t = 0.5,unit="cm"))+ylim(c(0,2))
    }else { 
      if(is.null(puritycol)){
        puritycol <- colset[1:length(unique(puritydata$Purity))]
        names(puritycol) <- unique(puritydata$Purity)
      }
      p_purity <- puritydata %>% ggplot(aes(Samples,1,fill=factor(Purity,levels = names(puritycol))))+geom_tile(col="black")+scale_fill_manual(values =puritycol)+theme_minimal()+theme(legend.position = "bottom",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"))+ylim(c(0,2))+guides(fill = guide_legend(nrow = 1))
      p_purity_legend <- as_ggplot(get_legend(p_purity+legend_size))
      #p_purity_legend <- p_purity_legend+theme(plot.margin = margin(b = 0))
      p_purity <- puritydata %>% ggplot(aes(Samples,1,fill=factor(Purity,levels = names(puritycol))))+geom_tile(col="black")+scale_fill_manual(purityname,values =puritycol)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"))+ylim(c(0,2))
    }
  }else {
    p_purity <- NULL
    p_purity_legend <- NULL
  }
  
  h_study <- if_else(is.null(p_study),0,0.14)
  h_purity <- if_else(is.null(p_study),0,0.11)
  
  p_legends <- plot_grid(p_study_legend+theme(plot.margin = margin(b = -2,t=0,unit = 'cm')), p_purity_legend+theme(plot.margin = margin(b = -2,t=0,unit = 'cm')), p_proportion_legend+theme(plot.margin = margin(t = -1,b=-0.5,unit = 'cm')), p_cosine_legend+theme(plot.margin = margin(t = -1,b=-0.5,unit = 'cm')),nrow = 2,ncol = 2)
  
  if(!is.null(p_purity) & !is.null(p_study) ){
    pall <- plot_grid(p_cluster,p_study,p_purity,p_mutation,p_cosine,p_proportion,p_legends+theme(plot.margin = margin(t=-3,unit = 'cm')),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,h_study,h_purity,4,0.05,7,1))
  }else if(is.null(p_purity) & is.null(p_study)){
    p_legends <- plot_grid(p_proportion_legend+theme(plot.margin = margin(t = 2,b=-0.5,unit = 'cm')), p_cosine_legend+theme(plot.margin = margin(t = -1,b=-0.5,unit = 'cm')),nrow=2,ncol = 1)
    pall <- plot_grid(p_cluster,p_mutation,p_cosine,p_proportion,p_legends+theme(plot.margin = margin(t=-3,unit = 'cm')),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,4,0.05,7,1))
  } else if(is.null(p_purity)){
    pall <- plot_grid(p_cluster,p_study,p_mutation,p_cosine,p_proportion,p_legends+theme(plot.margin = margin(t=-3,unit = 'cm')),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,h_study,4,0.05,7,1))
  }else{
    pall <- plot_grid(p_cluster,p_purity,p_mutation,p_cosine,p_proportion,p_legends+theme(plot.margin = margin(t=-3,unit = 'cm')),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,h_purity,4,0.05,7,1))
  }
  
  if(is.null(output_plot)){
    return(pall)
  }else{
    xleng <- 2.5+dim(sigdata)[1]*0.1
    xleng <- if_else(xleng>20,20,xleng)
    yleng <- 12
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = pall,width = plot_width,height = plot_height)
  }
  
}




## piechar common

piechart_plot <- function(data,colset=NULL,keep_legend = TRUE,legend_name=NULL, output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  ## data format: Type,Catelogy,Freq,Label
  colnames(data) <- c("Type","Catelogy",'Frequency','Label')
  
  if(is.null(colset)){
    uvalues <- unique(data$Catelogy)
    if((unique(uvalues %in% names(Subscolor))) == TRUE && length(unique(uvalues %in% names(Subscolor))) == 1){
      sigcolorindex <- as.character(Subscolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else if((unique(uvalues %in% names(SBScolor))) == TRUE && length(unique(uvalues %in% names(SBScolor))) == 1){
      sigcolorindex <- as.character(SBScolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else {
      sigcolorindex <- as.character(SBScolor[1:length(uvalues)])
      names(sigcolorindex) <- uvalues
    }
  }else{
    sigcolorindex <- colset
  }
  
  
  # convert nsig_data to another format for the piechart
  data_pie <- data %>%
    group_by(Type) %>% 
    arrange(Frequency) %>%
    mutate(
      end_angle = 2*pi*cumsum(Frequency)/sum(Frequency),   # ending angle for each pie slice
      start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
      mid_angle = 0.5*(start_angle + end_angle),   # middle of each pie slice, for the text label
      # horizontal and vertical justifications depend on whether we're to the left/right
      # or top/bottom of the pie
      hjust = ifelse(mid_angle > pi, 1, 0),
      vjust = ifelse(mid_angle < pi/2 | mid_angle > 3*pi/2, 0, 1)
    ) %>% 
    ungroup()
  
  # radius of the pie and radius for outside and inside labels
  rpie <- 1
  rlabel_out <- 1.05 * rpie
  rlabel_in <- 0.6 * rpie
  legend_name <- if_else(is.null(legend_name),'',legend_name)
  
  p <- ggplot(data_pie) +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie, start = start_angle, end = end_angle, fill = Catelogy)) +
    geom_text(aes(x = rlabel_in * sin(mid_angle), y = rlabel_in * cos(mid_angle), label = Label ), size = 4)+
    facet_wrap(~Type,nrow = 2)+
    coord_fixed()+
    labs(x="",y="")+
    scale_fill_manual(legend_name,values = sigcolorindex)+
    theme_ipsum_rc(axis = FALSE, grid = FALSE)+
    theme(axis.text.y = element_blank(),axis.text.x = element_blank(),strip.text = element_text(size = 14,hjust = 0.5,face = "bold"))
  
  if(!(keep_legend)){
    p <- p+theme(legend.position = 'none')
  }
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 7
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
} 



## bar plot 
barchart_plot <- function(data,colset=NULL,keep_legend = TRUE,legend_name=NULL, output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  ## data format: Type,Catelogy,Freq,Label
  colnames(data) <- c("Type","Catelogy",'Frequency')
  
  if(is.null(colset)){
    uvalues <- unique(data$Catelogy)
    if((unique(uvalues %in% names(Subscolor))) == TRUE && length(unique(uvalues %in% names(Subscolor))) == 1){
      sigcolorindex <- as.character(Subscolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else if((unique(uvalues %in% names(SBScolor))) == TRUE && length(unique(uvalues %in% names(SBScolor))) == 1){
      sigcolorindex <- as.character(SBScolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else {
      sigcolorindex <- as.character(SBScolor[1:length(uvalues)])
      names(sigcolorindex) <- uvalues
    }
  }else{
    sigcolorindex <- colset
  }
  legend_name <- if_else(is.null(legend_name),'',legend_name)
  
  p <- data %>% 
    mutate(Catelogy=fct_reorder(Catelogy,Frequency,.desc = TRUE)) %>% 
    ggplot(aes(Catelogy,Frequency,fill=Catelogy))+
    geom_col(width = 0.8)+
    geom_text(aes(label = percent(Frequency,accuracy = 0.1)), vjust = -.5,size=3.5)+
    facet_wrap(~Type,nrow = 2)+
    labs(x="",y="Frequency (%)")+
    scale_y_percent(breaks = pretty_breaks(),limits = c(0,1),expand = expand_scale(mult = 0.1))+
    scale_fill_manual(legend_name,values = sigcolorindex,breaks = names(sigcolorindex))+
    theme_ipsum_rc(axis_title_size = 12,axis_title_just = "m",axis = TRUE, grid = "Yy")+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),strip.text = element_text(size = 14,hjust = 0.5,face = "bold"),axis.line.x = element_line(colour = 'black',size = 0.25),axis.line.y = element_line(colour = 'black',size = 0.25),legend.position = 'bottom')+
    guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  
  if(!(keep_legend)){
    p <- p+theme(legend.position = 'none')
  }
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 2.5+length(uvalues)*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
} 

barchart_plot2 <- function(data, output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  ## data format: Type,Catelogy,Freq,Label
  colnames(data) <- c("Type","Catelogy",'Frequency')
  uvalues <- sort(unique(data$Catelogy))
  p <- data %>% 
    mutate(Catelogy=fct_reorder(Catelogy,Frequency,.desc = TRUE)) %>% 
    ggplot(aes(Catelogy,Frequency,fill=Frequency))+
    geom_col(width = 0.8)+
    #geom_text(aes(label = percent(Frequency,accuracy = 0.1)), vjust = -.5,size=3.5)+
    facet_wrap(~Type,nrow = 2)+
    labs(x="",y="Frequency")+
    scale_y_continuous(breaks = pretty_breaks(),expand = expansion(mult = c(0,0.1)))+
    #scale_fill_manual(legend_name,values = sigcolorindex,breaks = names(sigcolorindex))+
    theme_ipsum_rc(axis_title_size = 12,axis_title_just = "m",axis = TRUE, grid = "Yy",)+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = 8),strip.text = element_text(size = 14,hjust = 0.5,face = "bold"),axis.line.x = element_line(colour = 'black',size = 0.25),axis.line.y = element_line(colour = 'black',size = 0.25),legend.position = 'none')
  #    guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 2.5+length(uvalues)*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- 4
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
} 



prevalence_plot <- function(sigdata,nmutation = 0, legend_name="Sigantures", output_plot = NULL,plot_width=NULL, plot_height=NULL){
  pie_input <- sigdata %>% 
    summarise(across(where(is.numeric),~sum(.x,na.rm=TRUE))) %>% 
    mutate(Type="Prevalence by mutations") %>%
    select(Type,everything()) %>%
    adorn_percentages(denominator = 'row') %>% 
    pivot_longer(cols = -Type) %>% 
    mutate(lab=if_else(value>0.05,percent(value,accuracy = 0.1),''))
  
  bar_input <- sigdata %>% pivot_longer(cols=-Samples) %>% #,names_transform = list(key = forcats::fct_inorder)
    filter(value>nmutation) %>%
    count(name) %>%
    mutate(value=n/dim(sigdata)[1]) %>% 
    mutate(Type="Prevalence by samples") %>% 
    select(Type,Catelogy=name,Frequency=value) %>% 
    mutate(Catelogy=factor(Catelogy,levels = c(colnames(sigdata)[-1]))) %>% 
    arrange(Catelogy) %>% 
    mutate(Catelogy=as.character(Catelogy))
  
  p_piechart <- piechart_plot(data = pie_input,keep_legend = FALSE)
  p_barchart <- barchart_plot(data = bar_input,legend_name = legend_name,keep_legend = FALSE)
  
  uvalues <- sort(unique(bar_input$Catelogy))
  pall <- plot_grid(p_piechart+theme(plot.margin = margin(r = -2)),p_barchart,align = 'h',nrow = 1,rel_widths = c(3,length(uvalues)))
  
  if(is.null(output_plot)){
    return(pall)
  }else{
    uvalues <- sort(unique(bar_input$Catelogy))
    xleng <- 6+length(uvalues)*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- 5
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = pall,width = plot_width,height = plot_height)
  }
  
}



iupac <- function(base){
  base2 <- NA_character_
  if(base == "A"){base2 = "A"}
  if(base == "T"){base2 = "T"}
  if(base == "C"){base2 = "C"}
  if(base == "G"){base2 = "G"}
  if(base == "R"){base2 = c("A","G")}
  if(base == "Y"){base2 = c("C", "T")}
  if(base == "S"){base2 = c("G","C")}
  if(base == "W"){base2 = c("A","T")}
  if(base == "K"){base2 = c("G","T")}
  if(base == "M"){base2 = c("A","C")}
  if(base == "B"){base2 = c("C","G","T")}
  if(base == "D"){base2 = c("A","G","T")}
  if(base == "H"){base2 = c("A","C","T")}
  if(base == "V"){base2 = c("A","C","G")}
  if(base == "N"){base2 = c("A","T","C","G")}
  return(base2)
}

iupac_list <- c('A','T','C','G','R','Y','S','W','K','M','B','D','H','V','N')

context_plot <- function(data,pattern,data_return = FALSE,output_plot = NULL,plot_width=14, plot_height=9){
  # type <- 'C>T'
  # subtype1 <- 'N'
  # subtype2 <- 'G'
  #pattern <- paste0(subtype1,str_sub(type,1,1),subtype2,">",subtype1,str_sub(type,3,3),subtype2)
  type <- paste0(str_sub(pattern,2,2),str_sub(pattern,4,4),str_sub(pattern,6,6))
  subtype1 <- str_sub(pattern,1,1)
  subtype2 <- str_sub(pattern,3,3)
  pattern1 <- paste0("\n",subtype1,str_sub(type,1,1),subtype2,">",subtype1,str_sub(type,3,3),subtype2," context")
  pattern2 <- paste0(type," other contexts\n")
  
  tmpdata0 <- suppressMessages(data %>% group_by(Study,Sample) %>% summarise(Total=sum(Mutations,na.rm = TRUE)) %>% ungroup())
  tmpdata1 <- suppressMessages(data %>% filter(Type == type) %>% group_by(Study,Sample,Type) %>% summarise(N0=sum(Mutations,na.rm = TRUE)) %>% ungroup())
  tmpdata2 <- suppressMessages(data %>% filter(Type == type,SubType1 %in% iupac(subtype1),SubType2 %in% iupac(subtype2)) %>% group_by(Study,Sample,Type) %>% summarise(N1=sum(Mutations,na.rm = TRUE)) %>% ungroup())
  
  data <- suppressMessages(left_join(tmpdata0,tmpdata1) %>% 
                             left_join(tmpdata2) %>% 
                             mutate(N2=N0-N1) %>% 
                             mutate(N1=N1/Total,N2=N2/Total,Type=pattern) %>% 
                             rename(Pattern=Type) %>% 
                             filter(Total>200))
  
  if(data_return){
    return(data)
  }
  p <- data %>% 
    ggplot(aes(N1,N2,size=(Total),fill=Study))+
    geom_point(pch=21,stroke=0.25,col="black")+
    labs(x=pattern1,y=pattern2,size="Number of mutations")+
    #guides(fill=FALSE)+
    guides(fill = guide_legend(override.aes = list(size=4)))+
    xlim(c(0,1))+
    ylim(c(0,1))+
    scale_fill_manual(values = as.character(SBScolor))+
    theme_ipsum_rc(axis_title_size = 14,axis_title_just = "m",axis = TRUE, grid = "XxYy")+
    theme(axis.line.x = element_line(colour = 'black',size = 0),axis.line.y = element_line(colour = 'black',size = 0),legend.position = 'right')+
    scale_size_continuous(trans = "log10",labels = comma_format())
  
  if(is.null(output_plot)){
    return(p)
  }else{
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
}


content_extraction <- function(data){
  content_data_all <- tibble(Study=character(),Sample=character(),Total=numeric(),Pattern=character(),N0=numeric(),N1=numeric(),N2=numeric())
  for(i in c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")){
    for(j1 in iupac_list){
      for(j2 in iupac_list){
        if(j1 == j2 & j1== "N"){ next }
        pattern <- paste0(j1,str_sub(i,1,1),j2,">",j1,str_sub(i,3,3),j2)
        tmp <- context_plot(data = data,pattern = pattern,data_return = TRUE)
        content_data_all <- bind_rows(content_data_all,tmp)
      }
    }
  }
  return(content_data_all)
}



signature_association <- function(data,signature_name_input1="Siganture 1",signature_name_input2="Siganture 2", signature_both=FALSE,output_plot = NULL,plot_width=10, plot_height=10){
  # cancer_type_only
  # if(!is.null(cancer_type_input)){
  #   data <- data %>% filter(Cancer_Type==cancer_type_input)
  # }
  
  if(signature_both){
    data <- data %>% filter(Exposure1>0,Exposure2>0)
  }
  
  ntotal <- data %>% filter(!is.na(Exposure1),!is.na(Exposure2)) %>% dim() %>% .[[1]]
  
  if(ntotal==0){
    errinfo <- "ERROR: No sample assigned to both signatures in selected study"
    return(errinfo)
  }else {
    
    # data %>% 
    #   ggplot(aes(log10(Exposure1+1),log10(Exposure2+1)))+
    #   geom_point()+
    #   geom_smooth(method = "lm",se = TRUE)+
    #   labs(x=paste0('Nubmer of mutations in ',signature_name_input1, ' (log10)'),y=paste0('Nubmer of mutations in ',signature_name_input2,' (log10)'))+
    #   theme_ipsum_rc(axis_title_size = 12,axis_title_just = 'm',axis_col = 'black',ticks = T)
    #   
    
    ## generated another Rplot blank file???
    p <- ggstatsplot::ggscatterstats(
      data=data %>% mutate(Exposure1=log10(Exposure1+1),Exposure2=log10(Exposure2+1)),
      x=Exposure1,
      y=Exposure2,
      xlab=paste0('Nubmer of mutations in ',signature_name_input1, ' (log10)'),
      ylab=paste0('Nubmer of mutations in ',signature_name_input2,' (log10)'),
      marginal.type = "density",
      messages=FALSE,
    )
    
    if(is.null(output_plot)){
      return(p)
    }else{
      ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
    }
    
  }
  
  
}


# Genome2Size  ------------------------------------------------------------
genome2size <- function(genome){
  genomesize <- case_when(
    genome == "hg18" ~ 3080436051/10^6, 
    genome == "GRCh38" ~ 3217346917/10^6, 
    genome == "GRCh37" ~ 3101976562/10^6,
    genome == "hg19" ~ 3101976562/10^6,
    genome == "mmc10" ~ 2725537669/10^6,
    genome == "mmc9" ~ 2654911517/10^6,
    TRUE ~ NA_real_
  )
  return(genomesize)
}




# Rainfall_Plot -----------------------------------------------------------

# mutSNP ----------------------------------------------------------------
mutSNP.input <- function(mut.data, chr = "chr", pos = "pos", ref = "ref", 
                         alt = "alt", build = NULL, k = 10) {
  if (exists("mut.data", mode = "list")) {
    mut.data <- mut.data
  } else {
    if (file.exists(mut.data)) {
      mut.data <- utils::read.table(mut.data, sep = "\t", header = TRUE, 
                                    as.is = FALSE, check.names = FALSE)
    } else {
      stop("mut.data is neither a file nor a loaded data frame")
    }
  }
  mut.data <- mut.data[, c(chr, pos, ref, alt)]
  genome.opts = c("hg19", "hg18", "hg38")
  if (!build %in% genome.opts || is.null(build)) {
    stop("Available reference builds: hg18, hg19, hg38")
  }
  if (build == "hg19") {
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 
                 171115067, 159138663, 146364022, 141213431, 135534747, 
                 135006516, 133851895, 115169878, 107349540, 102531392, 
                 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 
                 51304566, 155270560, 59373566)
    bsg = BSgenome.Hsapiens.UCSC.hg19
  } else if (build == "hg18") {
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 
                 170899992, 158821424, 146274826, 140273252, 135374737, 
                 134452384, 132349534, 114142980, 106368585, 100338915, 
                 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 
                 49691432, 154913754, 57772954)
    bsg = BSgenome.Hsapiens.UCSC.hg18
  } else if (build == "hg38") {
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 
                 170805979, 159345973, 145138636, 138394717, 133797422, 
                 135086622, 133275309, 114364328, 107043718, 101991189, 
                 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 
                 50818468, 156040895, 57227415)
    bsg = BSgenome.Hsapiens.UCSC.hg38
  } else {
    stop("Available reference builds: hg18, hg19, hg38")
  }
  mut.data$build = build
  if (!all(mut.data$ref %in% DNA_BASES & mut.data$alt %in% DNA_BASES)) {
    stop("Only SNV substitutions are currently supported.")
  }
  ref_base = DNAStringSet(mut.data$ref)
  alt_base = DNAStringSet(mut.data$alt)
  conv.start = mut.data$pos - k
  conv.end = mut.data$pos + k
  context = getSeq(bsg, mut.data$chr, start = conv.start, end = conv.end)
  if (TRUE) {
    idx = mut.data$ref %in% c("A", "G")
    context[idx] = reverseComplement(context[idx])
    ref_base[idx] = reverseComplement(ref_base[idx])
    alt_base[idx] = reverseComplement(alt_base[idx])
  }
  mut.data$alteration = paste(ref_base, alt_base, sep = ">")
  mut.data$context = context
  # Replace chr x and y with numeric value (23 and 24) for better
  # ordering
  seq = gsub(pattern = "chr", replacement = "", x = mut.data$chr, 
             fixed = TRUE)
  seq = gsub(pattern = "X", replacement = "23", x = seq, fixed = TRUE)
  seq = gsub(pattern = "Y", replacement = "24", x = seq, fixed = TRUE)
  mut.data$seq = as.numeric(seq)
  mut.data = mut.data[order(mut.data$seq, mut.data$pos), ]
  chr.lens.sum = cumsum(chr.lens)
  chr.lens.sum = c(0, chr.lens.sum)
  mut.data$dis = c(mut.data$pos[1], diff(mut.data$pos + chr.lens.sum[mut.data$seq]))
  return(mut.data)
}




# Kataegis point ----------------------------------------------------------
katPoint <- function(data, sample = "sample", min.mut = 5, max.dis = 1000, 
                     txdb = NULL) {
  build = data$build[1]
  genome.opts = c("hg19", "hg18", "hg38")
  if (!build %in% genome.opts) {
    stop("Available reference builds: hg18, hg19, hg38")
  }
  if (build == "hg19") {
    chr.arm = c(1.25e+08, 93300000, 9.1e+07, 50400000, 48400000, 
                6.1e+07, 59900000, 45600000, 4.9e+07, 40200000, 53700000, 
                35800000, 17900000, 17600000, 1.9e+07, 36600000, 2.4e+07, 
                17200000, 26500000, 27500000, 13200000, 14700000, 60600000, 
                12500000)
  } else if (build == "hg18") {
    chr.arm = c(124300000, 93300000, 91700000, 50700000, 47700000, 
                60500000, 59100000, 45200000, 51800000, 40300000, 52900000, 
                35400000, 1.6e+07, 15600000, 1.7e+07, 38200000, 22200000, 
                16100000, 28500000, 27100000, 12300000, 11800000, 59500000, 
                11300000)
  } else if (build == "hg38") {
    chr.arm = c(123400000, 93900000, 90900000, 5e+07, 48800000, 
                59800000, 60100000, 45200000, 4.3e+07, 39800000, 53400000, 
                35500000, 17700000, 17200000, 1.9e+07, 36800000, 25100000, 
                18500000, 26200000, 28100000, 1.2e+07, 1.5e+07, 6.1e+07, 
                10400000)
  } else {
    stop("Available reference builds: hg18, hg19, hg38")
  }
  num = dim(data)[1] - 5
  katPoint <- matrix(nrow = num, ncol = 8)
  i = 1
  mutnum = 1
  Cmutnum = 0
  for (i in 1:num) {
    if (data$ref[i] %in% c("C", "G")){
      Cmutnum = Cmutnum + 1
    }
    if (data$dis[i + 1] <= max.dis) {
      mutnum = mutnum + 1
    } else {
      if (mutnum >= min.mut) {
        len = data$pos[i] - data$pos[i - mutnum + 1] + 1
        chr.n = gsub(pattern = "chr", replacement = "", x = data$chr[i], 
                     fixed = TRUE)
        chr.n = gsub(pattern = "X", replacement = "23", x = chr.n, 
                     fixed = TRUE)
        chr.n = gsub(pattern = "Y", replacement = "24", x = chr.n, 
                     fixed = TRUE)
        chr.n = as.numeric(chr.n)
        if (data$pos[i] <= chr.arm[chr.n]) {
          arm = paste(chr.n, "p", sep = "")
        } else if (data$pos[i - mutnum + 1] >= chr.arm[chr.n]) {
          arm = paste(chr.n, "q", sep = "")
        } else {
          arm = paste(chr.n, "p, ", chr.n, "q", sep = "")
        }
        katPoint[i, 1:8] = c(sample, data$chr[i], data$pos[i - mutnum + 1], data$pos[i], arm, len, mutnum, round(Cmutnum/mutnum,3))
      }
      mutnum = 1
      Cmutnum = 0
    }
  }
  katPoint.out = data.frame(na.omit(katPoint))
  names(katPoint.out) = c("sample", "chrom", "start", "end", "chrom.arm", "length", "number.mut", 
                          "weight.C>X")
  for (i in 1:dim(katPoint.out)[1]) {
    if (as.numeric(as.character(katPoint.out$"weight.C>X"[i])) < 0.8) {
      katPoint.out$confidence[i] = 0
    } else {
      katPoint.out$confidence[i] <- length(which(subset(katPoint.out,as.numeric(as.character(katPoint.out$"weight.C>X")) >= 0.8)$chrom == katPoint.out$chrom[i]))
      if (katPoint.out$confidence[i] > 3) {
        katPoint.out$confidence[i] = 3
      }
    }
  }
  if (!is.null(txdb)) {
    gr <- GRanges(seqnames = Rle(katPoint.out$chrom), ranges=IRanges(start = as.numeric(as.character(katPoint.out$start)), end =as.numeric(as.character(katPoint.out$end))))
    peakAnno <- ChIPseeker::annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
    katPoint.out$annotation <- peakAnno@anno$annotation
    katPoint.out$distanceToTSS <- peakAnno@anno$distanceToTSS
    katPoint.out$geneName <- peakAnno@anno$SYMBOL
    katPoint.out$geneID <- peakAnno@anno$geneId
  } 
  message(paste(dim(katPoint.out)[1], "potential kataegis events identified", 
                sep = " "))
  return(katPoint.out)
}


# Kataegis plot -----------------------------------------------------------

kataegis_rainfall_plot <- function(mutdata,sample_name="sample",genome_build = "hg19",reference_data_folder, chromsome=NULL,kataegis_highligh=FALSE,min.mut = 5,max.dis = 1000,filename=NULL){
  require(tidyverse)
  require(ggsci)
  require(hrbrthemes)
  require(cowplot)
  #require(ChIPseeker)
  
  genome_build <- str_to_lower(genome_build)
  
  if(genome_build %in% c('hg19','grch37')){
    require(BSgenome.Hsapiens.UCSC.hg19)
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    genome_build <- 'hg19'
    ref_file <- paste0(reference_data_folder,"/hg19_ref.RData")
    load(ref_file)
    hgref <- hg19
    TxDb.Hsapiens <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }
  
  if(genome_build %in% c('hg38','grch38')){
    require(BSgenome.Hsapiens.UCSC.hg38)
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    genome_build <- 'hg38'
    ref_file <- paste0(reference_data_folder,"/hg38_ref.RData")
    load(ref_file)
    hgref <- hg38
    TxDb.Hsapiens <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  
  
  #mutdata format:chr,pos,ref,alt
  mutdata <- mutdata %>% 
    mutate(chr=str_remove(chr,"^chr")) %>% 
    dplyr::filter(str_length(ref)==1,str_length(alt)==1,!str_detect(alt,"-"),!str_detect(ref,"-"),chr %in% c(1:22,"X","Y"))%>% 
    mutate(chr=paste0("chr",chr)) 
  
  if(!is.null(chromsome)){
    hgref <- hgref %>% dplyr::filter(chr==chromsome) %>% mutate(start2=1,end2=len)
    mutdata <- mutdata %>% dplyr::filter(chr==chromsome)
  }
  
  mutSNP = mutSNP.input(mut.data = as.data.frame(mutdata), chr = "chr", pos = "pos", ref = "ref", alt = "alt", build = genome_build)
  mutSNP$context <- ""
  katdata <- NA
  try({
    katdata <- katPoint(mutSNP,txdb = TxDb.Hsapiens,sample = sample_name,min.mut = min.mut,max.dis = max.dis)
  },silent = TRUE)
  #print(katdata)
  
  if(is.data.frame(katdata)){
    if(dim(katdata)[1]>0){
      katdata <- katdata %>% mutate(start=as.integer(start),end=as.integer(end)) %>% left_join(hgref,by=c('chrom'='chr'))
      for(i in 1:dim(katdata)[1]){
        kchr <- katdata$chrom[i]
        kstart <- katdata$start[i]
        kend <- katdata$end[i]
        mutSNP$context[mutSNP$chr==kchr & mutSNP$pos>=kstart & mutSNP$pos<=kend] <- 'Kataegsis'
      }
    }
  }
  
  #color themes for the mutation subtype
  SNVcolor <- pal_d3()(6)
  names(SNVcolor) = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  addcolor <- "gray80"
  names(addcolor) <- "Non-kataegis mutations"
  
  mutSNP <- mutSNP %>% 
    left_join(hgref %>% dplyr::select(chr,start2),by=c('chr'='chr')) %>% 
    mutate(pos=pos+start2-1) %>% 
    dplyr::rename(dist=dis) 
  
  if(kataegis_highligh){
    mutSNP <- mutSNP %>% mutate(SNV=if_else(context=="","Non-kataegis mutations",alteration))
    SNVcolor <- c(SNVcolor,addcolor)
  }else{
    mutSNP <- mutSNP %>% mutate(SNV=alteration) 
  }
  
  mutSNP <- mutSNP %>% mutate(SNV=factor(SNV,level=names(SNVcolor))) 
  
  p <- mutSNP %>% 
    ggplot(aes(pos,log10(dist),fill=SNV))+
    geom_vline(xintercept = c(hgref$start2,hgref$end2),col="#cccccc",lwd=0.4)+
    geom_hline(yintercept = c(2,3),col=c("red","blue"),lwd=0.4,linetype=2)+
    geom_point(shape=21,stroke=0.1,size=1.5)+
    labs(y="Intermutation distance (bp,log10)\n",x = "Chromosome")+
    scale_x_continuous(breaks = (hgref$start2+hgref$end2)/2,labels = str_remove(hgref$chr,'chr'),expand = c(0.005,0.005),limits = c(0,tail(hgref$end2,1)),guide = guide_axis(n.dodge=1))+
    scale_y_continuous(breaks = pretty_breaks())+
    scale_fill_manual("Subtype",values = SNVcolor,drop=FALSE,na.value="gray90")+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 14,grid="Y",ticks = T,axis = FALSE)+
    theme(plot.title = element_text(hjust = 0.5,face = "plain"),legend.text = element_text(size = 12),legend.position = "bottom")+
    guides(fill = guide_legend(override.aes = list(size = 3),nrow = 1))+
    panel_border(color = 'black',size = 0.5)+
    coord_cartesian(clip = 'off')
  
  if(is.data.frame(katdata)){
    if(dim(katdata)[1]>0){
      p <- p+annotate("segment", x = katdata$start2+(katdata$start+katdata$end)/2, xend = katdata$start2+(katdata$start+katdata$end)/2, y = 0.5, yend = 1,size=0.3, colour = "black",arrow=arrow(length = unit(0.25,"cm")))
      katdata <- katdata %>% dplyr::select(sample:geneID)
    }}
  ggsave(filename = filename,plot = p,width = 15,height = 7)
  
  return(katdata)
}






### Association function ### 


### other funcitons ###
rowAny <- function(x) rowSums(x) > 0
