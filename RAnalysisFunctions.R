library(ggplot2)
library(ggpubr)
library(Biostrings)
library(ggseqlogo)
library(seqinr)
library(dplyr)

###################################### Functions ###############################

get_commonSites <- function(common_signal_t1_path, common_signal_t2_path){
  
  common_signal_t1 <- read.table(common_signal_t1_path,header=T,sep=",")
  if(!("gene_name" %in% colnames(common_signal_t1))){
    common_signal_t1$gene_name <- "NONE"
  }
  common_signal_t1$diff_deletion_ratio <- common_signal_t1$treated_deletion_ratio - common_signal_t1$ctrl_deletion_ratio
  common_signal_t1$start <- unlist(lapply(strsplit(common_signal_t1$site,"-"),function(x) x[1]))
  common_signal_t1$end <- unlist(lapply(strsplit(common_signal_t1$site,"-"),function(x) x[2]))
  common_signal_t1[which(is.na(common_signal_t1$end)),"end"] <- common_signal_t1[which(is.na(common_signal_t1$end)),"start"]
  colnames(common_signal_t1) <- c("chr_name","site","treated_total_counts_t1",
                                  "treated_deletion_counts_t1","treated_deletion_ratio_t1",
                                  "ctrl_total_counts_t1","ctrl_deletion_counts_t1",
                                  "ctrl_deletion_ratio_t1","p_value_t1","gene_name",
                                  "diff_deletion_ratio_t1","start","end")
  
  common_signal_t2 <- read.table(common_signal_t2_path,header=T,sep=",")
  if(!("gene_name" %in% colnames(common_signal_t2))){
    common_signal_t2$gene_name <- "NONE"
  }
  common_signal_t2$diff_deletion_ratio <- common_signal_t2$treated_deletion_ratio - common_signal_t2$ctrl_deletion_ratio
  common_signal_t2$start <- unlist(lapply(strsplit(common_signal_t2$site,"-"),function(x) x[1]))
  common_signal_t2$end <- unlist(lapply(strsplit(common_signal_t2$site,"-"),function(x) x[2]))
  common_signal_t2[which(is.na(common_signal_t2$end)),"end"] <- common_signal_t2[which(is.na(common_signal_t2$end)),"start"]
  colnames(common_signal_t2) <- c("chr_name","site","treated_total_counts_t2",
                                  "treated_deletion_counts_t2","treated_deletion_ratio_t2",
                                  "ctrl_total_counts_t2","ctrl_deletion_counts_t2",
                                  "ctrl_deletion_ratio_t2","p_value_t2","gene_name",
                                  "diff_deletion_ratio_t2","start","end")
  
  common_signal <- inner_join(common_signal_t1, common_signal_t2,
                              by = c("chr_name","site","gene_name","start","end"))
  common_signal$diff_deletion_ratio_mean <- 
    (common_signal$diff_deletion_ratio_t1 + common_signal$diff_deletion_ratio_t2)/2
  
  common_signal_list <- common_signal[,c("chr_name","site","end","diff_deletion_ratio_mean")]
  colnames(common_signal_list) <- c("chr_name","site","chr_index","diff_deletion_ratio_mean.control")
  common_signal_list$chr_index <- as.numeric(common_signal_list$chr_index)
  common_signal_list$label <- "psiU"
  
  return(list(common_signal, common_signal_list))
  
}

get_totalSignal <- function(total_signal_ut_path, total_signal_t1_path, total_signal_t2_path){
  
  total_signal_ut <- read.table(total_signal_ut_path,header=T,sep="\t")
  total_signal_t1 <- read.table(total_signal_t1_path,header=T,sep="\t")
  total_signal_t2 <- read.table(total_signal_t2_path,header=T,sep="\t")
  
  total_signal_ut$total_count <-  total_signal_ut$A + total_signal_ut$G + 
    total_signal_ut$C + total_signal_ut$T + total_signal_ut$del_._count
  total_signal_ut$deletion_ratio <- total_signal_ut$del_._count/total_signal_ut$total_count
  colnames(total_signal_ut) <- c("chr_name","chr_index","ref_base","A_ut","G_ut","C_ut","T_ut",
                                 "del_start_count_ut","del_._count_ut","mis_start_count_ut",
                                 "total_count_ut","deletion_ratio_ut")
  
  total_signal_t1$total_count <-  total_signal_t1$A + total_signal_t1$G + 
    total_signal_t1$C + total_signal_t1$T + total_signal_t1$del_._count
  total_signal_t1$deletion_ratio <- total_signal_t1$del_._count/total_signal_t1$total_count
  colnames(total_signal_t1) <- c("chr_name","chr_index","ref_base","A_t1","G_t1","C_t1","T_t1",
                                 "del_start_count_t1","del_._count_t1","mis_start_count_t1",
                                 "total_count_t1","deletion_ratio_t1")
  
  total_signal_t2$total_count <-  total_signal_t2$A + total_signal_t2$G + 
    total_signal_t2$C + total_signal_t2$T + total_signal_t2$del_._count
  total_signal_t2$deletion_ratio <- total_signal_t2$del_._count/total_signal_t2$total_count
  colnames(total_signal_t2) <- c("chr_name","chr_index","ref_base","A_t2","G_t2","C_t2","T_t2",
                                 "del_start_count_t2","del_._count_t2","mis_start_count_t2",
                                 "total_count_t2","deletion_ratio_t2")
  
  total_signal_t1_merge <- inner_join(total_signal_ut,total_signal_t1,by=c("chr_name","chr_index","ref_base"))
  total_signal_t1_merge$diff_deletion_ratio_t1 <- total_signal_t1_merge$deletion_ratio_t1 - total_signal_t1_merge$deletion_ratio_ut
  total_signal_t2_merge <- inner_join(total_signal_ut,total_signal_t2,by=c("chr_name","chr_index","ref_base"))
  total_signal_t2_merge$diff_deletion_ratio_t2 <- total_signal_t2_merge$deletion_ratio_t2 - total_signal_t2_merge$deletion_ratio_ut
  
  total_signal_merge <- inner_join(total_signal_t1_merge,total_signal_t2_merge,
                                   by=c("chr_name","chr_index","ref_base",
                                        "A_ut","G_ut","C_ut","T_ut",
                                        "del_start_count_ut","del_._count_ut","mis_start_count_ut",
                                        "total_count_ut","deletion_ratio_ut"))
  total_signal_merge$diff_deletion_ratio_mean.treated <- 
    (total_signal_merge$diff_deletion_ratio_t1 + total_signal_merge$diff_deletion_ratio_t2)/2
  
  total_signal_list <- total_signal_merge[,c("chr_name","chr_index","diff_deletion_ratio_mean.treated")]
  
  return(list(total_signal_merge, total_signal_list))
  
}

get_changedSites <- function(common_signal_list, total_signal_list, 
                             absolute_threld=0.05, fold_threld=0.3){
  
  common_signal_withLabel <- inner_join(common_signal_list,total_signal_list,
                                        by =c("chr_name","chr_index"))
  common_signal_withLabel$ratioChange_afterWriterKO <- 
    common_signal_withLabel$diff_deletion_ratio_mean.treated - 
    common_signal_withLabel$diff_deletion_ratio_mean.control
  common_signal_withLabel <- common_signal_withLabel[order(common_signal_withLabel$ratioChange_afterWriterKO),]
  
  common_signal_withLabel$absoluteChangeValue_label <- "filter"
  common_signal_withLabel[which(abs(common_signal_withLabel$ratioChange_afterWriterKO) >= absolute_threld),"absoluteChangeValue_label"] <- "nonfilter"
  
  common_signal_withLabel$changeFold_afterWriterKO <- abs(common_signal_withLabel$ratioChange_afterWriterKO) / common_signal_withLabel$diff_deletion_ratio_mean.control
  common_signal_withLabel$relativeChangeValue_label <- "keep"
  down_label_vector_1 <- which(common_signal_withLabel$ratioChange_afterWriterKO<=0)
  down_label_vector_2 <- which(common_signal_withLabel$changeFold_afterWriterKO>=fold_threld)
  down_label_vector <- intersect(down_label_vector_1, down_label_vector_2)
  common_signal_withLabel[down_label_vector,"relativeChangeValue_label"] <- "down"
  up_label_vector_1 <- which(common_signal_withLabel$ratioChange_afterWriterKO>0)
  up_label_vector_2 <- which(common_signal_withLabel$changeFold_afterWriterKO>=fold_threld)
  up_label_vector <- intersect(up_label_vector_1, up_label_vector_2)
  common_signal_withLabel[up_label_vector,"relativeChangeValue_label"] <- "up"
  
  return(common_signal_withLabel)
}

get_dependentSites <- function(common_signal_list, total_signal_merge.writer_ko, 
                               absolute_threld = 0.05, relative_threld = 0.3){
  
  common_signal.writer_ko <- inner_join(common_signal_list, total_signal_merge.writer_ko, by = c("chr_name","chr_index"))
  common_signal.writer_ko$changeRatio_afterKO <- common_signal.writer_ko$diff_deletion_ratio_mean.treated - common_signal.writer_ko$diff_deletion_ratio_mean.control
  common_signal.writer_ko$changeFold_afterKO <- abs(common_signal.writer_ko$changeRatio_afterKO) / common_signal.writer_ko$diff_deletion_ratio_mean.control
  
  common_signal.writer_ko$dependentSitesLabel <- "non_dependentSites"
  dependentSites_index <- which(common_signal.writer_ko$changeRatio_afterKO <= -absolute_threld & common_signal.writer_ko$changeFold_afterKO >= relative_threld)
  common_signal.writer_ko[dependentSites_index,"dependentSitesLabel"] <- "dependentSites"
  
  common_signal.dependentSitesJudgement <- common_signal.writer_ko[,c("chr_name","chr_index",
                                                                      "diff_deletion_ratio_mean.control", "diff_deletion_ratio_mean.treated",
                                                                      "changeRatio_afterKO","changeFold_afterKO","dependentSitesLabel")]
  return(common_signal.dependentSitesJudgement)
  
}

changedSites_correlationPlot <- function(pus_dependent_signal, fold_threld = 0.3, absolute_threld = 0.05){
  
  cor_plot <- 
    ggplot(pus_dependent_signal,aes(x=diff_deletion_ratio_mean.control,
                                    y=diff_deletion_ratio_mean.treated,
                                    color=relativeChangeValue_label))+ 
    geom_point(size=1,shape=16,colour="lightgray",alpha=1)+
    geom_text(data = pus_dependent_signal[pus_dependent_signal$relativeChangeValue_label == "up", ], 
              aes(label = relativeChangeValue_label), vjust = 0, hjust = 0) + 
    # geom_smooth(method = "lm",se = F,color = "black")+
    # stat_cor(data = common_signal_withLabel,method = "pearson",digits=4,size=5) +
    theme_classic() +
    coord_equal() + 
    scale_x_continuous(expand = c(0.05, 0.05)) +  
    scale_y_continuous(expand = c(0.05, 0.05)) +  
    ggtitle("PsiU Ratio") + 
    theme(axis.title.x = element_text(size=15,colour = "black"),
          axis.title.y = element_text(size=15,colour = "black"),
          axis.text = element_text(size=15,colour="black"),
          plot.title = element_text(face = "bold",size=15), 
          panel.border = element_rect(fill = NA, colour = "black", size=1), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    labs(x="Control",y="WriterKO") +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) + 
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
    # geom_pointdensity(adjust = 4) +
    geom_abline(slope = 1, intercept = -absolute_threld, linetype = "dashed", color = "grey") + 
    geom_abline(slope = 1, intercept = absolute_threld, linetype = "dashed", color = "grey") +
    #geom_point(aes(color = ifelse(abs(diff_deletion_ratio_mean.treated-diff_deletion_ratio_mean.control) <= distance, "InArea", "OutArea"))) +
    #scale_color_manual(values = c("InArea" = "lightgrey", "OutArea" = "black")) +
    geom_point(aes(color = ifelse(abs(diff_deletion_ratio_mean.treated-diff_deletion_ratio_mean.control) <= absolute_threld, "InArea", 
                                  as.character(relativeChangeValue_label)))) +
    scale_color_manual(values = c("InArea" = "#F0F0F0", 
                                  "keep" = "grey", 
                                  "up" = "red", 
                                  "down" = "blue")) +
    geom_abline(slope = 1, intercept = 0,  color = "grey", size = 0.5,  linetype = "dashed") + 
    geom_abline(slope = 1 - fold_threld, intercept = 0,  color = "grey", size = 0.5,  linetype = "dashed") +
    geom_abline(slope = 1 + fold_threld, intercept = 0,  color = "grey", size = 0.5,  linetype = "dashed") 
  # (y-x)/x=0.5 y=1.5x 
  # (x-y)/x=0.5 y=0.5x 
  
  return(cor_plot)
  
}

changedSites_correlationPlot <- function(pus_dependent_signal, fold_threld = 0.3, absolute_threld = 0.05){
  
  cor_plot_1 <- 
    ggplot(pus_dependent_signal,aes(x=diff_deletion_ratio_mean.control,
                                    y=diff_deletion_ratio_mean.treated,
                                    color=relativeChangeValue_label,
                                    shape = relativeChangeValue_label))+ 
    geom_point(size=1,shape=16,colour="lightgray",alpha=1)+
    #geom_text(data = pus_dependent_signal[pus_dependent_signal$relativeChangeValue_label == "up", ], 
    #          aes(label = relativeChangeValue_label), vjust = 0, hjust = 0) + 
    # geom_smooth(method = "lm",se = F,color = "black")+
    # stat_cor(data = common_signal_withLabel,method = "pearson",digits=4,size=5) +
    theme_classic() +
    coord_equal() + 
    scale_x_continuous(expand = c(0.05, 0.05)) +  
    scale_y_continuous(expand = c(0.05, 0.05)) +
    ggtitle("PsiU Ratio") + 
    theme(axis.title.x = element_text(size=15,colour = "black"),
          axis.title.y = element_text(size=15,colour = "black"),
          axis.text = element_text(size=15,colour="black"),
          plot.title = element_text(face = "bold",size=15), 
          panel.border = element_rect(fill = NA, colour = "black", size=1),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    labs(x="Control",y="WriterKO") +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) + 
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
    # geom_pointdensity(adjust = 4) +
    geom_abline(slope = 1, intercept = -absolute_threld, linetype = "dashed", color = "grey") + 
    geom_abline(slope = 1, intercept = absolute_threld, linetype = "dashed", color = "grey") +
    #geom_point(aes(color = ifelse(abs(diff_deletion_ratio_mean.treated-diff_deletion_ratio_mean.control) <= distance, "InArea", "OutArea"))) +
    #scale_color_manual(values = c("InArea" = "lightgrey", "OutArea" = "black")) +
    geom_point(aes(color = ifelse(abs(diff_deletion_ratio_mean.treated-diff_deletion_ratio_mean.control) <= absolute_threld, "InArea", 
                                  as.character(relativeChangeValue_label)))) +
    scale_color_manual(values = c("InArea" = "#F0F0F0", 
                                  "keep" = "grey", 
                                  "up" = "red", 
                                  "down" = "blue")) +
    geom_abline(slope = 1, intercept = 0,  color = "grey", size = 0.5,  linetype = "dashed") + 
    geom_abline(slope = 1 - fold_threld, intercept = 0,  color = "grey", size = 0.5,  linetype = "dashed") +
    geom_abline(slope = 1 + fold_threld, intercept = 0,  color = "grey", size = 0.5,  linetype = "dashed") +
    scale_shape_manual(values = c("up" = 17, "down" = 24, "keep" = 16, "InArea" = 16))
  # (y-x)/x=0.5 y=1.5x 
  # (x-y)/x=0.5 y=0.5x 
  
  cor_plot_2 <- 
    ggplot(pus_dependent_signal,aes(x=diff_deletion_ratio_mean.control,
                                    y=diff_deletion_ratio_mean.treated,
                                    color=relativeChangeValue_label))+ 
    geom_point(data = subset(pus_dependent_signal, abs(diff_deletion_ratio_mean.treated - diff_deletion_ratio_mean.control) <= absolute_threld),
               size = 2, alpha = 1, shape = 16) +
    geom_point(data = subset(pus_dependent_signal, relativeChangeValue_label == "keep"),
               size = 2, alpha = 1, shape = 16) +
    #geom_point(data = subset(pus_dependent_signal, relativeChangeValue_label == "up"),
    #           size = 3, alpha = 1, shape = 17) + 
    geom_point(data = subset(pus_dependent_signal, relativeChangeValue_label == "up"),
               size = 2, alpha = 1, shape = 16) +  
    geom_point(data = subset(pus_dependent_signal, relativeChangeValue_label == "down"),
               size = 3, alpha = 1, shape = 25) +  
    # geom_text(data = pus_dependent_signal[pus_dependent_signal$relativeChangeValue_label == "up", ], 
    #          aes(label = relativeChangeValue_label), vjust = 0, hjust = 0) +
    # geom_smooth(method = "lm",se = F,color = "black")+
    # stat_cor(data = common_signal_withLabel,method = "pearson",digits=4,size=5) +
    theme_classic() +
    coord_equal() + 
    scale_x_continuous(expand = c(0.05, 0.05)) + 
    scale_y_continuous(expand = c(0.05, 0.05)) +  
    ggtitle("PsiU Ratio") + 
    theme(axis.title.x = element_text(size=15,colour = "black"),
          axis.title.y = element_text(size=15,colour = "black"),
          axis.text = element_text(size=15,colour="black"),
          plot.title = element_text(face = "bold",size=15), 
          panel.border = element_rect(fill = NA, colour = "black", size=1), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    labs(x="Control",y="WriterKO") +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) + 
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
    # geom_pointdensity(adjust = 4) +
    geom_abline(slope = 1, intercept = -absolute_threld, linetype = "dashed", color = "grey") + 
    geom_abline(slope = 1, intercept = absolute_threld, linetype = "dashed", color = "grey") +
    #geom_point(aes(color = ifelse(abs(diff_deletion_ratio_mean.treated-diff_deletion_ratio_mean.control) <= distance, "InArea", "OutArea"))) +
    #scale_color_manual(values = c("InArea" = "lightgrey", "OutArea" = "black")) +
    geom_point(aes(color = ifelse(abs(diff_deletion_ratio_mean.treated-diff_deletion_ratio_mean.control) <= absolute_threld, "InArea", 
                                  as.character(relativeChangeValue_label)))) +
    scale_color_manual(values = c("InArea" = "#F0F0F0", 
                                  #"keep" = "grey", # 
                                  "keep" = "#F0F0F0",
                                  #"up" = "#B22222", 
                                  "up" = "#F0F0F0",
                                  "down" = "#000080")) +
    geom_abline(slope = 1, intercept = 0,  color = "grey", size = 0.5,  linetype = "dashed") + 
    geom_abline(slope = 1 - fold_threld, intercept = 0,  color = "grey", size = 0.5,  linetype = "dashed") +
    geom_abline(slope = 1 + fold_threld, intercept = 0,  color = "grey", size = 0.5,  linetype = "dashed") 
  
  
  return(list(cor_plot_1,cor_plot_2))
  
}

###################################### Functions ###############################


