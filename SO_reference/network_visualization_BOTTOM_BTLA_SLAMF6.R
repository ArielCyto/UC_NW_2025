#https://github.com/Cytoreason/analysis-p13/blob/b88dd33ae68cdadef2ec38aca3153c8e2886d2a7/notes/Pso_CriteriaRanking/Graph_Criteria_NW_v2.R

devtools::load_all("~/code/p03-POC1-analysis/R")
dir.create("/home/coder/scratch", recursive = TRUE, showWarnings = FALSE)
options(
  cytoreason.cc_tempdir = "/home/coder/scratch/cytoreason_cache/cyto-cc",
  cytoreason.cc.client_tempdir = "/home/coder/scratch/cytoreason_cache/cyto-cc",
  # These ones are for good measures
  stringsAsFactors = FALSE, # always good to set to FALSE
  nwarnings = 10000 # capture all warnings,
)
library(cytoreason.ccm.pipeline)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(stringr)


#### Centrality ####

##### prepare data #####

# Define the Signature Table
SignatureTable <- read_asset("wf-b1136bff5e")

BulkAD_NW_cent <- readRDS(get_workflow_outputs("wf-7f3b0c0601"))
BulkAD_NW_cent_df <- BulkAD_NW_cent$FullParam
rm(BulkAD_NW_cent)
gc()

BulkAD_NW_cent_df_sub <- BulkAD_NW_cent_df[BulkAD_NW_cent_df$ListType %in% c("gene_list", paste("rand_iter", 1:1000, sep="_")),]

adj_AD_NW_cent <- readRDS(get_workflow_outputs("wf-cbd2fbb22f"))
adj_AD_NW_cent_df <- adj_AD_NW_cent$FullParam
rm(adj_AD_NW_cent)
gc()

adj_AD_NW_cent_df_sub <- adj_AD_NW_cent_df[adj_AD_NW_cent_df$ListType %in% c("gene_list", paste("rand_iter", 1:1000, sep="_")),]
rm(adj_AD_NW_cent_df)
gc()

BulkAD_NW_cent_df_sub$Type <- "bulk"
adj_AD_NW_cent_df_sub$Type <- "adjusted"

AD_NW_cent_df_sub <- rbind(BulkAD_NW_cent_df_sub, adj_AD_NW_cent_df_sub)

##### aggregate random #####

SignatureTable <- read_asset("wf-b1136bff5e")

#aggregate random and randomsmooth before scaling

AD_NW_cent_df_sub <-  AD_NW_cent_df_sub %>%
  dplyr::mutate(Identifier = ListName) %>%
  dplyr::left_join(SignatureTable)

mean_cent_df_random <- AD_NW_cent_df_sub %>%
  dplyr::mutate(Identifier = ListName) %>%
  dplyr::left_join(SignatureTable) %>%
  dplyr::filter(SubCategory %in% c("random", "smoothedRandom")) %>%
  dplyr::group_by(ListType,BroadCategory, SubCategory,Type) %>%
  dplyr::mutate(across(betweenness_mean:OverlapPercent, ~ ifelse(is.finite(.), ., NA))) %>%
  dplyr::summarise(across(betweenness_mean:OverlapPercent, mean, na.rm = TRUE))

mean_cent_df_random$ListName <- mean_cent_df_random$SubCategory
mean_cent_df_random$Identifier <- mean_cent_df_random$ListName
mean_cent_df_random$CleanName <- mean_cent_df_random$ListName

AD_NW_cent_df_sub$ModID <- NULL
AD_NW_cent_df_sub$Identifier_original <- NULL

##### combine dataframes #####

mean_cent_df_random <- mean_cent_df_random[,colnames(mean_cent_df_random) %in% colnames(AD_NW_cent_df_sub)]

AD_NW_cent_df_sub_combined <- rbind.fill(AD_NW_cent_df_sub[-grep("random|smoothedRandom", AD_NW_cent_df_sub$Identifier),] ,mean_cent_df_random )
cyto_cc_push_object(AD_NW_cent_df_sub_combined,"AD_NW_cent_df_sub_combined")
#wf-f463a895c7
#wf-c8c1f9bf0f #BOTTOM BTLA & SLAMF6

AD_NW_cent_df_sub_combined <- read_asset("wf-c8c1f9bf0f")
#remove largest component criterion

##### page rank gene_list vs random#####

AD_NW_cent_df_sub_combined_Target <- AD_NW_cent_df_sub_combined[!AD_NW_cent_df_sub_combined$Identifier %in% c("steadyArm__IL13.avgFDR","steadyArm__IL4.avgFDR"),]

AD_NW_cent_df_sub_combined <- AD_NW_cent_df_sub_combined[!grepl("3sd", AD_NW_cent_df_sub_combined$ListName),]
AD_NW_cent_df_sub_combined$Type <- factor(AD_NW_cent_df_sub_combined$Type)
AD_NW_cent_df_sub_combined <- AD_NW_cent_df_sub_combined[!grepl("avg", AD_NW_cent_df_sub_combined$SubCategory),]

#fix the factor (order) of the subcategories
filtered_data <- AD_NW_cent_df_sub_combined %>%
  dplyr::filter(ListType == "gene_list" & BroadCategory %in% c("Target", "negativeControls")) %>%
  dplyr::mutate(Type = factor(Type, levels = c("bulk", "adjusted"))) %>%
  # Ensure BroadCategory has "negativeControls" last
  dplyr::mutate(BroadCategory = factor(BroadCategory, levels = c("Target", "negativeControls"))) %>%
  # Order SubCategory by page_rank_median within each Type
  dplyr::group_by(Type, BroadCategory) %>%
  dplyr::arrange(Type, BroadCategory, desc(page_rank_median)) %>%
  dplyr::mutate(SubCategory = factor(SubCategory, levels = rev(unique(SubCategory)))) %>%
  dplyr::ungroup()

filtered_data$SubCategory <- as.character(filtered_data$SubCategory)
filtered_data[filtered_data$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"
filtered_data$SubCategory <- factor(filtered_data$SubCategory, levels=rev(unique(filtered_data$SubCategory)))

library(patchwork)

AD_NW_cent_df_sub_combined <- AD_NW_cent_df_sub_combined[AD_NW_cent_df_sub_combined$Type %in% "bulk",]

max_values <- AD_NW_cent_df_sub_combined[AD_NW_cent_df_sub_combined$ListType %in% "gene_list",] %>%
  dplyr::filter(BroadCategory %in% "negativeControls") %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(max_page_rank_median = max(page_rank_median, na.rm = TRUE)) 

AD_NW_cent_df_sub_combined[AD_NW_cent_df_sub_combined$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"

page_rank_median_g = AD_NW_cent_df_sub_combined %>%
  dplyr::filter(BroadCategory %in% c("Target", "negativeControls")) %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::select(c(ListName,page_rank_median,Type,SubCategory,CleanName,ListType)) %>%
  dplyr::mutate(SubCategory = factor(SubCategory,levels = rev(as.character(unique(filtered_data$SubCategory))))) %>%
  unique() %>%
  ggplot(aes(x = SubCategory,y = page_rank_median)) +
  facet_wrap(~ Type ,scales = "free",ncol = 1)+
  geom_point(data = function(x) subset(x,ListType != "gene_list"),aes(color = SubCategory), alpha=0.5 ,color = "grey") +
  geom_point(data = function(x) subset(x,ListType == "gene_list"),aes(color = SubCategory), size=3) +
  geom_hline(data = max_values, aes(yintercept = max_page_rank_median), linetype = "dashed", color = "black") +  # Add hline here
  coord_flip()+
  theme_minimal() + 
  theme(legend.position = "none")+
  #  theme(aspect.ratio = 1.5) + 
  border() +ylab("Median page rank") + xlab("") +
  #  scale_x_log10() +
  scale_color_manual(values = targetColors)

##### page rank tail probability#####

#nwCentrality <- read_asset("wf-0aa8f9bef8")
#nwCentrality <- read_asset("wf-f0233b864b")
nwCentrality <- read_asset("wf-69bcfe33cf") #bottom BTLA & SLAMF6

nwCentrality$Type <- factor(nwCentrality$Type, levels=c("bulk", "adjusted"))

#targetColors = read_asset("wf-8866a3b9f3") #new signature set

targetColors = read_asset("wf-c588727d1c") #bottom BTLA & SLAMF6

#targetColors['PD1'] <- "#FFE4C4"
names(targetColors)[names(targetColors) %in% "IL13.Dupi"] <- "IL13.minFDR"
names(targetColors)[names(targetColors) %in% "IL4.Dupi"] <- "IL4.minFDR"

nwCentrality <- nwCentrality[!grepl("3sd", nwCentrality$Identifier),]
nwCentrality <- nwCentrality[nwCentrality$BroadCategory %in% c("Target", "negativeControls"),]
nwCentrality <- nwCentrality[!grepl("avg", nwCentrality$Identifier),]

nwCentrality$log10pval <- -log10(nwCentrality$pval)

nwCentrality_page_rank <- nwCentrality[nwCentrality$Criteria %in% "page_rank_median",]
nwCentrality_page_rank[nwCentrality_page_rank$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"

#ZeroCenteredProbability = -(2*(0.05 - .5))
#NegLogPVAL_thr = -log10(1-abs(ZeroCenteredProbability))

nwCentrality_page_rank$SubCategory <- factor(nwCentrality_page_rank$SubCategory, levels=rev(unique(filtered_data$SubCategory)))

tail_probability__page_rank_g <- ggplot(nwCentrality_page_rank, aes(x=SubCategory, y=-log10(pval), fill=SubCategory)) + 
  facet_wrap(~Type, ncol=1)+
  geom_bar(position=position_dodge(), stat="identity",color = "black") + 
  scale_fill_manual(values = targetColors)+
  coord_flip()+
  border() +
  ylab("Permutation Tail \n Probability") + xlab("")+
  geom_hline(yintercept=1.3, linetype="dashed")+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"), 
                   legend.position = "none")


page_rank_median_g+tail_probability__page_rank_g
ggsave(paste0("~/code/p03-POC1-analysis/notes/data/page_rank_median_NW.tiff"),height = 8,width = 9)

##### page rank scaled score bar #####

#nwCentrality <- read_asset("wf-f0233b864b")
nwCentrality <- read_asset("wf-69bcfe33cf") #bottom BTLA & SLAMF6

nwCentrality$Type <- factor(nwCentrality$Type, levels=c("bulk", "adjusted"))
nwCentrality <- nwCentrality[nwCentrality$Type %in% "bulk",]

#targetColors = read_asset("wf-8866a3b9f3") #new signature set
targetColors = read_asset("wf-c588727d1c") #bottom BTLA & SLAMF6

#targetColors['PD1'] <- "#FFE4C4"
names(targetColors)[names(targetColors) %in% "IL13.Dupi"] <- "IL13.minFDR"
names(targetColors)[names(targetColors) %in% "IL4.Dupi"] <- "IL4.minFDR"

nwCentrality <- nwCentrality[!grepl("3sd", nwCentrality$Identifier),]
nwCentrality <- nwCentrality[nwCentrality$BroadCategory %in% c("Target", "negativeControls"),]
nwCentrality <- nwCentrality[!grepl("avg", nwCentrality$Identifier),]

nwCentrality$log10pval <- -log10(nwCentrality$pval)

nwCentrality_page_rank <- nwCentrality[nwCentrality$Criteria %in% "page_rank_median",]
nwCentrality_page_rank[nwCentrality_page_rank$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"

nwCentrality_page_rank$SubCategory <- factor(nwCentrality_page_rank$SubCategory, levels=rev(unique(filtered_data$SubCategory)))
nwCentrality_page_rank$Type <- factor(nwCentrality_page_rank$Type, levels=c("bulk", "adjusted"))


max_values <- nwCentrality_page_rank[nwCentrality_page_rank$BroadCategory %in% "negativeControls",] %>%
 # dplyr::filter(BroadCategory %in% "negativeControls") %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(max_SigScaled = max(SigScaled, na.rm = TRUE)) 

score__page_rank_g <- ggplot(nwCentrality_page_rank, aes(x=SubCategory, y=SigScaled, fill=SubCategory)) + 
  facet_wrap(~Type, ncol=1)+
  geom_bar(position=position_dodge(), stat="identity",color = "black") + 
  scale_fill_manual(values = targetColors)+
  coord_flip()+
  border() +
  ylab("Scaled score") + xlab("")+
  geom_hline(data = max_values, aes(yintercept = max_SigScaled), linetype = "dashed", color = "black") +  # Add hline here
    theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"), 
                   legend.position = "none")

page_rank_median_g+score__page_rank_g
ggsave(paste0("~/code/p03-POC1-analysis/notes/data/page_rank_median_scale_score_NW.tiff"),height = 8,width = 9)


##### overlap gene_list vs random#####

AD_NW_cent_df_sub_combined_Target <- AD_NW_cent_df_sub_combined[!AD_NW_cent_df_sub_combined$Identifier %in% c("steadyArm__IL13.avgFDR","steadyArm__IL4.avgFDR"),]

AD_NW_cent_df_sub_combined <- AD_NW_cent_df_sub_combined[!grepl("3sd", AD_NW_cent_df_sub_combined$ListName),]
AD_NW_cent_df_sub_combined$Type <- factor(AD_NW_cent_df_sub_combined$Type)
AD_NW_cent_df_sub_combined <- AD_NW_cent_df_sub_combined[!grepl("avg", AD_NW_cent_df_sub_combined$SubCategory),]
AD_NW_cent_df_sub_combined <- AD_NW_cent_df_sub_combined[AD_NW_cent_df_sub_combined$Type %in% "bulk",]

#fix the factor (order) of the subcategories
filtered_data <- AD_NW_cent_df_sub_combined %>%
  dplyr::filter(ListType == "gene_list" & BroadCategory %in% c("Target", "negativeControls")) %>%
  dplyr::mutate(Type = factor(Type, levels = c("bulk", "adjusted"))) %>%
  # Ensure BroadCategory has "negativeControls" last
  dplyr::mutate(BroadCategory = factor(BroadCategory, levels = c("Target", "negativeControls"))) %>%
  # Order SubCategory by page_rank_median within each Type
  dplyr::group_by(Type, BroadCategory) %>%
  dplyr::arrange(Type, BroadCategory, desc(OverlapPercent)) %>%
  dplyr::mutate(SubCategory = factor(SubCategory, levels = rev(unique(SubCategory)))) %>%
  dplyr::ungroup()
filtered_data[filtered_data$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"

library(patchwork)

max_values <- AD_NW_cent_df_sub_combined[AD_NW_cent_df_sub_combined$ListType %in% "gene_list",] %>%
  dplyr::filter(BroadCategory %in% "negativeControls") %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(max_OverlapPercent = max(OverlapPercent, na.rm = TRUE)) 

OverlapPercent_g = AD_NW_cent_df_sub_combined %>%
  dplyr::filter(BroadCategory %in% c("Target", "negativeControls")) %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::select(c(ListName,OverlapPercent,Type,SubCategory,CleanName,ListType)) %>%
  dplyr::mutate(SubCategory = factor(SubCategory,levels = rev(as.character(unique(filtered_data$SubCategory))))) %>%
  unique() %>%
  ggplot(aes(x = SubCategory,y = OverlapPercent)) +
  facet_wrap(~ Type ,scales = "free",ncol = 1)+
  geom_point(data = function(x) subset(x,ListType != "gene_list"),aes(color = SubCategory), alpha=0.5 ,color = "grey") +
  geom_point(data = function(x) subset(x,ListType == "gene_list"),aes(color = SubCategory), size=3) +
  geom_hline(data = max_values, aes(yintercept = max_OverlapPercent), linetype = "dashed", color = "black") +  # Add hline here
  coord_flip()+
  theme_minimal() + 
  theme(legend.position = "none")+
  #  theme(aspect.ratio = 1.5) + 
  border() +ylab("OverlapPercent") + xlab("") +
  #  scale_x_log10() +
  scale_color_manual(values = targetColors)

##### OverlapPercent tail probability#####

#nwCentrality <- read_asset("wf-0aa8f9bef8")
#nwCentrality <- read_asset("wf-f0233b864b")
nwCentrality <- read_asset("wf-69bcfe33cf") #bottom BTLA & SLAMF6

nwCentrality$Type <- factor(nwCentrality$Type, levels=c("bulk", "adjusted"))
nwCentrality <- nwCentrality[nwCentrality$Type %in% "bulk",]
#targetColors = read_asset("wf-8866a3b9f3") #new signature set

targetColors = read_asset("wf-c588727d1c") #bottom BTLA & SLAMF6

#targetColors['PD1'] <- "#FFE4C4"
names(targetColors)[names(targetColors) %in% "IL13.Dupi"] <- "IL13.minFDR"
names(targetColors)[names(targetColors) %in% "IL4.Dupi"] <- "IL4.minFDR"

nwCentrality <- nwCentrality[!grepl("3sd", nwCentrality$Identifier),]
nwCentrality <- nwCentrality[nwCentrality$BroadCategory %in% c("Target", "negativeControls"),]
nwCentrality <- nwCentrality[!grepl("avg", nwCentrality$Identifier),]

nwCentrality$log10pval <- -log10(nwCentrality$pval)

nwCentrality_OverlapPercent <- nwCentrality[nwCentrality$Criteria %in% "OverlapPercent",]

#ZeroCenteredProbability = -(2*(0.05 - .5))
#NegLogPVAL_thr = -log10(1-abs(ZeroCenteredProbability))

nwCentrality_OverlapPercent[nwCentrality_OverlapPercent$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"
nwCentrality_OverlapPercent$SubCategory <- factor(nwCentrality_OverlapPercent$SubCategory, levels=rev(unique(filtered_data$SubCategory)))

tail_probability__OverlapPercent_g <- ggplot(nwCentrality_OverlapPercent, aes(x=SubCategory, y=-log10(pval), fill=SubCategory)) + 
  facet_wrap(~Type, ncol=1)+
  geom_bar(position=position_dodge(), stat="identity",color = "black") + 
  scale_fill_manual(values = targetColors)+
  coord_flip()+
  border() +
  ylab("Permutation Tail \n Probability") + xlab("")+
  geom_hline(yintercept=1.3, linetype="dashed")+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"), 
                   legend.position = "none")


OverlapPercent_g+tail_probability__OverlapPercent_g
ggsave(paste0("~/code/p03-POC1-analysis/notes/data/OverlapPercent_NW.tiff"),height = 8,width = 9)

##### OverlapPercent scaled score bar #####

#nwCentrality <- read_asset("wf-f0233b864b")
nwCentrality <- read_asset("wf-69bcfe33cf") #bottom BTLA & SLAMF6

nwCentrality$Type <- factor(nwCentrality$Type, levels=c("bulk", "adjusted"))
nwCentrality <- nwCentrality[nwCentrality$Type %in% "bulk",]

#targetColors = read_asset("wf-8866a3b9f3") #new signature set
targetColors = read_asset("wf-c588727d1c") #bottom BTLA & SLAMF6

#targetColors['PD1'] <- "#FFE4C4"
names(targetColors)[names(targetColors) %in% "IL13.Dupi"] <- "IL13.minFDR"
names(targetColors)[names(targetColors) %in% "IL4.Dupi"] <- "IL4.minFDR"

nwCentrality <- nwCentrality[!grepl("3sd", nwCentrality$Identifier),]
nwCentrality <- nwCentrality[nwCentrality$BroadCategory %in% c("Target", "negativeControls"),]
nwCentrality <- nwCentrality[!grepl("avg", nwCentrality$Identifier),]

nwCentrality$log10pval <- -log10(nwCentrality$pval)

nwCentrality_OverlapPercent <- nwCentrality[nwCentrality$Criteria %in% "OverlapPercent",]

nwCentrality_OverlapPercent[nwCentrality_OverlapPercent$SubCategory %in% "PD1", "SubCategory"] <- "PDL1"
nwCentrality_OverlapPercent$SubCategory <- factor(nwCentrality_OverlapPercent$SubCategory, levels=rev(unique(filtered_data$SubCategory)))
nwCentrality_OverlapPercent$Type <- factor(nwCentrality_OverlapPercent$Type, levels=c("bulk", "adjusted"))


max_values <- nwCentrality_OverlapPercent[nwCentrality_OverlapPercent$BroadCategory %in% "negativeControls",] %>%
  # dplyr::filter(BroadCategory %in% "negativeControls") %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(max_SigScaled = max(SigScaled, na.rm = TRUE)) 

score__OverlapPercent_g <- ggplot(nwCentrality_OverlapPercent, aes(x=SubCategory, y=SigScaled, fill=SubCategory)) + 
  facet_wrap(~Type, ncol=1)+
  geom_bar(position=position_dodge(), stat="identity",color = "black") + 
  scale_fill_manual(values = targetColors)+
  coord_flip()+
  border() +
  ylab("Scaled score") + xlab("")+
  geom_hline(data = max_values, aes(yintercept = max_SigScaled), linetype = "dashed", color = "black") +  # Add hline here
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"), 
                   legend.position = "none")

OverlapPercent_g+score__OverlapPercent_g
ggsave(paste0("~/code/p03-POC1-analysis/notes/data/OverlapPercent_scale_score_NW.tiff"),height = 8,width = 9)

#### Topology ####

##### prepare data #####

SignatureTable <- read_asset("wf-b1136bff5e")

#BulkAD_NW_topology_df_sub <- read_asset("wf-07ac559790")
BulkAD_NW_topology_df_sub <- read_asset("wf-52a3f689de") #bottom BTLA & SLAMF6


#adj_AD_NW_topology_df_sub <- read_asset("wf-2d7b0b023e")
adj_AD_NW_topology_df_sub <- read_asset("wf-34716629eb") #bottom BTLA & SLAMF6

BulkAD_NW_topology_df_sub$Type <- "bulk"
adj_AD_NW_topology_df_sub$Type <- "adjusted"

AD_NW_topology_df_sub <- rbind(BulkAD_NW_topology_df_sub, adj_AD_NW_topology_df_sub)

##### aggregate random #####
#aggregate random and randomsmooth before scaling

SignatureTable <- read_asset("wf-b1136bff5e")

AD_NW_topology_df_sub <-  AD_NW_topology_df_sub %>%
  dplyr::mutate(Identifier = ListName) %>%
  dplyr::left_join(SignatureTable)

AD_NW_topology_df_sub$largest_connected_component_entrez <- NULL
AD_NW_topology_df_sub$all_components_frac <- NULL
AD_NW_topology_df_sub$largest_connected_component_symbols <- NULL


mean_topology_df_random <- AD_NW_topology_df_sub %>%
  dplyr::mutate(Identifier = ListName) %>%
  dplyr::left_join(SignatureTable) %>%
  dplyr::filter(SubCategory %in% c("random", "smoothedRandom")) %>%
  dplyr::group_by(ListType, BroadCategory, SubCategory, Type, geneset) %>%
  dplyr::mutate(across(number_of_genes:largest_connected_component_diameter_weighted_frac, 
                       ~ as.numeric(.))) %>%
  dplyr::mutate(across(number_of_genes:largest_connected_component_diameter_weighted_frac, 
                       ~ ifelse(is.nan(.) | !is.finite(.), NA, .))) %>%
  dplyr::summarise(across(number_of_genes:largest_connected_component_diameter_weighted_frac, 
                          mean, na.rm = TRUE))

mean_topology_df_random$ListName <- mean_topology_df_random$SubCategory
mean_topology_df_random$Identifier <- mean_topology_df_random$ListName
mean_topology_df_random$CleanName <- mean_topology_df_random$ListName


AD_NW_topology_df_sub$ModID <- NULL
AD_NW_topology_df_sub$Identifier_original <- NULL

##### combine dataframes #####

mean_topology_df_random <- mean_topology_df_random[,colnames(mean_topology_df_random) %in% colnames(AD_NW_topology_df_sub)]

AD_NW_topology_df_sub_combined <- rbind.fill(AD_NW_topology_df_sub[-grep("random|smoothedRandom", AD_NW_topology_df_sub$Identifier),] ,mean_topology_df_random )
cyto_cc_push_object(AD_NW_topology_df_sub_combined,"AD_NW_topology_df_sub_combined")
#wf-6a8264ced5
#wf-06a1afa239 #bottom BTLA & SLAMF6


##### density gene_list vs random#####
AD_NW_topology_df_sub_combined <- read_asset("wf-06a1afa239")
AD_NW_topology_df_sub_combined <- AD_NW_topology_df_sub_combined[!grepl("3sd", AD_NW_topology_df_sub_combined$ListName),]
AD_NW_topology_df_sub_combined$Type <- factor(AD_NW_topology_df_sub_combined$Type)
AD_NW_topology_df_sub_combined <- AD_NW_topology_df_sub_combined[!grepl("avg", AD_NW_topology_df_sub_combined$SubCategory),]
AD_NW_topology_df_sub_combined <- AD_NW_topology_df_sub_combined[AD_NW_topology_df_sub_combined$Type %in% "bulk",]

AD_NW_topology_df_sub_combined[AD_NW_topology_df_sub_combined$SubCategory %in% "PD1", 'SubCategory'] <- "PDL1"

#fix the factor (order) of the subcategories
filtered_data <- AD_NW_topology_df_sub_combined %>%
  dplyr::filter(ListType == "gene_list" & BroadCategory %in% c("Target", "negativeControls")) %>%
  dplyr::mutate(Type = factor(Type, levels = c("bulk", "adjusted"))) %>%
  dplyr::mutate(BroadCategory = factor(BroadCategory, levels = c("Target", "negativeControls"))) %>%
  dplyr::group_by(Type, BroadCategory) %>%
  dplyr::arrange(Type, BroadCategory, desc(density)) %>%
  dplyr::mutate(SubCategory = factor(SubCategory, levels = rev(unique(SubCategory)))) %>%
  dplyr::ungroup()


max_values <- AD_NW_topology_df_sub_combined[AD_NW_topology_df_sub_combined$ListType %in% "gene_list",] %>%
  dplyr::filter(BroadCategory %in% "negativeControls") %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(density = max(density, na.rm = TRUE)) 

density_g = AD_NW_topology_df_sub_combined %>%
  dplyr::filter(BroadCategory %in% c("Target", "negativeControls")) %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::select(c(ListName,density,Type,SubCategory,CleanName,ListType)) %>%
  dplyr::mutate(SubCategory = factor(SubCategory,levels = rev(as.character(unique(filtered_data$SubCategory))))) %>%
  unique() %>%
  ggplot(aes(x = SubCategory,y = density)) +
  facet_wrap(~ Type ,scales = "free",ncol = 1)+
  geom_point(data = function(x) subset(x,ListType != "gene_list"),aes(color = SubCategory), alpha=0.5 ,color = "grey") +
  geom_point(data = function(x) subset(x,ListType == "gene_list"),aes(color = SubCategory), size=3) +
  geom_hline(data = max_values, aes(yintercept = density), linetype = "dashed", color = "black") +  # Add hline here
  coord_flip()+
  theme_minimal() + 
  theme(legend.position = "none")+
  #  theme(aspect.ratio = 1.5) + 
  border() +ylab("Subgraph density") + xlab("") +
  #  scale_x_log10() +
  scale_color_manual(values = targetColors)


##### density tail probability #####

#nwTopology <- read_asset("wf-e5306fac89")
nwTopology <- read_asset("wf-68d58979bd") #bottom BTLA &SLAMF6

#targetColors = read_asset("wf-8866a3b9f3") #new signature set
targetColors = read_asset("wf-c588727d1c") 

#targetColors['PD1'] <- "#FFE4C4"
names(targetColors)[names(targetColors) %in% "IL13.Dupi"] <- "IL13.minFDR"
names(targetColors)[names(targetColors) %in% "IL4.Dupi"] <- "IL4.minFDR"

nwTopology <- nwTopology[!grepl("3sd", nwTopology$Identifier),]
nwTopology <- nwTopology[nwTopology$BroadCategory %in% c("Target", "negativeControls"),]
nwTopology <- nwTopology[!grepl("avg", nwTopology$Identifier),]
nwTopology <- nwTopology[nwTopology$Type %in% "bulk",]
nwTopology[nwTopology$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"

nwTopology$log10pval <- -log10(nwTopology$pval)

nwTopology_density <- nwTopology[nwTopology$Criteria %in% "density",]

nwTopology_density$SubCategory <- factor(nwTopology_density$SubCategory, levels=rev(unique(filtered_data$SubCategory)))

tail_probability__density_g <- ggplot(nwTopology_density, aes(x=SubCategory, y=-log10(pval), fill=SubCategory)) + 
  facet_wrap(~Type, ncol=1)+
  geom_bar(position=position_dodge(), stat="identity",color = "black") + 
  scale_fill_manual(values = targetColors)+
  coord_flip()+
  border() +
  ylab("Permutation Tail \n Probability") + xlab("")+
  geom_hline(yintercept=1.3, linetype="dashed")+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"), 
                   legend.position = "none")


density_g+tail_probability__density_g
ggsave(paste0("~/code/p03-POC1-analysis/notes/data/density_NW.tiff"),height = 8,width = 9)


##### density scaled score #####

#nwTopology <- read_asset("wf-e5306fac89")
nwTopology <- read_asset("wf-68d58979bd") #bottom BTLA & SLAMF6

#targetColors = read_asset("wf-8866a3b9f3") #new signature set
targetColors = read_asset("wf-c588727d1c") 


#targetColors['PD1'] <- "#FFE4C4"
names(targetColors)[names(targetColors) %in% "IL13.Dupi"] <- "IL13.minFDR"
names(targetColors)[names(targetColors) %in% "IL4.Dupi"] <- "IL4.minFDR"

nwTopology <- nwTopology[!grepl("3sd", nwTopology$Identifier),]
nwTopology <- nwTopology[nwTopology$BroadCategory %in% c("Target", "negativeControls"),]
nwTopology <- nwTopology[!grepl("avg", nwTopology$Identifier),]

nwTopology$log10pval <- -log10(nwTopology$pval)
nwTopology <- nwTopology[nwTopology$Type %in% "bulk",]
nwTopology[nwTopology$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"

nwTopology_density <- nwTopology[nwTopology$Criteria %in% "density",]

nwTopology_density$SubCategory <- factor(nwTopology_density$SubCategory, levels=rev(unique(filtered_data$SubCategory)))
nwTopology_density$Type <- factor(nwTopology_density$Type, levels=c("bulk", "adjusted"))

max_values <- nwTopology_density[nwTopology_density$BroadCategory %in% "negativeControls",] %>%
  # dplyr::filter(BroadCategory %in% "negativeControls") %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(max_SigScaled = max(SigScaled, na.rm = TRUE)) 


Sigscale_density_g <- ggplot(nwTopology_density, aes(x=SubCategory, y=SigScaled, fill=SubCategory)) + 
  facet_wrap(~Type, ncol=1)+
  geom_bar(position=position_dodge(), stat="identity",color = "black") + 
  scale_fill_manual(values = targetColors)+
  coord_flip()+
  border() +
  ylab("Scaled score") + xlab("")+
  geom_hline(data = max_values, aes(yintercept = max_SigScaled), linetype = "dashed", color = "black") +  # Add hline here
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"), 
                   legend.position = "none")

density_g+Sigscale_density_g
ggsave(paste0("~/code/p03-POC1-analysis/notes/data/density_scaled_score_NW.tiff"),height = 8,width = 9)



#### test co-correlation between selected centrality measures ####

#centrality
AD_NW_cent_df_sub_combined_gene_list <- AD_NW_cent_df_sub_combined[AD_NW_cent_df_sub_combined$ListType %in% "gene_list",] 
selected_param <- c("eigen_centrality_median", "degree_median", "page_rank_median", "closeness_median", "betweenness_median", "OverlapPercent", "ListName", "ListType", "Type", "Identifier", "CleanName", "BroadCategory", "SubCategory")
AD_NW_cent_df_sub_combined_gene_list_sub <- AD_NW_cent_df_sub_combined_gene_list[,colnames(AD_NW_cent_df_sub_combined_gene_list) %in% selected_param,]

AD_NW_cent_df_sub_combined_gene_list_sub <- reshape2::melt(AD_NW_cent_df_sub_combined_gene_list_sub, id.vars=c("ListName", "ListType", "Type", "Identifier", "CleanName", "BroadCategory", "SubCategory"))


#topology
AD_NW_topology_df_sub_combined_gene_list <- AD_NW_topology_df_sub_combined[AD_NW_topology_df_sub_combined$geneset %in% "gene_list",]
selected_param <- c("density", "ListName", "ListType", "Type", "Identifier", "CleanName", "BroadCategory", "SubCategory")
AD_NW_topology_df_sub_combined_gene_list_sub <- AD_NW_topology_df_sub_combined_gene_list[,colnames(AD_NW_topology_df_sub_combined_gene_list) %in% selected_param,]

AD_NW_topology_df_sub_combined_gene_list_sub <- reshape2::melt(AD_NW_topology_df_sub_combined_gene_list_sub, id.vars=c("ListName", "ListType", "Type", "Identifier", "CleanName", "BroadCategory", "SubCategory"))


AD_NW_combined_gene_list_sub <- rbind.fill(AD_NW_cent_df_sub_combined_gene_list_sub,AD_NW_topology_df_sub_combined_gene_list_sub)
AD_NW_combined_gene_list_sub <- AD_NW_combined_gene_list_sub[AD_NW_combined_gene_list_sub$BroadCategory %in% c("Target", "negativeControls"),]

AD_NW_combined_gene_list_sub_wide <- reshape2::dcast(AD_NW_combined_gene_list_sub, Identifier+Type ~variable, value.var = "value")

AD_NW_combined_gene_list_sub_wide$Identifier_type <- paste(AD_NW_combined_gene_list_sub_wide$Identifier, AD_NW_combined_gene_list_sub_wide$Type, sep="__")

row.names(AD_NW_combined_gene_list_sub_wide) <- AD_NW_combined_gene_list_sub_wide$Identifier_type
AD_NW_combined_gene_list_sub_wide$Identifier <- NULL
AD_NW_combined_gene_list_sub_wide$Type <- NULL
AD_NW_combined_gene_list_sub_wide$Identifier_type <- NULL


cor_criteria <- psych::corr.test(AD_NW_combined_gene_list_sub_wide, method="pearson")
cor_r <- cor_criteria$r


colors <- c("#f0e3ab", "#f49d77", "#ce5f5f", "#b23030")  # Example colors
breaks <- seq(0, 1, length.out = length(colors) )

ht_criteria_cor <- Heatmap(cor_r,
                   #split = Target_annot,
                   #rect_gp = gpar(col = "white", lwd = 2), row_title = NULL,
                  # col = colorRamp2(breaks, colors),
                   name="nw criteria measure value", 
                   row_names_gp = gpar(fontsize = 8), 
                   column_names_gp = gpar(fontsize = 7), show_row_names = TRUE, show_column_names = T,
                   row_title_gp = gpar(fontsize = 7),
                   row_title_rot=0,
                   column_names_rot = 90, column_names_max_height = unit(10,"cm"), 
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.2f", cor_r[i,j]), x, y, gp = gpar(fontsize = 7))})

draw(ht_criteria_cor, padding = unit(c(2, 2, 2, 90), "mm"), heatmap_legend_side = "left") #bottom, left, top, right paddings



