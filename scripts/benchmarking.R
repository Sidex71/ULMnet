library(decoupleR)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(caret)
library(org.Mm.eg.db)
library(CAMML)
library(Neighborseq)
#install.packages("gtools")
library(gtools)
library(data.table)
library(xgboost)
library(parallel)
library(multiROC)
library(RColorBrewer)
library(DoubletFinder)
#BiocManager::install("scDblFinder")
library(scDblFinder)
#################################ULM################################

##########invitro DC-T data#####################
Idents(invitro_obj) <- invitro_obj$sorting_scheme
###generate signature list from singlets
invit_mark <- subset(invitro_obj, idents = c("TCRb+","CD11c+"))
Markers <- FindAllMarkers(invit_mark, only.pos = T)
celltypes <- unique(Markers$cluster)
DEG_filt <-  as.list(rep(NA, length(celltypes)))

for (i in 1: length(celltypes)){
  filt <- Markers %>% filter(cluster==celltypes[i]) %>% filter(p_val_adj <0.05) %>%
    slice_max(., order_by = avg_log2FC, n = 100) %>% pull(gene)
  DEG_filt[[i]] <- filt
  names(DEG_filt)[i] <- as.character(celltypes)[i]
}
#saveRDS(DEG_filt, 'DEG_filt_invivo.rds')

mat <- as.matrix(invitro_obj@assays$RNA3@data)
sig <- data.frame(DEG_filt)
sig <- pivot_longer(sig, names_to = 'source', values_to = 'target', cols = everything())
sig$mor <- 1
acts <- run_ulm(mat=mat, net=sig, .source='source', .target='target',
                .mor='mor', minsize = 5)
acts_filt <- acts %>% filter(p_value <= 0.05 & score > 1) %>% mutate(count_ulm = 1)
saveRDS(acts, 'acts_invitro.rds')
cell_class <- acts_filt %>% 
  group_by(condition) %>% 
  mutate(count_ulm= sum(count_ulm),
         celltype_ulm= paste(source, collapse = '_'),
         avg_pvalue = mean(p_value),
         avg_score = mean(score))%>% 
  select(-c(source, p_value, score)) %>% distinct()
invit_meta <- invitro_obj@meta.data %>% left_join(cell_class, by = c('barcode' = 'condition'))
table(invit_meta$sorting_scheme, invit_meta$celltype_ulm)
invitro_obj@meta.data <- invit_meta
rownames(invitro_obj@meta.data) <- invit_meta$barcode
saveRDS(cell_class, 'cell_class_invitro.rds')
saveRDS(invitro_obj, 'invitro_obj.rds')
#invitro_obj@meta.data<- invitro_obj@meta.data %>% dplyr:: select(-c("count_ulm" , "celltype_ulm", "avg_pvalue", "avg_score", 'statistic'))

#############################confusion matrix

actual <- ifelse(invitro_obj$sorting_scheme=='CD11c+' | 
                   invitro_obj$sorting_scheme=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invitro_obj$celltype_ulm=='CD11c._TCRb.', 'Yes', 'No')

my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_invitro_ulm.rds')

#######################invivo DC-T################
#####get signature list from singlets
invivo_obj$cell_class <- str_replace(invivo_obj$sorting_scheme, 'Ag\\+ ', '') %>%
  str_replace('CD160\\+ ', '') %>% str_replace('Icos\\+ ', '') %>%
  str_replace('Tigit\\+ ', '')
Idents(invivo_obj) <- invivo_obj$cell_class
invit_mark <- subset(invivo_obj, idents = c("TCRb+","CD11c+"))
Markers <- FindAllMarkers(invit_mark, only.pos = T)
celltypes <- unique(Markers$cluster)
DEG_filt <-  as.list(rep(NA, length(celltypes)))

for (i in 1: length(celltypes)){
  filt <- Markers %>% filter(cluster==celltypes[i]) %>% filter(p_val_adj <0.05) %>%
    slice_max(., order_by = avg_log2FC, n = 100) %>% pull(gene)
  DEG_filt[[i]] <- filt
  names(DEG_filt)[i] <- as.character(celltypes)[i]
}
#saveRDS(DEG_filt, 'DEG_filt_invivo.rds')

mat <- as.matrix(invivo_obj@assays$RNA3@data)
sig <- data.frame(DEG_filt_invivo)
sig <- pivot_longer(sig, names_to = 'source', values_to = 'target', cols = everything())
sig$mor <- 1
acts <- run_ulm(mat=mat, net=sig, .source='source', .target='target',
                .mor='mor', minsize = 5)
acts_filt <- acts %>% filter(p_value <= 0.05 & score > 1) %>% mutate(count_ulm = 1)
saveRDS(acts, 'acts_invivo.rds')
cell_class <- acts_filt %>% 
  group_by(condition) %>% 
  mutate(count_ulm= sum(count_ulm),
         celltype_ulm= paste(source, collapse = '_'),
         avg_pvalue = mean(p_value),
         avg_score = mean(score))%>% 
  select(-c(source, p_value, score)) %>% distinct()
invit_meta <- invivo_obj@meta.data %>% left_join(cell_class, by = c('barcode' = 'condition'))
table(invit_meta$sorting_scheme, invit_meta$celltype_ulm)
invivo_obj@meta.data <- invit_meta
rownames(invivo_obj@meta.data) <- invit_meta$barcode
saveRDS(cell_class, 'cell_class_invivo.rds')
saveRDS(invivo_obj, 'invivo_obj.rds')

###################################confusion matrix
library(caret)
actual<- ifelse(invivo_obj$sorting_scheme=='Ag+ TCRb+ CD11c+' |
                  invivo_obj$sorting_scheme=='TCRb+ CD11c+', 'Yes', 'No')
predicted <- ifelse(invivo_obj$celltype_ulm=='CD11c._TCRb.', 'Yes', 'No')

my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_invivo_ulm.rds')

################################################Hepatocyte-endothelial data
library(tidyverse)
sig <- read_tsv('/mnt/8TB/users/shameed/shameed/Doublet predictions/liver/PanglaoDB_markers_27_Mar_2020.tsv')
sig <- sig %>% 
  filter(species != 'Hs' & (`cell type`== 'Hepatocytes' | `cell type`=="Endothelial cells")) %>%
  dplyr::select(c(`cell type`, `official gene symbol`)) %>% 
  mutate(cap = substr(`official gene symbol`, 1, 1),
         low=tolower(substr(`official gene symbol`, 2, nchar(`official gene symbol`)))) %>%
  mutate(gene = paste0(cap, low),
         mor = 1) %>% 
  dplyr:: select(-c(cap, low,`official gene symbol`))
#table(sig$`cell type`)
mat <- GetAssayData(PairedData, assay = 'RNA3', slot = 'data')

set.seed(101124)
acts <- run_ulm(mat=mat, net= sig, .source='cell type', .target='gene',
                .mor='mor', minsize = 5)
acts_filt <- acts %>% filter(p_value <= 0.05 & score > 1) %>% mutate(count_ulm = 1)

saveRDS(acts, 'acts_PairedData.rds')

cell_class <- acts_filt %>% 
  group_by(condition) %>% 
  mutate(count_ulm= sum(count_ulm),
         celltype_ulm= paste(source, collapse = '_'),
         avg_pvalue = mean(p_value),
         avg_score = mean(score))%>% 
  dplyr::select(-c(source, p_value, score)) %>% distinct()

invit_meta <- PairedData@meta.data %>% rownames_to_column('barcode_ulm') %>% 
  left_join(cell_class, by = c('barcode_ulm' = 'condition')) %>% column_to_rownames('barcode_ulm')
table(invit_meta$CellType, invit_meta$celltype_ulm) 
PairedData@meta.data <- invit_meta
saveRDS(cell_class, 'cell_class_PairedData.rds')
saveRDS(PairedData, 'PairedData.rds')

######################################confusion matrix
actual <- ifelse(MergedData$CellType=='Hep_Endo', 'Yes', 'No')
predicted <- MergedData$man_pred
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_marker.rds')
man_stat <- matcon$byClass[1:4]

#################################NEIGHBORSEQ############################################
#########################################invitro DC-T data
Idents(invitro_obj)<- invitro_obj$sorting_scheme
sing_obj <- subset(invitro_obj, idents= c('CD11c+', 'TCRb+'))
set.seed(23042024)
ns.data = prep_cell_mat(ge = sing_obj@assays$RNA3, celltypes =sing_obj$sorting_scheme, topn = 50)
#saveRDS(ns.data, 'ns.data.rds')
markers <- as.data.frame(ns.data$markers)
markers_adj <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n=50)
ns.data$markers <- markers_adj

mt = multiplet_types(ns.data$celltypes)
#saveRDS(mt, 'mt.rds')
am = artificial_multiplets(cell.mat = as.matrix(ns.data$cell.mat), 
                           celltypes = ns.data$celltypes, 
                           multiplet_classes = mt)
#saveRDS(am, 'am_invitro.rds')

rf = multiplet_rf(am)

#saveRDS(rf, 'rf_invitro.rds')
mroc = mroc_format(rf$test, rf$pred)
mroc_plot(mroc)
#saveRDS(mroc, 'mroc.rds')

#####create a new matrix from test data
my_cell.mat =invitro_obj@assays$RNA3@counts[unique(markers_adj$gene),]
pred = xgpred(rf, as.matrix(my_cell.mat))
#saveRDS(pred, 'pred_invitro.rds')
neigb_type<- colnames(pred)[apply(pred,1, which.max)]
neigb_pred <- ifelse(neigb_type=='CD11c+_CD11c+', 'CD11c+', 
                     ifelse(neigb_type=='TCRb+_TCRb+', 'TCRb+', neigb_type))
invitro_obj$neigb_pred <- neigb_pred

##########################confusion matrix and accuracy
actual <- ifelse(invitro_obj$sorting_scheme=='CD11c+' | 
                   invitro_obj$sorting_scheme=='TCRb+', 'No', 'Yes')
predicted <- ifelse(neigb_pred=='CD11c+_TCRb+', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
#saveRDS(matcon, 'confmat_invitro_neigb.rds')
#saveRDS(invitro_obj, 'invitro_obj.rds')

#######################################################invivo DC-T data
invivo_obj$cell_class <- str_replace(invivo_obj$sorting_scheme, 'Ag\\+ ', '') %>%
  str_replace('CD160\\+ ', '') %>% str_replace('Icos\\+ ', '') %>%
  str_replace('Tigit\\+ ', '')
Idents(invivo_obj) <- invivo_obj$cell_class
sing_obj <- subset(invivo_obj, idents= c('CD11c+', 'TCRb+'))
set.seed(230420242)
ns.data = prep_cell_mat(ge = sing_obj@assays$RNA3, celltypes =sing_obj$cell_class, topn = 50)
#saveRDS(ns.data, 'ns.data.rds')
markers <- as.data.frame(ns.data$markers)
markers_adj <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n=50)
ns.data$markers <- markers_adj

mt = multiplet_types(ns.data$celltypes)
#saveRDS(mt, 'mt.rds')
am = artificial_multiplets(cell.mat = as.matrix(ns.data$cell.mat), 
                           celltypes = ns.data$celltypes, 
                           multiplet_classes = mt)
saveRDS(am, 'am_invivo.rds')

rf = multiplet_rf(am)
#view(rf$test)
saveRDS(rf, 'rf_invivo.rds')
#> [1]  train-mlogloss:3.545750
mroc = mroc_format(rf$test, rf$pred)
mroc_plot(mroc)

###import all cells and predict
my_cell.mat = invivo_obj@assays$RNA3@counts[unique(markers_adj$gene),]
pred = xgpred(rf, as.matrix(my_cell.mat))
saveRDS(pred, 'pred_invivo.rds')
neigb_type<- colnames(pred)[apply(pred,1, which.max)]
neigb_pred <- ifelse(neigb_type=='CD11c+_CD11c+', 'CD11c+', 
                     ifelse(neigb_type=='TCRb+_TCRb+', 'TCRb+', neigb_type))
invivo_obj$neigb_pred <- neigb_pred

##########################################confusion matrix
actual <- ifelse(invivo_obj$cell_class=='CD11c+' | 
                   invivo_obj$cell_class=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invivo_obj$neigb_pred=='CD11c+_TCRb+', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(invivo_obj, 'invivo_obj.rds')
saveRDS(matcon, 'confmat_invivo_neigb.rds')

######################################################Hepatocyte-endothelial data
LiverData$lineage <- ifelse(LiverData$CellType=="Endothelial cells" | 
                              LiverData$CellType=="Hepatocytes"|
                              LiverData$CellType=="Hepatic progenitor cells", 'HE_Lineage', 'Others')
HepEndData <- subset(LiverData, subset = lineage == 'HE_Lineage')
HepEndData <- subset(HepEndData, subset = CellType ==c("Endothelial cells", "Hepatocytes")) 
#saveRDS(HepEndData, 'HepEndData.rds')
my_names<- rownames(PairedData)
my_genes <- unique(markers_adj$gene)
my_names <- my_names[my_names %in% my_genes]
set.seed(07042024)
ns.data$cell.mat <- HepEndData@assays$RNA3@counts[my_names,]
#install.packages("gtools")

mt = multiplet_types(ns.data$celltypes)
#saveRDS(mt, 'mt_pairedData.rds')
am = artificial_multiplets(cell.mat = as.matrix(ns.data$cell.mat), 
                           celltypes = ns.data$celltypes, 
                           multiplet_classes = mt)
#saveRDS(am, 'am_pairedData.rds')


rf = multiplet_rf(am)

#saveRDS(rf, 'rf_pairedData.rds')
mroc = mroc_format(rf$test, rf$pred)
mroc_plot(mroc)
#saveRDS(mroc, 'mroc.rds')

my_cell.mat = PairedData@assays$RNA3@counts[my_names,]
pred = xgpred(rf, as.matrix(my_cell.mat))
#saveRDS(pred, 'pred_pairedData.rds')
neigb_type<- colnames(pred)[apply(pred,1, which.max)]
table(neigb_type)
neigb_pred <- ifelse(neigb_type=='Endothelial cells_Endothelial cells', 'Endothelial cells', 
                     ifelse(neigb_type=='Hepatocytes_Hepatocytes', 'Hepatocytes', neigb_type))
table(neigb_pred)
PairedData$neigb <- neigb_pred
table(PairedData$CellType, PairedData$neigb)
#saveRDS(PairedData, 'PairedData.rds')


##########################################CICADA###############################################
#########################################################invitro DC-T
Idents(invitro_obj) <- invitro_obj$sorting_scheme
sing_obj <- subset(invitro_obj, idents= c('CD11c+', 'TCRb+'))
my_set <- CAMML::BuildGeneSets(as.matrix(sing_obj@assays$RNA3@data), species = 'Mm', labels = sing_obj$cell_class)
seu_score_3 <- CAMML(invitro_obj, my_set)
CAM_score <- as.data.frame(seu_score_3@assays$CAMML$data)

####function for cell classification
CAM_calling <- function(my_score, my_thresh, my_cellnames) {
  names(my_score) <- my_cellnames
  my_class <- my_score > my_thresh
  my_celltype <- ifelse(sum(my_class) == 0 | sum(my_class) == 1, 
                        names(which.max(my_score)),
                        paste(names(my_score)[my_class], collapse = '_'))
  return(my_celltype)
}
####classifying cells with the function
Cicada_out <- sapply(CAM_score, CAM_calling, 0.5, rownames(CAM_score))
table(Cicada_out)
Cicada_out <- as.data.frame(Cicada_out)
invitro_obj$Cicada <- Cicada_out
table(invitro_obj$sorting_scheme, invitro_obj$Cicada)

####################################confusion matrix
actual <- ifelse(invitro_obj$sorting_scheme=='CD11c+' | 
                   invitro_obj$sorting_scheme=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invitro_obj$Cicada=='CD11c+_TCRb+', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_Cicada_invitro.rds')


###########################################################invivo DC-T
Idents(invivo_obj) <- invivo_obj$sorting_scheme
sing_obj <- subset(invivo_obj, idents= c('CD11c+', 'TCRb+'))
library(org.Mm.eg.db)
my_set <- CAMML::BuildGeneSets(as.matrix(sing_obj@assays$RNA3@data), species = 'Mm', labels = sing_obj$cell_class)
seu_score_3 <- CAMML(invivo_obj, my_set)
CAM_score <- as.data.frame(seu_score_3@assays$CAMML$data)

####function for cell classification
CAM_calling <- function(my_score, my_thresh, my_cellnames) {
  names(my_score) <- my_cellnames
  my_class <- my_score > my_thresh
  my_celltype <- ifelse(sum(my_class) == 0 | sum(my_class) == 1, 
                        names(which.max(my_score)),
                        paste(names(my_score)[my_class], collapse = '_'))
  return(my_celltype)
}
####classifying cells with the function
Cicada_out <- sapply(CAM_score, CAM_calling, 0.5, rownames(CAM_score))
table(Cicada_out)
Cicada_out <- as.data.frame(Cicada_out)
invivo_obj$Cicada <- Cicada_out
table(invivo_obj$sorting_scheme, invivo_obj$Cicada)

####################################confusion matrix
actual <- ifelse(invivo_obj$cell_class=='CD11c+' | 
                   invivo_obj$cell_class=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invivo_obj$Cicada=='CD11c+_TCRb+', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon

saveRDS(matcon, 'confmat_Cicada_invivo.rds')

#######################################################Hepatocyte-endothelial data

liver_sub <- subset(LiverData, downsample = 1000)
my_set <- CAMML::BuildGeneSets(as.matrix(liver_sub@assays$RNA3@data), species = 'Mm', labels =liver_sub$CellType)
#saveRDS(my_set, 'my_set.rds')
hep_set <- my_set %>% filter(cell.type=='Endothelial cells'| cell.type=='Hepatocytes')
set.seed(111024)
seu_score <- CAMML(PairedData, gene.set.df= hep_set)
CAM_score <- as.data.frame(seu_score@assays$CAMML$data)

####function for cell classification
CAM_calling <- function(my_score, my_thresh, my_cellnames) {
  names(my_score) <- my_cellnames
  my_class <- my_score > my_thresh
  my_celltype <- ifelse(sum(my_class) == 0 | sum(my_class) == 1, 
                        names(which.max(my_score)),
                        paste(names(my_score)[my_class], collapse = '_'))
  return(my_celltype)
}
####classifying cells with the function
Cicada_out <- sapply(CAM_score, CAM_calling, 0.5, rownames(CAM_score))
table(Cicada_out)
PairedData$Cicada <- Cicada_out
#saveRDS(PairedData, 'PairedData.rds')

#########################################Model comparison and plotting########################################
################################################### DC-T
conf_df <- rbind(confmat_Cicada_invitro$byClass[1:4], 
                 confmat_invitro_neigb$byClass[1:4],
                 confmat_invitro_ulm$byClass[1:4])
my_meth <- c('CICADA', 'Neigborseq', 'ulm')
conf_df<- as.data.frame(conf_df)
conf_df$Method <- my_meth
conf_df$Tissue <- 'in vitro'
conf_df <- conf_df[, c('Tissue','Method', 'Sensitivity', "Specificity",  "Pos Pred Value", "Neg Pred Value")]
colnames(conf_df) <- c('Tissue','Method', 'Precision', "Neg pred",  "Sensitivity", "Specificity")
inv_df <- rbind(confmat_Cicada_invivo$byClass[1:4], 
                confmat_invivo_neigb$byClass[1:4],
                confmat_invivo_ulm$byClass[1:4])
my_meth <- c('CICADA', 'Neigborseq', 'ulm')
inv_df<- as.data.frame(inv_df)
inv_df$Method <- my_meth
inv_df$Tissue <- 'in vivo'
inv_df <- inv_df[, c('Tissue','Method', 'Sensitivity', "Specificity",  "Pos Pred Value", "Neg Pred Value")]
colnames(inv_df) <- c('Tissue','Method', 'Precision', "Neg pred",  "Sensitivity", "Specificity")

confmat_all <- rbind(conf_df, inv_df)
confmat_all$Overall_accuracy <- (confmat_all$Sensitivity + confmat_all$Specificity)/2
saveRDS(confmat_all, '/mnt/8TB/users/shameed/shameed/Doublet predictions/confmat_all.rds')

#############plot

#colnames(confmat_all)[5] <- '% doublet detected'
confmat_all$Method <- str_replace(confmat_all$Method, 'ulm', 'ULMnet')

data <- confmat_all %>% filter(Tissue=='in vitro') %>% 
  dplyr:: select(c(Method, Sensitivity, Specificity, Precision, Tissue)) %>%
  pivot_longer(cols = -c(Method, Tissue), names_to = 'metrics', values_to = 'percentage') %>%
  mutate(percentage = percentage * 100)
#data <- data %>% filter(!Method %in% c('Marker', 'AUC'))
p1<- ggplot(data, aes(x =factor(metrics, levels = c('Sensitivity', 'Precision', 'Specificity')), y = percentage, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.5) +  # Bar position dodge width should match text
  geom_text(aes(label = round(percentage, 1)), 
            position = position_dodge(width = 0.8),  # Ensure the same width in dodge for text
            vjust = -0.5, size = 4.5) +  # Adjust vjust for better alignment
  labs(title = "Performance Comparison by Method (in vitro)", y = "Percentage", x = '', fill = 'Methods') +
  theme_minimal() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, hjust = 0.5, face = 'bold'),
        legend.key.size = unit(1.5, "cm"),
        legend.title = element_text(size = 17, hjust = 0.5, face = 'bold'),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name= 'Method')

#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/stat_invitro.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()
########################
data <- confmat_all %>% filter(Tissue=='in vivo') %>% 
  dplyr:: select(c(Method, Sensitivity, Specificity, Precision, Tissue)) %>%
  pivot_longer(cols = -c(Method, Tissue), names_to = 'metrics', values_to = 'percentage') %>%
  mutate(percentage = percentage * 100)
#data <- data %>% filter(!Method %in% c('Marker', 'AUC'))
p2<- ggplot(data, aes(x =factor(metrics, levels = c('Sensitivity', 'Precision', 'Specificity')), y = percentage, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.5) +  # Bar position dodge width should match text
  geom_text(aes(label = round(percentage, 1)), 
            position = position_dodge(width = 0.8),  # Ensure the same width in dodge for text
            vjust = -0.5, size = 4.5) +  # Adjust vjust for better alignment
  labs(title = "Performance Comparison by Method (in vivo)", y = "Percentage", x = '', fill = 'Methods') +
  theme_minimal() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, hjust = 0.5, face = 'bold'),
        legend.key.size = unit(1.5, "cm"),
        legend.title = element_text(size = 17, hjust = 0.5, face = 'bold'),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name= 'Method')

#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/stat_invivo.png", width = 15, height = 10.5, units = 'in', res = 600)
p2
dev.off()


align_fig <-ggarrange(p1, p2, ncol = 2, nrow=1, common.legend = T, legend="right")
print(align_fig)

png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/stat_DC_T.png", width = 22, height = 10.5, units = 'in', res = 600)
align_fig
dev.off()
############################################Hepatocyte-endothelial
neigb_det <- table(PairedData$neigb)[1]
CICADA_det <- table(PairedData$Cicada)[2]
ulm_det <- table(PairedData$celltype_ulm)[1]
number_detected <- c(neigb_det, CICADA_det, ulm_det)
detection_rate <- number_detected *100/length(PairedData$CellType)
my_metrics <- data.frame(number_detected, detection_rate) 
rownames(my_metrics) <- c('Neiborseq', 'CICADA', 'ulm')
saveRDS(my_metrics, 'my_metrics.rds')

data <- my_metrics %>% rownames_to_column('Method')
colnames(data)[3] <- '% doublet detected'
data$Method <- str_replace(data$Method, 'ulm', 'ULMnet')

p1<-ggplot(data, aes(x = Method, y = `% doublet detected`, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.3) +  # Bar position dodge width should match text
  geom_text(aes(label = round(`% doublet detected`, 1)), 
            position = position_dodge(width = 0.5),  # Ensure the same width in dodge for text
            vjust = -0.5, size = 4.5) +  # Adjust vjust for better alignment
  labs(title = "Performance Comparison by Method", x = '', fill = 'Methods') +
  theme_minimal() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        #legend.text = element_text(size = 17, hjust = 0.5, face = 'bold'),
        #legend.key.size = unit(1.5, "cm"),
        #legend.title = element_text(size = 17, hjust = 0.5, face = 'bold'),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) + NoLegend() +
  scale_fill_discrete(name= 'Method')

#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/liver/plots/doub_detection_liver.png", width = 12, height = 8.5, units = 'in', res = 600)
p1
dev.off()

###################################Doublet finder#######################
Obj<- invitro_obj
sweep.res.list_Obj <-DoubletFinder:: paramSweep(Obj, PCs = 1:20, sct = FALSE)
sweep.stats_Obj <- summarizeSweep(sweep.res.list_Obj, GT = FALSE)
bcmvn_Obj <- find.pK(sweep.stats_Obj)
glimpse(bcmvn_Obj)
bcmvn_Obj %>% filter(BCmetric==max(BCmetric)) ## identify pk value (0.3)
homotypic.prop <- modelHomotypic(Obj@meta.data$sorting_scheme)
nExp_poi <- round(0.08*nrow(Obj@meta.data))  ## Assuming 8% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_Obj <- doubletFinder(Obj, PCs = 1:20, pN = 0.25, pK = 0.3, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
table(seu_Obj$DF.classifications_0.25_0.3_798)
DimPlot(seu_Obj, reduction = 'umap', label = T, group.by = 'DF.classifications_0.25_0.3_798')
invitro_obj$Doublet_finder <- seu_Obj$DF.classifications_0.25_0.3_798
invitro_obj$Doublet_finder_score <- seu_Obj$pANN_0.25_0.3_798
saveRDS(invitro_obj, '/mnt/8TB/users/shameed/shameed/Doublet predictions/invitro_obj.rds')


#################invivo
Obj<- invivo_obj
sweep.res.list_Obj <-DoubletFinder:: paramSweep(Obj, PCs = 1:20, sct = FALSE)
sweep.stats_Obj <- summarizeSweep(sweep.res.list_Obj, GT = FALSE)
bcmvn_Obj <- find.pK(sweep.stats_Obj)
glimpse(bcmvn_Obj)
bcmvn_Obj %>% filter(BCmetric==max(BCmetric)) ## identify pk value (0.19)
homotypic.prop <- modelHomotypic(Obj@meta.data$sorting_scheme)
nExp_poi <- round(0.07*nrow(Obj@meta.data))  ## Assuming 7% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_Obj <- doubletFinder(Obj, PCs = 1:20, pN = 0.25, pK = 0.19, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
table(seu_Obj$DF.classifications_0.25_0.19_542)
DimPlot(seu_Obj, reduction = 'umap', label = T, group.by = 'DF.classifications_0.25_0.19_542')
invivo_obj$Doublet_finder <- seu_Obj$DF.classifications_0.25_0.19_542
invivo_obj$Doublet_finder_score <- seu_Obj$pANN_0.25_0.19_542
saveRDS(invivo_obj, '/mnt/8TB/users/shameed/shameed/Doublet predictions/invivo_obj.rds')

#############confusion matrix
actual <- ifelse(invitro_obj$sorting_scheme=='CD11c+' | 
                   invitro_obj$sorting_scheme=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invitro_obj$Doublet_finder=='Doublet', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_DoubletFinder_invitro.rds')
############confusion matrix
actual <- ifelse(invivo_obj$cell_class=='CD11c+' | 
                   invivo_obj$cell_class=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invivo_obj$Doublet_finder=='Doublet', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_DoubletFinder_invivo.rds')


##############################scDblFinder#########################
sce <- as.SingleCellExperiment(invitro_obj)
sce <- scDblFinder(sce)
table(sce$scDblFinder.class)
invitro_obj$scbDB_score <- sce$scDblFinder.score
invitro_obj$scbDB_class <- sce$scDblFinder.class
saveRDS(invitro_obj, '/mnt/8TB/users/shameed/shameed/Doublet predictions/invitro_obj.rds')

###########confusion matrix
actual <- ifelse(invitro_obj$sorting_scheme=='CD11c+' | 
                   invitro_obj$sorting_scheme=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invitro_obj$scbDB_class=='doublet', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_scdb_invitro.rds')
################invivo
library(scDblFinder)
sce <- as.SingleCellExperiment(invivo_obj)
sce <- scDblFinder(sce)
table(sce$scDblFinder.class)
invivo_obj$scbDB_score <- sce$scDblFinder.score
invivo_obj$scbDB_class <- sce$scDblFinder.class

saveRDS(invivo_obj, '/mnt/8TB/users/shameed/shameed/Doublet predictions/invivo_obj.rds')
actual <- ifelse(invivo_obj$cell_class=='CD11c+' | 
                   invivo_obj$cell_class=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invivo_obj$scbDB_class=='doublet', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_scdb_invivo.rds')

#####################SCRUBLET###################################
mat <- t(invitro_obj@assays$RNA3$data)
write.csv(mat, '/mnt/8TB/users/shameed/shameed/Doublet predictions/mat_invitro.csv')
Scrublet_results_invitro <- read_csv("/mnt/8TB/users/shameed/shameed/Doublet predictions/Scrublet_results_invitro.csv")
table(Scrublet_results_invivo$predicted_doublet)
hist(Scrublet_results_invitro$doublet_score, breaks = 50,
     main = "Scrublet Doublet Scores",
     xlab = "Doublet Score",
     col = "lightblue")

#table(Scrublet_results_invitro$doublet_score > 0.2)

scrub_invt <- Scrublet_results_invitro

scrub_invt$scrublet_predic <- ifelse(scrub_invt$doublet_score > 0.2, 'doublet', 'singlet')
invitro_obj$scrublet_score <- scrub_invt$doublet_score
invitro_obj$scrublet_class <- scrub_invt$scrublet_predic

##########confusion matrix
actual <- ifelse(invitro_obj$sorting_scheme=='CD11c+' | 
                   invitro_obj$sorting_scheme=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invitro_obj$scrublet_class=='doublet', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_scrublet_invitro.rds')

#########################invivo
mat <- t(invivo_obj@assays$RNA3$data)
write.csv(mat, '/mnt/8TB/users/shameed/shameed/Doublet predictions/mat_invivo.csv')
Scrublet_results_invivo <- read_csv("/mnt/8TB/users/shameed/shameed/Doublet predictions/Scrublet_results_invivo.csv")

hist(Scrublet_results_invivo$doublet_score, breaks = 50,
     main = "Scrublet Doublet Scores",
     xlab = "Doublet Score",
     col = "lightblue")
summary(Scrublet_results_invivo$doublet_score)

#table(Scrublet_results_invivo$doublet_score > 0.2)

scrub_invt <- Scrublet_results_invivo

scrub_invt$scrublet_predic <- ifelse(scrub_invt$doublet_score > 0.2, 'doublet', 'singlet')
invivo_obj$scrublet_score <- scrub_invt$doublet_score
invivo_obj$scrublet_class <- scrub_invt$scrublet_predic

################
actual <- ifelse(invivo_obj$cell_class=='CD11c+' | 
                   invivo_obj$cell_class=='TCRb+', 'No', 'Yes')
predicted <- ifelse(invivo_obj$scrublet_class=='doublet', 'Yes', 'No')
my_tab <- table(actual, predicted)
matcon <-confusionMatrix(my_tab, positive = 'Yes')
matcon
saveRDS(matcon, 'confmat_scrublet_invivo.rds')

###############model comparison and plotting######################
################################################### DC-T
conf_df <- rbind(confmat_scrublet_invitro$byClass[1:4], 
                 confmat_DoubletFinder_invitro$byClass[1:4],
                 confmat_scdb_invitro$byClass[1:4],
                 confmat_invitro_ulm$byClass[1:4])
my_meth <- c('Scrublet', 'DoubletFinder', 'scDblFinder', 'ULMnet')
conf_df<- as.data.frame(conf_df)
conf_df$Method <- my_meth

conf_df <- conf_df[, c('Method', 'Sensitivity', "Specificity",  "Pos Pred Value", "Neg Pred Value")]
colnames(conf_df) <- c('Method', 'Precision', "Neg pred",  "Sensitivity", "Specificity")

############invivo
inv_df <- rbind(confmat_scrublet_invivo$byClass[1:4], 
                 confmat_DoubletFinder_invivo$byClass[1:4],
                 confmat_scdb_invivo$byClass[1:4],
                 confmat_invivo_ulm$byClass[1:4])
my_meth <- c('Scrublet', 'DoubletFinder', 'scDblFinder', 'ULMnet')
inv_df<- as.data.frame(inv_df)
inv_df$Method <- my_meth

inv_df <- inv_df[, c('Method', 'Sensitivity', "Specificity",  "Pos Pred Value", "Neg Pred Value")]
colnames(inv_df) <- c('Method', 'Precision', "Neg pred",  "Sensitivity", "Specificity")

 #############plots
library(tidyverse)
data <- conf_df %>% 
  dplyr:: select(c(Method, Sensitivity, Specificity, Precision)) %>%
  pivot_longer(cols = -c(Method), names_to = 'metrics', values_to = 'percentage') %>%
  mutate(percentage = percentage * 100)

p1<- ggplot(data, aes(x =factor(metrics, levels = c('Sensitivity', 'Precision', 'Specificity')), y = percentage, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.5) +  # Bar position dodge width should match text
  geom_text(aes(label = round(percentage, 1)), 
            position = position_dodge(width = 0.8),  # Ensure the same width in dodge for text
            vjust = -0.5, size = 4.5) +  # Adjust vjust for better alignment
  labs(title = "Performance Comparison by Method (in vitro)", y = "Percentage", x = '', fill = 'Methods') +
  theme_minimal() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, hjust = 0.5, face = 'bold'),
        legend.key.size = unit(1.5, "cm"),
        legend.title = element_text(size = 17, hjust = 0.5, face = 'bold'),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name= 'Method')

#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/stat_invitro_2.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()

##############invivo
data <- inv_df %>% 
  dplyr:: select(c(Method, Sensitivity, Specificity, Precision)) %>%
  pivot_longer(cols = -c(Method), names_to = 'metrics', values_to = 'percentage') %>%
  mutate(percentage = percentage * 100)

p2<- ggplot(data, aes(x =factor(metrics, levels = c('Sensitivity', 'Precision', 'Specificity')), y = percentage, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.5) +  # Bar position dodge width should match text
  geom_text(aes(label = round(percentage, 1)), 
            position = position_dodge(width = 0.8),  # Ensure the same width in dodge for text
            vjust = -0.5, size = 4.5) +  # Adjust vjust for better alignment
  labs(title = "Performance Comparison by Method (in vivo)", y = "Percentage", x = '', fill = 'Methods') +
  theme_minimal() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, hjust = 0.5, face = 'bold'),
        legend.key.size = unit(1.5, "cm"),
        legend.title = element_text(size = 17, hjust = 0.5, face = 'bold'),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name= 'Method')

png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/stat_invivo_2.png", width = 15, height = 10.5, units = 'in', res = 600)
p2
dev.off()
###########
align_fig <-ggarrange(p1, p2, ncol = 2, nrow=1, common.legend = T, legend="right")
print(align_fig)

png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/stat_DC_T_2.png", width = 20, height = 10.5, units = 'in', res = 600)
align_fig
dev.off()

##