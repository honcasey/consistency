library(PharmacoGx); library(dplyr); library(reshape2); library(survival); library(VennDiagram)
# downloadPSet("GDSC_2013"); downloadPSet("CCLE_2013")

# ==============================================================
# intersecting CCLE and GDSC for common drugs and cell lines
# ==============================================================
load("~/capsule/code/PSets/CCLE_2013.RData"); load("~/capsule/code/PSets/GDSC_2013.RData")
intersected <- intersectPSet(c(CCLE, GDSC), intersectOn = c("drugs", "cell.lines"))
# ^ 15 drugs and 514 cell lines in common!
ccle <- intersected$CCLE
gdsc <- intersected$GDSC

# venn diagram of common cell lines
dir.create("~/capsule/results/plots")
dir.create("~/capsule/results/plots/venn")
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@cell), area2=nrow(GDSC@cell), cross.area=nrow(ccle@cell), col = "black", cex=1.5, cat.cex=1, cat.col = c("black", "black"))
pdf("~/capsule/results/plots/venn/cell_intersection.pdf", height=4, width=4)
grid::grid.draw(venn.plot)
dev.off()
rm(venn.plot)

# venn diagram of common drugs
venn.plot2 <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@drug), area2=nrow(GDSC@drug), cross.area=nrow(ccle@drug), col = "black" ,cex=1.5, cat.cex=1, cat.col = c("black", "black"))
pdf("~/capsule/results/plots/venn/drug_intersection.pdf", height=4, width=4)
grid::grid.draw(venn.plot2)
dev.off()    
rm(venn.plot2)

rm(CCLE)
rm(GDSC)
rm(intersected)

coef <- c("pearson", "pearson p-value", "spearman", "spearman p-value")
common_drugs <- ccle@drug

# ====================================================
# loading manually curated list of drug information
# ====================================================
druglist <- as.data.frame(read.csv("PSets/drug_info.csv")) 
druglist[] <- lapply(druglist, as.character) # change class of dataframe columns from vector to character
rownames(druglist) <- rownames(ccle@drug)
# rownames(druglist) <- druglist$drug_name 

# ====================================================
# sort drugs into broad spectrum/targeted classes
# ====================================================
classes = data.frame(drug = rownames(ccle@drug), row.names = rownames(ccle@drug))
classes$BroadSpectrum_or_Targeted <- NA

for (drug in rownames(druglist)) {
  classes[drug, "BroadSpectrum_or_Targeted"] <- druglist[drug, "ChEMBL.drug.mechanism"]
}

# ====================================================
#       sort drugs into clue.io defined PCLs
# ====================================================
classes$PCL <- NA
for (drug in rownames(druglist)) {
  classes[drug, "PCL"] <- druglist[drug, "CMap.PCL"]
}

# ==============================================
#       get published ic50 and auc
# ==============================================
ccle_ic50 = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published"))

ccle_auc = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "auc_published"))

gdsc_ic50 = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published"))

gdsc_auc = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published"))

# ====================================================
# correlation of published IC50 between CCLE and GDSC 
# ====================================================
ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ccle_ic50))), row.names = colnames(ccle_ic50))
colnames(ic50.cor) <- coef
for (cell.line in rownames(ic50.cor)) {
  tryCatch(
    expr = {
      pearson.cor <- cor.test(x = ccle_ic50[, cell.line], y = gdsc_ic50[, cell.line], method = 'pearson', use = 'pairw')
  ic50.cor[cell.line, "pearson"] <- pearson.cor$estimate
  ic50.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ccle_ic50[, cell.line], y = gdsc_ic50[, cell.line], method = 'spearman', use = 'pairwise.complete.obs')
  ic50.cor[cell.line, "spearman"] <- spearman.cor$estimate
  ic50.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
    },
  error = {function(e){
    # message(cell.line, " does not have enough finite observations")
  }
  }
  )
}

# ====================================================
# correlation of published AUC between CCLE and GDSC (by cell line)
# ====================================================
auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ccle_auc))), row.names = colnames(ccle_auc))
colnames(auc.cor) <- coef

for (cell.line in rownames(auc.cor)) {
  tryCatch(
    expr = {
      pearson.cor <- cor.test(x = ccle_auc[, cell.line], y = gdsc_auc[, cell.line], method = 'pearson', use = 'pairw')
      auc.cor[cell.line, "pearson"] <- pearson.cor$estimate
      auc.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
      
      spearman.cor <- cor.test(x = ccle_auc[, cell.line], y = gdsc_auc[, cell.line], method = 'spearman', use = 'pairwise.complete.obs')
      auc.cor[cell.line, "spearman"] <- spearman.cor$estimate
      auc.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
    },
    error = {function(e){
      # message(cell.line, " does not have enough finite observations")
    }}
  )
}

# ====================================================
#    sort drugs into broad-spectrum vs. targeted 
# ====================================================
broad_spectrum_drugs <- classes[grep("broad-spectrum", classes$BroadSpectrum_or_Targeted), "drug"] 
targeted_drugs <- classes[grep("targeted", classes$BroadSpectrum_or_Targeted), "drug"] 

# ====================================================
#               sort drugs into PCLs 
# ====================================================
# make list of all PCLs and the drugs in each one
pcl_count <- dplyr::select(classes, PCL)
pcl_count$drug <- rownames(pcl_count)
pcls <- (reshape2::dcast(pcl_count, drug~PCL))
pcls <- pcls[, -1]
rm(pcl_count)
rm(classes)

pcl_list = list()
for (pcl in colnames(pcls)) {
  pcl_list[[pcl]] <- subset(pcls[, pcl], (!is.na(pcls[, pcl]))) #get the drugs in that PCL
}

# ====================================================
#   recompute sens measures after binary sorting 
# ====================================================
# ============= BROAD-SPECTRUM DRUGS ================
# cor of ic50 for each cell line + drugs that are broad-spectrum
ccle_ic50_brsp <- ccle_ic50[broad_spectrum_drugs, ] 
gdsc_ic50_brsp <- gdsc_ic50[broad_spectrum_drugs, ] 

ic50.brsp.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ccle_ic50))), row.names = colnames(ccle_ic50))
colnames(ic50.brsp.cor) <- coef

for (cell.line in rownames(ic50.brsp.cor)) {
  tryCatch(
    expr = {
      pearson.cor <- cor.test(x = ccle_ic50_brsp[, cell.line], y = gdsc_ic50_brsp[, cell.line], method = 'pearson', use = 'pairw')
      ic50.brsp.cor[cell.line, "pearson"] <- pearson.cor$estimate
      ic50.brsp.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
      
      spearman.cor <- cor.test(x = ccle_ic50_brsp[, cell.line], y = gdsc_ic50_brsp[, cell.line], method = 'spearman', use = 'pairw')
      ic50.brsp.cor[cell.line, "spearman"] <- spearman.cor$estimate
      ic50.brsp.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
    },
    error = {function(e){
      # message(cell.line, " does not have enough finite observations")
    }}
  )
}

# cor of auc for each cell line + drugs that are broad-spectrum
ccle_auc_brsp <- ccle_auc[broad_spectrum_drugs, ]
gdsc_auc_brsp <- gdsc_auc[broad_spectrum_drugs, ] 

auc.brsp.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ccle_ic50))), row.names = colnames(ccle_ic50))
colnames(auc.brsp.cor) <- coef

for (cell_line in rownames(auc.brsp.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = ccle_auc_brsp[, cell.line], y = gdsc_auc_brsp[, cell.line], method = 'pearson', use = 'pairw')
  auc.brsp.cor[cell.line, "pearson"] <- pearson.cor$estimate
  auc.brsp.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ccle_auc_brsp[, cell.line], y = gdsc_auc_brsp[, cell.line], method = 'spearman', use = 'pairwise.complete.obs')
  auc.brsp.cor[cell.line, "spearman"] <- spearman.cor$estimate
  auc.brsp.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
    },
  error = {function(e) {
    # message(cell.line, " does not have enough finite observations")
  }}
  )
}

# ============= TARGETED DRUGS ================
# cor of ic50 for each cell line + drugs that are targeted
ccle_ic50_targ <- ccle_ic50[targeted_drugs, ]
gdsc_ic50_targ <- gdsc_ic50[targeted_drugs, ] 

ic50.targ.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ccle_ic50))), row.names = colnames(ccle_ic50))
colnames(ic50.targ.cor) <- coef

for (cell.line in rownames(ic50.targ.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = ccle_ic50_targ[, cell.line], y = gdsc_ic50_targ[, cell.line], method = 'pearson', use = 'pairw')
  ic50.targ.cor[cell.line, "pearson"] <- pearson.cor$estimate
  ic50.targ.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ccle_ic50_targ[, cell.line], y = gdsc_ic50_targ[, cell.line], method = 'spearman', use = 'pairw')
  ic50.targ.cor[cell.line, "spearman"] <- spearman.cor$estimate
  ic50.targ.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
    },
  error = {function(e) {
    # message(cell.line, " does not have enough finite observations")
  }}
  )
}

# cor of auc for each cell line + drugs that are targeted
ccle_auc_targ <- ccle_auc[targeted_drugs, ]
gdsc_auc_targ <- gdsc_auc[targeted_drugs, ] 

auc.targ.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ccle_auc))), row.names = colnames(ccle_auc))
colnames(auc.targ.cor) <- coef

for (cell.line in rownames(auc.targ.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = ccle_auc_targ[, cell.line], y = gdsc_auc_targ[, cell.line], method = 'pearson', use = 'pairw')
  auc.targ.cor[cell.line, "pearson"] <- pearson.cor$estimate
  auc.targ.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ccle_auc_targ[, cell.line], y = gdsc_auc_targ[, cell.line], method = 'spearman', use = 'pairw')
  auc.targ.cor[cell.line, "spearman"] <- spearman.cor$estimate
  auc.targ.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
    },
  error = {function(e) {
    # message(cell.line, " does not have enough finite observations")
  }}
  )
}

# ====================================================
#   recompute sens measures after PCL sorting 
# ====================================================
min_drugs = 3 # keeping only pcls that contain at least 3 drugs per class
pcl_temp <- data.frame(matrix(NA, nrow = min_drugs)) 
for (pcl in colnames(pcls)) {
  if (!length(pcl_list[[pcl]]) == 1) {
    pcl_temp[, pcl] <- pcl_list[[pcl]]
  }
}
pcl_temp <- pcl_temp[, -1]

# correlation of auc in drugs by pcl for each cell line 
auc.pcl.cor = data.frame(matrix(NA, nrow = length(colnames(ccle_auc))), row.names = colnames(ccle_auc))

for (pcl in colnames(pcl_temp)) {
  auc.pcl.cor[, paste0(pcl, ".pearson")] <- NA
    auc.pcl.cor[, paste0(pcl, ".pearson p-value")] <- NA
  auc.pcl.cor[, paste0(pcl, ".spearman")] <- NA
    auc.pcl.cor[, paste0(pcl, ".spearman p-value")] <- NA
}
auc.pcl.cor <- auc.pcl.cor[, -1]

for (pcl in colnames(pcl_temp)) {
  ccle_auc_temp <- ccle_auc[pcl_list[[pcl]], ] # subset by the drugs in the PCL
  gdsc_auc_temp <- gdsc_auc[pcl_list[[pcl]], ]
  for (cell.line in rownames(auc.pcl.cor)) {
    tryCatch(
      expr = {
        pearson.cor <- cor.test(x = ccle_auc_temp[, cell.line], y = gdsc_auc_temp[, cell.line], method = 'pearson', use = 'pairw')
        auc.pcl.cor[cell.line, paste0(pcl, ".pearson")] <- pearson.cor$estimate
        auc.pcl.cor[cell.line, paste0(pcl, ".pearson p-value")] <- pearson.cor$p.value
  
        spearman.cor <- cor.test(x = ccle_auc_temp[, cell.line], y = gdsc_auc_temp[, cell.line], method = 'spearman', use = 'pairw')
        auc.pcl.cor[cell.line, paste0(pcl, ".spearman")] <- spearman.cor$estimate
        auc.pcl.cor[cell.line, paste0(pcl, ".spearman p-value")] <- spearman.cor$p.value
        },
      error = {function(e) {
        # message(cell.line, " does not have enough finite observations")
      }}
    )
  }
}
    

# correlation of ic50 in drugs by pcl for each cell line
ic50.pcl.cor = data.frame(matrix(NA, nrow = length(colnames(ccle_ic50))), row.names = colnames(ccle_ic50))

for (pcl in colnames(pcl_temp)) {
  ic50.pcl.cor[, paste0(pcl, ".pearson")] <- NA
    ic50.pcl.cor[, paste0(pcl, ".pearson p-value")] <- NA
  ic50.pcl.cor[, paste0(pcl, ".spearman")] <- NA
    ic50.pcl.cor[, paste0(pcl, ".spearman p-value")] <- NA
}
ic50.pcl.cor <- ic50.pcl.cor[, -1]

for (pcl in colnames(pcl_temp)) {
  ccle_ic50_temp <- ccle_ic50[pcl_list[[pcl]], ] # subset by the drugs in the PCL
  gdsc_ic50_temp <- gdsc_ic50[pcl_list[[pcl]], ]
  for (cell.line in rownames(ic50.pcl.cor)) {
    tryCatch(
      expr = {
        pearson.cor <- cor.test(x = ccle_ic50_temp[, cell.line], y = gdsc_ic50_temp[, cell.line], method = 'pearson', use = 'pairw')
        ic50.pcl.cor[cell.line, paste0(pcl, ".pearson")] <- pearson.cor$estimate
        ic50.pcl.cor[cell.line, paste0(pcl, ".pearson p-value")] <- pearson.cor$p.value
  
        spearman.cor <- cor.test(x = ccle_ic50_temp[, cell.line], y = gdsc_ic50_temp[, cell.line], method = 'spearman', use = 'pairw')
        ic50.pcl.cor[cell.line, paste0(pcl, ".spearman")] <- spearman.cor$estimate
        ic50.pcl.cor[cell.line, paste0(pcl, ".spearman p-value")] <- spearman.cor$p.value
      },
    error ={function(e) {
      # message(cell.line, " does not have enough finite observations")
    }}
  )
  }
}

# ====================================================
#       removing inconsistent cell lines
# ====================================================
min = 0.5
# ===== inconsistent based on binary filtering ====
rmv.ic50.binary.cor <- rbind(ic50.brsp.cor[which(ic50.brsp.cor$spearman > min), ], ic50.targ.cor[which(ic50.targ.cor$spearman > min), ])
bin.ic50.cc <- rownames(rmv.ic50.binary.cor)
rmv.auc.binary.cor <- rbind(auc.brsp.cor[which(auc.brsp.cor$spearman > min), ], auc.targ.cor[which(auc.targ.cor$spearman > min), ])
bin.auc.cc <- rownames(rmv.auc.binary.cor)

# ===== inconsistent based on pcls =====
rmv.ic50.pcl.cor <- rbind(ic50.pcl.cor[which(ic50.pcl.cor$`SRC inhibitor.spearman` > min), ], ic50.pcl.cor[which(ic50.targ.cor$`MEK inhibitor.spearman` > min), ])
pcl.ic50.cc <- rownames(rmv.ic50.pcl.cor)
rmv.auc.pcl.cor <- rbind(auc.pcl.cor[which(auc.pcl.cor$`MEK inhibitor.spearman` > min), ], auc.pcl.cor[which(auc.pcl.cor$`SRC inhibitor.spearman` > min), ])
pcl.auc.cc <- rownames(rmv.auc.pcl.cor)

# ====================================================
#    remove inconsistent cell lines within GDSC 
# ====================================================
# filtering based on binary classifying
gdsc_ic50_cc_binary = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published", cell.lines = bin.ic50.cc))
gdsc_auc_cc_binary = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published", cell.lines = bin.auc.cc))

# filtering based on PCL classifying
gdsc_ic50_cc_pcl = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published", cell.lines = pcl.ic50.cc))
gdsc_auc_cc_pcl = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published", cell.lines = pcl.auc.cc))

# ====================================================
#    remove inconsistent cell lines within CCLE 
# ====================================================
# filtering based on binary classifying
ccle_ic50_cc_binary = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published", cell.lines = bin.ic50.cc))
ccle_auc_cc_binary = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "auc_published", cell.lines = bin.auc.cc))

# filtering based on PCL classifying
ccle_ic50_cc_pcl = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published", cell.lines = pcl.ic50.cc))
ccle_auc_cc_pcl = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "auc_published", cell.lines = pcl.auc.cc))


# ======================================================================================
#       correlation of published sens measures between CCLE and GDSC (by drug)
# ======================================================================================
orig.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_ic50))), row.names = rownames(ccle_ic50))
colnames(orig.ic50.cor) <- coef

for (drug in rownames(orig.ic50.cor)) {
  pearson.cor <- cor.test(x = as.numeric(ccle_ic50[drug, ]), y = as.numeric(gdsc_ic50[drug, ]), method = 'pearson', use = 'pairw')
  orig.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
  orig.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_ic50[drug, ]), y = as.numeric(gdsc_ic50[drug, ]), method = 'spearman', use = 'pairw')
  orig.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
  orig.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

orig.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_auc))), row.names = rownames(ccle_auc))
colnames(orig.auc.cor) <- coef

for (drug in rownames(orig.auc.cor)) {
  pearson.cor <- cor.test(x = as.numeric(ccle_auc[drug, ]), y = as.numeric(gdsc_auc[drug, ]), method = 'pearson', use = 'pairw')
  orig.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  orig.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_auc[drug, ]), y = as.numeric(gdsc_auc[drug, ]), method = 'spearman', use = 'pairw')
  orig.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  orig.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}
marray::write.xls(orig.ic50.cor, "~/capsule/results/orig.ic50.cor.xls", row.names = TRUE, col.names = TRUE)
marray::write.xls(orig.auc.cor, "~/capsule/results/orig.auc.cor.xls", row.names = TRUE, col.names = TRUE)

# ======================================================================================
#   correlation between CCLE and GDSC after sorting into broad-spectrum/targeted
# ======================================================================================
# ==== correlation after binary sorting ====
bin.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_ic50_cc_binary))), row.names = rownames(ccle_ic50_cc_binary))
colnames(bin.ic50.cor) <- coef

for (drug in rownames(bin.ic50.cor)) {
  pearson.cor <- cor.test(x = as.numeric(ccle_ic50_cc_binary[drug, ]), y = as.numeric(gdsc_ic50_cc_binary[drug, ]), method = 'pearson', use = 'pairw')
  bin.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
  bin.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_ic50_cc_binary[drug, ]), y = as.numeric(gdsc_ic50_cc_binary[drug, ]), method = 'spearman', use = 'pairw')
  bin.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
  bin.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

bin.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_auc_cc_binary))), row.names = rownames(ccle_auc_cc_binary))
colnames(bin.ic50.cor) <- coef

for (drug in rownames(bin.auc.cor)) {
  pearson.cor <- cor.test(x = as.numeric(ccle_auc_cc_binary[drug, ]), y = as.numeric(gdsc_auc_cc_binary[drug, ]), method = 'pearson', use = 'pairw')
  bin.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  bin.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_auc_cc_binary[drug, ]), y = as.numeric(gdsc_auc_cc_binary[drug, ]), method = 'spearman', use = 'pairw')
  bin.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  bin.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# ======================================================================================
#           Harrell's CI after sorting into broad-spectrum/targeted
# ======================================================================================
binary.ic50.cor <- merge(orig.ic50.cor, bin.ic50.cor, by = 0, all = TRUE)
marray::write.xls(binary.ic50.cor, "~/capsule/results/binary.ic50.cor.xls", row.names = TRUE, col.names = TRUE)
binary.ic50.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = binary.ic50.cor)
saveRDS(binary.ic50.concordance, "~/capsule/results/binary.ic50.concordance.rds")

# AUC
binary.auc.cor <- merge(orig.auc.cor, bin.auc.cor, by = 0, all = TRUE)
marray::write.xls(binary.auc.cor, "~/capsule/results/binary.auc.cor.xls", row.names = TRUE, col.names = TRUE)
binary.auc.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = binary.auc.cor)
saveRDS(binary.auc.concordance, "~/capsule/results/binary.auc.concordance.rds")

# =============================================================================
#       correlation between CCLE and GDSC after sorting into PCLs
# =============================================================================
pcl.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_ic50_cc_pcl))), row.names = rownames(ccle_ic50_cc_pcl))
colnames(pcl.ic50.cor) <- coef

for (drug in rownames(pcl.ic50.cor)) {
  pearson.cor <- cor.test(x = as.numeric(ccle_ic50_cc_pcl[drug, ]), y = as.numeric(gdsc_ic50_cc_pcl[drug, ]), method = 'pearson', use = 'pairw')
  pcl.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
  pcl.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_ic50_cc_pcl[drug, ]), y = as.numeric(gdsc_ic50_cc_pcl[drug, ]), method = 'spearman', use = 'pairw')
  pcl.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
  pcl.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

pcl.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_auc_cc_pcl))), row.names = rownames(ccle_auc_cc_pcl))
colnames(pcl.auc.cor) <- coef

for (drug in rownames(pcl.auc.cor)) {
  pearson.cor <- cor.test(x = as.numeric(ccle_auc_cc_pcl[drug, ]), y = as.numeric(gdsc_auc_cc_pcl[drug, ]), method = 'pearson', use = 'pairw')
  pcl.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  pcl.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_auc_cc_pcl[drug, ]), y = as.numeric(gdsc_auc_cc_pcl[drug, ]), method = 'spearman', use = 'pairw')
  pcl.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  pcl.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# ====================================================================================
#               Harrell's CI after sorting into PCLs
# ====================================================================================
orig.pcl.ic50.cor <- merge(orig.ic50.cor, pcl.ic50.cor, by = 0, all = TRUE)
marray::write.xls(orig.pcl.ic50.cor, "~/capsule/results/pcl.ic50.cor.xls", row.names = TRUE, col.names = TRUE)
pcl.ic50.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.pcl.ic50.cor)
saveRDS(pcl.ic50.concordance, "~/capsule/results/pcl.ic50.concordance.rds")

orig.pcl.auc.cor <- merge(orig.auc.cor, pcl.auc.cor, by = 0, all = TRUE)
marray::write.xls(orig.pcl.auc.cor, "~/capsule/results/pcl.auc.cor.xls", row.names = TRUE, col.names = TRUE)
pcl.auc.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.pcl.auc.cor)
saveRDS(pcl.auc.concordance, "~/capsule/results/pcl.auc.concordance.rds")

# ======================================================
# 1) Bar plot representing the Spearman's rank correlation coefficient for IC50 drug sensitivity measures
ss <- as.matrix(t(cbind(orig.ic50.cor$spearman, bin.ic50.cor$spearman, pcl.ic50.cor$spearman)))
ss[!is.na(ss) & ss < 0] <- 0
names(ss) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/all_ic50_bar_plot.pdf")
mp <- barplot(ss, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(ss), v=0.9), each=2), ylab=expression("Spearman's rank correlation coefficient r"[s]), density=c(100,10,20) , angle=c(0,45,90), main="Spearman's rank correlation coefficient of IC50 of original values \nand after both filtering methods", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ", "pcl"), density=c(100,10,20), angle=c(0,45,90), bty="n", cex=0.75)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(ss)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# 2) Bar plot representing the Spearman's rank correlation coefficient for AUC drug sensitivity measures; significance is reported using an asterisk if two‐sided P value <0.05.
aa <- as.matrix(t(cbind(orig.auc.cor$spearman, orig.auc.cor$spearman, pcl.auc.cor$spearman)))
aa[!is.na(aa) & aa < 0] <- 0
names(aa) <- rownames(orig.auc.cor)

pdf("~/capsule/results/plots/all_auc_bar_plot.pdf")
bp <- barplot(aa, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(aa), v=0.9), each=2), ylab=expression("Spearman's rank correlation coefficient r"[s]), density=c(100,10,20) , angle=c(0,45,90), main="Spearman's rank correlation coefficient of AUC of original values \nand after both filtering methods", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ", "pcl"), density=c(100,10,20), angle=c(0,45,90), bty="n", cex=0.75)
text(x=apply(bp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=1, labels=toupper(names(aa)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# 3) scatter plot (for each drug) reporting AUC of all cell line
pdf("~/capsule/results/plots/drugs_auc_scatterplot.pdf", height = 14, width = 14)
par(mfrow=c(4,4))
for (drug in rownames(gdsc_auc)) {
  plot(as.numeric(ccle_auc[drug, ]), as.numeric(gdsc_auc[drug, ]), main = drug, xlab = "CCLE", ylab = "GDSC", type = "p", col = "royalblue2")
}
mtext("Comparison of original AUC values", outer=TRUE,  cex=1, line=-1.5)
dev.off()

#scatter plot for each drug reporting AUC with cell lines removed by binary filtering
pdf("~/capsule/results/plots/bin_auc_scatterplot.pdf", height = 14, width = 14)
par(mfrow=c(4,4))
for (drug in rownames(gdsc_auc_cc_binary)) {
  plot(as.numeric(ccle_auc_cc_binary[drug, ]), as.numeric(gdsc_auc_cc_binary[drug, ]), main = drug, xlab = "CCLE", ylab = "GDSC", type = "p", col = "royalblue2")
}
mtext("Comparison of AUC values by binary drug classification", outer=TRUE,  cex=1, line=-1.5)
dev.off()

#scatter plot for each drug reporting AUC with cell lines removed by pcl filtering
pdf("~/capsule/results/plots/pcl_auc_scatterplot.pdf", height = 14, width = 14)
par(mfrow=c(4,4))
for (drug in rownames(gdsc_auc_cc_pcl)) {
  plot(as.numeric(ccle_auc_cc_pcl[drug, ]), as.numeric(gdsc_auc_cc_pcl[drug, ]), main = drug, xlab = "CCLE", ylab = "GDSC", type = "p", col = "royalblue2")
}
mtext("Comparison of AUC values by PCL drug classification", outer=TRUE,  cex=1, line=-1.5)
dev.off()

# 4) scatter plot for each drug reporting IC50 of all cell lines
pdf("~/capsule/results/plots/drugs_ic50_scatterplot.pdf", height = 14, width = 14)
par(mfrow=c(4,4))
for (drug in rownames(gdsc_ic50)) {
  plot(as.numeric(ccle_ic50[drug, ]), as.numeric(gdsc_ic50[drug, ]), main = drug, xlab = "CCLE", ylab = "GDSC", type = "p", col = "royalblue2")
}
mtext("Comparison of reported IC50 values in all common cell lines", outer=TRUE,  cex=1, line=-1.5)
dev.off()

#scatter plot for each drug reporting ic50 with cell lines removed by broad-spectrum filtering
pdf("~/capsule/results/plots/bin_ic50_scatterplot.pdf", height = 14, width = 14)
par(mfrow=c(4,4))
for (drug in rownames(gdsc_ic50_cc_binary)) {
  plot(as.numeric(ccle_ic50_cc_binary[drug, ]), as.numeric(gdsc_ic50_cc_binary[drug, ]), main = drug, xlab = "CCLE", ylab = "GDSC", type = "p", col = "royalblue2")
}
mtext("Comparison of reported IC50 values after filtering cell lines by binary drug classification", outer=TRUE,  cex=1, line=-1.5)
dev.off()

#scatter plot for each drug reporting ic50 with cell lines removed by pcl filtering
pdf("~/capsule/results/plots/pcl_ic50_scatterplot.pdf", height = 14, width = 14)
par(mfrow=c(4,4))
for (drug in rownames(gdsc_ic50_cc_pcl)) {
  plot(as.numeric(ccle_ic50_cc_pcl[drug, ]), as.numeric(gdsc_ic50_cc_pcl[drug, ]), main = drug, xlab = "CCLE", ylab = "GDSC", type = "p", col = "royalblue2")
}
mtext("Comparison of reported IC50 values after filtering cell lines by PCL drug classification", outer=TRUE,  cex=1, line=-1.5)
dev.off()

# 5) Bar plot representing spearman's ranks for AUC drug measure comparing original with binary
tg <- as.matrix(t(cbind(orig.auc.cor$spearman, bin.auc.cor$spearman)))
tg[!is.na(tg) & tg < 0] <- 0
names(tg) <- rownames(orig.auc.cor)

pdf("~/capsule/results/plots/bin_auc_cor_bar_plot.pdf")
tgb <- barplot(tg, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(tg), v=0.9), each=2), ylab=expression("Spearman's rank correlation coefficient r"[s]), density=c(100,10) , angle=c(0,45), main="Spearman's rank correlation coefficient of AUC values \nafter filtering cell lines by binary drug classification", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ"), fill=c("black", "black"), density=c(100, 10), bty="n", cex=1)
text(x=apply(tgb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(tg)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 6) Bar plot representing spearman's ranks for AUC drug measure comparing original with pcls
mek <- as.matrix(t(cbind(orig.auc.cor$spearman, pcl.auc.cor$spearman)))
mek[!is.na(mek) & mek < 0] <- 0
names(mek) <- rownames(orig.auc.cor)

pdf("~/capsule/results/plots/pcl_auc_cor_bar_plot.pdf")
mb <- barplot(mek, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(mek), v=0.9), each=2), ylab=expression("Spearman's rank correlation coefficient r"[s]), density=c(100,10) , angle=c(0,45), main="Spearman's rank correlation coefficient for AUC values \nafter filtering cell lines by PCL drug classification", font.main = 1)
legend("topright", legend=c("orig", "PCL"), fill=c("black", "black"), density=c(100, 10), bty="n", cex=1)
text(x=apply(mb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(mek)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 8) Bar plot representing spearman's ranks for IC50 drug measure comparing original with binary
aa <- as.matrix(t(cbind(orig.ic50.cor$spearman, bin.ic50.cor$spearman)))
aa[!is.na(aa) & aa < 0] <- 0
names(aa) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/bin_ic50_cor_bar_plot.pdf")
aab <- barplot(aa, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(aa), v=0.9), each=2), ylab=expression("Spearman's rank correlation coefficient r"[s]), density=c(100,10) , angle=c(0,45), main="Spearman's rank correlation coefficient for IC50 values \nafter filtering cell lines by binary drug classification", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ"), fill=c("black", "black"), density=c(100, 10), bty="n", cex=1)
text(x=apply(aab, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(aa)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 9) Bar plot representing spearman's ranks for IC50 drug measure comparing original with pcl
bb <- as.matrix(t(cbind(orig.ic50.cor$spearman, pcl.ic50.cor$spearman)))
bb[!is.na(bb) & bb < 0] <- 0
names(bb) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/pcl_ic50_cor_bar_plot.pdf")
bbb <- barplot(bb, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(bb), v=0.9), each=2), ylab=expression("Spearman's rank correlation coefficient r"[s]), angle=c(0, 45), density=c(100, 10), main="Spearman's rank correlation coefficient for IC50 values \nafter filtering cell lines by PCL drug classification", font.main = 1)
legend("topright", legend=c("orig", "PCL"), fill=c("black", "black"), density=c(100, 10), bty="n", cex=1)
text(x=apply(bbb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(bb)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 12) bar plot representing all concordance indexes
all_ic50_concordances = c(binary.ic50.concordance$concordance, pcl.ic50.concordance$concordance)
names(all_ic50_concordances) <- c("Broad-Spectrum/Targeted", "PCL")

pdf("~/capsule/results/plots/ic50_concordances_bar_plot.pdf")
ylim <- c(0, 1.2*max(all_ic50_concordances))
ee <- barplot(all_ic50_concordances, main = "Concordance between Pearson correlations of IC50 values", width = 0.9, ylim = ylim, ylab = "concordance index")
text(x = ee, y = all_ic50_concordances, label = signif(all_ic50_concordances, 3), pos = 3, cex = 0.8, col = "black")
dev.off()

all_auc_concordances = c(binary.auc.concordance$concordance, pcl.auc.concordance$concordance)
names(all_auc_concordances) <- c("Broad-Spectrum/Targeted", "PCL")

pdf("~/capsule/results/plots/auc_concordances_bar_plot.pdf")
ylim <- c(0, 1.2*max(all_auc_concordances))
ff <- barplot(all_auc_concordances, main = "Concordance between Pearson correlations of AUC values", width = 0.9, ylim = ylim, ylab = "concordance index")
text(x = ff, y = all_auc_concordances, label = signif(all_auc_concordances, 3), pos = 3, cex = 0.8, col = "black")
dev.off()
