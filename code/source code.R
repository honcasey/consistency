library(PharmacoGx); library(dplyr); library(reshape2); library(survival); library(VennDiagram); library(readxl)
# downloadPSet("GDSC_2013"); downloadPSet("ctrp_2013")
CTRP <- downloadPSet("CTRPv2_2015")
GDSC <- downloadPSet("GDSC_2020(v2-8.2)")

# ==============================================================
# intersecting ctrp and GDSC for common drugs and cell lines
# ==============================================================
# load("~/capsule/code/PSets/ctrp_2013.RData"); load("~/capsule/code/PSets/GDSC_2013.RData")

intersected <- intersectPSet(c(CTRP, GDSC), intersectOn = c("drugs", "cell.lines"))
ctrp <- intersected$CTRP
gdsc <- intersected$GDSC

# venn diagram of common cell lines
dir.create("~/capsule/results/plots")
dir.create("~/capsule/results/plots/venn")
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CTRP@cell), area2=nrow(GDSC@cell), cross.area=nrow(ctrp@cell), category = c("CTRP", "GDSC"), col = "black", cex=1.5, cat.cex=1, cat.col = c("black", "black"),fontfamily = rep("sans", 3))
pdf("~/capsule/results/plots/venn/cell_intersection.pdf", height=4, width=5)
grid::grid.draw(venn.plot)
dev.off()
rm(venn.plot)

# venn diagram of common drugs
venn.plot2 <- VennDiagram::draw.pairwise.venn(area1=nrow(CTRP@drug), area2=nrow(GDSC@drug), cross.area=nrow(ctrp@drug), category = c("CTRP", "GDSC"), col = "black" ,cex=1.5, cat.cex=1, cat.col = c("black", "black"), fontfamily = rep("sans", 3))
pdf("~/capsule/results/plots/venn/drug_intersection.pdf", height=4, width=5)
grid::grid.draw(venn.plot2)
dev.off()    
rm(venn.plot2)

rm(CTRP)
rm(GDSC)
rm(intersected)

coef <- c("pearson", "pearson p-value", "spearman", "spearman p-value")
common_drugs <- ctrp@drug
marray::write.xls(common_drugs, "~/capsule/results/common_drugs.xls", row.names = TRUE, col.names = TRUE)

# ====================================================
# loading manually curated list of drug information
# ====================================================
druglist <- as.data.frame(read_excel("drug_info.xls"))
druglist[] <- lapply(druglist, as.character)
rownames(druglist) <- druglist$master_cpd_id

# ====================================================
# sort drugs into broad spectrum/targeted classes
# ====================================================
classes = data.frame(drug = rownames(ctrp@drug), row.names = rownames(ctrp@drug))
classes$BroadSpectrum_or_Targeted <- NA

for (drug in rownames(druglist)) {
  classes[drug, "BroadSpectrum_or_Targeted"] <- druglist[drug, "target type"]
}

# ====================================================
#       sort drugs into clue.io defined PCLs
# ====================================================
classes$PCL <- NA
for (drug in rownames(druglist)) {
  classes[drug, "PCL"] <- druglist[drug, "CMAP PCL"]
}

# ==============================================
#       get published ic50 and aac
# ==============================================
ctrp_ic50 = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "ic50_recomputed"))

ctrp_aac = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "aac_recomputed"))

gdsc_ic50 = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_recomputed"))

gdsc_aac = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed"))

# ====================================================
# correlation of published IC50 between CTRP and GDSC 
# ====================================================
ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ctrp_ic50))), row.names = colnames(ctrp_ic50))
colnames(ic50.cor) <- coef
for (cell.line in rownames(ic50.cor)) {
  tryCatch(
    expr = {
      pearson.cor <- cor.test(x = ctrp_ic50[, cell.line], y = gdsc_ic50[, cell.line], method = 'pearson', use = 'pairw')
  ic50.cor[cell.line, "pearson"] <- pearson.cor$estimate
  ic50.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ctrp_ic50[, cell.line], y = gdsc_ic50[, cell.line], method = 'spearman', use = 'pairwise.complete.obs')
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
# correlation of published AUC between ctrp and GDSC (by cell line)
# ====================================================
auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ctrp_aac))), row.names = colnames(ctrp_aac))
colnames(auc.cor) <- coef

for (cell.line in rownames(auc.cor)) {
  tryCatch(
    expr = {
      pearson.cor <- cor.test(x = ctrp_aac[, cell.line], y = gdsc_aac[, cell.line], method = 'pearson', use = 'pairw')
      auc.cor[cell.line, "pearson"] <- pearson.cor$estimate
      auc.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
      
      spearman.cor <- cor.test(x = ctrp_aac[, cell.line], y = gdsc_aac[, cell.line], method = 'spearman', use = 'pairwise.complete.obs')
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
# rm(pcl_count)
# rm(classes)

pcl_list = list()
for (pcl in colnames(pcls)) {
  pcl_list[[pcl]] <- subset(pcls[, pcl], (!is.na(pcls[, pcl]))) #get the drugs in that PCL
}

# ====================================================
#   recompute sens measures after binary sorting 
# ====================================================
# ============= BROAD-SPECTRUM DRUGS ================
# cor of ic50 for each cell line + drugs that are broad-spectrum
ctrp_ic50_brsp <- ctrp_ic50[broad_spectrum_drugs, ] 
gdsc_ic50_brsp <- gdsc_ic50[broad_spectrum_drugs, ] 

ic50.brsp.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ctrp_ic50))), row.names = colnames(ctrp_ic50))
colnames(ic50.brsp.cor) <- coef

for (cell.line in rownames(ic50.brsp.cor)) {
  tryCatch(
    expr = {
      pearson.cor <- cor.test(x = ctrp_ic50_brsp[, cell.line], y = gdsc_ic50_brsp[, cell.line], method = 'pearson', use = 'pairw')
      ic50.brsp.cor[cell.line, "pearson"] <- pearson.cor$estimate
      ic50.brsp.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
      
      spearman.cor <- cor.test(x = ctrp_ic50_brsp[, cell.line], y = gdsc_ic50_brsp[, cell.line], method = 'spearman', use = 'pairw')
      ic50.brsp.cor[cell.line, "spearman"] <- spearman.cor$estimate
      ic50.brsp.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
    },
    error = {function(e){
      # message(cell.line, " does not have enough finite observations")
    }}
  )
}

# cor of auc for each cell line + drugs that are broad-spectrum
ctrp_aac_brsp <- ctrp_aac[broad_spectrum_drugs, ]
gdsc_aac_brsp <- gdsc_aac[broad_spectrum_drugs, ] 

auc.brsp.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ctrp_ic50))), row.names = colnames(ctrp_ic50))
colnames(auc.brsp.cor) <- coef

for (cell_line in rownames(auc.brsp.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = ctrp_aac_brsp[, cell.line], y = gdsc_aac_brsp[, cell.line], method = 'pearson', use = 'pairw')
  auc.brsp.cor[cell.line, "pearson"] <- pearson.cor$estimate
  auc.brsp.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ctrp_aac_brsp[, cell.line], y = gdsc_aac_brsp[, cell.line], method = 'spearman', use = 'pairwise.complete.obs')
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
ctrp_ic50_targ <- ctrp_ic50[targeted_drugs, ]
gdsc_ic50_targ <- gdsc_ic50[targeted_drugs, ] 

ic50.targ.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ctrp_ic50))), row.names = colnames(ctrp_ic50))
colnames(ic50.targ.cor) <- coef

for (cell.line in rownames(ic50.targ.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = ctrp_ic50_targ[, cell.line], y = gdsc_ic50_targ[, cell.line], method = 'pearson', use = 'pairw')
  ic50.targ.cor[cell.line, "pearson"] <- pearson.cor$estimate
  ic50.targ.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ctrp_ic50_targ[, cell.line], y = gdsc_ic50_targ[, cell.line], method = 'spearman', use = 'pairw')
  ic50.targ.cor[cell.line, "spearman"] <- spearman.cor$estimate
  ic50.targ.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
    },
  error = {function(e) {
    # message(cell.line, " does not have enough finite observations")
  }}
  )
}

# cor of auc for each cell line + drugs that are targeted
ctrp_aac_targ <- ctrp_aac[targeted_drugs, ]
gdsc_aac_targ <- gdsc_aac[targeted_drugs, ] 

auc.targ.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ctrp_aac))), row.names = colnames(ctrp_aac))
colnames(auc.targ.cor) <- coef

for (cell.line in rownames(auc.targ.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = ctrp_aac_targ[, cell.line], y = gdsc_aac_targ[, cell.line], method = 'pearson', use = 'pairw')
  auc.targ.cor[cell.line, "pearson"] <- pearson.cor$estimate
  auc.targ.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ctrp_aac_targ[, cell.line], y = gdsc_aac_targ[, cell.line], method = 'spearman', use = 'pairw')
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
pcl_list <- pcl_list[names(pcl_list) != "NA"] #remove drugs with no PCL

for (pcl in colnames(pcls)) {
  if (length(pcl_list[[pcl]]) >= min_drugs) {
    pcl_temp[, pcl] <- unlist(pcl_list[[pcl]])
  }
}
pcl_temp <- pcl_temp[, -1]

# correlation of auc in drugs by pcl for each cell line 
auc.pcl.cor = data.frame(matrix(NA, nrow = length(colnames(ctrp_aac))), row.names = colnames(ctrp_aac))

for (pcl in colnames(pcl_temp)) {
  auc.pcl.cor[, paste0(pcl, ".pearson")] <- NA
    auc.pcl.cor[, paste0(pcl, ".pearson p-value")] <- NA
  auc.pcl.cor[, paste0(pcl, ".spearman")] <- NA
    auc.pcl.cor[, paste0(pcl, ".spearman p-value")] <- NA
}
auc.pcl.cor <- auc.pcl.cor[, -1]

for (pcl in colnames(pcl_temp)) {
  ctrp_aac_temp <- ctrp_aac[pcl_list[[pcl]], ] # subset by the drugs in the PCL
  gdsc_aac_temp <- gdsc_aac[pcl_list[[pcl]], ]
  for (cell.line in rownames(auc.pcl.cor)) {
    tryCatch(
      expr = {
        pearson.cor <- cor.test(x = ctrp_aac_temp[, cell.line], y = gdsc_aac_temp[, cell.line], method = 'pearson', use = 'pairw')
        auc.pcl.cor[cell.line, paste0(pcl, ".pearson")] <- pearson.cor$estimate
        auc.pcl.cor[cell.line, paste0(pcl, ".pearson p-value")] <- pearson.cor$p.value
  
        spearman.cor <- cor.test(x = ctrp_aac_temp[, cell.line], y = gdsc_aac_temp[, cell.line], method = 'spearman', use = 'pairw')
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
ic50.pcl.cor = data.frame(matrix(NA, nrow = length(colnames(ctrp_ic50))), row.names = colnames(ctrp_ic50))

for (pcl in colnames(pcl_temp)) {
  ic50.pcl.cor[, paste0(pcl, ".pearson")] <- NA
    ic50.pcl.cor[, paste0(pcl, ".pearson p-value")] <- NA
  ic50.pcl.cor[, paste0(pcl, ".spearman")] <- NA
    ic50.pcl.cor[, paste0(pcl, ".spearman p-value")] <- NA
}
ic50.pcl.cor <- ic50.pcl.cor[, -1]

for (pcl in colnames(pcl_temp)) {
  ctrp_ic50_temp <- ctrp_ic50[pcl_list[[pcl]], ] # subset by the drugs in the PCL
  gdsc_ic50_temp <- gdsc_ic50[pcl_list[[pcl]], ]
  for (cell.line in rownames(ic50.pcl.cor)) {
    tryCatch(
      expr = {
        pearson.cor <- cor.test(x = ctrp_ic50_temp[, cell.line], y = gdsc_ic50_temp[, cell.line], method = 'pearson', use = 'pairw')
        ic50.pcl.cor[cell.line, paste0(pcl, ".pearson")] <- pearson.cor$estimate
        ic50.pcl.cor[cell.line, paste0(pcl, ".pearson p-value")] <- pearson.cor$p.value
  
        spearman.cor <- cor.test(x = ctrp_ic50_temp[, cell.line], y = gdsc_ic50_temp[, cell.line], method = 'spearman', use = 'pairw')
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
rmv.ic50.binary.cor <- union(ic50.brsp.cor[which(ic50.brsp.cor$pearson > min), ], ic50.targ.cor[which(ic50.targ.cor$pearson > min), ])
bin.ic50.cc <- unique(rownames(rmv.ic50.binary.cor))
rmv.auc.binary.cor <- union(auc.brsp.cor[which(auc.brsp.cor$pearson > min), ], auc.targ.cor[which(auc.targ.cor$pearson > min), ])
bin.auc.cc <- rownames(rmv.auc.binary.cor)

# ===== inconsistent based on pcls =====
rmv.ic50.pcl.cor = data.frame()
for (pcl in colnames(pcl_temp)) {
  t <- paste0(pcl, ".pearson")
  rmv.ic50.pcl.cor <- rbind(rmv.ic50.pcl.cor, ic50.pcl.cor[which(ic50.pcl.cor[, t] > min), ])
}
pcl.ic50.cc <- rownames(rmv.ic50.pcl.cor)

rmv.auc.pcl.cor = data.frame()
for (pcl in colnames(pcl_temp)) {
  t <- paste0(pcl, ".pearson")
  rmv.auc.pcl.cor <- rbind(rmv.auc.pcl.cor, auc.pcl.cor[which(auc.pcl.cor[, t] > min), ])
}
pcl.auc.cc <- rownames(rmv.auc.pcl.cor)

# ====================================================
#    remove inconsistent cell lines within GDSC 
# ====================================================
# filtering based on binary classifying
gdsc_ic50_cc_binary = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_recomputed", cell.lines = bin.ic50.cc))
gdsc_aac_cc_binary = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed", cell.lines = bin.auc.cc))

# filtering based on PCL classifying
gdsc_ic50_cc_pcl = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_recomputed", cell.lines = pcl.ic50.cc))
gdsc_aac_cc_pcl = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed", cell.lines = pcl.auc.cc))

# ====================================================
#    remove inconsistent cell lines within ctrp 
# ====================================================
# filtering based on binary classifying
ctrp_ic50_cc_binary = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "ic50_recomputed", cell.lines = bin.ic50.cc))
ctrp_aac_cc_binary = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "aac_recomputed", cell.lines = bin.auc.cc))

# filtering based on PCL classifying
ctrp_ic50_cc_pcl = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "ic50_recomputed", cell.lines = pcl.ic50.cc))
ctrp_aac_cc_pcl = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "aac_recomputed", cell.lines = pcl.auc.cc))


# ======================================================================================
#       correlation of published sens measures between ctrp and GDSC (by drug)
# ======================================================================================
orig.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ctrp_ic50))), row.names = rownames(ctrp_ic50))
colnames(orig.ic50.cor) <- coef

for (drug in rownames(orig.ic50.cor)) {
  tryCatch(
    expr = {
      pearson.cor <- cor.test(x = as.numeric(ctrp_ic50[drug, ]), y = as.numeric(gdsc_ic50[drug, ]), method = 'pearson', use = 'pairw')
      orig.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
      orig.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
      
      spearman.cor <- cor.test(x = as.numeric(ctrp_ic50[drug, ]), y = as.numeric(gdsc_ic50[drug, ]), method = 'spearman', use = 'pairw')
      orig.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
      orig.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
    },
    error ={function(e) {
      # message(drug, " does not have enough finite observations")
    }}
  )
}

orig.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ctrp_aac))), row.names = rownames(ctrp_aac))
colnames(orig.auc.cor) <- coef

for (drug in rownames(orig.auc.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = as.numeric(ctrp_aac[drug, ]), y = as.numeric(gdsc_aac[drug, ]), method = 'pearson', use = 'pairw')
  orig.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  orig.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ctrp_aac[drug, ]), y = as.numeric(gdsc_aac[drug, ]), method = 'spearman', use = 'pairw')
  orig.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  orig.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
    },
  error ={function(e) {
    # message(drug, " does not have enough finite observations")
  }}
  )
}

marray::write.xls(orig.ic50.cor, "~/capsule/results/orig.ic50.cor.xls", row.names = TRUE, col.names = TRUE)
marray::write.xls(orig.auc.cor, "~/capsule/results/orig.auc.cor.xls", row.names = TRUE, col.names = TRUE)

# ======================================================================================
#   correlation between ctrp and GDSC after sorting into broad-spectrum/targeted
# ======================================================================================
# ==== correlation after binary sorting ====
bin.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ctrp_ic50_cc_binary))), row.names = rownames(ctrp_ic50_cc_binary))
colnames(bin.ic50.cor) <- coef

for (drug in rownames(bin.ic50.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = as.numeric(ctrp_ic50_cc_binary[drug, ]), y = as.numeric(gdsc_ic50_cc_binary[drug, ]), method = 'pearson', use = 'pairw')
  bin.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
  bin.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ctrp_ic50_cc_binary[drug, ]), y = as.numeric(gdsc_ic50_cc_binary[drug, ]), method = 'spearman', use = 'pairw')
  bin.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
  bin.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
    },
  error ={function(e) {
    # message(drug, " does not have enough finite observations")
  }}
  )
}

bin.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ctrp_aac_cc_binary))), row.names = rownames(ctrp_aac_cc_binary))
colnames(bin.auc.cor) <- coef

for (drug in rownames(bin.auc.cor)) {
    tryCatch(
      expr = {
  pearson.cor <- cor.test(x = as.numeric(ctrp_aac_cc_binary[drug, ]), y = as.numeric(gdsc_aac_cc_binary[drug, ]), method = 'pearson', use = 'pairw')
  bin.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  bin.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ctrp_aac_cc_binary[drug, ]), y = as.numeric(gdsc_aac_cc_binary[drug, ]), method = 'spearman', use = 'pairw')
  bin.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  bin.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
      },
  error ={function(e) {
    # message(drug, " does not have enough finite observations")
  }}
    )
  }

# ======================================================================================
#           Harrell's CI after sorting into broad-spectrum/targeted
# ======================================================================================
binary.ic50.cor <- merge(orig.ic50.cor, bin.ic50.cor, by = 0, all = TRUE)
binary.ic50.cor$pearson.difference <- bin.ic50.cor$pearson - orig.ic50.cor$pearson 
rownames(binary.ic50.cor) <- binary.ic50.cor$Row.names
marray::write.xls(binary.ic50.cor, "~/capsule/results/binary.ic50.cor.xls", row.names = TRUE, col.names = TRUE)

binary.ic50.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = binary.ic50.cor)
saveRDS(binary.ic50.concordance, "~/capsule/results/binary.ic50.concordance.rds")

# AUC
binary.auc.cor <- merge(orig.auc.cor, bin.auc.cor, by = 0, all = TRUE)
binary.auc.cor$pearson.difference <- bin.auc.cor$pearson - orig.auc.cor$pearson
rownames(binary.auc.cor) <- binary.auc.cor$Row.names
marray::write.xls(binary.auc.cor, "~/capsule/results/binary.auc.cor.xls", row.names = TRUE, col.names = TRUE)
binary.auc.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = binary.auc.cor)
saveRDS(binary.auc.concordance, "~/capsule/results/binary.auc.concordance.rds")

# =============================================================================
#       correlation between CTRP and GDSC after sorting into PCLs
# =============================================================================
pcl.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ctrp_ic50_cc_pcl))), row.names = rownames(ctrp_ic50_cc_pcl))
colnames(pcl.ic50.cor) <- coef

for (drug in rownames(pcl.ic50.cor)) {
    tryCatch(
      expr = {
  pearson.cor <- cor.test(x = as.numeric(ctrp_ic50_cc_pcl[drug, ]), y = as.numeric(gdsc_ic50_cc_pcl[drug, ]), method = 'pearson', use = 'pairw')
  pcl.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
  pcl.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ctrp_ic50_cc_pcl[drug, ]), y = as.numeric(gdsc_ic50_cc_pcl[drug, ]), method = 'spearman', use = 'pairw')
  pcl.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
  pcl.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
      },
  error ={function(e) {
    # message(drug, " does not have enough finite observations")
  }}
    )
}

pcl.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ctrp_aac_cc_pcl))), row.names = rownames(ctrp_aac_cc_pcl))
colnames(pcl.auc.cor) <- coef

for (drug in rownames(pcl.auc.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = as.numeric(ctrp_aac_cc_pcl[drug, ]), y = as.numeric(gdsc_aac_cc_pcl[drug, ]), method = 'pearson', use = 'pairw')
  pcl.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  pcl.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ctrp_aac_cc_pcl[drug, ]), y = as.numeric(gdsc_aac_cc_pcl[drug, ]), method = 'spearman', use = 'pairw')
  pcl.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  pcl.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
    },
  error ={function(e) {
    # message(drug, " does not have enough finite observations")
  }}
  )
}

# ====================================================================================
#               Harrell's CI after sorting into PCLs
# ====================================================================================
orig.pcl.ic50.cor <- merge(orig.ic50.cor, pcl.ic50.cor, by = 0, all = TRUE)
orig.pcl.ic50.cor$pearson.difference <- pcl.ic50.cor$pearson - orig.ic50.cor$pearson
rownames(orig.pcl.ic50.cor) <- orig.pcl.ic50.cor$Row.names
marray::write.xls(orig.pcl.ic50.cor, "~/capsule/results/pcl.ic50.cor.xls", row.names = TRUE, col.names = TRUE)
pcl.ic50.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.pcl.ic50.cor)
saveRDS(pcl.ic50.concordance, "~/capsule/results/pcl.ic50.concordance.rds")

orig.pcl.auc.cor <- merge(orig.auc.cor, pcl.auc.cor, by = 0, all = TRUE)
orig.pcl.auc.cor$pearson.difference <- pcl.auc.cor$pearson - orig.auc.cor$pearson
rownames(orig.pcl.auc.cor) <- orig.pcl.auc.cor$Row.names
marray::write.xls(orig.pcl.auc.cor, "~/capsule/results/pcl.auc.cor.xls", row.names = TRUE, col.names = TRUE)
pcl.auc.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.pcl.auc.cor)
saveRDS(pcl.auc.concordance, "~/capsule/results/pcl.auc.concordance.rds")

# ======================================================
# Bar plot of Pearson coefficients for IC50 of all drugs
aa <- as.matrix(t(cbind(orig.ic50.cor$pearson, bin.ic50.cor$pearson, pcl.ic50.cor$pearson)))
aa[!is.na(aa) & aa < 0] <- 0
names(aa) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/all_ic50_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
ab <- barplot(aa, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(aa), v=0.9), each=2), ylab="Pearson's correlation coefficient r", ylim = c(0, 1), density=c(100,10,0) , angle=c(0,45,90), main="Pearson's correlation coefficient of IC50 of original values and post filtering methods", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ", "pcl"), density=c(100,10,0), angle=c(0,45,90), bty="n", cex=0.75)
text(x=apply(ab, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(aa)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of difference in Pearson coefficients of IC50 after both filtering methods
hh <- as.matrix(t(cbind(binary.ic50.cor$pearson.difference, orig.pcl.ic50.cor$pearson.difference)))
names(hh) <- binary.ic50.cor$Row.names

pdf("~/capsule/results/plots/all_pearson_change_ic50_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
hb <- barplot(hh, beside = TRUE, space = c(0.1, 2), col=rep(rainbow(length(hh), v=0.9), each=2), ylab="Difference in Pearson's correlation coefficient", ylim=c(-1.06, 1.02), density=c(100, 10), angle=c(0, 45), main="Difference in Pearson's correlation coefficient of IC50 of original values and post filtering methods", font.main = 1)
legend("topleft", legend=c("binary", "PCL"), density=c(100, 10), angle=c(0, 45), bty = "n", cex=0.75)
text(x=apply(hb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(hh)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of difference in Pearson coefficients of IC50 after binary filtering methods
hh <- as.matrix(t(binary.ic50.cor$pearson.difference))
names(hh) <- binary.ic50.cor$Row.names

pdf("~/capsule/results/plots/bin_pearson_change_ic50_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
hb <- barplot(hh, beside = TRUE, space = c(0.1, 2), col=rep(rainbow(length(hh), v=0.9), each=2), ylab="Difference in Pearson's correlation coefficient", ylim=c(-1.06, 1.02), density=c(100), angle=c(0), main="Difference in Pearson's correlation coefficient of IC50 of original values and post binary drug filtering methods", font.main = 1)
legend("topleft", legend=c("binary"), density=c(100), angle=c(0), bty = "n", cex=0.75)
text(x=apply(hb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(hh)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of difference in Pearson coefficients of IC50 after PCL filtering methods
hh <- as.matrix(t(orig.pcl.ic50.cor$pearson.difference))
names(hh) <- orig.pcl.ic50.cor$Row.names

pdf("~/capsule/results/plots/pcl_pearson_change_ic50_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
hb <- barplot(hh, beside = TRUE, space = c(0.1, 2), col=rep(rainbow(length(hh), v=0.9), each=2), ylab="Difference in Pearson's correlation coefficient", ylim=c(-1.06, 1.02), density=c(100), angle=c(0), main="Difference in Pearson's correlation coefficient of IC50 of original values and post PCL drug filtering methods", font.main = 1)
legend("topleft", legend=c("PCL"), density=c(100), angle=c(0), bty = "n", cex=0.75)
text(x=apply(hb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(hh)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for IC50 (excluding drugs with pearson below 0)
dr <- paste(rownames(orig.ic50.cor[which(orig.ic50.cor$pearson > 0), ]))
dr <- append(dr, rownames(bin.ic50.cor[which(bin.ic50.cor$pearson > 0), ]))
dr <- append(dr, rownames(pcl.ic50.cor[which(pcl.ic50.cor$pearson > 0), ]))
dr <- unique(dr)

ss <- as.matrix(t(cbind(orig.ic50.cor[dr, "pearson"], bin.ic50.cor[dr, "pearson"], pcl.ic50.cor[dr, "pearson"])))
ss[!is.na(ss) & ss < 0] <- 0
names(ss) <- rownames(orig.ic50.cor[dr, ])

pdf("~/capsule/results/plots/positive_ic50_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
mp <- barplot(ss, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(ss), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,10,0) , angle=c(0,45,90), main="Pearson's correlation coefficient of IC50 of original values \nand post filtering methods", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ", "pcl"), density=c(100,10,0), angle=c(0,45,90), bty="n", cex=0.75)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(ss)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for IC50 for all drugs after binary sorting and PCL sorting
kk <- as.matrix(t(cbind(bin.ic50.cor$pearson, pcl.ic50.cor$pearson)))
kk[!is.na(kk) & kk < 0] <- 0
names(kk) <- rownames(bin.ic50.cor)

pdf("~/capsule/results/plots/bin_pcl_ic50_barplot.pdf", height = 8, width = 15)
par(mar=c(8,5,5,5))
kp <- barplot(kk, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(kk), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,0) , angle=c(0,45), main="Pearson's correlation coefficient of IC50 of values post filtering methods", font.main = 1)
legend("topright", legend=c("br-sp/targ", "pcl"), density=c(100,0), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(kp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(kk)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()


# Bar plot of Pearson coefficients for IC50 values for the 15 drugs with highest and lowest correlations after binary filtering
t <- na.omit(binary.ic50.cor[with(binary.ic50.cor,order(binary.ic50.cor$pearson.difference)), ])
highest <- rownames(tail(t, n=15))
lowest <- rownames(head(t, n=15))

cc <- as.matrix(t(cbind(binary.ic50.cor[highest, "pearson.x"], binary.ic50.cor[highest, "pearson.y"])))
names(cc) <- rownames(binary.ic50.cor[highest, ])

pdf("~/capsule/results/plots/top15_bin_ic50_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
cp <- barplot(cc, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(cc), v=0.9), each=2), ylim = c(-0.3, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of IC50 of original values \nand after binary drug classification for the 15 drugs with largest changes in correlation", font.main = 1)
legend("topleft", legend=c("orig", "binary"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(cp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(cc)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

dd <- as.matrix(t(cbind(binary.ic50.cor[lowest, "pearson.x"], binary.ic50.cor[lowest, "pearson.y"])))
names(dd) <- rownames(orig.ic50.cor[lowest, ])

pdf("~/capsule/results/plots/lowest15_bin_ic50_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
dp <- barplot(dd, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(cc), v=0.9), each=2), ylim = c(-0.2, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of IC50 of original values \nand after binary drug classification for the 15 drugs with lowest changes in correlation", font.main = 1)
legend("topright", legend=c("orig", "binary"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(dp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(dd)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for IC50 values for the 15 drugs with highest and lowest correlations after PCL filtering
w <- na.omit(orig.pcl.ic50.cor[with(orig.pcl.ic50.cor,order(orig.pcl.ic50.cor$pearson.difference)), ])
highest <- rownames(tail(w, n=15))
lowest <- rownames(head(w, n=15))

ee <- as.matrix(t(cbind(orig.pcl.ic50.cor[highest, "pearson.x"], orig.pcl.ic50.cor[highest, "pearson.y"])))
names(ee) <- rownames(orig.pcl.ic50.cor[highest, ])

pdf("~/capsule/results/plots/top15_pcl_ic50_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
ep <- barplot(ee, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(ee), v=0.9), each=2), ylim = c(-0.2, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of IC50 of original values \nand after PCL drug classification for the 15 drugs with largest changes in correlation", font.main = 1)
legend("topright", legend=c("orig", "PCL"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(ep, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(ee)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

ff <- as.matrix(t(cbind(orig.pcl.ic50.cor[lowest, "pearson.x"], orig.pcl.ic50.cor[lowest, "pearson.y"])))
names(ff) <- rownames(orig.pcl.ic50.cor[lowest, ])

pdf("~/capsule/results/plots/lowest15_pcl_ic50_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
fp <- barplot(ff, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(ff), v=0.9), each=2), ylim = c(-0.6, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of IC50 of original values \nand after PCL drug classification for the 15 drugs with lowest changes in correlation", font.main = 1)
legend("bottomright", legend=c("orig", "PCL"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(fp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(ff)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for AUC of all drugs
cc <- as.matrix(t(cbind(orig.auc.cor$pearson, orig.auc.cor$pearson, pcl.auc.cor$pearson)))
cc[!is.na(cc) & cc < 0] <- 0
names(cc) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/all_auc_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
bp <- barplot(cc, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(cc), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,10,0) , angle=c(0,45,90), main="Pearson's correlation coefficient of AUC of original values \nand post filtering methods", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ", "pcl"), density=c(100,10,0), angle=c(0,45,90), bty="n", cex=0.75)
text(x=apply(bp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(cc)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for AUC (excluding drugs with pearson below 0)
tr <- paste(rownames(orig.auc.cor[which(orig.auc.cor$pearson > 0), ]))
tr <- append(tr, rownames(bin.auc.cor[which(bin.auc.cor$pearson > 0), ]))
tr <- append(tr, rownames(pcl.auc.cor[which(pcl.auc.cor$pearson > 0), ]))
tr <- unique(tr)

tg <- as.matrix(t(cbind(orig.auc.cor[tr, "pearson"], bin.auc.cor[tr, "pearson"], pcl.auc.cor[tr, "pearson"])))
tg[!is.na(tg) & tg < 0] <- 0
names(tg) <- rownames(orig.auc.cor)

pdf("~/capsule/results/plots/positive_auc_bar_plot.pdf", height = 8, width = 15)
par(mar=c(8,5,5,5))
tgb <- barplot(tg, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(tg), v=0.9), each=2), ylim = c(0, 0.9), ylab="Pearson's correlation coefficient r", density=c(100,10, 0) , angle=c(0,45, 90), main="Pearson's correlation coefficient of AUC of original values \nand post filtering methods", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ", "pcl"), density=c(100, 10, 0),angle=c(0,45,90), bty="n", cex=0.75)
text(x=apply(tgb, 2, mean) + 1.5, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(tg)), srt=45, xpd=NA, font=0.75, cex = 0.75)
dev.off()

# Bar plot of difference in Pearson coefficients of AUC after both filtering methods
hh <- as.matrix(t(cbind(binary.auc.cor$pearson.difference, orig.pcl.auc.cor$pearson.difference)))
names(hh) <- binary.auc.cor$Row.names

pdf("~/capsule/results/plots/all_pearson_change_auc_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
hb <- barplot(hh, beside = TRUE, space = c(0.1, 2), col=rep(rainbow(length(hh), v=0.9), each=2), ylab="Difference in Pearson's correlation coefficient", ylim=c(-0.7, 0.5), density=c(100, 10), angle=c(0, 45), main="Difference in Pearson's correlation coefficient of AUC of original values and post filtering methods", font.main = 1)
legend("topleft", legend=c("binary", "PCL"), density=c(100, 10), angle=c(0, 45), bty = "n", cex=0.75)
text(x=apply(hb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(hh)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of difference in Pearson coefficients of AUC after binary filtering methods
hh <- as.matrix(t(binary.auc.cor$pearson.difference))
names(hh) <- binary.auc.cor$Row.names

pdf("~/capsule/results/plots/bin_pearson_change_auc_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
hb <- barplot(hh, beside = TRUE, space = c(0.1, 2), col=rep(rainbow(length(hh), v=0.9), each=2), ylab="Difference in Pearson's correlation coefficient", ylim=c(-0.5, 0.2), density=c(100), angle=c(0), main="Difference in Pearson's correlation coefficient of AUC of original values and post binary drug filtering method", font.main = 1)
legend("topleft", legend=c("binary"), density=c(100), angle=c(0), bty = "n", cex=0.75)
text(x=apply(hb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(hh)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of difference in Pearson coefficients of AUC after PCL filtering methods
hh <- as.matrix(t(orig.pcl.auc.cor$pearson.difference))
names(hh) <- orig.pcl.auc.cor$Row.names

pdf("~/capsule/results/plots/pcl_pearson_change_auc_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
hb <- barplot(hh, beside = TRUE, space = c(0.1, 2), col=rep(rainbow(length(hh), v=0.9), each=2), ylab="Difference in Pearson's correlation coefficient", ylim=c(-0.7, 0.5), density=c(100), angle=c(0), main="Difference in Pearson's correlation coefficient of AUC of original values and post binary drug filtering method", font.main = 1)
legend("topleft", legend=c("PCL"), density=c(100), angle=c(0), bty = "n", cex=0.75)
text(x=apply(hb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(hh)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()


# Bar plot of Pearson coefficients for AUC for all drugs after binary sorting and PCL sorting
ii <- as.matrix(t(cbind(bin.auc.cor$pearson, pcl.auc.cor$pearson)))
ii[!is.na(ii) & ii < 0] <- 0
names(ii) <- rownames(bin.auc.cor)

pdf("~/capsule/results/plots/bin_pcl_auc_barplot.pdf", height = 8, width = 15)
par(mar=c(8,5,5,5))
ip <- barplot(ii, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(ii), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,0) , angle=c(0,45), main="Pearson's correlation coefficient of AUC of values post filtering methods", font.main = 1)
legend("topright", legend=c("br-sp/targ", "pcl"), density=c(100,0), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(ip, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(ii)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for AUC values for the 15 drugs with highest and lowest correlations after binary filtering
q <- na.omit(binary.auc.cor[with(binary.auc.cor,order(binary.auc.cor$pearson.difference)), ])
highest <- rownames(tail(q, n=15))
lowest <- rownames(head(q, n=15))

qq <- as.matrix(t(cbind(binary.auc.cor[highest, "pearson.x"], binary.auc.cor[highest, "pearson.y"])))
names(qq) <- rownames(binary.auc.cor[highest, ])

pdf("~/capsule/results/plots/top15_bin_auc_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
qp <- barplot(qq, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(qq), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of AUC of original values \nand after binary drug classification for the 15 drugs with largest changes in correlation", font.main = 1)
legend("topright", legend=c("orig", "binary"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(qp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(qq)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

ss <- as.matrix(t(cbind(binary.auc.cor[lowest, "pearson.x"], binary.auc.cor[lowest, "pearson.y"])))
names(ss) <- rownames(binary.auc.cor[lowest, ])

pdf("~/capsule/results/plots/lowest15_bin_auc_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
sp <- barplot(ss, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(ss), v=0.9), each=2), ylim = c(-1, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of AUC of original values \nand after binary drug classification for the 15 drugs with lowest changes in correlation", font.main = 1)
legend("topleft", legend=c("orig", "binary"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(sp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(ss)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for AUC values for the 15 drugs with highest and lowest correlations after PCL filtering
y <- na.omit(pcl.auc.cor[with(pcl.auc.cor,order(pcl.auc.cor$pearson)), ])
highest <- rownames(tail(y, n=15))
lowest <- rownames(head(y, n=15))

yy <- as.matrix(t(cbind(orig.auc.cor[highest, "pearson"], pcl.auc.cor[highest, "pearson"])))
names(yy) <- rownames(orig.auc.cor[highest, ])

pdf("~/capsule/results/plots/top15_pcl_auc_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
yp <- barplot(yy, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(yy), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of AUC of original values \nand after PCL drug classification for the 15 drugs with highest correlations", font.main = 1)
legend("topleft", legend=c("orig", "PCL"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(yp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(yy)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

zz <- as.matrix(t(cbind(orig.auc.cor[lowest, "pearson"], pcl.auc.cor[lowest, "pearson"])))
names(zz) <- rownames(orig.auc.cor[lowest, ])

pdf("~/capsule/results/plots/lowest15_pcl_auc_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
zp <- barplot(zz, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(zz), v=0.9), each=2), ylim = c(-0.4, 0.8), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of AUC of original values \nand after PCL drug classification for the 15 drugs with lowest correlations", font.main = 1)
legend("topleft", legend=c("orig", "PCL"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(zp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(zz)), srt=45, xpd=NA, font=1, cex = 0.75)
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
