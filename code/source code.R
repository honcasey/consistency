library(PharmacoGx); library(dplyr); library(reshape2); library(survival); library(VennDiagram)
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
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CTRP@cell), area2=nrow(GDSC@cell), cross.area=nrow(ctrp@cell), col = "black", cex=1.5, cat.cex=1, cat.col = c("black", "black"))
pdf("~/capsule/results/plots/venn/cell_intersection.pdf", height=4, width=4)
grid::grid.draw(venn.plot)
dev.off()
rm(venn.plot)

# venn diagram of common drugs
venn.plot2 <- VennDiagram::draw.pairwise.venn(area1=nrow(CTRP@drug), area2=nrow(GDSC@drug), cross.area=nrow(ctrp@drug), col = "black" ,cex=1.5, cat.cex=1, cat.col = c("black", "black"))
pdf("~/capsule/results/plots/venn/drug_intersection.pdf", height=4, width=4)
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
# druglist <- as.data.frame(read.csv("PSets/drug_info.csv")) 
druglist <- as.data.frame(read_excel("drug_info.xls"))
druglist[] <- lapply(druglist, as.character) # change class of dataframe columns from vector to character
rownames(druglist) <- rownames(ctrp@drug)
# rownames(druglist) <- druglist$drug_name 

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
# ctrp_ic50 = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "ic50_published"))
ctrp_ic50 = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "ic50_recomputed"))

# ctrp_auc = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "auc_published"))
ctrp_aac = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "aac_recomputed"))

# gdsc_ic50 = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published"))
gdsc_ic50 = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_recomputed"))

# gdsc_auc = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published"))
gdsc_aac = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed"))

# ====================================================
# correlation of published IC50 between ctrp and GDSC 
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
rmv.ic50.binary.cor <- rbind(ic50.brsp.cor[which(ic50.brsp.cor$spearman > min), ], ic50.targ.cor[which(ic50.targ.cor$spearman > min), ])
bin.ic50.cc <- rownames(rmv.ic50.binary.cor)
rmv.auc.binary.cor <- rbind(auc.brsp.cor[which(auc.brsp.cor$spearman > min), ], auc.targ.cor[which(auc.targ.cor$spearman > min), ])
bin.auc.cc <- rownames(rmv.auc.binary.cor)

# ===== inconsistent based on pcls =====
rmv.ic50.pcl.cor = data.frame()
for (pcl in colnames(pcl_temp)) {
  t <- paste0(pcl, ".spearman")
  rmv.ic50.pcl.cor <- rbind(rmv.ic50.pcl.cor, ic50.pcl.cor[which(ic50.pcl.cor[, t] > min), ])
}
pcl.ic50.cc <- rownames(rmv.ic50.pcl.cor)

rmv.auc.pcl.cor = data.frame()
for (pcl in colnames(pcl_temp)) {
  t <- paste0(pcl, ".spearman")
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
marray::write.xls(binary.ic50.cor, "~/capsule/results/binary.ic50.cor.xls", row.names = TRUE, col.names = TRUE)
binary.ic50.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = binary.ic50.cor)
saveRDS(binary.ic50.concordance, "~/capsule/results/binary.ic50.concordance.rds")

# AUC
binary.auc.cor <- merge(orig.auc.cor, bin.auc.cor, by = 0, all = TRUE)
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
marray::write.xls(orig.pcl.ic50.cor, "~/capsule/results/pcl.ic50.cor.xls", row.names = TRUE, col.names = TRUE)
pcl.ic50.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.pcl.ic50.cor)
saveRDS(pcl.ic50.concordance, "~/capsule/results/pcl.ic50.concordance.rds")

orig.pcl.auc.cor <- merge(orig.auc.cor, pcl.auc.cor, by = 0, all = TRUE)
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

# 6) Bar plot representing spearman's ranks for AUC drug measure comparing original with pcls
mek <- as.matrix(t(cbind(orig.auc.cor$spearman, pcl.auc.cor$spearman)))
mek[!is.na(mek) & mek < 0] <- 0
names(mek) <- rownames(orig.auc.cor)

pdf("~/capsule/results/plots/pcl_auc_cor_bar_plot.pdf")
mb <- barplot(mek, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(mek), v=0.9), each=2), ylab=expression("Spearman's rank correlation coefficient r"[s]), density=c(100,10) , angle=c(0,45), main="Spearman's rank correlation coefficient for AUC values \nafter filtering cell lines by PCL drug classification", font.main = 1)
legend("topright", legend=c("orig", "PCL"), fill=c("black", "black"), density=c(100, 10), bty="n", cex=1)
text(x=apply(mb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.01), pos=1, labels=toupper(names(mek)), srt=50, xpd=NA, font=1, cex = 0.75)
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
