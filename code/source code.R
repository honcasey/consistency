library(PharmacoGx); library(dplyr); library(reshape2); library(survival); library(VennDiagram); library(readxl)
CTRP <- downloadPSet("CTRPv2_2015")
GDSC <- downloadPSet("GDSC_2020(v2-8.2)")

# ==============================================================
# intersecting CTRP and GDSC for common drugs and cell lines
# ==============================================================
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
# marray::write.xls(common_drugs, "~/capsule/results/common_drugs.xls", row.names = TRUE, col.names = TRUE)

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
# correlation of published AAC between ctrp and GDSC (by cell line)
# ====================================================
aac.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ctrp_aac))), row.names = colnames(ctrp_aac))
colnames(aac.cor) <- coef

for (cell.line in rownames(aac.cor)) {
  tryCatch(
    expr = {
      pearson.cor <- cor.test(x = ctrp_aac[, cell.line], y = gdsc_aac[, cell.line], method = 'pearson', use = 'pairw')
      aac.cor[cell.line, "pearson"] <- pearson.cor$estimate
      aac.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
      
      spearman.cor <- cor.test(x = ctrp_aac[, cell.line], y = gdsc_aac[, cell.line], method = 'spearman', use = 'pairwise.complete.obs')
      aac.cor[cell.line, "spearman"] <- spearman.cor$estimate
      aac.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
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

# cor of aac for each cell line + drugs that are broad-spectrum
ctrp_aac_brsp <- ctrp_aac[broad_spectrum_drugs, ]
gdsc_aac_brsp <- gdsc_aac[broad_spectrum_drugs, ] 

aac.brsp.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ctrp_ic50))), row.names = colnames(ctrp_ic50))
colnames(aac.brsp.cor) <- coef

for (cell_line in rownames(aac.brsp.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = ctrp_aac_brsp[, cell.line], y = gdsc_aac_brsp[, cell.line], method = 'pearson', use = 'pairw')
  aac.brsp.cor[cell.line, "pearson"] <- pearson.cor$estimate
  aac.brsp.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ctrp_aac_brsp[, cell.line], y = gdsc_aac_brsp[, cell.line], method = 'spearman', use = 'pairwise.complete.obs')
  aac.brsp.cor[cell.line, "spearman"] <- spearman.cor$estimate
  aac.brsp.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
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

# cor of aac for each cell line + drugs that are targeted
ctrp_aac_targ <- ctrp_aac[targeted_drugs, ]
gdsc_aac_targ <- gdsc_aac[targeted_drugs, ] 

aac.targ.cor = data.frame(matrix(ncol = length(coef), nrow = length(colnames(ctrp_aac))), row.names = colnames(ctrp_aac))
colnames(aac.targ.cor) <- coef

for (cell.line in rownames(aac.targ.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = ctrp_aac_targ[, cell.line], y = gdsc_aac_targ[, cell.line], method = 'pearson', use = 'pairw')
  aac.targ.cor[cell.line, "pearson"] <- pearson.cor$estimate
  aac.targ.cor[cell.line, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = ctrp_aac_targ[, cell.line], y = gdsc_aac_targ[, cell.line], method = 'spearman', use = 'pairw')
  aac.targ.cor[cell.line, "spearman"] <- spearman.cor$estimate
  aac.targ.cor[cell.line, "spearman p-value"] <- spearman.cor$p.value
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

# correlation of aac in drugs by pcl for each cell line 
aac.pcl.cor = data.frame(matrix(NA, nrow = length(colnames(ctrp_aac))), row.names = colnames(ctrp_aac))

for (pcl in colnames(pcl_temp)) {
  aac.pcl.cor[, paste0(pcl, ".pearson")] <- NA
    aac.pcl.cor[, paste0(pcl, ".pearson p-value")] <- NA
  aac.pcl.cor[, paste0(pcl, ".spearman")] <- NA
    aac.pcl.cor[, paste0(pcl, ".spearman p-value")] <- NA
}
aac.pcl.cor <- aac.pcl.cor[, -1]

for (pcl in colnames(pcl_temp)) {
  ctrp_aac_temp <- ctrp_aac[pcl_list[[pcl]], ] # subset by the drugs in the PCL
  gdsc_aac_temp <- gdsc_aac[pcl_list[[pcl]], ]
  for (cell.line in rownames(aac.pcl.cor)) {
    tryCatch(
      expr = {
        pearson.cor <- cor.test(x = ctrp_aac_temp[, cell.line], y = gdsc_aac_temp[, cell.line], method = 'pearson', use = 'pairw')
        aac.pcl.cor[cell.line, paste0(pcl, ".pearson")] <- pearson.cor$estimate
        aac.pcl.cor[cell.line, paste0(pcl, ".pearson p-value")] <- pearson.cor$p.value
  
        spearman.cor <- cor.test(x = ctrp_aac_temp[, cell.line], y = gdsc_aac_temp[, cell.line], method = 'spearman', use = 'pairw')
        aac.pcl.cor[cell.line, paste0(pcl, ".spearman")] <- spearman.cor$estimate
        aac.pcl.cor[cell.line, paste0(pcl, ".spearman p-value")] <- spearman.cor$p.value
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
rmv.aac.binary.cor <- union(aac.brsp.cor[which(aac.brsp.cor$pearson > min), ], aac.targ.cor[which(aac.targ.cor$pearson > min), ])
bin.aac.cc <- rownames(rmv.aac.binary.cor)

# ===== inconsistent based on pcls =====
rmv.ic50.pcl.cor = data.frame()
for (pcl in colnames(pcl_temp)) {
  t <- paste0(pcl, ".pearson")
  rmv.ic50.pcl.cor <- rbind(rmv.ic50.pcl.cor, ic50.pcl.cor[which(ic50.pcl.cor[, t] > min), ])
}
pcl.ic50.cc <- rownames(rmv.ic50.pcl.cor)

rmv.aac.pcl.cor = data.frame()
for (pcl in colnames(pcl_temp)) {
  t <- paste0(pcl, ".pearson")
  rmv.aac.pcl.cor <- rbind(rmv.aac.pcl.cor, aac.pcl.cor[which(aac.pcl.cor[, t] > min), ])
}
pcl.aac.cc <- rownames(rmv.aac.pcl.cor)

# ====================================================
#    remove inconsistent cell lines within GDSC 
# ====================================================
# filtering based on binary classifying
gdsc_ic50_cc_binary = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_recomputed", cell.lines = bin.ic50.cc))
gdsc_aac_cc_binary = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed", cell.lines = bin.aac.cc))

# filtering based on PCL classifying
gdsc_ic50_cc_pcl = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_recomputed", cell.lines = pcl.ic50.cc))
gdsc_aac_cc_pcl = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed", cell.lines = pcl.aac.cc))

# ====================================================
#    remove inconsistent cell lines within ctrp 
# ====================================================
# filtering based on binary classifying
ctrp_ic50_cc_binary = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "ic50_recomputed", cell.lines = bin.ic50.cc))
ctrp_aac_cc_binary = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "aac_recomputed", cell.lines = bin.aac.cc))

# filtering based on PCL classifying
ctrp_ic50_cc_pcl = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "ic50_recomputed", cell.lines = pcl.ic50.cc))
ctrp_aac_cc_pcl = as.data.frame(summarizeSensitivityProfiles(ctrp, sensitivity.measure = "aac_recomputed", cell.lines = pcl.aac.cc))


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

orig.aac.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ctrp_aac))), row.names = rownames(ctrp_aac))
colnames(orig.aac.cor) <- coef

for (drug in rownames(orig.aac.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = as.numeric(ctrp_aac[drug, ]), y = as.numeric(gdsc_aac[drug, ]), method = 'pearson', use = 'pairw')
  orig.aac.cor[drug, "pearson"] <- pearson.cor$estimate
  orig.aac.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ctrp_aac[drug, ]), y = as.numeric(gdsc_aac[drug, ]), method = 'spearman', use = 'pairw')
  orig.aac.cor[drug, "spearman"] <- spearman.cor$estimate
  orig.aac.cor[drug, "spearman p-value"] <- spearman.cor$p.value
    },
  error ={function(e) {
    # message(drug, " does not have enough finite observations")
  }}
  )
}

marray::write.xls(orig.ic50.cor, "~/capsule/results/orig.ic50.cor.xls", row.names = TRUE, col.names = TRUE)
marray::write.xls(orig.aac.cor, "~/capsule/results/orig.aac.cor.xls", row.names = TRUE, col.names = TRUE)

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

bin.aac.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ctrp_aac_cc_binary))), row.names = rownames(ctrp_aac_cc_binary))
colnames(bin.aac.cor) <- coef

for (drug in rownames(bin.aac.cor)) {
    tryCatch(
      expr = {
  pearson.cor <- cor.test(x = as.numeric(ctrp_aac_cc_binary[drug, ]), y = as.numeric(gdsc_aac_cc_binary[drug, ]), method = 'pearson', use = 'pairw')
  bin.aac.cor[drug, "pearson"] <- pearson.cor$estimate
  bin.aac.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ctrp_aac_cc_binary[drug, ]), y = as.numeric(gdsc_aac_cc_binary[drug, ]), method = 'spearman', use = 'pairw')
  bin.aac.cor[drug, "spearman"] <- spearman.cor$estimate
  bin.aac.cor[drug, "spearman p-value"] <- spearman.cor$p.value
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

# aac
binary.aac.cor <- merge(orig.aac.cor, bin.aac.cor, by = 0, all = TRUE)
binary.aac.cor$pearson.difference <- bin.aac.cor$pearson - orig.aac.cor$pearson
rownames(binary.aac.cor) <- binary.aac.cor$Row.names
marray::write.xls(binary.aac.cor, "~/capsule/results/binary.aac.cor.xls", row.names = TRUE, col.names = TRUE)
binary.aac.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = binary.aac.cor)
saveRDS(binary.aac.concordance, "~/capsule/results/binary.aac.concordance.rds")

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

pcl.aac.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ctrp_aac_cc_pcl))), row.names = rownames(ctrp_aac_cc_pcl))
colnames(pcl.aac.cor) <- coef

for (drug in rownames(pcl.aac.cor)) {
  tryCatch(
    expr = {
  pearson.cor <- cor.test(x = as.numeric(ctrp_aac_cc_pcl[drug, ]), y = as.numeric(gdsc_aac_cc_pcl[drug, ]), method = 'pearson', use = 'pairw')
  pcl.aac.cor[drug, "pearson"] <- pearson.cor$estimate
  pcl.aac.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ctrp_aac_cc_pcl[drug, ]), y = as.numeric(gdsc_aac_cc_pcl[drug, ]), method = 'spearman', use = 'pairw')
  pcl.aac.cor[drug, "spearman"] <- spearman.cor$estimate
  pcl.aac.cor[drug, "spearman p-value"] <- spearman.cor$p.value
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

orig.pcl.aac.cor <- merge(orig.aac.cor, pcl.aac.cor, by = 0, all = TRUE)
orig.pcl.aac.cor$pearson.difference <- pcl.aac.cor$pearson - orig.aac.cor$pearson
rownames(orig.pcl.aac.cor) <- orig.pcl.aac.cor$Row.names
marray::write.xls(orig.pcl.aac.cor, "~/capsule/results/pcl.aac.cor.xls", row.names = TRUE, col.names = TRUE)
pcl.aac.concordance <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.pcl.aac.cor)
saveRDS(pcl.aac.concordance, "~/capsule/results/pcl.aac.concordance.rds")

# ======================================================
# bar plot of original pearson coefficients for ic50
rr <- as.matrix(t(cbind(orig.ic50.cor$pearson)))
rr[!is.na(rr) & rr < 0] <- 0
names(rr) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/orig_ic50_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
rb <- barplot(rr, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(aa), v=0.9), each=2), ylab="Pearson's correlation coefficient r", ylim = c(0, 1), density=c(100) , angle=c(0), main="Pearson's correlation coefficient of IC50 of original values", font.main = 1)
legend("topright", legend=c("orig"), density=c(100), angle=c(0), bty="n", cex=0.75)
text(x=apply(rb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(rr)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# bar plot of original pearson coefficients for aac
rr <- as.matrix(t(cbind(orig.aac.cor$pearson)))
rr[!is.na(rr) & rr < 0] <- 0
names(rr) <- rownames(orig.aac.cor)

pdf("~/capsule/results/plots/orig_aac_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
rb <- barplot(rr, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(aa), v=0.9), each=2), ylab="Pearson's correlation coefficient r", ylim = c(0, 1), density=c(100) , angle=c(0), main="Pearson's correlation coefficient of AAC of original values", font.main = 1)
legend("topright", legend=c("orig"), density=c(100), angle=c(0), bty="n", cex=0.75)
text(x=apply(rb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(rr)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

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

# Bar plot of Pearson coefficients for aac of all drugs
cc <- as.matrix(t(cbind(orig.aac.cor$pearson, orig.aac.cor$pearson, pcl.aac.cor$pearson)))
cc[!is.na(cc) & cc < 0] <- 0
names(cc) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/all_aac_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
bp <- barplot(cc, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(cc), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,10,0) , angle=c(0,45,90), main="Pearson's correlation coefficient of aac of original values \nand post filtering methods", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ", "pcl"), density=c(100,10,0), angle=c(0,45,90), bty="n", cex=0.75)
text(x=apply(bp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(cc)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for aac (excluding drugs with pearson below 0)
tr <- paste(rownames(orig.aac.cor[which(orig.aac.cor$pearson > 0), ]))
tr <- append(tr, rownames(bin.aac.cor[which(bin.aac.cor$pearson > 0), ]))
tr <- append(tr, rownames(pcl.aac.cor[which(pcl.aac.cor$pearson > 0), ]))
tr <- unique(tr)

tg <- as.matrix(t(cbind(orig.aac.cor[tr, "pearson"], bin.aac.cor[tr, "pearson"], pcl.aac.cor[tr, "pearson"])))
tg[!is.na(tg) & tg < 0] <- 0
names(tg) <- rownames(orig.aac.cor)

pdf("~/capsule/results/plots/positive_aac_bar_plot.pdf", height = 8, width = 15)
par(mar=c(8,5,5,5))
tgb <- barplot(tg, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(tg), v=0.9), each=2), ylim = c(0, 0.9), ylab="Pearson's correlation coefficient r", density=c(100,10, 0) , angle=c(0,45, 90), main="Pearson's correlation coefficient of aac of original values \nand post filtering methods", font.main = 1)
legend("topright", legend=c("orig", "br-sp/targ", "pcl"), density=c(100, 10, 0),angle=c(0,45,90), bty="n", cex=0.75)
text(x=apply(tgb, 2, mean) + 1.5, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(tg)), srt=45, xpd=NA, font=0.75, cex = 0.75)
dev.off()

# Bar plot of difference in Pearson coefficients of aac after both filtering methods
hh <- as.matrix(t(cbind(binary.aac.cor$pearson.difference, orig.pcl.aac.cor$pearson.difference)))
names(hh) <- binary.aac.cor$Row.names

pdf("~/capsule/results/plots/all_pearson_change_aac_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
hb <- barplot(hh, beside = TRUE, space = c(0.1, 2), col=rep(rainbow(length(hh), v=0.9), each=2), ylab="Difference in Pearson's correlation coefficient", ylim=c(-0.7, 0.5), density=c(100, 10), angle=c(0, 45), main="Difference in Pearson's correlation coefficient of aac of original values and post filtering methods", font.main = 1)
legend("topleft", legend=c("binary", "PCL"), density=c(100, 10), angle=c(0, 45), bty = "n", cex=0.75)
text(x=apply(hb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(hh)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of difference in Pearson coefficients of aac after binary filtering methods
hh <- as.matrix(t(binary.aac.cor$pearson.difference))
names(hh) <- binary.aac.cor$Row.names

pdf("~/capsule/results/plots/bin_pearson_change_aac_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
hb <- barplot(hh, beside = TRUE, space = c(0.1, 2), col=rep(rainbow(length(hh), v=0.9), each=2), ylab="Difference in Pearson's correlation coefficient", ylim=c(-0.5, 0.2), density=c(100), angle=c(0), main="Difference in Pearson's correlation coefficient of aac of original values and post binary drug filtering method", font.main = 1)
legend("topleft", legend=c("binary"), density=c(100), angle=c(0), bty = "n", cex=0.75)
text(x=apply(hb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(hh)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of difference in Pearson coefficients of aac after PCL filtering methods
hh <- as.matrix(t(orig.pcl.aac.cor$pearson.difference))
names(hh) <- orig.pcl.aac.cor$Row.names

pdf("~/capsule/results/plots/pcl_pearson_change_aac_bar_plot.pdf", height = 9, width = 15)
par(mar=c(8,5,5,5))
hb <- barplot(hh, beside = TRUE, space = c(0.1, 2), col=rep(rainbow(length(hh), v=0.9), each=2), ylab="Difference in Pearson's correlation coefficient", ylim=c(-0.7, 0.5), density=c(100), angle=c(0), main="Difference in Pearson's correlation coefficient of AAC of original values and post PCL drug filtering method", font.main = 1)
legend("topleft", legend=c("PCL"), density=c(100), angle=c(0), bty = "n", cex=0.75)
text(x=apply(hb, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05) + 0.03, pos=2, labels=toupper(names(hh)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for aac for all drugs after binary sorting and PCL sorting
ii <- as.matrix(t(cbind(bin.aac.cor$pearson, pcl.aac.cor$pearson)))
ii[!is.na(ii) & ii < 0] <- 0
names(ii) <- rownames(bin.aac.cor)

pdf("~/capsule/results/plots/bin_pcl_aac_barplot.pdf", height = 8, width = 15)
par(mar=c(8,5,5,5))
ip <- barplot(ii, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(ii), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,0) , angle=c(0,45), main="Pearson's correlation coefficient of AAC of values post filtering methods", font.main = 1)
legend("topright", legend=c("br-sp/targ", "pcl"), density=c(100,0), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(ip, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(ii)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for aac values for the 15 drugs with highest and lowest correlations after binary filtering
q <- na.omit(binary.aac.cor[with(binary.aac.cor,order(binary.aac.cor$pearson.difference)), ])
highest <- rownames(tail(q, n=15))
lowest <- rownames(head(q, n=15))

qq <- as.matrix(t(cbind(binary.aac.cor[highest, "pearson.x"], binary.aac.cor[highest, "pearson.y"])))
names(qq) <- rownames(binary.aac.cor[highest, ])

pdf("~/capsule/results/plots/top15_bin_aac_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
qp <- barplot(qq, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(qq), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of AAC of original values \nand after binary drug classification for the 15 drugs with largest changes in correlation", font.main = 1)
legend("topright", legend=c("orig", "binary"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(qp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(qq)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

ss <- as.matrix(t(cbind(binary.aac.cor[lowest, "pearson.x"], binary.aac.cor[lowest, "pearson.y"])))
names(ss) <- rownames(binary.aac.cor[lowest, ])

pdf("~/capsule/results/plots/lowest15_bin_aac_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
sp <- barplot(ss, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(ss), v=0.9), each=2), ylim = c(-1, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of AAC of original values \nand after binary drug classification for the 15 drugs with lowest changes in correlation", font.main = 1)
legend("topleft", legend=c("orig", "binary"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(sp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(ss)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# Bar plot of Pearson coefficients for aac values for the 15 drugs with highest and lowest correlations after PCL filtering
y <- na.omit(pcl.aac.cor[with(pcl.aac.cor,order(pcl.aac.cor$pearson)), ])
highest <- rownames(tail(y, n=15))
lowest <- rownames(head(y, n=15))

yy <- as.matrix(t(cbind(orig.aac.cor[highest, "pearson"], pcl.aac.cor[highest, "pearson"])))
names(yy) <- rownames(orig.aac.cor[highest, ])

pdf("~/capsule/results/plots/top15_pcl_aac_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
yp <- barplot(yy, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(yy), v=0.9), each=2), ylim = c(0, 1), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of AAC of original values \nand after PCL drug classification for the 15 drugs with highest correlations", font.main = 1)
legend("topleft", legend=c("orig", "PCL"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(yp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(yy)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

zz <- as.matrix(t(cbind(orig.aac.cor[lowest, "pearson"], pcl.aac.cor[lowest, "pearson"])))
names(zz) <- rownames(orig.aac.cor[lowest, ])

pdf("~/capsule/results/plots/lowest15_pcl_aac_bar_plot.pdf", height = 8, width = 10)
par(mar=c(8,5,5,5))
zp <- barplot(zz, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(zz), v=0.9), each=2), ylim = c(-0.4, 0.8), ylab="Pearson's correlation coefficient r", density=c(100,10) , angle=c(0,45), main="Pearson's correlation coefficient of aac of original values \nand after PCL drug classification for the 15 drugs with lowest correlations", font.main = 1)
legend("topleft", legend=c("orig", "PCL"), density=c(100,10), angle=c(0,45), bty="n", cex=0.75)
text(x=apply(zp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.01), pos=2, labels=toupper(names(zz)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# 12) bar plot representing all concordance indexes
all_concordances = cbind("IC50" = rbind(binary.ic50.concordance$concordance, pcl.ic50.concordance$concordance), "AAC" = rbind(binary.aac.concordance$concordance, pcl.aac.concordance$concordance))
names(all_concordances) <- c("IC50", "AUC")

barplot(as.matrix(all_concordances))

pdf("~/capsule/results/plots/ic50_concordances_bar_plot.pdf")
ylim <- c(0, 1.2*max(all_ic50_concordances))
ee <- barplot(all_ic50_concordances, main = "Concordance between Pearson correlations of IC50 values", width = 0.9, ylim = ylim, ylab = "concordance index")
text(x = ee, y = all_ic50_concordances, label = signif(all_ic50_concordances, 3), pos = 3, cex = 0.8, col = "black")
dev.off()

all_aac_concordances = c(binary.aac.concordance$concordance, pcl.aac.concordance$concordance)
names(all_aac_concordances) <- c("Broad-Spectrum/Targeted", "PCL")

pdf("~/capsule/results/plots/aac_concordances_bar_plot.pdf")
ylim <- c(0, 1.2*max(all_aac_concordances))
ff <- barplot(all_aac_concordances, main = "Concordance between Pearson correlations of AAC values", width = 0.9, ylim = ylim, ylab = "concordance index")
text(x = ff, y = all_aac_concordances, label = signif(all_aac_concordances, 3), pos = 3, cex = 0.8, col = "black")
dev.off()
