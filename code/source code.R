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
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@cell), area2=nrow(GDSC@cell), cross.area=nrow(ccle@cell), col = "black", cex=1.5, cat.cex=1, cat.col = c("black", "black"))
pdf("~/capsule/results/cell_intersection.pdf", height=4, width=4)
grid::grid.draw(venn.plot)
dev.off()
rm(venn.plot)

# venn diagram of common drugs
venn.plot2 <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@drug), area2=nrow(GDSC@drug), cross.area=nrow(ccle@drug), col = "black" ,cex=1.5, cat.cex=1, cat.col = c("black", "black"))
pdf("~/capsule/results/drug_intersection.pdf", height=4, width=4)
grid::grid.draw(venn.plot)
dev.off()    
rm(venn.plot2)

rm(CCLE)
rm(GDSC)
rm(intersected)

coef <- c("pearson", "pearson p-value", "spearman", "spearman p-value")
common_drugs <- ccle@drug

# classes from PSet
# pcl_count <- plyr::count(ccle@drug$Mechanism.of.action) # there's a lot of classes here but each one only has 1 drug in it... probs don't use this
# class_count <- plyr::count(ccle@drug$Class) # there's only 3 classes and one of them just has 1 in it.. probs not gonna use this either

# ====================================================
# loading manually curated list of drug information
# ====================================================
druglist <- as.data.frame(read.csv("PSets/druglist330.csv")) 
# TO-DO: manually fix rownames of druglist so they match ccle@drug rownames (in excel), and make this csv actually presentable

druglist[] <- lapply(druglist, as.character) # change class of dataframe columns from vector to character
rownames(druglist) <- rownames(ccle@drug)
# rownames(druglist) <- druglist$drug_name 
# TO-DO:
# pcl_count2 = plyr::count(druglist$`clue.io PCL`)
# rownames(pcl_count2) <- pcl_count2$x

# ====================================================
# sort drugs into broad spectrum/targeted classes
# ====================================================
classes = data.frame(drug = rownames(ccle@drug), row.names = rownames(ccle@drug))
classes$BroadSpectrum_or_Targeted <- NA

for (drug in rownames(druglist)) {
  tryCatch(
    expr = {
      if (druglist[drug, "target..drug.mechanism.from.ChEMBL."] == "single protein") {
        classes[drug, "BroadSpectrum_or_Targeted"] <- "targeted"
      }
      else if (druglist[drug, "target..drug.mechanism.from.ChEMBL."] == "protein family") {
        classes[drug, "BroadSpectrum_or_Targeted"] <- "broad-spectrum"
      }
      else if (is.na(druglist[drug, "target..drug.mechanism.from.ChEMBL."])) {
        classes[drug, "BroadSpectrum_or_Targeted"] <- NA #TO-DO: fix this (maybe just redo the druglist column to be broad-spectrum)
      }
    },
    error = function(e){
      message(drug, " is unclassified")
    }
  )
}

# ====================================================
#       sort drugs into clue.io defined PCLs
# ====================================================
classes$PCL <- NA
for (drug in rownames(druglist)) {
  classes[drug, "PCL"] <- druglist[drug, "clue.io.PCL"]
}

# ==============================================
#       get published ic50 and auc
# ==============================================
ccle_ic50 = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published"))
# TO-DO: remove columns with all NA's?

ccle_auc = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "auc_published"))
# TO-DO: remove columns with all NA's?

gdsc_ic50 = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published"))
# TO-DO: remove columns with all NA's?

gdsc_auc = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published"))
# TO-DO: remove columns with all NA's?

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
# TO-DO: this doesn't work, intersect whatever drugs aren't in broad_spectrum_drugs and targeted_drugs
# unclassified_drugs <- classes[grep(NA, classes$BroadSpectrum_or_Targeted), "drug"] # drugs that aren't sorted

# ====================================================
#               sort drugs into PCLs 
# ====================================================
# make list of all PCLs and the drugs in each one
# TO-DO: make this nicer...
pcl_count <- dplyr::select(classes, PCL)
pcl_count$drug <- rownames(pcl_count)
pcls <- (reshape2::dcast(pcl_count, drug~PCL))
pcls <- pcls[, -1]

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
# ===== inconsistent based on broad-spectrum drugs =====
min = 0.5
rmv.ic50.brsp.cor <- ic50.brsp.cor[which(ic50.brsp.cor$spearman > min), ] # subset by removing cell lines that had low correlation (Spearman less than 0.5)
brsp.ic50.cc <- rownames(rmv.ic50.brsp.cor)
rmv.auc.brsp.cor <- auc.brsp.cor[which(auc.brsp.cor$spearman > min), ]
brsp.auc.cc <- rownames(rmv.auc.brsp.cor)
  
# ===== inconsistent based on targeted drugs =====
rmv.ic50.targ.cor <- ic50.targ.cor[which(ic50.targ.cor$spearman > min), ]
targ.ic50.cc <- rownames(rmv.ic50.targ.cor)
rmv.auc.targ.cor <- auc.targ.cor[which(auc.targ.cor$spearman > min), ]
targ.auc.cc <- rownames(rmv.auc.targ.cor)

# ===== inconsistent based on pcls =====
# MEK inhibitor
mek_cols = c("MEK inhibitor.pearson", "MEK inhibitor.pearson p-value", "MEK inhibitor.spearman", "MEK inhibitor.spearman p-value")
rmv.ic50.mek.cor <- ic50.pcl.cor[which(ic50.pcl.cor$`MEK inhibitor.spearman` > min), mek_cols]
mek.ic50.cc <- rownames(rmv.ic50.mek.cor)
rmv.auc.mek.cor <- auc.pcl.cor[which(auc.pcl.cor$`MEK inhibitor.spearman` > min), mek_cols]
mek.auc.cc <- rownames(rmv.auc.mek.cor)

# SRC inhibitor
src_cols = c("SRC inhibitor.pearson", "SRC inhibitor.pearson p-value", "SRC inhibitor.spearman", "SRC inhibitor.spearman p-value")
rmv.ic50.src.cor <- ic50.pcl.cor[which(ic50.pcl.cor$`SRC inhibitor.spearman` > min), src_cols]
src.ic50.cc <- rownames(rmv.ic50.src.cor)
rmv.auc.src.cor <- auc.pcl.cor[which(auc.pcl.cor$`SRC inhibitor.spearman` > min), src_cols]
src.auc.cc <- rownames(rmv.auc.src.cor)

# ====================================================
#    remove inconsistent cell lines within GDSC 
# ====================================================
# filtering based on broad-spectrum classification
gdsc_ic50_cc_brsp = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published", cell.lines = brsp.ic50.cc))
gdsc_auc_cc_brsp = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published", cell.lines = brsp.auc.cc))

# filtering based on targeted classification
gdsc_ic50_cc_targ = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published", cell.lines = targ.ic50.cc))
gdsc_auc_cc_targ = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published", cell.lines = targ.auc.cc))

# filtering based on MEK inhibitor 
gdsc_ic50_cc_mek = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published", cell.lines = mek.ic50.cc))
gdsc_auc_cc_mek = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published", cell.lines = mek.auc.cc))

# filtering based on SRC inhibitor
gdsc_ic50_cc_src = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published", cell.lines = src.ic50.cc))
gdsc_auc_cc_src = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published", cell.lines = src.auc.cc))

# ====================================================
#    remove inconsistent cell lines within CCLE 
# ====================================================
# filtering based on broad-spectrum classification
ccle_ic50_cc_brsp = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published", cell.lines = brsp.ic50.cc))
ccle_auc_cc_brsp = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "auc_published", cell.lines = brsp.auc.cc))

# filtering based on targeted classification
ccle_ic50_cc_targ = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published", cell.lines = targ.ic50.cc))
ccle_auc_cc_targ = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "auc_published", cell.lines = targ.auc.cc))

# filtering based on MEK inhibitor 
ccle_ic50_cc_mek = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published", cell.lines = mek.ic50.cc))
ccle_auc_cc_mek = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "auc_published", cell.lines = mek.auc.cc))

# filtering based on SRC inhibitor
ccle_ic50_cc_src = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published", cell.lines = src.ic50.cc))
ccle_auc_cc_src = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "auc_published", cell.lines = src.auc.cc))

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

# ======================================================================================
#   correlation between CCLE and GDSC after sorting into broad-spectrum/targeted
# ======================================================================================
# ==== correlation after binary sorting ====
# broad-spectrum 
# IC50
brsp.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_ic50_cc_brsp))), row.names = rownames(ccle_ic50_cc_brsp))
colnames(brsp.ic50.cor) <- coef

for (drug in rownames(brsp.ic50.cor)) {
  pearson.cor <- cor.test(x = as.numeric(ccle_ic50_cc_brsp[drug, ]), y = as.numeric(gdsc_ic50_cc_brsp[drug, ]), method = 'pearson', use = 'pairw')
  brsp.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
  brsp.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_ic50_cc_brsp[drug, ]), y = as.numeric(gdsc_ic50_cc_brsp[drug, ]), method = 'spearman', use = 'pairw')
  brsp.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
  brsp.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# AUC
brsp.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_auc_cc_brsp))), row.names = rownames(ccle_auc_cc_brsp))
colnames(brsp.auc.cor) <- coef

for (drug in rownames(brsp.auc.cor)) {  # did not have enough obsrvations as only one cell line had enough sens measures, so all are NA
  pearson.cor <- cor.test(x = as.numeric(ccle_auc_cc_brsp[drug, ]), y = as.numeric(gdsc_auc_cc_brsp[drug, ]), method = 'pearson', use = 'pairw')
  brsp.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  brsp.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_auc_cc_brsp[drug, ]), y = as.numeric(gdsc_auc_cc_brsp[drug, ]), method = 'spearman', use = 'pairw')
  brsp.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  brsp.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# targeted
# IC50
targ.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_ic50_cc_targ))), row.names = rownames(ccle_ic50_cc_targ))
colnames(targ.ic50.cor) <- coef

for (drug in rownames(brsp.ic50.cor)) {
  pearson.cor <- cor.test(x = as.numeric(ccle_ic50_cc_targ[drug, ]), y = as.numeric(gdsc_ic50_cc_targ[drug, ]), method = 'pearson', use = 'pairw')
  targ.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
  targ.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_ic50_cc_targ[drug, ]), y = as.numeric(gdsc_ic50_cc_targ[drug, ]), method = 'spearman', use = 'pairw')
  targ.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
  targ.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# AUC
targ.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_auc_cc_targ))), row.names = rownames(ccle_auc_cc_targ))
colnames(targ.auc.cor) <- coef

for (drug in rownames(targ.auc.cor)) {
  pearson.cor <- cor.test(x = as.numeric(ccle_auc_cc_targ[drug, ]), y = as.numeric(gdsc_auc_cc_targ[drug, ]), method = 'pearson', use = 'pairw')
  targ.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  targ.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_auc_cc_targ[drug, ]), y = as.numeric(gdsc_auc_cc_targ[drug, ]), method = 'spearman', use = 'pairw')
  targ.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  targ.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# ======================================================================================
#           Harrell's CI after sorting into broad-spectrum/targeted
# ======================================================================================
# IC50
# broad-spectrum:
orig.brsp.ic50.cor <- merge(orig.ic50.cor, brsp.ic50.cor, by = 0, all = TRUE)
marray::write.xls(orig.brsp.ic50.cor, "~/capsule/results/broad_sp_ic50_corr.xls", row.names = TRUE, col.names = TRUE)

orig.brsp.ic50.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.brsp.ic50.cor)
# should i do concordance between spearman?

# n= 15 
# Concordance= 0.7714 se= 0.08645
# concordant discordant     tied.x     tied.y    tied.xy 
# 81         24          0          0          0 

saveRDS(orig.brsp.ic50.conc, "~/capsule/results/orig.brsp.ic50.conc.rds")

# targeted:
orig.targ.ic50.cor <- merge(orig.ic50.cor, targ.ic50.cor, by = 0, all = TRUE)
# TO-DO: make plot

marray::write.xls(orig.targ.ic50.cor, "~/capsule/results.targeted_ic50_corr.xls", row.names = TRUE, col.names = TRUE)

orig.targ.ic50.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.targ.ic50.cor)

# n= 15 
# Concordance= 0.8476 se= 0.0464
# concordant discordant     tied.x     tied.y    tied.xy 
# 89         16          0          0          0 

saveRDS(orig.targ.ic50.conc, "~/capsule/results/orig.targ.ic50.conc.rds")

# AUC
# broad-spectrum:
orig.brsp.auc.cor <- merge(orig.auc.cor, brsp.auc.cor, by = 0, all = TRUE)
marray::write.xls(orig.brsp.auc.cor, "~/capsule/results/orig.brsp.auc.cor.xls", row.names = TRUE, col.names = TRUE)

# orig.brsp.auc.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.brsp.auc.cor)
# no concordance since brsp.auc.cor has no observations
# TO-DO: make plot
# saveRDS(orig.brsp.auc.conc, "~/capsule/results/orig.brsp.auc.conc.rds")

# targeted:
orig.targ.auc.cor <- merge(orig.auc.cor, targ.auc.cor, by = 0, all = TRUE)
marray::write.xls(orig.targ.auc.cor, "~/capsule/results/orig.targ.auc.cor.xls", row.names = TRUE, col.names = TRUE)

orig.targ.auc.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.targ.auc.cor)

# n= 15 
# Concordance= 0.7429 se= 0.08497
# concordant discordant     tied.x     tied.y    tied.xy 
# 78         27          0          0          0 

# TO-DO: make plot

saveRDS(orig.targ.auc.conc, "~/capsule/results/orig.targ.auc.conc.rds")

# =============================================================================
#       correlation between CCLE and GDSC after sorting into PCLs
# =============================================================================
# ================ MEK inhibitor class ================
# IC50
mek.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_ic50_cc_mek))), row.names = rownames(ccle_ic50_cc_mek))
colnames(mek.ic50.cor) <- coef

for (drug in rownames(mek.ic50.cor)) {
    pearson.cor <- cor.test(x = as.numeric(ccle_ic50_cc_mek[drug, ]), y = as.numeric(gdsc_ic50_cc_mek[drug, ]), method = 'pearson', use = 'pairw')
  mek.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
  mek.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_ic50_cc_mek[drug, ]), y = as.numeric(gdsc_ic50_cc_mek[drug, ]), method = 'spearman', use = 'pairw')
  mek.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
  mek.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# AUC
mek.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_auc_cc_mek))), row.names = rownames(ccle_auc_cc_mek))
colnames(mek.auc.cor) <- coef

for (drug in rownames(mek.auc.cor)) {
    pearson.cor <- cor.test(x = as.numeric(ccle_auc_cc_mek[drug, ]), y = as.numeric(gdsc_auc_cc_mek[drug, ]), method = 'pearson', use = 'pairw')
  mek.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  mek.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_auc_cc_mek[drug, ]), y = as.numeric(gdsc_auc_cc_mek[drug, ]), method = 'spearman', use = 'pairw')
  mek.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  mek.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# ================ SRC inhibitor class ================
# IC50
src.ic50.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_ic50_cc_src))), row.names = rownames(ccle_ic50_cc_src))
colnames(src.ic50.cor) <- coef

for (drug in rownames(src.ic50.cor)) {
    pearson.cor <- cor.test(x = as.numeric(ccle_ic50_cc_src[drug, ]), y = as.numeric(gdsc_ic50_cc_src[drug, ]), method = 'pearson', use = 'pairw')
  src.ic50.cor[drug, "pearson"] <- pearson.cor$estimate
  src.ic50.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_ic50_cc_src[drug, ]), y = as.numeric(gdsc_ic50_cc_src[drug, ]), method = 'spearman', use = 'pairw')
  src.ic50.cor[drug, "spearman"] <- spearman.cor$estimate
  src.ic50.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# AUC
src.auc.cor = data.frame(matrix(ncol = length(coef), nrow = length(rownames(ccle_auc_cc_src))), row.names = rownames(ccle_auc_cc_src))
colnames(mek.auc.cor) <- coef

for (drug in rownames(mek.auc.cor)) {
    pearson.cor <- cor.test(x = as.numeric(ccle_auc_cc_src[drug, ]), y = as.numeric(gdsc_auc_cc_src[drug, ]), method = 'pearson', use = 'pairw')
  src.auc.cor[drug, "pearson"] <- pearson.cor$estimate
  src.auc.cor[drug, "pearson p-value"] <- pearson.cor$p.value
  
  spearman.cor <- cor.test(x = as.numeric(ccle_auc_cc_src[drug, ]), y = as.numeric(gdsc_auc_cc_src[drug, ]), method = 'spearman', use = 'pairw')
  src.auc.cor[drug, "spearman"] <- spearman.cor$estimate
  src.auc.cor[drug, "spearman p-value"] <- spearman.cor$p.value
}

# ====================================================================================
#               Harrell's CI after sorting into PCLs
# ====================================================================================
# ================ MEK inhibitor class ================
# IC50
orig.mek.ic50.cor <- merge(orig.ic50.cor, mek.ic50.cor, by = 0, all = TRUE)
# TO-DO: make plot
marray::write.xls(orig.mek.ic50.cor, "~/capsule/results/orig.mek.ic50.cor.xls", row.names = TRUE, col.names = TRUE)

orig.mek.ic50.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.mek.ic50.cor)

# n= 15 
# Concordance= 0.7905 se= 0.06929
# concordant discordant     tied.x     tied.y    tied.xy 
# 83         22          0          0          0 

saveRDS(orig.mek.ic50.conc, "~/capsule/results/orig.mek.ic50.conc.rds")

# AUC
orig.mek.auc.cor <- merge(orig.auc.cor, mek.auc.cor, by = 0, all = TRUE)
marray::write.xls(orig.mek.auc.cor, "~/capsule/results/orig.mek.auc.cor.xls", row.names = TRUE, col.names = TRUE)

orig.mek.auc.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.mek.auc.cor)

# n= 15 
# Concordance= 0.7905 se= 0.08568
# concordant discordant     tied.x     tied.y    tied.xy 
# 83         22          0          0          0 

# TO-DO: make plot

saveRDS(orig.mek.auc.conc, "~/capsule/results/orig.mek.auc.conc.rds")

# ================ SRC inhibitor class ================
# IC50
orig.src.ic50.cor <- merge(orig.ic50.cor, src.ic50.cor, by = 0, all = TRUE)
marray::write.xls(orig.src.ic50.cor, "~/capsule/results/orig.src.ic50.cor.xls", row.names = TRUE, col.names = TRUE)

orig.src.ic50.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.src.ic50.cor)

# n=14 (1 observation deleted due to missingness)
# Concordance= 0.7692 se= 0.06774
# concordant discordant     tied.x     tied.y    tied.xy 
# 70         21          0          0          0 

# TO-DO: make plot

saveRDS(orig.src.ic50.conc, "~/capsule/results/orig.src.ic50.conc.rds")

# AUC
orig.src.auc.cor <- merge(orig.auc.cor, src.auc.cor, by = 0, all = TRUE)
marray::write.xls(orig.src.auc.cor, "~/capsule/results/orig.src.auc.cor.xls", row.names = TRUE, col.names = TRUE)

orig.src.auc.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.src.auc.cor)

#n= 15 
# Concordance= 0.6571 se= 0.09447
# concordant discordant     tied.x     tied.y    tied.xy 
# 69         36          0          0          0 

# TO-DO: make plot

saveRDS(orig.src.auc.conc, "~/capsule/results/orig.src.auc.conc.rds")

# ======================================================
# 1) Bar plot representing the Spearman's rank correlation coefficient for IC50 drug sensitivity measures
ss <- as.matrix(t(cbind(orig.ic50.cor$spearman, brsp.ic50.cor$spearman, targ.ic50.cor$spearman, src.ic50.cor$spearman, mek.ic50.cor$spearman)))
ss[!is.na(ss) & ss < 0] <- 0
names(ss) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/all_ic50_bar_plot.pdf")
mp <- barplot(ss, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(ss), v=0.9), each=2), ylab=expression("r"[s]), density=c(0,10,20,30,7) , angle=c(0,45,90,11,36), main="Spearman's rank correlation coefficient of IC50", font.main = 1)
legend("topright", legend=c("orig", "br-sp", "targ", "src", "mek"), density=c(0,10,20,30,7), angle=c(0,45,90,11,36), bty="n", cex=0.75)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(ss)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# 2) Bar plot representing the Spearman's rank correlation coefficient for AUC drug sensitivity measures; significance is reported using an asterisk if two‐sided P value <0.05.
# x axis = the drugs (which are the rownames of brsp.auc.cor), y axis = the spearman rank
# each drug has two bars, one with original correlation, one with new correlation when broad spectrum, one with new corr when targ, one with new corr when mek, one with new corr when src
aa <- as.matrix(t(cbind(orig.auc.cor$spearman, brsp.auc.cor$spearman, targ.auc.cor$spearman, src.auc.cor$spearman, mek.auc.cor$spearman)))
aa[!is.na(aa) & aa < 0] <- 0
names(aa) <- rownames(orig.auc.cor)

pdf("~/capsule/results/plots/all_auc_bar_plot.pdf")
bp <- barplot(aa, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(aa), v=0.9), each=2), ylab=expression("r"[s]), density=c(0,10,20,30,7) , angle=c(0,45,90,11,36), main="Spearman's rank correlation coefficient of AUC", font.main = 1)
legend("topright", legend=c("orig", "br-sp", "targ", "src", "mek"), density=c(0,10,20,30,7), angle=c(0,45,90,11,36), bty="n", cex=0.75)
text(x=apply(bp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=1, labels=toupper(names(aa)), srt=45, xpd=NA, font=1, cex = 0.75)
dev.off()

# 3) scatter plot (for each drug) reporting drug sensitivity of all cell line
# x axis = CCLE reported measures, y axis = GDSC reported measures
pdf("~/capsule/results/plots/drugs_auc_cor_barplot.pdf", height = 14, width = 14)
par(mfrow=c(4,4))
for (drug in rownames(gdsc_auc)) {
  plot(as.numeric(ccle_auc[drug, ]), as.numeric(gdsc_auc[drug, ]), main = drug, xlab = "ccle", ylab = "gdsc", type = "p")
}
dev.off()

# 5) Bar plot representing spearman's ranks for AUC drug measure comparing original with targeted
tg <- as.matrix(t(cbind(orig.auc.cor$spearman, targ.auc.cor$spearman)))
tg[!is.na(tg) & tg < 0] <- 0
names(tg) <- rownames(orig.auc.cor)

pdf("~/capsule/results/plots/targeted_auc_cor_bar_plot.pdf")
tgb <- barplot(tg, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(tg), v=0.9), each=2), ylab=expression("r"[s]), angle=c(45, -45), density=c(100, 40), main="Spearman's rank correlation coefficient of AUC", font.main = 1)
legend("topright", legend=c("orig", "targeted"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
text(x=apply(tgb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(tg)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 6) Bar plot representing spearman's ranks for AUC drug measure comparing original with mek drugs
mek <- as.matrix(t(cbind(orig.auc.cor$spearman, mek.auc.cor$spearman)))
mek[!is.na(mek) & mek < 0] <- 0
names(mek) <- rownames(orig.auc.cor)

pdf("~/capsule/results/plots/MEK_auc_cor_bar_plot.pdf")
mb <- barplot(mek, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(mek), v=0.9), each=2), ylab=expression("r"[s]), angle=c(45, -45), density=c(100, 40), main="Spearman's rank correlation coefficient for AUC", font.main = 1)
legend("topright", legend=c("orig", "targeted"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
text(x=apply(mb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(mek)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 7) Bar plot representing spearman's ranks for AUC drug measure comparing original with src drugs
src <- as.matrix(t(cbind(orig.auc.cor$spearman, src.auc.cor$spearman)))
src[!is.na(src) & src < 0] <- 0
names(src) <- rownames(orig.auc.cor)

pdf("~/capsule/results/plots/SRC_auc_cor_bar_plot.pdf")
sb <- barplot(src, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(mek), v=0.9), each=2), ylab=expression("r"[s]), angle=c(45, -45), density=c(100, 40), main="Spearman's rank correlation coefficient for AUC", font.main = 1)
legend("topright", legend=c("orig", "targeted"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
text(x=apply(sb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(src)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 8) Bar plot representing spearman's ranks for IC50 drug measure comparing original with broad spectrum
aa <- as.matrix(t(cbind(orig.ic50.cor$spearman, src.ic50.cor$spearman)))
aa[!is.na(aa) & aa < 0] <- 0
names(aa) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/brsp_ic50_cor_bar_plot.pdf")
aab <- barplot(aa, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(aa), v=0.9), each=2), ylab=expression("r"[s]), angle=c(45, -45), density=c(100, 40), main="Spearman's rank correlation coefficient for IC50", font.main = 1)
legend("topright", legend=c("orig", "targeted"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
text(x=apply(aab, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(aa)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 9) Bar plot representing spearman's ranks for IC50 drug measure comparing original with targeted
bb <- as.matrix(t(cbind(orig.ic50.cor$spearman, targ.ic50.cor$spearman)))
bb[!is.na(bb) & bb < 0] <- 0
names(bb) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/targ_ic50_cor_bar_plot.pdf")
bbb <- barplot(bb, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(bb), v=0.9), each=2), ylab=expression("r"[s]), angle=c(45, -45), density=c(100, 40), main="Spearman's rank correlation coefficient for IC50", font.main = 1)
legend("topright", legend=c("orig", "targeted"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
text(x=apply(bbb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(bb)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 10) Bar plot representing spearman's ranks for IC50 drug measure comparing original with mek drugs
cc <- as.matrix(t(cbind(orig.ic50.cor$spearman, mek.ic50.cor$spearman)))
cc[!is.na(cc) & cc < 0] <- 0
names(cc) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/mek_ic50_cor_bar_plot.pdf")
ccb <- barplot(cc, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(cc), v=0.9), each=2), ylab=expression("r"[s]), angle=c(45, -45), density=c(100, 40), main="Spearman's rank correlation coefficient for IC50", font.main = 1)
legend("topright", legend=c("orig", "targeted"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
text(x=apply(ccb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(cc)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 11) Bar plot representing spearman's ranks for IC50 drug measure comparing original with src drugs
dd <- as.matrix(t(cbind(orig.ic50.cor$spearman, src.ic50.cor$spearman)))
dd[!is.na(dd) & dd < 0] <- 0
names(dd) <- rownames(orig.ic50.cor)

pdf("~/capsule/results/plots/src_ic50_cor_bar_plot.pdf")
ddb <- barplot(dd, beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(dd), v=0.9), each=2), ylab=expression("r"[s]), angle=c(45, -45), density=c(100, 40), main="Spearman's rank correlation coefficient for IC50", font.main = 1)
legend("topright", legend=c("orig", "targeted"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
text(x=apply(ddb, 2, mean), y=par("usr")[3] - (par("usr")[4]*0.05), pos=1, labels=toupper(names(dd)), srt=50, xpd=NA, font=1, cex = 0.75)
dev.off()

# 12) bar plot representing all concordance indexes
all_ic50_concordances = c(orig.brsp.ic50.conc$concordance, orig.targ.ic50.conc$concordance, orig.mek.ic50.conc$concordance, orig.src.ic50.conc$concordance)
names(all_ic50_concordances) <- c("Broad", "Targ","MEK Inhib", "SRC Inhib")

pdf("~/capsule/results/plots/ic50_concordances_bar_plot.pdf")
ee <- barplot(all_ic50_concordances, main = "Concordance between Pearson correlations of IC50")
dev.off()

all_auc_concordances = c(orig.targ.auc.conc$concordance, orig.mek.auc.conc$concordance, orig.src.auc.conc$concordance)
names(all_auc_concordances) <- c("Targ","MEK Inhib", "SRC Inhib")

pdf("~/capsule/results/plots/auc_concordances_bar_plot.pdf")
ff <- barplot(all_auc_concordances, main = "Concordance between Pearson correlations of AUC")
dev.off()
