library(PharmacoGx); library(dplyr); library(reshape2); library(survival); library(WriteXLS); library(VennDiagram) #TO-DO: do the required thing to install packages
# downloadPSet("GDSC_2013"); downloadPSet("CCLE_2013")

# ==== loading PSets and intersecting CCLE and GDSC for common drugs and cell lines ====
load("~/capsule/code/PSets/CCLE_2013.RData"); load("~/capsule/code/PSets/GDSC_2013.RData")
intersected <- intersectPSet(c(CCLE, GDSC), intersectOn = c("drugs", "cell.lines"))
# ^ 15 drugs and 514 cell lines in common!
ccle <- intersected$CCLE
gdsc <- intersected$GDSC

# venn diagram of common cell lines
pdf("~/capsule/results/cell_intersection.pdf", height=4, width=4)
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@cell), area2=nrow(GDSC@cell), cross.area=nrow(ccle@cell), fill=c(mycol[1], mycol[2]), lty="blank",cex=1.5, cat.cex=1, cat.col = c("black", "black"))
dev.off()

# venn diagram of common drugs
pdf("~/capsule/results/drug_intersection.pdf"), height=4, width=4)
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@drug), area2=nrow(GDSC@drug), cross.area=nrow(ccle@drug), fill=c(mycol[1], mycol[2]), lty="blank",cex=1.5, cat.cex=1, cat.col = c("black", "black"))
dev.off()      

rm(CCLE)
rm(GDSC)
rm(intersected)

# classes from PSet
# pcl_count <- plyr::count(ccle@drug$Mechanism.of.action) # there's a lot of classes here but each one only has 1 drug in it... probs don't use this
# class_count <- plyr::count(ccle@drug$Class) # there's only 3 classes and one of them just has 1 in it.. probs not gonna use this either

# ==== loading manually curated list of drug information ====
druglist <- as.data.frame(read.csv("PSets/druglist330.csv")) 
# TO-DO: manually fix rownames of druglist so they match ccle@drug rownames (in excel), and make this csv actually presentable

druglist[] <- lapply(druglist, as.character) # change class of dataframe columns from vector to character
rownames(druglist) <- rownames(ccle@drug)
# rownames(druglist) <- druglist$drug_name # TO-DO:
# pcl_count2 = plyr::count(druglist$`clue.io PCL`)
# rownames(pcl_count2) <- pcl_count2$x

# ==== sort drugs into broad spectrum/targeted classes ====
classes = data.frame(drug = rownames(ccle@drug), row.names = rownames(ccle@drug))
classes$BroadSpectrum_or_Targeted <- NA

for (drug in rownames(druglist[14:15,])) { # to-do: put into try catch
  if (druglist[drug, "target..drug.mechanism.from.ChEMBL."] == "single protein") {
    classes[drug, "BroadSpectrum_or_Targeted"] <- "targeted"
  }
  else if (druglist[drug, "target..drug.mechanism.from.ChEMBL."] == "protein family") {
    classes[drug, "BroadSpectrum_or_Targeted"] <- "broad-spectrum"
  }
  else if (is.na(druglist[drug, "target..drug.mechanism.from.ChEMBL."])) {
    classes[drug, "BroadSpectrum_or_Targeted"] <- NA #TO-DO: fix this (maybe just redo the druglist column to be broad-spectrum)
  }
}

# ==== sort drugs into clue.io defined PCLs ====
classes$PCL <- NA
for (drug in rownames(druglist)) {
  classes[drug, "PCL"] <- druglist[drug, "clue.io.PCL"]
}

ccle_ic50 = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published"))
# TO-DO: remove columns with all NA's?

ccle_auc = as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "auc_published"))
# TO-DO: remove columns with all NA's?

gdsc_ic50 = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "ic50_published"))
# TO-DO: remove columns with all NA's?

gdsc_auc = as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "auc_published"))
# TO-DO: remove columns with all NA's?

# ==== correlation of published IC50 between CCLE and GDSC ====
ic50.cor = data.frame(cell_line = colnames(ccle_ic50), row.names = colnames(ccle_ic50))
ic50.cor$pearson <- NA
ic50.cor$spearman <- NA

# TO-DO: simplify for-loops with lapply
for (cell_line in rownames(ic50.cor)) {
  ic50.cor[cell_line, "pearson"] <- cor(ccle_ic50[, cell_line], gdsc_ic50[, cell_line], method = "pearson", use = "pairwise.complete.obs")
  ic50.cor[cell_line, "spearman"] <- cor(ccle_ic50[, cell_line], gdsc_ic50[, cell_line], method = "spearman", use = "pairwise.complete.obs")
}

# ==== correlation of published AUC between CCLE and GDSC ====
# by cell line
auc.cor = data.frame(cell_line = colnames(ccle_auc), row.names = colnames(ccle_auc))
auc.cor$pearson <- NA
auc.cor$spearman <- NA
for (cell_line in rownames(auc.cor)) {
  auc.cor[cell_line, "pearson"] <- cor(ccle_auc[, cell_line], gdsc_auc[, cell_line], method = "pearson", use = "pairwise.complete.obs")
  auc.cor[cell_line, "spearman"] <- cor(ccle_auc[, cell_line], gdsc_auc[, cell_line], method = "spearman", use = "pairwise.complete.obs")
}

# ==== sort drugs into broad-spectrum vs. targeted ====
broad_spectrum_drugs <- classes[grep("broad-spectrum", classes$BroadSpectrum_or_Targeted), "drug"] 
targeted_drugs <- classes[grep("targeted", classes$BroadSpectrum_or_Targeted), "drug"] 
# TO-DO: this doesn't work, intersect whatever drugs aren't in broad_spectrum_drugs and targeted_drugs
# unclassified_drugs <- classes[grep(NA, classes$BroadSpectrum_or_Targeted), "drug"] # drugs that aren't sorted

# ==== sort drugs into PCLs ====
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

# ==== recompute sens measures after binary sorting ====
# === BROAD-SPECTRUM DRUGS ===
# correlation of IC50 in broad-spectrum drugs
ccle_ic50_brsp <- ccle_ic50[broad_spectrum_drugs, ] # ic50 for each cell lines and drugs that are broad-spectrum
gdsc_ic50_brsp <- gdsc_ic50[broad_spectrum_drugs, ] 
ic50.brsp.cor = data.frame(cell_line = colnames(ccle_ic50), row.names = colnames(ccle_ic50))
ic50.brsp.cor$pearson <- NA
ic50.brsp.cor$spearman <- NA
for (cell_line in rownames(ic50.brsp.cor)) {
  ic50.brsp.cor[cell_line, "pearson"] <- cor(ccle_ic50_brsp[, cell_line], gdsc_ic50_brsp[, cell_line], method = "pearson", use = "pairwise.complete.obs")
  ic50.brsp.cor[cell_line, "spearman"] <- cor(ccle_ic50_brsp[, cell_line], gdsc_ic50_brsp[, cell_line], method = "spearman", use = "pairwise.complete.obs")
}

# correlation of AUC in broad-spectrum drugs
ccle_auc_brsp <- ccle_auc[broad_spectrum_drugs, ] # auc for each cell lines and drugs that are broad-spectrum
gdsc_auc_brsp <- gdsc_auc[broad_spectrum_drugs, ] 
auc.brsp.cor = data.frame(cell_line = colnames(ccle_auc), row.names = colnames(ccle_auc))
auc.brsp.cor$pearson <- NA
auc.brsp.cor$spearman <- NA
for (cell_line in rownames(auc.brsp.cor)) {
  auc.brsp.cor[cell_line, "pearson"] <- cor(ccle_auc_brsp[, cell_line], gdsc_auc_brsp[, cell_line], method = "pearson", use = "pairwise.complete.obs")
  auc.brsp.cor[cell_line, "spearman"] <- cor(ccle_auc_brsp[, cell_line], gdsc_auc_brsp[, cell_line], method = "spearman", use = "pairwise.complete.obs")
}

# === TARGETED DRUGS ===
# correlation of IC50 in targeted drugs
ccle_ic50_targ <- ccle_ic50[targeted_drugs, ] # ic50 for each cell lines and drugs that are targeted
gdsc_ic50_targ <- gdsc_ic50[targeted_drugs, ] 
ic50.targ.cor = data.frame(cell_line = colnames(ccle_ic50), row.names = colnames(ccle_ic50))
ic50.targ.cor$pearson <- NA
ic50.targ.cor$spearman <- NA
for (cell_line in rownames(ic50.targ.cor)) {
  ic50.targ.cor[cell_line, "pearson"] <- cor(ccle_ic50_targ[, cell_line], gdsc_ic50_targ[, cell_line], method = "pearson", use = "pairwise.complete.obs")
  ic50.targ.cor[cell_line, "spearman"] <- cor(ccle_ic50_targ[, cell_line], gdsc_ic50_targ[, cell_line], method = "spearman", use = "pairwise.complete.obs")
}

# correlation of auc in targeted drugs
ccle_auc_targ <- ccle_auc[targeted_drugs, ] # auc for each cell lines and drugs that are targeted
gdsc_auc_targ <- gdsc_auc[targeted_drugs, ] 
auc.targ.cor = data.frame(cell_line = colnames(ccle_auc), row.names = colnames(ccle_auc))
auc.targ.cor$pearson <- NA
auc.targ.cor$spearman <- NA
for (cell_line in rownames(auc.targ.cor)) {
  auc.targ.cor[cell_line, "pearson"] <- cor(ccle_auc_targ[, cell_line], gdsc_auc_targ[, cell_line], method = "pearson", use = "pairwise.complete.obs")
  auc.targ.cor[cell_line, "spearman"] <- cor(ccle_auc_targ[, cell_line], gdsc_auc_targ[, cell_line], method = "spearman", use = "pairwise.complete.obs")
}

# ==== recompute sens measures after sorting into PCLs ====
# keeping only pcls that contain at least 3 drugs 
pcl_temp <- data.frame(matrix(NA, nrow = 3)) #TO-DO: remove all the hard-coded numbers/stuff
for (pcl in colnames(pcls)) {
  if (!length(pcl_list[[pcl]]) == 1) {
    pcl_temp[, pcl] <- pcl_list[[pcl]]
  }
}
pcl_temp <- pcl_temp[, -1]

# correlation of auc in drugs by pcl for each cell line 
auc.pcl.cor = data.frame(cell_line = colnames(ccle_auc), row.names = colnames(ccle_auc))

for (pcl in colnames(pcl_temp)) {
  auc.pcl.cor[, paste0(pcl, ".pearson")] <- NA
  auc.pcl.cor[, paste0(pcl, ".spearman")] <- NA
}

for (pcl in colnames(pcl_temp)) {
  ccle_auc_temp <- ccle_auc[pcl_list[[pcl]], ] # subset by the drugs in the PCL
  gdsc_auc_temp <- gdsc_auc[pcl_list[[pcl]], ]

  for (cell_line in rownames(auc.pcl.cor)) {
    auc.pcl.cor[cell_line, paste0(pcl, ".pearson")] <- cor(ccle_auc_temp[, cell_line], gdsc_auc_temp[, cell_line], method = "pearson", use = "pairwise.complete.obs")
    auc.pcl.cor[cell_line, paste0(pcl, ".spearman")] <- cor(ccle_auc_temp[, cell_line], gdsc_auc_temp[, cell_line], method = "spearman", use = "pairwise.complete.obs")
  }
}

# correlation of ic50 in drugs by pcl for each cell line
ic50.pcl.cor = data.frame(cell_line = colnames(ccle_auc), row.names = colnames(ccle_auc))
for (pcl in colnames(pcl_temp)) {
  ic50.pcl.cor[, paste0(pcl, ".pearson")] <- NA
  ic50.pcl.cor[, paste0(pcl, ".spearman")] <- NA
}

for (pcl in colnames(pcl_temp)) {
  ccle_ic50_temp <- ccle_ic50[pcl_list[[pcl]], ] # subset by the drugs in the PCL
  gdsc_ic50_temp <- gdsc_ic50[pcl_list[[pcl]], ]
  
  for (cell_line in rownames(ic50.pcl.cor)) {
    ic50.pcl.cor[cell_line, paste0(pcl, ".pearson")] <- cor(ccle_ic50_temp[, cell_line], gdsc_ic50_temp[, cell_line], method = "pearson", use = "pairwise.complete.obs")
    ic50.pcl.cor[cell_line, paste0(pcl, ".spearman")] <- cor(ccle_ic50_temp[, cell_line], gdsc_ic50_temp[, cell_line], method = "spearman", use = "pairwise.complete.obs")
  }
}

# ==== removing inconsistent cell lines ====
# inconsistent based on broad-spectrum drugs:
rmv.ic50.brsp.cor <- ic50.brsp.cor[which(ic50.brsp.cor$spearman > 0.5), ] # subset by removing cell lines that had low correlation (Spearman less than 0.5)
brsp.ic50.cc <- rownames(rmv.ic50.brsp.cor)
rmv.auc.brsp.cor <- auc.brsp.cor[which(auc.brsp.cor$spearman > 0.5), ]
brsp.auc.cc <- rownames(rmv.auc.brsp.cor)

# inconsistent based on targeted drugs:
rmv.ic50.targ.cor <- ic50.targ.cor[which(ic50.targ.cor$spearman > 0.5), ]
targ.ic50.cc <- rownames(rmv.ic50.targ.cor)
rmv.auc.targ.cor <- auc.targ.cor[which(auc.targ.cor$spearman > 0.5), ]
targ.auc.cc <- rownames(rmv.auc.targ.cor)

# inconsistent based on pcls:
# MEK inhibitor
rmv.ic50.mek.cor <- ic50.pcl.cor[which(ic50.pcl.cor$`MEK inhibitor.spearman` > 0.5), c("cell_line", "MEK inhibitor.pearson", "MEK inhibitor.spearman")]
mek.ic50.cc <- rownames(rmv.ic50.mek.cor)
rmv.auc.mek.cor <- auc.pcl.cor[which(auc.pcl.cor$`MEK inhibitor.spearman` > 0.5), c("cell_line", "MEK inhibitor.pearson", "MEK inhibitor.spearman")]
mek.auc.cc <- rownames(rmv.auc.mek.cor)

# SRC inhibitor
rmv.ic50.src.cor <- ic50.pcl.cor[which(ic50.pcl.cor$`SRC inhibitor.spearman` > 0.5), c("cell_line", "SRC inhibitor.pearson", "SRC inhibitor.spearman")]
src.ic50.cc <- rownames(rmv.ic50.src.cor)
rmv.auc.src.cor <- auc.pcl.cor[which(auc.pcl.cor$`SRC inhibitor.spearman` > 0.5), c("cell_line", "SRC inhibitor.pearson", "SRC inhibitor.spearman")]
src.auc.cc <- rownames(rmv.auc.src.cor)

# ==== remove inconsistent cell lines within GDSC ====
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

# ==== remove inconsistent cell lines within CCLE ====
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

# ==== original Pearson & Spearman correlation of sens measures of each common cell line between CCLE and GDSC ====
orig.ic50.cor = data.frame(drug = rownames(ccle_ic50), row.names = rownames(ccle_ic50))
orig.ic50.cor$pearson <- NA; orig.ic50.cor$spearman <- NA 
for (drug in rownames(orig.ic50.cor)) {
  orig.ic50.cor[drug, "pearson"] <- cor(as.numeric(ccle_ic50[drug, ]), as.numeric(gdsc_ic50[drug, ]), method = "pearson", use = "pairwise.complete.obs")
  orig.ic50.cor[drug, "spearman"] <- cor(as.numeric(ccle_ic50[drug, ]), as.numeric(gdsc_ic50[drug, ]), method = "spearman", use = "pairwise.complete.obs")
}

orig.auc.cor = data.frame(drug = rownames(ccle_auc), row.names = rownames(ccle_auc))
orig.auc.cor$pearson <- NA; orig.auc.cor$spearman <- NA
for (drug in rownames(orig.auc.cor)) {
  orig.auc.cor[drug, "pearson"] <- cor(as.numeric(ccle_auc[drug,]), as.numeric(gdsc_auc[drug,]), method = "pearson", use = "pairwise.complete.obs")
  orig.auc.cor[drug, "spearman"] <- cor(as.numeric(ccle_auc[drug,]), as.numeric(gdsc_auc[drug,]), method = "spearman", use = "pairwise.complete.obs")
}
# ==== Pearson & Spearman correlation between CCLE and GDSC after binary sorting ====
# broad-spectrum 
brsp.ic50.cor = data.frame(drug = rownames(ccle_ic50_cc_brsp), row.names = rownames(ccle_ic50_cc_brsp))
brsp.ic50.cor$pearson <- NA
brsp.ic50.cor$spearman <- NA
for (drug in rownames(brsp.ic50.cor)) {
  brsp.ic50.cor[drug, "pearson"] <- cor(as.numeric(ccle_ic50_cc_brsp[drug, ]), as.numeric(gdsc_ic50_cc_brsp[drug, ]), method = "pearson", use = "pairwise.complete.obs")
  brsp.ic50.cor[drug, "spearman"] <- cor(as.numeric(ccle_ic50_cc_brsp[drug, ]), as.numeric(gdsc_ic50_cc_brsp[drug, ]), method = "spearman", use = "pairwise.complete.obs")
}

brsp.auc.cor = data.frame(drug = rownames(ccle_auc_cc_brsp), row.names = rownames(ccle_auc_cc_brsp))
brsp.auc.cor$pearson <- NA
brsp.auc.cor$spearman <- NA
for (drug in rownames(brsp.auc.cor)) {
  brsp.auc.cor[drug, "pearson"] <- cor(as.numeric(ccle_auc_cc_brsp[drug, ]), as.numeric(gdsc_auc_cc_brsp[drug, ]), method = "pearson", use = "pairwise.complete.obs")
  brsp.auc.cor[drug, "spearman"] <- cor(as.numeric(ccle_auc_cc_brsp[drug, ]), as.numeric(gdsc_auc_cc_brsp[drug, ]), method = "spearman", use = "pairwise.complete.obs")
}

# targeted
targ.ic50.cor = data.frame(drug = rownames(ccle_ic50_cc_targ), row.names = rownames(ccle_ic50_cc_targ))
targ.ic50.cor$pearson <- NA
targ.ic50.cor$spearman <- NA
for (drug in rownames(targ.ic50.cor)) {
  targ.ic50.cor[drug, "pearson"] <- cor(as.numeric(ccle_ic50_cc_targ[drug, ]), as.numeric(gdsc_ic50_cc_targ[drug, ]), method = "pearson", use = "pairwise.complete.obs")
  targ.ic50.cor[drug, "spearman"] <- cor(as.numeric(ccle_ic50_cc_targ[drug, ]), as.numeric(gdsc_ic50_cc_targ[drug, ]), method = "spearman", use = "pairwise.complete.obs")
}

targ.auc.cor = data.frame(drug = rownames(ccle_auc_cc_targ), row.names = rownames(ccle_auc_cc_targ))
targ.auc.cor$pearson <- NA; targ.auc.cor$spearman <- NA
for (drug in rownames(targ.auc.cor)) {
  targ.auc.cor[drug, "pearson"] <- cor(as.numeric(ccle_auc_cc_targ[drug, ]), as.numeric(gdsc_auc_cc_targ[drug, ]), method = "pearson", use = "pairwise.complete.obs")
  targ.auc.cor[drug, "spearman"] <- cor(as.numeric(ccle_auc_cc_targ[drug, ]), as.numeric(gdsc_auc_cc_targ[drug, ]), method = "spearman", use = "pairwise.complete.obs")
}

# ==== Harrell's Concordance Index between CCLE and GDSC after binary sorting ====
# IC50:
# broad-spectrum filtering:
orig.brsp.ic50 <- merge(orig.ic50.cor, brsp.ic50.cor, by = 0, all = TRUE)
WriteXLS::WriteXLS(orig.brsp.ic50.cor, "~/capsule/results.broad_sp_ic50_corr.xlsx", row.names = TRUE, col.names = TRUE, AdjWidth = TRUE, BoldHeaderRow = TRUE)

orig.brsp.ic50.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.brsp)
saveRDS(orig.brsp.ic50.conc, "~/capsule/results/orig.brsp.ic50.conc.rds")

# targeted:
orig.targ.ic50 <- merge(orig.ic50.cor, targ.ic50.cor, by = 0, all = TRUE)
WriteXLS::WriteXLS(orig.targ.ic50, "~/capsule/results.targeted_ic50_corr.xlsx", row.names = TRUE, col.names = TRUE, AdjWidth = TRUE, BoldHeaderRow = TRUE)

orig.targ.ic50.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.targ)
saveRDS(orig.targ.ic50.conc, "~/capsule/results/orig.targ.ic50.conc.rds")

# AUC:
#brsp
orig.brsp.auc <- merge(orig.auc.cor, brsp.auc.cor, by = 0, all = TRUE)
WriteXLS::WriteXLS(orig.brsp.auc, "~/capsule/results.broad_sp_auc_corr.xlsx", row.names = TRUE, col.names = TRUE, AdjWidth = TRUE, BoldHeaderRow = TRUE)

orig.brsp.auc.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.brsp.auc)
saveRDS(orig.brsp.auc.conc, "~/capsule/results/orig.brsp.auc.conc.rds")

#targ
orig.targ.auc <- merge(orig.auc.cor, targ.auc.cor, by = 0, all = TRUE)
WriteXLS::WriteXLS(orig.targ.auc, "~/capsule/results.targ_auc_corr.xlsx", row.names = TRUE, col.names = TRUE, AdjWidth = TRUE, BoldHeaderRow = TRUE)

orig.targ.auc.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.targ.auc)
saveRDS(orig.targ.auc.conc, "~/capsule/results/orig.targ.auc.conc.rds")

# ==== Pearson & Spearman correlation between CCLE and GDSC after sorting into PCLs ====
# MEK inhibitor
mek.ic50.cor = data.frame(drug = rownames(ccle_ic50_cc_mek), row.names = rownames(ccle_ic50_cc_mek))
mek.ic50.cor$pearson <- NA; mek.ic50.cor$spearman <- NA
for (drug in rownames(mek.ic50.cor)) {
  mek.ic50.cor[drug, "pearson"] <- cor(as.numeric(ccle_ic50_cc_mek[drug, ]), as.numeric(gdsc_ic50_cc_mek[drug, ]), method = "pearson", use = "pairwise.complete.obs")
  mek.ic50.cor[drug, "spearman"] <- cor(as.numeric(ccle_ic50_cc_mek[drug, ]), as.numeric(gdsc_ic50_cc_mek[drug, ]), method = "spearman", use = "pairwise.complete.obs")
}

mek.auc.cor = data.frame(drug = rownames(ccle_auc_cc_mek), row.names = rownames(ccle_auc_cc_mek))
mek.auc.cor$pearson <- NA; mek.auc.cor$spearman <- NA
for (drug in rownames(mek.auc.cor)) {
  mek.auc.cor[drug, "pearson"] <- cor(as.numeric(ccle_auc_cc_mek[drug, ]), as.numeric(gdsc_auc_cc_mek[drug, ]), method = "pearson", use = "pairwise.complete.obs")
  mek.auc.cor[drug, "spearman"] <- cor(as.numeric(ccle_auc_cc_mek[drug, ]), as.numeric(gdsc_auc_cc_mek[drug, ]), method = "spearman", use = "pairwise.complete.obs")
}

# SRC inhibitor
src.ic50.cor = data.frame(cell_line = rownames(ccle_ic50_cc_src), row.names = rownames(ccle_ic50_cc_src))
src.ic50.cor$pearson <- NA; src.ic50.cor$spearman <- NA
for (drug in rownames(src.ic50.cor)) {
  src.ic50.cor[drug, "pearson"] <- cor(as.numeric(ccle_ic50_cc_src[drug, ]), as.numeric(gdsc_ic50_cc_src[drug, ]), method = "pearson", use = "pairwise.complete.obs")
  src.ic50.cor[drug, "spearman"] <- cor(as.numeric(ccle_ic50_cc_src[drug, ]), as.numeric(gdsc_ic50_cc_src[drug, ]), method = "spearman", use = "pairwise.complete.obs")
}

src.auc.cor = data.frame(cell_line = rownames(ccle_auc_cc_src), row.names = rownames(ccle_auc_cc_src))
src.auc.cor$pearson <- NA; src.auc.cor$spearman <- NA
for (drug in rownames(src.auc.cor)) {
  src.auc.cor[drug, "pearson"] <- cor(as.numeric(ccle_auc_cc_src[drug, ]), as.numeric(gdsc_auc_cc_src[drug, ]), method = "pearson", use = "pairwise.complete.obs")
  src.auc.cor[drug, "spearman"] <- cor(as.numeric(ccle_auc_cc_src[drug, ]), as.numeric(gdsc_auc_cc_src[drug, ]), method = "spearman", use = "pairwise.complete.obs")
}

# ==== Harrell's Concordance Index between CCLE and GDSC after sorting into PCLs ====
# MEK
# IC50
orig.mek.ic50 <- merge(orig.ic50.cor, mek.ic50.cor, by = 0, all = TRUE)
WriteXLS::WriteXLS(orig.mek.ic50, "~/capsule/results.mek_ic50_corr.xlsx", row.names = TRUE, col.names = TRUE, AdjWidth = TRUE, BoldHeaderRow = TRUE)
orig.mek.ic50.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.mek.ic50)
saveRDS(orig.mek.ic50.conc, "~/capsule/results/orig.mek.ic50.conc.rds")


# AUC
orig.mek.auc <- merge(orig.auc.cor, mek.auc.cor, by = 0, all = TRUE)
WriteXLS::WriteXLS(orig.mek.auc, "~/capsule/results.mek_auc_corr.xlsx", row.names = TRUE, col.names = TRUE, AdjWidth = TRUE, BoldHeaderRow = TRUE)
orig.mek.auc.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.mek.auc)
saveRDS(orig.mek.auc.conc, "~/capsule/results/orig.mek.auc.conc.rds")

# SRC
# IC50
orig.src.ic50 <- merge(orig.ic50.cor, src.ic50.cor, by = 0, all = TRUE)
WriteXLS::WriteXLS(orig.src.ic50, "~/capsule/results.src_ic50_corr.xlsx", row.names = TRUE, col.names = TRUE, AdjWidth = TRUE, BoldHeaderRow = TRUE)

orig.src.ic50.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.src.ic50)
saveRDS(orig.src.ic50.conc, "~/capsule/results/orig.src.ic50.conc.rds")

# AUC
orig.src.auc <- merge(orig.auc.cor, src.auc.cor, by = 0, all = TRUE)
WriteXLS::WriteXLS(orig.src.auc, "~/capsule/results.src_auc_corr.xlsx", row.names = TRUE, col.names = TRUE, AdjWidth = TRUE, BoldHeaderRow = TRUE)

orig.src.auc.conc <- survival::concordance(object = pearson.x ~ pearson.y, data = orig.src.auc)
saveRDS(orig.src.auc.conc, "~/capsule/results/orig.src.auc.conc.rds")
