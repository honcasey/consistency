library(PharmacoGx); library(dplyr); library(reshape2)
# downloadPSet("GDSC_2013"); downloadPSet("CCLE_2013")

# ==== loading PSets and intersecting CCLE and GDSC for common drugs and cell lines ====
load("~/capsule/code/PSets/CCLE_2013.RData"); load("~/capsule/code/PSets/GDSC_2013.RData")
intersected <- intersectPSet(c(CCLE, GDSC), intersectOn = c("drugs", "cell.lines"))
# ^ 15 drugs and 514 cell lines in common!
ccle <- intersected$CCLE
gdsc <- intersected$GDSC
rm(CCLE)
rm(GDSC)
#rm(intersected)

# classes from PSet
# pcl_count <- plyr::count(ccle@drug$Mechanism.of.action) # there's a lot of classes here but each one only has 1 drug in it... probs don't use this
# class_count <- plyr::count(ccle@drug$Class) # there's only 3 classes and one of them just has 1 in it.. probs not gonna use this either

# ==== loading manually curated list of drug information ====
druglist <- as.data.frame(read.csv("PSets/druglist330.csv")) 
# TO-DO: manually fix rownames of druglist so they match ccle@drug rownames (in excel), and make this csv actually presentable

druglist[] <- lapply(druglist, as.character) # change class of dataframe columns from vector to character
rownames(druglist) <- druglist$drug_name # rownames(ccle@drug)
# pcl_count2 = plyr::count(druglist$`clue.io PCL`)
# rownames(pcl_count2) <- pcl_count2$x

# ==== sort drugs into broad spectrum/targeted classes ====
classes = data.frame(drug = rownames(ccle@drug), row.names = rownames(ccle@drug))
classes$BroadSpectrum_or_Targeted <- NA

for (drug in rownames(druglist)) {
  if (druglist[drug, "target..drug.mechanism.from.ChEMBL."] == "single protein") {
    classes[drug, "BroadSpectrum_or_Targeted"] <- "targeted"
  }
  else if (druglist[drug, "target..drug.mechanism.from.ChEMBL."] == "protein family") {
    classes[drug, "BroadSpectrum_or_Targeted"] <- "broad-spectrum"
  }
  else if (is.na(druglist[drug, "target..drug.mechanism.from.ChEMBL."])) {
    classes[drug, "BroadSpectrum_or_Targeted"] <- NA
  }
}

# ==== sort drugs into clue.io defined PCLs ====
classes$PCL <- NA
for (drug in rownames(druglist)) {
  classes[drug, "PCL"] <- druglist[drug, "clue.io.PCL"]
}

ccle_ic50 <- as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "ic50_published"))
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

# TO-DO: simplify for-loops with lapply?
for (cell_line in rownames(ic50.cor)) {
  ic50.cor[cell_line, "pearson"] <- cor(ccle_ic50[, cell_line], gdsc_ic50[, cell_line], method = "pearson", use = "pairwise.complete.obs")
  ic50.cor[cell_line, "spearman"] <- cor(ccle_ic50[, cell_line], gdsc_ic50[, cell_line], method = "spearman", use = "pairwise.complete.obs")
}

# ==== correlation of published AUC between CCLE and GDSC ====
auc.cor = data.frame(cell_line = colnames(ccle_auc), row.names = colnames(ccle_auc))
auc.cor$pearson <- NA
auc.cor$spearman <- NA
for (cell_line in rownames(auc.cor)) {
  auc.cor[cell_line, "pearson"] <- cor(ccle_auc[, cell_line], gdsc_auc[, cell_line], method = "pearson", use = "pairwise.complete.obs")
  auc.cor[cell_line, "spearman"] <- cor(ccle_auc[, cell_line], gdsc_auc[, cell_line], method = "spearman", use = "pairwise.complete.obs")
}

# ==== sort drugs into broad-spectrum vs. targeted ====
broad_spectrum_drugs <- classes[grep("broad-spectrum", classes$BroadSpectrum_or_Targeted), "drug"] # drugs that are broad-spectrum
targeted_drugs <- classes[grep("targeted", classes$BroadSpectrum_or_Targeted), "drug"] # drugs that are targeted
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

# ==== remove inconsistent cell lines within GDSC ====
# inconsistent as in removing cell lines that had low correlation between measures (inconsistent defined as Spearman rank less than 0.5)

# ==== remove inconsistent cell lines within CCLE ====

# ==== recompute sens measures after sorting into PCLs ====
#how to organize all these objects into one?
auc.pcl.cor = data.frame(cell_line = colnames(ccle_auc), row.names = colnames(ccle_auc))
for (pcl in colnames(pcls)) {
  auc.pcl.cor[, paste0(pcl, ".pearson")] <- NA
  auc.pcl.cor[, paste0(pcl, ".spearman")] <- NA
}

for (pcl in pcl_list) {
  temp <- assign(pcl, data.frame(pcl_list[[pcl]])) # name of the PCL 
  ccle_auc_temp <- ccle_auc[pcl_list[[pcl]], ] # subset by the drugs in the PCL
  gdsc_auc_temp <- gdsc_auc[pcl_list[[pcl]], ]
  auc.temp.cor = data.frame(cell_line = colnames(ccle_auc), row.names = colnames(ccle_auc))
  auc.temp.cor$pearson <- NA
  auc.temp.cor$spearman <- NA
  for (cell_line in rownames(auc.temp.cor)) {
    auc.temp.cor[cell_line, "pearson"] <- cor(ccle_auc_temp[, cell_line], gdsc_auc_temp[, cell_line], method = "pearson", use = "pairwise.complete.obs")
    auc.temp.cor[cell_line, "spearman"] <- cor(ccle_auc_temp[, cell_line], gdsc_auc_temp[, cell_line], method = "spearman", use = "pairwise.complete.obs")
  }
}


# ==== remove inconsistent cell lines within GDSC ====

# ==== remove inconsistent cell lines within CCLE ====

# ==== Pearson correlation between CCLE and GDSC after binary sorting ====

# ==== Harrell's CI correlation between CCLE and GDSC after binary sorting ====

# ==== Pearson correlation between CCLE and GDSC after sorting into PCLs ====

# ==== Harrell's CI correlation between CCLE and GDSC after sorting into PCLs ====
