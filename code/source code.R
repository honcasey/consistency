library(PharmacoGx); library(CePa); library(xml2)
downloadPSet("GDSC_2013"); downloadPSet("CCLE_2013")
load("~/capsule/code/PSets/CCLE_2013.RData"); load("~/capsule/code/PSets/GDSC_2013.RData")
intersected <- intersectPSet(c(CCLE, GDSC), intersectOn = c("drugs", "cell.lines"))
common_cells <- as.list(intersected[["CCLE"]]@cell)

drug_list <- as.list(intersected$CCLE@drug[["drugid"]])

# TO-DO: get list of all drugs
# TO-DO: sort into broad spectrum/targeted classes
# TO-DO: sort into more specific PCLs

PCL_list = list()
PCL_list$erlotinib <- CePa::read.gct("~/capsule/code/erlotinib_output.gct")
PCL_list$lapatinib <- CePa::read.gct("~/capsule/code/lapatinib.gct")
PCL_list$nilotinib <- CePa::read.gct("~/capsule/code/nilotinib.gct")
PCL_list$paclitaxel <- CePa::read.gct("~/capsule/code/paclitaxel.gct")
PCL_list$PD0325901 <- CePa::read.gct("~/capsule/code/PD-0325901.gct")
PCL_list$PD0332991 <- CePa::read.gct("~/capsule/code/PD-0332991.gct")
PCL_list$PHA665752 <- CePa::read.gct("~/capsule/code/PHA-665752.gct")
PCL_list$sorafenib <- CePa::read.gct("~/capsule/code/Sorafenib.gct")
PCL_list$PLX4720 <- CePa::read.gct("~/capsule/code/PLX-4720.gct")
PCL_list$saracatinib <- CePa::read.gct("~/capsule/code/saracatinib.gct")
PCL_list$selumetinib <- CePa::read.gct("~/capsule/code/selumetinib.gct")
PCL_list$crizotinib <- CePa::read.gct("~/capsule/code/crizotinib.gct")

# classes from PSet
pcl_count <- plyr::count(intersected[["CCLE"]]@drug$Mechanism.of.action) # there's a lot of classes here but each one only has 1 drug in it... probs don't use this
class_count <- plyr::count(intersected[["CCLE"]]@drug$Class) # there's only 3 classes and one of them just has 1 in it.. probs not gonna use this either

# classes from CMAP
druglist <- as.data.frame(read_csv("PSets/druglist330.csv"))
rownames(druglist) <- druglist$drug_name
pcl_count2 = plyr::count(druglist$`clue.io PCL`)
rownames(pcl_count2) <- pcl_count2$x

ccle_ic50 = data.frame(row.names = rownames(intersected[["CCLE"]]@sensitivity[["profiles"]]))
ccle_ic50$ic50_published <- intersected[["CCLE"]]@sensitivity[["profiles"]][, "ic50_published"]
ccle_ic50$drug <- NA
ccle_ic50$pcl <- NA

for (drug in rownames(ccle_ic50)) {
  which_drug <- stringr::str_split(drug, "_")[[1]][[2]]
  ccle_ic50[drug, "drug"] <- which_drug # get the drug name
  ccle_ic50[drug, "pcl"] <- druglist[which_drug, 'clue.io PCL'] # which pcl is this drug in
}


ccle_auc = data.frame(row.names = rownames(intersected[["CCLE"]]@sensitivity[["profiles"]]))
ccle_auc$auc_published <- intersected[["CCLE"]]@sensitivity[["profiles"]][, "auc_published"]
ccle_auc$drug <- NA
ccle_auc$pcl <- NA
for (drug in rownames(ccle_auc)) {
  which_drug <- stringr::str_split(drug, "_")[[1]][[2]]
  ccle_auc[drug, "drug"] <- which_drug # get the drug name
  ccle_auc[drug, "pcl"] <- druglist[which_drug, 'clue.io PCL'] # which pcl is this drug in
}


gdsc_ic50 = data.frame(row.names = rownames(intersected[["GDSC"]]@sensitivity[["profiles"]]))
gdsc_ic50$ic50_published <- intersected[["GDSC"]]@sensitivity[["profiles"]][, "ic50_published"]
gdsc_ic50$drug <- NA
gdsc_ic50$pcl <- NA
for (drug in rownames(gdsc_ic50)) {
  drugid <- stringr::str_split(drug, "_")[[1]][[2]]
  gdsc_drugid <- intersected[["GDSC"]]@drug[which(intersected[["GDSC"]]@drug$drugid == drugid), "drug.name"]
  unique_drugid <- intersected[["GDSC"]]@curation[["drug"]][which(intersected[["GDSC"]]@curation[["drug"]]$GDSC.drugid == gdsc_drugid), "unique.drugid"]# switch to unique drug id
  gdsc_ic50[drug, "drug"] <- unique_drugid # get the drug name
  gdsc_ic50[drug, "pcl"] <- druglist[unique_drugid, 'clue.io PCL'] # which pcl is this drug in
}
# TO-DO: manually fix drugs that are differently named in gdsc_ic50 & gdsc_auc

gdsc_auc = data.frame(row.names = rownames(intersected[["GDSC"]]@sensitivity[["profiles"]]))
gdsc_auc$auc_published <- intersected[["GDSC"]]@sensitivity[["profiles"]][, "auc_published"]
gdsc_auc$drug <- NA
gdsc_auc$pcl <- NA
for (drug in rownames(gdsc_auc)) {
  drugid <- stringr::str_split(drug, "_")[[1]][[2]]
  gdsc_drugid <- intersected[["GDSC"]]@drug[which(intersected[["GDSC"]]@drug$drugid == drugid), "drug.name"]
  unique_drugid <- intersected[["GDSC"]]@curation[["drug"]][which(intersected[["GDSC"]]@curation[["drug"]]$GDSC.drugid == gdsc_drugid), "unique.drugid"]# switch to unique drug id
  gdsc_auc[drug, "drug"] <- unique_drugid # get the drug name
  gdsc_auc[drug, "pcl"] <- druglist[unique_drugid, 'clue.io PCL'] # which pcl is this drug in
}





# could be used later
for (drug in rownames(druglist)) {
  for (cell in intersected[["CCLE"]]@cell[["cellid"]]) {
    rownm <- paste0("drugid_", drug, "_", cell) # which row to get in the sensitivity profiles
    drug_ic50 <- intersected[["CCLE"]]@sensitivity[["profiles"]][rownm, "ic50_published"] # get the drug's ic50
  }
}
