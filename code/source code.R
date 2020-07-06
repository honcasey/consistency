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


# 2. TO-DO: 

