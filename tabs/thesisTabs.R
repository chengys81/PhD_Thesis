#########################################
# tab1 primer list  
#########################################
library(Hmisc)
setwd("/Users/yongsheng/Documents/PhD/TeX/PPP2R5C CYS PhD/tabs/")

wb <- loadWorkbook("./tabs.xlsx")
sheets <- getSheets(wb) # generate a list of sheets
sheetNames <- names(sheets) # get names of each sheets
list <- list()
# both loops works
for (i in seq_along(sheetNames)) {
    list[[i]] <- read.xlsx("./tabs.xlsx", sheetName=sheetNames[i])
}

for (sheet in sheetNames) {
    list[[i]] <-read.xlsx("./tabs.xlsx", sheetName=sheet)
}

names(list) <- sheetNames

latex.list(list,longtable=TRUE,rowname=NULL)


primerList <- read.xlsx("./tabs.xlsx", sheetName = "primer")
colnames(primerList) <- c("Oligo No.", "Primer Sequence", "Name")

latex(primerList,longtable=TRUE,rowname=NULL)

instrumentList <- read.xlsx("./tabs.xlsx", sheetName = "instrument")
latex(instrumentList,longtable=TRUE,rowname=NULL)


antibodyList <- read.xlsx("./tabs.xlsx", sheetName = "antibody")
latex(antibodyList,longtable=TRUE,rowname=NULL)

chemicalList <- read.xlsx("./tabs.xlsx", sheetName = "chemical")
latex(chemicalList,longtable=TRUE,rowname=NULL)

kit <- read.xlsx("./tabs.xlsx", sheetName = "kit")
latex(kit,longtable=TRUE,rowname=NULL)

soft <- read.xlsx("./tabs.xlsx", sheetName = "software")
latex(soft,longtable=TRUE,rowname=NULL)

buffer <- read.xlsx("./tabs.xlsx", sheetName = "buffers")
latex(buffer,longtable=TRUE,rowname=NULL)

tab1.1 <- read.xlsx("./tabsChap1.xlsx", sheetName = "tab1.1")
latex(tab1.1,rowname=NULL)

tab1.2 <- read.xlsx("./tabsChap1.xlsx", sheetName = "tab1.2")
latex(tab1.2,rowname=NULL)


tab2.1 <- read.table("./refed_1.txt", sep="\t", skip=1, nrows=31,header=TRUE)
latex(tab2.1[,c(1,2,5)],rowname=NULL,booktabs=TRUE)

tab2.2 <- read.table("./fast_1.txt", sep="\t", skip=1, nrows=29,header=TRUE)
latex(tab2.2[,c(1,2,5)],rowname=NULL,booktabs=TRUE)

tab2.3 <- read.table("./random_1.txt", sep="\t", skip=1, nrows=24,header=TRUE)
latex(tab2.3[,c(1,2,5)],rowname=NULL,booktabs=TRUE)

tab2.4 <- read.table("./activated TFs.txt", sep="\t", skip=1, nrows=18,header=TRUE)
latex(tab2.4[,c(1,2,5)],rowname=NULL,booktabs=TRUE)
