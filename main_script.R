# LOAD PACKAGES AND FUNCTIONS --------------------------------------------------

source("./codes/my_packages.R")
source("./codes/my_functions.R")

# DIFFERENTIAL EXPRESSION ANALYSIS PARAMETERS ----------------------------------

## GSE25001_Whole
gse <- "GSE25001"
gset <- getGEO(gse, GSEMatrix =TRUE, getGPL = F, AnnotGPL = F)
gset <- gset[[1]]
show(gset)
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
LogC

### Phases: 53 early (0-3 day of illness), 100 late (4-8d), 34 convalescent (follow up)
gsmsP <- paste0("2012XX20101X22211011002222222011X0111011101111X211",
               "10011121110101010211210X1121X0111012101X1X2011101X",
               "01011X211011X1222202010XX21211X011101010X1110X1101",
               "0120121X11010010XX0011011011011X01021X011211011111",
               "201020101")
groupsP <- make.names(c("Early","Late","Convalescent"))

### Severity: 122 nonDSS (acute and convalescent), 65 DSS (acute and convalescent)
gsmsS <- paste0("0001XX10011X10111111110000001001X0001001100000X100",
               "00001000001111000100100X0000X0000111001X1X1000000X",
               "00001X100001X1000010110XX11111X111000000X0000X0000",
               "0011110X00001100XX0100111000000X00110X001100000000",
               "000110011")
groupsS <- make.names(c("nonDSS","DSS"))

## GSE28405_Genome_Tolfvenstam
gse <- "GSE28405"
gset <- getGEO(gse, GSEMatrix =TRUE, getGPL = T, AnnotGPL = T)
gset <- gset[[1]]
show(gset)
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
LogC

### Phases: 31 early (less 72h days after symptoms), 31 late (3-7d), 
# 31 convalescent (2-5wk), 26 control (less 72h days after symptoms)
gsmsP <- paste0("33333333333333333333333333000000000000000000000000",
               "00000001111111111111111111111111111111222222222222",
               "2222222222222222222")
groupsP <- make.names(c("Early","Late","Convalescent","Control"))
### Severity: 93 DF, 26 Control
gsmsS <- paste0("11111111111111111111111111000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "0000000000000000000")
groupsS <- make.names(c("DF","Control"))

## GSE28991_Genome_11
gse <- "GSE28991"
gset <- getGEO(gse, GSEMatrix =TRUE, getGPL = F, AnnotGPL = F)
gset <- gset[[1]]
show(gset)
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
LogC
### Phases: 11 early (0-3 days after symptoms), 11 late (4-7d), 11 convalescent (15-25d)
gsmsP <- "012012012012012012012120120012012"
groupsP <- make.names(c("Early","Late","Convalescent"))
gsmsS <- NULL
groupsS <- NULL

## GSE28988_Global
gse <- "GSE28988"
gset <- getGEO(gse, GSEMatrix =TRUE, getGPL = F, AnnotGPL = F)
gset <- gset[[1]]
show(gset)
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
LogC

### Phases: 114 early (less 72h days after symptoms), 96 late (3-7d), 96 convalescent (2-5wk)
gsmsP <- paste0("01201201201200120120120120101201201201201201201201",
               "20120120120012020201200120012010120120120201200201",
               "01201201201201201201201201201200201201201201201012",
               "01201201201201201200010001202012012012012012012012",
               "01201201201201201201201201201201201201201201201201",
               "20120120120120120120120120102012012012010120120120",
               "120120")
groupsP <- make.names(c("Early","Late","Convalescent"))
gsmsS <- NULL
groupsS <- NULL

## GSE40628_Dengue
gse <- "GSE40628"
gset <- getGEO(gse, GSEMatrix =TRUE, getGPL = F, AnnotGPL = F)
gset <- gset[[1]]
show(gset)
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
LogC

### Phases: 1 early (0-3d), 17 Late (4-8d), 10 Convalescent (15+d), 4 healthy control 
gsmsP <- "X111212110312X11123111212122133122"
groupsP <- make.names(c("Early","Late","Convalescent","Control"))
### Severity: 17 nonDSS (WHO stage 1-2 or indet), 13 DSS (WHO stage 3-4), 4 healthy control
gsmsS <- "0000110110201111012000010001022101"
groupsS <- make.names(c("NonDSS","DSS","Control"))

## GSE43777_Sequential
gse <- "GSE43777"
### GPL570 	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
gset <- getGEO(gse, GSEMatrix =TRUE, getGPL = F, AnnotGPL = F)
gset <- gset[[2]]
show(gset)
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
LogC

### Phases: 26 early (0-3d), 51 late (4-10d), 24 convalescent
gsmsP <- paste0("00211221111121110121011121112200112111111121201112",
               "01120011200120001200112001120012011201120112112001",
               "2")
groupsP <- make.names(c("Early","Late","Convalescent"))
### Severity: 52 DF, 49 DHF
gsmsS <- paste0("00011100001110000001000000000000000000011110011111",
               "11111111100000000011111000001111111111111111111000",
               "0")
groupsS <- make.names(c("DF","DHF"))

## GSE51808_Systems
gse <- "GSE51808"
gset <- getGEO(gse, GSEMatrix =TRUE, getGPL = F, AnnotGPL = F)
gset <- gset[[1]]
show(gset)
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
LogC

### Phases: 28 Acute, 19 Convalescent, 9 Control
gsmsP <- "00000000000000000000000000001111111111111111111222222222"
groupsP <- make.names(c("Acute","Convalescent","Control"))
### Severity: 31 DF, 16 DHF, 9 Control
gsmsS <- "00001010010111000010110010000000100011001011000222222222"
groupsS <- make.names(c("DF","DHF","Control"))

# DATA FILTERING ---------------------------------------------------------------

# list all tables
original_data <- list.files(path = "./data", pattern = "_vs_", 
                            full.names = TRUE, recursive = TRUE)
# example: original_name <- "./data/GSE25001_Whole/GSE25001.100_4to8d_vs_34_convalescent.tsv"

# apply the DEG_filter function to each table in list
lapply(original_data, DEG_filter)

# DE ANALYSIS METADATA TABLE ---------------------------------------------------

# list paths
outputs1 <- list.dirs(path = "./output", full.names = T, recursive = F)
GSEs <- grep(x = basename(outputs1), pattern = "^GSE")
outputs2 <- outputs1[GSEs]
folders <- list.dirs(path = outputs2, full.names = T, recursive = F)

# turn to df and manipulate columns
outputs_table <- as.data.frame(folders)
outputs_table <- dplyr::mutate(outputs_table, folders = (gsub(x = outputs_table$folders, pattern = "./output/", "")))
outputs_table <- dplyr::mutate(outputs_table, folder1 = (dirname(outputs_table$folders)))
outputs_table <- dplyr::mutate(outputs_table, GSE = (gsub(x = outputs_table$folder1, pattern = "_.*", "")))
outputs_table <- dplyr::mutate(outputs_table, studytitle = (gsub(x = outputs_table$folder1, pattern = "^[[:alnum:]]+_", "")))
outputs_table <- dplyr::mutate(outputs_table, folder2 = (basename(outputs_table$folders)))
outputs_table <- dplyr::mutate(outputs_table, temp1 = (gsub(x = outputs_table$folder2, pattern = "^[[:alnum:]]+\\.", "")))
outputs_table <- dplyr::mutate(outputs_table, g1N = (gsub(x = outputs_table$temp1, pattern = "_.*", "")))
outputs_table <- dplyr::mutate(outputs_table, g1cond = (gsub(x = outputs_table$temp1, pattern = "_vs.*", "")))
outputs_table <- dplyr::mutate(outputs_table, g1cond = (gsub(x = outputs_table$g1cond, pattern = "^[[:alnum:]]+_", "")))
outputs_table <- dplyr::mutate(outputs_table, temp2 = (gsub(x = outputs_table$temp1, pattern = ".*vs_", "")))
outputs_table <- dplyr::mutate(outputs_table, g2N = (gsub(x = outputs_table$temp2, pattern = "_.*", "")))
outputs_table <- dplyr::mutate(outputs_table, g2cond = (sub(x = outputs_table$temp2, pattern = "^[[:alnum:]]+_", "")))

exclude <- grep(colnames(outputs_table), pattern = "temp.")
outputs_table <- outputs_table[,-exclude]

# Number of DEGs (total, up, down) found in each comparison
for (a in c("UP","DOWN","DEG")) {
  outputs_files <- list.files(path = "./output", pattern = paste0(a, "_.*\\.tsv"), full.names = T, recursive = T)
  rownumbers <- lapply(outputs_files, function(i){
    nrow(read.table(i, header = TRUE, sep = "\t", dec = ".", quote = "")) })
  names(rownumbers) <- gsub(x = outputs_files, pattern = paste0("./output/|/", a, "_.*"),replacement = "")
  row_df <- melt(rownumbers)
  colnames(row_df) <- c(a, "folders")
  outputs_table <- merge(x = outputs_table, y = row_df, by = "folders", all.x = TRUE)
}

write.table(outputs_table, file = "clipboard", sep="\t", dec=".", 
            quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

# Add groups' time and severity classification manually.
# Saved on './output/tabela_comparisons.txt'

# DEGs PLOT --------------------------------------------------------------------

# Number of DEGs found by comparison, grouped by GSE.

# read comparisons table, turn number of downregulated genes negative,
# and turn df into long format
COMPARISONS <- read.table(file = "./output/tabela_comparisons.txt", 
                          header = TRUE, sep = "\t", dec = ".", quote = "")
my_table1 <- dplyr::mutate(COMPARISONS, DOWNNeg = -DOWN)
my_table2 <- melt.data.frame(my_table1, measure.vars = c("UP", "DOWNNeg"))

# diverging bars plot:
ggplot(my_table2, aes( y=value, 
                       x=gsub(x = folder2, pattern = ".+\\.", replacement = "")))+
  geom_bar(stat='identity', aes(fill=variable), width = 0.8)+
  scale_fill_brewer(name = element_blank(), palette = "Set1",
                    labels = c("up regulated", "down regulated")) +
  labs(x = element_blank(), y = "DEGs") +
  coord_flip() +
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        panel.background = element_rect( fill = "white", color = "grey"),
        panel.grid = element_blank(),
        strip.text.y = element_text(vjust = 1, size = 8), 
        strip.background = element_rect(colour="grey", fill="white")) +
  facet_grid(rows = vars(GSE), scale = "free_y", space = "free_y")

#ggsave(filename = "DEGs.png", path = "./output/figures", width = 2500, height = 3300, 
#       units = "px", dpi = 300, device = "png", scale = 1)

# INTERFEROME ------------------------------------------------------------------

## SEARCH HELPER ---------------------------------------------------------------

# gets the gene list to clipboard. 
# Paste it in the interferome and download files manually.
# moves downloaded files from interferome to folders in output.
# repeat for all comparisons.

# list paths
outputs1 <- list.dirs(path = "./output", full.names = T, recursive = F)
GSEs <- grep(x = basename(outputs1), pattern = "^GSE")
outputs2 <- outputs1[GSEs]
folders <- list.dirs(path = outputs2, full.names = T, recursive = F)

# run the next line only once before starting
i <- 1

### UP 

# UP gene symbol list to clipboard:
print(c(i, folders[i]))
my_file <- list.files(path = folders[i], pattern = "UP_.+\\.tsv", full.names = TRUE, recursive = F)
table <- read.table(file = my_file, header = TRUE, sep = "\t", dec = ".", quote = "")
write.table(table$Gene.symbol, file = "clipboard", 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = FALSE)

# submit to interferome search, download files

# UP downloaded files moving:
# if downloaded files exist, create UP folder, move downloaded files to UP folder
dl_files <- list.files(path = "C:/Users/julia/Downloads", pattern = ".+SearchResults.txt",
                       full.names = TRUE, recursive = FALSE)
if(length(dl_files) != 0) {
  new_path_interferome <- paste0(folders[i], "/Interferome_hs_UP")
  fs::dir_create(path = new_path_interferome)
  fs::file_move(new_path = new_path_interferome, path = dl_files)}

### DOWN

# DOWN gene symbol list to clipboard:
print(c(i, folders[i]))
my_file <- list.files(path = folders[i], pattern = "DOWN_.+\\.tsv",
                      full.names = TRUE, recursive = F)
table <- read.table(file = my_file, header = TRUE, sep = "\t", dec = ".", quote = "")
write.table(table$Gene.symbol, file = "clipboard", 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = FALSE)

# submit to interferome search, download files

# DOWN downloaded files moving:
# if downloaded files exist, create DOWN folder, move downloaded files to DOWN folder
dl_files <- list.files(path = "C:/Users/julia/Downloads", pattern = ".+SearchResults.txt",
                       full.names = TRUE, recursive = FALSE)
if(length(dl_files) != 0) {
  new_path_interferome <- paste0(folders[i], "/Interferome_hs_DOWN")
  fs::dir_create(path = new_path_interferome)
  fs::file_move(new_path = new_path_interferome, path = dl_files)}
i <- i + 1
# go back to ### UP and repeat with all folders

## RESULTS SUMMARY TABLE -------------------------------------------------------

# Summarizes the data downloaded from the interferome
# from the data search results (experimental data), as it is more accurate.

# list paths
outputs1 <- list.dirs(path = "./output", full.names = T, recursive = F)
GSEs <- grep(x = basename(outputs1), pattern = "^GSE")
outputs2 <- outputs1[GSEs]
folders <- list.dirs(path = outputs2, full.names = T, recursive = F)

# load the comparisons table
COMPARISONS <- read.table(file = "./output/tabela_comparisons.txt", 
                          header = T, sep = "\t", dec = ".", quote = "")
comparisons_table <- COMPARISONS
  
# list all 'data search results' files
DataSearchRes <- list.files(path = paste0("./output"), pattern = 
                              ".+_DataSearchResults\\.txt", full.names = T, recursive = T)

# get number of unique genes from 'Data search results'
DataN <- lapply(DataSearchRes, function(i){
  table <- read.table(file = i, header = TRUE, sep = "\t", dec = ".", quote = "", skip = 19)
  table1 <- table[!duplicated(table$Gene.Name), ]
  nrow(table1)})
# add file names and manipulate
names(DataN) <- gsub(x = DataSearchRes, pattern = "/[[:alnum:]]+_[[:alnum:]]+\\.txt|\\./output/", replacement = "") 
DataNdf <- melt(DataN)
DataNdf <- dplyr::mutate(DataNdf, folder1 = (gsub(x = DataNdf$L1, pattern = "/.*", "")))
DataNdf <- dplyr::mutate(DataNdf, temp1 = (gsub(x = DataNdf$L1, pattern = "^[[:alnum:]_]+/", "")))
DataNdf <- dplyr::mutate(DataNdf, folder2 = (gsub(x = DataNdf$temp1, pattern = "/.*", "")))
DataNdf <- dplyr::mutate(DataNdf, UPDOWN = (gsub(x = DataNdf$L1, pattern = ".*_", "")))
exclude <- grep(colnames(DataNdf), pattern = "temp.")
DataNdf <- DataNdf[,-exclude]


# get number of UP genes and merge tables
UPDataN <- subset(x = DataNdf, DataNdf$UPDOWN == "UP")
comparisons_table <- merge(x = comparisons_table, y = UPDataN[,c(1,4)], by = "folder2", all.x = TRUE)
colnames(comparisons_table)[colnames(comparisons_table)== "value"] <- "UP_interferome"
# get number of DOWN genes and merge tables
DOWNDataN <- subset(x = DataNdf, DataNdf$UPDOWN == "DOWN")
comparisons_table <- merge(x = comparisons_table, y = DOWNDataN[,c(1,4)], by = "folder2", all.x = TRUE)
colnames(comparisons_table)[colnames(comparisons_table)== "value"] <- "DOWN_interferome"

### number of unique genes by IFN type from 'data search results"

# for each IFN type,
for (IFN_type in c("I","II","III")) {
# for each table
  # subset by select IFN type; remove duplicates; get number of genes;
  DataNType <- lapply(DataSearchRes, function(i){
    table <- read.table(file = i, header = TRUE, sep = "\t", dec = ".", quote = "", skip = 19)
    table1 <- subset(x = table, table$Inteferome.Type == IFN_type)
    table2 <- table1[!duplicated(table1$Gene.Name), ]
    nrow(table2)})
# turn to df, manipulate names
  names(DataNType) <- gsub(x = dirname(DataSearchRes), pattern = "\\./output/", replacement = "") 
  DataNTypedf <- melt(DataNType)
  DataNTypedf <- dplyr::mutate(DataNTypedf, folder1 = (gsub(x = DataNTypedf$L1, pattern = "/.*", "")))
  DataNTypedf <- dplyr::mutate(DataNTypedf, temp1 = (gsub(x = DataNTypedf$L1, pattern = "^[[:alnum:]_]+/", "")))
  DataNTypedf <- dplyr::mutate(DataNTypedf, folder2 = (gsub(x = DataNTypedf$temp1, pattern = "/.*", "")))
  DataNTypedf <- dplyr::mutate(DataNTypedf, UPDOWN = (gsub(x = DataNTypedf$L1, pattern = ".*_", "")))
  exclude <- grep(colnames(DataNTypedf), pattern = "temp.")
  DataNTypedf <- DataNTypedf[,-exclude]

# get number of UP unique genes by selected IFN type
  UPDataNType <- subset(x= DataNTypedf, DataNTypedf$UPDOWN == "UP")
  comparisons_table <- merge(x=comparisons_table, y =UPDataNType[,c(1,4)], by = "folder2", all.x = TRUE)
  colnames(comparisons_table)[colnames(comparisons_table)== "value"] <- paste0("UP_Type_",IFN_type)

# get number of DOWN unique genes by selected IFN type
  DOWNDataNType <- subset(x= DataNTypedf, DataNTypedf$UPDOWN == "DOWN")
  comparisons_table <- merge(x=comparisons_table, y =DOWNDataNType[,c(1,4)], by = "folder2", all.x = TRUE)
  colnames(comparisons_table)[colnames(comparisons_table)== "value"] <- paste0("DOWN_Type_",IFN_type)
}

# arrange columns order
# turn NAs to 0
comparisons_table <- comparisons_table %>% dplyr::relocate(DEG, UP, .after = g2sev)
comparisons_table <- comparisons_table %>% dplyr::relocate(DOWN, DOWN_interferome, 
                     DOWN_Type_I, DOWN_Type_II ,DOWN_Type_III, .after = UP_Type_III)
comparisons_table <- comparisons_table %>% dplyr::relocate(folder2, .after = studytitle)
comparisons_table[is.na(comparisons_table)] <- 0

# write table with interferome results summarized
write.table(comparisons_table, file = "./output/tabela_interferome_datasearchres.tsv", 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)

#### PLOT ----------------------------------------------------------------------

# diverging and overlapped bars plot that represents the summary table of the DEGs
# and IFN-regulated DEGs, by comparisons used in further analyses
# check rascunho for other versions

# get the interferome results summary table
tabela <- read.table(file = "./output/tabela_interferome_datasearchres.tsv", 
                     header = T, sep = "\t", dec = ".", quote = "")

# get the comparisons used in the heatmaps, and rbind them
  tabelaTime <- subset(tabela, 
                       ((tabela$g1time == "early_acute" | tabela$g1time == "late_acute")
                        & tabela$g2time == "convalescent"))
  tabelaTime <- dplyr::mutate(tabelaTime, 
                              Names = (paste0(g1N, " ", 
                                              gsub(g1time, pattern = "_", replacement = " "), 
                                              " vs. ", g2N, " ", g2time, " (", GSE, ")")))
  tabelaTime$Classif <- "Phase"
  tabelaTime <- dplyr::relocate(tabelaTime, Names, Classif, .after = folder2)
  
  tabela1 <- tabelaTime
  tabelaUP <- tabela1[,c(1:21)]
  tabelaUP$UPDOWN <- "UP"
  colnames(tabelaUP)[colnames(tabelaUP)%in%(colnames(tabelaUP)[16:21])] <- c("Total_DEGs", "DEGs", "IFNreg_DEGs", "Type_I", "Type_II", "Type_III")
  tabelaDOWN <- tabela1[,-c(17:21)]
  tabelaDOWN$UPDOWN <- "DOWN"
  colnames(tabelaDOWN)[colnames(tabelaDOWN)%in% (colnames(tabelaDOWN)[16:21])] <- c("Total_DEGs", "DEGs", "IFNreg_DEGs", "Type_I", "Type_II", "Type_III")
  tabela2 <- rbind.data.frame(tabelaUP, tabelaDOWN)
  tabela4 <- tidyr::pivot_longer(data = tabela2[,-(19:21)], cols = c("DEGs", "IFNreg_DEGs"),names_to = "Type", values_to = "Number")
  tabelaTime <- tabela4
  
  orderTime <- c("early_acute", "late_acute")
  tabelaTime <- arrange(transform(tabelaTime, g1time=factor(g1time,levels=orderTime)),g1time,GSE)
  
  tabelaSev <- subset(tabela, 
                      tabela$g1time == "all_acute" & 
                        (tabela$g1sev == "DF" | tabela$g1sev == "DHF" | tabela$g1sev == "DSS")
                      & tabela$g2time == "convalescent"
                      & !(tabela$g2sev == "DF" | tabela$g2sev == "DHF" | tabela$g2sev == "DSS"))
  tabelaSev <- dplyr::mutate(tabelaSev, 
                             Names = (paste0(g1N," ", g1sev, 
                                             " acute vs. ", g2N, " ", g2time, " (", GSE, ")")))
  tabelaSev$Classif <- "Severity"
  tabelaSev <- dplyr::relocate(tabelaSev, Names, Classif, .after = folder2)
  
  tabela1 <- tabelaSev
  tabelaUP <- tabela1[,c(1:21)]
  tabelaUP$UPDOWN <- "UP"
  colnames(tabelaUP)[colnames(tabelaUP)%in%(colnames(tabelaUP)[16:21])] <- c("Total_DEGs", "DEGs", "IFNreg_DEGs", "Type_I", "Type_II", "Type_III")
  tabelaDOWN <- tabela1[,-c(17:21)]
  tabelaDOWN$UPDOWN <- "DOWN"
  colnames(tabelaDOWN)[colnames(tabelaDOWN)%in% (colnames(tabelaDOWN)[16:21])] <- c("Total_DEGs", "DEGs", "IFNreg_DEGs", "Type_I", "Type_II", "Type_III")
  tabela2 <- rbind.data.frame(tabelaUP, tabelaDOWN)
  tabela4 <- tidyr::pivot_longer(data = tabela2[,-(19:21)], cols = c("DEGs", "IFNreg_DEGs"),names_to = "Type", values_to = "Number")
  tabelaSev <- tabela4
  
  orderSev <- c("DF", "DHF", "DSS")
  tabelaSev <-arrange(transform(tabelaSev, g1sev=factor(g1sev,levels=orderSev)),g1sev,GSE)

  tabela16 <- rbind.data.frame(tabelaTime, tabelaSev)
  neworder <- rev(unique(tabela16$Names))
  tabela16 <- arrange(transform(tabela16, Names=factor(Names,levels=neworder)),Names)
  
  
## actual plot:
ggplot(data = tabela16) +
  #downregulated total DEGs and downregulated IFN-regulated DEGs, and labels
  geom_bar(data = subset(tabela16, (UPDOWN == "DOWN")), 
           position = position_dodge(), stat = "identity", 
           aes(fill = paste(UPDOWN, Type), y = -Number, x = Names),
           alpha = 1, width = 0.8)+
  geom_text(data = subset(tabela16, (UPDOWN == "DOWN")),
            aes(y = -Number, x = Names, label = Number, hjust = 1.1, 
                vjust =  ifelse(test = (Type == "IFNreg_DEGs"), yes = -0.4, no = 1.3)),
            size = 3) +
  #upregulated total DEGs and upregulated IFN-regulated DEGs, and labels
  geom_bar(data = subset(tabela16, (UPDOWN == "UP")), 
           position = position_dodge(), stat = "identity", 
           aes(fill = paste(UPDOWN, Type), y = Number, x = Names), 
           alpha = 1, width = 0.8)+
  geom_text(data = subset(tabela16, (UPDOWN == "UP")),
            aes(y = Number, x = Names, label = Number, hjust = -0.1, 
                vjust =  ifelse(test = (Type == "IFNreg_DEGs"), yes = -0.4, no = 1.3)),
            size = 3) +
  theme(legend.position = "bottom", 
        legend.title =  element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(), 
#        axis.text.y = element_text(hjust = 0.5),
        panel.background = element_rect( fill = "white", color = "grey"),
        panel.grid = element_blank(),
        strip.text.y = element_text(vjust = 1, hjust = 0.5), 
        strip.background = element_rect(colour="grey", fill="gray95")) +
  labs(x = element_blank(), y = element_blank()) +
  coord_flip() +
  facet_grid(rows = vars(Classif), scale = "free_y", space = "free_y") +
  scale_fill_manual(values=c("#5782bb","#a2bad9", "#db443c", "#eb9793"),
                    labels = c("Down-regulated DEGs", "Down-regulated interferon-regulated DEGs", 
                               "Up-regulated DEGs", "Up-regulated interferon-regulated DEGs")) +
  guides(fill = guide_legend(nrow = 2, byrow = F)) 

#ggsave(filename = "IFN_DEGs_classif6.png", path = "./output/figures", width = 2700, height = 1800, 
#       units = "px", dpi = 300, device = "png", scale = 1)

#ggsave(filename = "IFN_DEGs_classif2.svg", path = "./output/figures", width = 2900, height = 1800, 
#       units = "px", dpi = 300, device = "svg", scale = 1)

# GENE LISTS BY IFN TYPE -------------------------------------------------------

# for each IFN type:
  # gets gene lists results from the interferome, separated by:
  # up or downregulated genes
# adds them to the comparisons table.

# list all 'data search results' files
DataSearchRes <- list.files(path = "./output", pattern = 
                 ".+_DataSearchResults\\.txt", full.names = T, recursive = T)
# get the interferome results summary table
tabela <- read.table(file = "./output/tabela_interferome_datasearchres.tsv", 
                     header = T, sep = "\t", dec = ".", quote = "")


# 'melt' UP and DOWN, then IFN type
tabelaUP <- tabela[,c(1:15, 17:19)]
tabelaUP$UPDOWN <- "UP"
colnames(tabelaUP)[colnames(tabelaUP)%in%(colnames(tabelaUP)[14:18])] <- c("Total_DEGs", "DEGs", "I", "II", "III")

tabelaDOWN <- tabela[,-c(15:19, 21)]
tabelaDOWN$UPDOWN <- "DOWN"
colnames(tabelaDOWN)[colnames(tabelaDOWN)%in% (colnames(tabelaDOWN)[14:18])] <- c("Total_DEGs", "DEGs", "I", "II", "III")

tabela1 <- rbind.data.frame(tabelaUP, tabelaDOWN)
tabela2 <- tidyr::pivot_longer(data = tabela1, cols = c("I", "II", "III"),names_to = "IFN_TYPE", values_to = "Type_N" )


# for each IFN type, for each comparison, for UP and DOWN, get list of unique genes
TypeGenesList <- lapply(c("I","II","III"), genes_by_type)

# rbind the gene list df's of each type
TypeGenesAll <- rbind.fill(TypeGenesList)
# merge the gene lists with the comparisons table
tabela3 <- merge(x=tabela2, y =TypeGenesAll[,-(1:3)], by = c("folder2", "UPDOWN", "IFN_TYPE"), all.x = TRUE)
#arrange columns order
tabela3 <- tabela3 %>% dplyr::relocate(folder2, .after = studytitle)
tabela3 <- tabela3 %>% dplyr::relocate(UPDOWN, .after = DEGs)
tabela3 <- tabela3 %>% dplyr::relocate(IFN_TYPE, .after = UPDOWN)

write.table(tabela3, file = "./output/tabela_genelists_IFNtype.tsv", 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

# circos tables ----------------------------------------------------------------
tabcircos <- read.table(file = "./degs_ifn_late.txt", header = F, sep = "\t", dec = ".")
nomes1 <- as.character(tabcircos[1,])
nomes2 <- as.character(tabcircos[2,])
nomes3 <- paste(nomes1, nomes2)
nomes4 <- gsub(nomes3, pattern = " ", replacement = "_")

E_circos <- list()
for (i in 1:ncol(tabcircos)) {
  list1 <- tabcircos[-c(1,2),i]
  empties <- which(list1 == "")
  if (is_empty(empties)== F) {
  list2 <-list1[-empties]
  } else {
    list2 <- list1
  }
  E_circos[[i]] <- list2
}
names(E_circos) <- nomes4

E_circos_matrix <- ldply(E_circos, function(i){
  teste <- lapply(E_circos, intersect, i)
  as.numeric(lapply(teste, length))
})
colnames(E_circos_matrix) <- c("data", nomes4)


write.table(E_circos_matrix, file = "./output/circos_tables/Circos_matrix_late.txt", 
            quote = T, sep = "\t", row.names = F, dec = ".", col.names = T)

GSE_merged <- list()
for (GSE in c("GSE25001", "GSE28405", "GSE28988", "GSE28991", "GSE43777")) {
  GSE_table <- tabcircos[-c(1,2),which(tabcircos[1,] == GSE)]
  GSE_I <- GSE_table[,1]
  GSE_II <- GSE_table[,2]
  GSE_III <- GSE_table[,3]
  GSE_intersec <- c(GSE_I, GSE_II, GSE_III)
  GSE_intersec <- GSE_intersec[!duplicated(GSE_intersec)]
  GSE_intersec <- GSE_intersec[-which(GSE_intersec == "")]
  pos <- grep(c("GSE25001", "GSE28405", "GSE28988", "GSE28991", "GSE43777"), pattern = GSE)
  GSE_merged[[pos]] <- GSE_intersec
  names(GSE_merged)[pos] <- GSE
}

E_circos_matrix_GSE <- ldply(E_circos, function(i){
  teste <- lapply(GSE_merged, intersect, i)
  as.numeric(lapply(teste, length))
})

colnames(E_circos_matrix_GSE) <- c("data", names(GSE_merged))

write.table(E_circos_matrix_GSE, file = "./output/circos_tables/Circos_matrix_GSE_late.txt", 
            quote = T, sep = "\t", row.names = F, dec = ".", col.names = T)



E_circos_Type <- E_circos[grep(names(E_circos), pattern = "_I$")]
E_Type_1x1_intersection <- list()
for (i in 1:length(E_circos_Type)) {
  for (j in 1:length(E_circos_Type)){
    if (names(E_circos_Type[i]) != names(E_circos_Type[j])) {
      E_Type_1x1_intersection[[i]] <- intersect(E_circos_Type[[i]],E_circos_Type[[j]])
    }
  }
}

  lapply(E_circos_Type, function(i){
  teste <- lapply(E_circos_Type, function(j){
    if (i != j) {intersect(i,j)}
    })
})
which

E_Type_1x1_intersection <- lapply(E_circos_Type, function(i){
  teste <- lapply(E_circos_Type, intersect, i)
})
E_Type_1x1_intersection_merged <- list()
for (i in 1:length(E_Type_1x1_intersection)) {
  same <- which(names(E_Type_1x1_intersection)[i] == names(E_Type_1x1_intersection[[i]]))
  E_Type_1x1_intersection[[i]] <- E_Type_1x1_intersection[[i]][-same]
  merged <- c(E_Type_1x1_intersection[[i]][1],E_Type_1x1_intersection[[i]][2],
    E_Type_1x1_intersection[[i]][3],E_Type_1x1_intersection[[i]][4])
  merged_unique <- unique(merged)
  E_Type_1x1_intersection_merged[[i]] <- merged_unique
    }

list_intersec <- Reduce(intersect, genelists)
interseccao <- as.data.frame(interseccao)
colnames(interseccao) <- "V1"

names(E_Type_1x1_intersection[1])
names(E_Type_1x1_intersection[[1]])


## upset tables ----------------------------------------------------------------

    for (EL in c("early_acute","late_acute")) {
    for (IFN_Type in c("IFN Type I","IFN Type II","IFN Type III")) {
      tabupset <- read.table(file = "./upset_table_phases.txt", header = F, sep = "\t", dec = ".")
      tabupsetType <- tabupset[which(tabupset[3,] == IFN_Type)]
      tabupset <- tabupsetType[which(tabupsetType[2,] == EL)]
      
      nomes1 <- as.character(tabupset[1,])
      nomes2 <- as.character(tabupset[2,])
      nomes21 <- gsub(nomes2, pattern = "_acute", replacement = "")
      nomes3 <- as.character(tabupset[3,])
      nomes31 <- gsub(nomes3, pattern = " ", replacement = "_")
      nomes4 <- paste(nomes1, nomes21, nomes31, sep= "_")
      
      tabupset2 <- tabupset[-(1:4),]
      colnames(tabupset2) <- nomes4
      colnames(tabupset2)
      name <- paste0("./output/upset_tables/upset_table_", 
                     gsub(x = EL, pattern = "acute", replacement = ""),
                     gsub(x = IFN_Type, pattern = "IFN Type ", replacement = ""),
                     ".txt")
      name
      write.table(tabupset2, file = name, 
                  quote = T, sep = "\t", row.names = F, dec = ".", col.names = T)
    }
  }

# HEATMAP BY IFN TYPE ----------------------------------------------------------

## GENE INTERSECTIONS ----------------------------------------------------------

# Get intersections of genes in groups of interest

# load the comparisons table
COMPARISONS <- read.table(file = "./output/tabela_comparisons.txt", 
                          header = T, sep = "\t", dec = ".", quote = "")

#set name of the group of comparisons (for the output tables' names)
group_name <- "DF_vs_conv"

# filter by desired group of comparisons (specify desired groups conditions)
tabela_filtered <- filter(COMPARISONS, COMPARISONS$g1sev == "DF" &
                            COMPARISONS$g2sev != "DF" &
                            COMPARISONS$g2time == "convalescent")
tabela_filtered <- filter(COMPARISONS, COMPARISONS$g1time == "early_acute" & 
                            COMPARISONS$g2time == "convalescent")
compar_names <- tabela_filtered$folders

# get the intersection of the genes and write 3 tables, one for each IFN type
genes_intersection(compar_names = compar_names, group_name = group_name)

### Heatmap Matrix -------------------------------------------------------------

# For selected comparison groups, get the intersections, join intersection genes,
# and get logFC from tables downloaded from GEO2R (before filtering)

# list all intersection tables
intersec_tables_list <- list.files(path = "./output/Intersection_Tables/", pattern = "\\.tsv", 
                                   full.names = TRUE, recursive = F)
intersec_tables <- lapply(intersec_tables_list, function(i){
  read.table(file = i, header = TRUE, sep="\t", dec=".", quote = "")})
names(intersec_tables) <- gsub(x = intersec_tables_list, pattern = ".+/|\\.tsv", replacement = "") 

# set desired intersections and final table name
group_names <- c("DF_vs_conv", "DHF_vs_conv", "DSS_vs_conv")
group_names <- c("early_vs_conv", "late_vs_conv")
table_name <- "merged_sev"

# get matrices (3 tables, one for each IFN type)
heatmap_matrix(group_names = group_names, table_name = table_name)

a <- grep(names(intersec_tables), pattern = paste0(group_names, collapse = "|"))
intersec_tables2 <- intersec_tables[a]

intersec_tables3 <- lapply(intersec_tables2, function(i){
  colnames(i) <- gsub(colnames(i), pattern = "_.*", replacement = "")
  i
})

merged_table <- reduce(.x = intersec_tables3, .f = full_join,
                         by = "Gene.symbol")
intersec_genes <- merged_table[1]

for (i in 1:length(intersec_tables3)) {
  nome <- gsub(names(intersec_tables3[i]), pattern = "early_vs_conv_tipo", replacement = "E")
  nome <- gsub(nome, pattern = "late_vs_conv_tipo", replacement = "L")
  intersec_tables3[[i]] <- mutate(intersec_tables3[[i]], type = nome)
}

merged_table_e <- reduce(.x = intersec_tables3[1:3], .f = full_join,
                         by = c("Gene.symbol", "GSE25001", "GSE28405", "GSE28988", "GSE28991", "GSE43777"))
colnames(merged_table_e) <- c("Gene.symbol", "GSE25001_E", "GSE28405_E", 
                              "GSE28988_E", "GSE28991_E", "GSE43777_E", "type_E_I", "type_E_II", "type_E_III")
merged_table_e <- merged_table_e %>% dplyr::relocate(type_E_I, type_E_II, type_E_III, .before = Gene.symbol)

merged_table_l <- reduce(.x = intersec_tables3[4:6], .f = full_join,
                         by = c("Gene.symbol", "GSE25001", "GSE28405", "GSE28988", "GSE28991", "GSE43777"))
colnames(merged_table_l) <- c("Gene.symbol", "GSE25001_L", "GSE28405_L", 
                              "GSE28988_L", "GSE28991_L", "GSE43777_L", "type_L_I", "type_L_II", "type_L_III")
merged_table_l <- merged_table_l %>% dplyr::relocate(type_L_I, type_L_II, type_L_III, .before = Gene.symbol)

merged_table_2 <- full_join(x = merged_table_e, y =  merged_table_l, by = "Gene.symbol")
merged_table_2 <- merged_table_2 %>% dplyr::relocate(type_L_I, type_L_II, type_L_III, .before = Gene.symbol)
intersec_labels <- merged_table_2[,1:7]

#(run heatmap matrix function)
tabelas <- intersec_tables2
folders <- unique(folders)

merged_table_3 <- full_join(x = intersec_labels, y =  FC_table, by = "Gene.symbol")

write.table(merged_table_3, file = paste0("./output/Intersection_Tables/Intersection_phases_morpheus.tsv"), 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

write.table(FC_table, file = paste0("./output/Intersection_Tables/Intersection_sev_morpheus.tsv"), 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

### Alluvial plot --------------------------------------------------------------
# Alluvial plot of GO BP for heatmap 
sankey_table <- read.table(file = "clipboard", sep = "\t")

sankey_list <- list()
sankey_dfs <- data.frame()
for (i in 1:(nrow(sankey_table))) {
  GO_id <- sankey_table[i,2]
  sankey_list[i] <- strsplit(sankey_table[i,8], split = ";")
  names(sankey_list)[[i]] <- GO_id
  df <- as.data.frame(sankey_list[i])
  df <- df %>% dplyr::mutate(go_id = GO_id)
  colnames(df)[1] <- "gene"
  sankey_dfs <- rbind(sankey_dfs, df)
} 

sankey_dfs <- dplyr::mutate(sankey_dfs, connection = 1)

intersec_genes <- read.table(file = "clipboard", sep = "\t")
colnames(intersec_genes)[1] <- "gene"

#merged <- merge.data.frame(intersec_genes, sankey_dfs, by.x = "V1", by.y = "gene", all = T, sort = F)
#merged2 <- full_join(intersec_genes, sankey_dfs, by = "gene")
merged3 <- right_join(intersec_genes, sankey_dfs, by = "gene")
merged3 <- merged3 %>% dplyr::mutate(gene = factor(gene, levels=unique(merged3$gene)))

sankey_table2 <- sankey_table %>% dplyr::mutate(bp_name = gsub(V1, pattern = "\\(.*\\)", replacement = ""))
colnames(sankey_table2)[2] <- "go_id"

merged4 <- full_join(merged3, sankey_table2[,c(2,9,1)], by = "go_id")

#go_id_levels <- read.table(file = "clipboard", sep = "\t", header = F)
#go_id_levels <- as.character(go_id_levels[,1])

#merged4 <- merged4 %>% dplyr::mutate(go_id = factor(go_id, levels=go_id_levels))

ggplot(merged3,
       aes(y = connection, axis1 = gene, axis2 = go_id)) +
  geom_alluvium(aes(fill = go_id), curve_type = "cubic") +
  geom_stratum(aes(fill = go_id)) +
  #geom_flow(aes(fill = go_id)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  #scale_fill_viridis_d() +
  theme_void()
#600 x 1500
  
  # matrix degs x BPs
  
  sankey_matrix <- ldply(as.list(intersec_genes[,1]), function(i){
    teste <- lapply(sankey_list, function(j){
      if (i %in% j) {1} else {0} 
    })
    teste <- as.numeric(teste)
  })
  
  colnames(sankey_matrix) <- names(sankey_list)
  row.names(sankey_matrix) <- intersec_genes[,1]
  
  write.table(sankey_matrix, file = "clipboard", sep = "\t", col.names = T, row.names = T, 
              quote = F)

# ENRICHMENT ANALYSIS ----------------------------------------------------------

# Clusters were selected visually on the heatmaps
# The genes were submitted to Enrichr for the enrichment analysis
# The GO Biological Process tables were downloaded
# To select the top 5 GO BP results:

Enrichrtxt <- list.files(path = "./data/Enrichr_GO_BP", pattern = "\\.txt", 
                            full.names = TRUE)

# export tables in the list as sheets in a xlsx
lapply(Enrichrtxt[-3], function(i){
  tabela <- read.table(file = i, header = TRUE, sep="\t", dec=".", quote = "")
  exclude <- grep(x = colnames(tabela), pattern = "Old")
  tabela2 <- tabela[-(exclude)]
  tabela2 <- tabela2 %>% dplyr::mutate(TermID = (gsub(x = Term, pattern = ".*\\(|\\)", "")))
  tabela2 <- tabela2 %>% dplyr::relocate(TermID, .after = Term)
  tabela3 <- subset(tabela2, tabela2$Adjusted.P.value < 0.05)
  xlname <- gsub(x = i, pattern = ".*/|\\.txt", replacement = "")
  xlname <- gsub(x = xlname, pattern = "Phases", replacement = "P")
  xlname <- gsub(x = xlname, pattern = "Sev", replacement = "S")
  write.xlsx(x = tabela3, file = "./output/Enrichr_filtered2.xlsx", 
             sheetName = xlname, row.names = F, append = TRUE)
})

#after submitting to revigo, results (revigo.csv) were downloaded
Revigo <- list.files(path = "./data/Enrichr_GO_BP", pattern = "revigo\\.csv", full.names = TRUE)
  tabela <- read.csv(file = i, header = TRUE, sep=",", dec=".")

# dotplots ---------------------------------------------------------------------
Enrichrtxt <- list.files(path = "./data/Enrichr_GO_BP", pattern = "\\.txt", 
                         full.names = TRUE, recursive = TRUE)
lapply(Enrichrtxt, function(i){
  nome <- gsub(x = i, pattern = ".*/|\\.txt", replacement = "")
  
  df <- read.table(file = i, header = TRUE, sep="\t", dec=".", quote = "")
  df <- df %>% filter(Adjusted.P.value < 0.05)
  df <- arrange(df, desc(Combined.Score))
  df$Term <- factor(df$Term, levels = df$Term[order(df$Combined.Score)])
  df <- df[1:15,]
  
  # create the plot
  ggplot(df, aes(x = Combined.Score, y = Term, color = Adjusted.P.value)) +
    geom_point(mapping = aes(size = Combined.Score)) +
    # scale_x_continuous(limits = c(0, 70)) +
    theme_dotplot +
    xlab ("Combined.Score") +
    ylab("GO category") +
    ggtitle(paste("DotPlot by GO -", nome))+
    scale_color_gradient(low="green", high="royalblue")
  
  ggsave(filename = paste0(nome, "_dotplot.png"), path = "./output/figures/dotplots", 
         #         width = 3700, height = 2100, 
         units = "px", dpi = 300, device = "png", scale = 1)
})


df$Term <- factor(df$Term, levels = df$Term[order(df$Adjusted.P.value)])

ggplot(df, aes(x = Adjusted.P.value, y = Term, color = Combined.Score)) +
  geom_point(mapping = aes(size = Adjusted.P.value)) +
  # scale_x_continuous(limits = c(0, 70)) +
  theme_dotplot +
  xlab ("Adjusted.P.value") +
  ylab("GO category") +
  ggtitle("DotPlot by GO")+
  scale_color_gradient(low="green", high="royalblue")

# cluster profiler parameters --------------------------------------------------
  
intersec_tables_list <- list.files(path = "./output/Intersection_Tables/", pattern = "\\.txt", 
                                     full.names = TRUE, recursive = F)
intersec_tables_list <- intersec_tables_list[-4]

intersec_tables <- lapply(intersec_tables_list, function(i){
  tabela <- read.table(file = i, header = TRUE, sep="\t", dec=".")
  geneList <- tabela[,2:3]
  colnames(geneList) <- c("ID", "LogFC")
  geneList
  })
names(intersec_tables) <- gsub(x = intersec_tables_list, pattern = ".+/|\\.txt", replacement = "") 

i <- 3
geneList <- intersec_tables[[i]]
nome <- names(intersec_tables)[[i]]

dotplot(egoBP, showCategory=20) + ggtitle("dotplot for BP")+ scale_color_viridis()
dotplot(egoBP, showCategory=15) + ggtitle("dotplot for BP")+ scale_color_viridis()

#svg: 550 x 500
ggsave(filename = "L_II_4.png", path = "./output/figures/dotplots", width = 1800, height = 2000, 
       units = "px", dpi = 300, device = "png", scale = 1)

egoxBP <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')

#early
goidBP <- c(  "GO:0009615", "GO:0060337", 
              "GO:0050792", #"GO:0048525", "GO:0019058","GO:1903900",
              "GO:0060333", "GO:0002697")
#late
goidBP <- c("GO:0006260", "GO:0007059", "GO:0000280", "GO:0000075", "GO:0044839"#, "GO:0071897"
            )

egoxBP@result <- egoxBP@result[goidBP, ]

ggsave(filename = "E_III_network1.png", path = "./output/figures/networks", width = 700, height = 600, 
       units = "px", dpi = 300, device = "png", scale = 1)

#BPs table
BPs <- as.data.frame(egoxBP@result)
#filtered <- dplyr::filter(BPs, p.adjust < 0.05) #p.adj or q val?
nome2 <- gsub(nome, pattern = "_vs_conv_tipo|_entrez", replacement = "")
write.xlsx(BPs, file = "./output/clusterprofiler_BPs.xlsx", sheetName = nome2, 
           col.names = T, row.names = F, append = T)
#write.table(BPs, file = "./output/clusterprofiler_BP_E_III.tsv", sep = "\t", 
#            col.names = T, row.names = F, quote = T)
  
# PCA --------------------------------------------------------------------------  

## GSE25001 expression matrix
  nome <- "./data/GSE25001_Whole/GSE25001_Whole_NetAn_SM_Entrez_Normalized_Filtered.txt" 
  #tabela1 <- read.table(file = nome, header = F , 
  #                      sep = "\t", dec = ".", row.names = 1)
  tabela2 <- read.table(file = nome, header = F, 
                        sep = "\t", dec = ".", row.names = 1, skip = 3)
  header <- read.table(file = nome, header = T, 
                        sep = "\t", dec = ".", row.names = 1,  nrows = 2)
  #GSE25001 entrez degs early vs conv
  #degs_table <- read_excel(path = "./output/GSE25001_Whole/GSE25001.53_72h_vs_34_convalescent/DEG_GSE25001.53_72h_vs_34_convalescent_SynGO.xlsx")
  degs_table <- read_excel(path = "./data/GSE25001_Whole/SynGO_25001_early_intersec.xlsx")
  #GSE25001 entrez degs late vs conv
  #degs_table <- read_excel(path = "./output/GSE25001_Whole/GSE25001.100_4to8d_vs_34_convalescent/DEG_GSE25001.100_4to8d_vs_34_convalescent_SynGO.xlsx")
  degs_table <- read_excel(path = "./data/GSE25001_Whole/SynGO_25001_late_intersec.xlsx")

  tabela_DEGs <- merge.data.frame(y= degs_table[,c(1,2)], x = tabela2, by.y = "entrezgene", by.x = "row.names")
  tabela_DEGs <- dplyr::relocate(tabela_DEGs, query, Row.names, .before = 1)
  tabela_DEGs <- tabela_DEGs[,-2]
  tabela_DEGs <- tabela_DEGs %>% remove_rownames %>% column_to_rownames(var = "query")
  
  # select samples, change to remove early or late samples accordingly
  nas <- which(is.na(header[2,]))
  header2 <- header[,-nas]
  excl <- which(header2[1,] == "Early")
  header3 <- header2[,-excl]
  
  tabela_DEGs2 <- tabela_DEGs[,-nas]
  tabela_DEGs2 <- tabela_DEGs2[,-excl]
  # skip to 'common'
  
## GSE43777 expression matrix
  nome2 <- "./data/GSE43777_Sequential/GSE43777_filtered/data_normalized.csv"
  #GSE43777 entrez degs early vs conv
  degs_table <- read_excel(path = "./data/GSE43777_Sequential/SynGO_43777_early_intersec.xlsx")
  #GSE43777 entrez degs late vs conv
  degs_table <- read_excel(path = "./data/GSE43777_Sequential/SynGO_43777_late_intersec.xlsx")
  
  tabela2 <- read.table(file = nome2, header = T , 
                         sep = ",", dec = ".", row.names = 1)
  header2 <- read.table(file = "./data/GSE43777_Sequential/GSE43777_filtered/GSE43777_NetAn_SM.txt", header = T, 
                         sep = "\t", dec = ".", row.names = 1,  nrows = 2)
  #degs_table <- degs_table[-(which(duplicated(degs_table[,2]))),]

  tabela_DEGs <- merge.data.frame(y= degs_table[,c(1,2)], x = tabela2, by.y = "entrezgene", by.x = "row.names")
  tabela_DEGs <- dplyr::relocate(tabela_DEGs, query, Row.names, .before = 1)
  tabela_DEGs <- tabela_DEGs[,-2]
  tabela_DEGs <- tabela_DEGs %>% remove_rownames %>% column_to_rownames(var = "query")

  # select samples, change to remove early or late samples accordingly
  excl <- which(header2[1,] == "Early")
  header3 <- header2[,-excl]
  tabela_DEGs2 <- tabela_DEGs
  tabela_DEGs2 <- tabela_DEGs2[,-excl]

  #remove outliers from GSE43777
  excl2 <- grep(colnames(header3), pattern = "GSM1071100|GSM1071104|GSM1071105|GSM1071108")
  header3 <- header3[,-excl2]
  tabela_DEGs2 <- tabela_DEGs2[,-excl2]
  
  
# common
header3[3,] <- paste(header3[1,],header3[2,], sep = "_")
groups <- as.factor(header3[3,])
levels(groups)
tabela4 <- t(tabela_DEGs2)
tabela4 <- as.data.frame(tabela4)
tabela5 <- cbind(as.data.frame(groups), tabela4)
#write.csv(tabela5, file = "./output/43777_Late_PCA.csv",  row.names = F)

PCAf <- tabela5
#PCAf <- filter(PCAf, groups > 2 )
groups <- as.factor(PCAf$groups)
PCAf$groups <- NULL  

res.pca <- prcomp(PCAf, scale = TRUE)

res.pca4 <- res.pca
groups4 <- groups
# run PCA script

# RF matrix only: 
excl2 <- which(header3[1,] == "Convalescent")
header3 <- header3[,-excl2]
tabela_DEGs2 <- tabela_DEGs2[,-excl2]

tabela5 <- cbind(tabela4, as.data.frame(as.numeric(groups)))
colnames(tabela5)[ncol(tabela5)] <- "Group"
write.csv(tabela5, file = "./25001_Early_RF.csv",  row.names = F)

# reorder groups
groups11 <- factor(groups1, levels = c("Convalescent_nonDSS", "Convalescent_DSS", "Early_DSS", "Early_nonDSS"))
groups12 <- factor(groups2, levels = c("Convalescent_nonDSS", "Convalescent_DSS", "Late_DSS", "Late_nonDSS"))
groups13 <- factor(groups3, levels = c("Convalescent_DF", "Convalescent_DHF",  "Early_DHF", "Early_DF"))
groups14 <- factor(groups4, levels = c("Convalescent_DF", "Convalescent_DHF", "Late_DHF", "Late_DF"))


# SEVERITY METANALYSIS --------------------------------------------------------- 

# expression matrix ------------------------------------------------------------
nome <- "./data/Severity_metanalysis/meta_late_0/ExpressAnalyst_merged_data.txt" 
#tabela1 <- read.table(file = nome, header = T , 
#                      sep = "\t", dec = ".", row.names = 1)
tabela2 <- read.table(file = nome, header = F, 
                      sep = "\t", dec = ".", row.names = 1, skip = 3)
header <- read.table(file = nome, header = T, 
                     sep = "\t", dec = ".", row.names = 1,  nrows = 2, comment.char = "")

# entrez and gene name table
entrez_table <- read_excel(path = "./data/Severity_metanalysis/meta_late_0/meta_synGO.xlsx")

# intersection IFNome genes 
# (metassignificant and genes common to at least 3 datasets in late vs conv)
intersec <- read.table(file = "clipboard", header = F) 
colnames(intersec) <- "symbol"
tabela_DEGs <- merge.data.frame(x = intersec, y = entrez_table, by = "symbol")
tabela_DEGs <- drop_na(data = tabela_DEGs, "name")

# merge with expression matrix
tabela_DEGs2 <- merge.data.frame(y= tabela_DEGs[,c(1,2)], x = tabela2, by.y = "query", by.x = "row.names")

tabela_DEGs2 <- dplyr::relocate(tabela_DEGs2, symbol, Row.names, .before = 1)
tabela_DEGs2 <- tabela_DEGs2[,-2]
tabela_DEGs2 <- tabela_DEGs2 %>% remove_rownames %>% column_to_rownames(var = "symbol")

replace <- which(header[2,] == "GSE25001_non-normalized_NetAn_late_log2_0.txt")
header[2,replace] <- "GSE25001"
replace <- which(header[2,] == "GSE43777_GPL570_NetAn_SM_late.txt")
header[2,replace] <- "GSE43777"

header[3,] <- paste(header[1,],header[2,], sep = "_")
groups <- as.factor(header[3,])
levels(groups)
tabela4 <- t(tabela_DEGs2)
tabela4 <- as.data.frame(tabela4)
tabela5 <- cbind(as.data.frame(groups), tabela4)
write.csv(tabela5, file = "./data/Severity_metanalysis/meta_PCA.csv",  row.names = F)

PCAf <- tabela5
#PCAf <- filter(PCAf, groups > 2 )
groups <- as.factor(PCAf$groups)
PCAf$groups <- NULL  

res.pca <- prcomp(PCAf, scale = TRUE)

# combat 
BiocManager::install("sva")
library("sva")

tabela_DEGs3 <- tabela_DEGs2
colnames(tabela_DEGs3) <- colnames(header)
matriz <- as.matrix(tabela_DEGs3)
lote <- header[2,]
lote <- as.factor(lote)
lote <- as.numeric(lote)

# parametric adjustment
combat1 = ComBat(dat=matriz, batch=lote, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

tabela6 <- t(combat1)
tabela6 <- as.data.frame(tabela6)
tabela7 <- cbind(as.data.frame(groups), tabela6)
write.csv(tabela7, file = "./data/Severity_metanalysis/meta_combat_PCA.csv",  row.names = F)
PCAf <- tabela7


# RF matrix 
groups <- as.factor(header[1,])
levels(groups)

tabela8 <- cbind(tabela6, as.data.frame(as.numeric(groups)))
colnames(tabela8)[ncol(tabela8)] <- "Group"
write.csv(tabela8, file = "./meta_Late_RF.csv",  row.names = F)


# combat table for Express Analyst
tabela9 <- as.data.frame(combat1)
header2 <- read.table(file = nome, header = F, 
                      sep = "\t", dec = ".", row.names = 1,  nrows = 3, comment.char = "")
header2[3,] <- header[2,]
colnames(header2) <- colnames(header)
tabela10 <- rbind.data.frame(header2, tabela9)
write.table(tabela10, file = "./data/Severity_metanalysis/meta_combat_NetAn.csv",
            row.names = T, col.names = F, sep = ",")


# after Differential Expression Analysis in ExpressAnalyst
# get the DEA results for the 10 genes selected in the random forest
tabela11 <- read.table(file = "./data/Severity_metanalysis/SigGene_DSS_vs_nonDSS.csv", header = T, 
                     sep = ",", dec = ".", quote = "\"\"\"")

# the 10 genes selected in the random forest
RF_res <- read.table(file = "clipboard", header = F) 
colnames(RF_res) <- "Symbols"
tabela12 <- inner_join(x = RF_res, y = tabela11, by = "Symbols")
write.table(tabela12, file = "clipboard",
            row.names = F, col.names = T, sep = ",")


## heatmap - expression of 10 RF genes, expression after combat

headerA <- read.table(file = "./data/GSE25001_Whole/GSE25001_NetAn_SM.txt", header = T, 
                       sep = "\t", dec = ".", row.names = 1, nrows = 2)
headerE <- read.table(file = "./data/GSE43777_Sequential/GSE43777_GPL570_NetAn_SM.txt", header = T, 
                      sep = "\t", dec = ".", row.names = 1, nrows = 2)
header_joined <- bind_cols(headerA, headerE)

RF_genes <- c('CA1', 'DYRK2', 'SLPI', 'NFE2L3', 'SELENBP1', 'ARL4A', 'MS4A4A',
              'LAG3', 'RUNDC3A', 'SESN3')
RF_genes <- as.data.frame(RF_genes)

tabela_combat <- read.table(file = "./data/Severity_metanalysis/meta_combat_NetAn.csv",
                            header = T,sep = ",", dec = ".",  comment.char = "")
#tabela_combat2 <- filter(tabela_combat, X.NAME %in% RF_genes)
tabela_combat3 <- merge.data.frame(x = RF_genes, y = tabela_combat, 
                                   by.x = 'RF_genes', by.y = "row.names", sort = F)
tabela_combat3 <- column_to_rownames(tabela_combat3, var = 'RF_genes')
tabela_combat4 <- bind_rows(tabela_combat3, header_joined)
write.table(tabela_combat4, file = "clipboard",
            row.names = F, col.names = T, sep = ",")
