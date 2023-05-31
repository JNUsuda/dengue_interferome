# DEG FILTER-EXPORT FUNCTION ---------------------------------------------------

DEG_filter <- function(original_name) {
  # file and path names 
  out_path <- gsub(x = original_name, pattern =".tsv|.txt", replacement = "")
  out_path <- gsub(x = out_path, pattern ="/data/", replacement = "/output/")
  new_folder <- gsub(x = out_path, pattern =".*/", replacement = "")
  
  # create directory for GSE, if it doesn't already exists
  GSE_path <- file.path(dirname(out_path))
  if (!(GSE_path %in% list.dirs(path = "./output"))){
    dir.create(GSE_path)} 
  # create directory for the comparison, inside the GSE folder
  if (!(out_path %in% list.dirs(path = "./output"))){
    dir.create(out_path)} 
  
  # table import
  table <- read.table(file = original_name, header = TRUE, sep = "\t", dec = ".", quote = "" )
  
  # if Dengue dataset, get gene symbols:
  if(basename(GSE_path) == "GSE40628_Dengue"){
    adf_table <- read.table(file = "./data/GSE40628_Dengue/GPL16021_Print_980_edited.txt",
                            header = TRUE, sep = "\t", quote = "", fill = T)
    adf_table2 <- dplyr::mutate(adf_table, printID = (gsub(adf_table$Feature.Identifier, 
                  pattern = "URN:LSID:smd.stanford.edu:feature.feature\\|feature.printid:", replacement = "")))
    adf_table2 <- adf_table2 %>% dplyr::relocate(printID, .before = adf_table2$Block.Column)
    adf_table3 <- adf_table2[,c("printID","Reporter.Database.Entry.image.", 
                                "Composite.Element.Database.Entry.Gene.Symbol.", 
                                "Composite.Element.Database.Entry.Geneid..Locusid..")]
    colnames(adf_table3) <- c("ID", "IMAGE", "Gene.symbol", "Entrez_ID")
    
    GEO2R_table <- dplyr::mutate(table, image = (gsub(table$CLONE_ID, pattern = "IMAGE:", replacement = "")))
    GEO2R_table <- GEO2R_table %>% dplyr::relocate(image, .after = CLONE_ID)
    
    table <- merge(x=GEO2R_table, y = adf_table3, by = "ID",  all.x = T) 
    write.table(table, file = (paste0(out_path, "/GeneSymbol_", new_folder, ".tsv")), 
                sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  
  #if RNA-seq dataset, rename columns
  if(basename(GSE_path) == "GSE94892_Transcriptome_RNAseq"){
    colnames(table)[colnames(table)%in% c("gene", "log2.fold_change.", 
    "p_value", "q_value")] <- c("Gene.symbol", "logFC",  "P.Value", "adj.P.Val")}
  
  # subset by adjusted p-value < 0.05
  p_filtered <- subset(x = table, table$adj.P.Val < 0.05)
  # subset by log FC < -1 OR > 1
  p_FC_filtered <- subset(x = p_filtered, subset = (p_filtered$logFC < -1 | p_filtered$logFC > 1 ) )
  # remove duplicates by gene symbol, after filtering
  DEGs <- p_FC_filtered[!duplicated(p_FC_filtered$Gene.symbol), ]

  # remove empty values
  if (any(nchar(DEGs$Gene.symbol) == 0) == TRUE){
    to_remove <- which(nchar(DEGs$Gene.symbol) == 0)
    DEGs <- DEGs[-to_remove,]}
  if (any(DEGs$Gene.symbol == " ") == TRUE){
    to_remove <- which(DEGs$Gene.symbol == " ")
    DEGs <- DEGs[-to_remove,]}
  
  # subset upregulated and downregulated genes
  upreg <- subset(x = DEGs, subset = (DEGs$logFC > 1))
  downreg <- subset(x = DEGs, subset = (DEGs$logFC < -1))
  
  # Export tables as .tsv:
  # all filtered DEGs
  write.table(DEGs, file = (paste0(out_path, "/DEG_", new_folder, ".tsv")), 
              sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)
  # upregulated DEGs
  write.table(upreg, file = (paste0(out_path, "/UP_", new_folder, ".tsv")), 
              sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)
  # downregulated DEGs
  write.table(downreg, file = (paste0(out_path, "/DOWN_", new_folder, ".tsv")), 
              sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # export excel spreadsheet with the all, upregulated and downregulated DEGs as 'tabs'
  if (!(is.null(DEGs) | nrow(DEGs) == 0)) { 
    write.xlsx(x = DEGs, file = (paste0(out_path, "/excel_", new_folder, ".xlsx")), 
               sheetName = "DEG", row.names = FALSE)}
  if (!(is.null(upreg)| nrow(upreg) == 0)) {
    write.xlsx(x = upreg, file = (paste0(out_path, "/excel_", new_folder, ".xlsx")), 
               sheetName = "UP", row.names = FALSE, append = TRUE)}
  if (!(is.null(downreg)| nrow(downreg) == 0)) {
    write.xlsx(x = downreg, file = (paste0(out_path, "/excel_", new_folder, ".xlsx")), 
               sheetName = "DOWN", row.names = FALSE, append = TRUE)}
  
}

# GENE LIST BY IFN TYPE FUNCTION -----------------------------------------------

genes_by_type <- function(IFN_type) {
  # get gene lists
  TypeGenes <- lapply(DataSearchRes, function(i){
    table <- read.table(file = i, header = TRUE, sep = "\t", dec = ".", quote = "", skip = 19)
    table1 <- subset(x = table, table$Inteferome.Type == IFN_type)
    table2 <- table1[!duplicated(table1$Gene.Name), ]
    table2$Gene.Name})
  # turn into df
  names(TypeGenes) <- gsub(x = DataSearchRes, pattern = "/[[:alnum:]]+_[[:alnum:]]+\\.txt|\\./output/", replacement = "") 
  TypeGenesdf <- ldply(TypeGenes ,rbind)
  # get comparisons names, add col with IFN type, arrange cols
  TypeGenesdf <- dplyr::mutate(TypeGenesdf, 
                  folder1 = (gsub(x = TypeGenesdf$.id, pattern = "/.*", "")))
  TypeGenesdf <- TypeGenesdf %>% dplyr::relocate(folder1, .after = .id)
  TypeGenesdf <- dplyr::mutate(TypeGenesdf, 
                  temp1 = (gsub(x = TypeGenesdf$.id, pattern = "^[[:alnum:]_]+/", "")))
  TypeGenesdf <- TypeGenesdf %>% dplyr::relocate(temp1, .after = folder1)
  TypeGenesdf <- dplyr::mutate(TypeGenesdf, 
                  folder2 = (gsub(x = TypeGenesdf$temp1, pattern = "/.*", "")))
  TypeGenesdf <- TypeGenesdf %>% dplyr::relocate(folder2, .after = temp1)
  TypeGenesdf <- dplyr::mutate(TypeGenesdf, 
                  UPDOWN = (gsub(x = TypeGenesdf$.id, pattern = ".*_", "")))
  TypeGenesdf <- TypeGenesdf %>% dplyr::relocate(UPDOWN, .after = folder2)
  TypeGenesdf$IFN_TYPE <- IFN_type
  TypeGenesdf <- TypeGenesdf %>% dplyr::relocate(IFN_TYPE, .after = folder2)
  TypeGenesdf
}

# GENE INTERSECTION FUNCTION ---------------------------------------------------

genes_intersection <- function(compar_names, group_name) {
  for (IFN_type in c("I", "II", "III")) {
    # for each comparison:
    # load the interferome UP DATA table, filter by IFN type, remove gene symbol duplicates
    # load the interferome DOWN DATA table, filter by IFN type, remove gene symbol duplicates
    # join gene symbols from UP and DOWN, save result, add names
    genelists <- lapply(compar_names, function(i){
      my_fileUP <- list.files(path = paste0("./output/",i ,"/Interferome_hs_UP"), pattern = ".+_DataSearchResults\\.txt", 
                              full.names = TRUE, recursive = F)
      tableUP <- read.table(file = my_fileUP, header = TRUE, sep = "\t", dec = ".", quote = "", skip = 19)
      tableUP <- subset(x = tableUP, tableUP$Inteferome.Type == IFN_type)
      tableUP <- tableUP[!duplicated(tableUP$Gene.Name), ]
      my_fileDOWN <- list.files(path = paste0("./output/",i ,"/Interferome_hs_DOWN"), pattern = ".+_DataSearchResults\\.txt", 
                                full.names = TRUE, recursive = F)
      tableDOWN <- read.table(file = my_fileDOWN, header = TRUE, sep = "\t", dec = ".", quote = "", skip = 19)
      tableDOWN <- subset(x = tableDOWN, tableDOWN$Inteferome.Type == IFN_type)
      tableDOWN <- tableDOWN[!duplicated(tableDOWN$Gene.Name), ]
      genes <- c(tableUP$Gene.Name, tableDOWN$Gene.Name)
      genes})
    # add the folders names to the gene lists
    names(genelists) <- compar_names
    
    # get intersection of the comparisons
    # create df from intersect genes as column, rename column
    interseccao <- Reduce(intersect, genelists)
    interseccao <- as.data.frame(interseccao)
    colnames(interseccao) <- "Gene.symbol"
    
    # make df with intersection genes and logFC:
    # for each comparison, create new col with logFC, based on gene symbol
    FC_list <- lapply(compar_names, function(i){
      my_file <- list.files(path = paste0("./output/",i ), pattern = "DEG.+\\.tsv", 
                            full.names = TRUE, recursive = F)
      my_table <- read.table(file = my_file, header = TRUE, sep = "\t", dec = ".", quote = "")
      colunas <- match(c("logFC","Gene.symbol"), colnames(my_table))
      intersec_table <- merge.data.frame(x = interseccao, y = my_table[, colunas], by = "Gene.symbol" , all.x = TRUE)
      colnames(intersec_table)[colnames(intersec_table) == "logFC"] <- i
      intersec_table})
    # join all in one table
    FC_table <- Reduce(merge.data.frame, FC_list)
    
    #write table
    FC_table_name <- paste0("./output/Intersection_Tables/", group_name, "_tipo_", IFN_type, ".tsv")
    write.table(FC_table, file = FC_table_name, sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }}

# HEATMAP MATRIX FUNCTION ------------------------------------------------------

heatmap_matrix <- function(group_names, table_name) {
  for (IFN_type in c("I", "II", "III")) {
    
    # get intersection tables
    names1 <- paste0(group_names,"_tipo_", IFN_type, "$", collapse = "|")
    pos <- grep(names(intersec_tables), pattern = names1)
    tabelas <- intersec_tables[pos]
    
    # get names of included comparisons
    compari_names <- lapply(tabelas, function(i){
      colnames(i[-1])})
    folders <- unlist(compari_names, use.names = FALSE)
    folders <- sub(x = folders, pattern = "\\.", replacement = "/")
    
    # merge tables by gene symbol, keep all genes from all
    # actually get only the merged gene.symbols column 
    if (any(lapply(tabelas, nrow) == 0)) {
      ignore <- which(lapply(tabelas, nrow) == 0)
      merged_table <- reduce(.x = tabelas[-ignore], .f = full_join, by = "Gene.symbol")
    } else {
      merged_table <- reduce(.x = tabelas, .f = full_join, by = "Gene.symbol")}
    intersec_genes <- merged_table[1]
    
    # make df with intersection genes and logFC:
    # for each comparison, 
    # get original data table
    FC_list <- lapply(folders, function(i){
      if(grepl(pattern = "GSE94892_Transcriptome_RNAseq", x = i, fixed = T)){
        my_table <- read.table(file = paste0("./data/", i, ".txt" ), 
                               header = TRUE, sep = "\t", dec = ".", quote = "")
        colnames(my_table)[colnames(my_table)%in% c("gene", "log2.fold_change.",
        "p_value", "q_value")] <- c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")
      }
      if(grepl(pattern = "GSE40628_Dengue", x = i, fixed = T)){
        my_table <- read.table(file = paste0("./output/",i, "/GeneSymbol_", basename(i), ".tsv" ),
                               header = TRUE, sep = "\t", dec = ".", quote = "")
      } else {
        my_table <- read.table(file = paste0("./data/",i, ".tsv" ), 
                               header = TRUE, sep = "\t", dec = ".", quote = "")
      }
    # filter by p-value and remove duplicates
      my_table1 <- subset(x = my_table, my_table$adj.P.Val < 0.05)
      my_table2 <- my_table1[!duplicated(my_table1$Gene.symbol), ]
    # create new col with logFC, based on gene symbol
      colunas <- match(c("logFC","Gene.symbol"), colnames(my_table2))
      intersec_table <- merge.data.frame(x = intersec_genes, 
                        y = my_table2[, colunas], by = "Gene.symbol", all.x = TRUE)
      colnames(intersec_table)[colnames(intersec_table) == "logFC"] <- gsub(i, pattern = ".+/", replacement = "")
      intersec_table
    })
    # join all in one table
    FC_table <- reduce(.x = FC_list, .f = full_join, by = "Gene.symbol")

    # replace NA's with 0's, write table
    FC_table[is.na(FC_table)] <- 0
    write.table(FC_table, file = paste0("./output/Intersection_Tables/", table_name, "_tipo_", IFN_type, ".tsv"), 
                sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)
}}

#NETWORKANALYST FUNCTIONS ------------------------------------------------------

#Function to turn gsms and group from geo2R to group classification
NetAn_groups <- function(gsms, groups) {
  sml <- strsplit(gsms, split="")[[1]]
  if (any(sml == "X")) {
    Xs <- which(sml == "X")
    sml[Xs] <- (max(as.double(sml[-Xs]))+ 1)
    Xgroups <- make.names(c(groups,"X"))
    gs <- factor(sml)
    levels(gs) <- Xgroups
    grupos <- as.character(gs)
    grupos <- gsub(grupos, pattern = "X", replacement = "NA")
    grupos
  } else {
    gs <- factor(sml)
    levels(gs) <- groups
    grupos <- as.character(gs)
    grupos
  }
}

#Function to export the series matrix (SM) table with sample classification 
#on a format suited to ExpressAnalyst 
NetAn_SM <- function(gset, gsmsS, groupsS, gsmsP, groupsP, GPLannot) {
  df <- as.data.frame(ex)
  
  if (GPLannot == TRUE) {
    probes <- fData(gset)
    cols <- which(colnames(probes) %in% c("ID","Gene ID"))
    
    Genedf <- merge(x = probes[,cols], y = df, by ='row.names', all = T)
    Genedf <- Genedf[,-c(1,2)]
    colnames(Genedf)[1] <- "Gene.id"
    Genedf <- Genedf %>% filter(Genedf$Gene.id != "")
    Genedf <- Genedf %>% filter(!is.na(Genedf$Gene.id))

    if (is.null(gsmsS)== TRUE) {
      df2 <- add_row(.data = Genedf, .before = 1)
      df2[1,1] <- "#CLASS:Severity"
    } else {
      gruposS <- NetAn_groups(gsms = gsmsS, groups = groupsS)
      gruposS <- c("#CLASS:Severity", gruposS)
      df2 <- rbind(gruposS, Genedf)
    }

    if (is.null(gsmsP)== TRUE) {
      df3 <- add_row(.data = df2, .before = 1)
      df2[1,1] <- "#CLASS:Phases"
    } else {
      gruposP <- NetAn_groups(gsms = gsmsP, groups = groupsP)
      gruposP <- c("#CLASS:Phases", gruposP)
      df3 <- rbind (gruposP, df2)
    }
    colnames(df3)[1] <- "#NAME"
    
    folders <- list.dirs(path = "./data", recursive = F)
    name <- paste0(folders[grep(folders, pattern = gse)],"/", gse, "_NetAn_SM.txt")
    write.table(df3, file = name, quote = T, sep = "\t", row.names = F, dec = ".", col.names = T)
  }
  else {
    if (is.null(gsmsS)== TRUE) {
      df2 <- add_row(.data = df, .before = 1)
    } else {
      gruposS <- NetAn_groups(gsms = gsmsS, groups = groupsS)
      df2 <- rbind(gruposS, df)
    }
    rownames(df2)[1] <- "#CLASS:Severity"
    
    if (is.null(gsmsP)== TRUE) {
      df3 <- add_row(.data = df2, .before = 1)
    } else {
      gruposP <- NetAn_groups(gsms = gsmsP, groups = groupsP)
      df3 <- rbind (gruposP, df2)
    }
    rownames(df3)[1] <- "#CLASS:Phases"
    
    df4 <- tibble::rownames_to_column(df3, "row_names")
    colnames(df4)[1] <- "#NAME"
    
    folders <- list.dirs(path = "./data", recursive = F)
    name <- paste0(folders[grep(folders, pattern = gse)],"/", gse, "_NetAn_SM.txt")
    write.table(df4, file = name, quote = T, sep = "\t", row.names = F, dec = ".", col.names = T)
  }
}

