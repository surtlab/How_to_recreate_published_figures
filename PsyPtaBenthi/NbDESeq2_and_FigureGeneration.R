## RNA-seq of Nicotiana benthamiana responses to Pta and Psy
#Morgan Carter 2021-2022
#Anything within the prepped directory is an input file, including genomes, count data, and phenotype csv for DESeq2

# Packages and setup ------------------------------------------------------
install.packages('DESeq2')
install.packages('ballgown')
install.packages('dplyr')
install.packages('goseq')
install.packages('ggplot2')
install.packages('pheatmap')
install.packages('viridis')
install.packages('VennDiagram')
install.packages('patchwork')
install.packages('ggforce')

library(DESeq2)
library(ballgown)
library(dplyr)
library(goseq)
library(ggplot2)
library(pheatmap)
library(viridis)
library(VennDiagram)
library(patchwork)
library(ggforce)
setwd("C:/Users/sumwh/Box/BaltrusLab/Ptab_B728a_RNAseq3") #laptop
setwd("E:/Box/BaltrusLab/Ptab_B728a_RNAseq3") #Desktop


# DESeq2 Analysis ---------------------------------------------------------

#reading in gene count files and formatting gene names
counts_data_Nb <- read.csv("prepped/Nb101_gene_count_matrix.csv", header=T)
counts_data_Nb$gene_id <- lapply(counts_data_Nb$gene_id, function(x) gsub("STRG.[0-9]*\\|(NbD......$)","\\1.1", x)) #regex to change STRG to .1 transcripts
counts_data_Nb$gene_id <- lapply(counts_data_Nb$gene_id, function(x) gsub(" CDS","", x))#removing CDS
counts_data_Nb %>% #combining counts from duplicate gene rows in weird cases
        group_by(gene_id) %>%
        summarise_all(sum) %>%
        data.frame() -> counts_data_Nb

#formatting it for DESeq count table
row.names(counts_data_Nb) <- counts_data_Nb$gene_id
counts_data_Nb <-counts_data_Nb[,-1]
counts_data_Nb <-counts_data_Nb[-nrow(counts_data_Nb),] #removing final row with NAs from stringtie output

#reading in phenotype table and turning integers into factors
coldata_Nb <- read.csv("prepped/Nb_Pheno_combined.csv", row.names = "ids")

#checking that rownames in the phenotype and column names in gene count files match
all(rownames(coldata_Nb) %in% colnames(counts_data_Nb))
all(rownames(coldata_Nb) == colnames(counts_data_Nb))

#the main function as outlined by DESeq manual to make data set
dds_Nb_combo <- DESeqDataSetFromMatrix(countData = counts_data_Nb, 
                                 colData = coldata_Nb, 
                                 design = ~ group)

#running DEseq
dds.dataNb <- DESeq(dds_Nb_combo)

#creating the output table showing all genes and whether they are DEGs
contrast_comboB <- c("group", "5B", "0B")
contrast_comboT <- c("group", "5T", "0T")

#Test to change Cook's distance cut off
m=12 #number of samples
b=4 #number of conditions
cooks_cut=qf(.97, b, m - b) #determining the cooksCutoff to set if we want to use a 97% quantile instead of 99%
res_Nb_comboB<-results(dds.dataNb,
                       contrast=contrast_comboB,
                       independentFiltering = T,
                       alpha = 0.01,
                       cooksCutoff = cooks_cut)
summary(res_Nb_comboB)
#graphs to look into that
plot(mcols(dds.dataNb)$baseMean, mcols(dds.dataNb)$maxCooks) #plot to see Cooks distance for each gene plotted against base mean
stopifnot(all.equal(rownames(dds.dataNb), rownames(res_Nb_comboB)))
plot(res_Nb_comboB$log2FoldChange, mcols(dds.dataNb)$maxCooks) #plot to see Cooks distance for each gene plotted against log2FoldChange

#pulling results from the dds object based on parameters
res_Nb_comboB<-results(dds.dataNb,
                     contrast=contrast_comboB,
                     independentFiltering = T,
                     alpha = 0.01,
                     cooksCutoff = cooks_cut)
summary(res_Nb_comboB)
write.csv(as.data.frame(res_Nb_comboB),file="DEseq_NbComboB_out.csv")

res_Nb_comboT<-results(dds.dataNb,
                       contrast=contrast_comboT,
                       independentFiltering = T,
                       alpha = 0.01,
                       cooksCutoff = cooks_cut)
summary(res_Nb_comboT)
write.csv(as.data.frame(res_Nb_comboT),file="DEseq_NbComboT_out.csv")

# Creating plots ----------------------------------------------------------
#note some plots were further modified for aesthetics (text and colors) in Illustrator
#rlog transformation
rld_Nb_combo <- rlog(dds_Nb_combo, blind=F)

#PCA plots
plotPCA(rld_Nb_combo, intgroup=c("group"))

#Heatmap
degs_namesB <-
        rownames(subset(res_Nb_comboB,
                        padj < 0.01))
degs_namesT <-
        rownames(subset(res_Nb_comboT,
                        padj < 0.01))
degs_names <- union(degs_namesB, degs_namesT)

ann_colors = list(
        group = c("0B" ="#5BA8A0", "0T" = "#CBe54E", "5B" = "#3B5284", "5T"= "#94B447")
)

pheatmap(assay(rld_Nb_combo)[degs_names,],
         show_rownames = F,
         filename = "Plots/Nb_heatmap.pdf", 
         annotation_col = coldata_Nb,
         scale="row",
         color=magma(15),
         annotation_colors = ann_colors,
         show_colnames = F
         )

# Adding annotations from gff file to DESeq2 output ---------------
#converts gff file to useful info from last column for manipulating in R
gffdata_mRNA = gffRead('prepped/Nb101.gff')
gff_stuff <- c("Name", "gene", 'Codon_start', "product", "protein_id", "NCBI Feature Key", "NCBI Join Type", "ID") #will be column headers for all of attributes
GeneInfo <- as.data.frame(matrix(ncol=length(gff_stuff), nrow=0, dimnames=list(NULL,gff_stuff))) #making a new, empty data frame
gffdata_gene_mod <- gffdata_mRNA[9]
gffdata_gene_mod <- distinct(gffdata_mRNA[9]) #there are multiple lines for each gene for different segments of CDS, this removes duplicates
num <- dim(gffdata_gene_mod)
for (i in 1:num[1]) { #for every line in the gff file (this takes a while because so many genes in genome)
        GeneInfo[i,] <- NA
        teststr <- unlist(strsplit(as.character(gffdata_gene_mod[i,1]), ";", fixed = T))
        for (j in 1:length(teststr)) {
                teststr2 <-unlist(strsplit(teststr[j], "=", fixed = T))
                GeneInfo[i,teststr2[1]] <- teststr2[2]
        }
}
GeneInfo2 <- as.data.frame(matrix(ncol=2, nrow=0, dimnames=list(NULL, c("product", "protein_id")))) #reducing to just the protein info and ids from the product col
for (i in 1:num[1]) { 
        GeneInfo2[i,] <- NA
        teststr <- unlist(strsplit(as.character(GeneInfo$product[i]), "  (", fixed = T))
        GeneInfo2[i,1] <- teststr[1]
        GeneInfo2[i,2] <- gsub(")", "",teststr[2])
}
row.names(GeneInfo2) <- GeneInfo$protein_id
write.csv(GeneInfo2, file="Nb101_gff_info.csv")

#Adding GO terms from a separate csv file
Nb101_GO <- read.csv("Prepped/Nb101_GO.csv", header=T, fill=T)
Nb101_GO <- Nb101_GO[,1:12]
colnames(Nb101_GO)[1] <- "Name"
gene_info <- read.csv("Nb101_gff_info.csv", header=T)
DeSeq_Info_0 <- merge.data.frame(gene_info, Nb101_GO, by.x="X", by.y="Name", all=T)
write.csv(DeSeq_Info_0, file="Nb101_gff_info_GO.csv")

#combining with the DESeq2 output - Ptab and B728a
DE_data1 <- read.csv("DEseq_NbComboT_out.csv", header=T)
colnames(DE_data1) <- c("X", "baseMean.Ptab", "L2FC.Ptab", "lfcSE.Ptab", "stat.Ptab", "pvalue.Ptab", "padj.Ptab")
DE_data2 <- read.csv("DEseq_NbComboB_out.csv", header=T)
colnames(DE_data2) <- c("X", "baseMean.B728a", "L2FC.B728a", "lfcSE.B728a", "stat.B728a", "pvalue.B728a", "padj.B728a")
gene_info <- read.csv("Nb101_gff_info_GO.csv", header=T)
gene_info <- gene_info[,-1]
DeSeq_Info_T <- merge.data.frame(gene_info, DE_data1, by="X")
DeSeq_Info_TB <- merge.data.frame(DeSeq_Info_T, DE_data2, by="X")

write.csv(DeSeq_Info_TB, file="Annotated_Out/Nb101_DESeq_All_Annotated_GO.csv")


# Subsetting DE Genes -----------------------------------------------------
#Taking subsets for different groups of DE genes between conditions for Ptab/B728a
DESeq_GO<- read.csv("Annotated_Out/Nb101_DESeq_All_Annotated_GO.csv", header=T)

#one comparison only, exclude shared with other comparisons (for two I had to include is.na for padj because the numbers werent lining up with the venn diagram)
DESeq_GO.B728a <- filter(DESeq_GO, padj.B728a<0.01 & abs(L2FC.B728a)>2 & (padj.Ptab>=0.01|abs(L2FC.Ptab)<=2|is.na(padj.Ptab)==T))
write.csv(DESeq_GO.B728a, file="Annotated_Out/Nb101_DESeq_B728a_only_DE_Annotated_GO.csv")
DESeq_GO.Ptab  <- filter(DESeq_GO, padj.Ptab<0.01  & abs(L2FC.Ptab)>2 & (padj.B728a>=0.01|abs(L2FC.B728a)<=2|is.na(padj.B728a)==T))
write.csv(DESeq_GO.Ptab, file="Annotated_Out/Nb101_DESeq_Ptab_only_DE_Annotated_GO.csv")

#two comparisons, exclude shared with other comparisons
DESeq_GO.Time <- filter(DESeq_GO, padj.B728a<0.01 & padj.Ptab<0.01 & abs(L2FC.B728a)>2 & abs(L2FC.Ptab)>2)
write.csv(DESeq_GO.Time, file="Annotated_Out/Nb101_DESeq_Time_Annotated_GO.csv")

# Creating a Venn Diagram -------------------------------------------------
#reading in data
NbComboB <- read.csv("DEseq_NbComboB_out.csv")
NbComboT <- read.csv("DEseq_NbComboT_out.csv")
#filtering for DEGs
NbComboB.degs <- filter(NbComboB, padj<0.01 & abs(log2FoldChange)>2)[,1] 
NbComboT.degs <- filter(NbComboT, padj<0.01 & abs(log2FoldChange)>2)[,1]
#drawing the plot
venn.plot <- venn.diagram(x=list("Psy"= NbComboB.degs, "Pta"= NbComboT.degs),
                          filename = "Plots/Venn_0vs5_pretty_correct.tiff",
                          imagetype="tiff",
                          col = c("#3B5284", "#94B447"),
                          lwd =0.5,
                          cex = 0.5,
                          fontfamily = "sans",
                          height = 720, 
                          width = 720, 
                          resolution = 300,
                          cat.cex = 0.4,
                          cat.fontface = "bold",
                          cat.default.pos = "outer",
                          cat.fontfamily = "sans",
                          fill = c("#5BA8A0", "#CBe54E"),
                          cat.col = c("#3B5284", "#94B447"),
                          cat.dist = c(0.055, 0.055),
                          margin=0.1
                          )

# GOSeq Analysis ----------------------------------------------------------
#creating lists for each type of GO term
Nb101_GO <- read.csv("Prepped/Nb101_GO.csv", header=T, fill=T)
colnames(Nb101_GO)[1]<- "Name"
list.P <- list()
for(i in 1:length(Nb101_GO[,1])){
   list.P[Nb101_GO$Name[i]] <- strsplit(as.character(Nb101_GO$GO.P.ID[i]), ";", fixed = T)
}
list.F <- list()
for(i in 1:length(Nb101_GO[,1])){
        list.F[Nb101_GO$Name[i]] <- strsplit(as.character(Nb101_GO$GO.F.ID[i]), ";", fixed = T)
}
list.C <- list()
for(i in 1:length(Nb101_GO[,1])){
        list.C[Nb101_GO$Name[i]] <- strsplit(as.character(Nb101_GO$GO.C.ID[i]), ";", fixed = T)
}

#figuring out transcript lengths by adding up the CDSs in the gff file
Nb_gff <- gffRead("prepped/Nb101.gff")
Nb_gff_att <- Nb_gff[9]
current_id <- c("Name=NbD000418.1 CDS") #have to start with the first name attribute - manually pull from file
j=1 #just a counter since the rows wont match up between input and output data frame
transc_length=0
LengthDF <- as.data.frame(matrix(ncol=2, nrow=0, dimnames=list(NULL,c("Name","Length")))) #making a new, empty data frame
for(i in 1:length(Nb_gff[,1])){
        teststr <- unlist(strsplit(as.character(Nb_gff_att[i,1]), ";", fixed = T)) #split the attribute column up
        new_id <- teststr[1]#pull the Name from the attribute string
        if(new_id%in%current_id==T){ #see if it matches the row before
                transc_length <- transc_length+Nb_gff$end[i]-Nb_gff$start[i]+1 #add the transcript length
        }else{ 
                LengthDF[j,] <- NA #create a row
                LengthDF$Name[j] <- gsub("Name=(NbD........) CDS", "\\1", current_id) #log the current ID name
                LengthDF$Length[j] <- transc_length #log the final transcript length
                current_id <- new_id #set up new_id as current id
                transc_length <- Nb_gff$end[i]-Nb_gff$start[i]+1 #start a new transcript
                j=j+1 #count up one
        }
}
LengthDF[j,] <- NA #need to have this one more time outside of the loop to catch the last one
LengthDF$Name[j] <- gsub("Name=(NbD........) CDS", "\\1", current_id) 
LengthDF$Length[j] <- transc_length
LengthDF <- LengthDF[order(LengthDF$Name),]
write.csv(LengthDF, file="Nb_transc_lengths.csv")

#need to subset for ones that had >0 base mean and ones that came out of DE file (some were lost)
Nb_DE <- read.csv(file="Annotated_out/Nb101_DESeq_All_Annotated_GO.csv", header=T)
Nb_DE  <- filter(Nb_DE, baseMean.Ptab>0)
Nb_DE  <- Nb_DE[,-1]
write.csv(Nb_DE, file="Annotated_out/Nb101_DESeq_All_Annotated_GO_mod.csv")

#Getting the DE genes and vectors for GOSeq
LengthDF <- read.csv("Nb_transc_lengths.csv", header=T)[,2:3]
allGenes <- Nb_DE[,1]

Length_DF_filter<- filter(LengthDF, Name%in%allGenes)
geneName <- Length_DF_filter$Name
geneLength <- Length_DF_filter$Length #creating the vector needed for GOseq

DEGenes.TimeOnly <- read.csv(file="Annotated_out/Nb101_DESeq_Time_Annotated_GO.csv", header=T)
DEGenes.TimeOnly <- DEGenes.TimeOnly$X
gene.vector.TimeOnly=as.integer(allGenes%in%DEGenes.TimeOnly)
names(gene.vector.TimeOnly)=allGenes

DEGenes.B728aOnly <- read.csv(file="Annotated_out/Nb101_DESeq_B728a_only_DE_Annotated_GO.csv", header=T)
DEGenes.B728aOnly <- DEGenes.B728aOnly$X
gene.vector.B728aOnly=as.integer(allGenes%in%DEGenes.B728aOnly)
names(gene.vector.B728aOnly)=allGenes

DEGenes.PtabOnly <- read.csv(file="Annotated_out/Nb101_DESeq_Ptab_only_DE_Annotated_GO.csv", header=T)
DEGenes.PtabOnly <- DEGenes.PtabOnly$X
gene.vector.PtabOnly=as.integer(allGenes%in%DEGenes.PtabOnly)
names(gene.vector.PtabOnly)=allGenes 

#Actual analysis now that I have all the inputs correct and graphing
pwf.TimeOnlyC=nullp(gene.vector.TimeOnly, bias.data = geneLength)
GO.wall.TimeOnlyC <- goseq(pwf.TimeOnlyC, gene2cat=list.C, test.cats=c("GO:CC"))
GO.wall.TimeOnlyC$padj_over=p.adjust(GO.wall.TimeOnlyC$over_represented_pvalue, method="BH")
write.csv(GO.wall.TimeOnlyC, file= "Annotated_out/GO/GO_CC_Nb_TimeOnly.csv")

pwf.TimeOnlyM=nullp(gene.vector.TimeOnly, bias.data = geneLength)
GO.wall.TimeOnlyM <- goseq(pwf.TimeOnlyM, gene2cat=list.F, test.cats=c("GO:MF"))
GO.wall.TimeOnlyM$padj_over=p.adjust(GO.wall.TimeOnlyM$over_represented_pvalue, method="BH")
write.csv(GO.wall.TimeOnlyM, file= "Annotated_out/GO/GO_MF_Nb_TimeOnly.csv")

pwf.TimeOnlyB=nullp(gene.vector.TimeOnly, bias.data = geneLength)
GO.wall.TimeOnlyB <- goseq(pwf.TimeOnlyB, gene2cat=list.P, test.cats=c("GO:BP"))
GO.wall.TimeOnlyB$padj_over=p.adjust(GO.wall.TimeOnlyB$over_represented_pvalue, method="BH")
write.csv(GO.wall.TimeOnlyB, file= "Annotated_out/GO/GO_BP_Nb_TimeOnly.csv")

GO.wall.TimeOnly <- rbind(GO.wall.TimeOnlyC, GO.wall.TimeOnlyM, GO.wall.TimeOnlyB) #combining sig terms
GO.wall.TimeOnly <- filter(GO.wall.TimeOnly, padj_over<0.05)
GO.wall.TimeOnly <- filter(GO.wall.TimeOnly, is.na(ontology)==F) #remove NA obsolete category
write.csv(GO.wall.TimeOnly, file= "Annotated_out/GO/GO_allSig_Nb_TimeOnly.csv")

pwf.B728aOnlyC=nullp(gene.vector.B728aOnly, bias.data = geneLength)
GO.wall.B728aOnlyC <- goseq(pwf.B728aOnlyC, gene2cat=list.C, test.cats=c("GO:CC"))
GO.wall.B728aOnlyC$padj_over=p.adjust(GO.wall.B728aOnlyC$over_represented_pvalue, method="BH")
write.csv(GO.wall.B728aOnlyC, file= "Annotated_out/GO/GO_CC_Nb_B728aOnly.csv")

pwf.B728aOnlyM=nullp(gene.vector.B728aOnly, bias.data = geneLength)
GO.wall.B728aOnlyM <- goseq(pwf.B728aOnlyM, gene2cat=list.F, test.cats=c("GO:MF"))
GO.wall.B728aOnlyM$padj_over=p.adjust(GO.wall.B728aOnlyM$over_represented_pvalue, method="BH")
write.csv(GO.wall.B728aOnlyM, file= "Annotated_out/GO/GO_MF_Nb_B728aOnly.csv")

pwf.B728aOnlyB=nullp(gene.vector.B728aOnly, bias.data = geneLength)
GO.wall.B728aOnlyB <- goseq(pwf.B728aOnlyB, gene2cat=list.P, test.cats=c("GO:BP"))
GO.wall.B728aOnlyB$padj_over=p.adjust(GO.wall.B728aOnlyB$over_represented_pvalue, method="BH")
write.csv(GO.wall.B728aOnlyB, file= "Annotated_out/GO/GO_BP_Nb_B728aOnly.csv")

GO.wall.B728aOnly <- rbind(GO.wall.B728aOnlyC, GO.wall.B728aOnlyM, GO.wall.B728aOnlyB) #combining sig terms
GO.wall.B728aOnly <- filter(GO.wall.B728aOnly, padj_over<0.05)
write.csv(GO.wall.B728aOnly, file= "Annotated_out/GO/GO_allSig_Nb_B728aOnly.csv")

pwf.PtabOnlyC=nullp(gene.vector.PtabOnly, bias.data = geneLength)
GO.wall.PtabOnlyC <- goseq(pwf.PtabOnlyC, gene2cat=list.C, test.cats=c("GO:CC"))
GO.wall.PtabOnlyC$padj_over=p.adjust(GO.wall.PtabOnlyC$over_represented_pvalue, method="BH")
write.csv(GO.wall.PtabOnlyC, file= "Annotated_out/GO/GO_CC_Nb_PtabOnly.csv")

pwf.PtabOnlyM=nullp(gene.vector.PtabOnly, bias.data = geneLength)
GO.wall.PtabOnlyM <- goseq(pwf.PtabOnlyM, gene2cat=list.F, test.cats=c("GO:MF"))
GO.wall.PtabOnlyM$padj_over=p.adjust(GO.wall.PtabOnlyM$over_represented_pvalue, method="BH")
write.csv(GO.wall.PtabOnlyM, file= "Annotated_out/GO/GO_MF_Nb_PtabOnly.csv")

pwf.PtabOnlyB=nullp(gene.vector.PtabOnly, bias.data = geneLength)
GO.wall.PtabOnlyB <- goseq(pwf.PtabOnlyB, gene2cat=list.P, test.cats=c("GO:BP"))
GO.wall.PtabOnlyB$padj_over=p.adjust(GO.wall.PtabOnlyB$over_represented_pvalue, method="BH")
write.csv(GO.wall.PtabOnlyB, file= "Annotated_out/GO/GO_BP_Nb_PtabOnly.csv")

GO.wall.PtabOnly <- rbind(GO.wall.PtabOnlyC, GO.wall.PtabOnlyM, GO.wall.PtabOnlyB) #combining sig terms
GO.wall.PtabOnly <- filter(GO.wall.PtabOnly, padj_over<0.05)
write.csv(GO.wall.PtabOnly, file= "Annotated_out/GO/GO_allSig_Nb_PtabOnly.csv")

# Comparison GO term Visualization ---------------------------------------------
GO.wall.B728aOnly <- read.csv(file= "Annotated_out/GO/GO_allSig_Nb_B728aOnly.csv", header=T)
GO.wall.B728aOnly$categoryAll <- paste(GO.wall.B728aOnly$category, GO.wall.B728aOnly$term) #setting up for graphing
GO.wall.B728aOnly$categoryAll <- substr(GO.wall.B728aOnly$categoryAll, 1, 40)
GO.wall.B728aOnly$categoryAll <- factor(GO.wall.B728aOnly$categoryAll, levels=rev(GO.wall.B728aOnly$categoryAll))

GO.wall.PtabOnly <- read.csv(file= "Annotated_out/GO/GO_allSig_Nb_PtabOnly.csv", header=T)
GO.wall.PtabOnly$categoryAll <- paste(GO.wall.PtabOnly$category, GO.wall.PtabOnly$term) #setting up for graphing
GO.wall.PtabOnly$categoryAll <- substr(GO.wall.PtabOnly$categoryAll, 1, 40)
GO.wall.PtabOnly$categoryAll <- factor(GO.wall.PtabOnly$categoryAll, levels=rev(GO.wall.PtabOnly$categoryAll))

GO.wall.All <- rbind(GO.wall.B728aOnly, GO.wall.PtabOnly)
GO.wall.All$categoryAll <- factor(GO.wall.All$categoryAll, levels=rev(GO.wall.All$categoryAll))

DEGenes.B728aOnly <- read.csv(file="Annotated_out/Nb101_DESeq_B728a_only_DE_Annotated_GO.csv", header=T)[,-1:-2]
DEGenes.B728aOnly$AllGO <- paste(DEGenes.B728aOnly$GO.P.ID, DEGenes.B728aOnly$GO.F.ID, DEGenes.B728aOnly$GO.C.ID, sep=";")
DEGenes.PtabOnly <- read.csv(file="Annotated_out/Nb101_DESeq_Ptab_only_DE_Annotated_GO.csv", header=T)[,-1:-2]
DEGenes.PtabOnly$AllGO <- paste(DEGenes.PtabOnly$GO.P.ID, DEGenes.PtabOnly$GO.F.ID, DEGenes.PtabOnly$GO.C.ID, sep=";")

GO_terms_int <- c(GO.wall.All$category) #create vector of the enriched GO terms
GO_pd_allBB <- as.data.frame(matrix(ncol=4, nrow=0)) #create new data frame to fill in
for (GO in GO_terms_int) { #pulls genes that have enriched GO terms
        GO_genes <-subset(DEGenes.B728aOnly, grepl(GO, DEGenes.B728aOnly$AllGO)) #subset for ones with each GO term
        GO_plot_data <- GO_genes[, c("X", "L2FC.B728a")] #subset to just two columns
        GO_plot_data$term <- GO #add in the GO term they were pulled for
        GO_plot_data$ontology <- GO.wall.All$ontology[which(GO_terms_int %in% GO)] #add ontology category
        GO_pd_allBB <- rbind(GO_pd_allBB, GO_plot_data) #add to dataframe for this GO term
}
GO_pd_allBB$term <- factor(GO_pd_allBB$term, levels=rev(GO.wall.All$category)) #changing factor level to match other graph

GO_terms_int <- c(GO.wall.All$category) #create vector of the enriched GO terms
GO_pd_allBT <- as.data.frame(matrix(ncol=4, nrow=0)) #create new data frame to fill in
for (GO in GO_terms_int) { #pulls genes that have enriched GO terms
        GO_genes <-subset(DEGenes.PtabOnly, grepl(GO, DEGenes.PtabOnly$AllGO)) #subset for ones with each GO term
        if(dim(GO_genes)[1]!=0){
                GO_plot_data <- GO_genes[, c("X", "L2FC.Ptab")] #subset to just two columns
                GO_plot_data$term <- GO #add in the GO term they were pulled for
                GO_plot_data$ontology <- GO.wall.All$ontology[which(GO_terms_int %in% GO)] #add ontology category
                GO_pd_allBT <- rbind(GO_pd_allBT, GO_plot_data) #add to dataframe for this GO term
         }else{ #this is just creating a row with a nonsense L2FC so that the category doesn't get skipped
                 new_df <- as.data.frame(matrix(ncol=4, nrow=1))
                 colnames(new_df) <- c("X", "L2FC.Ptab", "term", "ontology")
                 new_df[1,2] <- 50
                 new_df[1,3] <- GO #add in the GO term they were pulled for
                 new_df[1,4] <- GO.wall.All$ontology[which(GO_terms_int %in% GO)] #add ontology category
                 GO_pd_allBT <- rbind(GO_pd_allBT, new_df)
       }
}
GO_pd_allBT$term <- factor(GO_pd_allBT$term, levels=rev(GO.wall.All$category))

enrich_plot <- ggplot(GO.wall.All, aes(x=categoryAll, y=-log10(padj_over), fill=ontology)) +
        geom_rect(aes(ymin=-Inf, ymax=Inf, xmin="GO:0016117 carotenoid biosynthetic proce", xmax=Inf), 
                  fill="#d3daeb",
                  alpha=0.02) + #adding in a rectangle to shade the in -2 to 2
        geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-Inf, xmax="GO:0003700 DNA-binding transcription fac"), 
                  fill="#dde8c4",
                  alpha=0.02) +
        geom_bar(stat="identity") +
        scale_fill_manual(name="GO Category", 
                          values=c("#feb078", "#b73779", "#2c115f"), 
                          breaks=c('CC', 'MF', 'BP')) +
        xlab("GO Term") +
        ylab("Enrichment\nScore") +
        #ggtitle("Enriched GO Terms (Ptab)") +
        ylim(0,15) +
        coord_flip() +
        theme_light() +
        theme(
                #plot.title=element_text(angle=0, size=14, face="bold", vjust=1, hjust=0.5),
                axis.text.x=element_text(angle=0, size=8),
                axis.text.y=element_text(angle=0, size=8),
                axis.title=element_text(size=10, face="bold"),
                title=element_text(size=10),
                legend.background=element_rect(),
                legend.key=element_blank(),     #removes the border
                legend.text=element_text(size=8),
                legend.title=element_text(size=10, face="bold"),
                legend.position = c(-1, -0.06), 
                legend.direction = "horizontal",
                legend.title.align = 0) +
        annotate("text", x = "GO:0042644 chloroplast nucleoid", y = 14, label = "Psy-Response Enriched", angle=90, color="#3B5284", fontface="bold") +
        annotate("text", x = "GO:0010182 sugar mediated signaling path", y = 14, label = "Pta-Response Enriched", angle=90, color="#94B447", fontface="bold")


gene_plotB <- ggplot(GO_pd_allBB, aes(y=L2FC.B728a, x=term)) + #creating the plot
        geom_rect(aes(ymin=-2, ymax=2, xmin=-Inf, xmax=Inf), fill="grey86") + #adding in a rectangle to shade the in -2 to 2
        geom_rect(aes(ymin=-Inf, ymax=Inf, xmin="GO:0016117", xmax=Inf), 
                  fill="#d3daeb",
                  alpha=0.02) + #adding in a rectangle to shade the in -2 to 2
        geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-Inf, xmax="GO:0003700"), 
                  fill="#dde8c4",
                  alpha=0.02) +
        geom_jitter(aes(color = ontology), 
                    width=0.2,
                    alpha = .5,
                    show.legend = FALSE)+
        scale_color_manual(name="GO Category",
                           values=c("#feb078", "#b73779", "#2c115f"), 
                           breaks=c('CC', 'MF', 'BP')) +
        ylim(-12, 27) +
        ylab("L2FC of Psy \nDE Genes") +
        #ggtitle("Corresponding Genes") +
        coord_flip() +
        theme_light() +
        theme(
                plot.title=element_text(angle=0, size=12, face="bold", vjust=1, hjust=0.5),
                axis.text.y=element_blank(), #get rid of y axis
                axis.title.y=element_blank(), #get rid of y axis
                axis.ticks.y=element_blank(), #get rid of y axis
                axis.text.x=element_text(angle=0, size=8),
                axis.title=element_text(size=10, face="bold"),
                title=element_text(size=10),
        ) +
        guides(color = guide_legend(override.aes = list(size=5, alpha=1))) + #resizing the legend dots
        geom_hline(yintercept = c(-2, 2), linetype = "dotted", lwd=0.5) #adding a line to 2 and -2 #resizing the legend dots

gene_plotT <- ggplot(GO_pd_allBT, aes(y=L2FC.Ptab, x=term)) + #creating the plot
        geom_rect(aes(ymin=-2, ymax=2, xmin=-Inf, xmax=Inf), fill="grey86") + #adding in a rectangle to shade the in -2 to 2
        geom_rect(aes(ymin=-Inf, ymax=Inf, xmin="GO:0016117", xmax=Inf), 
                  fill="#d3daeb",
                  alpha=0.02) + #adding in a rectangle to shade the in -2 to 2
        geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-Inf, xmax="GO:0003700"), 
                  fill="#dde8c4",
                  alpha=0.02) +
        geom_jitter(aes(color = ontology), 
                    width=0.2,
                    alpha = .5,
                    show.legend = FALSE)+
        scale_color_manual(name="GO Category",
                           values=c("#feb078", "#b73779", "#2c115f"), 
                           breaks=c('CC', 'MF', 'BP')) +
        ylim(-12, 27) +
        ylab("L2FC of Pta \nDE Genes") +
        #ggtitle("Corresponding Genes") +
        coord_flip() +
        theme_light() +
        theme(
                plot.title=element_text(angle=0, size=12, face="bold", vjust=1, hjust=0.5),
                axis.text.y=element_blank(), #get rid of y axis
                axis.title.y=element_blank(), #get rid of y axis
                axis.ticks.y=element_blank(), #get rid of y axis
                axis.text.x=element_text(angle=0, size=8),
                axis.title=element_text(size=10, face="bold"),
                title=element_text(size=10),
        ) +
        guides(color = guide_legend(override.aes = list(size=5, alpha=1))) + #resizing the legend dots
        geom_hline(yintercept = c(-2, 2), linetype = "dotted", lwd=0.5) #adding a line to 2 and -2 #resizing the legend dots

total_plot <- enrich_plot+gene_plotB+gene_plotT
total_plot

ggsave("Plots/Figure4.tiff",
       device="tiff",
       plot=total_plot, 
       height = 135, 
       width = 172, 
       units="mm")
ggsave("Plots/Figure4.png",
       device="png",
       plot=total_plot, 
       height = 135, 
       width = 172, 
       units="mm")

# Overlap GO term visualization ----------------------------------------------
#All GO terms
GO.wall.Time <- read.csv(file= "Annotated_out/GO/GO_allSig_Nb_TimeOnly.csv", header=T)
DEGenes.Time <- read.csv(file="Annotated_out/Nb101_DESeq_Time_Annotated_GO.csv", header=T)[,-1:-2]
DEGenes.Time$AllGO <- paste(DEGenes.Time$GO.P.ID, DEGenes.Time$GO.F.ID, DEGenes.Time$GO.C.ID, sep=";")

GO.wall.Time$category <- paste(GO.wall.Time$category, GO.wall.Time$term) #setting up for graphing
GO.wall.Time$category <- substr(GO.wall.Time$category, 1, 40)
GO.wall.Time$category <- factor(GO.wall.Time$category, levels=rev(GO.wall.Time$category)) 

enrich_plot <- ggplot(GO.wall.Time, aes(x=category, y=-log10(padj_over), fill=ontology)) +
        geom_bar(stat="identity",
        )+
        scale_fill_manual(name="GO Category",
                          values=c("#feb078", "#b73779", "#2c115f"), 
                          breaks=c('CC', 'MF', 'BP')) +
        xlab("GO Term") +
        ylab("Enrichment Score") +
        theme_light()+
        theme(
                legend.background=element_rect(),
                #plot.title=element_text(angle=0, size=14, face="bold", vjust=1, hjust=0.5),
                axis.text.x=element_text(angle=0, size=8),
                axis.text.y=element_text(angle=0, size=6),
                axis.title=element_text(size=10, face="bold"),
                legend.key=element_blank(),     #removes the border
                legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
                legend.text=element_text(size=10),  #Text size
                title=element_text(size=12)) +
        coord_flip()

ggsave("Plots/FigureS2.tiff",
       device="tiff",
       plot=enrich_plot, 
       height = 240, 
       width = 172, 
       units="mm")
ggsave("Plots/FigureS2.png",
       device="png",
       plot=enrich_plot, 
       height = 240, 
       width = 172, 
       units="mm")

#BP only
GO.wall.TimeOnlyBP <- read.csv(file= "Annotated_out/GO/GO_BP_Nb_TimeOnly.csv", header=T)

GO.wall.TimeOnlyBP <- filter(GO.wall.TimeOnlyBP, padj_over<0.01)

DEGenes.Time <- read.csv(file="Annotated_out/Nb101_DESeq_Time_Annotated_GO.csv", header=T)[,-1:-2]
DEGenes.Time$AllGO <- paste(DEGenes.Time$GO.P.ID, DEGenes.Time$GO.F.ID, DEGenes.Time$GO.C.ID, sep=";")

GO_terms_int <- c(GO.wall.TimeOnlyBP$category) #create vector of the enriched GO terms
GO_pd_allBB <- as.data.frame(matrix(ncol=4, nrow=0)) #create new data frame to fill in
for (GO in GO_terms_int) { #pulls genes that have enriched GO terms
        GO_genes <-subset(DEGenes.Time, grepl(GO, DEGenes.Time$AllGO)) #subset for ones with each GO term
        if(dim(GO_genes)[1]!=0){
                GO_plot_data <- GO_genes[, c("X", "L2FC.B728a")] #subset to just two columns
                GO_plot_data$term <- GO #add in the GO term they were pulled for
                GO_plot_data$ontology <- GO.wall.TimeOnlyBP$ontology[which(GO_terms_int %in% GO)] #add ontology category
                GO_pd_allBB <- rbind(GO_pd_allBB, GO_plot_data) #add to dataframe for this GO term
        }else{ #this is just creating a row with a nonsense L2FC so that the category doesn't get skipped
                new_df <- as.data.frame(matrix(ncol=4, nrow=1))
                colnames(new_df) <- c("X", "L2FC.Ptab", "term", "ontology")
                new_df[1,2] <- 50
                new_df[1,3] <- GO #add in the GO term they were pulled for
                new_df[1,4] <- GO.wall.TimeOnlyBP$ontology[which(GO_terms_int %in% GO)] #add ontology category
                GO_pd_allBT <- rbind(GO_pd_allBT, new_df)
        }
}
GO_pd_allBB$term <- factor(GO_pd_allBB$term, levels=rev(GO.wall.TimeOnlyBP$category)) #changing factor level to match other graph

GO_terms_int <- c(GO.wall.TimeOnlyBP$category) #create vector of the enriched GO terms
GO_pd_allBT <- as.data.frame(matrix(ncol=4, nrow=0)) #create new data frame to fill in
for (GO in GO_terms_int) { #pulls genes that have enriched GO terms
        GO_genes <-subset(DEGenes.Time, grepl(GO, DEGenes.Time$AllGO)) #subset for ones with each GO term
        if(dim(GO_genes)[1]!=0){
                GO_plot_data <- GO_genes[, c("X", "L2FC.Ptab")] #subset to just two columns
                GO_plot_data$term <- GO #add in the GO term they were pulled for
                GO_plot_data$ontology <- GO.wall.TimeOnlyBP$ontology[which(GO_terms_int %in% GO)] #add ontology category
                GO_pd_allBT <- rbind(GO_pd_allBT, GO_plot_data) #add to dataframe for this GO term
        }else{ #this is just creating a row with a nonsense L2FC so that the category doesn't get skipped
                new_df <- as.data.frame(matrix(ncol=4, nrow=1))
                colnames(new_df) <- c("X", "L2FC.Ptab", "term", "ontology")
                new_df[1,2] <- 50
                new_df[1,3] <- GO #add in the GO term they were pulled for
                new_df[1,4] <- GO.wall.TimeOnlyBP$ontology[which(GO_terms_int %in% GO)] #add ontology category
                GO_pd_allBT <- rbind(GO_pd_allBT, new_df)
        }
}
GO_pd_allBT$term <- factor(GO_pd_allBT$term, levels=rev(GO.wall.TimeOnlyBP$category))

GO.wall.TimeOnlyBP$category <- paste(GO.wall.TimeOnlyBP$category, GO.wall.TimeOnlyBP$term) #setting up for graphing
GO.wall.TimeOnlyBP$category <- substr(GO.wall.TimeOnlyBP$category, 1, 40)
GO.wall.TimeOnlyBP$category <- factor(GO.wall.TimeOnlyBP$category, levels=rev(GO.wall.TimeOnlyBP$category)) 

enrich_plot <- ggplot(GO.wall.TimeOnlyBP, aes(x=category, y=-log10(padj_over), fill=ontology)) +
        geom_bar(stat="identity",
                 show.legend = FALSE
                 )+
        scale_fill_manual(values=c("#feb078", "#b73779", "#2c115f"), breaks=c('CC', 'MF', 'BP')) +
        xlab("GO Term") +
        ylab("Enrichment\nScore") +
        theme_light()+
        theme(
                legend.background=element_rect(),
                axis.text.x=element_text(angle=0, size=8),
                axis.text.y=element_text(angle=0, size=8),
                axis.title=element_text(size=10, face="bold"),
                title=element_text(size=12)) +
        coord_flip()

gene_plotB <- ggplot(GO_pd_allBB, aes(y=L2FC.B728a, x=term)) + #creating the plot
        geom_rect(aes(ymin=-2, ymax=2, xmin=-Inf, xmax=Inf), fill="grey86") + #adding in a rectangle to shade the in -2 to 2
        geom_jitter(aes(color = ontology), 
                    width=0.2,
                    alpha = .5,
                    show.legend = FALSE)+
        scale_color_manual(name="GO Category",
                           values=c("#feb078", "#b73779", "#2c115f"), 
                           breaks=c('CC', 'MF', 'BP')) +
        ylim(-25, 16) +
        ylab("L2FC of Psy \nDE Genes") +
        coord_flip() +
        theme_light() +
        theme(
                plot.title=element_text(angle=0, size=12, face="bold", vjust=1, hjust=0.5),
                axis.text.y=element_blank(), #get rid of y axis
                axis.title.y=element_blank(), #get rid of y axis
                axis.ticks.y=element_blank(), #get rid of y axis
                axis.text.x=element_text(angle=0, size=8),
                axis.title=element_text(size=10, face="bold"),
                title=element_text(size=12),
        ) +
        guides(color = guide_legend(override.aes = list(size=5, alpha=1))) + #resizing the legend dots
        geom_hline(yintercept = c(-2, 2), linetype = "dotted", lwd=0.5) #adding a line to 2 and -2 #resizing the legend dots

gene_plotT <- ggplot(GO_pd_allBT, aes(y=L2FC.Ptab, x=term)) + #creating the plot
        geom_rect(aes(ymin=-2, ymax=2, xmin=-Inf, xmax=Inf), fill="grey86") + #adding in a rectangle to shade the in -2 to 2
        geom_jitter(aes(color = ontology), 
                    width=0.2,
                    alpha = .5,
                        show.legend = FALSE)+
        scale_color_manual(name="GO Category",
                           values=c("#feb078", "#b73779", "#2c115f"), 
                           breaks=c('CC', 'MF', 'BP')) +
        ylim(-25, 16) +
        ylab("L2FC of Pta \nDE Genes") +
        coord_flip() +
        theme_light() +
        theme(
                plot.title=element_text(angle=0, size=12, face="bold", vjust=1, hjust=0.5),
                axis.text.y=element_blank(), #get rid of y axis
                axis.title.y=element_blank(), #get rid of y axis
                axis.ticks.y=element_blank(), #get rid of y axis
                axis.text.x=element_text(angle=0, size=8),
                axis.title=element_text(size=10, face="bold"),
                title=element_text(size=10),
        ) +
        guides(color = guide_legend(override.aes = list(size=5, alpha=1))) + #resizing the legend dots
        geom_hline(yintercept = c(-2, 2), linetype = "dotted", lwd=0.5) #adding a line to 2 and -2 #resizing the legend dots

total_plot <- enrich_plot+gene_plotB+gene_plotT
total_plot
ggsave("Plots/Figure2.tiff",
       device="tiff",
       plot=total_plot, 
       height = 135, 
       width = 172, 
       units="mm")
ggsave("Plots/Figure2.png",
       device="png",
       plot=total_plot, 
       height = 135, 
       width = 172, 
       units="mm")

# L2FC visualization ------------------------------------------------------
DESeq_GO<- read.csv("Annotated_Out/Nb101_DESeq_All_Annotated_GO_mod.csv", header=T)
DESeq_GO <- DESeq_GO[,-1]
colnames(DESeq_GO)[1]<- "X" 
allGenes <- DESeq_GO[,1]

#one comparison only, exclude shared with other comparisons (for two I had to include is.na for padj because the numbers werent lining up with the venn diagram)
DESeq_GO.B728a <- filter(DESeq_GO, padj.B728a<0.01 & abs(L2FC.B728a)>2 & (padj.Ptab>=0.01|abs(L2FC.Ptab)<=2|is.na(padj.Ptab)==T))
DESeq_GO.Ptab  <- filter(DESeq_GO, padj.Ptab<0.01  & abs(L2FC.Ptab)>2 & (padj.B728a>=0.01|abs(L2FC.B728a)<=2|is.na(padj.B728a)==T))
#two comparisons, exclude shared with other comparisons
DESeq_GO.Time <- filter(DESeq_GO, padj.B728a<0.01 & padj.Ptab<0.01 & abs(L2FC.B728a)>2 & abs(L2FC.Ptab)>2)

#make a file with relevant info for plot
plot_data <- DESeq_GO[, c("X", "L2FC.B728a", "L2FC.Ptab")] #subset to just two columns
DEGenes.TimeOnly <- DESeq_GO.Time$X
gene.vector.TimeOnly=as.integer(allGenes%in%DEGenes.TimeOnly)
plot_data$BothSig <- gene.vector.TimeOnly
DEGenes.Ptab <- DESeq_GO.Ptab$X
gene.vector.Ptab=as.integer(allGenes%in%DEGenes.Ptab)
plot_data$PtabSig <- gene.vector.Ptab
DEGenes.B728a <- DESeq_GO.B728a$X
gene.vector.B728a=as.integer(allGenes%in%DEGenes.B728a)
plot_data$B728aSig <- gene.vector.B728a
plot_data$Sig <- as.character(paste(plot_data$BothSig, plot_data$PtabSig, plot_data$B728aSig, sep=""))
plot_data2 <- subset(plot_data, plot_data$Sig!="000") #only the genes that are significant


L2FC_plot <- ggplot(plot_data2, aes(y=L2FC.B728a, x=L2FC.Ptab)) + #creating the plot
        geom_rect(aes(ymin=-2, ymax=2, xmin=-Inf, xmax=Inf), fill="grey96") + #adding in a rectangle to shade the in -2 to 2
        geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=-2, xmax=2), fill="grey96") + #adding in a rectangle to shade the in -2 to 2
        geom_jitter(aes(color=Sig), 
                    width=0.2,
                    alpha = .5)+
        scale_color_manual(name="DE in Response to",
                           values=c("#b73779", "#94B447", "#3B5284"),
                           breaks=c('100', '010', '001'),
                           labels = c("Both", "Pta", "Psy")) +
        ylab("L2FC for Response to Psy") +
        xlab("L2FC for Response to Pta") +
        ylim(-29, 29)+
        xlim(-33, 33)+
        theme_light() +
        theme(
                plot.title=element_text(angle=0, size=12, face="bold", vjust=1, hjust=0.5),
                axis.text.y=element_text(angle=0, size=10),
                axis.title.y=element_text(size=10, face="bold"), 
                axis.text.x=element_text(angle=0, size=10),
                axis.title=element_text(size=10, face="bold"),
                legend.background=element_rect(),
                legend.key=element_blank(),     #removes the border
                legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
                legend.text=element_text(size=10),
                legend.title=element_text(size=10, face="bold"),
                title=element_text(size=10),
        ) +
        guides(color = guide_legend(override.aes = list(size=5, alpha=1)))+#resizing the legend dots
        geom_hline(yintercept = c(-2, 2), linetype = "dotted", lwd=0.5)+ #adding a line to 2 and -2 #resizing the legend dots
        geom_vline(xintercept = c(-2, 2), linetype = "dotted", lwd=0.5)+#adding a line to 2 and -2 #resizing the legend dots
        geom_mark_ellipse(data = plot_data2 %>% 
                                filter(L2FC.B728a>18&L2FC.Ptab<10),
                          aes(color = Sig),
                          expand=unit(2, "mm"),
                          show.legend = FALSE)+
        geom_mark_ellipse(data = plot_data2 %>% 
                                filter(L2FC.B728a<5&L2FC.Ptab>15),
                          aes(color = Sig),
                          expand=unit(2, "mm"),
                          show.legend = FALSE)+
        geom_mark_ellipse(data = plot_data2 %>% 
                                filter(L2FC.B728a < -18 & L2FC.Ptab > -3),
                          aes(color = Sig),
                          expand=unit(2, "mm"),
                          show.legend = FALSE)
        
ggsave("Plots/Figure3.tiff",
       device="tiff",
       plot=L2FC_plot, 
       height = 90, 
       width = 134, 
       units="mm")
ggsave("Plots/Figure3.png",
       device="png",
       plot=L2FC_plot, 
       height = 90, 
       width = 134, 
       units="mm")
