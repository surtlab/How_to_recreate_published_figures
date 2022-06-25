## RNA-seq of Nicotiana benthamiana responses to Pta and Psy
#Morgan Carter 2021-2022
#Anything within the prepped directory is an input file, including genomes, count data, and phenotype csv for DESeq2
#Important note is that the B728a count data relies on time zero points from JGI (accessions Gp0060698, Gp0060702, and Gp0060705) which are first added as new samples into our B728a count data output from featureCounts

# Packages and setup ------------------------------------------------------
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
library(ggpattern)
setwd("E:/Box/BaltrusLab/Ptab_B728a_RNAseq3") #Desktop


# B728a Count Data --------------------------------------------------------
R1 <- read.delim('prepped/B728a_T0_JGI/20C_R1.txt', header=T) #load in the three text files with count data
R2 <- read.delim('prepped/B728a_T0_JGI/20C_R2.txt', header=T)
R3 <- read.delim('prepped/B728a_T0_JGI/20C_R3.txt', header=T)
R1 <- R1[, c("Locus.Tag", "Reads.Count")] #subset to just counts and locus tag
R2 <- R2[, c("Locus.Tag", "Reads.Count")] 
R3 <- R3[, c("Locus.Tag", "Reads.Count")] 
B728a <- merge.data.frame(R1, R2, by="Locus.Tag") #Merging them by the locus tag
B728a <- merge.data.frame(B728a, R3, by="Locus.Tag")
colnames(B728a) <- c("Locus.Tag", "B0_B728a_1", "B0_B728a_2", "B0_B728a_3") #renaming to match other data

#need to remove first line from featureCounts output and simplify column names or next steps will yield error
count_B728a <- read.delim("prepped/featureCounts_mapped_results_B728a_4DE.txt", header=T) #read in other count data
count_B728a <- count_B728a[, c("Geneid", "B5_B728a_1", "B5_B728a_2", "B5_B728a_3")] #subset to just two columns
count_B728a <- merge.data.frame(B728a, count_B728a, by.x="Locus.Tag", by.y="Geneid") #Merging them by the locus tag

write.table(count_B728a, file="prepped/B728a_JGIB0_NbB5_count.txt", sep="\t") #writing out resulting count data which will be used for DESeq2

# DESeq2 Analysis ---------------------------------------------------------
#reading in gene count files
counts_data_B728a <- read.table("prepped/B728a_JGIB0_NbB5_count.txt", header=T, sep="")
counts_data_Ptab <- read.table("prepped/featureCounts_mapped_results_Ptab_4DE.txt", header=T, sep="\t")
counts_data_Ptab <- subset(counts_data_Ptab, select=-c(Chr, Start, End, Strand, Length))

#formatting it for DESeq count table
row.names(counts_data_B728a) <- counts_data_B728a$Locus.Tag
counts_data_B728a <-counts_data_B728a[,-1]
counts_data_B728a[,1:6] <- lapply(counts_data_B728a, function(x) ceiling(x))
row.names(counts_data_Ptab) <- counts_data_Ptab$Geneid
counts_data_Ptab <-counts_data_Ptab[,-1]
counts_data_Ptab[,1:6] <- lapply(counts_data_Ptab, function(x) ceiling(x))

#reading in phenotype table and turning integers into factors
coldata_B728a <- read.csv("prepped/B728a_Pheno.csv", row.names = "ids")
coldata_B728a[,1] <- as.factor(coldata_B728a[,1])
coldata_Ptab <- read.csv("prepped/Ptab_Pheno.csv", row.names = "ids")
coldata_Ptab[,1] <- as.factor(coldata_Ptab[,1])

#checking that rownames in the phenotype and column names in gene count files match
all(rownames(coldata_B728a) %in% colnames(counts_data_B728a))
all(rownames(coldata_B728a) == colnames(counts_data_B728a))
all(rownames(coldata_Ptab) %in% colnames(counts_data_Ptab))
all(rownames(coldata_Ptab) == colnames(counts_data_Ptab))

#the main function as outlined by DESeq manual to make data set
dds_B728a <- DESeqDataSetFromMatrix(countData = counts_data_B728a, 
                                    colData = coldata_B728a, 
                                    design = ~ time)
dds_Ptab <- DESeqDataSetFromMatrix(countData = counts_data_Ptab, 
                                   colData = coldata_Ptab, 
                                   design = ~ time)

#running DEseq
dds.dataB728a <- DESeq(dds_B728a)
dds.dataPtab <- DESeq(dds_Ptab)

#creating the output table showing all genes and whether they are DEGs
contrast_time <- c("time", "5", "0") # The name provided in the second element (in this case 0) is the level that is used as baseline

res_B728a_time<-results(dds.dataB728a,
                        contrast=contrast_time,
                        independentFiltering = T,
                        alpha = 0.01)
summary(res_B728a_time)
write.csv(as.data.frame(res_B728a_time),file="DEseq_B728aTime_out.csv")

res_Ptab_time<-results(dds.dataPtab,
                       contrast=contrast_time,
                       independentFiltering = T,
                       alpha = 0.01)
summary(res_Ptab_time)
write.csv(as.data.frame(res_Ptab_time),file="DEseq_PtabTime_out.csv")

# Creating plots ----------------------------------------------------------
#note some plots were further modified for aesthetics (text and colors) in Illustrator
rld_B728a <- rlog(dds_B728a, blind=F)
rld_Ptab <- rlog(dds_Ptab, blind=F)

plotPCA(rld_B728a, intgroup=c("time"))
plotPCA(rld_Ptab, intgroup=c("time"))

degs_namesB <-
  rownames(subset(res_B728a_time,
                  padj < 0.01 & abs(log2FoldChange)>2))
degs_namesT <-
  rownames(subset(res_Ptab_time,
                  padj < 0.01 & abs(log2FoldChange)>2))

ann_colors_Ptab = list(
  group = c("0" = "#CBe54E", "5" = "#94B447")
)
ann_colors_B728a = list(
  group = c("0" ="#5BA8A0", "5" = "#3B5284")
)

pheatmap(assay(rld_B728a)[degs_namesB,],
         show_rownames = F,
         annotation_col = coldata_B728a,
         scale="row",
         color=magma(15),
         filename = "Plots/B728a_heatmap.pdf", 
         show_colnames = F
)

pheatmap(assay(rld_Ptab)[degs_namesT,],
         show_rownames = F,
         annotation_col = coldata_Ptab,
         scale="row",
         color=magma(15),
         filename = "Plots/Ptab_heatmap.pdf", 
         show_colnames = F
)

# Adding annotations from gff file to DESeq2 output - bacteria -----------------------
##B728a file
gffdata_mRNA = gffRead('prepped/Pseudomonas_syringae_B728a_112.gff')
gffdata_mRNA <- gffdata_mRNA[gffdata_mRNA$feature == "CDS",] #limiting it to one entry per transcript
gff_stuff <- c("start", "end", "strand", "ID", "Parent", 'locus', "name") #will be column headers for all of attributes
GeneInfo <- as.data.frame(matrix(ncol=length(gff_stuff), nrow=0, dimnames=list(NULL,gff_stuff))) #making a new, empty data frame
gffdata_gene_mod <- gffdata_mRNA[9]
num <- dim(gffdata_gene_mod)
for (i in 1:num[1]) { #for every line in the gff file
  GeneInfo[i,] <- NA
  GeneInfo$start[i] <- gffdata_mRNA$start[i]
  GeneInfo$end[i] <- gffdata_mRNA$end[i]
  GeneInfo$strand[i] <- gffdata_mRNA$strand[i]
  teststr <- unlist(strsplit(as.character(gffdata_gene_mod[i,1]), ";", fixed = T)) #split up that attribute column
  for (j in 1:length(teststr)) {
    teststr2 <-unlist(strsplit(teststr[j], "=", fixed = T))
    GeneInfo[i, j+3] <- teststr2[2]
  }
}
row.names(GeneInfo) <- gsub("gene", "", GeneInfo$Parent)
write.csv(GeneInfo, file="B728a_gff_info.csv")

DE_data <- read.csv("DEseq_B728aTime_out.csv", header=T)
gene_info <- read.csv("B728a_gff_info.csv", header=T)
gene_info <- gene_info[,c("locus","start","end","strand","name")]
DeSeq_Info_B728a <- merge.data.frame(gene_info, DE_data, by.x="locus", by.y="X")
write.csv(DeSeq_Info_B728a, file="Annotated_Out/B728a_DESeq_Time_Annotated.csv")

##Ptab file
gffdata_mRNA = gffRead('prepped/Ptab_ncbi_all.gff3') #had to delete the chromosome sequence from the gff file Dave provided
gffdata_mRNA <- gffdata_mRNA[gffdata_mRNA$feature == "CDS",] #limiting it to one entry per transcript
gff_stuff <- c("start", "end", "strand", "Parent", "EC_number", "locus_tag", "codon_start", "translation", "transl_table", "inference", "product", "note") #will be column headers for all of attributes
GeneInfo <- as.data.frame(matrix(ncol=length(gff_stuff), nrow=0, dimnames=list(NULL,gff_stuff))) #making a new, empty data frame
gffdata_gene_mod <- gffdata_mRNA[9]
num <- dim(gffdata_gene_mod)
for (i in 1:num[1]) { #for every line in the gff file
  GeneInfo[i,] <- NA
  GeneInfo$start[i] <- gffdata_mRNA$start[i]
  GeneInfo$end[i] <- gffdata_mRNA$end[i]
  GeneInfo$strand[i] <- gffdata_mRNA$strand[i]
  teststr <- unlist(strsplit(as.character(gffdata_gene_mod[i,1]), ";", fixed = T))
  for (j in 1:length(teststr)) {
    teststr2 <-unlist(strsplit(teststr[j], "=", fixed = T))
    GeneInfo[i, teststr2[1]] <- teststr2[2]
  }
  if (is.na(GeneInfo$locus_tag[i])==T) {
    GeneInfo$locus_tag[i] <- gsub(".t00", "", GeneInfo$Parent[i]) #this if statement is added in to make the locus_tag column have all the PSYTB #s
  }
}

write.csv(GeneInfo, file="Ptab_gff_info.csv")

DE_data <- read.csv("DEseq_PtabTime_out.csv", header=T)
gene_info <- read.csv("Ptab_gff_info.csv", header=T)
gene_info <- gene_info[,c("locus_tag","start","end","strand","product", "EC_number")]
DeSeq_Info_Ptab <- merge.data.frame(gene_info, DE_data, by.x = "locus_tag", by.y="X", all=T)
DeSeq_Info_Ptab <- DeSeq_Info_Ptab[3:dim(DeSeq_Info_Ptab)[1],] #getting rid of two not PSYTB locus tag names
write.csv(DeSeq_Info_Ptab, file="Annotated_Out/Ptab_DESeq_Time_Annotated.csv")

# Determining orthologs between Ptab and B728a ----------------------------
#Data initially came from the output of OMA, run on the HPC.
OMA <- read.table("prepped/Psyringae_B728a_proteins-Ptab_protein_all.txt", sep='\t', quote="", fill=T)
gene_num <- dim(OMA)[1]
col_head <- c("locus_tag_B728a", "locus_tag_Ptab")
GeneInfo <- as.data.frame(matrix(ncol=length(col_head), nrow=gene_num, dimnames=list(NULL,col_head)))
for(i in 1:gene_num){
  teststr <- unlist(strsplit(as.character(OMA[i,3]), " ", fixed = T))
  lt1 <- grep("locus_tag", teststr, value = T)
  lt2 <- gsub("\\[locus_tag=", "", lt1)
  GeneInfo$locus_tag_B728a[i] <- gsub("\\]", "", lt2)
  teststr <- unlist(strsplit(as.character(OMA[i,4]), " ", fixed = T))
  lt3 <- grep("locus_tag", teststr, value = T)
  lt4 <- gsub("\\[locus_tag=", "", lt3)
  GeneInfo$locus_tag_Ptab[i] <- gsub("\\]", "", lt4)
}
write.csv(GeneInfo, file="Pseudomonas_ortho_loci.csv")

#Merging DE Data (renaming columns to be unique) and the ortholog names - had to change all PSYTB genes to the same prefix
DE_data1 <- read.csv("Annotated_Out/Ptab_DEseq_Time_Annotated.csv", header=T, row.names = 1)
colnames(DE_data1) <- c("locus.Ptab", "start.Ptab", "end.Ptab", "strand.Ptab", "product.Ptab", "EC_number.Ptab", "baseMean.Ptab", "L2FC.Ptab", "lfcSE.Ptab", "stat.Ptab", "pvalue.Ptab", "padj.Ptab")
DE_data2 <- read.csv("Annotated_Out/B728a_DEseq_Time_Annotated.csv", header=T, row.names = 1)
colnames(DE_data2) <- c("locus.B728a", "start.B728a", "end.B728a", "strand.B728a", "name.B728a", "baseMean.B728a", "L2FC.B728a", "lfcSE.B728a", "stat.B728a", "pvalue.B728a", "padj.B728a")
ortho_list <- read.csv("Pseudomonas_ortho_loci.csv", header=T, row.names = 1)
#DeSeq_Ortho_Ptab <- merge.data.frame(ortho_list, DE_data1, by.x = "locus_tag_Ptab", by.y="locus_tag")
DeSeq_Ortho_B728a <- merge.data.frame(ortho_list, DE_data2, by.x = "locus_tag_B728a", by.y="locus.B728a", all=T)
DeSeq_Ortho_B728a_Ptab <- merge.data.frame(DeSeq_Ortho_B728a, DE_data1, by.x = "locus_tag_Ptab", by.y="locus.Ptab", all=T)
write.csv(DeSeq_Ortho_B728a_Ptab, file="Annotated_Out/Pseudomonas_DESeq_Ortho_Annotated_All.csv")

#just checking to see how often there are duplicate orthologs
n_occur <- data.frame(table(ortho_list$locus_tag_Ptab))
n_occur_mult <- n_occur[n_occur$Freq > 1,]
sum(n_occur_mult[,2])
n_occur2 <- data.frame(table(ortho_list$locus_tag_B728a))
n_occur2_mult <- n_occur2[n_occur2$Freq > 1,]
sum(n_occur2_mult[,2])

# Creating a Venn Diagram -------------------------------------------------
#reading in data
Ptab_B728a <- read.csv("Annotated_Out/Pseudomonas_DESeq_Ortho_Annotated_All.csv", row.names = 1)

#filtering for DEGs that have orthologs
Ptab.degs <- filter(Ptab_B728a, padj.Ptab<0.05 & abs(L2FC.Ptab)>2 & is.na(locus_tag_B728a)==F)[,1] 
B728a.degs <- filter(Ptab_B728a, padj.B728a<0.05 & abs(L2FC.B728a)>2 & is.na(locus_tag_Ptab)==F)[,1]

#drawing the plot
venn.plot <- venn.diagram(x=list("Psy"= B728a.degs, "Pta"= Ptab.degs),
                          filename = "Plots/Ortho_Ptab_B728a_Venn.svg",
                          imagetype="svg",
                          col = c("#3B5284", "#94B447"),
                          lwd =0.5,
                          cex = 1,
                          fontfamily = "sans",
                          height = 10, 
                          width = 10, 
                          resolution = 300,
                          cat.cex = 1,
                          cat.fontface = "bold",
                          cat.default.pos = "outer",
                          cat.fontfamily = "sans",
                          fill = c("#5BA8A0", "#CBe54E"),
                          cat.col = c("#3B5284", "#94B447"),
                          #cat.dist = c(0.055, 0.055),
                          margin=0.1
            )

#filtering for DEGs that don't have orthologs to get a venn diagram with no overlap for those genes
Ptab.degs2 <- filter(Ptab_B728a, padj.Ptab<0.01 & abs(L2FC.Ptab)>2 & is.na(locus_tag_B728a)==T)[,1] 
B728a.degs2 <- filter(Ptab_B728a, padj.B728a<0.01 & abs(L2FC.B728a)>2 & is.na(locus_tag_Ptab)==T)[,2]

venn.plot <- venn.diagram(x=list("Psy"= B728a.degs2, "Pta"= Ptab.degs2),
                          filename = "Plots/NoOrtho_Ptab_B728a_Venn.svg",
                          imagetype="svg",
                          col = c("#3B5284", "#94B447"),
                          lwd =0.5,
                          cex = 1,
                          fontfamily = "sans",
                          height = 10, 
                          width = 10, 
                          resolution = 300,
                          cat.cex = 1,
                          cat.fontface = "bold",
                          cat.default.pos = "outer",
                          cat.fontfamily = "sans",
                          fill = c("#5BA8A0", "#CBe54E"),
                          cat.col = c("#3B5284", "#94B447"),
                          cat.dist = c(0.055, 0.055),
                          margin=0.1
)

# Figure Heat Maps --------------------------------------------------------
#have to create the data tables manually based on genes of interest and L2FC values and I had them in their own folder (ForBacterialFigures)

T3_data <- read.csv("ForBacterialFigures/T3v2.csv", header=T)

T3SS_data <- subset(T3_data,
                   Effector=="N")
  
T3E_data <- subset(T3_data,
                   Effector=="Y")

T3SS <- ggplot(T3SS_data, aes(x=Strain, y=Name, fill= L2FC))+ #sets data for plot, divide by strain and name, color based on L2FC
  geom_tile_pattern( #Make boxes for the plot
    color = "black", #Line color between boxes
    aes(pattern = Threshold), #Adds in stripes for Threshold
    show.legend = F) + #removes legend
  scale_pattern_manual(values=c('none', 'stripe'), na.value = "none") + #sets while Threshold value is stripes or not
  geom_text(aes(label = Locus), color = "black", size = 3.5)+ #Adds in locus names in black font over boxes
  scale_fill_viridis(option="magma", na.value = "white", limits = c(-6, 9))+ #sets color palette and scale of box fill
  theme_minimal()+ #removes grey
  scale_y_discrete(expand=c(0, 0))+ #removes distance between edge of graph and labels
  scale_x_discrete(expand=c(0, 0))+ #removes distance between edge of graph and labels
  #ylab("Gene") + #y axis label
  ggtitle("Type III Secretion System")+ #plot title
  theme(
    axis.text.x=element_text(size=12, face="bold", color = "gray16"), #size and angle of x axis values
    axis.text.y=element_text(angle=0, size=10), #size and angle of y axis values
    axis.title.x=element_blank(), #removes x axis title
    axis.title.y=element_blank(), #size of y axis title text
    title=element_text(size=12)) #size of title text

T3E <- ggplot(T3E_data, aes(x=Strain, y=Name, fill= L2FC))+
  geom_tile_pattern(color = "black", aes(pattern = Threshold), show.legend = F)+
  scale_pattern_manual(values=c('stripe', 'none'), na.value = "none")+ 
  geom_text(aes(label = Locus), color = "black", size = 3.5)+
  scale_fill_viridis(option="magma", na.value = "white", limits = c(-6, 9))+
  theme_minimal()+
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  #ylab("Gene") +
  ggtitle("Type III Effectors")+
  theme(
    axis.text.x=element_text(size=12, face="bold", color = "gray16"),
    axis.title.x=element_blank(),
    axis.text.y=element_text(angle=0, size=10),
    axis.title.y=element_blank(),
    title=element_text(size=12))


TabTox_data <- read.csv("ForBacterialFigures/TabTox.csv", header=T)

TabTox <- ggplot(TabTox_data, aes(x=Strain, y=Name, fill= L2FC))+
  geom_tile_pattern(color = "black", aes(pattern = Threshold), show.legend = F)+
  scale_pattern_manual(values=c('stripe', 'none'))+ 
  geom_text(aes(label = Locus), color="white", size = 3.5)+
  scale_fill_viridis(option="magma", na.value = "white", limits = c(-6, 9))+
  theme_minimal()+
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  ylab("Gene") +
  ggtitle("Tabtoxin and\nPhevamine")+
  theme(
    axis.text.x=element_text(size=12, face="bold", color = "gray16"),
    axis.title.x=element_blank(),
    axis.text.y=element_text(angle=0, size=10),
    axis.title.y=element_text(size=12, face="bold"),
    title=element_text(size=12))

SyrSyp_data <- read.csv("ForBacterialFigures/SyrSyp.csv", header=T)

SyrSyp <- ggplot(SyrSyp_data, aes(x=Strain, y=Name, fill= L2FC))+
  geom_tile_pattern(color = "black", aes(pattern = Threshold), show.legend = F)+
  scale_pattern_manual(values=c('stripe', 'none'))+ 
  geom_text(aes(label = Locus), color = "black", size = 3.5)+
  scale_fill_viridis(option="magma", na.value = "white", limits = c(-6, 9))+
  theme_minimal()+
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  ylab("Gene") +
  ggtitle("Syringomycin,\nSyringolin,\nand Syringopeptin")+
  theme(
    axis.text.x=element_text(size=12, face="bold", color = "gray16"),
    axis.title.x=element_blank(),
    axis.text.y=element_text(angle=0, size=10),
    axis.title.y=element_text(size=12, face="bold"),
    title=element_text(size=12))

Pil_data <- read.csv("ForBacterialFigures/Pilv2.csv", header=T)

Pil<- ggplot(Pil_data, aes(x=Strain, y=Name, fill= L2FC))+
  geom_tile_pattern(color = "black", aes(pattern = Threshold), show.legend = F)+
  scale_pattern_manual(values=c('stripe', 'none', "none"), na.value = "none")+
  geom_text(aes(label = Locus), color = "black", size = 3.5)+
  scale_fill_viridis(option="magma", na.value = "white", limits = c(-6, 9))+
  theme_minimal()+
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  #ylab("Gene") +
  ggtitle("Type IV Pili")+
  theme(
    axis.text.x=element_text(size=12, face="bold", color = "gray16"),
    axis.text.y=element_text(angle=0, size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    title=element_text(size=12),
    legend.key=element_blank())

Alg_data <- read.csv("ForBacterialFigures/alg.csv", header=T)

Alg<- ggplot(Alg_data, aes(x=Strain, y=Name, fill= L2FC))+
  geom_tile_pattern(color = "black", aes(pattern = Threshold), show.legend = F)+
  scale_pattern_manual(values=c('stripe', 'none', "none"), na.value = "none")+
  geom_text(aes(label = Locus), color = "black", size = 3.5)+
  scale_fill_viridis(option="magma", na.value = "white", limits = c(-6, 9))+
  theme_minimal()+
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  #ylab("Gene") +
  ggtitle("Alginate")+
  theme(
    axis.text.x=element_text(size=12, face="bold", color = "gray16"),
    axis.title.x=element_blank(),
    #axis.text.y=element_text(angle=0, size=10),
    axis.title.y=element_blank(),
    title=element_text(size=12))

Flg_data <- read.csv("ForBacterialFigures/Flagella2.csv", header=T)

Flg<- ggplot(Flg_data, aes(x=Strain, y=Name, fill= L2FC))+
  geom_tile_pattern(color = "black", aes(pattern = Threshold))+
  scale_pattern_manual(values=c('stripe', 'none'), na.value = "none")+ 
  geom_text(aes(label = Locus), color = "black", size = 3.5)+
  scale_fill_viridis(option="magma", na.value = "white", limits = c(-6, 9))+
  theme_minimal()+
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  #ylab("Gene") +
  ggtitle("Flagella")+
  theme(
    axis.text.x=element_text(size=12, face="bold"),
    axis.title.x=element_blank(),
    #axis.text.y=element_text(angle=0, size=10),
    axis.title.y=element_blank(),
    title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.text=element_text(size=10),
    legend.title=element_text(size=12, face="bold"),
    legend.title.align = 0)+
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 20),
         pattern = "none")
         
layout <- "
CAABBGG
CAABBGG
CAABBGG
CAABBGG
DEEBBGG
DEEFFGG
DEEFFGG
DEEFFGG
"
#using patchwork to put them together
final <-
  T3SS + T3E + TabTox + SyrSyp + Pil + Alg + Flg + plot_layout(design=layout, guides="collect")# + plot_annotation(tag_levels = 'A')
#exported it as 13 x 8 and modified from there in Illustrator

#Supp for tailocin data
Tail_data <- read.csv("ForBacterialFigures/Tailocins.csv", header=T)

Tailocin <- ggplot(Tail_data , aes(x=Strain, y=Name, fill= L2FC))+ #sets data for plot, divide by strain and name, color based on L2FC
  geom_tile_pattern( #Make boxes for the plot
    color = "black", #Line color between boxes
    aes(pattern = Threshold)) + #Adds in stripes for Threshold
  scale_pattern_manual(values=c('stripe', 'none'), na.value = "none") + #sets while Threshold value is stripes or not
  geom_text(aes(label = Locus), color = "black", size = 3.5)+ #Adds in locus names in black font over boxes
  scale_fill_viridis(option="magma", na.value = "white", limits = c(-6, 9))+ #sets color palette and scale of box fill
  theme_minimal()+ #removes grey
  scale_y_discrete(expand=c(0, 0))+ #removes distance between edge of graph and labels
  scale_x_discrete(expand=c(0, 0))+ #removes distance between edge of graph and labels
  ggtitle("Tailocin Operons")+ #plot title
  theme(
    axis.text.x=element_text(size=12, face="bold"),
    axis.title.x=element_blank(),
    #axis.text.y=element_text(angle=0, size=10),
    axis.title.y=element_blank(),
    title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.text=element_text(size=10),
    legend.title=element_text(size=12, face="bold"),
    legend.title.align = 0)+
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 20),
         pattern = "none")
