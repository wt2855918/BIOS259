sdev -c 1 -m 128GB -p howchang,normal,sfgf -t 12:00:00

module load hdf5/1.12.0  
module load R/4.0.2

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Annotate loops
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
.libPaths(c("/oak/stanford/groups/howchang/users/ydzhao/R_Library",.libPaths()[2]))
library(Gviz)
library(GenomicInteractions)
library(GenomicRanges)
library(InteractionSet)
library(rtracklayer)
library(biomaRt)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#FYI to check the avaible TXdb and the tables
#ucscGenomes()[ , "db"]
#supportedUCSCtables()

#hg38.refseq.db <- makeTxDbFromUCSC(genome="hg38", table="refGene")

hg38.refseq.db <- TxDb.Hsapiens.UCSC.hg38.knownGene
refseq.genes = genes(hg38.refseq.db)
refseq.transcripts = transcriptsBy(hg38.refseq.db, by="gene")
non_pseudogene = names(refseq.transcripts) %in% unlist(refseq.genes$gene_id) 
refseq.transcripts = refseq.transcripts[non_pseudogene] 
refseq.promoters = promoters(refseq.transcripts, upstream=2500, downstream=2500)
refseq.promoters = unlist(refseq.promoters)
refseq.promoters$tx_name = sapply(refseq.promoters$tx_name,function(x) strsplit(x,".",fix=T)[[1]][1])

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id"), filter = "ensembl_transcript_id",
               values = refseq.promoters$tx_name, mart = mart)
               
refseq.promoters$geneSymbol <- genes$hgnc_symbol[match(refseq.promoters$tx_name, genes$ensembl_transcript_id)]

names(refseq.promoters) <- refseq.promoters$geneSymbol
na.symbol <- is.na(names(refseq.promoters))
names(refseq.promoters)[na.symbol] <- refseq.promoters$tx_name[na.symbol]

# unlist object so "strand" is one vector
refseq.transcripts.ul = unlist(refseq.transcripts) 
# terminators can be called as promoters with the strand reversed
strand(refseq.transcripts.ul) = ifelse(strand(refseq.transcripts.ul) == "+", "-", "+") 
refseq.terminators.ul = promoters(refseq.transcripts.ul, upstream=1000, downstream=1000) 
# change back to original strand
strand(refseq.terminators.ul) = ifelse(strand(refseq.terminators.ul) == "+", "-", "+") 
# `relist' maintains the original names and structure of the list
refseq.terminators = relist(refseq.terminators.ul, refseq.transcripts)

annotation.features = list(promoter=refseq.promoters, gene.body=refseq.transcripts)

#####@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#define the function
#####@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#define the function
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}

#####@@@@@@@@@@@@@@@@@@@@
#Begin to annotate loops
####@@@@@@@@@@@@@@@@@@@@@
tmpinf = "/oak/stanford/groups/howchang/users/ydzhao/Resources/Pub_Dat/TCGA/Merged.FitHiChIP.interactions_Q0.1_MergeNearContacts_normcounts.txt"
data = read.table(tmpinf,sep="\t",quote=NULL,header=T)
	
anchor_1 = data[,c(1,2,3)]
anchor_2 = data[,c(4,5,6)]
	
anchor_1 = peakDF2GRanges(anchor_1)
anchor_2 = peakDF2GRanges(anchor_2)
	
data = GenomicInteractions(anchor_1, anchor_2)

annotateInteractions(data, annotation.features)
annotationFeatures(data)
		
res = data
tmp1 = anchorOne(res)
tmp2 = anchorTwo(res)

tmp1 = as.data.frame(tmp1)
tmp2 = as.data.frame(tmp2)

colnames(tmp1) = paste0("anchor_1_",colnames(tmp1))
colnames(tmp2) = paste0("anchor_2_",colnames(tmp2))

tmp = cbind(tmp1,tmp2)
tmp = as.data.frame(tmp)
row.names(tmp) = paste0(tmp$anchor_1_seqnames,"_",tmp$anchor_1_start,"_",tmp$anchor_1_end,"_",tmp$anchor_2_seqnames,"_",tmp$anchor_2_start,"_",tmp$anchor_2_end)
res =tmp
myoutf = "/oak/stanford/groups/howchang/users/ydzhao/Resources/Pub_Dat/TCGA/Merge_loops_annotated_knownGene.Rda"
save(res,file=myoutf)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#Loop anchor annotation
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
rm(list=ls())
.libPaths(c("/oak/stanford/groups/howchang/users/ydzhao/R_Library",.libPaths()[2]))

library(GenomicInteractions)
library(GenomicRanges)
library(InteractionSet)
library(ChIPpeakAnno)
library(ArchR)
###################################################
#define the convert function
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}

tmpinf = "/oak/stanford/groups/howchang/users/ydzhao/Resources/Pub_Dat/TCGA/Merged.FitHiChIP.interactions_Q0.1_MergeNearContacts_normcounts.txt"
data = read.table(tmpinf,sep="\t",quote=NULL,header=T)
row.names(data) = paste0(data$chr1,"_",data$s1,"_",data$e1,"_",data$chr2,"_",data$s2,"_",data$e2)

anchor_1 = data[,c("chr1","s1","e1")]
anchor_2 = data[,c("chr2","s2","e2")]

anchor_1 = peakDF2GRanges(anchor_1)
anchor_2 = peakDF2GRanges(anchor_2)

res = matrix(0,nrow(data),6)
row.names(res) = row.names(data)
colnames(res) = c("C1_op","C2_op","C3_op","C4_op","C5_op","C6_op")
res = as.data.frame(res)

tmpinf1 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_immune_final/PeakCalls/C1-reproduciblePeaks.gr.rds"
tmpinf2 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_immune_final/PeakCalls/C2-reproduciblePeaks.gr.rds"
tmpinf3 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_immune_final/PeakCalls/C3-reproduciblePeaks.gr.rds"
tmpinf4 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_immune_final/PeakCalls/C4-reproduciblePeaks.gr.rds"
tmpinf5 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_immune_final/PeakCalls/C5-reproduciblePeaks.gr.rds"
tmpinf6 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_immune_final/PeakCalls/C6-reproduciblePeaks.gr.rds"

#@@@@@@@@@@@@@@@@@@@
#Plasma cell
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf1)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"C1_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#B cell
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf2)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"C2_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#T/NK cell
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf3)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"C3_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#Stromal cell
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf4)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"C4_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#Macrophage cell
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf5)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"C5_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#Microgrila
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf6)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"C6_op"] = 1
Value = apply(res,1,sum)
res[,"Non_Tumor"] = Value

Value = apply(res[,c("C1_op","C2_op","C3_op","C5_op")],1,sum)
res[,"Immune"] = Value
myoutf = "/oak/stanford/groups/howchang/users/ydzhao/Resources/Pub_Dat/TCGA/Union_HiChIP_loop_non_tumor_annotation.txt"
write.table(res,myoutf,sep="\t",quote=F)


tmpinf1 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_cancer/PeakCalls/BLCA-reproduciblePeaks.gr.rds"
tmpinf2 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_cancer/PeakCalls/BRCA-reproduciblePeaks.gr.rds"
tmpinf3 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_cancer/PeakCalls/COAD-reproduciblePeaks.gr.rds"
tmpinf4 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_cancer/PeakCalls/GBMx-reproduciblePeaks.gr.rds"
tmpinf5 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_cancer/PeakCalls/KIRC-reproduciblePeaks.gr.rds"
tmpinf6 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_cancer/PeakCalls/KIRP-reproduciblePeaks.gr.rds"
tmpinf7 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_cancer/PeakCalls/LUAD-reproduciblePeaks.gr.rds"
tmpinf8 = "/scratch/users/ydzhao/scATAC/Data/Benchmark/TCGA/pan_cancer/PeakCalls/SKCM-reproduciblePeaks.gr.rds"

res = matrix(0,nrow(data),8)
row.names(res) = row.names(data)
colnames(res) = c("BLCA_op","BRCA_op","COAD_op","GBM_op","KIRC_op","KIRP_op","LUAD_op","SKCM_op")
res = as.data.frame(res)

#@@@@@@@@@@@@@@@@@@@
#BLCA
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf1)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"BLCA_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#BRCA
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf2)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"BRCA_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#COAD
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf3)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"COAD_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#GBM
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf4)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"GBM_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#KIRC
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf5)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"KIRC_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#KIRP
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf6)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"KIRP_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#LUAD
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf7)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"LUAD_op"] = 1

#@@@@@@@@@@@@@@@@@@@
#SKCM
#@@@@@@@@@@@@@@@@@@@
info = readRDS(tmpinf8)
op1 = findOverlaps(anchor_1,info)
op2 = findOverlaps(anchor_2,info)

op1 = as.data.frame(op1)
op2 = as.data.frame(op2)

tag1 = unique(op1$queryHits)
tag2 = unique(op2$queryHits)
tag = intersect(tag1,tag2)
res[tag,"SKCM_op"] = 1

Value = apply(res,1,sum)
res[,"Tumor"] = Value

myoutf = "/oak/stanford/groups/howchang/users/ydzhao/Resources/Pub_Dat/TCGA/Union_HiChIP_loop_tumor_annotation.txt"
write.table(res,myoutf,sep="\t",quote=F)

####################################################################
myinf1 = "/oak/stanford/groups/howchang/users/ydzhao/Resources/Pub_Dat/TCGA/Union_HiChIP_loop_non_tumor_annotation.txt"
myinf2 = "/oak/stanford/groups/howchang/users/ydzhao/Resources/Pub_Dat/TCGA/Union_HiChIP_loop_tumor_annotation.txt"
myinf3 = "/scratch/users/ydzhao/scATAC/Data/Loop_infiltration_correlation.txt"
myinf4 = "/scratch/users/ydzhao/scATAC/Data/Loop_immune_correlation_Spearman.txt"

res1 = read.table(myinf1,sep="\t",quote=NULL)
res2 = read.table(myinf2,sep="\t",quote=NULL)
res3 = read.table(myinf3,sep="\t",quote=NULL)
res4 = read.table(myinf4,sep="\t",quote=NULL)

tag1 = res1$Immune > 0
tag2 = res2$Tumor > 0

xx1 = row.names(res1)[tag1]
xx2 = row.names(res2)[tag2]

com = intersect(xx1,xx2)
tag1 = which(xx1 %in% com)
tag2 = which(xx2 %in% com)

xx1 = xx1[-tag1]
xx2 = xx2[-tag2]

yy1 = res3[xx1,"Leukocyte.Fraction"]
yy2 = res3[xx2, "Leukocyte.Fraction"]

tag1 = res3[,"Leukocyte.Fraction"] > 0.4
tag2 = res4[,"Immune"] > 0.4

xx1 = row.names(res3)[tag1]
xx2 = row.names(res4)[tag2]

tag = res2$Tumor > 0
tumor = row.names(res2)[tag]
tag = res4$C7 > 0.3
tumor_loop = row.names(res4)[tag]

tag = res1$Immune > 0
non_tumor = row.names(res1)[tag]
tag = res1$Immune > 0.3
non_tumor_loop = row.names(res1)[tag]
non_tumor_loop = intersect(non_tumor_loop, non_tumor)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Plasma loops
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tag = res1$C1_op == 1 & res1$C5_op == 0 & res1$C2_op == 0 & res1$C3_op == 0 & res1$C4_op == 0 
loop = row.names(res1)[tag]
com = intersect(loop,tumor)
tag = which(loop %in% com)
loop = loop[-tag]

Immune_correlation = res3[loop,]
ATAC_correlation = res4[loop,]
tag = ATAC_correlation[,"C1"] >= 0.3
target = row.names(ATAC_correlation)[tag]
shared = row.names(ATAC_correlation)[tag]
shared = length(shared)
take = length(loop)

white = sum(res4[,"C1"] >= 0.3,na.rm=T)
black = nrow(res4) - white
ER = (shared/take)/(white/(white+black))
myp = sum(dhyper(shared:take, white, black, take))

> ER
[1] 1.490828
> myp
[1] 6.161813e-09

Plasma_correlation = res3[loop,"Leukocyte.Fraction"]
Tumor_correlation = res3[tumor_loop,"Leukocyte.Fraction"]

t.test(Plasma_correlation,Tumor_correlation)

	Welch Two Sample t-test

data:  Plasma_correlation and Tumor_correlation
t = 21.205, df = 3135.4, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.05755216 0.06927940
sample estimates:
  mean of x   mean of y 
 0.03550687 -0.02790891 

tag = res3[target,"Leukocyte.Fraction" ] > 0
target = target[tag]
#159 out of 219 out of 2981

myoutf = "/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C1_General_loop.Rda"
C1_loop_final_pan = target
save(C1_loop_final_pan,file=myoutf)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#B loops
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tag = res1$C2_op == 1 & res1$C1_op == 0 & res1$C5_op == 0 & res1$C3_op == 0 & res1$C4_op == 0 
loop = row.names(res1)[tag]
com = intersect(loop,tumor)
tag = which(loop %in% com)
loop = loop[-tag]

Immune_correlation = res3[loop,]
ATAC_correlation = res4[loop,]
tag = ATAC_correlation[,"C2"] >= 0.3
target = row.names(ATAC_correlation)[tag]
shared = row.names(ATAC_correlation)[tag]
shared = length(shared)
take = length(loop)

white = sum(res4[,"C2"] >= 0.3,na.rm=T)
black = nrow(res4) - white
ER = (shared/take)/(white/(white+black))
myp = sum(dhyper(shared:take, white, black, take))

B_correlation = res3[loop,"Leukocyte.Fraction"]
Tumor_correlation = res3[tumor_loop,"Leukocyte.Fraction"]

t.test(B_correlation, Tumor_correlation)

	Welch Two Sample t-test

data:  B_correlation and Tumor_correlation
t = 17.937, df = 2499.2, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.05189486 0.06463388
sample estimates:
  mean of x   mean of y 
 0.03035546 -0.02790891 

tmp = res3[target,"Leukocyte.Fraction" ]

187 out of 243 out of 2395

tag = res3[target,"Leukocyte.Fraction" ] > 0
target = target[tag]
myoutf = "/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C2_General_loop.Rda"
C2_loop_final_pan = target
save(C2_loop_final_pan,file=myoutf)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#T/NK loops
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tag = res1$C3_op == 1 & res1$C1_op == 0 & res1$C5_op == 0 & res1$C2_op == 0 & res1$C4_op == 0 
loop = row.names(res1)[tag]
com = intersect(loop,tumor)
tag = which(loop %in% com)
loop = loop[-tag]

Immune_correlation = res3[loop,]
ATAC_correlation = res4[loop,]
tag = ATAC_correlation[,"C3"] >= 0.3
target = row.names(ATAC_correlation)[tag]
shared = row.names(ATAC_correlation)[tag]
shared = length(shared)
take = length(loop)

white = sum(res4[,"C3"] >= 0.3,na.rm=T)
black = nrow(res4) - white
ER = (shared/take)/(white/(white+black))
myp = sum(dhyper(shared:take, white, black, take))

T_NK_correlation = res3[loop,"Leukocyte.Fraction"]
Tumor_correlation = res3[tumor_loop,"Leukocyte.Fraction"]

t.test(T_NK_correlation,Tumor_correlation)

	Welch Two Sample t-test

data:  T_NK_correlation and Tumor_correlation
t = 12.758, df = 595.79, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.08059295 0.10991914
sample estimates:
  mean of x   mean of y 
 0.06734714 -0.02790891 

tmp = res3[target,"Leukocyte.Fraction"]
35 out of 40 out of 492


tag = res3[target,"Leukocyte.Fraction" ] > 0
target = target[tag]
myoutf = "/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C3_General_loop.Rda"
C3_loop_final_pan = target
save(C3_loop_final_pan,file=myoutf)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Macrophage loops
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tag = res1$C5_op == 1 & res1$C1_op == 0 & res1$C2_op == 0 & res1$C3_op == 0 & res1$C4_op == 0 
loop = row.names(res1)[tag]
com = intersect(loop,tumor)
tag = which(loop %in% com)
loop = loop[-tag]

Immune_correlation = res3[loop,]
ATAC_correlation = res4[loop,]
tag = ATAC_correlation[,"C5"] >= 0.3
target = row.names(ATAC_correlation)[tag]
shared = row.names(ATAC_correlation)[tag]
shared = length(shared)
take = length(loop)

white = sum(res4[,"C5"] >= 0.3,na.rm=T)
black = nrow(res4) - white
ER = (shared/take)/(white/(white+black))
myp = sum(dhyper(shared:take, white, black, take))

> ER
[1] 2.348261
> myp
[1] 2.526444e-60

Macrophage_correlation = res3[loop,"Leukocyte.Fraction"]
Tumor_correlation = res3[tumor_loop,"Leukocyte.Fraction"]

t.test(Macrophage_correlation,Tumor_correlation)

	Welch Two Sample t-test

data:  Macrophage_correlation and Tumor_correlation
t = 24.522, df = 5996.5, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.04868705 0.05714791
sample estimates:
  mean of x   mean of y 
 0.02500857 -0.02790891 

tmp = res3[target,"Leukocyte.Fraction"]

#331 out of 450 out of 5434

tag = res3[target,"Leukocyte.Fraction" ] > 0
target = target[tag]
myoutf = "/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C5_General_loop.Rda"
C5_loop_final_pan = target
save(C5_loop_final_pan,file=myoutf)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Stromal loops
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tag = res1$C4_op == 1 & res1$C1_op == 0 & res1$C5_op == 0 & res1$C2_op == 0 & res1$C3_op == 0 
loop = row.names(res1)[tag]
com = intersect(loop,tumor)
tag = which(loop %in% com)
loop = loop[-tag]

Immune_correlation = res3[loop,]
ATAC_correlation = res4[loop,]
tag = ATAC_correlation[,"C4"] >= 0.3
target = row.names(ATAC_correlation)[tag]
shared = row.names(ATAC_correlation)[tag]
shared = length(shared)
take = length(loop)

white = sum(res4[,"C4"] >= 0.3,na.rm=T)
black = nrow(res4) - white
ER = (shared/take)/(white/(white+black))
myp = sum(dhyper(shared:take, white, black, take))

> ER
[1] 1.271948
> myp
[1] 9.924529e-05
> 

tag = res3[target,"Leukocyte.Fraction" ] > 0
target = target[tag]
myoutf = "/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C4_General_loop.Rda"
C4_loop_final_pan = target
save(C4_loop_final_pan,file=myoutf)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Tumor loops
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tag = res2$Tumor > 0
loop = row.names(res2)[tag]
com = intersect(loop,non_tumor)
tag = which(loop %in% com)
loop = loop[-tag]

Immune_correlation = res3[loop,]
ATAC_correlation = res4[loop,]
tag = ATAC_correlation[,"C7"] >= 0.3
target = row.names(ATAC_correlation)[tag]
shared = row.names(ATAC_correlation)[tag]
shared = length(shared)
take = length(loop)

white = sum(res4[,"C7"] >= 0.3,na.rm=T)
black = nrow(res4) - white
ER = (shared/take)/(white/(white+black))
myp = sum(dhyper(shared:take, white, black, take))

tmpinf = "/scratch/users/ydzhao/scATAC/Data/Loop_purity_correlation.txt"
res5 = read.table(tmpinf,sep="\t",quote=NULL)

Tumor_correlation = res5[loop,"IHC"]
non_tumor_correlation = res5[non_tumor_loop,"IHC"]

tmp = res3[target,"Leukocyte.Fraction"]
#10216 out of 18143 out of 140796


tag = res3[target,"Leukocyte.Fraction" ] < 0
target = target[tag]
myoutf = "/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C7_General_loop.Rda"
C7_loop_final_pan = target
save(C7_loop_final_pan,file=myoutf)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Check the Genes
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
load("/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C1_General_loop.Rda")
load("/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C2_General_loop.Rda")
load("/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C3_General_loop.Rda")
load("/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C4_General_loop.Rda")
load("/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C5_General_loop.Rda")
load("/scratch/users/ydzhao/scATAC/Data/Specific_loop_General_Final_Pan/C7_General_loop.Rda")

load("/oak/stanford/groups/howchang/users/ydzhao/Resources/Pub_Dat/TCGA/Merge_loops_annotated_knownGene.Rda")
res1 = res[C1_loop_final_pan,]
res2 = res[C2_loop_final_pan,]
res3 = res[C3_loop_final_pan,]
res4 = res[C4_loop_final_pan,]
res5 = res[C5_loop_final_pan,]
res7 = res[C7_loop_final_pan,]

C1_Gene = unique(as.vector(unlist(c(res1$anchor_1_promoter.id, res1$anchor_2_promoter.id))))
C2_Gene = unique(as.vector(unlist(c(res2$anchor_1_promoter.id, res2$anchor_2_promoter.id))))
C3_Gene = unique(as.vector(unlist(c(res3$anchor_1_promoter.id, res3$anchor_2_promoter.id))))
C4_Gene = unique(as.vector(unlist(c(res4$anchor_1_promoter.id, res4$anchor_2_promoter.id))))
C5_Gene = unique(as.vector(unlist(c(res5$anchor_1_promoter.id, res5$anchor_2_promoter.id))))
C7_Gene = unique(as.vector(unlist(c(res7$anchor_1_promoter.id, res7$anchor_2_promoter.id))))




files
 [1] "scATAC_BRCA_14AD76EE_12F9_40B3_8DCD_4A256E02CF8D_X003_S02_B1_T1_immune_working"
 [2] "scATAC_BRCA_7C6A3AE4_E2EA_42B3_B3F1_81C19E6F2170_X005_S02_B1_T1_immune_working"
 [3] "scATAC_BRCA_8D1E6006_85CB_484A_8B5C_30766D90137B_X001_S01_B1_T1_immune_working"
 
 
 
 [4] "scATAC_BRCA_8D1E6006_85CB_484A_8B5C_30766D90137B_X003_S03_B1_T2_immune_working"
 [5] "scATAC_BRCA_94AF19F0_1F2A_41EC_8CB6_96C76227811F_X013_S01_B1_T1_immune_working"
 [6] "scATAC_BRCA_C147AAD5_A8F1_41D5_8709_21820BE50902_X008_S02_B1_T1_immune_working"
 [7] "scATAC_BRCA_C9C8D426_A3FD_4455_89A9_768BC01D66A9_X009_S02_B1_T1_immune_working"
 [8] "scATAC_BRCA_DD69EDE9_142D_46E2_AA06_58D07D3230FB_X014_S08_B1_T1_immune_working"
 [9] "scATAC_COAD_0914606C_2CA1_4287_B530_DB70EA93ED6C_X006_S03_B1_T1_immune_working"
[10] "scATAC_GBMx_09C0DCE7_D669_4D28_980D_BF71179116A4_X005_S04_B1_T1_immune_working"
[11] "scATAC_GBMx_6BEE2CB6_9AFD_42A6_9C26_9C4428FBABFA_X004_S04_B1_T1_immune_working"
[12] "scATAC_GBMx_ED12A6C9_D96E_49C4_B882_6382A5FB4538_X003_S05_B1_T1_immune_working"
[13] "scATAC_KIRC_7D6A394E_01EC_4C58_A010_1C63147A376A_X005_S05_B1_T1_immune_working"
[14] "scATAC_LUAD_9B9C5C6D_1755_41CD_9BC1_BFBC5D383DE0_X003_S06_B1_T1_immune_working"
[15] "scATAC_LUAD_A1ABC0B2_D7F8_45B5_B431_AE23B19B8DED_X009_S06_B1_T1_immune_working"
[16] "scATAC_LUAD_CA6F245B_30E0_48DE_AB15_23736E21D312_X006_S06_B1_T1_immune_working"
[17] "scATAC_SKCM_211D9CF4_3348_4DCD_8A01_6827435DDB3D_X004_S07_B1_T1_immune_working"
[18] "scATAC_SKCM_FDA487D2_5293_4315_9212_3836856CCFFB_X008_S06_B1_T1_immune_working"
