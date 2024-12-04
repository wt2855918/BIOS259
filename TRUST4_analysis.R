#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#[1]CESC fastq2trust Changed because of Github
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
.libPaths(c("/oak/stanford/groups/howchang/users/ydzhao/R_Library_4.2",.libPaths()[2]))

mydir  = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/Scripts/RNA/fastq2trust/"
dir.create(mydir,recursive=T)
setwd(mydir)

mydir = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/fastq/"
files = list.files(mydir)

#######################
#separate name group
#######################
sam = gsub("_1.fq","",files)
sam = gsub("_2.fq","",sam)
sam = unique(sam)

for(k in 1 : length(sam))
{
	cat("\r",k)
	
	tmpinf = paste0(mydir,sam[k])
	
	myoutf1 = paste("job_", sam[k], ".sh", sep="")
	conOut = file(myoutf1, "w")
	
	curLine = c("#!/bin/bash")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = paste0("#SBATCH --job-name=",myoutf1)
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --nodes=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --mem-per-cpu=64GB")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --ntasks-per-node=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --time=48:00:00")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --output=%x.%j.out")
	writeLines(curLine, conOut)
	curLine = c("#SBATCH --error=%x.%j.err")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --partition=howchang,owners,normal,sfgf")
	writeLines(curLine, conOut)
	writeLines("", conOut)
		
	curLine = c("source activate /oak/stanford/groups/howchang/users/ydzhao/Conda_Library/NGS")
	writeLines(curLine, conOut)
	
	curLine = c("module load biology")
	writeLines(curLine, conOut)
	
	curLine = c("module load samtools")
	writeLines(curLine, conOut)
	
	input_1 = paste0(mydir,sam[k],"_1.fq")
	input_2 = paste0(mydir,sam[k],"_2.fq")
	
	output = paste0("/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/trust/",sam[k],"/")
	dir.create(output,recursive=T)

	curLine = paste0("run-trust4 --abnormalUnmapFlag -1 ",input_1," -2 ",input_2," -f /oak/stanford/groups/howchang/users/ydzhao/Softwares/TRUST4/hg38_bcrtcr.fa --ref /oak/stanford/groups/howchang/users/ydzhao/Softwares/TRUST4/human_IMGT+C.fa -o ",output)
	writeLines(curLine, conOut)	
	close(conOut)

}

sbatch job_Sample_SQ24058578-WX556-1-WX556-1.sh

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#[2]OV fastq2bam
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
.libPaths(c("/oak/stanford/groups/howchang/users/ydzhao/R_Library_4.2",.libPaths()[2]))

mydir  = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/Scripts/RNA/fastq2bam/"
dir.create(mydir,recursive = T)
setwd(mydir)

mydir = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/fastq/"
files = list.files(mydir)
files = gsub("_1.fq","",files)
files = gsub("_2.fq","",files)
files = unique(files)
sam = files

for(k in 1 : length(sam))
{
	cat("\r",k)

	myoutf1 = paste("job_", sam[k], ".sh", sep="")
	conOut = file(myoutf1, "w")
	
	curLine = c("#!/bin/bash")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = paste0("#SBATCH --job-name=",myoutf1)
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --nodes=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --mem-per-cpu=64GB")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --ntasks-per-node=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --time=30:00:00")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --output=%x.%j.out")
	writeLines(curLine, conOut)
	curLine = c("#SBATCH --error=%x.%j.err")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --partition=howchang,owners,normal,sfgf")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("module load biology")
	writeLines(curLine, conOut)

	curLine = c("module load star/2.5.4b")
	writeLines(curLine, conOut)
	
	curLine = c("cd /oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/fastq/")
	writeLines(curLine, conOut)
	
	genome_file = "/oak/stanford/groups/howchang/users/ydzhao/Resources/Genome/Hg38/star_genome_d1_vd1_gtfv22/"
	myoutf = paste0("/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/STAR/",sam[k],"_")
	
	input_file_1 = paste0(sam[k],"_1.fq")
	input_file_2 = paste0(sam[k],"_2.fq")
	
	curLine = paste0("STAR --genomeDir ",genome_file," --readFilesIn ",input_file_1," ",input_file_2," --outFileNamePrefix ",myoutf," --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMattributes Standard")
		
	writeLines(curLine, conOut)	
	close(conOut)
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Summary the library size
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
.libPaths(c("/oak/stanford/groups/howchang/users/ydzhao/R_Library_4.2",.libPaths()[2]))

mydir  = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/Scripts/RNA/bam2lib/"
dir.create(mydir)
setwd(mydir)

mydir = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/STAR/"
files = list.files(mydir)
tag = grep(".bam",files)
files = files[tag]
#tag = grep("gz",files)
#files = files[-tag]
sam = gsub("_Aligned.sortedByCoord.out.bam","",files)

for(k in 1 : length(sam))
{
	cat("\r",k)
	
	myoutf1 = paste("job_", sam[k], ".sh", sep="")
	conOut = file(myoutf1, "w")
	
	curLine = c("#!/bin/bash")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = paste0("#SBATCH --job-name=",myoutf1)
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --nodes=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --mem-per-cpu=64GB")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --ntasks-per-node=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --time=10:00:00")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --output=%x.%j.out")
	writeLines(curLine, conOut)
	curLine = c("#SBATCH --error=%x.%j.err")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --partition=howchang,owners,normal,sfgf")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	input = paste0(mydir,files[k])
	output = paste0("/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/Lib/",sam[k],"_BAM_summary.txt")
	
	curLine = c("module load biology")
	writeLines(curLine, conOut)	
	
	curLine = c("module load samtools/1.8")
	writeLines(curLine, conOut)	
	
	curLine = paste0("samtools flagstat ",input," > ",output)
	writeLines(curLine, conOut)	
	close(conOut)

}

#samtools flagstat your.bam > your_bam_summary.txt

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#CESC TCR/BCR calling
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
.libPaths(c("/oak/stanford/groups/howchang/users/ydzhao/R_Library_4.2",.libPaths()[2]))
library(dplyr)
library(ineq)

#source("/oak/stanford/groups/howchang/users/ydzhao/Softwares/trust4_metric_functions.R")
mydir = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/trust/"
folders = list.files(mydir)

data = matrix(0,length(folders),2)
row.names(data) = folders
colnames(data) = c("BCR","TCR")
data = as.data.frame(data)

for(i in 1 : length(folders))
{
	cat("\r",i)
	
	tmpinf = paste0(mydir,folders[i],"/_report.tsv")
	res = read.csv(tmpinf,sep="\t",quote=NULL)
	colnames(res) = c("count","frequency","CDR3nt","CDR3aa","V","D","J","C","cid","cid_full_length")

	tag = grep("partial",res$CDR3aa)
	if(sum(tag)>0)
	{
		res = res[-tag,]
	}
	
	
	tag = grep("out_of_frame",res$CDR3aa)
	if(sum(tag)>0)
	{
		res = res[-tag,]
	}
	
	tag1 = grep("IG",res$V)
	tag2 = grep("TR",res$V)
	
	res1 = res[tag1,]
	res2 = res[tag2,]
	
	BCR = ineq(res1[,1],type="Gini")
	TCR = ineq(res2[,1],type="Gini")
	
	data[i,1] = BCR
	data[i,2] = TCR
	
}

myoutf = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/CESC_TRUST_Gini_index.txt"
write.table(data,myoutf,sep="\t",quote=F)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Calculate the CPM for TCR/BCR
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
.libPaths(c("/oak/stanford/groups/howchang/users/ydzhao/R_Library_4.2",.libPaths()[2]))
library(dplyr)
library(ineq)

#source("/oak/stanford/groups/howchang/users/ydzhao/Softwares/trust4_metric_functions.R")
mydir = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/trust/"
folders = list.files(mydir)

data = matrix(0,length(folders),2)
row.names(data) = folders
colnames(data) = c("BCR","TCR")
data = as.data.frame(data)

for(i in 1 : length(folders))
{
	cat("\r",i)
	
	tmpinf = paste0(mydir,folders[i],"/_report.tsv")
	res = read.csv(tmpinf,sep="\t",quote=NULL)
	colnames(res) = c("count","frequency","CDR3nt","CDR3aa","V","D","J","C","cid","cid_full_length")

	tag = grep("partial",res$CDR3aa)
	if(sum(tag)>0)
	{
		res = res[-tag,]
	}
	
	
	tag = grep("out_of_frame",res$CDR3aa)
	if(sum(tag)>0)
	{
		res = res[-tag,]
	}
	
	tag1 = grep("IG",res$V)
	tag2 = grep("TR",res$V)
	
	res1 = res[tag1,]
	res2 = res[tag2,]
	
	myinf = paste0("/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/CESC_WGS/data/RNA/Lib/",folders[i],"_BAM_summary.txt")
	lib_size = read.table(myinf,sep="\t",quote=NULL)
	value = lib_size[1,1]
	value = strsplit(value," ")[[1]][1]
	value = as.numeric(value)
	
	
	BCR = sum(res1[,1])/value
	TCR = sum(res2[,1])/value
	
	BCR = log2(BCR*1000000+1)
	TCR = log2(TCR*1000000+1)
	
	data[i,1] = BCR
	data[i,2] = TCR
	
}

myoutf = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/OV_ATAC_RNA_WGS/Data/Data/RNA/OV_TRUST_CPM.txt"
write.table(data,myoutf,sep="\t",quote=F)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Calculate the counts for TCR/BCR
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
.libPaths(c("/oak/stanford/groups/howchang/users/ydzhao/R_Library",.libPaths()[2]))
library(dplyr)
library(ineq)

#source("/oak/stanford/groups/howchang/users/ydzhao/Softwares/trust4_metric_functions.R")
mydir = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/OV_ATAC_RNA_WGS/Data/Data/RNA/trust/"
folders = list.files(mydir)

data = matrix(0,length(folders),2)
row.names(data) = folders
colnames(data) = c("BCR","TCR")
data = as.data.frame(data)

for(i in 1 : length(folders))
{
	cat("\r",i)
	
	tmpinf = paste0(mydir,folders[i],"/_report.tsv")
	res = read.csv(tmpinf,sep="\t",quote=NULL)
	colnames(res) = c("count","frequency","CDR3nt","CDR3aa","V","D","J","C","cid","cid_full_length")

	tag = grep("partial",res$CDR3aa)
	if(sum(tag)>0)
	{
		res = res[-tag,]
	}
	
	
	tag = grep("out_of_frame",res$CDR3aa)
	if(sum(tag)>0)
	{
		res = res[-tag,]
	}
	
	tag1 = grep("IG",res$V)
	tag2 = grep("TR",res$V)
	
	res1 = res[tag1,]
	res2 = res[tag2,]
	
	#BCR = sum(res1[,1])/value
	#TCR = sum(res2[,1])/value
	
	#BCR = log2(BCR*1000000+1)
	#TCR = log2(TCR*1000000+1)
	
	BCR = sum(res1[,1])
	TCR = sum(res2[,1])
	
	data[i,1] = BCR
	data[i,2] = TCR
	
}

myoutf = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/OV_ATAC_RNA_WGS/Data/Data/RNA/OV_TRUST_count.txt"
write.table(data,myoutf,sep="\t",quote=F)




