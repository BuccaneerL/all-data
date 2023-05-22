###############1.GTEX###########################
#setwd("C:\\Rlangrange\\Disulfidptosis\\new\\1.GTEx")
library(data.table) #载入数据用
rm(list=ls())

#表达矩阵
exp_gtex.=fread("gtex_gene_expected_count.gz",header = T, sep = '\t',data.table = F)
rownames(exp_gtex.)=exp_gtex.[,1]
exp_gtex.=exp_gtex.[,-1]

#样本信息
data_cl=fread("GTEX_phenotype.gz",header = T, sep = '\t',data.table = F)
data_cl=data_cl[,c(1,3)]
names(data_cl)=c('Barcode','Tissue')
data_cl=data_cl[data_cl$Tissue == 'Brain',] #筛选出Prostate的数据

#注释信息
annotat=fread("gencode.v23.annotation.gene.probemap",header = T, sep = '\t',data.table = F)
annotat=annotat[,c(1,2)]
rownames(annotat)=annotat[,1] #这里没有选择删去id这一列

#筛选，筛选之后还剩100个barcode
exp_gtex.=exp_gtex.[,colnames(exp_gtex.) %in% data_cl$Barcode]
#还原为count
exp_gtex.=2^exp_gtex.-1

#基因注释
exp_gtex.=as.matrix(exp_gtex.)
t_index=intersect(rownames(exp_gtex.),rownames(annotat)) #行名取交集，t_index中是能够进行注释的probe_id
exp_gtex.=exp_gtex.[t_index,]
annotat=annotat[t_index,]
rownames(exp_gtex.)=annotat$gene

#去除重复基因名
t_index1=order(rowMeans(exp_gtex.),decreasing = T)
t_data_order=exp_gtex.[t_index1,]
keep=!duplicated(rownames(t_data_order))#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
exp_gtex.=t_data_order[keep,]#得到最后处理之后的表达谱矩阵

#读出
write.csv(exp_gtex.,file = "exp_gtex.count.csv",quote = FALSE)


################2.tcga############################
##第一步.TCGA数据下载

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("coin",force = TRUE)

BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data",force = T)

BiocManager::install("BioinformaticsFMRP/TCGAbiolinks",force = T)

BiocManager::install("SummarizedExperiment",force = TRUE)

library(SummarizedExperiment)
library(TCGAbiolinks)



####使用TCGAbiolink这个包下载，三个步骤，查询，下载，整理
#查询
query <- GDCquery(project = "TCGA-LGG",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",   
                  workflow.type = "STAR - Counts", 
                  legacy = FALSE)

####获取TCGA中各种癌症的project
projects <- TCGAbiolinks::getGDCprojects()$project_id 
projects <- projects[grepl('^TCGA', projects, perl=TRUE)]



#下载，将文件下载到工作目录中，并在目录下创建GDCdata的文件夹
#下次使用时，不再需要运行该句代码，重复下载
GDCdownload(query = query)


#整理，将下载好的文件整理成SE格式
SE=GDCprepare(query = query)

#保存该数据，以后可以直接加载使用
save(SE,file = 'TCGA_LGG_SE.Rdata')

load('TCGA_LGG_SE.Rdata')

names(SE@assays)    ###assays包含的是打包的各种数据类型的基因表达量的信息
names(SE@rowRanges) ###rowRanges包含基因注释信息

###获取我们所需要的基因注释信息
mydf <- as.data.frame(rowRanges(SE))
###查看基因类型
table(mydf$gene_type)

###获取所需的基因表达量的数据类型：tpm_unstrand，fpkm_unstrand，unstranded（counts）
gene_expr<- as.data.frame(assay(SE,i='tpm_unstrand'))
# gene_expr<- as.data.frame(assay(SE,i='unstranded'))

###Ensemble ID转换
gene_expr <- cbind(type=mydf$gene_type,ID=mydf$gene_name,gene_expr)

library(limma)
test4<-as.matrix(gene_expr )
rownames(test4)<-test4[,2]
exp<-test4[,3:ncol(test4)]#去第一行
dimnames<-list(rownames(exp),colnames(exp))
data<-matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

t_index=order(rowMeans(data),decreasing = T)#计算所有行平均值，按降序排列
t_data_order=data[t_index,]#调整表达谱的基因顺序
keep=!duplicated(rownames(t_data_order))#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
exp_tcga=t_data_order[keep,]#得到最后处理之后的表达谱矩阵

write.csv(exp_tcga,file = "exp_tcga.LGG.count.csv",quote = FALSE) 

gene_expr<- as.data.frame(assay(SE,i='tpm_unstrand'))
gene_expr <- cbind(type=mydf$gene_type,ID=mydf$gene_name,gene_expr)
library(limma)
test4<-as.matrix(gene_expr )
rownames(test4)<-test4[,2]
exp<-test4[,3:ncol(test4)]#去第一行
dimnames<-list(rownames(exp),colnames(exp))
data<-matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
t_index=order(rowMeans(data),decreasing = T)#计算所有行平均值，按降序排列
t_data_order=data[t_index,]#调整表达谱的基因顺序
keep=!duplicated(rownames(t_data_order))#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
exp_tcga=t_data_order[keep,]#得到最后处理之后的表达谱矩阵
write.csv(exp_tcga,file = "exp_tcga.LGG.TPM.csv",quote = FALSE) 






###############3.数据合并#############
#setwd("C:\\Rlangrange\\Disulfidptosis\\new\\3.HE QUPICI")

exp_tcga<- read.csv(file = "exp_tcga.LGG.count.csv", header=T, row.names=1,check.names=FALSE)
Tumor <- substr(colnames(exp_tcga),14,16)=='01A'
Tumor <- exp_tcga[,Tumor]

Normal <- substr(colnames(exp_tcga),14,16)=='11A'
Normal <- exp_tcga[,Normal]

ncol(Tumor)  ##查看肿瘤样本
ncol(Normal ) ##查看正常样本
exp_tcga <- cbind(Normal,Tumor)#合并505

exp_gtex <- read.csv(file = "exp_gtex.count.csv", header=T, row.names=1,check.names=FALSE)

t_index=intersect(rownames(exp_gtex),rownames(exp_tcga))
exp_m = cbind(exp_tcga[t_index,],exp_gtex[t_index,])
rm(exp_tcga,exp_gtex)
save(exp_m,file="exp_m_Rdata")


################################提取出死亡相关的基因##########################################
#与死亡相关基因取交集

library(limma)

load("exp_m_Rdata")
data=exp_m#表达输入文件
geneFile="DRGs.txt"       #基因列表文件
gene=read.table(geneFile, header=F, check.names=F, sep="\t")
sameGene=intersect(gene[,1], rownames(data))
same<- as.data.frame(sameGene)
geneExp=data[sameGene,]
exprSetX<-data[sameGene,]
save(exprSetX,file= "exprsetx.Rdata")
########################4.limma差异分析 ###########################################  
# setwd("C:\\Rlangrange\\3")
library(limma)

rm(list=ls())
load( "exprsetx.Rdata")

# 需要进行log2数据转换
exprSet <- as.data.frame(exprSetX)
exprSet <- log2(exprSet+1)
exprSetX<-exprSet


#exprSetX<-combat_Expr
#表达矩阵(exprSet)
dat1 <- exprSetX[,c(1:505)]
dat2 <- exprSetX[,c(506:1653)]
dat <- cbind(dat2,dat1)
exprSet=dat

#分组矩阵(designtable)
group <- rep(c("N","T"),c(1148,505))
designtable<-model.matrix(~0+factor(group))#将分组的向量类型转化为因子
colnames(designtable)<-levels(factor(group))#改列名
designtable<-as.data.frame(designtable)
exprSet<-as.data.frame(exprSet)
rownames(designtable)<-colnames(exprSet)#改行名

#比较矩阵(contrast)
contrast.matrix=makeContrasts("T - N",levels=designtable)#将control和HHT进行比较，定义水平值即排序的优先值按什么比较（这里只有两个组所以没意义，多组有意义）

#limma DEG
fit<-lmFit(exprSet,designtable)#线性拟合模型构建【需要两个东西：exprSet和designtable】。得到的结果再和contrast一起导入
fit2<-contrasts.fit(fit,contrast.matrix)#引入比较矩阵
fit2<-eBayes(fit2)#贝斯检验：利用上一步contrasts.fit（）的结果
allDiff<-topTable(fit2,coef="T - N",n=Inf)#输出：根据比较的组修改coef值；INF输出前面的多少

#########挑选DEGs##
allDiff$group<-ifelse(allDiff$P.Value>0.05,"no_change",#建立一个GROUP,筛选P值大于0.05的就没有显著性差异（P值是根据前面的数据出来的，而校正后的P值是根据P值出来的算双重保险，但是删选出来的基因较少故一般用P值），
                      ifelse(allDiff$logFC>1,"up",#对于小于0.05的就用logFC进行比较，等于HHT减去CONTRAST后取了一个log2,相当于是看两组表达量相差了多少倍
                             ifelse(allDiff$logFC< -1,"down","no_change")))#logFC大于1表示上调，小于-1表示下调，之间表示没有差异
table(allDiff$group)
allDiffDEGs <- select(allDiff,c(logFC,P.Value )) %>% arrange(.,P.Value) 
write.csv(allDiffDEGs,file = "allDiffDEGs.csv")


dif<-allDiff[allDiff[,"P.Value"]<0.05&abs(allDiff[,"logFC"])>1,]#取出P值小于0.05而且绝对值大于1的便于后续做KEGG等
library(dplyr)
DEGs <- select(dif,c(logFC,P.Value )) %>% arrange(.,P.Value)  
head(DEGs) 

save(DEGs,group,exprSet, allDiff, file = "DEGs.Rdata")
write.csv(DEGs,file = "DEGs.csv")


###############################5.箱线图##################################

library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)

rm(list=ls())
load("DEGs.Rdata")
##基因表达
Exp<-t(exprSet)

gene <- c("FLNA","FLNB","MYH9","TLN1","ACTB","MYL6","MYH10","CAPZB","DSTN","IQGAP1","ACTN4","PDLIM1","CD2AP","INF2","SLC7A11")#这里我们只选择这几个基因做数据
gene <- as.vector(gene)
Exp_plot <- Exp[,gene]#提取需要作图得基因表达信息

##样本
Type <- rep(c("N","T"),c(1148,505))
info<-as.data.frame(Type)
info$Sample<-rownames(Exp_plot)
#加载样本信息
Exp_plot<-as.data.frame(Exp_plot)
Exp_plot<- Exp_plot[info$Sample,]
Exp_plot$sam=info$Type
Exp_plot$sam <- factor(Exp_plot$sam,levels=c("N","T"))
#设置颜色
col <-c("#5CB85C","#D9534F")
#循环作图
plist2<-list()
for (i in 1:length(gene)){
  bar_tmp<-Exp_plot[,c(gene[i],"sam")]
  colnames(bar_tmp)<-c("Expression","sam")
  my_comparisons1 <- list(c("N", "T")) 
  
  pb1<-ggboxplot(bar_tmp,
                 x="sam",
                 y="Expression",
                 color="sam",
                 fill=NULL,
                 add = "jitter",
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size=0.01,
                 font.label = list(size=30), 
                 palette = col)+theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb1<-pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 0.5))
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gene[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA")#
  pb1<-pb1+stat_compare_means(method="t.test",hide.ns = F,
                              comparisons =c(my_comparisons1),
                              label="p.signif")
  plist2[[i]]<-pb1 
}

plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],
          plist2[[4]],plist2[[5]],plist2[[6]],
          plist2[[7]],plist2[[8]],plist2[[9]],
          plist2[[10]],plist2[[11]],plist2[[12]],plist2[[13]],
          plist2[[14]],plist2[[15]],ncol=5)#ncol=5表示图片排为几列



# -------------------------------------------------------------------------

# 火山图
# 将文件改为allDiff,列名logFC，PValue
rm(list=ls())
load("DEGs.Rdata")
head(allDiff)
colnames(allDiff)[4] = "PValue" 

cutoff_logFC = 1
cutoff_PValue = 0.05

allDiff$change = as.factor(ifelse(allDiff$PValue < cutoff_PValue & abs(allDiff$logFC) > cutoff_logFC,
                                  ifelse(allDiff$logFC >cutoff_logFC,'UP','DOWN'),'NOT') )

head(allDiff)
table(allDiff$change)

# 将"change"列转换为字符向量
allDiff$change <- as.character(allDiff$change)

# 将"ACTN4"和"IQGAP1"对应的行的"change"值修改为"risk"
allDiff["ACTN4", "change"] <- "risk"
allDiff["IQGAP1", "change"] <- "risk"

head(allDiff)
table(allDiff$change)
library(ggplot2)

# 设置基因名标签的阈值，只给"risk"类别的点添加标签
label_threshold <- -log10(0.05)

ggplot(data = allDiff, aes(x = logFC, y = -log10(PValue), color = change)) +     
  geom_hline(yintercept = -log10(cutoff_PValue), lty = 4, lwd = 0.6, alpha = 0.8) + 
  geom_vline(xintercept = c(cutoff_logFC, -cutoff_logFC), lty = 4, lwd = 0.6, alpha = 0.8) + 
  geom_point() +
  geom_text(aes(label = ifelse(-log10(PValue) >= label_threshold & change == "risk", rownames(allDiff), "")),
            vjust = 0.5, hjust = 1.25, size = 5, color = "black") +
  scale_color_manual(values = c("grey", "green", "red")) +
  labs(x = "log2 (fold change)",
       y = "-log10 (P.Value)",
       title = "Volcano Of DEGs") +
  theme_bw(base_size = 25)





#差异基因相关性分析
library(corrplot)
load("DEGs.Rdata")
data=exprSet#表达输入文件
M = t(data)
matrix <- cor (M,method = "pearson")#计算样本间的相关系数
#直接画图，不设置其它参数
corrplot(matrix, 
         order = "hclust", 
         addrect = 1,
         tl.pos="lt",# tl.*参数依次设置变量标签的位置为左侧和右侧，为1，颜色为黑色，倾斜45度，距离热图边界0.5
         tl.cex=1.5, #字符大小
         tl.col="red",#颜色为黑色，，
         tl.srt =45,#倾斜45度
         tl.offset=0.25,#距离热图边界0.5
         cl.pos = "r",#cl.pos参数设置图例放置位置，可选参数：r(图右侧)，b(下侧)，n(不显示)。
         cl.length = 5,#cl.length设置图例的颜色数量，设置了col参数，则默认使用col设置的颜色数量。
         cl.cex = 1.5,# cl.cex参数设置图例数值字符大小
         cl.ratio = 0.4,#cl.ratio设置颜色图例的宽度
         cl.offset=0.01#cl.offset设置数值与颜色图例距离。
         )

write.csv(matrix,file = "matrix.csv")


