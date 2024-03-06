
########

gencode.pcg$gene_id=gsub("\\..*","",gencode.pcg$gene_id)
##################
mycolors=c("#7F3C8D","#11A579","#3969AC","#80BA5A","#F2B701","#E73F74","#E68310","#008695","#CF1C90","#F97B72","#4B4B8F","#A5AA99")

clin.color=mycolors[4:5]
#################################
load("raw_datas/TCGA/luad.tcga.exp.RData")
dim(luad.tcga.exp) ## 25554*574

luad.tcga.exp=luad.tcga.exp[rownames(luad.tcga.exp) %in% gencode.pcg$gene_name,]
dim(luad.tcga.exp) ## 19037 * 574
range(luad.tcga.exp)
## 01: 513 肿瘤样本
## 02: 2
## 11:59 癌旁正常样本
table(substr(colnames(luad.tcga.exp),14,15))
####################
tcga.t.exp=luad.tcga.exp[,substr(colnames(luad.tcga.exp),14,15)=="01"]
tcga.n.exp=luad.tcga.exp[,substr(colnames(luad.tcga.exp),14,15)=="11"]
dim(tcga.t.exp) ## 513
dim(tcga.n.exp) ## 59

################################ WGCNA analysis
tcga.exp=cbind(tcga.t.exp,tcga.n.exp)
tcga.pheno=rbind(data.frame(Sample=colnames(tcga.t.exp),tissue="Tumor"),
                  data.frame(Sample=colnames(tcga.n.exp),tissue="Normal"))
##### cluster gene
library(WGCNA)
mads=apply(tcga.exp,1,mad)
pcg.selected=names(sort(mads,decreasing = T))[1:(0.7*length(mads))]
length(pcg.selected) ## 13325
###############
edata.wgcna=tcga.exp[pcg.selected,]
edata.wgcna=t(edata.wgcna)
dim(edata.wgcna) ## 572*13325
range(edata.wgcna)
## 0.00000 16.87588
########
mg_wgcna_get_power=function(exp,RsquaredCut=0.85,height=0,net_type='unsigned',maxpowers=30,blockSize = 7000){
  library(WGCNA)
  logs=c()
  #print('=====')
  filter.exp=exp
  sampleTree = hclust(dist(filter.exp), method = "average")
  #pdf(file = paste0(MG_GENE_FOLDER,'/SampleCluster.pdf'),width = 12,height = 6)
  #plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  #dev.off()
  #height=140000
  #print(height)
  #print(min(sampleTree$height))
  if(height>min(sampleTree$height)){
    cte=cutree(tree = sampleTree,h = height)
    #plot(sampleTree)
    #abline(h=height)
    s.inds=match(names(which(cte==names(which.max(table(cte))))),row.names(filter.exp))
    #sum(sampleTree$height<height)
    logs=c(logs,paste0('remove samples:',(nrow(filter.exp)-length(s.inds)),' row count')) 
    filter.exp=filter.exp[s.inds,]
    #print('==0')
  }
  #print('===1')
  powers = c(c(1:10), seq(from = 12, to=maxpowers, by=2))
  #print('===2')
  #print(head(filter.exp))
  sft = pickSoftThreshold(filter.exp, powerVector=powers, 
                          RsquaredCut=RsquaredCut,
                          networkType=net_type, verbose=5,blockSize = blockSize)
  #print('===3')
  cutPower=sft$powerEstimate
  if(is.na(cutPower)){
    cutPower=0
  }
  print(paste0('power succed:',cutPower,',starting plot'))
  logs=c(logs,paste0('power succed:',cutPower,',starting plot')) 
  
  #library(customLayout)
  #lay1 <- lay_new(
  #  matrix(1:1, nc = 1),
  #  widths = c(5),
  #  heights = c(4))
  #lay2 <- lay_new(
  #  matrix(1:2, nc = 2),
  #  widths = c(2.5,2.5),
  #  heights = c(4))
  
  #cl = lay_bind_row(lay1, lay2,heights = c(0.8,1))
  #lay_set(cl)
  
  layout(matrix(c(1,1,2,3),2,2,byrow=T),widths=c(1,1),heights=c(0.8,1))
  
  mai=par('mai')
  
  #pdf(file = paste0(MG_GENE_FOLDER,'/powers.pdf'),width = 9,height = 8)
  mai1=mai
  mai1[1]=0
  par(mai=mai1)
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  if(height>min(sampleTree$height)){
    abline(h=height,col='red')
  }
  par(mai=mai)
  
  #lay_show(cl)
  #par(mfrow = c(1,2))
  
  cex1 = rep(0.9,length(sft$fitIndices[,1]))
  col1=rep('red',length(sft$fitIndices[,1]))
  if(cutPower>0){
    cex1[which(sft$fitIndices[,1]==cutPower)]=1.5
    col1[which(sft$fitIndices[,1]==cutPower)]='blue'
  }
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col=col1)
  # abline(h=0.85,col="red")
  abline(h=RsquaredCut ,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
       cex=cex1, col=col1)
  arg=paste0('height=',height,',net_type=',net_type,',blockSize = ',blockSize)
  return(list(cutPower=cutPower,log=logs,Exp=filter.exp,arg=arg))  
}

pdf('analysis/01_WGCNA/Fig1ABC.pdf',width = 8,height = 8)
power=mg_wgcna_get_power(edata.wgcna,RsquaredCut=0.85)
dev.off()

power$cutPower ## 7
net=mg_WGCNA_getModule(edata.wgcna,power = power$cutPower
                       , deepSplit = 2, mergeCutHeight = 0.2
                       , minModuleSize = 90)
length(table(net$Modules[,2])) ## 15

pdf('analysis/01_WGCNA/Fig1D.pdf',height = 4,width = 6)
plotDendroAndColors(net$Tree, net$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

fig1e=mg_barplot_point(labels = names(table(net$Modules[,2]))
                       ,values = as.numeric(table(net$Modules[,2]))
                       ,point_sizes = 2
                       ,point_cols = names(table(net$Modules[,2]))
                       ,xlab = 'Number of Genes',legend.pos = NULL)
fig1e=fig1e+theme(axis.text.y = element_text(family = "Times", face = "plain"),
                  axis.title.y = element_text(family = "Times", face = "plain"),
                  panel.background = element_rect(fill = "white",
                                                  colour = "black"),
                  plot.background = element_rect(fill = "white",
                                                 colour = "black")
)
fig1e 

savePDF('analysis/01_WGCNA/Fig1E.pdf',fig1e,height = 6,width = 4)

#### 模块特征向量聚类
# Calculate eigengenes
MEs = net$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('analysis/01_WGCNA/Fig1F.pdf',height = 6,width = 10,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()

####################################
datTraits0 = data.frame(tcga.pheno,stringsAsFactors = F)
datTraits0$tissue=factor(datTraits0$tissue,c("Normal","Tumor"))
str(datTraits0)

######
# datTraits=datTraits0[,c("tissue",'age','ROS')]
datTraits=datTraits0
datTraits[,c("tissue")]=sapply(datTraits[, c("tissue")], function(x)as.numeric(as.factor(x)))
head(datTraits)
dim(datTraits) ## 572*2

####### Calculate module eigengenes
MEs<-net$MEs
dim(MEs) ## 572*15
## Define numbers of genes and samples
nGenes = ncol(edata.wgcna)
nSamples = nrow(edata.wgcna)
## 计算module和clinical data的相关性.
modTraitCor = WGCNA::cor(MEs[,rownames(MEDiss)[METree$order]]
                         , datTraits$tissue
                         , use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, nSamples)
## 对correlation进行可视化:
textMatrix = paste(signif(modTraitCor, 2), "\n(", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

dev.off()
pdf('analysis/01_WGCNA/Fig1G.pdf',width = 6,height = 7)
labeledHeatmap(Matrix = data.frame(modTraitCor), 
               # xLabels = colnames(modTraitCor),
               xLabels='Tissue Type',
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = rownames(modTraitCor), 
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix), setStdMargins = FALSE,
               cex.text = 0.5, zlim = c(-1,1),xLabelsAngle =0,xLabelsAdj =0.5,
               main = paste("Module-trait relationships"))
dev.off()
rownames(modTraitCor)

dim(edata.wgcna) ## 572*13325
dim(net$MEs) ## 572 * 15
geneModuleMembership <- as.data.frame(signedKME(edata.wgcna, data.frame(net$MEs), outputColumnName = ""))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

#计算基因和clinical data的相关性
all(rownames(edata.wgcna)==(datTraits$Sample))
geneTraitSignificance <- as.data.frame(cor(edata.wgcna, data.frame(tissue=datTraits$tissue), use = 'pairwise.complete.obs'))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
head(geneTraitSignificance)

###模块内分析：鉴定具有高GS和MM的基因
### Intramodular analysis: identifying genes with high GS and MM
############
get_module_hub_genes=function(net=NULL,module=NULL,trait=NULL,output=NULL,MM=0.6,GS=0.5,pval=0.05){
  modNames = substring(names(net$MEs), 3)
  moduleColors = unname(net$Modules[,2])
  #########
  module = module
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  pdf(output,height = 6,width = 6)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, trait]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",trait),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                     , col = module,lwd=2)
  abline(v = MM, col = "red", lwd = 2, lty = 2)
  abline(h = GS, col = "red", lwd = 2, lty = 2)
  dev.off()
  print(colnames(geneModuleMembership)[column])
  # inds=(abs(geneModuleMembership[moduleGenes, column])>MM & MMPvalue[moduleGenes, column]<pval) & (abs(geneTraitSignificance[moduleGenes, trait])>GS & GSPvalue[moduleGenes, trait]<pval)
  inds=((geneModuleMembership[moduleGenes, column])>MM & MMPvalue[moduleGenes, column]<pval) & ((geneTraitSignificance[moduleGenes, trait])>GS & GSPvalue[moduleGenes, trait]<pval)
  hub.genes=(rownames(geneModuleMembership)[moduleGenes])[inds]
  print(length(hub.genes))
  return(hub.genes)
}
table(net$Modules[,2])
rownames(modTraitCor)
#### module-GENE
blue.module.genes=rownames(net$Modules)[net$Modules[,2]=='blue']
length(blue.module.genes) ## 2264

blue.hub.genes=get_module_hub_genes(net = net,module = "blue",trait = "tissue"
                                    ,output = "analysis/01_WGCNA/tcga.wgcna.blue.scatterplot.pdf"
                                    ,MM=0.7,GS=0.4,pval=0.05)
length(blue.hub.genes) ## 214 相关
dev.off()
################
blue.module.enrich.res=mg_clusterProfiler(blue.hub.genes)

dt=blue.module.enrich.res$GO_BP@result
dt=dt[which(dt$p.adjust<0.05),]
dim(dt) ## 349*9

# dt[grep("oxygen species",dt$Description),1:2]
# dt[grep("oxida",dt$Description),1:2]
# inds=c(dt[grep("oxygen species",dt$Description),1],dt[grep("oxida",dt$Description),1])
# length(inds)
# 
# x=blue.module.enrich.res$GO_BP@result %>% top_n(n=-8,wt = p.adjust) %>% rownames()
# indexs=x
# indexs=c(inds,x)

p1=barplot(blue.module.enrich.res$KEGG,width=0.5,showCategory=20)+labs(title = 'KEGG')
p2=barplot(blue.module.enrich.res$GO_BP,width=0.5,showCategory=20)+labs(title = 'GO_BP')
p3=barplot(blue.module.enrich.res$GO_CC,width=0.5,showCategory=20)+labs(title = 'GO_CC')
p4=barplot(blue.module.enrich.res$GO_MF,width=0.5,showCategory=20)+labs(title = 'GO_MF')

figure=mg_merge_plot(p1,p2,p3,p4,nrow = 2,ncol = 2,labels = LETTERS[3:6])
figure
savePDF('analysis/01_WGCNA/tcga.wgcna.blue.enrich.barplot.pdf',figure,height = 12,width = 14)
write.table(data.frame(Gene=blue.hub.genes,module="blue"),'analysis/01_WGCNA/tcga.wgcna.blue.hub.genes.txt',sep = "\t",quote = F)
write.table(blue.module.enrich.res$Enrich_tab,'analysis/01_WGCNA/tcga.wgcna.blue.enrich.res.txt',quote = F)

write.table(data.frame(Gene=blue.hub.genes,module="blue"),'results/Files/tcga.wgcna.blue.hub.genes.txt',sep = "\t",quote = F)
write.table(blue.module.enrich.res$Enrich_tab,'results/Files/tcga.wgcna.blue.enrich.res.txt',quote = F)

##############################################################
############################### DEGs analysis
##############################################################
################  DEGs analysis
table(tcga.pheno$tissue)
# Tumor Normal 
# 513     59 
res=mg_limma_DEG(exp = tcga.exp,group = tcga.pheno$tissue,ulab = "Tumor",dlab = "Normal")
res$Summary
# 1.2-fold    1.3-fold    1.5-fold    2-fold     
# p<0.05   "5533|4580" "4229|3832" "2575|2812" "1007|1593"
# p<0.01   "5355|4456" "4175|3799" "2560|2804" "1005|1593"
# FDR<0.05 "5496|4559" "4220|3829" "2575|2812" "1007|1593"
# FDR<0.01 "5283|4414" "4142|3784" "2557|2802" "1005|1593"
tcga.degs=res$DEG
head(tcga.degs)

library(dplyr)
tcga.up=tcga.degs %>% filter(logFC > log2(2) & adj.P.Val<0.05) %>% rownames()
tcga.down=tcga.degs %>% filter(logFC < (-log2(2)) & adj.P.Val<0.05) %>% rownames()
length(tcga.up) ## 1007
length(tcga.down) ## 1593

df1=data.frame(Gene=tcga.up,Type='Up')
df2=data.frame(Gene=tcga.down,Type='Down')
df=rbind(df1,df2)
write.table(df,file = 'analysis/02_DEGs/tcga.degs.txt',sep = "\t",quote = F)
write.table(df,file = 'results/Files/tcga.degs.txt',sep = "\t",quote = F)

tcga.up.enrich.res=mg_clusterProfiler(tcga.up)
tcga.down.enrich.res=mg_clusterProfiler(tcga.down)
write.table(tcga.up.enrich.res$GO_BP@result,file = 'analysis/02_DEGs/tcga.up.enrich.res.txt',sep = "\t",quote = F)
write.table(tcga.down.enrich.res$GO_BP@result,file = 'analysis/02_DEGs/tcga.down.enrich.res.txt',sep = "\t",quote = F)

write.table(tcga.up.enrich.res$GO_BP@result,file = 'results/Files/tcga.up.enrich.res.txt',sep = "\t",quote = F)
write.table(tcga.down.enrich.res$GO_BP@result,file = 'results/Files/tcga.down.enrich.res.txt',sep = "\t",quote = F)

##############
library(ggplot2)
p1=barplot(tcga.up.enrich.res$GO_BP,width=0.5,showCategory=20,font.size=10)
p2=barplot(tcga.down.enrich.res$GO_BP,width=0.5,showCategory=20,font.size=10)
# p1=p1+  scale_x_discrete(
#   labels = function(x) stringr::str_wrap(x, width = 50),
#   drop = FALSE
# )

figure=mg_merge_plot(p1,p2,nrow = 1,ncol = 2,labels = LETTERS[1:2])
figure
savePDF('analysis/02_DEGs/tcga.degs.enrich.barplot.pdf',figure,height = 8,width = 16)
############### 火山图
tcga.degs.volcano=mg_volcano(logfc = tcga.degs$logFC,pvalue = tcga.degs$adj.P.Val,
                             colors = c(mycolors[3], 'grey', mycolors[6]),
                             cutFC = log2(2),
                             cutPvalue = 0.05,
                             symbol = rownames(tcga.degs),
                             leg = 'TCGA (Tumor vs Normal)',
                             legend.pos = 'tl',
                             ylab = '-log10(FDR)',
                             xlab = 'log2(FC)')
tcga.degs.volcano
savePDF('analysis/02_DEGs/tcga.degs.volcano.pdf',tcga.degs.volcano,height = 8,width = 8)

########################## 热图
library(ComplexHeatmap)
dt=tcga.exp[c(tcga.up,tcga.down),tcga.pheno$Sample]
pdf('analysis/02_DEGs/tcga.degs.heatmap.pdf',height = 8,width = 8)
Heatmap(as.matrix(t(scale(t(dt))))
        , name = "Expres"
        , row_split = c(rep('Up',length(tcga.up)),rep("Down",length(tcga.down)))
        , cluster_rows = T
        , cluster_row_slices = T
        , row_title_gp = gpar(fill = mycolors[c(6,3)])
        , show_row_dend = F
        , column_split = tcga.pheno$tissue
        , column_title_gp = gpar(fill = clin.color[c(2,1)])
        , cluster_columns = F
        , cluster_column_slices = T
        , show_column_dend = F
        , show_column_names = F
        , show_row_names = F
        , row_names_gp = gpar(fontsize = 10)
        , col = circlize::colorRamp2(c(-4, 0, 4), c('#3B4992FF', 'white', '#EE0000FF'))
        , border = TRUE)

dev.off()

##################
res=mg_limma_DEG(exp = tcga.exp,group = tcga.pheno$tissue,ulab = "Tumor",dlab = "Normal")
res$Summary
degs.res=res$DEG
degs.res$SYMBOL=rownames(degs.res)
head(degs.res)
degs.res=degs.res[,c('SYMBOL','logFC')]
head(degs.res)

library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
gene<-str_trim(degs.res$SYMBOL,"both") #定义gene
#开始ID转换
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- merge(degs.res,gene,by="SYMBOL")
head(gene_df)

geneList<-gene_df$logFC #第二列可以是folodchange，也可以是logFC
names(geneList)=gene_df$ENTREZID #使用转换好的ID
geneList=sort(geneList,decreasing = T) #从高到低排序
head(geneList)

geneList<-gene_df$logFC #第二列可以是folodchange，也可以是logFC
names(geneList)=gene_df$SYMBOL #使用转换好的ID
geneList=sort(geneList,decreasing = T) #从高到低排序
head(geneList)
########## gmt 文件
kegmt<-read.gmt("raw_datas/h.all.v2023.1.Hs.symbols.gmt") #读gmt文件
# selected=unique(kegmt$ont)[grep("OXIDATIVE_STRESS",unique(kegmt$ont))]
# kegmt=kegmt[kegmt$ont %in% selected, ]
unique(kegmt$ont) ## 50 个通路
head(kegmt)

KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
write.table(KEGG@result,'analysis/03_GSEA/tcga.GSEA.results.txt',sep = "\t",quote = F)

dt=KEGG@result
table(dt$p.adjust<0.05)
inds=dt$Description[which(dt$NES>0)]
inds

library(enrichplot)
#特定通路作图
p=gseaplot2(KEGG,inds,color=mycolors[1:11],base_size =10,pvalue_table = F) # 按第一个做二维码图，并显示p值
p
savePDF('analysis/03_GSEA/tcga.GSEA.pdf',p,height = 8,width = 14)

#####################
length(tcga.up) ## 1007
length(blue.hub.genes) ## 214

dev.off()
mg_venn_plot(data_list = list('Blue module'=blue.hub.genes,
                              'Up-regulated'=tcga.up),fill = c("blue",mycolors[6]))

alist=list('Blue module'=blue.hub.genes,"Up-regulated"=tcga.up)
p=plot(eulerr::venn(alist),labels = list(col = "gray20", font = 2), 
       edges = list(col = "gray60", lex = 1),
       fills = list(fill = c("blue",mycolors[6]), alpha = 0.6),
       quantities = list(cex = .8, col = 'gray20'))
p
savePDF('analysis/04_diagnostic_model/tcga.wgcna.degs.venn.plot.pdf',p,height = 5,width = 5)

####
hub.genes=Reduce(intersect,list(blue.hub.genes,c(tcga.up)))
length(hub.genes) ## 192

################################# 诊断模型
head(tcga.pheno)
all(colnames(tcga.exp)==(tcga.pheno$Sample))
expr_mat=t(tcga.exp[hub.genes,tcga.pheno$Sample])
dim(expr_mat) ## 572*192
tcga.pheno$tissue=factor(tcga.pheno$tissue,levels = c('Normal','Tumor'))

############################################################
####################### lasso analysis
library(glmnet)
set.seed(123)
fit=glmnet(x = expr_mat,y = tcga.pheno$tissue,family = "binomial"
           , nlambda = 100
           , alpha = 1) 
set.seed(123)
cv.fit<-cv.glmnet(x = expr_mat,y = tcga.pheno$tissue,family = "binomial"
                  , nlambda = 100
                  , nfolds = 10
                  , alpha = 1)
cv.fit$lambda.min #0.002637551

options(ggrepel.max.overlaps = Inf)
pdf('analysis/04_diagnostic_model/lasso.pdf',height = 5,width = 10,onefile = F)
par(mfrow=c(1,2))
plot(fit,xvar="lambda",label=T)
plot(cv.fit)
dev.off()

lambda=cv.fit$lambda.min
coefficients<-coef(fit,s=lambda)
Active.Index<-which(coefficients[,1]!=0)
genes=row.names(coefficients)[Active.Index]
Active.coefficients<-coefficients[Active.Index]  
###
genes=genes[-1]
Active.coefficients=Active.coefficients[-1]

LASSO_genes=genes
LASSO_genes
# [1] "SLC2A1"   "IQGAP3"   "UBE2T"    "E2F8"     "TEDC2"    "RMI2"     "ARHGEF39"
# [8] "RCC1"     "SRPK1"    "SPATS2"   "FAM136A"

############################################################
########### 支持向量机方法(SVM)
library(caret)
library(kernlab)
######## Controlling the Feature Selection Algorithms
ctrl = rfeControl(functions = caretFuncs, method ="cv", number= 10, verbose = FALSE)
######## Backwards Feature Selection
sizes= c(2:(ncol(expr_mat)-1))
sizes
sizes=c(2:50)
set.seed(123)
svm_model = caret::rfe(x = expr_mat,y =tcga.pheno$tissue, sizes =sizes,
                       rfeControl = ctrl, method ="svmLinear")
svm_model
# The top 5 variables (out of 18):
#   SAPCD2, ERCC6L, TEDC2, UBE2T, RCC1

optsize=svm_model$optsize
results=svm_model$results
accuracy=results$Accuracy[which(results$Variables==optsize)]
accuracy=round(accuracy,digits = 3)
optsize ## 18
accuracy ## 0.991

# plot the results
plot(svm_model, type=c("g", "o"))

pdf('analysis/04_diagnostic_model/SVM.pdf',height = 5,width = 5,onefile = F)
plot(svm_model$results$Variables,svm_model$results$Accuracy,type="l",col="blue",xlim=c(1,50),
     xlab="Number of features",ylab="10x CV accuracy")
points(optsize, accuracy, cex = 2, pch = 1, col ="red")
text(optsize,0.99,paste(optsize,accuracy,sep = " - "),col="red",cex=1)
dev.off()

SVM_genes=svm_model$optVariables
SVM_genes
# [1] "SAPCD2"   "ERCC6L"   "TEDC2"    "UBE2T"    "RCC1"     "FAM136A"  "PAICS"   
# [8] "TOP2A"    "SPAG5"    "CCNB1"    "HJURP"    "IQGAP3"   "CDC25C"   "ARHGEF39"
# [15] "KIF2C"    "CDCA8"    "ORC6"     "GINS1"

x=intersect(SVM_genes,LASSO_genes)
x ## "TEDC2"    "UBE2T"    "RCC1"     "FAM136A"  "IQGAP3"   "ARHGEF39" 

############################################################
############ 随机森林方法
library(randomForest)
# 1.寻找最优参数mtry，即指定节点中用于二叉树的最佳变量个数
n = ncol(expr_mat)
print(n) ## 192
err=as.numeric()     # 设置模型误判率向量初始值
for(i in 1:(n-1)){
  set.seed(123456)
  fit = randomForest(expr_mat,tcga.pheno$tissue,mtry=i)
  err[i] = mean(fit$err.rate)
}
plot(err,ylab="10x CV accuracy")

mtry=which.min(err)
mtry ## 137

# 2.寻找最佳参数ntree，即指定随机森林所包含的最佳决策树数目
set.seed(123456)
fit = randomForest(expr_mat,tcga.pheno$tissue,mtry=mtry,ntree=1000)
plot(fit)    #绘制模型误差与决策树数量关系图 
abline(v = 150, col = "blue")
abline(v = 600, col = "blue")
abline(v = 610, col = "blue")
abline(h = 0.01048951, col = "red")
abline(h = 0.06779661, col = "red")
abline(h = 0.003898635, col = "red")
summary(fit)

which.min(apply(fit$err.rate,1,function(x){min(x)}))

ntree=610

# 3.随机森林模型搭建
set.seed(123456)
RF_model = randomForest(expr_mat,tcga.pheno$tissue,mtry=mtry,ntree=ntree,importance=TRUE,proximity=TRUE)    
RF_model
##输出变量重要性:分别从精确度递减和均方误差递减的角度来衡量重要程度。
importance=RF_model$importance
importance=data.frame(importance)
dim(importance) ## 192*4
varImpPlot(RF_model)

pdf("analysis/04_diagnostic_model/RF.top10.pdf",height = 5,width = 10)
varImpPlot(RF_model, n.var = 10, main = 'Top 10 - variable importance')
dev.off()

RF_genes_1=importance %>% dplyr::arrange(desc(MeanDecreaseAccuracy))%>% head(n=10) %>% rownames()
RF_genes_2=importance %>% dplyr::arrange(desc(MeanDecreaseGini))%>% head(n=10) %>% rownames()

RF_genes=intersect(RF_genes_1,RF_genes_2)
RF_genes
######### ntree=200
# [1] "FAM136A" "UBE2T"   "XRCC2"   "SRPK1"   "MYBL2"   "UBE2C"   "CDC20"   "SLC2A1" 
# [9] "BIRC5"
######### ntree=610
# "SAPCD2"  "TEDC2"   "FAM136A" "RCC1"    "UBE2T"   "HJURP"
###############
dev.off()
alist=list('LASSO'=LASSO_genes,"SVM"=SVM_genes,'RF'=RF_genes)
model.hub.venn=plot(eulerr::venn(alist),labels = list(col = "gray20", font = 2), 
       edges = list(col = "gray60", lex = 1),
       fills = list(fill = mycolors[1:3], alpha = 0.6),
       quantities = list(cex = .8, col = 'gray20'))
model.hub.venn

########
final.genes=Reduce(intersect,list('LASSO'=LASSO_genes,"SVM"=SVM_genes,"RF"=RF_genes))
final.genes
# "UBE2T"   "TEDC2"   "RCC1"    "FAM136A"

################# 转录组诊断模型
###############################
library(pROC)
final.genes
gene.colors=mycolors[1:4]
names(gene.colors)=final.genes

T_model_data=data.frame(expr_mat[,final.genes],tissue=tcga.pheno$tissue)
dim(T_model_data) ## 572*5

p.all <- list()
p.all <- c(p.all,list(model.hub.venn))
for (gene in final.genes) {
  dt <- T_model_data[, c('tissue', gene)]
  colnames(dt) <- c('disease', 'Gene')
  roc.res <- roc(disease ~ Gene,
                 data=dt,
                 aur=TRUE,
                 ci=TRUE, 
                 smooth=F)
  # roc 
  p <- ggroc(roc.res,color=gene.colors[gene], legacy.axes = TRUE )+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color="darkgrey", linetype=8)+
    # theme_gray() +
    theme_bw()+
    annotate("text",x=0.75,y=0.25,label=paste("AUC = ", round(roc.res$auc,3))) +
    ggtitle(paste0(gene, ' ROC'))
  p.all[[gene]] <- p
  
}
length(p.all)
hub.genes.roc <- cowplot::plot_grid(plotlist = p.all,ncol = 3)
hub.genes.roc

######################## 合并到一起
# This is equivalent to using roc.formula:
head(T_model_data)
final.genes
roc.list <- roc(tissue ~ UBE2T+TEDC2+RCC1+FAM136A, data = T_model_data)

auc.res=sapply(roc.list,function(x){round(x$auc,3)},simplify=T)
auc.res
data.labels=paste0(names(auc.res)," , AUC = ",paste(round(auc.res,2)))

ggroc(roc.list, legacy.axes = TRUE )+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="darkgrey", linetype=8)+
  scale_colour_manual(labels=data.labels,values = gene.colors)+
  theme_bw()+labs(title = "TCGA")+
  theme(legend.position = c(0.8, 0.2))


# ################### SVM model
# library(pROC)
# library(e1071)
# paste0(final.genes,collapse = '+')
# dim(T_model_data)
# fit =svm(tissue ~ . ,T_model_data, probability=TRUE)
# summary(fit)
# pred <- predict(fit)
# pred <- as.ordered(pred)
# model.roc <- roc(T_model_data$tissue,pred)
# ########################### 
# pdf('analysis/06_Model/svm.model.roc.pdf',height =5,width = 5)
# plot(model.roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), 
#      grid.col=c("blue", "red"),
#      max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=TRUE)
# dev.off()

######### logit模型
dim(T_model_data)
fit2 <- caret::train(tissue ~ . , data = T_model_data,
                     method = 'glmnet',
                     family = 'binomial')
fit2
pred <- predict(fit2, T_model_data[,1:4],'raw')
pred <- as.ordered(pred)
pred
model.roc <- roc(T_model_data$tissue,pred)
########################### 
pdf('analysis/04_diagnostic_model/logit.model.roc.pdf',height =5,width = 5)
plot(model.roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), 
     grid.col=c("blue", "red"),
     max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()

######################
library(ggpubr)
p.all=list()
for(gene in final.genes){
  dt=T_model_data[,c("tissue",gene)]
  dt$tissue=factor(dt$tissue,levels = c("Normal","Tumor"))
  p=ggboxplot(dt,x='tissue',y=colnames(dt)[2],fill = "tissue", palette = clin.color,main=gene,ylab="Expression",xlab="")+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none",text = element_text(family = 'Times'))+
    stat_compare_means(method = "wilcox.test",label.x=1.5,label="p.signif")
  p
  p.all=c(p.all,list(p))
}
hub.genes.boxplot <- cowplot::plot_grid(plotlist = p.all,nrow = 1)
hub.genes.boxplot

figure=mg_merge_plot(hub.genes.roc,hub.genes.boxplot,nrow = 2,heights = c(2,1))
figure

savePDF('analysis/04_diagnostic_model/tcga.T_model.hub.genes.roc.pdf',figure,height = 15,width = 15)


################# 外部数据集-验证
load('raw_datas/GEO/gse30219.exp.RData')
gse30219.exp=gse30219.exp[rownames(gse30219.exp) %in% gencode.pcg$gene_name,]
dim(gse30219.exp) ## 17283*307
########## 临床信息整理
load("raw_datas/GEO/GSE30219.RData")
gse30219.cli<-GSE30219$Sample
rownames(gse30219.cli)=gse30219.cli$Acc
gse30219.cli$Title=gsub("(.*?)\\d+","\\1",gse30219.cli$Title,perl = TRUE)
gse30219.cli$Title=gsub("(.*?) ADK$","\\1 ADC",gse30219.cli$Title)
table(gse30219.cli$Title,gse30219.cli$Source)
dim(gse30219.cli)
table(gse30219.cli$tissue,gse30219.cli$histology)
table(gse30219.cli$Source,gse30219.cli$histology)
gse30219.cli=reshape::rename(gse30219.cli,c("Source"='tissue'))
gse30219.cli$tissue[which(gse30219.cli$tissue=="Lung Tumour")]="Tumor"
gse30219.cli$tissue[which(gse30219.cli$tissue=="Non Tumoral Lung")]="Normal"
gse30219.cli$tissue=factor(gse30219.cli$tissue,levels = c("Normal","Tumor"))
table(gse30219.cli$tissue)

test_data1=data.frame(tissue=gse30219.cli$tissue,t(gse30219.exp[final.genes,]))
p.all <- list()
for (gene in final.genes) {
  dt <- test_data1[, c('tissue', gene)]
  colnames(dt) <- c('disease', 'Gene')
  roc.res <- roc(disease ~ Gene,
                 data=dt,
                 aur=TRUE,
                 ci=TRUE, 
                 smooth=F)
  # roc 
  p <- ggroc(roc.res,color=gene.colors[gene], legacy.axes = TRUE )+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color="darkgrey", linetype=8)+
    # theme_gray() +
    theme_bw()+
    annotate("text",x=0.75,y=0.25,label=paste("AUC = ", round(roc.res$auc,3))) +
    ggtitle(paste0(gene, ' ROC'))
  p.all[[gene]] <- p
  
}
length(p.all)
cowplot::plot_grid(plotlist = p.all,ncol = 4)

#######
final.genes
roc.list <- roc(tissue ~ UBE2T+TEDC2+RCC1+FAM136A, data = test_data1)

auc.res=sapply(roc.list,function(x){round(x$auc,3)},simplify=T)
auc.res
data.labels=paste0(names(auc.res)," , AUC = ",paste(round(auc.res,2)))

gse30219.hub.genes.roc=ggroc(roc.list, legacy.axes = TRUE )+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="darkgrey", linetype=8)+
  scale_colour_manual(labels=data.labels,values = gene.colors)+
  theme_bw()+labs(title = "GSE30219",color="Gene")+
  theme(legend.position = c(0.8, 0.2))
gse30219.hub.genes.roc


p.all=list()
for(gene in final.genes){
  dt=test_data1[,c("tissue",gene)]
  p=ggboxplot(dt,x='tissue',y=colnames(dt)[2],fill = "tissue", palette = clin.color,main=gene,ylab="Expression",xlab="")+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none",text = element_text(family = 'Times'))+
    stat_compare_means(method = "wilcox.test",label.x=1.5,label="p.signif")
  p
  p.all=c(p.all,list(p))
}
length(p.all)
gse30219.hub.genes.boxplot <- cowplot::plot_grid(plotlist = p.all,ncol = 4)
gse30219.hub.genes.boxplot
head(test_data1)
gse30219.hub.genes.boxplot=mg_PlotMutiBoxplot(test_data1[,-1]
                                              , group = test_data1$tissue
                                              , legend.pos = 'top'
                                              , add = 'boxplot'
                                              , xangle = 45
                                              , ylab = 'Expression'
                                              , group_cols = clin.color
                                              , fill = T
                                              , test_method = 'wilcox.test') + labs(title = 'GSE30219',fill = 'Tissue')
                     
gse30219.hub.genes.boxplot

###########
load('raw_datas/GEO/gse31210.exp.RData')
gse31210.exp=gse31210.exp[rownames(gse31210.exp) %in% gencode.pcg$gene_name,]
dim(gse31210.exp) ## 17283*246
########## 临床信息整理
load('raw_datas/GEO/GSE31210.RData')
gse31210.cli<-GSE31210$Sample
rownames(gse31210.cli)=gse31210.cli$Acc
table(gse31210.cli$tissue)
gse31210.cli$tissue[which(gse31210.cli$tissue=="primary lung tumor")]="Tumor"
gse31210.cli$tissue[which(gse31210.cli$tissue=="normal lung")]="Normal"
gse31210.cli$tissue=factor(gse31210.cli$tissue,levels = c("Normal","Tumor"))
table(gse31210.cli$tissue)
# Normal  Tumor 
# 20    226
###########
test_data2=data.frame(tissue=gse31210.cli$tissue,t(gse31210.exp[final.genes,]))
head(test_data2)
p.all <- list()
for (gene in final.genes) {
  dt <- test_data2[, c('tissue', gene)]
  colnames(dt) <- c('disease', 'Gene')
  roc.res <- roc(disease ~ Gene,
                 data=dt,
                 aur=TRUE,
                 ci=TRUE, 
                 smooth=F)
  # roc 
  p <- ggroc(roc.res,color=gene.colors[gene], legacy.axes = TRUE )+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color="darkgrey", linetype=8)+
    # theme_gray() +
    theme_bw()+
    annotate("text",x=0.75,y=0.25,label=paste("AUC = ", round(roc.res$auc,3))) +
    ggtitle(paste0(gene, ' ROC'))
  p.all[[gene]] <- p
  
}
length(p.all)
cowplot::plot_grid(plotlist = p.all,ncol = 4)


final.genes
roc.list <- roc(tissue ~ UBE2T+TEDC2+RCC1+FAM136A, data = test_data2)

auc.res=sapply(roc.list,function(x){round(x$auc,3)},simplify=T)
auc.res
data.labels=paste0(names(auc.res)," , AUC = ",paste(round(auc.res,2)))

gse31210.hub.genes.roc=ggroc(roc.list, legacy.axes = TRUE )+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="darkgrey", linetype=8)+
  scale_colour_manual(labels=data.labels,values = gene.colors)+
  theme_bw()+labs(title = "GSE31210",color="Gene")+
  theme(legend.position = c(0.8, 0.2))
gse31210.hub.genes.roc


p.all=list()
for(gene in final.genes){
  dt=test_data2[,c('tissue',gene)]
  p=ggboxplot(dt,x='tissue',y=colnames(dt)[2],fill = "tissue", palette = clin.color,main=gene,ylab="Expression",xlab="")+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none",text = element_text(family = 'Times'))+
    stat_compare_means(method = "wilcox.test",label.x=1.5,label="p.signif")
  p
  p.all=c(p.all,list(p))
}
length(p.all)
gse31210.hub.genes.boxplot <- cowplot::plot_grid(plotlist = p.all,ncol = 4)
gse31210.hub.genes.boxplot

gse31210.hub.genes.boxplot=mg_PlotMutiBoxplot(test_data2[,-1]
                                              , group = test_data2$tissue
                                              , legend.pos = 'top'
                                              , add = 'boxplot'
                                              , xangle = 45
                                              , ylab = 'Expression'
                                              , group_cols = clin.color
                                              , fill = T
                                              , test_method = 'wilcox.test') + labs(title = 'GSE31210',fill = 'Tissue')

gse31210.hub.genes.boxplot

plot1=mg_merge_plot(gse30219.hub.genes.roc,gse30219.hub.genes.boxplot,nrow = 2,ncol = 1,widths = c(1,1),labels = LETTERS[1:2])
plot2=mg_merge_plot(gse31210.hub.genes.roc,gse31210.hub.genes.boxplot,nrow = 2,ncol = 1,widths = c(1,1),labels = LETTERS[3:4])

figure=mg_merge_plot(plot1,plot2,nrow = 1,ncol =2)
figure
savePDF('analysis/05_model_validation/T_model.validation.roc.pdf',figure,height = 10,width = 10)

#####################################################
ggplotTimeROC=function(time,status,score,mks=c(1,3,5),pal=NULL){
  #time=g.os
  #status=g.ev
  #score=as.numeric(cpm.score)
  #cx=coxRun(data.frame(time,status,score))
  #if(cx[1]<=1){
  #  score=-1*score
  #}
  roc.tm=mg_surv_pROC(time,status,score,mks)
  print('roc.tm')
  print((roc.tm))
  library(survival)
  library(ggplot2)
  mks=mg_predict_time_ymd(time,mks)
  print(mks)  
  ROC.DSST=timeROC::timeROC(T=time,
                            delta=status
                            ,marker=score,
                            cause=1,weighting="marginal",
                            times=mks,
                            iid=TRUE)
  print(ROC.DSST)
  mks=mks[which(!is.na(ROC.DSST$AUC)&ROC.DSST$AUC>0)]
  print(mks)
  if(length(mks)>0){
    if(max(ROC.DSST$AUC)<0.5){
      score=-1*score
    }
    ROC.DSST=timeROC::timeROC(T=time,
                              delta=status
                              ,marker=score,
                              cause=1,weighting="marginal",
                              times=mks,
                              iid=TRUE)
    print(ROC.DSST$times)
    if(max(ROC.DSST$times)<20){
      lb=paste0(ROC.DSST$times,'-Years')
    }else if(max(ROC.DSST$times)<365){
      lb=paste0(round(ROC.DSST$times/12,0),'-Years')
    }else{
      lb=paste0(round(ROC.DSST$times/365,0),'-Years')
    }
    
    lbs=paste0(lb,',AUC=',round(ROC.DSST$AUC,2),',95%CI(',paste0(round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,1]/100,2),'-',
                                                                 round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,2]/100,2)),')')
    ########
    p.dat=rbind()
    print(length(roc.tm))
    for(i in 1:length(roc.tm)){
      #print(i)
      r1=roc.tm[[i]]
      x1=1-r1$specificities
      y1=r1$sensitivities
      #print(cbind(1-r1$specificities,r1$sensitivities))
      nx1=unique(x1)
      ny1=c()
      for(x in unique(x1)){
        x.inds=which(x1==x)
        if(length(x.inds)>0&x<0.5){
          ny1=c(ny1,min(y1[x.inds]))
        }else if(length(x.inds)>0){
          ny1=c(ny1,max(y1[x.inds]))
        }else{
          ny1=c(ny1,y1[x.inds][1])
        }
      }
      #print(cbind(nx1,ny1))
      p.dat=rbind(p.dat,data.frame(x=nx1, y=ny1,rep(lbs[i],length(nx1)),stringsAsFactors = F))
    }
    colnames(p.dat)=c('V1','V2','Type')
    p.dat=as.data.frame(p.dat)
    
    ###############################
    if(is.null(pal)){
      pal=c(pal_npg(alpha =0.8)(9)[c(3,4,1,9)]) ## 4种
    }
    
    custome_theme=function (base_size = 12, base_family = "", border = FALSE, 
                            margin = TRUE, legend = c("top", "bottom", "left", 
                                                      "right", "none"), x.text.angle = 0) 
    {
      half_line <- base_size/2
      if (!is.numeric(legend)) 
        legend <- match.arg(legend)
      if (x.text.angle > 5) 
        xhjust <- 1
      else xhjust <- NULL
      if (border) {
        panel.border <- element_rect(fill = NA, colour = "black", 
                                     size = 0.7)
        axis.line <- element_blank()
      }
      else {
        panel.border <- element_blank()
        axis.line = element_line(colour = "black", size = 0.5)
      }
      if (margin) 
        plot.margin <- ggplot2::margin(half_line, half_line, half_line, 
                              half_line)
      else plot.margin <- unit(c(0.5, 0.3, 0.3, 0.3), "mm")
      .theme <- theme_bw(base_size = base_size, base_family = base_family) %+replace% 
        theme(panel.border = panel.border, panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), axis.line = axis.line, 
              axis.text = element_text(color = "black"), 
              legend.key = element_blank(), strip.background = element_rect(fill = "#F2F2F2", 
                                                                            colour = "black", size = 0.7), plot.margin = plot.margin, 
              legend.position = legend, complete = TRUE)
      if (x.text.angle != 0) 
        .theme <- .theme + theme(axis.text.x = element_text(angle = x.text.angle, 
                                                            hjust = xhjust))
      .theme
    }
    
    p1=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))
    p1=p1+geom_line(aes(colour=Type),lwd=1.1)+custome_theme(base_size=12,base_family="Times",border=T,legend=c(0.7,0.2))
    p1=p1+xlab('False positive fraction')+ylab('True positive fraction')
    # p1=p1+stat_smooth(aes(colour=Type),se = FALSE, size = 1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction')
    p1=p1+theme(legend.background = element_rect(fill = NA, colour = NA))
    p1=p1+scale_colour_manual(values = pal)
    return(p1)
  }else{
    return(mg_getplot_bank('No data plot by ROC!'))
  }
}

## 无虚线
ggplotKMCox=function(dat,title=NULL,legend.title='Groups',labs=NULL,add_text=NULL,pal=NULL){
  library(ggplot2)
  library(ggsci)
  library(survival)
  library(ggpubr)
  library(survminer)
  ######
  colnames(dat)=c('time','status','groups')
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  if(is.null(pal)){
    pal=pal_lancet()(9)[c(2,4,3,1,5:6,9)]
  }
  surv=survminer::ggsurvplot(sf, data = dat
                             , palette = pal 
                             , pval = TRUE
                             , surv.median.line = 'hv'
                             , conf.int = T
                             # , linetype = "strata"
                             , title=title
                             , legend.title=legend.title
                             , xlab = "Time(years)"
                             , conf.int.style = 'step'
                             , pval.coord = c(0, 0.2)#Add p-value
                             , risk.table = TRUE
                             , ggtheme = theme_pubr(base_size=12,base_family='Times')
                             , font.family='Times'
                             , risk.table.y.text = FALSE
                             , legend.labs = labs)
  p1=surv$plot
  p2=surv$table
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(0.9,0.3),align = "v")
  return(g2)
}
####################
final.genes
# "UBE2T"    "TEDC2"   "RCC1"     "FAM136A" 
load("raw_datas/TCGA/tcga.t.cli_use.RData")
dim(tcga.t.cli_use) ## 513*13
dim(tcga.t.exp) ## 19037*513
all(rownames(tcga.t.cli_use)==colnames(tcga.t.exp))
head(tcga.t.cli_use)

y=apply(t(tcga.t.exp[final.genes,]),2,function(x){ifelse(x>median(x),"High","Low")})
dt=cbind(tcga.t.cli_use,y)

p.all=list()
for(gene in final.genes){
  p=ggplotKMCox(data.frame(time = dt$OS.time/365
                            , event = dt$OS
                            , groups=dt[,gene])
                 , legend.title = gene
                 , title = 'TCGA'
                 , pal = clin.color
                 , labs = c('High','Low')
                 , add_text = '')
  p.all=c(p.all,list(p))
  
}
length(p.all)
hub.genes.km <- cowplot::plot_grid(plotlist = p.all,ncol = 4)
hub.genes.km
#############
dt=cbind(tcga.t.cli_use,t(tcga.t.exp[c(final.genes),]))
dt=dt[!is.na(dt$OS),]
dt=dt[which(dt$OS.time>30),]
head(dt)

p.all=list()
for(gene in final.genes){
  p=ggplotTimeROC(time = dt$OS.time,
                  status = dt$OS,
                  score = dt[,gene])+labs(title = gene)
  p.all=c(p.all,list(p))
  
}

hub.genes.timeROC <- cowplot::plot_grid(plotlist = p.all,ncol = 4)
hub.genes.timeROC

############################# 与临床信息的关系
######## 可视化
mg_Forestplot=function(df_m,outFile,width=6,height=3){
  colnames(df_m)=c('HR','HR.95L','HR.95H','pvalue')
  gene=rownames(df_m)
  # hr=sprintf("%.3f",df_m$"HR")
  # hrLow=sprintf("%.3f",df_m$"HR.95L")
  # hrHigh=sprintf("%.3f",df_m$"HR.95H")
  hr=format(df_m$"HR", digits = 3)
  hrLow=format(df_m$"HR.95L", digits = 3)
  hrHigh=format(df_m$"HR.95H",digits =3)
  Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
  # pVal=ifelse(df_m$pvalue<1e-8, "<1e-8", sprintf("%.3f", df_m$pvalue))
  pVal=format(df_m$pvalue, digits = 3)
  #########
  pdf(file=outFile, width = width, height =height,onefile = FALSE)
  n=nrow(df_m)
  nRow=n+1
  ylim=c(1,nRow)
  layout.show(layout(matrix(c(1,2),nc=2),width=c(2,1.2)))
  #森林图左边的基因信息
  xlim = c(0,3)
  par(mar=c(4,1,2,1),mpg=c(2,0.5,0))
  # par(mar=c(3,2,1.5,1.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mpg=c(2,0.5,0))
  # par(mar=c(3,1,1.5,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="black",lwd=2.5,lty=5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
  points(as.numeric(hr), n:1, pch = 20, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}

unicox<-function(vars=c("T","N"),time=NULL,event=NULL,data=LUAD_clinical){
  library(survival)
  require(survminer)
  y<-Surv(as.numeric(time),
          as.numeric(event))
  
  ## for循环
  pvalue<-c()
  HR<-c()
  lower<-c()
  upper<-c()
  varname<-c()
  
  for(i in vars ){
    cox.fit_uni<-coxph(y~data[[i]])
    uni_res <-summary(cox.fit_uni)
    pvalue[[i]]<-uni_res$waldtest[[3]]#pvalue
    HR[[i]]<-uni_res$conf.int[[1]]# HR提取
    lower[[i]] <-uni_res$conf.int[[3]]#lower
    upper[[i]] <-uni_res$conf.int[[4]]# upper
  }
  ## 保存结果到数据框中
  univar_res<-data.frame(
    # varname<-vars,
    HR<-as.numeric(HR),
    # CI=paste(as.numeric(lower),"-",as.numeric(upper),sep = ""),
    HR.95L<-as.numeric(lower),
    HR.95H<-as.numeric(upper),
    pvalue<-as.numeric(pvalue)
  )
  colnames(univar_res)<-c("HR","HR.95L","HR.95H","pvalue")
  rownames(univar_res)=vars
  univar_res## 最后的函数运行结果
}
multicox<-function(vars=c("T","N","M","age"),time=NULL,event=NULL,data=LUAD_clinical,forest=T){
  library(survival)
  require(survminer)
  y<-Surv(as.numeric(time),
          as.numeric(event))
  ## 构建公式
  FM<-as.formula(paste0("y~",paste(vars,collapse = "+")))
  cox.fit_multi <- coxph(FM,data = data)
  munivar_res<-summary(cox.fit_multi)#cox的结果
  pvalue<-munivar_res$coefficients[,"Pr(>|z|)"]#pvalue
  HR<-munivar_res$coefficients[,"exp(coef)"]# HR
  lower<-munivar_res$conf.int[,3]
  upper<-munivar_res$conf.int[,4]
  ## 保存结果到数据框中
  munivar_res<-data.frame(
    # varname<-vars,
    HR<-as.numeric(HR),
    # CI=paste(as.numeric(lower),"-",as.numeric(upper),sep = ""),
    HR.95L<-as.numeric(lower),
    HR.95H<-as.numeric(upper),
    pvalue<-as.numeric(pvalue)
  )
  colnames(munivar_res)<-c("HR","HR.95L","HR.95H","pvalue")
  rownames(munivar_res)=vars
  
  ## 保存pdf到工作目录
  if (forest==T){
    ggforest(cox.fit_multi,data = data)
    ggsave(filename = "mult_cox.pdf")}
  munivar_res## 最后的函数运行结果
  
}

## 批量进行单因素cox回归分析
tcga_clinical=cbind(tcga.t.cli_use,t(tcga.t.exp[final.genes,]))
table(tcga_clinical$`T Stage`) # T1:163,T2:263,T3:43,T4:18
table(tcga_clinical$`N Stage`) # N0:317,N1:92,N2:68,N3:2
table(tcga_clinical$`M Stage`) # M0:324,M1:24
table(tcga_clinical$Stage) # I:263,II:115,III:79,IV:25
table(tcga_clinical$Gender) # F:262,M:228
table(tcga_clinical$Smoker) # Never:68,Ever:292,Current:116

tcga_clinical$T.Stage=ifelse(!is.na(tcga_clinical$`T Stage`),ifelse(tcga_clinical$`T Stage` %in% c('T1','T2'),'T1+T2','T3+T4'),NA)
tcga_clinical$N.Stage=ifelse(!is.na(tcga_clinical$`N Stage`),ifelse(tcga_clinical$`N Stage` %in% c('N0'),'N0','N1+N2+N3+N4'),NA)
tcga_clinical$M.Stage=ifelse(!is.na(tcga_clinical$`M Stage`),ifelse(tcga_clinical$`M Stage` %in% c('M0'),'M0','M1'),NA)
tcga_clinical$Stage_2=ifelse(!is.na(tcga_clinical$Stage),ifelse(tcga_clinical$Stage %in% c('I','II'),'I+II','III+IV'),NA)
tcga_clinical$Smoker_2=ifelse(!is.na(tcga_clinical$Smoker),ifelse(tcga_clinical$Smoker %in% c('Ever','Current'),'Yes','No'),NA)

tcga_clinical$Gender=factor(tcga_clinical$Gender,levels = c('FEMALE','MALE'))
tcga_clinical$Smoker_2=factor(tcga_clinical$Smoker_2,levels = c("No","Yes"))

##################### OS
## 单因素
univar_res<-unicox(vars=c("T.Stage","N.Stage","M.Stage","Stage_2","Age","Gender","Smoker_2",final.genes),time = tcga_clinical$OS.time,event = tcga_clinical$OS,data=tcga_clinical)
univar_res
rownames(univar_res)[4]='Stage'
rownames(univar_res)[7]='Smoker'
table(univar_res$pvalue<0.05)
rownames(univar_res)[which(univar_res$pvalue<0.05)]
# "T.Stage" "N.Stage" "M.Stage" "Stage"   "UBE2T"   "TEDC2"   "FAM136A"
## 多因素
univar.sig=c("UBE2T","TEDC2","FAM136A")
mutivar_res<-multicox(vars=c("T.Stage","N.Stage","M.Stage","Stage_2",univar.sig),time = tcga_clinical$OS.time,event = tcga_clinical$OS,data=tcga_clinical,forest = F)
mutivar_res
rownames(mutivar_res)[4]='Stage'
rownames(mutivar_res)[mutivar_res$pvalue<0.05]
# "T.Stage"  "N.Stage"
mg_Forestplot(df_m = univar_res,outFile = 'analysis/06_Clinical/TCGA.univar.forestplot.pdf',height = 4,width = 6)
mg_Forestplot(df_m = mutivar_res,outFile = 'analysis/06_Clinical/TCGA.mutivar.forestplot.pdf',height = 4,width = 6)
##################### PFS
univar_res<-unicox(vars=c("T.Stage","N.Stage","M.Stage","Stage_2","Age","Gender","Smoker_2",final.genes),time = tcga_clinical$PFI.time,event = tcga_clinical$PFI,data=tcga_clinical)
univar_res
rownames(univar_res)[4]='Stage'
rownames(univar_res)[7]='Smoker'
rownames(univar_res)[which(univar_res$pvalue<0.05)]
# "T.Stage" "N.Stage" "Stage"   "UBE2T"
## 多因素
univar.sig=c("UBE2T")
mutivar_res<-multicox(vars=c("T.Stage","N.Stage","Stage_2",univar.sig),time = tcga_clinical$PFI.time,event = tcga_clinical$PFI,data=tcga_clinical,forest = F)
mutivar_res
rownames(mutivar_res)[3]='Stage'
rownames(mutivar_res)[which(mutivar_res$pvalue<0.05)]
# "T.Stage" "N.Stage"
mg_Forestplot(df_m = univar_res,outFile = 'analysis/06_Clinical/TCGA.univar.forestplot.PFI.pdf',height = 4,width = 6)
mg_Forestplot(df_m = mutivar_res,outFile = 'analysis/06_Clinical/TCGA.mutivar.forestplot.PFI.pdf',height = 4,width = 6)

################# 跟临床信息的相关性分析
final.genes
tcga_clinical_2=cbind(tcga.t.cli_use,t(tcga.t.exp[final.genes,]))
head(tcga_clinical_2)
tcga_clinical_2$T.Stage=tcga_clinical_2$`T Stage`
tcga_clinical_2$N.Stage=tcga_clinical_2$`N Stage`
tcga_clinical_2$M.Stage=tcga_clinical_2$`M Stage`
tcga_clinical_2$Stage_2=tcga_clinical_2$Stage
tcga_clinical_2$Smoker_2=tcga_clinical_2$Smoker

tcga_clinical_2$T.Stage=factor(tcga_clinical_2$T.Stage,levels = c('T1','T2','T3','T4'))
tcga_clinical_2$N.Stage=factor(tcga_clinical_2$N.Stage,levels = c('N0','N1','N2','N3'))
tcga_clinical_2$M.Stage=factor(tcga_clinical_2$M.Stage,levels = c('M0','M1'))
tcga_clinical_2$Stage_2=factor(tcga_clinical_2$Stage_2,levels = c('I','II','III','IV'))
tcga_clinical_2$Gender=factor(tcga_clinical_2$Gender,levels = c('FEMALE','MALE'))
tcga_clinical_2$Smoker_2=factor(tcga_clinical_2$Smoker_2,levels = c("Never","Ever","Current"))

tcga_clinical_2$T.Stage=as.numeric(tcga_clinical_2$T.Stage)
tcga_clinical_2$N.Stage=as.numeric(tcga_clinical_2$N.Stage)
tcga_clinical_2$M.Stage=as.numeric(tcga_clinical_2$M.Stage)
tcga_clinical_2$Stage_2=as.numeric(tcga_clinical_2$Stage_2)
tcga_clinical_2$Gender=as.numeric(tcga_clinical_2$Gender)
tcga_clinical_2$Smoker_2=as.numeric(tcga_clinical_2$Smoker_2)

####################
dt=tcga_clinical_2[,c(final.genes,'T.Stage','N.Stage','M.Stage','Stage_2','Gender','Age','Smoker_2')]
dt=reshape::rename(dt,c('Stage_2'='Stage','Smoker_2'='Smoker'))
res <- psych::corr.test(as.matrix(dt))
cor_M <- res$r
cor_P <- res$p

mycolor <- pal_d3(alpha =1)(7)
library(ggcorrplot)
p <- ggcorrplot(cor_M, 
                hc.order = TRUE,
                colors = c(mycolor[1], "white", mycolor[2]),
                ggtheme = ggplot2::theme_bw,
                # type = c("lower"),
                p.mat = cor_P)
p

library(corrplot)
cor.colors=pal_d3(alpha =1)(7)[c(1,3,2,4)]
scales::show_col(cor.colors)
cor.colors

pdf('analysis/06_Clinical/TCGA.markers.clin.corrplot.pdf', width = 6, height = 6, onefile = F)
col <- colorRampPalette(c(cor.colors[1:2], "white", cor.colors[3:4]))
corrplot(cor_M,
         order = 'original', 
         type="lower",
         addCoef.col = 'black', # 相关性系数文字颜色
         tl.col = "black",
         # tl.pos = 'd',
         cl.pos = 'n',
         diag=F,
         p.mat=cor_P,
         # insig="label_sig",
         col = col(200)[1:200])
dev.off()

######################### 与TME的关系
get_correlation=function(x,y,wider=T,longer=F){
  cor.res=psych::corr.test(x=x,y = y)
  df_cor=cor.res$r
  df_pval=cor.res$p
  dim(df_cor)
  ########
  library(tidyverse)
  g = pivot_longer(data=rownames_to_column(data.frame(df_cor,check.names = F),var = "from"),
                   cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                   names_to = "to",
                   values_to = "cor")
  gp = pivot_longer(data=rownames_to_column(data.frame(df_pval,check.names = F)),
                    cols = 2:(ncol(df_pval)+1),
                    names_to = "gene",
                    values_to = "pvalue")
  all(g$from==gp$rowname & g$to==gp$gene)
  head(g)
  g$pvalue=gp$pvalue
  g <- g %>%mutate(p.signif = cut(pvalue, breaks = c(0,0.0001, 0.001, 0.01, 0.05,1),
                                  labels = c("****", "***", "**","*",""),
                                  right = FALSE, include.lowest = TRUE))
  head(g)
  g$cor1 <- cut(abs(g$cor),# 绝对值
                breaks = c(0, 0.3, 0.5, 0.7, 0.9, 1),
                labels = c("< 0.3","0.3 - 0.5","0.5 - 0.7","0.7 - 0.9","> 0.9"),
                right=FALSE) # right=FALSE表示表示区间为左闭右开
  g$pvalue1 <- cut(g$pvalue,
                   breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                   labels = c("< 0.0001","< 0.001","< 0.01","< 0.05","> 0.05"),
                   right=FALSE) 
  
  cor.R <- pivot_wider(data=g[,c(1:2,3)],names_from ='to',values_from ='cor')
  cor.R=data.frame(cor.R,check.names = F)
  rownames(cor.R)=cor.R$from
  cor.R=cor.R[,-1]
  cor.R=as.matrix(cor.R)
  
  cor.Pval <- pivot_wider(data=g[,c(1:2,4)],names_from ='to',values_from ='pvalue')
  cor.Pval=data.frame(cor.Pval,check.names = F)
  rownames(cor.Pval)=cor.Pval$from
  cor.Pval=cor.Pval[,-1]
  cor.Pval=as.matrix(cor.Pval)
  all(rownames(cor.Pval)==rownames(cor.R))
  all(colnames(cor.Pval)==colnames(cor.R))
  if(wider){
    return(list(cor.R=cor.R,cor.Pval=cor.Pval))
  } else if(longer){
    return(g)
  }
}

plotMutiHeatmap1=function(up,down,up_pval,down_pval,up_break,down_break,up_colors,down_colors,row_split=NULL,column_split=NULL,title=''){
  library("ComplexHeatmap")
  library("circlize")
  UpColor <- circlize::colorRamp2(breaks = up_break, colors = up_colors) ## Up 相关性的颜色
  DnColor <- circlize::colorRamp2(breaks = down_break, colors = down_colors) ## Down 相关性的颜色
  
  
  DiagFunc <- function(up, down,up_pval,down_pval){
    function(j, i, x, y, width, height, fill){
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                   unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
      
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                   unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                   gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
      
      if(up_pval[i, j]<0.05){
        txt="****"
        if(up_pval[i, j]>=0.01 & up_pval[i, j]<0.05){
          txt='*'
        }else if(up_pval[i, j]>=0.001 & up_pval[i, j]<0.01){
          txt='**'
        }else if(up_pval[i, j]>=0.0001 & up_pval[i, j]<0.001){
          txt='***'
        }
        grid.text(label=txt,x=(x + 0.5*width),
                  y=(y+ 0.5*height),just = c('right','top'))
      }
      
      if(down_pval[i, j]<0.05){
        txt="****"
        if(down_pval[i, j]>=0.01 & down_pval[i, j]<0.05){
          txt='*'
        }else if(down_pval[i, j]>=0.001 & down_pval[i, j]<0.01){
          txt='**'
        }else if(down_pval[i, j]>=0.0001 & down_pval[i, j]<0.001){
          txt='***'
        }
        grid.text(label=txt,x=(x - 0.5*width),
                  y=(y - 0.5*height),just = c('left','bottom'))
      }
      
    }
  }
  
  p1 <- Heatmap(up, column_title = title
                , rect_gp = gpar(type = "none")
                , show_heatmap_legend = F
                , cluster_rows = T
                , cluster_columns = T
                , row_split=row_split
                , column_split=column_split
                , cell_fun = DiagFunc(up = up, down = down,up_pval = up_pval,down_pval = down_pval) 
  ) 
  p1
  ########### down legend
  col_fun = colorRamp2(down_break, down_colors) 
  lgd <- Legend(title = "Correlation", 
                col_fun = col_fun, 
                at = c(-1,0,1), 
                labels = c("-1","0","1"),  
                direction = "horizontal" 
  )
  ########### up legend
  col_fun2 = colorRamp2(up_break, up_colors) 
  lgd2 <- Legend(title = "-log10(p value)", 
                 col_fun = col_fun2, 
                 at = up_break, 
                 labels = c('0',"1","2","3","4",">5"),  
                 direction = "horizontal"
  )
  
  # draw(p1, annotation_legend_list = list(lgd,lgd2), annotation_legend_side = "bottom"
  #      ,heatmap_legend_side = "bottom", merge_legend = TRUE)
  draw(p1, annotation_legend_list = list(lgd), annotation_legend_side = "bottom"
       ,heatmap_legend_side = "bottom", merge_legend = TRUE)
}

############
get.IOBR.immu.format=function(tcga.t.exp.cibersort){
  tcga.t.exp.cibersort = data.frame(tcga.t.exp.cibersort)
  rownames(tcga.t.exp.cibersort) = tcga.t.exp.cibersort$ID
  tcga.t.exp.cibersort = tcga.t.exp.cibersort[, -1]
  colnames(tcga.t.exp.cibersort) = gsub('(.*)_.*', "\\1", colnames(tcga.t.exp.cibersort))
  return(tcga.t.exp.cibersort)
}
load("raw_datas/Pathway_ssGSEA/tcga.exp.cibersort.RData")
tcga.exp.cibersort=get.IOBR.immu.format(tcga.exp.cibersort)
tcga.t.cibersort=tcga.exp.cibersort[colnames(tcga.t.exp),1:22]
dim(tcga.t.cibersort) ### 513*22

################# 正常样本
dim(tcga.n.exp)
tcga.n.cibersort=tcga.exp.cibersort[colnames(tcga.n.exp),1:22]
dim(tcga.n.cibersort) ## 59*22

############
dt1=data.frame(t(tcga.t.exp[final.genes,]),check.names = F)
dt2=data.frame(tcga.t.cibersort,check.names = F)
dim(dt1) ## 513*2
dim(dt2) ## 513*22

alist=get_correlation(x = dt1,y = dt2,wider = T)
down=alist[[1]]
down_pval=alist[[2]]
down[which(is.na(down),arr.ind = T)]=0
down_pval[which(is.na(down_pval),arr.ind = T)]=1

########################### control
dt1=data.frame(t(tcga.n.exp[final.genes,]),check.names = F)
dt2=data.frame(tcga.n.cibersort,check.names = F)
dim(dt1) ## 59*2
dim(dt2) ## 59*22

alist=get_correlation(x = dt1,y = dt2,wider = T)
up=alist[[1]]
up_pval=alist[[2]]

########## correlation
up=up
down=down
########## correlation pvalue
up_pval=up_pval
down_pval=down_pval
########## correlation colors
cor.colors=pal_d3(alpha =1)(7)[c(1,3,2,4)]
scales::show_col(cor.colors)
cor.colors
up_break=c(-1,-0.5,0,0.5,1)
down_break=c(-1,-0.5,0,0.5,1)
up_colors=c(cor.colors[1:2],"white",cor.colors[3:4])
down_colors=c(cor.colors[1:2],"white",cor.colors[3:4])

########## group colors
# column_group = paths$Type[match(colnames(up), paths$Pathways)]
# row_group = c(rep("OSR",13),"ROS")
##########
pdf('analysis/07_TME/tcga.markers.cibersort.cor.heatmap.pdf',height = 8,width = 12)
plotMutiHeatmap1(up,down,up_pval,down_pval,up_break,down_break,up_colors,down_colors,column_split=NULL,row_split=NULL,title = '')
dev.off()

##############################
dt1=data.frame(t(tcga.t.exp[final.genes,]),check.names = F)
dt2=data.frame(tcga.t.cibersort,check.names = F)
dim(dt1) ## 513*2
dim(dt2) ## 513*22

df=get_correlation(x = dt1,y = dt2,wider = F,longer = T)
head(df)
max(abs(df$cor))

p.all=list()
for(i in final.genes){
  print(i)
  dat=df %>% dplyr::filter(from==i) %>% dplyr::arrange((cor))
  dat$to = factor(dat$to, levels = dat$to)
  dat=dat[!is.na(dat$cor),]
  dim(dat)
  print(range(dat$cor,na.rm = T))
  p = dat %>% ggplot(aes(x = cor, y = to, color = pvalue1)) +
    scale_color_manual(name="pvalue",values = c("#7F3C8D","#E73F74", "#F2B701", "#80BA5A", "gray"),breaks = c("< 0.0001","< 0.001","< 0.01","< 0.05","> 0.05"))+
    geom_segment(aes(x = 0, y = to, xend = cor, yend = to),size = 1) +
    geom_point(aes(size = cor1))+
    theme_test()+
    geom_vline(xintercept = c(0.0),size=0.25)+
    labs(size = "Cor")+
    labs(x = NULL, y = "")+
    theme(axis.line = element_line(size = 0.25),
          plot.margin =  unit(c(0.1,0.1,0.1,0.1),'cm'),
          axis.ticks = element_line(colour = "black",size = 0.25),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8,face = "plain",color = "black"), 
          legend.text = element_text(size = 6,face = "plain",color = "black"),
          legend.title = element_text(size = 6,face = "plain",color = "black"),
          legend.box.spacing = unit(1.2,'cm'))+
    coord_cartesian(clip = 'off',xlim = c(-0.35,0.35))+
    annotate(geom="text",x=0.4,y=dat$to,color="black", size=3,label=dat$pvalue1)+
    theme(aspect.ratio = 1.5,
          legend.position = c(0.95,0),
          legend.justification = c(0.95,0),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA))
  p=p+ggtitle(i)+labs(x='Correlation')
  p
  p.all=c(p.all,list(p))
}
length(p.all)
p.all[[1]]
plot1=mg_merge_plot(p.all,nrow = 1,ncol = 4,common.legend = T,labels = LETTERS[1:4],font.label = list(size = 14, color = "black", face ="bold" , family = 'Times'))
plot1

######### 与ESTIMATE软件的相关性
load("raw_datas/Pathway_ssGSEA/tcga.exp.estimate.RData")
tcga.exp.estimate=get.IOBR.immu.format(tcga.exp.estimate)
tcga.t.estimate=tcga.exp.estimate[colnames(tcga.t.exp),1:3]
dim(tcga.t.estimate)

df=data.frame(tcga.t.estimate,t(tcga.t.exp[final.genes,]))
head(df)

p.all=list()
for(i in c("StromalScore","ImmuneScore","ESTIMATEScore")){
  for(gene in final.genes){
    dat=df[,c(i,gene)]
    # dat=tidyr::pivot_longer(data = dt,cols=c(2:6),names_to="Gene",values_to ="Expression")
    # dat=data.frame(dat)
    head(dat)
    p <- ggscatter(dat, x = gene, y = i,
                   color =gene.colors[gene],
                   add = "reg.line", 
                   conf.int = TRUE)
    p=p + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),digits = 3)
    p
    p.all=c(p.all,list(p))
  }
}

length(p.all)
plot2=mg_merge_plot(p.all,nrow = 3,ncol = 4,labels = LETTERS[2:4],font.label = list(size = 14, color = "black", face ="bold" , family = 'Times'),common.legend = T)
plot2

figure=mg_merge_plot(plot1,plot2,nrow = 2,ncol = 1,heights = c(1,2))
figure
ggsave('analysis/07_TME/TCGA.marker.TME.corrplot.pdf',figure,width = 20,height = 20)

#################
dt1=data.frame(t(tcga.t.exp[final.genes,]),check.names = F)
dt2=data.frame(tcga.t.estimate,check.names = F)
dim(dt1) ## 513*5
dim(dt2) ## 513*3

alist=get_correlation(x = dt1,y = dt2,wider = T)
cor_mat=alist[[1]]
cor_pval=alist[[2]]

library(ggplot2)
library(ggradar)
# save(cor_mat,file = 'model.estimate.cor.RData')
# setwd("Z:/users/lishuang/Work1/Work/20231116_LUAD_Radiomics")
# load("model.estimate.cor.RData")
range(cor_mat)
dt=cor_mat %>% t() %>% data.frame() %>% add_rownames(var = "group")
dt

p=ggradar(dt,values.radar = c("-0.5", "-0.35", "-0.2"),
          grid.min = -0.5,
          grid.mid = -0.35,
          grid.max = -0.2,
          group.line.width = 1,
          group.point.size = 2,
          background.circle.colour = "white",
          gridline.mid.colour = "grey",
          base.size=12,
          legend.text.size=12,
          legend.position = "bottom")
p
ggsave('analysis/07_TME/TCGA.marker.ESTIMATE.corrplot.pdf',p,width = 4,height = 4)
save.image(file = '20231116_LUAD_Transcriptome_final.RData')

