### 路径

## 导入数据
load('raw_datas/GEO/gse72094.exp.RData')
load('raw_datas/GEO/gse31210.exp.RData')
load('raw_datas/GEO/gse30219.exp.RData')
load("raw_datas/GEO/gse50081.exp.RData")
load("raw_datas/GEO/gse37745.exp.RData")

######################### 去除批次效应
library("FactoMineR")
library("factoextra")
pca.plot = function(dat,col,title ="PCA - Biplot"){
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               addEllipses = TRUE,
               legend.title = "Groups",
               title=title)
}
###### PCG
comm.genes=Reduce(intersect,list(rownames(gse72094.exp),
                                 rownames(gse31210.exp),
                                 rownames(gse30219.exp),
                                 rownames(gse50081.exp)
                                 # rownames(gse37745.exp)
                                 ))
length(comm.genes) ## 18685
########
geo.batch = c(
  rep("GSE72094", ncol(gse72094.exp)),
  rep("GSE31210", ncol(gse31210.exp)),
  rep("GSE30219", ncol(gse30219.exp)),
  rep("GSE50081", ncol(gse50081.exp))
  # rep("GSE37745", ncol(gse50081.exp))
)
######## PCG
geo.exp = cbind(gse72094.exp[comm.genes, ],
                gse31210.exp[comm.genes, ],
                gse30219.exp[comm.genes, ],
                gse50081.exp[comm.genes, ]
                # gse37745.exp[comm.genes, ]
)

#################
geo.exp.all=geo.exp
dim(geo.exp.all) ## PCG
######### 导入临床信息数据
load("raw_datas/GEO/geo.t.cli.RData")
table(geo.t.cli$batch,geo.t.cli$Histology)
geo.t.cli=geo.t.cli[which(geo.t.cli$Histology=="LUAD"),]
dim(geo.t.cli) ## 940例LUAD
####### selected LUAD
dim(geo.exp.all)
geo.exp.all=geo.exp.all[,colnames(geo.exp.all) %in% rownames(geo.t.cli)]
dim(geo.exp.all) ## 834
geo.t.cli=geo.t.cli[colnames(geo.exp.all),]
dim(geo.t.cli)
######## 批次校正之前
pdf('raw_datas/GEO/geo.exprs.pca.pdf',height = 4,width = 4)
pca.plot(dat = geo.exp.all,col = geo.t.cli$batch)
dev.off()

######## 批次校正
library(sva)
geo.exp.combat=ComBat(dat = as.matrix(geo.exp.all)
                      , batch = geo.t.cli$batch
                      , par.prior = T)

##################
pdf('raw_datas/GEO/geo.exprs.adjusted.pca.pdf',height = 4,width = 4)
pca.plot(dat = geo.exp.combat,col = geo.t.cli$batch)
dev.off()

######### normalize.quantiles 标准化
geo.exp.norm=preprocessCore::normalize.quantiles(geo.exp.combat)
rownames(geo.exp.norm)=rownames(geo.exp.combat)
colnames(geo.exp.norm)=colnames(geo.exp.combat)
range(geo.exp.norm)
##################
geo.exp.pcg=geo.exp.norm[comm.genes,]
save(geo.exp.pcg,file = 'raw_datas/GEO/geo.exp.pcg.RData')
################ END