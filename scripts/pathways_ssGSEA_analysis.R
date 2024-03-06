setwd('Z:/users/lishuang/Work1/Work/20231116_LUAD_Radiomics/')
source('Z:/projects/codes/mg_base.R')
load("raw_datas/TCGA/luad.tcga.exp.RData")
# load("raw_datas/GEO/gse31210.exp.RData")
# load("raw_datas/GEO/gse72094.exp.RData")
# load("raw_datas/GEO/gse50081.exp.RData")
# load("raw_datas/GEO/geo.exp.pcg.RData")

####
h <- clusterProfiler::read.gmt("raw_datas/h.all.v2023.1.Hs.symbols.gmt")
geneSets<-lapply(unique(h$term),function(x){h$gene[h$term==x]})
names(geneSets) <- unique(h$term)

########## TCGA
tcga.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(luad.tcga.exp)
                                                , genelist = geneSets)
save(tcga.h.all.ssgsea,file = 'raw_datas/Pathway_ssGSEA/tcga.h.all.ssgsea.RData')
# ########## GSE31210
# gse31210.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse31210.exp)
#                                                     , genelist = geneSets)
# save(gse31210.h.all.ssgsea,file = 'raw_datas/Pathway_ssGSEA/gse31210.h.all.ssgsea.RData')
# ########## GSE72094
# gse72094.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse72094.exp)
#                                                       , genelist = geneSets)
# save(gse72094.h.all.ssgsea,file = 'raw_datas/Pathway_ssGSEA/gse72094.h.all.ssgsea.RData')
# ########## GSE50081
# gse50081.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse50081.exp)
#                                                       , genelist = geneSets)
# save(gse50081.h.all.ssgsea,file = 'raw_datas/Pathway_ssGSEA/gse50081.h.all.ssgsea.RData')
# 
# ########## GEO
# geo.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(geo.exp.pcg)
#                                                       , genelist = geneSets)
# save(geo.h.all.ssgsea,file = 'raw_datas/Pathway_ssGSEA/geo.h.all.ssgsea.RData')
# 
# ########### metabolism pathway
# load("raw_datas/metabolism.paths.genelist.RData")
# length(paths.genelist)
# tcga.metabolism.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(luad.tcga.exp)
#                                                   , genelist = paths.genelist)
# save(tcga.metabolism.ssgsea,file = 'raw_datas/Pathway_ssGSEA/tcga.metabolism.ssgsea.RData')
# 
# geo.metabolism.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(geo.exp.pcg)
#                                                  , genelist = paths.genelist)
# save(geo.metabolism.ssgsea,file = 'raw_datas/Pathway_ssGSEA/geo.metabolism.ssgsea.RData')


##################
library(IOBR)
#### ESTIMATE
tcga.exp.estimate<-deconvo_estimate(eset=luad.tcga.exp)
save(tcga.exp.estimate,file='raw_datas/Pathway_ssGSEA/tcga.exp.estimate.RData')
### CIBERSORT
tcga.exp.cibersort<-deconvo_cibersort(eset=luad.tcga.exp,arrays=T)
save(tcga.exp.cibersort,file='raw_datas/Pathway_ssGSEA/tcga.exp.cibersort.RData')
############ TIMER 
tcga.exp.timer<-deconvo_timer(eset=as.matrix(luad.tcga.exp),indications=rep('hnsc',ncol(luad.tcga.exp)))
save(tcga.exp.timer,file='raw_datas/Pathway_ssGSEA/tcga.exp.timer.RData')
############ EPIC 
tcga.exp.epic<-deconvo_epic(eset=as.matrix(luad.tcga.exp),tumor = TRUE)
save(tcga.exp.epic,file='raw_datas/Pathway_ssGSEA/tcga.exp.epic.RData')
############ MCP-counter 
tcga.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(luad.tcga.exp))
save(tcga.exp.mcp,file='raw_datas/Pathway_ssGSEA/tcga.exp.mcp.RData')

