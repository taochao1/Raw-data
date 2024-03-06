
gencode.pcg$gene_id=gsub("\\..*","",gencode.pcg$gene_id)

########### GDC 下载数据处理
### TCGA-临床信息
luad.tcga.cli=read.table('raw_datas/TCGA/TCGA-LUAD_clinical.txt',sep = "\t",header = T,stringsAsFactors = F,check.names = F)
luad.tcga.cli$project_id="TCGA-LUAD"
####
tcga.cli=luad.tcga.cli
###########
tcga.cli$Tumor_Sample_Barcode=tcga.cli$A0_Samples
tcga.cli$os=tcga.cli$A1_OS
tcga.cli$os_status=tcga.cli$A2_Event
dim(tcga.cli)
tcga.cli=tcga.cli[,c(179:ncol(tcga.cli),1:178)]
names(table(tcga.cli$os_status))

tcga.cli=data.frame(SampleID=paste0(tcga.cli$Tumor_Sample_Barcode,"-01"),tcga.cli,stringsAsFactors = F,check.names = F)
rownames(tcga.cli)=tcga.cli$SampleID

table(tcga.cli$os_status)
tcga.cli$os_status[tcga.cli$os_status=='Alive']=0
tcga.cli$os_status[tcga.cli$os_status=='Dead']=1
tcga.cli$os_status[tcga.cli$os_status=='']=NA
table(tcga.cli$os_status)
#################
### TCGA-表达谱
###TCGA-LUAD FPKM数据
# luad.tcga.exp.fpkm=read.table("raw_datas/TCGA/TCGA-LUAD_FPKM.txt"
#                          ,row.names = 1,sep="\t",header=T,as.is=T,quote="\""
#                          ,fill=T,check.names = F,stringsAsFactors = F)
# dim(luad.tcga.exp.fpkm)##60483*574
# head(rownames(luad.tcga.exp.fpkm))
# luad.tcga.exp.tpm=mg_FPKM2TPMs(luad.tcga.exp.fpkm)
# dim(luad.tcga.exp.tpm)##60483*574
# luad.tcga.exp.log=log2(luad.tcga.exp.tpm+1)
# luad.tcga.exp=exp_ensg2symbol(luad.tcga.exp.log)
# dim(luad.tcga.exp)## 25554*574
# #######保存
# save(luad.tcga.exp,file='raw_datas/TCGA/luad.tcga.exp.RData')
load("raw_datas/TCGA/luad.tcga.exp.RData")
luad.tcga.exp=luad.tcga.exp[rownames(luad.tcga.exp) %in% gencode.pcg$gene_name,]
dim(luad.tcga.exp) ## 19037 * 574
range(luad.tcga.exp)
## 01: 513
## 02: 2
## 11:59
table(substr(colnames(luad.tcga.exp),14,15))
####################
tcga.t.exp=luad.tcga.exp[,substr(colnames(luad.tcga.exp),14,15)=="01"]
tcga.n.exp=luad.tcga.exp[,substr(colnames(luad.tcga.exp),14,15)=="11"]
dim(tcga.t.exp) ## 516
dim(tcga.n.exp) ## 59

#######
tcga.t.cli=tcga.cli[na.omit(match(colnames(tcga.t.exp),tcga.cli$SampleID)),]
tcga.t.cli=dplyr::rename(tcga.t.cli,c('OS.time'='A1_OS'
                                      , 'OS' = 'A2_Event'))
tcga.t.cli$OS[which(tcga.t.cli$OS=='Alive')]=0
tcga.t.cli$OS[which(tcga.t.cli$OS=='Dead')]=1
tcga.t.cli$OS[which(tcga.t.cli$OS=='')]=NA
table(tcga.t.cli$OS)
tcga.t.cli=crbind2DataFrame(tcga.t.cli)

table(tcga.t.cli$OS.time>=30,tcga.t.cli$OS)
table(tcga.t.cli$OS.time<365*15,tcga.t.cli$OS)
# tcga.t.cli=tcga.t.cli[which(tcga.t.cli$os>=30),]
####
table(tcga.t.cli$project_id) ## 490 例LUAD
tcga.t.cli=tcga.t.cli[tcga.t.cli$project_id=='TCGA-LUAD',]
tcga.t.exp=luad.tcga.exp[,tcga.t.cli$SampleID]
dim(tcga.t.exp) ## 19037*490
########## 临床信息整理
clin.selected=c('A17_Age','A18_Sex','A3_T','A4_N','A5_M','A6_Stage')

setdiff(clin.selected,colnames(tcga.t.cli))

luad.tcga.t.cli=tcga.t.cli[,clin.selected]

luad.tcga.t.cli.tnm=clean_TNMStage(st=luad.tcga.t.cli$A3_T
                                   , sn = luad.tcga.t.cli$A4_N
                                   , sm = luad.tcga.t.cli$A5_M
                                   , ss = luad.tcga.t.cli$A6_Stage
                                   , sex = luad.tcga.t.cli$A18_Sex
                                   , age = luad.tcga.t.cli$A17_Age
                                   , age_cut = NULL)
row.names(luad.tcga.t.cli.tnm)=rownames(tcga.t.cli)
dim(luad.tcga.t.cli.tnm)

colnames(luad.tcga.t.cli.tnm)

#############
tcga.t.cli_use=cbind(tcga.t.cli[,c('OS.time','OS','A8_New_Event_Time','A8_New_Event','tobacco_smoking_history')],
                     luad.tcga.t.cli.tnm)
colnames(tcga.t.cli_use)
median(tcga.t.cli_use$Age,na.rm = T)

tcga.t.cli_use$Age1=ifelse(tcga.t.cli_use$Age>60,'>60','<=60')
tcga.t.cli_use$Status=ifelse(tcga.t.cli_use$OS==0,'Alive','Dead')

tcga.t.cli_use=dplyr::rename(tcga.t.cli_use,c('T Stage'='Clinical_T'
                                              ,'N Stage'='Clinical_N'
                                              ,'M Stage'='Clinical_M'
                                              ,'Stage'='Clinical_Stage'
                                              ,'PFI.time'='A8_New_Event_Time'
                                              ,'PFI'='A8_New_Event'
                                              ,'Smoker'='tobacco_smoking_history'
))
dim(tcga.t.cli_use)
tcga.t.cli_use$Stage=gsub("Stage ","",tcga.t.cli_use$Stage)
tcga.t.cli_use$Smoker[which(tcga.t.cli_use$Smoker==1)]='Never'
tcga.t.cli_use$Smoker[which(tcga.t.cli_use$Smoker==2)]='Current'
tcga.t.cli_use$Smoker[which(tcga.t.cli_use$Smoker==3)]='Ever'
tcga.t.cli_use$Smoker[which(tcga.t.cli_use$Smoker==4)]='Ever'
tcga.t.cli_use$Smoker[which(tcga.t.cli_use$Smoker==5)]='Ever'

save(tcga.t.cli_use,file = 'raw_datas/TCGA/tcga.t.cli_use.RData')
save(tcga.t.exp,file = 'raw_datas/TCGA/tcga.t.exp.RData')
