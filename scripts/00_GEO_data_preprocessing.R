##################### GEO 数据


head(gencode.pcg)
## GSE72094 dataset
# GSE72094=getGEOExpDataByCel('GSE72094')
# save(GSE72094,file='raw_datas/GEO/GSE72094.RData')
load('raw_datas/GEO/GSE72094.RData')
########## 基因表达谱整理
# range(GSE72094$Exp$GPL15048_Data)
# gse72094.exp<-exp_probe2symbol_v2(GSE72094$Exp$GPL15048_Data,GPL = 'GPL15048')
# save(gse72094.exp,file='raw_datas/GEO/gse72094.exp.RData')
load('raw_datas/GEO/gse72094.exp.RData')
gse72094.exp=gse72094.exp[rownames(gse72094.exp) %in% gencode.pcg$gene_name,]

########## 临床信息整理
gse72094.cli<-GSE72094$Sample
rownames(gse72094.cli)=gse72094.cli$Acc
table(gse72094.cli$Source) ## 442例肺腺癌

######### TNM
gse72094.cli.tnm<-clean_TNMStage(ss = gse72094.cli$Stage
                                 , sex = gse72094.cli$gender
                                 , age = gse72094.cli$age_at_diagnosis
                                 , age_cut = NULL)
rownames(gse72094.cli.tnm)<-gse72094.cli$Acc

gse72094.cli.os<-data.frame(OS=gse72094.cli$vital_status
                            ,OS.time=gse72094.cli$survival_time_in_days)
rownames(gse72094.cli.os)<-gse72094.cli$Acc
gse72094.cli.os<-crbind2DataFrame(gse72094.cli.os)

table(gse72094.cli.os$OS)
gse72094.cli.os$OS[gse72094.cli.os$OS=='Dead']<-1
gse72094.cli.os$OS[gse72094.cli.os$OS=='Alive']<-0
gse72094.cli.os$OS[gse72094.cli.os$OS=='NA']<-NA
gse72094.cli.os<-crbind2DataFrame(gse72094.cli.os)

#############
gse72094.cli_df=data.frame(gse72094.cli.os,gse72094.cli.tnm
                           ,Smoking_status=gse72094.cli$smoking_status,stringsAsFactors = F,check.names = F)
head(gse72094.cli_df)
table(gse72094.cli_df$Smoking_status)
gse72094.cli_df$Smoking_status[which(gse72094.cli_df$Smoking_status=="Missing")]=NA
gse72094.cli_df$Histology="LUAD"

############################
gse72094.t.cli=gse72094.cli_df[!is.na(gse72094.cli_df$OS)&!is.na(gse72094.cli_df$OS.time),]
######### OS.time>0 样本保留
gse72094.t.cli<-gse72094.t.cli[which(gse72094.t.cli$OS.time>0),]
dim(gse72094.t.cli)
table(gse72094.t.cli$OS.time>30)
table(gse72094.t.cli$OS.time<365*15)

########
gse72094.t.exp=gse72094.exp[,rownames(gse72094.t.cli)]
dim(gse72094.t.exp)
range(gse72094.t.exp)
########################################
## GSE31210 dataset
# GSE31210=getGEOExpDataByCel('GSE31210')
# save(GSE31210,file='raw_datas/GEO/GSE31210.RData')
load('raw_datas/GEO/GSE31210.RData')
########## 基因表达谱整理
# gse31210.exp<-exp_probe2symbol_v2(GSE31210$Exp$GPL570_Data,GPL='GPL570')
# save(gse31210.exp,file='raw_datas/GEO/gse31210.exp.RData')
load('raw_datas/GEO/gse31210.exp.RData')
gse31210.exp=gse31210.exp[rownames(gse31210.exp) %in% gencode.pcg$gene_name,]

########## 临床信息整理
gse31210.cli<-GSE31210$Sample
rownames(gse31210.cli)=gse31210.cli$Acc
###########
table(gse31210.cli$tissue)##226例肺腺癌+20例癌旁组织样本
gse31210.cli=gse31210.cli[which(gse31210.cli$tissue=="primary lung tumor"),]
gse31210.exp=gse31210.exp[,rownames(gse31210.cli)]
dim(gse31210.exp)
#############
gse31210.cli.tnm<-clean_TNMStage(ss = gse31210.cli$`pathological stage`
                                 , sex = gse31210.cli$gender
                                 , age = gse31210.cli$`age (years)`
                                 , age_cut = NULL)
rownames(gse31210.cli.tnm)<-gse31210.cli$Acc
head(gse31210.cli.tnm)
######### 生存时间
gse31210.cli.os<-data.frame(OS=gse31210.cli$death
                            ,OS.time=gse31210.cli$`days before death/censor`
                            ,PFI=gse31210.cli$relapse
                            ,PFI.time=gse31210.cli$`days before relapse/censor`)
rownames(gse31210.cli.os)<-gse31210.cli$Acc

gse31210.cli.os<-crbind2DataFrame(gse31210.cli.os)
table(gse31210.cli.os$OS)
gse31210.cli.os$OS[gse31210.cli.os$OS=='dead']<-1
gse31210.cli.os$OS[gse31210.cli.os$OS=='alive']<-0

table(gse31210.cli.os$PFI)
gse31210.cli.os$PFI[gse31210.cli.os$PFI=='relapsed']<-1
gse31210.cli.os$PFI[gse31210.cli.os$PFI=='not relapsed']<-0

gse31210.cli.os<-crbind2DataFrame(gse31210.cli.os)
##########
gse31210.cli_df=data.frame(gse31210.cli.os,gse31210.cli.tnm
                           ,Smoking_status=gse31210.cli$`smoking status`,stringsAsFactors = F,check.names = F)
table(gse31210.cli_df$Smoking_status)
gse31210.cli_df$Smoking_status=gsub("(.*)\\-.*","\\1",gse31210.cli_df$Smoking_status)
gse31210.cli_df$Histology="LUAD"

######### 样本筛选
gse31210.t.cli=gse31210.cli_df[!is.na(gse31210.cli_df$OS)&!is.na(gse31210.cli_df$OS.time),]
######### OS.time>=30 样本保留
gse31210.t.cli<-gse31210.t.cli[which(gse31210.t.cli$OS.time>0),]
dim(gse31210.t.cli)
table(gse31210.t.cli$OS.time>30)
table(gse31210.t.cli$OS.time<365*15)

########
gse31210.t.exp=gse31210.exp[,rownames(gse31210.t.cli)]
dim(gse31210.t.exp)
range(gse31210.t.exp)

########################################
## GSE30219 dataset
# GSE30219=getGEOExpDataByCel('GSE30219')
# save(GSE30219,file='raw_datas/GEO/GSE30219.RData')
load('raw_datas/GEO/GSE30219.RData')
########## 基因表达谱整理
# gse30219.exp<-exp_probe2symbol_v2(GSE30219$Exp$GPL570_Data,GPL='GPL570')
# save(gse30219.exp,file='raw_datas/GEO/gse30219.exp.RData')
load('raw_datas/GEO/gse30219.exp.RData')
gse30219.exp=gse30219.exp[rownames(gse30219.exp) %in% gencode.pcg$gene_name,]

########## 临床信息整理
gse30219.cli<-GSE30219$Sample
rownames(gse30219.cli)=gse30219.cli$Acc
gse30219.cli$Title=gsub("(.*?)\\d+","\\1",gse30219.cli$Title,perl = TRUE)
gse30219.cli$Title=gsub("(.*?) ADK$","\\1 ADC",gse30219.cli$Title)
table(gse30219.cli$Title,gse30219.cli$Source)
dim(gse30219.cli)
table(gse30219.cli$tissue,gse30219.cli$histology)
table(gse30219.cli$Source,gse30219.cli$histology)

##307例样本，其中包含了14例非肿瘤肺和293例肺癌，在肺癌样本中包含了21例小细胞肺癌
###### 去掉Non tumoral lung
gse30219.cli<-gse30219.cli[-which(gse30219.cli$histology=='NTL'),]
dim(gse30219.cli)
######
## Lung cancer ADC: ADC
## Lung cancer Basaloid: BAS【肺的基底细胞样癌】
## Lung cancer Carcinoid: CARCI
## Lung cancer Large Cell Carcinoma: LCC
## Lung cancer Large Cell Neuroendocrine: LCNE
## Lung cancer_Other: Other
## Lung Squamous Cell Carcinoma: SQC
## Small Cell Lung Carcinoma: SCC
df=gse30219.cli[,c("Title","histology")]
df=df[!duplicated(df),]
df
###### 小细胞肺癌：SCLC是一种高级别神经内分泌癌，与其他神经内分泌起源的肿瘤归为一类
gse30219.cli<-gse30219.cli[-which(gse30219.cli$histology=='SCC'),]
dim(gse30219.cli)

###############
gse30219.exp=gse30219.exp[,rownames(gse30219.cli)]
dim(gse30219.exp)##272例NSCLC

###########
table(gse30219.cli$Title,gse30219.cli$tissue)
gse30219.cli.tnm<-clean_TNMStage(st = gse30219.cli$`pt stage`
                                 , sn = gse30219.cli$`pn stage`
                                 , sm=gse30219.cli$`pm stage`
                                 , sex = gse30219.cli$gender
                                 , age = gse30219.cli$`age at surgery`
                                 , age_cut = NULL)
rownames(gse30219.cli.tnm)<-gse30219.cli$Acc
head(gse30219.cli.tnm)

table(gse30219.cli.tnm$Clinical_T,gse30219.cli.tnm$Clinical_M)

gse30219.cli.tnm$Clinical_Stage[which(gse30219.cli.tnm$Clinical_N=="N0" & (gse30219.cli.tnm$Clinical_T %in% c('T1')) & (gse30219.cli.tnm$Clinical_M=="M0"))]='Stage I'
gse30219.cli.tnm$Clinical_Stage[which(gse30219.cli.tnm$Clinical_N=="N0" & (gse30219.cli.tnm$Clinical_T %in% c('T2')) & (gse30219.cli.tnm$Clinical_M=="M0"))]='Stage II'
gse30219.cli.tnm$Clinical_Stage[which(gse30219.cli.tnm$Clinical_N=="N0" & (gse30219.cli.tnm$Clinical_T %in% c('T3')) & (gse30219.cli.tnm$Clinical_M=="M0"))]='Stage II'
gse30219.cli.tnm$Clinical_Stage[which(gse30219.cli.tnm$Clinical_N=="N0" & (gse30219.cli.tnm$Clinical_T %in% c('T4')) & (gse30219.cli.tnm$Clinical_M=="M0"))]='Stage III'
gse30219.cli.tnm$Clinical_Stage[which(gse30219.cli.tnm$Clinical_N=="N1" & (gse30219.cli.tnm$Clinical_T %in% c('T1','T2') & (gse30219.cli.tnm$Clinical_M=="M0")))]='Stage II'
gse30219.cli.tnm$Clinical_Stage[which(gse30219.cli.tnm$Clinical_N=="N1" & (gse30219.cli.tnm$Clinical_T %in% c('T3','T4')) & (gse30219.cli.tnm$Clinical_M=="M0"))]='Stage III'
gse30219.cli.tnm$Clinical_Stage[which(gse30219.cli.tnm$Clinical_N=="N2" & (gse30219.cli.tnm$Clinical_M=="M0"))]='Stage III'
gse30219.cli.tnm$Clinical_Stage[which(gse30219.cli.tnm$Clinical_N=="N3" & (gse30219.cli.tnm$Clinical_M=="M0"))]='Stage III'
gse30219.cli.tnm$Clinical_Stage[which((gse30219.cli.tnm$Clinical_M=="M1"))]='Stage IV'
table(gse30219.cli.tnm$Clinical_Stage)

######### 生存时间
gse30219.cli.os<-data.frame(OS=gse30219.cli$status
                            ,OS.time=gse30219.cli$`follow-up time (months)`
                            ,PFI=gse30219.cli$`relapse (event=1; no event=0)`
                            ,PFI.time=gse30219.cli$`disease free survival in months`)
rownames(gse30219.cli.os)<-gse30219.cli$Acc

gse30219.cli.os<-crbind2DataFrame(gse30219.cli.os)
table(gse30219.cli.os$OS)
gse30219.cli.os$OS[which(gse30219.cli.os$OS=='DEAD')]<-1
gse30219.cli.os$OS[which(gse30219.cli.os$OS=='ALIVE')]<-0
table(gse30219.cli.os$PFI)

gse30219.cli.os$OS.time=gse30219.cli.os$OS.time*30
gse30219.cli.os$PFI.time=gse30219.cli.os$PFI.time*30
gse30219.cli.os<-crbind2DataFrame(gse30219.cli.os)
##########
gse30219.cli_df=data.frame(gse30219.cli.os,gse30219.cli.tnm)
gse30219.cli_df$Histology=gse30219.cli$histology
table(gse30219.cli_df$Histology)
dim(gse30219.cli_df)

gse30219.cli_df$Histology[gse30219.cli_df$Histology %in% "ADC"]="LUAD"
gse30219.cli_df$Histology[gse30219.cli_df$Histology %in% "SQC"]="LUSC"
gse30219.cli_df$Histology[gse30219.cli_df$Histology %in% c("BAS","CARCI","Other")]="Other"

######### 样本筛选
gse30219.t.cli=gse30219.cli_df[!is.na(gse30219.cli_df$OS)&!is.na(gse30219.cli_df$OS.time),]
table(gse30219.t.cli$Histology,gse30219.t.cli$OS.time>=30)

######### OS.time>0 样本保留
gse30219.t.cli<-gse30219.t.cli[which(gse30219.t.cli$OS.time>0),]
dim(gse30219.t.cli)
table(gse30219.t.cli$OS.time>30)
table(gse30219.t.cli$OS.time<365*15)

########
gse30219.t.exp=gse30219.exp[,rownames(gse30219.t.cli)]
dim(gse30219.t.exp)
range(gse30219.t.exp)

########################################
## GSE50081 dataset
# GSE50081=getGEOExpDataByCel('GSE50081')
# save(GSE50081,file='raw_datas/GEO/GSE50081.RData')
load('raw_datas/GEO/GSE50081.RData')
########## 基因表达谱整理
# gse50081.exp<-exp_probe2symbol_v2(GSE50081$Exp$GPL570_Data,GPL='GPL570')
# save(gse50081.exp,file='raw_datas/GEO/gse50081.exp.RData')
load('raw_datas/GEO/gse50081.exp.RData')
gse50081.exp=gse50081.exp[rownames(gse50081.exp) %in% gencode.pcg$gene_name,]

########## 临床信息整理
gse50081.cli<-GSE50081$Sample
rownames(gse50081.cli)=gse50081.cli$Acc
###########
gse50081.cli<-gse50081.cli[which(gse50081.cli$Source=='Primary NSCLC'),]
gse50081.cli.tnm<-clean_TNMStage(ss = gse50081.cli$Stage
                                 ,st=gse50081.cli$`t-stage`
                                 ,sn=gse50081.cli$`n-stage`
                                 ,sm=gse50081.cli$`m-stage`
                                 , sex = gse50081.cli$Sex
                                 , age = gse50081.cli$age
                                 , age_cut = NULL)
rownames(gse50081.cli.tnm)<-gse50081.cli$Acc
head(gse50081.cli.tnm)
######### 生存时间
gse50081.cli.os<-data.frame(OS=gse50081.cli$status
                            ,OS.time=gse50081.cli$`survival time`
                            ,PFI=gse50081.cli$recurrence
                            ,PFI.time=gse50081.cli$`disease-free survival time`)
rownames(gse50081.cli.os)<-gse50081.cli$Acc
gse50081.cli.os<-crbind2DataFrame(gse50081.cli.os)

table(gse50081.cli.os$OS)
gse50081.cli.os$OS[which(gse50081.cli.os$OS=='dead')]<-1
gse50081.cli.os$OS[which(gse50081.cli.os$OS=='alive')]<-0

table(gse50081.cli.os$PFI)
gse50081.cli.os$PFI[which(gse50081.cli.os$PFI=='Y')]<-1
gse50081.cli.os$PFI[which(gse50081.cli.os$PFI=='N')]<-0
gse50081.cli.os$PFI[which(gse50081.cli.os$PFI=='U')]<-NA

gse50081.cli.os<-crbind2DataFrame(gse50081.cli.os)
gse50081.cli.os$OS.time=gse50081.cli.os$OS.time*365
gse50081.cli.os$PFI.time=gse50081.cli.os$PFI.time*365

gse50081.cli_df=data.frame(gse50081.cli.os,gse50081.cli.tnm
                           ,Smoking_status=gse50081.cli$smoking,stringsAsFactors = F,check.names = F)
head(gse50081.cli_df)
table(gse50081.cli_df$Smoking_status)
gse50081.cli_df$Smoking_status[which(gse50081.cli_df$Smoking_status=="Unable to determine")]=NA
gse50081.cli_df$Smoking_status[which(gse50081.cli_df$Smoking_status=="Ex-smoker")]='Ever'
gse50081.cli_df$Histology=gse50081.cli$histology
table(gse50081.cli_df$Histology)

gse50081.cli_df$Histology[gse50081.cli_df$Histology %in% c("adenocarcinoma")]="LUAD"
gse50081.cli_df$Histology[gse50081.cli_df$Histology %in% c("squamous cell carcinoma","squamous cell carcinoma X2")]="LUSC"
gse50081.cli_df$Histology[gse50081.cli_df$Histology %in% c("large cell carcinoma")]="LCC"
## Include large cell and adeno-squamous carcinoma.
gse50081.cli_df$Histology[gse50081.cli_df$Histology %in% c("adenosquamous carcinoma","NSClarge cell carcinoma-mixed","NSCLC-favor adenocarcinoma")]="Other"
table(gse50081.cli_df$Histology)

######### 样本筛选
gse50081.t.cli=gse50081.cli_df[!is.na(gse50081.cli_df$OS)&!is.na(gse50081.cli_df$OS.time),]
######### OS.time>0 样本保留
gse50081.t.cli<-gse50081.t.cli[which(gse50081.t.cli$OS.time>0),]
dim(gse50081.t.cli)
########
gse50081.t.exp=gse50081.exp[,rownames(gse50081.t.cli)]
dim(gse50081.t.exp)
range(gse50081.t.exp)

## GSE37745 dataset
# GSE37745=getGEOExpDataByCel('GSE37745')
# save(GSE37745,file='raw_datas/GEO/GSE37745.RData')
load('raw_datas/GEO/GSE37745.RData')
########## 基因表达谱整理
# gse37745.exp<-exp_probe2symbol_v2(GSE37745$Exp$GPL570_Data,GPL='GPL570')
# save(gse37745.exp,file='raw_datas/GEO/gse37745.exp.RData')
load('raw_datas/GEO/gse37745.exp.RData')
gse37745.exp=gse37745.exp[rownames(gse37745.exp) %in% gencode.pcg$gene_name,]

########## 临床信息整理
gse37745.cli<-GSE37745$Sample
rownames(gse37745.cli)=gse37745.cli$Acc
###########
table(gsub("(.*)patient.*","\\1",gse37745.cli$Source))

gse37745.cli.tnm<-clean_TNMStage(ss = gse37745.cli$`tumor stage`
                                 , sex = gse37745.cli$gender
                                 , age = gse37745.cli$age
                                 , age_cut = NULL)
rownames(gse37745.cli.tnm)<-gse37745.cli$Acc
head(gse37745.cli.tnm)
table(gse37745.cli.tnm$Clinical_Stage)
table(gse37745.cli$`tumor stage`)

######### 生存时间
gse37745.cli.os<-data.frame(OS=gse37745.cli$dead
                            ,OS.time=gse37745.cli$`days to determined death status`
                            ,PFI=gse37745.cli$recurrence
                            ,PFI.time=gse37745.cli$`days to recurrence / to last visit`)
rownames(gse37745.cli.os)<-gse37745.cli$Acc
gse37745.cli.os<-crbind2DataFrame(gse37745.cli.os)

table(gse37745.cli.os$OS)
gse37745.cli.os$OS[which(gse37745.cli.os$OS=='yes')]<-1
gse37745.cli.os$OS[which(gse37745.cli.os$OS=='no')]<-0

table(gse37745.cli.os$PFI)
gse37745.cli.os$PFI[which(gse37745.cli.os$PFI=='yes')]<-1
gse37745.cli.os$PFI[which(gse37745.cli.os$PFI=='no')]<-0
gse37745.cli.os$PFI[which(gse37745.cli.os$PFI=='not known')]<-NA
gse37745.cli.os$PFI.time[which(gse37745.cli.os$PFI.time=="not known")]=NA
gse37745.cli.os<-crbind2DataFrame(gse37745.cli.os)

gse37745.cli_df=data.frame(gse37745.cli.os,gse37745.cli.tnm
                           ,stringsAsFactors = F,check.names = F)
head(gse37745.cli_df)
gse37745.cli_df$Histology=gse37745.cli$histology
table(gse37745.cli_df$Histology,gse37745.cli_df$OS.time>30)

## Adenocarcinoma:LUAD
## Squamous cell ca:LUSC
## Large cell ca:LCC
gse37745.cli_df$Histology[which(gse37745.cli_df$Histology=="adeno")]="LUAD"
gse37745.cli_df$Histology[which(gse37745.cli_df$Histology=="squamous")]="LUSC"
gse37745.cli_df$Histology[which(gse37745.cli_df$Histology=="large")]="LCC"
table(gse37745.cli_df$Histology,gse37745.cli_df$OS.time>30)

######### 样本筛选
gse37745.t.cli=gse37745.cli_df[!is.na(gse37745.cli_df$OS)&!is.na(gse37745.cli_df$OS.time),]
######### OS.time>0 样本保留
gse37745.t.cli<-gse37745.t.cli[which(gse37745.t.cli$OS.time>0),]
dim(gse37745.t.cli)
########
gse37745.t.exp=gse37745.exp[,rownames(gse37745.t.cli)]
dim(gse37745.t.exp)
range(gse37745.t.exp)

################################ 临床信息的整理
table(gse72094.t.cli$Histology)## LUAD：398
table(gse31210.t.cli$Histology)## LUAD：226
table(gse50081.t.cli$Histology)## LUAD,LUSC,LCC,Other
# LCC  LUAD  LUSC Other 
# 7   127    43     4 
table(gse30219.t.cli$Histology)## LUAD,LUSC,LCC,LCNE,Other
# LCC  LCNE  LUAD  LUSC Other 
# 3    56    83    61    65 
table(gse37745.t.cli$Histology)## LUAD,LUSC,LCC
# LCC LUAD LUSC 
# 24  106   66 

####################
gse72094.t.cli=dplyr::rename(gse72094.t.cli,c('Stage'='Clinical_Stage','Smoker'='Smoking_status'))
gse31210.t.cli=dplyr::rename(gse31210.t.cli,c('Stage'='Clinical_Stage','Smoker'='Smoking_status'))
gse50081.t.cli=dplyr::rename(gse50081.t.cli,c('T Stage'='Clinical_T'
                                              ,'N Stage'='Clinical_N'
                                              ,'M Stage'='Clinical_M'
                                              ,'Stage'='Clinical_Stage'
                                              ,'Smoker'='Smoking_status'))
gse30219.t.cli=dplyr::rename(gse30219.t.cli,c('T Stage'='Clinical_T'
                                              ,'N Stage'='Clinical_N'
                                              ,'M Stage'='Clinical_M'
                                              ,'Stage'='Clinical_Stage'))
gse37745.t.cli=dplyr::rename(gse37745.t.cli,c('Stage'='Clinical_Stage'))

gse72094.t.cli$Stage=gsub("Stage ","",gse72094.t.cli$Stage)
gse31210.t.cli$Stage=gsub("Stage ","",gse31210.t.cli$Stage)
gse50081.t.cli$Stage=gsub("Stage ","",gse50081.t.cli$Stage)
gse30219.t.cli$Stage=gsub("Stage ","",gse30219.t.cli$Stage)
gse37745.t.cli$Stage=gsub("Stage ","",gse37745.t.cli$Stage)
#############
geo.t.cli=dplyr::bind_rows(gse72094.t.cli,gse31210.t.cli,
                           gse50081.t.cli,gse30219.t.cli,
                           gse37745.t.cli)
head(geo.t.cli)

geo.t.cli$batch=c(rep("GSE72094",nrow(gse72094.t.cli)),
                  rep("GSE31210",nrow(gse31210.t.cli)),
                  rep("GSE50081",nrow(gse50081.t.cli)),
                  rep("GSE30219",nrow(gse30219.t.cli)),
                  rep("GSE37745",nrow(gse37745.t.cli)))

table(geo.t.cli$batch,geo.t.cli$Histology)
head(geo.t.cli)

table(geo.t.cli$batch,geo.t.cli$Histology)
# LCC LCNE LUAD LUSC Other
# GSE30219   3   56   83   61    65
# GSE31210   0    0  226    0     0
# GSE37745  24    0  106   66     0
# GSE50081   7    0  127   43     4
# GSE72094   0    0  398    0     0
###
save(geo.t.cli,file = 'raw_datas/GEO/geo.t.cli.RData')

write.table(geo.t.cli,file = 'raw_datas/Processed/geo.t.cli.txt',sep = "\t",quote = F,row.names = T,col.names = T)
write.table(gse72094.t.cli,file = 'raw_datas/Processed/gse72094.t.cli.txt',sep = "\t",quote = F,row.names = T,col.names = T)
write.table(gse31210.t.cli,file = 'raw_datas/Processed/gse31210.t.cli.txt',sep = "\t",quote = F,row.names = T,col.names = T)
write.table(gse50081.t.cli,file = 'raw_datas/Processed/gse50081.t.cli.txt',sep = "\t",quote = F,row.names = T,col.names = T)
write.table(gse30219.t.cli,file = 'raw_datas/Processed/gse30219.t.cli.txt',sep = "\t",quote = F,row.names = T,col.names = T)
write.table(gse37745.t.cli,file = 'raw_datas/Processed/gse37745.t.cli.txt',sep = "\t",quote = F,row.names = T,col.names = T)

save(gse72094.t.cli,file = 'raw_datas/GEO/gse72094.t.cli.RData')
save(gse31210.t.cli,file = 'raw_datas/GEO/gse31210.t.cli.RData')
save(gse50081.t.cli,file = 'raw_datas/GEO/gse50081.t.cli.RData')
save(gse30219.t.cli,file = 'raw_datas/GEO/gse30219.t.cli.RData')
save(gse37745.t.cli,file = 'raw_datas/GEO/gse37745.t.cli.RData')

save(gse72094.t.exp,file = 'raw_datas/GEO/gse72094.t.exp.RData')
save(gse31210.t.exp,file = 'raw_datas/GEO/gse31210.t.exp.RData')
save(gse50081.t.exp,file = 'raw_datas/GEO/gse50081.t.exp.RData')
save(gse30219.t.exp,file = 'raw_datas/GEO/gse30219.t.exp.RData')
save(gse37745.t.exp,file = 'raw_datas/GEO/gse37745.t.exp.RData')

########## END
save.image(file = 'proj_GEO_data_preprocessed.RData')
