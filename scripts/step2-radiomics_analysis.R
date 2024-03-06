
id=openxlsx::read.xlsx('raw_datas/TCGA-LUAD_radiomics_2.xlsx',colNames =F,rowNames =F)
id=data.frame(id,check.names = F,stringsAsFactors = F)
head(id)
data.frame(table(id$X3))
# Var1 Freq
# 1    N    6
# 2    Y   63
####
luad.radiomics.dat=openxlsx::read.xlsx('raw_datas/TCGA-LUAD_radiomics_1.xlsx',colNames =T,rowNames =T)
luad.radiomics.dat=data.frame(luad.radiomics.dat,stringsAsFactors = F,check.names = F)
dim(luad.radiomics.dat)
head(colnames(luad.radiomics.dat))
luad.radiomics.dat=luad.radiomics.dat[,-grep("diagnostics",colnames(luad.radiomics.dat))]
rownames(luad.radiomics.dat)=id$X2[match(rownames(luad.radiomics.dat),id$X1)]
dim(luad.radiomics.dat) ## 63*1130
luad.radiomics.dat=luad.radiomics.dat[,grep("^original_",colnames(luad.radiomics.dat))]
dim(luad.radiomics.dat) ## 63*107
luad.radiomics.dat[1:4,1:4]
table(apply(luad.radiomics.dat,2,sd)==0)
# FALSE  TRUE 
# 105     2 
colnames(luad.radiomics.dat)[which(apply(luad.radiomics.dat,2,sd)==0)]
(luad.radiomics.dat)[,which(apply(luad.radiomics.dat,2,sd)==0)]

luad.radiomics.dat=luad.radiomics.dat[,-which(apply(luad.radiomics.dat,2,sd)==0)]
luad.radiomics.dat=apply(luad.radiomics.dat, c(1,2), function(x){as.numeric(x)})
rownames(luad.radiomics.dat)=paste0(rownames(luad.radiomics.dat),"-01")
dim(luad.radiomics.dat) ## 63*105

####################
mg_lasso_cox_use=function(dat,time,event,nfolds=3,lambda.min=T,show_text=T,figLabels=c('A','B')){
  library("glmnet") 
  library('survival')
  t.inds=which(!is.na(time)&!is.na(event)&time>0)
  dat=dat[t.inds,]
  time=as.numeric(time[t.inds])
  event=as.numeric(event[t.inds])
  y=Surv(time,event)
  set.seed(89898989)
  fit1_cv = cv.glmnet(as.matrix(dat), y, family = "cox", nfolds=nfolds
                      ,nlambda=10000, alpha=1
  )
  fit<-glmnet(dat, y, family = "cox")
  if(lambda.min){
    lambda=fit1_cv$lambda.min
  }else{
    lambda=fit1_cv$lambda.1se
  }
  coefficients<-coef(fit,s=lambda)
  Active.Index<-which(coefficients[,1]!=0)
  genes=row.names(coefficients)[Active.Index]
  print(genes)
  Active.coefficients<-coefficients[Active.Index]  
  g=mg_plot_lasso(fit,fit1_cv,lambda = lambda,show_text=show_text,figLabels=figLabels)
  return(list(Mode1=fit,Model2=fit1_cv,Genes=genes,Coef=Active.coefficients,lambda=lambda,plot=g))
}

mg_plot_lasso_use <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
  if(is.null(lambda)){
    lmda=cv_fit$lambda.min
  }else{
    lmda=lambda
  }
  fit.coef=fit$beta[(apply(fit$beta,1,function(x){
    return(sum(x!=0))
  })>0),]
  
  fit.coef=as.matrix(fit.coef)
  colnames(fit.coef)=fit$lambda
  #fit$lambda==cv_fit$lambda
  library(ggplot2)
  dat=data.table::melt(t(as.matrix(fit.coef)))
  dat_z=dat[which(dat$value==0),]
  dat=dat[which(dat$value!=0),]
  dat.sv=rbind()
  for (u in unique(dat_z[,2])) {
    t.z=dat_z[which(dat_z[,2]==u),1]
    t.zx=max(t.z)
    dat.sv=rbind(dat.sv,c(t.zx,u,0))
    t.zn=min(t.z)
    if(t.zx!=t.zn){
      dat.sv=rbind(dat.sv,c(t.zn,u,0))
    }
  }
  colnames(dat.sv)=colnames(dat_z)
  #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
  dat=crbind2DataFrame(rbind(dat,dat.sv))
  mn=min(-log(dat$Var1))
  mx=max(-log(dat$Var1))
  if(show_text){
    mx=(mx-mn)*0.1+mx
  }
  p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
  p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
  if(show_text){
    fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
    for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
    p=p+ggrepel::geom_label_repel(
      aes(label = Var2,color=Var2),
      data = for_label,hjust = 0
    )
  }
  p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
  p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
  tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                 ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
  p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
    geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
    geom_point(aes(colour=col))
  p1=p1+theme_bw()+theme(legend.position = "none")
  # gal=p+p1
  gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                        ,align = "hv"
                        ,labels = figLabels)
  return(gal)
}

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

################################################
load("raw_datas/TCGA/luad.tcga.exp.RData")
load("raw_datas/TCGA/tcga.t.cli_use.RData")
intersect(rownames(luad.radiomics.dat),substr(colnames(luad.tcga.exp),1,15)) ##29
comm=intersect(rownames(luad.radiomics.dat),substr(rownames(tcga.t.cli_use),1,15)) ##29
luad.radiomics.dat=luad.radiomics.dat[comm,]
dim(luad.radiomics.dat)
luad.radiomics.cli=tcga.t.cli_use[comm,]
################################################  单因素cox分析
tcga.rad.cox.res=cox_batch(t((luad.radiomics.dat)),time = luad.radiomics.cli$OS.time,event = luad.radiomics.cli$OS)
tcga.rad.cox.res
table(tcga.rad.cox.res$p.value<0.05)
# FALSE  TRUE 
# 101     4 
################################################
R_model_data=cbind(data.frame(ID=rownames(luad.radiomics.cli),luad.radiomics.cli[,c("OS.time","OS")],stringsAsFactors = F,check.names = F),luad.radiomics.dat)
##################
# "UBE2T","TEDC2","RCC1","FAM136A"
final.genes=c("UBE2T","TEDC2","RCC1","FAM136A")
# UBE2T: 0个特征
# TEDC2:7个特征
# RCC1: 1个特征：original_glszm_SmallAreaEmphasis
# FAM136A:0个特征
x=as.numeric(tcga.t.exp['TEDC2',rownames(luad.radiomics.cli)])
tcga.rad.pheno=ifelse(x>median(x),"High","Low")
tcga.rad.pheno=factor(tcga.rad.pheno,levels = c("Low","High"))
tcga.rad.pheno

library(glmnet)
set.seed(123)
fit=glmnet(x = scale(luad.radiomics.dat),y = tcga.rad.pheno,
           family = binomial(link = 'logit'),
           nlambda = 1000,
           alpha = 1) 
set.seed(123)
cv.fit<-cv.glmnet(x = scale(luad.radiomics.dat),
                  y = tcga.rad.pheno, 
                  family = binomial(link = 'logit'),
                  nlambda = 1000,
                  nfolds = 5,
                  type.measure = 'deviance',
                  alpha = 1)
cv.fit$lambda.min ### 0.06387721

dev.off()
pdf('analysis/08_Radiomics_model/tcga.radi.lasso.pdf',height = 5,width = 10,onefile = F)
par(mfrow=c(1,2))
plot(fit,xvar="lambda",label=T)
plot(cv.fit)
dev.off()

lambda=cv.fit$lambda.min
coefficients<-coef(fit,s=lambda)
Active.Index<-which(coefficients[,1]!=0)
features=row.names(coefficients)[Active.Index]
Active.coefficients<-coefficients[Active.Index]  
###
features=features[-1]
Active.coefficients=Active.coefficients[-1]

LASSO_features=features
LASSO_features
# [1] "original_shape_Elongation"                   
# [2] "original_firstorder_Median"                  
# [3] "original_firstorder_TotalEnergy"             
# [4] "original_glrlm_RunLengthNonUniformity"       
# [5] "original_glszm_LargeAreaEmphasis"            
# [6] "original_glszm_SmallAreaEmphasis"            
# [7] "original_glszm_SmallAreaLowGrayLevelEmphasis"

mg_violin_use=function(data,group.col=NULL,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL){
  library(ggplot2)
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }else if(legend.pos=='none'){
    pos='none'
  }
  
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  if(!is.null(ylim)){
    data_m=data_m[data_m[,2]<=max(ylim),]
    data_m=data_m[which(!is.na(data_m[,1])),]
    ylim[2]=1.2*ylim[2]
  }
  ######## violin
  p1<-ggplot(data_m,aes(x=Group,y=value,color=Group,fill=Group))+
    geom_violin(width=0.5, alpha=0.7)
  
  if(ct<=8){
    if(!is.null(group.col)){
      p1=p1+scale_color_manual(values = group.col)
      p1=p1+scale_fill_manual(values = group.col)
    }else{
      p1=p1+scale_fill_manual(values=ggsci::pal_d3("category20", alpha = 0.6)(9))
    }
  }else if(ct<=10){
    # p1=p1+ggsci::scale_fill_npg(name=leg.title)
    if(!is.null(group.col)){
      p1=p1+scale_fill_manual(values = group.col)
    }else{
      p1=p1+scale_fill_manual(values = c(pal_lancet('lanonc',alpha =0.8)(9)[c(2,1,4,3,5:6)]))
    }
    
  }else if(ct<=20){
    p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  }else if(ct<=30){
    cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }else if(ct<=38){
    cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10)
                ,ggsci::pal_d3("category20", alpha = 0.6)(20)
                ,ggsci::pal_nejm("default", alpha = 0.6)(8))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }
  
  ####### boxplot
  p1=p1+geom_boxplot(width=0.2,fill="white",outlier.shape = NA)
  #######
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.6,shape=21, size=1,show.legend=FALSE,width = 0.15)
    }else{
      p1<-p1+geom_jitter(alpha=0.6,shape=21, size=1,show.legend=FALSE,width = 0.15,size=point_size)
    }
  }
  
  p1=p1+theme_bw()
  p1=p1+theme(axis.text.x=tx, #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
              axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
              axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
              #panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
              legend.text=element_text(face="plain", family="Times", colour="black"  #设置图例的子标题的字体属性
              ),
              legend.title=element_text(face="plain", family="Times", colour="black" #设置图例的总标题的字体属性
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
              ,panel.grid.major = element_blank(),   #不显示网格线
              panel.grid.minor = element_blank()
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- anova(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til)
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    if(length(comps)<7){
      p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
    }
  }
  return(p1)
}

dt=data.frame(luad.radiomics.dat[,LASSO_features],class=tcga.rad.pheno)
table(dt$class)

p.all=list()
for(i in LASSO_features){
  p=mg_violin_use(data = data.frame(dt[,c('class',i)])
                  ,melt = T
                  ,group.col = clin.color
                  ,xlab = ''
                  ,ylab = i
                  ,jitter=F
                  ,test_method = 'kruskal.test'
                  ,cmp_test_method = 'wilcox.test'
                  # ,legend.pos = NULL
                  ,show_compare = T)+labs(fill="TEDC2",color="TEDC2")
  p.all=c(p.all,list(p))
  
}
length(p.all)
figure=mg_merge_plot(p.all,nrow = 2,ncol = 4,common.legend = T)
figure
savePDF('analysis/08_Radiomics_model/tcga.radi.features.violin.pdf',figure,height = 6,width = 10)
################ R_model_data
R_model_data=data.frame(scale(luad.radiomics.dat[,LASSO_features]),class=tcga.rad.pheno)
R_model_data$class=factor(R_model_data$class,levels = c('Low','High'))
dim(R_model_data) ### 29*13
###############
library(caret)
set.seed(20231212)
trainIndex <- createDataPartition(R_model_data$class, p = 0.6,
                                  list = FALSE,
                                  times = 1)

head(trainIndex)
TRAIN=R_model_data[trainIndex,]
TEST=R_model_data[-trainIndex,]
table(TRAIN$class)
table(TEST$class)

###################
trControl <- trainControl(method = 'repeatedcv',
                          number = 3,
                          repeats =  5,
                          classProbs = T,
                          search = 'random')

logit.CV <- caret::train(class ~ . , data = TRAIN,
                         method = 'glmnet',
                         trControl = trControl,
                         family = 'binomial' )
logit.CV

plot(logit.CV)

dim(TRAIN)
class <- predict(logit.CV, TRAIN[,-ncol(TRAIN)],'raw')
probs <- predict(logit.CV, TRAIN[,-ncol(TRAIN)],'prob')
TRAIN$class

pred=class
truth=TRAIN$class
xtab <- table(pred, truth)
confusionMatrix(xtab,positive ="High")

###################
library(ROCR)
pred <-prediction(probs[,1],TRAIN$class)
############ ROC curve
perf <- performance(pred,"tpr","fpr")
pdf('analysis/08_Radiomics_model/TRAIN.radi.ROC.pdf',height = 5,width = 5)
plot(perf, col='blue',lty=2)
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
auc=round(auc,digits = 2)
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

############ PR curve
perf <- performance(pred,"prec","rec")
pdf('analysis/08_Radiomics_model/TRAIN.radi.PR.pdf',height = 5,width = 5)
plot(perf, col='blue',lty=2)
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
auc=round(auc,digits = 2)
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("PR curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(1,-1)
dev.off()
########################## 验证集
class <- predict(logit.CV, TEST)
probs <- predict(logit.CV, TEST,'prob')

threshold <- 0.6
class <- ifelse(probs[,1] >= threshold, "Low", "High")
class <- factor(class,levels = c("Low","High"))

pred=class
truth=TEST$class
xtab <- table(pred, truth)
xtab
confusionMatrix(xtab,positive ="High")

library(ROCR)
pred <-prediction(probs[,1],TEST$class)
############ ROC curve
perf <- performance(pred,"tpr","fpr")
plot(perf, col='blue',lty=2)
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
auc=round(auc,digits = 2)
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
############ PR curve
perf <- performance(pred,"prec","rec")
plot(perf, col='blue',lty=2)
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
auc=round(auc,digits = 2)
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("PR curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(1,-1)

###################DCA################################
library(rms)
library(rmda)

#临床决策曲线DCA
traindata=R_model_data[trainIndex,]
table(traindata$class)
traindata$class=factor(traindata$class,levels = c("High","Low"),labels = c(1,0))
traindata$class=as.numeric(as.character(traindata$class))
#校准曲线
dd <- datadist(traindata)
options(datadist="dd")
fmla <- as.formula(paste0("class ~",paste0(colnames(traindata)[-ncol(traindata)],collapse = '+')))
fmla

# 进行hosmer-lemeshow 检验
library(ResourceSelection)
model_glm <- glm(fmla,data = traindata, family = binomial)

# hosmer-lemeshow 检验
p.hoslem <- hoslem.test(model_glm$y, fitted(model_glm), g=10)$p.value
p.hoslem

fit1 <- lrm(fmla,data = traindata,x=T,y=T)
cal1 <- calibrate(fit1, method='boot', B=500)

pdf("analysis/08_Radiomics_model/tcga.Calibration.TRAIN.pdf",width = 5,height = 5)
plot(cal1,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "Prediced Probability",
     ylab = "Observed Probability",
     legend = FALSE)
lines(cal1[,c("predy","calibrated.corrected")], 
      type = 'l', #连线的类型，可以是"p","b","o"
      lwd = 3, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = "#2166AC") #连线的颜色
lines(cal1[,c("predy","calibrated.orig")],type="l",pch=16,lwd=3,col="tomato")
abline(0,1,
       lty = 2, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = "#224444")#对角线的颜色
legend(0.6,0.2,
       c("Apparent","Bias-corrected","Ideal"), 
       lty = c(2,1,1), 
       lwd = c(2,3,3), 
       col = c("black","#2166AC","tomato"), 
       bty = "n"
)
text(0,0,bquote("Hosmer-Lemeshow "~italic(P)~" = "~.(round(p.hoslem,3))),adj = 0)
dev.off()

#### 测试集的校准曲线
testdata=TEST
testdata$class=factor(testdata$class,levels = c("High","Low"),labels = c(1,0))
testdata$class=as.numeric(as.character(testdata$class))
# 首先获取测试集的预测结果
probability <- predict(fit1, testdata, type = 'fitted')
# 直接使用val.prob即可实现
val.prob(probability, testdata$class,statloc = F,cex = 1)
############################### 临床决策曲线分析
fmla
fit2<- decision_curve(formula = fmla,
                       data = traindata, 
                       family = binomial(link ='logit'))
pdf('analysis/08_Radiomics_model/tcga.DCA.TRAIN.pdf',height = 5,width = 5)
plot_decision_curve(fit2,
                    curve.names="DCA model",xlab="Threshold probability",
                    cost.benefit.axis =FALSE,col= "#E73F74",
                    confidence.intervals=FALSE,
                    standardize = FALSE)

dev.off()

############################################
logit.CV
fit1 = glmnet::glmnet(x = as.matrix(TRAIN[,-ncol(TRAIN)]), y=TRAIN$class, family=binomial(link='logit'),alpha=logit.CV$bestTune$alpha,lambda=logit.CV$bestTune$lambda)
rs = predict.glmnet(fit1, newx = as.matrix(R_model_data[,-ncol(R_model_data)]),s='lambda.min',type='link')
rs

x=names(coef(fit1)[,1])[2:8]
score.manual=as.matrix(R_model_data[,x]) %*% as.numeric(coef(fit1)[,1][2:8])+as.numeric(coef(fit1)[,1][1])
head(cbind(score.manual,score.TCGA))

# mosaic::zscore()
score.TCGA=data.frame(luad.radiomics.cli[,c('OS.time','OS')],riskscore=(rs[,1]))
library(survminer)
res.cut <- surv_cutpoint(score.TCGA, time = "OS.time", event = "OS",
                         variables = c('riskscore'))
print(res.cut$cutpoint$cutpoint) # -1.965349
score.TCGA$risktype=ifelse(score.TCGA$riskscore>res.cut$cutpoint$cutpoint,"High","Low")

head(score.TCGA)

p1=ggplotTimeROC(time = score.TCGA$OS.time,
                 status = score.TCGA$OS,
                 score = score.TCGA$riskscore,mks = c(1,3,5))+labs(title = 'TCGA')
p1


p2=ggplotKMCox(data.frame(time = score.TCGA$OS.time/365
                          , event = score.TCGA$OS
                          , groups=score.TCGA$risktype)
               , legend.title = 'Radiomic Score'
               , title = 'TCGA'
               , pal = clin.color
               , labs = c('Low','High')
               , add_text = '')
p2

plot1=mg_merge_plot(p1,p2,nrow = 1,ncol = 2)
plot1
# savePDF('analysis/08_Radiomics_model/radi.model.timeROC.pdf',plot1,height = 5,width = 5)
###################
dt=cbind(R_model_data,score.TCGA)
p1=mg_violin_use(data = data.frame(dt[trainIndex,c('class','riskscore')])
                ,melt = T
                ,group.col = clin.color
                ,xlab = ''
                ,ylab = 'Radiomic Score'
                ,jitter=F
                ,test_method = 'kruskal.test'
                ,cmp_test_method = 'wilcox.test'
                ,legend.pos = NULL
                ,show_compare = T)
p1

p2=mg_violin_use(data = data.frame(dt[-trainIndex,c('class','riskscore')])
                 ,melt = T
                 ,group.col = clin.color
                 ,xlab = ''
                 ,ylab = 'Radiomic Score'
                 ,jitter=F
                 ,test_method = 'kruskal.test'
                 ,cmp_test_method = 'wilcox.test'
                 ,legend.pos = NULL
                 ,show_compare = T)
p2
plot2=mg_merge_plot(p1,p2,nrow = 1,ncol = 2)
# savePDF('analysis/08_Radiomics_model/radi.model.marker.violin.pdf',plot2,height = 5,width = 10)

figure=mg_merge_plot(plot1,plot2,nrow = 2,ncol = 1)
savePDF('analysis/08_Radiomics_model/radi.model.prognosis.pdf',figure,height = 10,width = 10)

############################### 预后模型
R_clinical=data.frame(tcga_clinical[rownames(score.TCGA),],score.TCGA[,c("riskscore","risktype")])
dim(R_clinical)

univar_res<-unicox(vars=c("T.Stage","N.Stage","M.Stage","Stage_2","Age","Gender","Smoker_2",'riskscore'),time = R_clinical$OS.time,event = R_clinical$OS,data=R_clinical)
univar_res
rownames(univar_res)[4]='Stage'
rownames(univar_res)[7]='Smoker'
table(univar_res$pvalue<0.05)
rownames(univar_res)[which(univar_res$pvalue<0.05)]
# "N.Stage" "M.Stage" "Stage" 
## 多因素
mutivar_res<-multicox(vars=c("Stage_2","Age","Gender","Smoker_2",'riskscore'),time = R_clinical$OS.time,event = R_clinical$OS,data=R_clinical,forest = F)
mutivar_res
rownames(mutivar_res)[mutivar_res$pvalue<0.05]
# "Stage_2"
save.image(file = '20231116_LUAD_Radiomics_all_final.RData')

