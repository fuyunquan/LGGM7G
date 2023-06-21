#筛选与巨噬细胞相关的m7G基因做候选基因集，进行模型筛选

#最后选择第一次批次矫正后TCGA+GEO合并的修正模型

# 2022.08.17 一个惨痛的教训，修正集的预测因子漏了一个指标，算重复了一个。
# 后续分析所有手算lp的图都要重画了 46张图，删图只需要5秒，画图需要两个月
# 所以，建议让代码自动换行！！！！！检查好代码。太痛了


rm(list=ls()) #删除工作历史的所有变量
ls()


setwd("E:/A_sudakeyan/LGGm7G")#指定工作目录
getwd()


library(glmnet)
library(survival)
library(rms)
library(MASS)

library(car)
library(ggplot2)
library(glmnet)
library(ggrisk)

library(survminer)
library(ggridges)
library(pROC)
library(plotROC)
library(riskRegression)
library(survivalROC)
library(pec)


#合并GEO的数据集后少了几个基因，但是筛出来的结果不影响
load("./rdata/tcga_data1.RData")

#TCGA 71 gene
#load("./rdata/tcga_all.RData")

#巨噬细胞相关
load("./rdata/tcag_m7g_mdeg.RData")

genenames <- tcag_m7g_mdeg[,1]

# 问了CGGA数据官方，微阵列数据已经做了标准化处理，不建议合并其他数据集分析
#usedata1 <- tcga_data1[genenames]#报错，CGGA和GEO的基因没对上！！！

#先简单对一下吧
genenames2 <- colnames(tcga_data1)[-c(1:12)]

genecommon <- intersect(genenames,genenames2)
#19个，就少了一个
#不含CGGA微阵列301后，20个全在 TCGA GEO CGGA693+325
#21 无GEO

usedata1 <- tcga_data1[genecommon]

usedata1 <- cbind(tcga_data1[1:12],usedata1)


##############lasso###################
###########    OS  #############

x <- as.matrix(usedata1[,-c(1:12)])
y <- Surv(usedata1$OS.time,usedata1$OS)

fit <-glmnet(x, y, family = "cox")
plot(fit, xvar = "lambda", label = TRUE)
print(fit)
#交叉验证
set.seed(1234)
cv.fit <- cv.glmnet(x, y, family="cox")
abline(v=log(c(cv.fit$lambda.min,cv.fit$lambda.1se)),lty=2)
plot(cv.fit)

#取最小值
cv.fit$lambda.min
Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Index
Active.Coefficients
row.names(Coefficients)[Active.Index]
# 9个
# 5个  ……?
# 还是5个 一模一样 TCGA CGGA

#[1] "ENSG00000035403" "ENSG00000104894" "ENSG00000107581" "ENSG00000114867"
#[5] "ENSG00000128595" "ENSG00000142864" "ENSG00000147133" "ENSG00000167257"
#[9] "ENSG00000215301"

#[1] "ENSG00000104894" "ENSG00000107581" "ENSG00000128595"
#[4] "ENSG00000130309" "ENSG00000215301"

#拟合生存回归模型 这里是原模型
model1.norm <- cph(Surv(usedata1$OS.time,usedata1$OS) ~ ENSG00000104894 + 
                   ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
                   ENSG00000215301,
                 data = usedata1,x=T,y=T,surv=T)


###  重新导入代码9.3 中去除批次效应的TCGA数据 

#load("./rdata/tcga_data2.RData") 这里是只看去除批次效应后TCGA的结果 不影响模型内部验证性能

#第二次去除批次效应，不带GSE16011
#load("./rdata/tcga_all2.RData")

#usedata1 <- tcga_all2

#拟合生存回归模型
#model1.norm<-cph(Surv(usedata1$OS.time,usedata1$OS) ~ ENSG00000104894 + 
 #                  ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
  #                 ENSG00000215301,
   #              data = usedata1,x=T,y=T,surv=T)

#c-index.bootstrap
# set.seed(1234)
# model1_c<-validate(model1.norm,B=1000,dxy = T)
# slope 0.9399
#(model1_c[1,1]+1)/2
#0.7939374
#0.7868049
#0.7868049 2022.06.20 去除批次效应后不影响

#(model1_c[1,5]+1)/2
#0.7793071
#0.7792708 ……区别竟然不大
#0.7794153

#C指数置信区间
library(survcomp)

cindex <- concordance.index(predict(model1.norm),surv.time = usedata1$OS.time,
                            surv.event = usedata1$OS,method = "noether")

cindex$c.index; cindex$lower; cindex$upper
# [1] 0.786984
# [1] 0.7415874
# [1] 0.8323806

###### 内外验证集都要去 去除批次效应后再来看结果


#### 2022.06.20 第二次去除批次效应 模型第一次重校准：CGGA校准斜率*模型系数
### 废弃不用 下面那个方法结果更好
#提取系数：
#coef <- as.data.frame(model1.norm$coefficients)

#调整系数：slope看外部验证
#adjustcoef <- slope*coef

#usedata1$lp <- 0.1446081*usedata1$ENSG00000104894 - 0.4797180*usedata1$ENSG00000107581 + 
 # 0.2905599*usedata1$ENSG00000128595 + 0.2400222*usedata1$ENSG00000128595 + 0.3123324*usedata1$ENSG00000215301

#model1.adj <- cph(Surv(usedata1$OS.time,usedata1$OS)~lp,
 #                 data = usedata1,x=T,y=T,surv=T)

####  废弃不用↑ 


####2022.06.21 09:52  第一次去除批次效应 模型修正 : 合并TCGA和GSE16011,重新评估系数
 
### 最后决定用这个结果，因为GEO、CGGA校准都不错

# #导入第一次的俩个数据集
# load("./rdata/tcga_data2_2.RData")
# load("./rdata/gse16011_m7ggene_surv_no2.RData")
# 
# 
# use_tcga <- tcga_data2[,c(1,3,4,13:17)]
# colnames(gse16011_m7ggene_surv_no2)[1] <- "sample"
# gse16011_m7ggene_surv_no2$OS.time <- gse16011_m7ggene_surv_no2$OS.time*365
# 
# usedata1 <- rbind(use_tcga,gse16011_m7ggene_surv_no2)
# 
# GEO_revisionset <- usedata1
# save(GEO_revisionset,file = "./rdata/GEO_revisionset2.RData")

load("./rdata/GEO_revisionset2.RData")

usedata1 <- GEO_revisionset

#拟合生存回归模型
model1.adj <- cph(Surv(usedata1$OS.time,usedata1$OS) ~ ENSG00000104894 + 
                   ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
                   ENSG00000215301,
                 data = usedata1,x=T,y=T,surv=T)

model1.adj$coefficients # 这是重新估计出来的系数

#  重新在推导集中拟合模型
usedata3 <- tcga_data2

#提取系数，手算lp
#usedata3$lp <- usedata3$ENSG00000104894*0.01511283 - usedata3$ENSG00000107581*0.49321889 + 
 # usedata3$ENSG00000128595*0.28688898 + usedata3$ENSG00000130309*0.44763118 + usedata3$ENSG00000215301*0.20994054

#GEO的非分位数标准化后，系数都不一样了，救命
usedata3$lp <- -usedata3$ENSG00000104894*0.009889326 - 
  usedata3$ENSG00000107581*0.494768762 + 
  usedata3$ENSG00000128595*0.361307666 + 
  usedata3$ENSG00000130309*0.346389452 + 
  usedata3$ENSG00000215301*0.131503000


# 
model1.adj <- cph(Surv(usedata3$OS.time,usedata3$OS)~lp,
                  data = usedata3,x=T,y=T,surv=T)
# 
# #c-index.bootstrap
# set.seed(1234)
# model1_c<-validate(model1.adj,B=1000,dxy = T)
# 
# (model1_c[1,1]+1)/2
# # 0.7808983 TCGA
# # 0.7851332 GEO NO quanli 算错了x
# # 0.7915227 GEO NO quanli 正确的lp √
# 
# (model1_c[1,5]+1)/2
# # 0.7805978 TCGA
# # 0.7850061 GEO NO quanli 算错了x
# # 0.7915079 GEO NO quanli 正确的lp √
# 
# model1_c
# # 所以这里C指数应该是指纯TCGA数据集
# # 这个结果的校准曲线看起来还不错，先在9.4中验证一下

### 

###############   下面这段也没用   ###############################
###  2022.07.24 看文献，试一下其他修正方法(花里胡哨的)
### 一个非常头疼的问题，为什么2000,2004,2008三篇文献讲的方法大同小异
## 但是这个不同的点很要命啊，我到底按哪个来？哦，还有篇2012年的
## 2000：保留原来的系数，model做lp，在验证集里筛新变量
## 结合2000年和2008

# 先导入批次矫正后的数据 TCGA
# load("./rdata/tcga_data3.RData")
# 
# usedata1 <- tcga_data3
# # 先在这个数据集里拟合原模型，求lp
# model1 <- cph(Surv(usedata1$OS.time,usedata1$OS) ~ ENSG00000104894 + 
#                     ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
#                     ENSG00000215301,
#                   data = usedata1,x=T,y=T,surv=T)
# 
# model1$coefficients # 和仅有五个基因的矫正后系数不一样了
# 
# usedata1$lp <- 0.1749335*usedata1$ENSG00000104894 - 0.5999493*usedata1$ENSG00000107581 + 
#   0.3594843*usedata1$ENSG00000128595 + 0.3007320 *usedata1$ENSG00000130309 + 
#   0.3900521*usedata1$ENSG00000215301
# 
# # 再在修正集里代入模型lp
# 
# # 在修正集中用lp+剩余预后因子
# # 啊糟了，没有留其他的预后因子做批次矫正。尴尬了，要重新整理数据。去提取那20个预后因子
# load("./rdata/gse16011_m7ggene_surv_no3.RData")
# 
# use_tcga <- tcga_data3[,c(1,3,4,13:32)]
# colnames(gse16011_m7ggene_surv_no3)[1] <- "sample"
# gse16011_m7ggene_surv_no3$OS.time <- gse16011_m7ggene_surv_no3$OS.time*365
# 
# usedata1 <- rbind(use_tcga,gse16011_m7ggene_surv_no3)
# 
# GEO_revisionset3 <- usedata1
# save(GEO_revisionset3,file = "./rdata/GEO_revisionset3.RData")
# 
# 
# # 在修正集中用lp+剩余预后因子
# model1 <- cph(Surv(usedata1$OS.time,usedata1$OS) ~ ENSG00000104894 + 
#                 ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
#                 ENSG00000215301,
#               data = usedata1,x=T,y=T,surv=T)
# 
# model1$coefficients
# 
# usedata1$lp <- 0.03433276*usedata1$ENSG00000104894 - 0.22970824*usedata1$ENSG00000107581 + 
#   0.36686929*usedata1$ENSG00000128595 + 0.27377293 *usedata1$ENSG00000130309 + 
#   0.27377293*usedata1$ENSG00000215301
# 
# usedata2 <- subset(usedata1, select=-c(ENSG00000104894,ENSG00000107581,ENSG00000128595,
#                                        ENSG00000130309,ENSG00000215301))
# # 不要用带引号的列名 ！！！！
# 
# 
# # 构建新的备选模型
# model2 <- cph(Surv(OS.time,OS) ~ lp + ENSG00000009307 + ENSG00000035403 +
#                 ENSG00000068650 + ENSG00000108424 + ENSG00000108654 +
#                 ENSG00000111142 + ENSG00000114867 + ENSG00000124151 + 
#                 ENSG00000142864 + ENSG00000147133 + ENSG00000167257 + 
#                 ENSG00000169813 + ENSG00000181904 + ENSG00000205155,
#               data = usedata2,x=T,y=T,surv=T)
# 
# # 用逐步回归：
# model2_step <- stepAIC(model2,direction = "both")
# 
# model2_step <- cph(formula = Surv(OS.time, OS) ~ lp + ENSG00000009307 + ENSG00000108424 + 
#                      ENSG00000114867 + ENSG00000167257 + ENSG00000169813, data = usedata2, x = T, y = T, surv = T)
# 
# 
# 
# #c-index.bootstrap
# set.seed(1234)
# model2_c <- validate(model2_step,B=1000,dxy = T)
# 
# (model2_c[1,1]+1)/2
# # 0.7325667
# 
# (model2_c[1,5]+1)/2
# # 0.7255729  结果没上一个好啊
# 
# usedata3 <- tcga_data3
# 
# usedata3$lp <- 0.03433276*usedata3$ENSG00000104894 - 0.22970824*usedata3$ENSG00000107581 + 
#   0.36686929*usedata3$ENSG00000128595 + 0.27377293 *usedata3$ENSG00000130309 + 
#   0.27377293*usedata3$ENSG00000215301
# 
# model3 <- cph(formula = Surv(OS.time, OS) ~ lp ,
#               data = usedata3, x = T, y = T, surv = T)
# 
# set.seed(1234)
# model3_c <- validate(model3,B=1000,dxy = T)
# 
# (model3_c[1,1]+1)/2
# # 0.78309
# 
# (model3_c[1,5]+1)/2
# 0.7828776  ???怎么这个高 哦，和上一个没什么差别

###########  那这个最后也没什么用啊 ##################



### 2022.06.21 10:45  第二次批次效应去掉GSE16011 TCGA+CGGA 325 合并，重新评估系数
### 废弃不用 上面那个结果就挺好
# 导入两个数据
#load("./rdata/tcga_all2.RData")
#load("./rdata/cgga_all2.RData")

#提取tcga数据集
#use_tcga <- tcga_all2[,c(1,3,4,13:17)]

# 提取325数据集
#rownames(cgga_all2) <- cgga_all2$CGGA_ID

#load("./rdata/cgga325_all.RData")
#cgga_sample <- cgga325_all$CGGA_ID

#use_cgga <- cgga_all2[cgga_sample,c(1,6,5,12:16)]
#colnames(use_cgga)[1] <- "sample"


# 合并数据集
#usedata1 <- rbind(use_tcga,use_cgga)


#拟合生存回归模型
#model1.adj <- cph(Surv(usedata1$OS.time,usedata1$OS) ~ ENSG00000104894 + 
 #                   ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
  #                  ENSG00000215301,
   #               data = usedata1,x=T,y=T,surv=T)

#model1.adj$coefficients

#  重新在推导集中拟合模型
#usedata3 <- tcga_all2

#提取系数，手算lp
#usedata3$lp <- usedata3$ENSG00000104894*0.1484534 - usedata3$ENSG00000107581*0.4583745 + 
 # usedata3$ENSG00000128595*0.3223749 + usedata3$ENSG00000130309*0.3047288 + usedata3$ENSG00000130309*0.2308864


#model1.adj <- cph(Surv(usedata3$OS.time,usedata3$OS)~lp,
 #                 data = usedata3,x=T,y=T,surv=T)

#c-index.bootstrap
#model1_c<-validate(model1.adj,B=1000,dxy = T)

#(model1_c[1,1]+1)/2
# 0.7754374
#(model1_c[1,5]+1)/2
# 0.7746758

####   废弃不用↑



#############  模型粗表现 survival ROC  ################################
#nobs <- NROW(usedata1)
cutoff1 <- 365 #这里就是时间节点的设置
cutoff2 <-1095
cutoff3 <-1825

# 原模型
# usedata3 <- tcga_data2
# model1.norm <- cph(Surv(OS.time,OS) ~ ENSG00000104894 + 
#                      ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
#                      ENSG00000215301,
#                    data = usedata3,x=T,y=T,surv=T)
# 
# usedata3$lp <- predict(model1.norm, newdata=usedata3, type="lp")
# 


#修正模型，上面手算lp  后面看model1.adj
#usedata1$lp <- predict(model1.norm, newdata=usedata1, type="lp") 

#修正模型  TCGA数据集存的是usedata3
summary(usedata3$lp)


#1-year
Mayo4.1= survivalROC(Stime=usedata3$OS.time,##生存时间
                     status=usedata3$OS,## 终止事件    
                     marker = usedata3$lp, ## marker value    
                     predict.time = cutoff1,## 预测时间截点
                     method = 'KM')##span,NNE法的namda

str(Mayo4.1)## list结构,是每一个marker的cutoff值都计算出相应的TP,FP

#3-year
Mayo4.2 <- survivalROC(Stime=usedata3$OS.time,##生存时间
                     status=usedata3$OS,## 终止事件    
                     marker = usedata3$lp, ## marker value    
                     predict.time = cutoff2,## 预测时间截点
                     method = 'KM')

#5-year
Mayo4.3 <- survivalROC(Stime=usedata3$OS.time,##生存时间
                     status=usedata3$OS,## 终止事件    
                     marker = usedata3$lp, ## marker value    
                     predict.time = cutoff3,## 预测时间截点
                     method = 'KM')

#画图
plot(Mayo4.1$FP, Mayo4.1$TP, ## x=FP,y=TP
     type="l",col="#BC3C29FF", ##线条设置
     xlim=c(0,1), ylim=c(0,1),  lwd = 1.5,
     xlab=("1 - Specificity"), ##连接
     ylab="Sensitivity",
     main="Time dependent ROC")## \n换行符
abline(0,1,col="gray", lwd = 1.5,lty=2)##线条颜色

#lines函数在原有基础上继续绘图 #E18727FF
#legend函数增加legend
lines(Mayo4.2$FP, Mayo4.2$TP, type="l", lwd = 1.5,col="#0072B5FF",xlim=c(0,1), ylim=c(0,1))
lines(Mayo4.3$FP, Mayo4.3$TP, type="l", lwd = 1.5,col="#E18727FF",xlim=c(0,1), ylim=c(0,1))
legend(0.45,0.3,c(paste("AUC of 1-year =",round(Mayo4.1$AUC,3)),
                  paste("AUC of 3-year =",round(Mayo4.2$AUC,3)),
                  paste("AUC of 5-year =",round(Mayo4.3$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 3,col=c("#BC3C29FF","#0072B5FF",'#E18727FF'),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)# 



# 校准图  在后面统一成手算的格式了

# usedata3 <- tcga_data2
# model1.norm <- cph(Surv(OS.time,OS) ~ ENSG00000104894 + 
#                      ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
#                      ENSG00000215301,
#                    data = usedata3,x=T,y=T,surv=T)
# 
# usedata3$lp <- predict(model1.norm, newdata=usedata3, type="lp") 

#1 year
model1.norm_1year<-cph(Surv(usedata3$OS.time,usedata3$OS)~lp,
                       data = usedata3,x=T,y=T,surv=T,time.inc=cutoff1)

p1<- calibrate(model1.norm_1year,#模型名称
               cmethod='KM',
               method='boot',#检测方法
               u=cutoff1,#评估的时间，注：一定要与模型的时间一致
               m=120, #每次抽样的样本量，
               B=500)#抽样次数
#注，m值的确定：m=数据总数/3-4,即你想让最终的校准曲线有3个点，那就是m=数据总数/3
#B值一般1000，电脑配置不好可以选500,300,100等

plot(p1,
     add=F,#增加第二条线
     conf.int=T,#95%CI
     subtitles = F,#关闭副标题
     cex.subtitles=0.8, #副标题大小
     lwd=2,#95%CI粗细
     lty=1,#95%CI实线，2=虚线
     errbar.col="#76549A",#95%CI颜色
     xlim=c(0.0,1),#x轴范围
     ylim=c(0.0,1),
     xlab="Predicted(%)",
     ylab="Observed(%)",
     col="#76549A",#曲线颜色
     bty="l",#L表示只画左边和底部，象形字(看字知意)
     par.corrected=list(col="#76549A", lty=1, lwd=2, pch=4))

lines(p1[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细  要的就是这个
      pch = 16, #点的形状，可以是0-20
      col = c("#76549A")) #连线的颜色


#3years
model1.norm_3year <- cph(Surv(usedata3$OS.time,usedata3$OS)~lp,
                         data = usedata3,x=T,y=T,surv=T,time.inc=cutoff2)

p2<- calibrate(model1.norm_3year,
               cmethod='KM',
               method='boot',
               u=cutoff2,
               m=120, 
               B=500)

plot(p2,
     add=T,
     conf.int=T,
     subtitles = F,
     #cex.subtitles=0.8, 
     lwd=2,
     lty=1,
     errbar.col="#DF7861",
     xlim=c(0.0,1),
     ylim=c(0.0,1),
     col="#DF7861",
     par.corrected=list(col="#DF7861", lty=1, lwd=2, pch=4))

lines(p2[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细  要的就是这个
      pch = 16, #点的形状，可以是0-20
      col = c("#DF7861")) #连线的颜色


#5 years
model1.norm_5year <- cph(Surv(usedata3$OS.time,usedata3$OS)~lp,
                         data = usedata3,x=T,y=T,surv=T,time.inc=cutoff3)

p3<- calibrate(model1.norm_5year,
               cmethod='KM',
               method='boot',
               u=cutoff3,
               m=120,
               B=500)

plot(p3,
     add=T,
     conf.int=T,
     subtitles = F,
     #cex.subtitles=0.8, 
     lwd=2,
     lty=1,
     errbar.col="#94B49F",
     xlim=c(0.0,1),
     ylim=c(0.0,1),col="#94B49F",
     par.corrected=list(col="#94B49F", lty=1, lwd=2, pch=4))


lines(p3[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细  要的就是这个
      pch = 16, #点的形状，可以是0-20
      col = c("#94B49F")) #连线的颜色



#加上图例
legend("bottomright", cex = 1,
       legend=c("1-year","3-year", "5-year"),
       col=c("#76549A","#DF7861", "#94B49F"), lwd=3,bty="n")
#调整对角线
abline(0,1,lty=3,lwd=1,col="grey")



#计算一下Brier Score.bootstrap
m <- nrow(usedata3)
# 这里是原模型的
model1.norm <- cph(Surv(OS.time,OS) ~ ENSG00000104894 + 
                    ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
                    ENSG00000215301,
                  data = usedata3,x=T,y=T,surv=T)

bs1 <- Score(list("1-year"=model1.norm),
             Surv(OS.time,OS)~1,
             data=usedata3,times=cutoff1,
             plots="Calibration",metrics=c("Brier"),split.method="bootcv",
             B=100,N=m)$Brier$score$Brier[-1]
bs1
# 0.05123262 model1.norm

bs3 <- Score(list("3-year"=model1.norm),
             Surv(OS.time,OS)~1,
             data=usedata3,times=cutoff2,
             plots="Calibration",metrics=c("Brier"),split.method="bootcv",
             B=100,N=m)$Brier$score$Brier[-1]

bs3
# 0.1285301 model1.norm

bs5 <- Score(list("5-year"=model1.norm),
             Surv(OS.time,OS)~1,
             data=usedata3,times=cutoff3,
             plots="Calibration",metrics=c("Brier"),split.method="bootcv",
             B=100,N=m)$Brier$score$Brier[-1]
bs5
# 0.1884266 model1.norm








###################################
###################################
############################# 外部验证 ##################################

#load("./rdata/cgga_data2.RData")
#load("./rdata/gse16011_m7ggene_surv_no2.RData")

#是第二次去除批次效应，不带GEO
#load("./rdata/cgga_all2.RData")

#rownames(cgga_data2) <- cgga_data2$CGGA_ID
#cgga_sample <- cgga693_all$CGGA_ID

#load("./rdata/cgga325_all.RData")
#cgga_sample <- cgga325_all$CGGA_ID


#usedata2 <- cgga_data2[cgga_sample,]
####

## 第一次去除批次效应 TCGA+GEO 合并数据集 重估系数 模型修正
load("./rdata/GEO_revisionset2.RData")
load("./rdata/cgga_data2_2.RData")

usedata2 <- GEO_revisionset
summary(usedata2$OS.time)
#usedata2$OS.time <- usedata2$OS.time*365


#usedata2$lp <- predict(model1.norm, newdata=usedata2, type="lp")

#usedata2$lp <- usedata2$ENSG00000104894*0.01511283 - usedata2$ENSG00000107581*0.49321889 + 
#  usedata2$ENSG00000128595*0.28688898 + usedata2$ENSG00000130309*0.44763118 + usedata2$ENSG00000130309*0.20994054

# 2022.06.24 GEO NO quantil
usedata2$lp <- -usedata2$ENSG00000104894*0.009889326 - 
  usedata2$ENSG00000107581*0.494768762 + 
  usedata2$ENSG00000128595*0.361307666 + 
  usedata2$ENSG00000130309*0.346389452 + 
  usedata2$ENSG00000215301*0.131503000 

# 这重校准一点变化都没有，放弃了
#usedata2$lp <- 0.9053482*usedata2$lp
#usedata2$lp <- 0.9382395*usedata2$lp

summary(usedata2$lp)

c_harrell_external <- (cph(Surv(OS.time,OS)~lp, data=usedata2,x=TRUE,y=TRUE)$stats["Dxy"]+1)/2
c_harrell_external
# 0.6912347 CGGA total
# 0.7169662 GSE16011
# 0.6912127 CGGA total 2

# 0.7588311 TCGA+GEO 0.7643094 no quantil
# 0.7048975 CGGA     0.7107278

# 2022.08.19 重算正确的lp
# 0.7640502 revision(geo no quantil)
# 0.7000132 CGGA

#校准斜率
slope <- cph(Surv(OS.time,OS)~lp, data=usedata2,x=TRUE,y=TRUE)$coef
slope
# 0.8308 CGGA total 
# 0.6669 GSE16011 心累啊
# 0.8308 CGGA total2  那基本就没怎么变

# 0.839722 改了跟没改一样啊气死我了  0.9053482 geo no quantil   TCGA+GEO
# 0.9011807                          0.9382395 CGGA


# 2022.08.19 重算正确的lp
# 0.999842 revision(geo no quantil)
# 0.9978944 CGGA ……这么好啊？

#C指数置信区间
library(survcomp)
fit <- coxph(formula = Surv(usedata2$OS.time, usedata2$OS) ~ lp,
             data = usedata2)
cindex <- concordance.index(predict(fit),surv.time = usedata2$OS.time,
                            surv.event = usedata2$OS,method = "noether")

cindex$c.index; cindex$lower; cindex$upper
# 0.6912024 0.6482143 0.7341904
# 0.7169103 0.6634669 0.7703537

# 0.7588318 0.7233283 0.7943352  
# 0.7048675 0.6637778 0.7459573

# geo no quantil
# 0.7643133 0.7643133 0.7988929
# 0.7106989 0.6704558 0.750942


# 2022.08.19 重算正确的lp
# 0.7640539 0.7298940 0.7982137 REVISION
# 0.6999824 0.6584457 0.7415191 CGGA

# 计算一下Brier Score.bootstrap
m <- nrow(usedata2)

# TCGA
model1.adj <- cph(formula = Surv(OS.time, OS) ~ lp, data = usedata3, 
    x = T, y = T, surv = T)

model1.CGGA <- cph(formula = Surv(OS.time, OS) ~ lp, data = usedata2, 
                  x = T, y = T, surv = T)

bs1 <- Score(list("1-year"=model1.CGGA),
             Surv(OS.time,OS)~1,
             data=usedata2,times=cutoff1,
             plots="Calibration",metrics=c("Brier"),split.method="bootcv",
             B=100,N=m)$Brier$score$Brier[-1]
bs1
# 0.06160649 model1.adj CGGA
# 0.05942758 model1.CGGA

bs3 <- Score(list("3-year"=model1.CGGA),
             Surv(OS.time,OS)~1,
             data=usedata2,times=cutoff2,
             plots="Calibration",metrics=c("Brier"),split.method="bootcv",
             B=100,N=m)$Brier$score$Brier[-1]

bs3
# 0.1388182 model1.adj CGGA
# 0.142122  model1.CGGA

bs5 <- Score(list("5-year"=model1.CGGA),
             Surv(OS.time,OS)~1,
             data=usedata2,times=cutoff3,
             plots="Calibration",metrics=c("Brier"),split.method="bootcv",
             B=100,N=m)$Brier$score$Brier[-1]
bs5
# 0.1816846 model1.adj CGGA
# 0.180023  model1.CGGA

##############survival ROC
#nobs <- NROW(usedata2)
cutoff1 <- 365 #这里就是时间节点的设置
cutoff2 <-1095
cutoff3 <-1825



Mayo4.1= survivalROC(Stime=usedata2$OS.time,##生存时间
                     status=usedata2$OS,## 终止事件    
                     marker = usedata2$lp, ## marker value    
                     predict.time = cutoff1,## 预测时间截点
                     method = 'KM')##span,NNE法的namda

str(Mayo4.1)## list结构,是每一个marker的cutoff值都计算出相应的TP,FP


Mayo4.2= survivalROC(Stime=usedata2$OS.time,##生存时间
                     status=usedata2$OS,## 终止事件    
                     marker = usedata2$lp, ## marker value    
                     predict.time = cutoff2,## 预测时间截点
                     method = 'KM')

Mayo4.3= survivalROC(Stime=usedata2$OS.time,##生存时间
                     status=usedata2$OS,## 终止事件    
                     marker = usedata2$lp, ## marker value    
                     predict.time = cutoff3,## 预测时间截点
                     method = 'KM')


plot(Mayo4.1$FP, Mayo4.1$TP, ## x=FP,y=TP
     type="l",col="#BC3C29FF", ##线条设置
     xlim=c(0,1), ylim=c(0,1), lwd=1.5,
     xlab=("1 - Specificity"), ##连接
     ylab="Sensitivity",
     main="Time dependent ROC")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色

#lines函数在原有基础上继续绘图 #E18727FF
#legend函数增加legend
lines(Mayo4.2$FP, Mayo4.2$TP, type="l",lwd=1.5,col="#0072B5FF",xlim=c(0,1), ylim=c(0,1))
lines(Mayo4.3$FP, Mayo4.3$TP, type="l",lwd=1.5,col="#E18727FF",xlim=c(0,1), ylim=c(0,1))
legend(0.45,0.3,c(paste("AUC of 1-year =",round(Mayo4.1$AUC,3)),
                  paste("AUC of 3-year =",round(Mayo4.2$AUC,3)),
                  paste("AUC of 5-year =",round(Mayo4.3$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("#BC3C29FF","#0072B5FF",'#E18727FF'),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)# 



#########  校准图 后面为了统一画风，用手算的代码，这里就注释掉了
#1 years

## 正确的外部验证校准图画法：
# 
# 
# plotCalibration(Score(list("1-year"=model1.adj),
#                       Surv(OS.time,OS)~1,
#                       data=usedata2,times=cutoff1,
#                       plots="Calibration"),method = "quantile",
#                 cens.method="local",legend = F,
#                 col = "#BC3C29FF",
#                 lwd=2,q=5)
# 
# 
# plotCalibration(Score(list("3-year"=model1.adj),
#                       Surv(OS.time,OS)~1,
#                       data=usedata2,times=cutoff2,
#                       plots="cal"),method = "quantile",
#                 cens.method="local",add=T,
#                 col = "#0072B5FF",
#                 lwd=2,q=5)
# 
# plotCalibration(Score(list("5-year"=model1.adj),
#                       Surv(OS.time,OS)~1,
#                       data=usedata2,times=cutoff3,
#                       plots="cal"),method = "quantile",
#                 cens.method="local",add=T,
#                 col = "#E18727FF",
#                 lwd=2,q=5)
# 
# #加上图例
# legend("topleft", cex = 1,bty="n",
#        legend=c("1-year",
#                 "3-year",
#                 "5-year"), 
#        col=c("#BC3C29FF","#0072B5FF",'#E18727FF'), lwd=2)

# GSE16011 烂得令人发指 
# 天哪，我要从头
#CGGA 数据集拆开后，验证不错有改进的空间。但是校准斜率一个693＞1，一个325＜1

# 想像不到 TCGA+GEO的修正结果不错
# GEO NO quantil 校准度在下降


## 2022.08.16 想要替换S0提高验证集的校准度，就只能手算。这样画出来的图长得就不一样了哇
load("./rdata/tcga_data2_2.RData")
load("./rdata/gse16011_m7ggene_surv_no2.RData")

load("./rdata/GEO_revisionset2.RData") #修正集
load("./rdata/cgga_data2_2.RData")

# 先补两个未修正模型的校准图 用的是TCGA的S0
# 再补充修正模型替换S0后的校准曲线  用的是CGGA的S0

# 我这里先拟合一下原模型
usedata1 <- tcga_data2
model1.norm <- cph(Surv(usedata1$OS.time,usedata1$OS) ~ ENSG00000104894 + 
                     ENSG00000107581 + ENSG00000128595 + ENSG00000130309 + 
                     ENSG00000215301,
                   data = usedata1,x=T,y=T,surv=T)

usedata1$lp <- predict(model1.norm, newdata=usedata1, type="lp")
# usedata2 <- usedata1


# usedata2 <- tcga_data2 #TCGA_S0 数据集修正模型的校准图
# usedata2 <- gse16011_m7ggene_surv_no2
usedata2 <- cgga_data2
# usedata2 <- GEO_revisionset
# usedata2$OS.time <- usedata2$OS.time*365

## 验证集 二选一跑结果
# 原始模型 未修正系数 中心化
# usedata2$lp <- predict(model1.norm, newdata=usedata2, type="lp")
# usedata2$lp <- usedata2$lp - mean(usedata2$lp)

# 这里是修正后的模型系数 未中心化
usedata2$lp <- -usedata2$ENSG00000104894*0.009889326 - 
  usedata2$ENSG00000107581*0.494768762 + 
  usedata2$ENSG00000128595*0.361307666 + 
  usedata2$ENSG00000130309*0.346389452 + 
  usedata2$ENSG00000215301*0.131503000


# usedata2$lp <- 0.9382395*usedata2$lp 这里处理之后校准曲线也不怎么样

# TCGA:496
# GSE16011: 106
# CGGA:408
# Revision:602
sample_size <- 408
S1_nomogram <- seq(1:sample_size)
S3_nomogram <- seq(1:sample_size)
S5_nomogram <- seq(1:sample_size)

# 这里是分组
group <- 4 #GSE16011人少，分2组；CGGA分四组,REVISION 5GROUP
survival_predicted_1_combine <- seq(1:group)
survival_predicted_3_combine <- seq(1:group)
survival_predicted_5_combine <- seq(1:group)

survival_observed_1_combine <- seq(1:group)
survival_observed_3_combine <- seq(1:group)
survival_observed_5_combine <- seq(1:group)

survival_observed_1_var_combine <- seq(1:group)
survival_observed_3_var_combine <- seq(1:group)
survival_observed_5_var_combine <- seq(1:group)



# 注意核对好时哪个模型的基线生存率哟(修正的 10.4 基线生存率 TCGA/CGGA) 
# 原始模型的基线生存率要另算(未修正的 9.2 基线生存率)
usedata2$S1_nomogram <- 0.9938915^exp(usedata2$lp)
usedata2$S3_nomogram <- 0.9756272^exp(usedata2$lp)
usedata2$S5_nomogram <- 0.9578915^exp(usedata2$lp)

# evenly divide patients into 5 groups 我这边分五个组好了，人少
# GSE16011 人更少，只能分个3组试试

usedata2$group5 <- cut(usedata2$lp,
                       quantile(usedata2$lp, seq(0,1,0.25)), right=FALSE, labels=c(1:group))
usedata2$group5[usedata2$lp==max(usedata2$lp)] <- group

survival_predicted_1 <- 0
survival_predicted_3 <- 0
survival_predicted_5 <- 0

# 1-year predicted survival
survival_predicted_1 <- aggregate(usedata2$S1_nomogram, list(usedata2$group5), mean)
survival_predicted_1_combine <- data.frame(survival_predicted_1_combine,survival_predicted_1$x)

# 3-year predicted survival
survival_predicted_3 <- aggregate(usedata2$S3_nomogram, list(usedata2$group5), mean)
survival_predicted_3_combine <- data.frame(survival_predicted_3_combine,survival_predicted_3$x)

# 5-year predicted survival
survival_predicted_5 <- aggregate(usedata2$S5_nomogram, list(usedata2$group5), mean)
survival_predicted_5_combine <- data.frame(survival_predicted_5_combine,survival_predicted_5$x)

# observed survival
survival_observed_1 <- 0
survival_observed_1_var <- 0

survival_observed_3 <- 0
survival_observed_3_var <- 0

survival_observed_5 <- 0
survival_observed_5_var <- 0

for (j in 1:group) {
  
  data_temp <- subset(usedata2,usedata2$group5==j)
  
  fit_calibration <- survfit(Surv(OS.time,OS) ~ 1, data=data_temp)
  
  survival_observed_1[j] <- min(fit_calibration$surv[fit_calibration$time <= 365])
  survival_observed_1_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 365],1))^2
  
  survival_observed_3[j] <- min(fit_calibration$surv[fit_calibration$time <= 1095])
  survival_observed_3_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 1095],1))^2
  
  survival_observed_5[j] <- min(fit_calibration$surv[fit_calibration$time <= 1825])
  survival_observed_5_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 1825],1))^2
}


survival_observed_1_combine <- data.frame(survival_observed_1_combine,survival_observed_1)
survival_observed_3_combine <- data.frame(survival_observed_3_combine,survival_observed_3)
survival_observed_5_combine <- data.frame(survival_observed_5_combine,survival_observed_5)

survival_observed_1_var_combine <- data.frame(survival_observed_1_var_combine,survival_observed_1_var)
survival_observed_3_var_combine <- data.frame(survival_observed_3_var_combine,survival_observed_3_var)
survival_observed_5_var_combine <- data.frame(survival_observed_5_var_combine,survival_observed_5_var)


# plot calibration plot

# 1-year

survival_predicted_1_final <- exp(log(survival_predicted_1_combine[,-1]))
survival_observed_1_final <- exp(log(survival_observed_1_combine[,-1]))
survival_observed_1_var_final <- (survival_observed_1_var_combine[,-1])

survival_lower1_final <- exp(log(survival_observed_1_final) - qnorm(0.975)*survival_observed_1_var_final^0.5)
survival_upper1_final <- exp(log(survival_observed_1_final) + qnorm(0.975)*survival_observed_1_var_final^0.5)
survival_comparison1 <- data.frame(survival_predicted_1_final, survival_observed_1_final,
                                   survival_lower1_final,survival_upper1_final)

survival_comparison1$survival_upper1_final <- ifelse(survival_comparison1$survival_upper1_final>1, 1,survival_comparison1$survival_upper1_final) 
# 没有survival_upper1_final这个变量啊 所以把他改成2了

survival_comparison1$underestimate <- (survival_comparison1$survival_observed_1_final-survival_comparison1$survival_predicted_1_final)/survival_comparison1$survival_observed_1_final

c1 <- ggplot(data=survival_comparison1, aes(x=survival_predicted_1_final, y=survival_observed_1_final)) +
  geom_line(size=1, colour="#76549A")+
  geom_errorbar(data=survival_comparison1, mapping=aes(x=survival_predicted_1_final, ymin=survival_lower1_final,
                                                       ymax=survival_upper1_final),
                colour="#76549A", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="#76549A")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )
c1


# 3-year

survival_predicted_3_final <- exp(log(survival_predicted_3_combine[,-1]))
survival_observed_3_final <- exp(log(survival_observed_3_combine[,-1]))
survival_observed_3_var_final <- (survival_observed_3_var_combine[,-1])

survival_lower3_final<-exp(log(survival_observed_3_final) - qnorm(0.975)*survival_observed_3_var_final^0.5)
survival_upper3_final<-exp(log(survival_observed_3_final) + qnorm(0.975)*survival_observed_3_var_final^0.5)
survival_comparison3 <- data.frame(survival_predicted_3_final, survival_observed_3_final,
                                   survival_lower3_final,survival_upper3_final)

survival_comparison3$survival_upper3_final<-ifelse(survival_comparison3$survival_upper3_final>1, 1,survival_comparison3$survival_upper3_final)
survival_comparison3$underestimate<-(survival_comparison3$survival_observed_3_final-survival_comparison3$survival_predicted_3_final)/survival_comparison3$survival_observed_3_final

c2 <- ggplot(data=survival_comparison3, aes(x=survival_predicted_3_final, y=survival_observed_3_final)) +
  geom_line(size=1, colour="#DF7861")+
  geom_errorbar(data=survival_comparison3, mapping=aes(x=survival_predicted_3_final, ymin=survival_lower3_final,
                                                       ymax=survival_upper3_final),
                colour="#DF7861", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="#DF7861")+
  xlim(0,1)+ #坐标轴改一下，从0开始
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )
c2


# 5-year
survival_predicted_5_final <- exp(log(survival_predicted_5_combine[,-1]))
survival_observed_5_final <- exp(log(survival_observed_5_combine[,-1]))
survival_observed_5_var_final <- (survival_observed_5_var_combine[,-1])

survival_lower5_final<-exp(log(survival_observed_5_final) - qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_upper5_final<-exp(log(survival_observed_5_final) + qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_comparison5 <- data.frame(survival_predicted_5_final, survival_observed_5_final,
                                   survival_lower5_final,survival_upper5_final)

survival_comparison5$survival_upper5_final<-ifelse(survival_comparison5$survival_upper5_final>1, 1,survival_comparison5$survival_upper5_final)
survival_comparison5$underestimate<-(survival_comparison5$survival_observed_5_final-survival_comparison5$survival_predicted_5_final)/survival_comparison5$survival_observed_5_final

c3 <- ggplot(data=survival_comparison5, aes(x=survival_predicted_5_final, y=survival_observed_5_final)) +
  geom_line(size=1, colour="#94B49F")+
  geom_errorbar(data=survival_comparison5, mapping=aes(x=survival_predicted_5_final, ymin=survival_lower5_final,
                                                       ymax=survival_upper5_final),
                colour="#94B49F", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="#94B49F")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.line = element_line(color='black'),
  )

c3


## 先合并数据集

survival_comparison1 <- survival_comparison1

colnames(survival_comparison1) <- c("survival_predicted_final",
                                    "survival_observed_final","survival_lower_final",
                                    "survival_upper_final","underestimate")

survival_comparison1$group <- 1


colnames(survival_comparison3) <- c("survival_predicted_final",
                                    "survival_observed_final","survival_lower_final",
                                    "survival_upper_final","underestimate")

survival_comparison3$group <- 3


colnames(survival_comparison5) <- c("survival_predicted_final",
                                    "survival_observed_final","survival_lower_final",
                                    "survival_upper_final","underestimate")

survival_comparison5$group <- 5


survival_comparison <- rbind(survival_comparison1,survival_comparison3,survival_comparison5)
survival_comparison$group <- as.factor(survival_comparison$group)

ctotal <- ggplot(data=survival_comparison,
                 aes(x=survival_predicted_final, 
                     y=survival_observed_final,
                     group=group,color=group)) +
  geom_line(size=0.95,alpha=0.95)+
  geom_errorbar(data=survival_comparison,
                mapping=aes(x=survival_predicted_final, ymin=survival_lower_final,
                            ymax=survival_upper_final),
                size=0.95,alpha=0.7, linetype=1,width = 0.035)+
  geom_point(size=2)+
  #加分组颜色
  scale_colour_manual(values = c("#76549A" ,"#DF7861","#94B49F"),
                      labels = c("1-year","3-year", "5-year"))+
  
  
  #设置坐标轴
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.line = element_line(color='black')
  )

ctotal



