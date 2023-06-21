rm(list=ls()) #删除工作历史的所有变量
ls()


setwd("E:/A_sudakeyan/LGGm7G")#指定工作目录
getwd()


library(glmnet)
library(survival)
library(rms)
library(MASS)

#查看数据集缺失情况
library(VIM)
library(mde)
library(tidyverse)

library(mice)

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
library(compareC)

#C指数置信区间
library(survcomp)

# 载入数据
#load("./rdata/tcga_all.RData")


#合并GEO的数据集后少了几个基因，但是筛出来的结果不影响
#load("./rdata/tcga_data1.RData")

# 先导入最后确定下来的数据集，和GEO去除批次效应的那个

load("./rdata/tcga_data2_2.RData")

####################################

usedata1 <- tcga_data2


#手算gene lp
usedata1$lp_gene <- -usedata1$ENSG00000104894*0.009889326 - 
  usedata1$ENSG00000107581*0.494768762 + 
  usedata1$ENSG00000128595*0.361307666 + 
  usedata1$ENSG00000130309*0.346389452 + 
  usedata1$ENSG00000215301*0.131503000


dev.new()
aggr(usedata1[,c(5:12)],numbers=TRUE,
     cex.axis = 1,cex.numbers=1, las=2,
     ylab=c("Histogram of missing data","Pattern")) 
#numbers 和 ylab属于函数plot(  )中，numbers=	TRUE用于指定图形显示相关数据，ylab指定图形的纵坐标名称



# 查看 数据集的缺失情况
# library(tidyverse)
# library(mde)
usedata1[,3:12] %>% na_summary()
# 最高7%  那就插补十次好了

usedata1[,3:12] %>% 
  na_summary() %>% 
  ggplot(aes(x=variable,y = percent_complete, fill = variable)) + 
  geom_col() + 
  theme_bw(base_size = 18)



###插补

#要插补的数据另加一个新的数据集
impu_data <- usedata1[,c(3:12,18)]

#改名字，不能以数字开头，里面不能有反斜杠，不然插补都会报错
colnames(impu_data)[8] <- "IDH_mutation_status"
colnames(impu_data)[9] <- "X1p19q_codeletion_status"
rownames(impu_data) <- usedata1$sample
#以累积危险率(风险值)作为生存分析多重插补时间的结局变量
impu_data$haz_os <- nelsonaalen(impu_data,OS.time,OS)

#以OS作为结局

#去掉几个不需要被利用其补充缺失值信息的变量,因为有了风险值，所以不要生存时间了
pred <- quickpred(impu_data,exclude = c("OS.time"))

#开始插补imputation

#插补了10个数据集
imp <- mice(impu_data,m=10,seed=1,pred=pred)


#考虑到如果并行分析，计算在原始模型中平均S0后，模型的展示太过复杂，还是堆叠吧
#堆叠数据集
usedata1_imp <- complete(imp,action =1)
for (i in 2:10){  #从2开始
  usedata1_imp <- rbind(usedata1_imp,complete(imp,action =i))
}

#记得加权1/10
usedata1_imp$w <- 1/10 

#给基因lp改个名字，不然后面算模型的lp重名了
#colnames(usedata1_imp)[11] <- "lp_gene"

tcga_stack <- usedata1_imp

save(tcga_stack,file = "./rdata/tcga_stack.RData")



## 存了个列表
n_impu <- 10

tcga_impu <- vector(n_impu,mode="list")

for (i in 1:10) {
  tcga_impu[[i]] <- mice::complete(imp, i)
}

#保存数据
save(tcga_impu,file = "./rdata/tcga_impu.RData")




########  开始变量筛选
load("./rdata/tcga_stack.RData")

usedata1_imp <- tcga_stack
model1_imp <- cph(Surv(usedata1_imp$OS.time, usedata1_imp$OS) ~ Radiation_Therapy +
                    Grade + Age + Gender + Pharmaceutical_Therapy + IDH_mutation_status +
                    X1p19q_codeletion_status + MGMT.promoter.status + lp_gene,
                  data = usedata1_imp, x = T, y = T, surv = T,
                  weights = w)

str(usedata1_imp)

#select model
model1_imp_step<-stepAIC(model1_imp,direction = "backward") # 根据文献，还是后退法好，因为会优先考虑到全模型

model1_imp_step <- cph(formula = Surv(OS.time, OS) ~ Grade + 
      Age + Gender + IDH_mutation_status + 
        X1p19q_codeletion_status + lp_gene, 
      data = usedata1_imp, 
      weights = w, x = T, y = T, surv = T)
# 不过后退法筛出来的变量也没区别
model1_imp_step$coefficients

# 嘿我就不信这个邪了，我要去用lasso (代码10.0 lasso筛变量)
# ……还是不变的六个

# [1] "Grade"             "Age"               "Gender"           
# [4] "IDH_mutation_status"        "X1p19q_codeletion_status" "lp_gene" 

#c-index.bootstrap
set.seed(1234) # 设个随机数种子，免得bootstrap结果老是不一样
model1_c <- validate(model1_imp_step,B=1000,dxy = T)

(model1_c[1,1]+1)/2
# 0.8509231
# 0.8444222 2022.08.02 虽然很郁闷，但是lasso结果也是如此
# 0.8486942 2022.08.21 发现lp_gene算错后重跑的，结果C指数也没升

(model1_c[1,5]+1)/2
# 0.8501702
# 0.8438163  lasso
# 0.8481161  2022.08.21 发现lp_gene算错后重跑的，结果C指数也没升


model1_c[3,3] # 0.9938204

#计算置信区间
# 
# cindex <- concordance.index(predict(model1_imp_step),
#                             surv.time = usedata1_imp$OS.time, 
#                             surv.event = usedata1_imp$OS,
#                             method = "noether")
# cindex$c.index; cindex$lower; cindex$upper
#[1] 0.8509111 [1] 0.8407392 [1] 0.8610831
#[1] 0.8444031 [1] 0.8339304 [1] 0.8548759 lasso 我要emo好久了
#[1] 0.8487170 [1] 0.8385001 [1] 0.8589339 2022.08.21 发现lp_gene算错后重跑的

#内部验证的C指数与原始表现相差不多，值得庆幸

# 这里这个置信区间要重新算啊,因为堆叠数据集的SE会偏小

load("./rdata/tcga_impu.RData")
# 备用个新数据集
data_impu <- tcga_impu


#插补的十次 / 30次 GSE16011
n_impu <- 10

# C statistics
c_harrell_apparent <- 0 #粗表现
c_harrell_apparent_var<-0 #方差

for (i in 1:n_impu) {
  
  data_impu[[i]]$lp <- predict(model1_imp_step, newdata=data_impu[[i]], type="lp")   # linear predictor
  
  final_model <- cph(Surv(OS.time, OS)~lp, data=data_impu[[i]],x=TRUE,y=TRUE) ##
  c_harrell_apparent[i] <- (final_model$stats["Dxy"]+1)/2
  c_harrell_apparent_var[i] <- vardiffC(data_impu[[i]]$OS.time, data_impu[[i]]$OS, data_impu[[i]]$lp, data_impu[[i]]$lp)$est.varCxy
}


# C statistics and its 95%CI
C_statistics_final <- mean(c_harrell_apparent)
C_statistics_var_final <- mean(c_harrell_apparent_var) + (1+1/n_impu)*var(c_harrell_apparent_var)
C_statistics_final-qnorm(0.975)*C_statistics_var_final^0.5 # lower limit of 95%
C_statistics_final+qnorm(0.975)*C_statistics_var_final^0.5 # upper limit of 95%


C_statistics_final
# 0.8486942(0.8165096-0.8808789)

##########模型粗表现(在推导集中的表现)
####ROC

nobs <- NROW(usedata1_imp)
cutoff1 <- 365 #这里就是时间节点的设置
cutoff2 <-1095
cutoff3 <-1825

#这里是算最终模型的线性预测值
usedata1_imp$lp <- predict(model1_imp_step, newdata=usedata1_imp, type="lp")   # linear predictor

#1-year
Mayo4.1= survivalROC(Stime=usedata1_imp$OS.time,##生存时间
                     status=usedata1_imp$OS,## 终止事件    
                     marker = usedata1_imp$lp, ## marker value    
                     predict.time = cutoff1,## 预测时间截点
                     method = 'KM')##span,NNE法的namda

str(Mayo4.1)## list结构,是每一个marker的cutoff值都计算出相应的TP,FP

#3-year
Mayo4.2= survivalROC(Stime=usedata1_imp$OS.time,##生存时间
                     status=usedata1_imp$OS,## 终止事件    
                     marker = usedata1_imp$lp, ## marker value    
                     predict.time = cutoff2,## 预测时间截点
                     method = 'KM')

#5-year
Mayo4.3= survivalROC(Stime=usedata1_imp$OS.time,##生存时间
                     status=usedata1_imp$OS,## 终止事件    
                     marker = usedata1_imp$lp, ## marker value    
                     predict.time = cutoff3,## 预测时间截点
                     method = 'KM')

#画图
plot(Mayo4.1$FP, Mayo4.1$TP, ## x=FP,y=TP
     type="l",col="#BC3C29FF", ##线条设置
     xlim=c(0,1), ylim=c(0,1), lwd = 1.5,
     xlab=("1 - Specificity"), ##连接
     ylab="Sensitivity",
     main="Time dependent ROC")## \n换行符
abline(0,1,col="gray",lty=2, lwd = 1.5,)##线条颜色

#lines函数在原有基础上继续绘图 #E18727FF
#legend函数增加legend
lines(Mayo4.2$FP, Mayo4.2$TP, type="l", lwd = 1.5,col="#0072B5FF",xlim=c(0,1), ylim=c(0,1))
lines(Mayo4.3$FP, Mayo4.3$TP, type="l", lwd = 1.5,col="#E18727FF",xlim=c(0,1), ylim=c(0,1))
legend(0.45,0.3,c(paste("AUC of 1-year =",round(Mayo4.1$AUC,3)),
                  paste("AUC of 3-year =",round(Mayo4.2$AUC,3)),
                  paste("AUC of 5-year =",round(Mayo4.3$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("#BC3C29FF","#0072B5FF",'#E18727FF'),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)# 


########  模型校准图  堆叠的数据计算出来的置信区间会偏小，最好手算 ##############

# #1 years
# model1.1_year<- cph(formula = Surv(usedata1_imp$OS.time, usedata1_imp$OS) ~ lp,
#                     data = usedata1_imp, 
#                     weights = w, x = T, y = T, surv = T,
#                     time.inc=cutoff1)
# 
# p1<- calibrate(model1.1_year,#模型名称
#                cmethod='KM',
#                method='boot',#检测方法
#                u=cutoff1,#评估的时间，注：一定要与模型的时间一致
#                m=1240, #每次抽样的样本量，
#                B=500)#抽样次数
# #注，m值的确定：m=数据总数/3-4,即你想让最终的校准曲线有3个点，那就是m=数据总数/3
# #B值一般1000，电脑配置不好可以选500,300,100等
# 
# plot(p1,
#      add=F,#增加第二条线
#      conf.int=T,#95%CI
#      subtitles = F ,#关闭副标题
#      cex.subtitles=0.85, #副标题大小
#      lwd=1.2,#95%CI粗细
#      lty=1,#95%CI实线，2=虚线
#      errbar.col="#76549A",#95%CI颜色
#      xlim=c(0.0,1),#x轴范围
#      ylim=c(0.0,1),
#      xlab="Predicted(%)",
#      ylab="Observed(%)",
#      col="#76549A",#曲线颜色
#      bty="l",#L表示只画左边和底部，象形字(看字知意)
#      par.corrected=list(col="#76549A", lty=1, lwd=1.2, pch=4))
# lines(p1[,c('mean.predicted',"KM")], 
#       type = 'b', #连线的类型，可以是"p","b","o"
#       lwd = 1.5, #连线的粗细  要的就是这个
#       pch = 16, #点的形状，可以是0-20
#       col = c("#76549A")) #连线的颜色
# 
# 
# #3years
# model1.3_year<- cph(formula = Surv(usedata1_imp$OS.time, usedata1_imp$OS) ~lp,
#                     data = usedata1_imp, 
#                     weights = w, x = T, y = T, surv = T,
#                     time.inc=cutoff2)
# 
# p2<- calibrate(model1.3_year,
#                cmethod='KM',
#                method='boot',
#                u=cutoff2,
#                m=1240, 
#                B=500)
# 
# plot(p2,
#      add=T,
#      conf.int=T,
#      subtitles = F,
#      #cex.subtitles=0.8, 
#      lwd=1.2,
#      lty=1,
#      errbar.col="#DF7861",
#      xlim=c(0.0,1),
#      ylim=c(0.0,1),
#      col="#DF7861",
#      par.corrected=list(col="#DF7861", lty=1, lwd=1.2, pch=4))
# 
# lines(p2[,c('mean.predicted',"KM")], 
#       type = 'b', #连线的类型，可以是"p","b","o"
#       lwd = 1.5, #连线的粗细  要的就是这个
#       pch = 16, #点的形状，可以是0-20
#       col = c("#DF7861")) #连线的颜色
# 
# 
# #5 years
# model1.5_year<- cph(formula = Surv(usedata1_imp$OS.time, usedata1_imp$OS) ~lp,
#                     data = usedata1_imp, 
#                     weights = w, x = T, y = T, surv = T,
#                     time.inc=cutoff3)
# 
# p3<- calibrate(model1.5_year,
#                cmethod='KM',
#                method='boot',
#                u=cutoff3,
#                m=1240, 
#                B=500)
# 
# plot(p3,
#      add=T,
#      conf.int=T,
#      subtitles = F,
#      #cex.subtitles=0.8, 
#      lwd=1.2,
#      lty=1,
#      errbar.col="#94B49F",
#      xlim=c(0.0,1),
#      ylim=c(0.0,1),col="#94B49F",
#      par.corrected=list(col="#94B49F", lty=1, lwd=1.2, pch=4))
# 
# lines(p3[,c('mean.predicted',"KM")], 
#       type = 'b', #连线的类型，可以是"p","b","o"
#       lwd = 1.5, #连线的粗细  要的就是这个
#       pch = 16, #点的形状，可以是0-20
#       col = c("#94B49F")) #连线的颜色
# 
# #加上图例
# legend("bottomright", cex = 1,
#        legend=c("1-year","3-year", "5-year"),
#        col=c("#76549A","#DF7861", "#94B49F"), lwd=2,bty="n")
# #调整对角线
# abline(0,1,lty=3,lwd=1,col="grey")



#### ！！！忘记个重要事情，堆叠数据集的置信区间会偏小，换句话说就是不准
### 要重新手画一遍
## Calibration plot

## 这里的插补数据集用列表的
load("./rdata/tcga_impu.RData")
data_impu <- tcga_impu

n_impu <- 10
sample_size <- 496
S1_nomogram <- seq(1:sample_size)
S3_nomogram <- seq(1:sample_size)
S5_nomogram <- seq(1:sample_size)

# 这里是分组
group <- 4
survival_predicted_1_combine <- seq(1:group)
survival_predicted_3_combine <- seq(1:group)
survival_predicted_5_combine <- seq(1:group)

survival_observed_1_combine <- seq(1:group)
survival_observed_3_combine <- seq(1:group)
survival_observed_5_combine <- seq(1:group)

survival_observed_1_var_combine <- seq(1:group)
survival_observed_3_var_combine <- seq(1:group)
survival_observed_5_var_combine <- seq(1:group)


for (i in 1:n_impu) {
  
  data_impu[[i]]$lp <- predict(model1_imp_step, newdata=data_impu[[i]], type="lp")   # linear predictor
  
  # calculate predicted survival probability 
  # use the Baseline survival probability of the prediction model 要先去算一下基线生存率
  data_impu[[i]]$S1_nomogram <- 0.9729659^exp(data_impu[[i]]$lp)
  data_impu[[i]]$S3_nomogram <- 0.8332288^exp(data_impu[[i]]$lp)
  data_impu[[i]]$S5_nomogram <- 0.6481703^exp(data_impu[[i]]$lp)
  
  # evenly divide patients into 4 groups 我这边分4个组好了，人少
  
  data_impu[[i]]$group5<-cut(data_impu[[i]]$lp, quantile(data_impu[[i]]$lp, seq(0,1,0.25)), right=FALSE, labels=c(1:group))
  data_impu[[i]]$group5[data_impu[[i]]$lp==max(data_impu[[i]]$lp)] <- group
  
  survival_predicted_1 <- 0
  survival_predicted_3 <- 0
  survival_predicted_5 <- 0
  
  # 1-year predicted survival
  survival_predicted_1 <- aggregate(data_impu[[i]]$S1_nomogram, list(data_impu[[i]]$group5), mean)
  survival_predicted_1_combine <- data.frame(survival_predicted_1_combine,survival_predicted_1$x)
  
  # 3-year predicted survival
  survival_predicted_3 <- aggregate(data_impu[[i]]$S3_nomogram, list(data_impu[[i]]$group5), mean)
  survival_predicted_3_combine <- data.frame(survival_predicted_3_combine,survival_predicted_3$x)
  
  # 5-year predicted survival
  survival_predicted_5 <- aggregate(data_impu[[i]]$S5_nomogram, list(data_impu[[i]]$group5), mean)
  survival_predicted_5_combine <- data.frame(survival_predicted_5_combine,survival_predicted_5$x)
  
  # observed survival
  survival_observed_1 <- 0
  survival_observed_1_var <- 0
  
  survival_observed_3 <- 0
  survival_observed_3_var <- 0
  
  survival_observed_5 <- 0
  survival_observed_5_var <- 0
  
  for (j in 1:group) {
    
    data_temp <- subset(data_impu[[i]],data_impu[[i]]$group5==j)
    
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
  
}

# plot calibration plot

# 1-year

survival_predicted_1_final <- exp(rowMeans(log(survival_predicted_1_combine[,-1])))
survival_observed_1_final <- exp(rowMeans(log(survival_observed_1_combine[,-1])))
survival_observed_1_var_final <- rowMeans(survival_observed_1_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_1_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)
# 注意这里，多重插补之后的方差计算是不一样的，要根据rubin法则重新算

survival_lower2_final <- exp(log(survival_observed_1_final) - qnorm(0.975)*survival_observed_1_var_final^0.5)
survival_upper2_final <- exp(log(survival_observed_1_final) + qnorm(0.975)*survival_observed_1_var_final^0.5)
survival_comparison2 <- data.frame(survival_predicted_1_final, survival_observed_1_final,
                                   survival_lower2_final,survival_upper2_final)

survival_comparison2$survival_upper2_final <- ifelse(survival_comparison2$survival_upper2_final>1, 1,survival_comparison2$survival_upper2_final) 
# 没有survival_upper1_final这个变量啊 所以把他改成2了

survival_comparison2$underestimate <- (survival_comparison2$survival_observed_1_final-survival_comparison2$survival_predicted_1_final)/survival_comparison2$survival_observed_1_final

c1 <- ggplot(data=survival_comparison2, aes(x=survival_predicted_1_final, y=survival_observed_1_final)) +
  geom_line(size=1, colour="#A574B0")+
  geom_errorbar(data=survival_comparison2, mapping=aes(x=survival_predicted_1_final, ymin=survival_lower2_final,
                                                       ymax=survival_upper2_final),
                colour="#A574B0", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="#A574B0")+
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

survival_predicted_3_final <- exp(rowMeans(log(survival_predicted_3_combine[,-1])))
survival_observed_3_final <- exp(rowMeans(log(survival_observed_3_combine[,-1])))
survival_observed_3_var_final <- rowMeans(survival_observed_3_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_3_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)
# 注意这里，多重插补之后的方差计算是不一样的，要根据rubin法则重新算

survival_lower3_final<-exp(log(survival_observed_3_final) - qnorm(0.975)*survival_observed_3_var_final^0.5)
survival_upper3_final<-exp(log(survival_observed_3_final) + qnorm(0.975)*survival_observed_3_var_final^0.5)
survival_comparison3 <- data.frame(survival_predicted_3_final, survival_observed_3_final,
                                   survival_lower3_final,survival_upper3_final)

survival_comparison3$survival_upper3_final<-ifelse(survival_comparison3$survival_upper3_final>1, 1,survival_comparison3$survival_upper3_final)
survival_comparison3$underestimate<-(survival_comparison3$survival_observed_3_final-survival_comparison3$survival_predicted_3_final)/survival_comparison3$survival_observed_3_final

c2 <- ggplot(data=survival_comparison3, aes(x=survival_predicted_3_final, y=survival_observed_3_final)) +
  geom_line(size=1, colour="#ffa631")+
  geom_errorbar(data=survival_comparison3, mapping=aes(x=survival_predicted_3_final, ymin=survival_lower3_final,
                                                       ymax=survival_upper3_final),
                colour="#ffa631", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="#ffa631")+
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
survival_predicted_5_final <- exp(rowMeans(log(survival_predicted_5_combine[,-1])))
survival_observed_5_final <- exp(rowMeans(log(survival_observed_5_combine[,-1])))
survival_observed_5_var_final <- rowMeans(survival_observed_5_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_5_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)
# 注意这里，多重插补之后的方差计算是不一样的，要根据rubin法则重新算

survival_lower5_final<-exp(log(survival_observed_5_final) - qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_upper5_final<-exp(log(survival_observed_5_final) + qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_comparison5 <- data.frame(survival_predicted_5_final, survival_observed_5_final,
                                   survival_lower5_final,survival_upper5_final)

survival_comparison5$survival_upper5_final<-ifelse(survival_comparison5$survival_upper5_final>1, 1,survival_comparison5$survival_upper5_final)
survival_comparison5$underestimate<-(survival_comparison5$survival_observed_5_final-survival_comparison5$survival_predicted_5_final)/survival_comparison5$survival_observed_5_final

c3 <- ggplot(data=survival_comparison5, aes(x=survival_predicted_5_final, y=survival_observed_5_final)) +
  geom_line(size=1, colour="#66B884")+
  geom_errorbar(data=survival_comparison5, mapping=aes(x=survival_predicted_5_final, ymin=survival_lower5_final,
                                                       ymax=survival_upper5_final),
                colour="#66B884", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="#66B884")+
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


# 想办法改一改，把三条线叠加到一张图上

## 先合并数据集

survival_comparison1 <- survival_comparison2

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
                size=0.95,alpha=0.7, linetype=1,width = 0.025)+
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



##############################
##############################
### 外部验证 

#会有点麻烦，涉及外部验证的插补，合并结果。
#先另开代码把插补搞定，还有改变量名

#关于多重插补后我该用堆叠还是并行分析，有待思考
#在外部验证应该用并行分析，但我没找到文献可咋办
#内外部数据集做法不一样，被问起来不好解释啊真头痛
#所以只能都用并行分析
#……好麻烦啊
# 想了想还是只在验证集中用并行分析吧
# 开发集太————————麻烦了

## 外部验证的分析另开代码写
