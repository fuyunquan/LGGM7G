# 接代码9.2 修正模型后的验证
# 这里后来又重画了一次GEO数据集没有分位数标准化的图
# 总的来讲，校准度有点下降
# 想说既然基线生存率相差有点大，反正后面的临床都更新了，那不如这边也替换掉算了


rm(list=ls()) #删除工作历史的所有变量
ls()


setwd("E:/A_sudakeyan/LGGm7G")#指定工作目录
getwd()


library(glmnet)
library(survival)
library(rms)
library(MASS)
library(pec)

library(ggplot2)
library(survminer)

library(reshape)
library(data.table)
library(dplyr)

## 第一次去除批次效应 TCGA+GEO 合并数据集 重估系数 模型修正
load("./rdata/tcga_data2_2.RData")
load("./rdata/gse16011_m7ggene_surv_no2.RData")
load("./rdata/cgga_data2_2.RData") # test set

load("./rdata/GEO_revisionset2.RData") #修正集


# 修正模型在改基线生存率前，用的还是TCGA的基线生存率。
# 所以在这里先跑出TCGA的基线生存率，然后再在后面吧修正集结果当验证集步骤跑
usedata1 <- tcga_data2


#改成 年
usedata1$OS.time <- usedata1$OS.time/365

summary(usedata1$OS.time) 


# # 先算一次原始模型的
# model1.norm <- cph(Surv(usedata1$OS.time,usedata1$OS) ~ ENSG00000104894 +
#                      ENSG00000107581 + ENSG00000128595 + ENSG00000130309 +
#                      ENSG00000215301,
#                    data = usedata1,x=T,y=T,surv=T)
# 
# usedata1$lp <- predict(model1.norm, newdata=usedata1, type="lp")
# summary(usedata1$lp)

## updating model
# GEO的非分位数标准化后
usedata1$lp <- -usedata1$ENSG00000104894*0.009889326 -
  usedata1$ENSG00000107581*0.494768762 +
  usedata1$ENSG00000128595*0.361307666 +
  usedata1$ENSG00000130309*0.346389452 +
  usedata1$ENSG00000215301*0.131503000

# 手动中心化一次
mean1 <- mean(usedata1$lp) #注意这里一定要对推导集的均数重新赋值，不然后面手算验证集减的这个均数就是0了
usedata1$lp <- usedata1$lp-mean1
summary(usedata1$lp)

## 计算预测概率

#中位数做截断值
cutpoint<-median(usedata1$lp)
cutpoint

usedata1$riskgroup <- ifelse(usedata1$lp>=cutpoint,'2','1')
summary(as.factor(usedata1$riskgroup))

#对于每个t，求所有人的平均，然后再画

final_model <- cph(Surv(OS.time, OS)~lp, data=usedata1,init=1, iter.max=0,x=TRUE,y=TRUE)
basehazard<-basehaz(final_model)
#baseriskav<-as.data.frame(basehazard$hazard)
#basehazardtime<-as.data.frame(basehazard$time)



#Baseline hazard average
#baserisk<-baseriskav
#basehazardtime<-basehazardtime
baserisk<-basehazard  #建模和验证用同一个
baserisk$'Baseline Survival'<-exp(-baserisk$hazard)
colnames(baserisk)[2]<-'Year'

calipre3<-subset(usedata1, select = c("lp", "riskgroup"))

for (i in 1:496) {
  calipre3[[paste0(baserisk[i,2])]] <-exp(-baserisk[i,1])^(exp(calipre3$lp))
}



predictc3<-calipre3 %>% group_by(riskgroup) %>% summarise_all(funs(mean))
predictc3$lp<-NULL
predictcl3 <- melt(setDT(predictc3), id.vars = c("riskgroup"), variable.name = "time")
colnames(predictcl3)[3]<-'predictedsurv'
colnames(predictcl3)[2]<-'timepredict'

predictcl3$timepredict <- as.numeric(as.character(predictcl3$timepredict))


#Calculate observed survival

# group "low"
data_temp1 <- subset(usedata1, usedata1$riskgroup == 1)

fit_observe1<-survfit(Surv(OS.time,OS)~1,data=data_temp1)
fit_observe1$lower

survival_observe1<- data.frame(as.data.frame(fit_observe1$surv),
                               as.data.frame(fit_observe1$lower),
                               as.data.frame(fit_observe1$upper),
                               as.data.frame(fit_observe1$time))
colnames(survival_observe1)[1]<-'observesur'
colnames(survival_observe1)[2]<-'oblower'
colnames(survival_observe1)[3]<-'obupper'
colnames(survival_observe1)[4]<-'observetime'
survival_observe1$group<-1


#group "high"
data_temp2 <- subset(usedata1, usedata1$riskgroup == 2)

fit_observe2 <- survfit(Surv(OS.time,OS) ~ 1, data=data_temp2)

survival_observe2<- data.frame(as.data.frame(fit_observe2$surv),
                               as.data.frame(fit_observe2$lower),
                               as.data.frame(fit_observe2$upper),
                               as.data.frame(fit_observe2$time))

colnames(survival_observe2)[1]<-'observesur'
colnames(survival_observe2)[2]<-'oblower'
colnames(survival_observe2)[3]<-'obupper'
colnames(survival_observe2)[4]<-'observetime'
survival_observe2$group<-2


survival_observe<-rbind(survival_observe1, survival_observe2)

survival_observe$group<-as.factor(survival_observe$group)

predictcl3$riskgroup <- as.factor(predictcl3$riskgroup)


colnames(survival_observe)[5]<-'Group'
colnames(predictcl3)[1]<-'Group'

predictcl3$Group <- factor(predictcl3$Group,levels = c(1,2),
                           labels = c("Low","High"))

survival_observe$Group <- factor(survival_observe$Group,levels = c(1,2),
                                 labels = c("Low","High"))
# Plot calibration graph
four3<- ggplot() +
  geom_line(data=predictcl3, aes(x=timepredict, y=predictedsurv, group=Group, color=Group), linetype = "dashed", size=1)+
  geom_line(data=survival_observe, aes(x=observetime, y=observesur, group=Group, color=Group))+
  geom_ribbon(data=survival_observe, aes(ymin=oblower, ymax=obupper, x=observetime, y=observesur, group=Group,
                                         fill=Group), alpha = 0.15, show.legend=F)+
  scale_x_continuous(limits = c(0,10),name='Follow-up Year', )+
  scale_y_continuous(limits = c(0,1), name='Survival Probability')+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color='black'),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))+
  scale_fill_manual(values = c("#00A8CC","#FFA41B"))+
  scale_color_manual(values = c("#00A8CC","#FFA41B"))

four3


########## 外部验证
# load("./rdata/gse16011_m7ggene_surv_no2.RData")
# load("./rdata/cgga_data2_2.RData")

# 在这里，把修正集当验证集跑
usedata2 <- GEO_revisionset

summary(usedata2$OS.time)

#
usedata2$OS.time<-usedata2$OS.time/365 
summary(usedata2$OS.time)

# # 原模型的
# usedata2$lp <- predict(model1.norm, newdata=usedata2, type="lp")
# 


#GEO的非分位数标准化后 修正后模型的
usedata2$lp <- -usedata2$ENSG00000104894*0.009889326 - 
  usedata2$ENSG00000107581*0.494768762 +
  usedata2$ENSG00000128595*0.361307666 + 
  usedata2$ENSG00000130309*0.346389452 + 
  usedata2$ENSG00000215301*0.131503000



summary(usedata2$lp)


usedata2$OS<-as.numeric(usedata2$OS)
usedata2$OS.time<-as.numeric(usedata2$OS.time)

# 2022.12.08  错误纠正：basehazard用的是推导集中中心化的结果，验证集的lp也是经历推导集的中心化。所以原始模型不用再自己手动中心化一次了
#手动中心化一次 修正模型 这里中心化用的是TCGA数据集的。因为用的是他的基线生存率
usedata2$lp <- usedata2$lp-mean1 
summary(usedata2$lp)

# 中位数做阈值分组
cutpoint <- median(usedata2$lp)
cutpoint

usedata2$riskgroup <- ifelse(usedata2$lp>=cutpoint,'2','1')

# 先跑不改基线生存率的
# 然后用CGGA的基线生存率替换，要重算
final_model <- cph(Surv(OS.time, OS)~lp, data=usedata1,init=1, iter.max=0,x=TRUE,y=TRUE)
basehazard <- basehaz(final_model)

baserisk <- basehazard  # 用CGGA的基线生存率替换
baserisk$'Baseline Survival'<- exp(-baserisk$hazard)
colnames(baserisk)[2] <- 'Year'


calipre3 <- subset(usedata2, select = c("lp", "riskgroup"))

for (i in 1:419) {
  calipre3[[paste0(baserisk[i,2])]] <-exp(-baserisk[i,1])^(exp(calipre3$lp))
}



predictc3<-calipre3 %>% group_by(riskgroup) %>% summarise_all(funs(mean))
predictc3$lp<-NULL
predictcl3 <- melt(setDT(predictc3), id.vars = c("riskgroup"), variable.name = "time")
colnames(predictcl3)[3]<-'predictedsurv'
colnames(predictcl3)[2]<-'timepredict'

predictcl3$timepredict <- as.numeric(as.character(predictcl3$timepredict))


#Calculate observed survival

# group "low"
data_temp1 <- subset(usedata2, usedata2$riskgroup == 1)

fit_observe1<-survfit(Surv(OS.time,OS)~1,data=data_temp1)
fit_observe1$upper

survival_observe1<- data.frame(as.data.frame(fit_observe1$surv),
                               as.data.frame(fit_observe1$lower),
                               as.data.frame(fit_observe1$upper),
                               as.data.frame(fit_observe1$time))
colnames(survival_observe1)[1]<-'observesur'
colnames(survival_observe1)[2]<-'oblower'
colnames(survival_observe1)[3]<-'obupper'
colnames(survival_observe1)[4]<-'observetime'
survival_observe1$group<-1


#group "high"
data_temp2 <- subset(usedata2, usedata2$riskgroup == 2)

fit_observe2 <- survfit(Surv(OS.time,OS) ~ 1, data=data_temp2)

survival_observe2<- data.frame(as.data.frame(fit_observe2$surv),
                               as.data.frame(fit_observe2$lower),
                               as.data.frame(fit_observe2$upper),
                               as.data.frame(fit_observe2$time))

colnames(survival_observe2)[1]<-'observesur'
colnames(survival_observe2)[2]<-'oblower'
colnames(survival_observe2)[3]<-'obupper'
colnames(survival_observe2)[4]<-'observetime'
survival_observe2$group<-2


survival_observe<-rbind(survival_observe1, survival_observe2)

survival_observe$group<-as.factor(survival_observe$group)

predictcl3$riskgroup <- as.factor(predictcl3$riskgroup)


colnames(survival_observe)[5]<-'Group'
colnames(predictcl3)[1]<-'Group'

predictcl3$Group <- factor(predictcl3$Group,levels = c(1,2),
                           labels = c("Low","High"))

survival_observe$Group <- factor(survival_observe$Group,levels = c(1,2),
                                 labels = c("Low","High"))
# Plot calibration graph
four3<- ggplot() +
  geom_line(data=predictcl3, aes(x=timepredict, y=predictedsurv, group=Group, color=Group), linetype = "dashed", size=1)+
  geom_line(data=survival_observe, aes(x=observetime, y=observesur, group=Group, color=Group))+
  geom_ribbon(data=survival_observe, aes(ymin=oblower, ymax=obupper, x=observetime, y=observesur, group=Group,
                                         fill=Group), alpha = 0.15, show.legend=F)+
  scale_x_continuous(limits = c(0,10),name='Follow-up Year', )+
  scale_y_continuous(limits = c(0,1), name='Survival Probability')+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color='black'),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))+
  scale_fill_manual(values = c("#00A8CC","#FFA41B"))+
  scale_color_manual(values = c("#00A8CC","#FFA41B"))

four3
