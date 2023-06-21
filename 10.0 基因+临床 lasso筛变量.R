# 想要尝试用LASSO筛选变量，除了堆叠的数据集外，还需要个插补数据集的列表
# 《clinical prediction model》23.5.2 Impact of Imputation
# 两条路子：1.选平均λ值 2.选统一保留的变量
# 然后再在堆叠数据集里拟合


rm(list=ls()) #删除工作历史的所有变量
ls()


setwd("E:/A_sudakeyan/LGGm7G")#指定工作目录
getwd()


library(glmnet)
library(survival)
library(rms)
library(MASS)

#查看数据集缺失情况
library(mde)
library(tidyverse)

library(VIM)
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

#C指数置信区间
library(survcomp)


# 先导入最后确定下来的数据集，和GEO去除批次效应的那个

load("./rdata/tcga_data2_2.RData")

####################################

# 存个副本做分析
usedata1 <- tcga_data2

#手算gene lp
#GEO的非分位数标准化后，系数都不一样了，校准度有点降
usedata1$lp_gene <- -usedata1$ENSG00000104894*0.009889326 - 
  usedata1$ENSG00000107581*0.494768762 + 
  usedata1$ENSG00000128595*0.361307666 + 
  usedata1$ENSG00000130309*0.346389452 + 
  usedata1$ENSG00000215301*0.131503000 

#改名字，不能以数字开头，里面不能有反斜杠，不然插补都会报错
colnames(usedata1)[11] <- "X1p19q_codeletion"

# 查看 数据集的缺失情况
#library(tidyverse)
usedata1[,3:12] %>% na_summary()
# 最高7.2580645%  那就插补十次好了



usedata1[,3:12] %>% 
  na_summary() %>% 
  ggplot(aes(x=variable,y = percent_complete, fill = variable)) + 
  geom_col() + 
  # 改字体大小
  theme_bw(base_size = 12)+ 
  #去掉X轴标签
  scale_x_discrete(labels = NULL) +
  theme(legend.key.size = unit(0.1, "inches"))


# 看随机缺失的图
windows()
marginmatrix(usedata1[,3:12]) #变量太多，不好看


#######  插补

#要插补的数据另加一个新的数据集
impu_data <- usedata1[,c(3:12,18)]

rownames(impu_data) <- usedata1$sample

#以累积危险率(风险值)作为生存分析多重插补时间的结局变量
impu_data$haz_os <- nelsonaalen(impu_data,OS.time,OS)


#去掉几个不需要被利用其补充缺失值信息的变量,因为有了风险值，所以不要生存时间了
pred <- quickpred(impu_data,exclude = c("OS.time"))


#开始插补imputation

#插补了10个数据集
imp <- mice(impu_data,m=10,seed=1,pred=pred)


# 考虑到如果并行分析，计算在原始模型中平均S0后，模型的展示太过复杂，还是堆叠吧
## 堆叠数据集
usedata1_imp <- complete(imp,action =1)
for (i in 2:10){  #从2开始
  usedata1_imp <- rbind(usedata1_imp,complete(imp,action =i))
}

#记得加权1/10
usedata1_imp$w <- 1/10 

# 保存数据
save(usedata1_imp,file = "./rdata/usedata1_imp.RData")

load("./rdata/usedata1_imp.RData")

## 但是要用每个插补后的数据集走LASSO的变量筛选

# 存了个列表
n_impu <- 10

tcga_impu <- vector(n_impu,mode="list")

for (i in 1:10) {
  tcga_impu[[i]] <- mice::complete(imp, i)
}


#保存数据
save(tcga_impu,file = "./rdata/tcga_impu.RData")


############## 开始LASSO 

## 1. 先试试那个找平均λ的

lambda2 <- 0

for (i in 1:10){
  
  # 提取单个插补的数据集
  data <- as.data.frame(tcga_impu[i])
  x <- data.matrix(data[,c(3:11)]) # 协变量
  y <- Surv(data$OS.time,data$OS) # 生存
  
  #交叉验证
  set.seed(1234)
  cv.fit <- cv.glmnet(x, y, family="cox")
  
  
  #取最小值
  lambda <- cv.fit$lambda.min
  lambda2 <- rbind(lambda2,lambda)
  
}

# 看了一眼lambda，感觉相差不大啊

lambda_mean <- mean(lambda2[-1,]) #0.01685655

save(lambda2,file = "./rdata/lambda2.RData")

load("./rdata/lambda2.RData")

## 用λ平均值筛变量

data2 <- usedata1_imp
x <- data.matrix(data2[,c(3:11)]) # 协变量
y <- Surv(data2$OS.time,data2$OS) # 生存

fit <-glmnet(x, y, family = "cox")
plot(fit, xvar = "lambda", label = TRUE)
print(fit)
abline(v=log(c(lambda_mean)),lty=2)

# 筛变量
Coefficients <- coef(fit, s = lambda_mean)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Index
Active.Coefficients
row.names(Coefficients)[Active.Index]
# [1] "Grade"                    "Age"                     
# [3] "Gender"                   "IDH_mutation_status"     
# [5] "X1p19q_codeletion_status" "lp_gene"  

# 还是这样 

# 无语了 

# 给λ画个柱状图

data <- data.frame(value=lambda2[-1,],name=c(1:10))

# Specific color for each bar? Use a well known palette
library(RColorBrewer)
coul <- brewer.pal(10, "Set3") 

barplot(height=data$value, names=data$name, 
        col=coul,
        xlab="lambda", 
        ylab="times",
        horiz=T, las=1)

data$name <- as.factor(data$name)
plt <- ggplot(data,aes(x=name,y=value,fill = name)) + 
  geom_bar(stat="identity")+
  coord_flip()+
  labs(y="lambda",x="times")+
  scale_fill_brewer(palette = 'Spectral')+
  theme_bw()+
  theme(legend.position = 'none')
  

plt
