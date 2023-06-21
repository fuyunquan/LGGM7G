# 血泪教训史
# 反正把建模和验证的数据集一起去除批次效应，结果会好些
# 第二次，决定删掉GSE16011 校准不好真是烦死了

#虽然重做一次GEO尝试合并TCGA 
#但是Rank-in 文件里没有分位数标准化  又要重做了


rm(list=ls()) #删除工作历史的所有变量
ls()


setwd("E:/A_sudakeyan/LGGm7G")#指定工作目录
getwd()


##  导入数据
load("./rdata/tcag_m7g_mdeg.RData")

genenames <- tcag_m7g_mdeg[,1]


#GEO
#load("./rdata/gse16011_m7ggene_surv.RData")
load("./rdata/gse16011_m7ggene_surv_no.RData")

names <- colnames(gse16011_m7ggene_surv)
gse16011_m7ggene_surv_no <- gse16011_m7ggene_surv_no[names]
save(gse16011_m7ggene_surv_no,file = "./rdata/gse16011_m7ggene_surv_no.RData")

# 20 genes

gse16011_m7ggene_surv_no2 <- gse16011_m7ggene_surv_no[genenames]

gse16011_m7ggene_surv_no2 <- cbind(gse16011_m7ggene_surv_no[,1:3],gse16011_m7ggene_surv_no2)



#CGGA
load("./rdata/cgga_all.RData")

#TCGA 返璞归真了属于是
load("./rdata/tcga_all.RData")


#之前的结果告诉我，全部基因一起去除批次效应，去除不了
#所以决定单独把预后基因摘出来，去除批次效应


#取出预后基因

progenes <- c("ENSG00000104894","ENSG00000107581","ENSG00000128595",
              "ENSG00000130309","ENSG00000215301")
# progenes <- genenames

cgga_gene <- cgga_all[progenes]
rownames(cgga_gene) <- cgga_all$CGGA_ID

gse16011_gene <- gse16011_m7ggene_surv_no[progenes]
rownames(gse16011_gene) <- gse16011_m7ggene_surv_no$ID_REF

tcga_gene <- tcga_all[progenes]
rownames(tcga_gene) <- tcga_all$sample

#### 合并

#使用dplyr包中的bind_rows()函数按列合并，
#不要求合并字段的名称必须相同，这个函数会自己做判断
library(dplyr)
df3 <-  dplyr::bind_rows(tcga_gene,cgga_gene,gse16011_gene)


#转置
usedata1 <- as.data.frame(t(df3))



#PCA查批次效应
library(ggfortify)

# apply PCA - scale. = TRUE is highly advisable, but default is FALSE. rma data have done log2,so no scale here
#par(mfcol=c(2,2))
out_pca <- prcomp(df3)
plot(out_pca,type="l")
autoplot(out_pca,data=df3,size=0.1,label=TRUE,label.size=3)
#这回是很明显的3堆
# 2堆

library(sva)
# 注意这里对应的是样本数
batch1 <- c(rep('TCGA',nrow(tcga_all)),rep('CGGA',nrow(cgga_all)),rep('GSE16011',nrow(gse16011_gene)))

x_merge2 <- as.matrix(usedata1)#关键的一步，转换为有向量、数值或者矩阵

#4.校正其实就一步
combat_edata <- ComBat(dat = x_merge2, batch = batch1)


#再看一次批次效应
# apply PCA - scale. = TRUE is highly advisable, but default is FALSE. rma data have done log2,so no scale here
#par(mfcol=c(2,2))
out_pca <- prcomp(t(combat_edata))
plot(out_pca,type="l")  #这里为什么一直都没降下去
autoplot(out_pca,data=df3,size=0.1,label=TRUE,label.size=3)
# 这里倒是聚成了一团


#然后再拆出来
combat_edata <- as.data.frame(t(combat_edata))

combat_edata$sample_ID <- rownames(combat_edata)


##  挨个拆数据集

## TCGA
tcga_gene2 <- combat_edata[grep(pattern = "TCGA",combat_edata[,21]),] #最后一列
colnames(tcga_gene2)[21] <- "sample"
tcga_clin <- tcga_all[,1:12]
tcga_data2 <- merge(tcga_clin,tcga_gene2,by="sample")
tcga_data3 <- tcga_data2
save(tcga_data3,file = "./rdata/tcga_data3.RData") # 20个预后基因一起


# 这是GEO组内批次效应矫正(分位数标准化)后
save(tcga_data2,file = "./rdata/tcga_data2.RData")

# 这里是保留GEO不分位数标准化后的结果
save(tcga_data2,file = "./rdata/tcga_data2_2.RData")

# 这里是第二次去除批次效应，不带GEO
tcga_all2 <- tcga_data2
save(tcga_all2,file = "./rdata/tcga_all2.RData")


## CGGA
cgga_gene2 <- combat_edata[grep(pattern = "CGGA",combat_edata[,21]),]
colnames(cgga_gene2)[21] <- "CGGA_ID"
cgga_clin <- cgga_all[,1:11]
cgga_data2 <- merge(cgga_clin,cgga_gene2,by="CGGA_ID")
cgga_data3 <- cgga_data2
save(cgga_data3,file = "./rdata/cgga_data3.RData") # 20 genes


# GEO 分位数标准化
save(cgga_data2,file = "./rdata/cgga_data2.RData")

# 这里是保留GEO不分位数标准化后的结果
save(cgga_data2,file = "./rdata/cgga_data2_2.RData")

# 这里是第二次去除批次效应，不带GEO
cgga_all2 <- cgga_data2
save(cgga_all2,file = "./rdata/cgga_all2.RData")


## GSE16011
gse16011_gene2 <- combat_edata[grep(pattern = "GSM",combat_edata[,21]),]
colnames(gse16011_gene2)[21] <- "ID_REF"
gse16011_clin <- gse16011_m7ggene_surv_no[,1:3]
gse16011_m7ggene_surv_no2 <- merge(gse16011_clin,gse16011_gene2,by="ID_REF")
gse16011_m7ggene_surv_no3 <- gse16011_m7ggene_surv_no2
save(gse16011_m7ggene_surv_no3,file = "./rdata/gse16011_m7ggene_surv_no3.RData")
# 20 genes

save(gse16011_m7ggene_surv_no2,file = "./rdata/gse16011_m7ggene_surv_no2.RData")
# 这里是保留GEO不分位数标准化后的结果




###  到此为止，三个数据集去除批次效应完毕，并且分别合并临床数据后保存好了 2022.06.19

### 去掉GEO数据集，做第二次批次矫正 2022.06.20 23:48 

### 然后返回代码9.2 去做验证






