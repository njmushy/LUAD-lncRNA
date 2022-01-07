#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
#install.packages("gplots")


setwd("D:\\original data")
rt=read.table("all_checkpoint_group.txt",row.names="geneNames",sep="\t",header=T)
#加载limma包，用于校正和比较差异
library(limma)

#normalize以下为对数据的校正
pdf(file="rawBox.pdf")
#画箱线图
boxplot(rt,col = "blue",xaxt = "n",outline = F)
dev.off()
#校正
rt=normalizeBetweenArrays(as.matrix(rt),method="scale")
pdf(file="normalBox.pdf")
#画箱线图
boxplot(rt,col = "red",xaxt = "n",outline = F)
dev.off()

#differential差异分析
class <- c(rep("con",26),rep("treat",26))  #前26列对照组，后26列为处理组
design <- model.matrix(~factor(class))
colnames(design) <- c("con","treat")
#算方差
fit <- lmFit(rt,design)
#贝叶斯检验
fit2 <- eBayes(fit)
#输出基因
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000) #选择20万是为了输出所有基因
#写入表格
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#找出差异两倍以上，pvalue小于0.05，写入表格
diffLab <- allDiff[with(allDiff, ((logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05 )), ]
write.table(diffLab,file="diffExp.xls",sep="\t",quote=F)

#差异基因表达水平，用于共表达
diffExpLevel <- rt[rownames(diffLab),]
qvalue=allDiff[rownames(diffLab),]$adj.P.Val
diffExpQvalue=cbind(qvalue,diffExpLevel)
write.table(diffExpQvalue,file="diffExpLevel.xls",sep="\t",quote=F)

#heatmap热图
#后面0.00001防止出现0而报错
hmExp=log10(diffExpLevel+0.00001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",height=150,width=30)
#图边距下左上右
par(oma=c(3,3,3,5))
heatmap.2(hmMat,col='greenred',trace="none",cexCol=1)
dev.off()

#volcano火山图
pdf(file="vol.pdf")
#定义横坐标和纵坐标最大值
xMax=max(-log10(allDiff$adj.P.Val))
yMax=max(abs(allDiff$logFC))
#图形化，打出所有点，都是黑色
plot(-log10(allDiff$adj.P.Val), allDiff$logFC, xlab="-log10(adj.P.Val)",ylab="logFC",main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
#找出差异点，打成红色
diffSub=subset(allDiff, allDiff$adj.P.Val<0.05 & abs(allDiff$logFC)>1)
points(-log10(diffSub$adj.P.Val), diffSub$logFC, pch=20, col="red",cex=0.4)
#加个虚线
abline(h=0,lty=2,lwd=3)
dev.off()
