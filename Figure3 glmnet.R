#install.packages('survival')
#install.packages('survminer')
#install.packages("glmnet")


#引用包
library(survival)
library(survminer)
library(glmnet)
coxPfilter=0.01         #单因素cox方法显著性过滤标准
setwd("D:\\original data")         #设置工作目录
rt=read.table("lncRNA_CPM.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
rt$futime=rt$futime/365

#单因素cox分析
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
	if(sd(rt[,i])<0.001){next}
	#单因素cox分析
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	#保留显著性基因
	if(coxP<coxPfilter){
	        sigGenes=c(sigGenes,i)
			outTab=rbind(outTab,
			             cbind(id=i,
			             HR=coxSummary$conf.int[,"exp(coef)"],
			             HR.95L=coxSummary$conf.int[,"lower .95"],
			             HR.95H=coxSummary$conf.int[,"upper .95"],
			             pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			             )
	}
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

###lasso筛选基因
rt=read.table("uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)       #读取文件
rt$futime[rt$futime<=0]=0.003
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)


#COX模型构建
rt=read.table("lasso.SigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")    #多因素cox根据AIC值进行过滤，如果剩余的基因数目多，可以把这行前面的#号去掉
multiCoxSum=summary(multiCox)

#输出模型参数
outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)

#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
riskOut=cbind(rt[,outCol],riskScore,risk)
riskOut=cbind(id=rownames(riskOut),riskOut)
write.table(riskOut,file="geneRisk.txt",sep="\t",quote=F,row.names=F)

#绘制森林图
pdf(file="multiforest.pdf",width = 10,height = 6,onefile = FALSE)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

