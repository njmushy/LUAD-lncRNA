######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

install.packages("survivalROC")

library(survivalROC)
setwd("D:\\original data")      #???ù???Ŀ
rt=read.table("sample_risk.txt",header=T,sep="\t",check.names=F,row.names=1)    #??ȡcox?ع??????ļ?
rt$futime=rt$futime/365
rocCol=rainbow(ncol(rt)-2)
aucText=c()

#????risk score??ROC????
pdf(file="multiROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
  xlab="False positive rate", ylab="True positive rate",
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#?????????ٴ???״??ROC????
j=1
for(i in colnames(rt[,3:(ncol(rt)-1)])){
	roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")
	j=j+1
	aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
	lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056
