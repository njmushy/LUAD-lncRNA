######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("survival")
#install.packages("survminer")

setwd("D:\\original data")              #???ù???Ŀ
library(survival)
library("survminer")
rt=read.table("sample_risk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

#????????????
pdf(file="survival.pdf",onefile = FALSE,
       width = 6,             #ͼƬ?Ŀ???
       height =5)             #ͼƬ?ĸ߶?
ggsurvplot(fit, 
           data=rt,
           conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 2,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()

summary(fit)    #?�???????????

######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056
