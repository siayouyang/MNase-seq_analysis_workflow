setwd("C:/Users/sia_c/Desktop")
mwtest<-read.csv("wilcoxon.csv", header = T)
a=na.exclude(mwtest[,5])
a = as.double(a)
b=na.exclude(mwtest[,21])
wilcox.test(x= a, y=b,alternative = 'two.sided')$p.value
t.test(x=a, y=b, alternative = "two.sided", var.equal = FALSE)$p.value
print(a)

  