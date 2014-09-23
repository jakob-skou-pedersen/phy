#Simulate data for testing DiscContFactor in phy-library

set.seed(6)
O1 <- c(rnorm(100, 0.2, 0.2),rnorm(100, 1.1, 0.1))

df <- data.frame(NAME=1:200, O1=O1, O2=c(rep('a',100),rep('b',100)),stringsAsFactors = F)
df <- df[sample(nrow(df)),]
df$NAME <- 1:200
mean(subset(df, O2 == 'a')$O1)
mean(subset(df, O2 == 'b')$O1)
var(subset(df, O2 == 'a')$O1)
var(subset(df, O2 == 'b')$O1)

library(ggplot2)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(df, aes(x=O1,fill=O2)) + geom_density() + theme_bw() + scale_fill_manual(values=cbbPalette[c(2,4)])

setwd('~/git/phyGithub/src/test/dfgTrainMixEx/')
write.table(df,"test1VarData.tab",sep="\t",row.names=F, quote=F)

#Remove some of the labels
df2 <- df
df2[sample(1:200,100),3] <- "*"
sum(df2[,3] == 'a')
write.table(df2,"test2VarData.tab",sep="\t",row.names=F, quote=F)

#Remove all labels
df3 <- df
df3[,3] <- "*"
write.table(df3,"test3VarData.tab",sep="\t",row.names=F, quote = F)

