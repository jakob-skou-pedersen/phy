#Simulate data for testing DiscContFactor in phy-library
set.seed(1)
O1 <- rnorm(100, 0.2, 0.2)
O2 <- rnorm(100, 0.2,0.2) + O1
cor(O1,O2)
df <- data.frame(NAME=1:100, O1=O1, O2=O2,stringsAsFactors = F)

library(ggplot2)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(df, aes(x=O1, y=O2)) + geom_point() + theme_bw() + geom_smooth(method=lm,se=FALSE)

setwd("~/git/phyGithub/src/test/dfgTrain2DNorm/")
write.table(df,"test1VarData.tab",sep="\t",row.names=F, quote=F)
lm(O2 ~ O1, df)
