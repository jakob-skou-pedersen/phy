#Simulate data for testing DiscContFactor in phy-library
O1 <- rnorm(100, 0.2, 0.2)
O2 <- rnorm(100, 0.2,0.2) + O1
cor(O1,O2)
df <- data.frame(NAME=1:100, O1=O1, O2=O2,stringsAsFactors = F)
write.table(df,"test1VarData.tab",sep="\t",row.names=F, quote=F)
lm(O2 ~ O1, df)