#Simulate data
O1 <- rnorm(1000,0.2,0.3)

df <- data.frame(NAME=1:1000, O1=O1, stringsAsFactors = F)

write.table(df,"~/git/phyGithub/src/test/dfgTrainNormEx/test1VarData.tab",sep="\t",row.names=F, quote=F)
