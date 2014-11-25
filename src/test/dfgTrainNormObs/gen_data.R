#Data for normobs
set.seed(1)
df <- data.frame(X1 = rnorm(10,2,sqrt(2)), X3 = rnorm(10, 0,sqrt(2)))
write.table(format(df,digits=4),file="~/git/phy/src/test/dfgTrainNormObs/varData.tab",quote = F ,sep='\t')
