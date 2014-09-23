# Examples training a mixture of normal distributions
# test1VarData.tab: Fully observed 2D normally distributed data

# Simulated data 100 samples 
# O1 normally distributed with mean 0.2 and var 0.04
# O2 normally distributed with mean O1+0.2 and var 0.04
../../dfgTrain --emTrain test1VarData.tab --writeInfo
