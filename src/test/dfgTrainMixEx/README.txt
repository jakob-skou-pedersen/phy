# Examples training a mixture of normal distributions
# test1VarData.tab: Class labels known
# test2VarData.tab: Some class labels known
# test3VarData.tab: Unknown labels

# Simulated data 100 samples from each class 
# class a: mean 0.2 sd 0.2(var 0.04)
# class b: mean 1.1 sd 0.1(var 0.01)

../../dfgTrain --emTrain test1VarData.tab --writeInfo
