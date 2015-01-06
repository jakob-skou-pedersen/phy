#Essentially binomial sample n=15
#Bg p=0.5, fg p=0.25
#Evaluate score log(L_fg/L_bg)

#Run
./scoreEval -s ./test/scoreEvalEx1/dfgSpec/ -a 1.0 -n 10000
#Try also different values of a, remember a=0.0 is naive sampling

#True Values are
MAX:  	6.08198
99.99%:	4.98336
99.9%:	3.88475
99.0%:	2.78614
95.0%:	1.68753

Normalized score mean: -2.157616


