check_PROGRAMS = testRunner testParser
TESTS = $(check_PROGRAMS)

AM_CPPFLAGS = -I$(top_srcdir) 

#tmpRun_LDADD = ../libphy.la -llapack -lntl -lopt -lnewmat -lm
#tmpRun_SOURCES = tmp.cpp


# on athena, test lib has diff name. TO DO: This should be encoded in the configure script.
#testRunner_LDADD = ../libphy.la -lboost_unit_test_framework-gcc41-mt -llapack -lntl -lopt -lnewmat -lm
# on mac:
#testRunner_LDADD = ../libphy.la -lboost_unit_test_framework-xgcc40-mt -llapack -lntl -lopt -lnewmat -lm
testRunner_LDADD = ../libphy.la -lboost_unit_test_framework -llapack -lntl -lopt -lnewmat -lm

testRunner_SOURCES = testRunner.cpp PhyTestDef.cpp PhyIO_test.cpp TRCTMCTree_test.cpp DiscreteFactorGraph_test.cpp	\
	PhyloTreeIO_test.cpp Observations_test.cpp PhyDef_test.cpp							\
	RateMatrix_test.cpp Factors_test.cpp ama_test.cpp								\
	utils_test.cpp Grammar_test.cpp StatDist_test.cpp EmitWrap_test.cpp

#testRunner_SOURCES = testRunner.cpp PhyTestDef.cpp Grammar_test.cpp 

testParser_LDADD = ../libphy.la -llapack -lntl -lopt -lnewmat -lm
testParser_SOURCES = testParser.cpp

