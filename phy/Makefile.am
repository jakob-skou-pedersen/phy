## Source directory

SUBDIRS = test

h_sources =  PhyDef.h GrammarDef.h NewickNode.h NewickParserBase.h		\
             NewickParser.h Tree.h PhyIO.h DiscreteFactorGraph.h		\
             DfgIO.h TRCTMCTree.h						\
             PhyloTree.h PhyloTreeIO.h FactorGraph.h utils.h			\
             utilsLinAlg.h Observations.h ObservationsIO.h			\
             DfgDataWrap.h xdoubleMod.h Factors.h FactorsIO.h			\
             RateMatrix.h ama.h SeqIO.h DataIO.h optimizexx.h Emit.h		\
             Grammar.h GrammarTrain.h GrammarWrap.h GrammarIO.h EmitWrap.h	\
             EmitWrapIO.h EmitWrapPhylo.h EmitWrapDfg.h StatDist.h

ccp_sources = NewickNode.cpp NewickLexer.cpp NewickParserBase.cpp NewickParser.cpp PhyIO.cpp	\
              DiscreteFactorGraph.cpp DfgIO.cpp TRCTMCTree.cpp					\
              PhyloTree.cpp PhyloTreeIO.cpp FactorGraph.cpp					\
              utils.cpp Observations.cpp							\
              ObservationsIO.cpp DfgDataWrap.cpp						\
              Factors.cpp FactorsIO.cpp RateMatrix.cpp ama.cpp SeqIO.cpp			\
              DataIO.cpp optimizexx.cpp Grammar.cpp GrammarTrain.cpp				\
              GrammarWrap.cpp GrammarIO.cpp EmitWrap.cpp EmitWrapIO.cpp EmitWrapDfg.cpp 

library_includedir=$(includedir)/$(PACKAGE)-$(API_VERSION)/$(PACKAGE)
library_include_HEADERS = $(h_sources)

#  compiler options
AM_CPPFLAGS = -I$(top_srcdir)

lib_LTLIBRARIES= libphy.la
libphy_la_SOURCES= $(h_sources) $(ccp_sources)
#libphy_la_LIBADD= sub/libsub.la
# the library flags are necessary on mac. To do: find more elegant way of including these
libphy_la_LDFLAGS= -version-info $(LIB_VERSION) -lntl -lopt -lnewmat
