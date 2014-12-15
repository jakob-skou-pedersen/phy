/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#define BOOST_TEST_DYN_LINK

// bost test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// tested code (the defines allow access to private data)
#include "phy/test/PhyTestDef.h"

#define private public
#define protected public

#include "phy/DiscreteFactorGraph.h"
#include "phy/DfgIO.h"
#include "phy/TRCTMCTree.h"
#include "phy/PhyIO.h"

#undef protected
#undef private 

// input/output not included by above
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <algorithm>

// for simulation
#include <boost/random/mersenne_twister.hpp>

// other stuff
#include <boost/foreach.hpp>
#include "boost/tuple/tuple.hpp"

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;

BOOST_AUTO_TEST_CASE(DFGNode_1) 
{
  // test variable properties
  unsigned dim = 4;
  DFGNode varNode(dim);
  
  BOOST_CHECK( varNode.dimension == dim );
  BOOST_CHECK( varNode.isFactor == false );
  BOOST_CHECK( varNode.potential.size1() == 0);
  BOOST_CHECK( varNode.potential.size2() == 0);
}


BOOST_AUTO_TEST_CASE(DFGNode_2) 
{
  // test 1D-factor properties
  unsigned size = 4;
  matrix_t potential(1, size);

  DFGNode varNode(potential);
  
  BOOST_CHECK( varNode.dimension == 1 );
  BOOST_CHECK( varNode.isFactor == true );
  BOOST_CHECK( varNode.potential.size1() == 1);
  BOOST_CHECK( varNode.potential.size2() == size);
}


BOOST_AUTO_TEST_CASE(DFGNode_3) 
{
  // test 2D-factor properties
  unsigned size = 4;
  matrix_t potential(size, size,0);

  DFGNode varNode(potential);
  
  BOOST_CHECK( varNode.dimension == 2 );
  BOOST_CHECK( varNode.isFactor == true );
  BOOST_CHECK( varNode.potential.size1() == size);
  BOOST_CHECK( varNode.potential.size2() == size);
}

DFG mkSimpleFactorGraph()
{
  // strategy: define simple factor graph and check properties
  // factor graph is a phylo tree with one root and two leaves (and one initial distribution factor).
  vector<unsigned> varDimensions;
  vector<matrix_t> facPotentials;
  vector<vector<unsigned> > varNeighbors(3);
  vector<vector<unsigned> > facNeighbors(3);
  
  // three variables
  unsigned dim = 4;
  varDimensions.push_back(dim);
  varDimensions.push_back(dim);
  varDimensions.push_back(dim);
  
  // 1D-potential
  matrix_t pot1D(1, dim);
  for (unsigned i = 0; i < dim; i++)
    pot1D(0, i) = 0.25;
  
  // 2D-potential
  matrix_t Q = jukesCantorRM();
  matrix_t pot2D(dim, dim);
  number_t t = 0.5;
  initTransitionMatrix(Q, pot2D, t);

  // three factors
  facPotentials.push_back(pot2D);  // pot2D = P(x,y) = exp(Qt)
  facPotentials.push_back(pot2D);
  facPotentials.push_back(pot1D);    // pot1D = equiFreq = [0.25, 0.25, 0.25, 0.25]
  
  // setting up the graph (tree) structure
  // this is done purely in terms of factor links

  // fac 0
  facNeighbors[0].push_back(0);  // neighbor order defines variable order for factor potentials!
  facNeighbors[0].push_back(1);

  // fac 1
  facNeighbors[1].push_back(0);  // neighbor order defines variable order for factor potentials!
  facNeighbors[1].push_back(2);

  // fac 2
  facNeighbors[2].push_back(0);
  
  return DFG(varDimensions, facPotentials, facNeighbors);
}


BOOST_AUTO_TEST_CASE(DFG_1) 
{
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)

  // check that links between nodes are defined both ways. 
  for (unsigned i = 0; i < fg.neighbors.size(); i++)
    for (unsigned j = 0; j < fg.neighbors[i].size(); j++) {
      unsigned other = fg.neighbors[i][j];
      vector<unsigned> const & nb = fg.neighbors[other];
      BOOST_CHECK( find( nb.begin(), nb.end(), i) < nb.end() ); 
    }
  
  // check factor indexing and neighbors based on known links
  vector<unsigned> nb = fg.neighbors[ fg.factors[0] ];
  BOOST_CHECK( find( nb.begin(), nb.end(), fg.variables[0] ) < nb.end() ); // fac 0 links var 0

  nb = fg.neighbors[ fg.factors[2] ];
  BOOST_CHECK( find( nb.begin(), nb.end(), fg.variables[0] ) < nb.end() ); // fac 2 links var 0

  nb = fg.neighbors[ fg.factors[1] ];
  BOOST_CHECK( not ( find( nb.begin(), nb.end(), fg.variables[1] ) < nb.end() ) ); // fac 1 does not link var 1

  // check variable indexing and neighbors based on known links
  nb = fg.neighbors[ fg.variables[0] ];
  BOOST_CHECK( find( nb.begin(), nb.end(), fg.factors[0] ) < nb.end() ); // var 0 links fac 0

  nb = fg.neighbors[ fg.variables[2] ];
  BOOST_CHECK( find( nb.begin(), nb.end(), fg.factors[1] ) < nb.end() ); // var 2 links fac 1

  nb = fg.neighbors[ fg.variables[1] ];
  BOOST_CHECK( not ( find( nb.begin(), nb.end(), fg.factors[1] ) < nb.end() ) ); // var 1 does not link fac 1

  //check potentials
  BOOST_CHECK( fg.getFactor(0).potential.size1() == 4 );
  BOOST_CHECK( fg.getFactor(0).potential.size2() == 4 );
  BOOST_CHECK( fg.getFactor(1).potential.size1() == 4 );
  BOOST_CHECK( fg.getFactor(1).potential.size2() == 4 );
  BOOST_CHECK( fg.getFactor(2).potential.size1() == 1 );
  BOOST_CHECK( fg.getFactor(2).potential.size2() == 4 );
}


stateMaskVec_t const mkStateMasks(string const & obsStr, StateMaskMap const & stateMaskMap, StateMap const & metaNucMap)
{
  vector<string> obsVec = stringToVectorOfStrings(obsStr);
  stateMaskVec_t stateMasks(obsStr.size(), NULL);
  for (unsigned i = 0; i < obsStr.size(); i++)
    stateMasks[i] = & stateMaskMap.metaState2StateMask( metaNucMap.symbol2State(obsVec[i]) );
  return stateMasks;
}


BOOST_AUTO_TEST_CASE(writeDot_1) 
{
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  string const fileName("fgSimplePhylo.dot");
  fg.writeDot(fileName);
  BOOST_TEST_MESSAGE("Wrote factor graph to dot file: fgSimplePhylo.dot");
}


BOOST_AUTO_TEST_CASE(initMessages_1) 
{
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  fg.initMessages();

  unsigned i = 0;
  unsigned j = 1;
  unsigned nb = fg.neighbors[i][j];
  unsigned k  = getIndex(fg.neighbors[nb], i);
  BOOST_CHECK(fg.inMessages_[i][j]  == & fg.outMessages_[nb][k]);
  BOOST_CHECK(fg.inMessages_[nb][k] == & fg.outMessages_[i][j]);
}


BOOST_AUTO_TEST_CASE(sumProduct_1) 
{
  // setup
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  StateMap metaNucMap = mkMetaNucStateMap();
  StateMaskMap stateMaskMap( metaNucMap ); 
  stateMaskVec_t stateMasks = mkStateMasks(string("---"), stateMaskMap, metaNucMap);
  xnumber_t p = fg.calcNormConst(stateMasks);

  BOOST_CHECK_CLOSE(toNumber(p), 1.0, EPS);

  // output norm const
  //  fg.runSumProduct(stateMasks);
  //   cout << "inMessages: " << endl;
  //   for (unsigned i = 0; i < fg.inMessages_.size(); i ++) {
  //     cout << "state: " << i << endl;
  //     for (unsigned j = 0; j < fg.inMessages_[i].size(); j ++)
  //       cout << "from: " << fg.neighbors[i][j] << * fg.inMessages_[i][j] << endl;
  //     cout << endl;
  //   }
  // 
  //   cout << "outMessages: " << endl;
  //   for (unsigned i = 0; i < fg.outMessages_.size(); i ++) {
  //     cout << "state: " << i << endl;
  //     for (unsigned j = 0; j < fg.outMessages_[i].size(); j ++) 
  //       cout << "to: " << fg.neighbors[i][j] << fg.outMessages_[i][j] << endl;
  //     cout << endl;
  //   }
 
  // output norm const
  //  cout << "DFG::calcNormConst: " << fg.calcNormConst(0, stateMasks[0], fg.inMessages_[0]) << endl;
}

BOOST_AUTO_TEST_CASE(sumProduct_2)
{
  //Test disconnected graphs
  //Read factor graph
  string const inPrefix                = "./data/dfgSpecUnconnected/";
  string const stateMapsFile           = "stateMaps.txt";
  string const factorPotentialsFile    = "factorPotentials.txt";
  string const variablesFile           = "variables.txt";
  string const factorGraphFile         = "factorGraph.txt";
  string const varDataFile             = "test1VarData.txt";

  //Read dfgInfo
  DfgInfo dfgInfo = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile, inPrefix + variablesFile, inPrefix + factorGraphFile);

  //Setup data input
  stateMaskVec_t stateMasks( dfgInfo.varNames.size() );

  dfgInfo.dfg.runSumProduct(stateMasks);
  BOOST_CHECK_CLOSE( toNumber(dfgInfo.dfg.calcNormConst(stateMasks)), 1.0, EPS );

  vector<string> varNames;
  vector<symbol_t> varVec;
  vector<unsigned> maxVarStates( dfgInfo.varNames.size() );

  varNames.push_back("O1"); varNames.push_back("O3");
  varVec.push_back("B"); varVec.push_back("B");
  vector<unsigned> varMap = mkSubsetMap( dfgInfo.varNames, varNames);
  
  dfgInfo.stateMaskMapSet.symbols2StateMasks(stateMasks, varVec, varMap);
  dfgInfo.dfg.runSumProduct(stateMasks);
  BOOST_CHECK_CLOSE( toNumber(dfgInfo.dfg.calcNormConst(stateMasks)), 0.198, EPS );
  BOOST_CHECK_CLOSE( toNumber(dfgInfo.dfg.runMaxSum(stateMasks, maxVarStates)), 0.05832, EPS);

  BOOST_CHECK( maxVarStates.at(0) == 1);
  BOOST_CHECK( maxVarStates.at(1) == 1);
  BOOST_CHECK( maxVarStates.at(2) == 1);
  BOOST_CHECK( maxVarStates.at(3) == 1);
  BOOST_CHECK( maxVarStates.at(4) == 0);
}


BOOST_AUTO_TEST_CASE(calcVariableMarginals_1) 
{
  // setup
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  StateMap metaNucMap = mkMetaNucStateMap();
  StateMaskMap stateMaskMap( metaNucMap ); 
  stateMaskVec_t stateMasks = mkStateMasks(string("--A"), stateMaskMap, metaNucMap);
  fg.runSumProduct(stateMasks);

  // calcVariableMarginals
  vector<xvector_t> variableMarginals;
  fg.calcVariableMarginals(variableMarginals, stateMasks);

  // test that marginals sum to one
  for (unsigned i = 0; i < variableMarginals.size(); i ++)
    BOOST_CHECK_CLOSE( toNumber( sum(variableMarginals[i]) ), 1.0, EPS );

  // output variableMarginals
  //   cout << "variableMarginals " << endl;
  //   for (unsigned i = 0; i < variableMarginals.size(); i ++)
  //     cout << "state: " << i << " marginals: " << variableMarginals[i] << endl;
}


BOOST_AUTO_TEST_CASE(calcFactorMarginals_1) 
{
  // setup
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  StateMap metaNucMap = mkMetaNucStateMap();
  StateMaskMap stateMaskMap( metaNucMap ); 
  stateMaskVec_t stateMasks = mkStateMasks(string("--A"), stateMaskMap, metaNucMap);
  fg.runSumProduct(stateMasks);

  // calcFactorMarginals
  vector<xmatrix_t> factorMarginals;
  fg.calcFactorMarginals(factorMarginals);

  // test that marginals sum to one
  for (unsigned i = 0; i < factorMarginals.size(); i ++)
    BOOST_CHECK_CLOSE( toNumber( sumMatrix(factorMarginals[i]) ), 1.0, EPS );

  // output calcFactorMarginals
  //   cout << endl << "factorMarginals " << endl;
  //   for (unsigned i = 0; i < factorMarginals.size(); i ++)
  //     cout << "state: " << i << " marginals: " << factorMarginals[i] << endl;
}


string const varVecToObsStr(vector<unsigned> const & varVec, StateMap const & metaNucMap)
{
  string obsStr;
  BOOST_FOREACH(unsigned state, varVec) {
    obsStr += metaNucMap.state2Symbol(state);
  }
  return obsStr;
}


BOOST_AUTO_TEST_CASE(maxSum_1) 
{
  // setup
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  StateMap metaNucMap = mkMetaNucStateMap();
  StateMaskMap stateMaskMap( metaNucMap ); 

  string inObsStr("TTC");
  stateMaskVec_t stateMasks = mkStateMasks(inObsStr, stateMaskMap, metaNucMap);
  vector<unsigned> maxVar;
  fg.initMaxVariables(maxVar);
  xnumber_t p = fg.runMaxSum(stateMasks, maxVar);
  string outObsStr = varVecToObsStr(maxVar, metaNucMap);
  BOOST_CHECK(inObsStr == outObsStr);

  BOOST_CHECK_CLOSE(toNumber(p), 0.019313169120593811, EPS);
}


vector<string> mkObsStringTable()
{
  vector<string> strVec;
  strVec.push_back("--AAACTG..");
  strVec.push_back("-AAGCCTG-A");

  return strVec;
}


BOOST_AUTO_TEST_CASE(calcNormConsMultObs_1) 
{
  // setup
  DFG dfg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  StateMap metaNucMap = mkMetaNucStateMap();
  StateMaskMapSet stateMaskMapSet(metaNucMap, 3);
  vector<string> const strVec = mkObsStringTable();
  vector<unsigned> varMap;
  varMap.push_back(1);
  varMap.push_back(2);
  long seqLength = strVec[0].size();
  unsigned ranVarCount = dfg.variables.size();
  stateMask2DVec_t stateMask2DVec(seqLength, stateMaskVec_t(ranVarCount) );

  mkStateMask2DVec(strVec, stateMask2DVec, varMap, stateMaskMapSet);
  vector<xnumber_t> result = calcNormConsMultObs(stateMask2DVec, dfg);

  BOOST_CHECK_CLOSE(toNumber( result[0] ), 1.0, EPS);
  BOOST_CHECK_CLOSE(toNumber( result[1] ), 0.25, EPS);
  BOOST_CHECK_CLOSE(toNumber( result[2] ), toNumber( result[5] ), EPS);
  BOOST_CHECK_CLOSE(toNumber( result[8] ), 1.0, EPS);
}


BOOST_AUTO_TEST_CASE(calcVarAccMarMultObs_1) 
{
  // setup
  DFG dfg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  StateMap metaNucMap = mkMetaNucStateMap();
  StateMaskMapSet stateMaskMapSet(metaNucMap, 3);
  vector<string> const strVec = mkObsStringTable();
  vector<unsigned> varMap;
  varMap.push_back(1);
  varMap.push_back(2);
  long seqLength = strVec[0].size();
  unsigned ranVarCount = dfg.variables.size();
  stateMask2DVec_t stateMask2DVec(seqLength, stateMaskVec_t(ranVarCount) );

  mkStateMask2DVec(strVec, stateMask2DVec, varMap, stateMaskMapSet);
  vector<vector_t> result = calcVarAccMarMultObs(stateMask2DVec, dfg);

  // output result
  //   BOOST_FOREACH(vector_t & v, result) {
  //     BOOST_FOREACH(number_t x, v)
  //       cout << x << " ";
  //     cout << endl;
  //   }

  BOOST_CHECK_CLOSE(sum(result[0]), 10.0, EPS);
  BOOST_CHECK_CLOSE(sum(result[1]), 10.0, EPS);
  BOOST_CHECK_CLOSE(sum(result[2]), 10.0, EPS);
}


BOOST_AUTO_TEST_CASE(readVariables_1) 
{
  string const fileName1 = "./data/dfgSpec/test1Variables.txt";
  string const fileName2 = "./output/test1Variables.out.txt";

  // setup
  map<string, string> varMap = readVariables(fileName1);

  // test
   BOOST_CHECK(varMap.size() == 3);
   BOOST_CHECK(varMap["H"]  == "twoState");
   BOOST_CHECK(varMap["O1"] == "twoState");
   BOOST_CHECK(varMap["O2"] == "twoState");
}


BOOST_AUTO_TEST_CASE(readAndWriteFactorGraph_1) 
{
  string const fileName1 = "./data/dfgSpec/test1FactorGraph.txt";
  string const fileName2 = "./output/test1FactorGraph.out.txt";

  // data structures holding factor graphs
  vector<string> varNames, varNames2;
  vector<string> facNames, facNames2;
  vector<string> potNames, potNames2; 
  vector<vector<unsigned> > facNeighbors, facNeighbors2;

   // setup
  readFactorGraph(fileName1, varNames, facNames, potNames, facNeighbors);
  writeFactorGraph(fileName2, varNames, facNames, potNames, facNeighbors);
  readFactorGraph(fileName2, varNames2, facNames2, potNames2, facNeighbors2);

   BOOST_CHECK(varNames  == varNames2);
   BOOST_CHECK(facNames  == facNames2);
   BOOST_CHECK(potNames  == potNames2);
   BOOST_CHECK(facNeighbors  == facNeighbors2);
}


BOOST_AUTO_TEST_CASE(readAndWriteDfgInfo_1) 
{
  string const inPrefix    = "./data/dfgSpec/";
  string const outPrefix    = "./output/";

  string const stateMapsFile           = "test1StateMaps.txt";
  string const factorPotentialsFile    = "test1Potentials.txt";
  string const variablesFile           = "test1Variables.txt";
  string const factorGraphFile         = "test1FactorGraph.txt";

  DfgInfo dfgInfo = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile, inPrefix + variablesFile, inPrefix + factorGraphFile);
  writeDfgInfo(dfgInfo, outPrefix + stateMapsFile, outPrefix + factorPotentialsFile, outPrefix + variablesFile, outPrefix + factorGraphFile);
  DfgInfo dfgInfo2 = readDfgInfo(outPrefix + stateMapsFile, outPrefix + factorPotentialsFile, outPrefix + variablesFile, outPrefix + factorGraphFile);

  // check that the original data structures equal the reread ones
  BOOST_CHECK(dfgInfo.varNames         == dfgInfo2.varNames);
  BOOST_CHECK(dfgInfo.facNames         == dfgInfo2.facNames);
  BOOST_CHECK(dfgInfo.dfg.factors      == dfgInfo2.dfg.factors);
  BOOST_CHECK(dfgInfo.dfg.variables    == dfgInfo2.dfg.variables);
  BOOST_CHECK(dfgInfo.dfg.neighbors    == dfgInfo2.dfg.neighbors);
}

BOOST_AUTO_TEST_CASE(readDfgInfo_2){
  string const inPrefix = "./data/dfgSpecNorm/";
  string const stateMapsFile           = "test1StateMaps.txt";
  string const factorPotentialsFile    = "test1Potentials.txt";
  string const variablesFile           = "test1Variables.txt";
  string const factorGraphFile         = "test1FactorGraph.txt";

  vector<string> varNames;
  vector<string> facNames;
  vector<string> potNames; 
  vector<vector<unsigned> > facNeighbors;

  //Check readFactorGraph
  readFactorGraph( inPrefix + factorGraphFile, varNames, facNames, potNames, facNeighbors);
  BOOST_CHECK(varNames.at(0) == "O1");
  BOOST_CHECK(varNames.size() == 1);
  BOOST_CHECK(facNames.at(0) == "O1");
  BOOST_CHECK(potNames.at(0) == "norm");
  BOOST_CHECK(potNames.size() == 1);
  BOOST_CHECK(facNeighbors.at(0).at(0) == 0);

  //Check readStateMapFile
  map<string, StateMapPtr_t> smMap = readStateMapFile( inPrefix + stateMapsFile );
  BOOST_CHECK( smMap["realMap"]->stateCount() == 202);

  //Check readFactorFile
  map<string, AbsBasFacPtr_t> facMap = readFactorFile( inPrefix + factorPotentialsFile);
  BOOST_CHECK( facMap["norm"]->type() == "discCont");
  BOOST_CHECK( facMap["norm"]->size1_ == 1);
  BOOST_CHECK( facMap["norm"]->size2_ == 202);

  //Check readVariables
  map<string, string> var2smMap = readVariables( inPrefix + variablesFile);
  BOOST_CHECK( var2smMap["O1"] == "realMap" );
   
  //Checking DfgInfo constructor and the methods it calls
  vector<AbsBasFacPtr_t> facVec( mkFacVec(potNames, facMap));
  CompositeFactorSet facSet(facVec);
  vector<StateMapPtr_t> stateMapVec( mkStateMapVec(varNames, var2smMap , smMap) );
  BOOST_CHECK( stateMapVec.size() == 1) ;
  BOOST_CHECK( stateMapVec.at(0)->name_ == "realMap");

  StateMaskMapPtr_t smm_ptr = StateMaskMapPtr_t( new StateMaskMap(*stateMapVec[0]) );
  StateMaskMapSet stateMaksMapSet( stateMapVec);
}

BOOST_AUTO_TEST_CASE(readDfgInfo_3)
{
  //Test subscribed variables
  //Test subscription factors
  string const inPrefix = "./data/dfgSpecSubscribe/";
  string const stateMapsFile           = "stateMaps.txt";
  string const factorPotentialsFile    = "factorPotentials.txt";
  string const variablesFile           = "variables.txt";
  string const factorGraphFile         = "factorGraph.txt";

  vector<string> varNames;
  vector<string> facNames;
  vector<string> potNames; 
  vector<vector<unsigned> > facNeighbors;

  //Read dfgInfo
  DfgInfo dfgInfo = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile, inPrefix + variablesFile, inPrefix + factorGraphFile);

  //Check variables
  BOOST_CHECK_EQUAL( dfgInfo.varNames.at(0), "O1");
  BOOST_CHECK_EQUAL( dfgInfo.varNames.at(1), "O2");

  //Check factors
  BOOST_CHECK( dfgInfo.facNames.at(0) == "O1Prior");
  BOOST_CHECK( dfgInfo.facNames.at(1) == "O2Prior");

  //Check subscriptions
  BOOST_CHECK( dfgInfo.subNames.at( dfgInfo.subscribedVars.at(0).at(0) ) == "N");
  BOOST_CHECK( dfgInfo.subNames.at( dfgInfo.subscribedVars.at(0).at(1) ) == "P1");
  BOOST_CHECK( dfgInfo.subNames.at( dfgInfo.subscribedVars.at(1).at(0) ) == "N");
  BOOST_CHECK( dfgInfo.subNames.at( dfgInfo.subscribedVars.at(1).at(1) ) == "P2");
  BOOST_CHECK( dfgInfo.subNames.size() == 3);

}

BOOST_AUTO_TEST_CASE(readDfgInfo_4){
  //Test subscribed variables
  //Test subscription factors
  string const inPrefix = "./data/dfgSpecMichal/";
  string const stateMapsFile           = "stateMaps.txt";
  string const factorPotentialsFile    = "factorPotentials.txt";
  string const variablesFile           = "variables.txt";
  string const factorGraphFile         = "factorGraph.txt";

  vector<string> varNames;
  vector<string> facNames;
  vector<string> potNames; 
  vector<vector<unsigned> > facNeighbors;

  //Read dfgInfo
  DfgInfo dfgInfo = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile, inPrefix + variablesFile, inPrefix + factorGraphFile);

  //Check variables(order may be irrelevant)
  BOOST_CHECK_EQUAL( dfgInfo.varNames.at(0), "prob");
  BOOST_CHECK_EQUAL( dfgInfo.varNames.at(1), "M_p");
  BOOST_CHECK_EQUAL( dfgInfo.varNames.at(2), "M_gb");

  //Check factors(order may be irrelevant)
  BOOST_CHECK( dfgInfo.facNames.at(0) == "POST_prob");
  BOOST_CHECK( dfgInfo.facNames.at(1) == "prob.M_p");
  BOOST_CHECK( dfgInfo.facNames.at(2) == "prob.M_gb");
  BOOST_CHECK( dfgInfo.facNames.at(3) == "POST_M_p");
  BOOST_CHECK( dfgInfo.facNames.at(4) == "POST_M_gb");

  //Check Subscriptions(order may irrelevant)
  BOOST_CHECK( dfgInfo.subNames.size() == 11);

}

BOOST_AUTO_TEST_CASE(VarData_1) 
{
  vector<string> varNames;
  varNames.push_back("H");
  varNames.push_back("O1");
  varNames.push_back("O2");

  string fileName = "./data/test1VarData.tab";

  VarData varData(fileName, varNames);
  string id;
  vector<string> data( varData.names().size() );
  vector< vector<string> > dataVec;


  while ( varData.next(id, data) ) {
    dataVec.push_back(data);
  }

  BOOST_CHECK(dataVec.size() == 2);
  BOOST_CHECK(dataVec[0].size() == 2);
  BOOST_CHECK(dataVec[1].size() == 2);
  BOOST_CHECK(dataVec[0][0] == "A");
}

BOOST_AUTO_TEST_CASE(DFG_components_1)
{
  string const inPrefix                = "./data/dfgSpecUnconnected/";
  string const stateMapsFile           = "stateMaps.txt";
  string const factorPotentialsFile    = "factorPotentials.txt";
  string const variablesFile           = "variables.txt";
  string const factorGraphFile         = "factorGraph.txt";

  //Read dfgInfo
  DfgInfo dfgInfo = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile, inPrefix + variablesFile, inPrefix + factorGraphFile);

  dfgInfo.dfg.initComponents();

  //Check components
  //TODO: remove implicit assumption that variables are the first nodes, followed by factors
  int HvarIdx = std::find(dfgInfo.varNames.begin(), dfgInfo.varNames.end(), "H")-dfgInfo.varNames.begin();
  int O1varIdx = std::find(dfgInfo.varNames.begin(), dfgInfo.varNames.end(), "O1")-dfgInfo.varNames.begin();
  int O2varIdx = std::find(dfgInfo.varNames.begin(), dfgInfo.varNames.end(), "O2")-dfgInfo.varNames.begin();
  int O3varIdx = std::find(dfgInfo.varNames.begin(), dfgInfo.varNames.end(), "O3")-dfgInfo.varNames.begin();
  int O4varIdx = std::find(dfgInfo.varNames.begin(), dfgInfo.varNames.end(), "O4")-dfgInfo.varNames.begin();
  
  int HfacIdx = std::find(dfgInfo.facNames.begin(), dfgInfo.facNames.end(), "H")-dfgInfo.facNames.begin();
  int HO1facIdx = std::find(dfgInfo.facNames.begin(), dfgInfo.facNames.end(), "H.O1")-dfgInfo.facNames.begin();
  int HO2facIdx = std::find(dfgInfo.facNames.begin(), dfgInfo.facNames.end(), "H.O2")-dfgInfo.facNames.begin();
  int O3facIdx = std::find(dfgInfo.facNames.begin(), dfgInfo.facNames.end(), "O3")-dfgInfo.facNames.begin();
  int O4facIdx = std::find(dfgInfo.facNames.begin(), dfgInfo.facNames.end(), "O4")-dfgInfo.facNames.begin();

  int varSize = dfgInfo.varNames.size();

  //Check same components
  BOOST_CHECK( dfgInfo.dfg.components.at(HvarIdx) == dfgInfo.dfg.components.at(HfacIdx+varSize));
  BOOST_CHECK( dfgInfo.dfg.components.at(HvarIdx) == dfgInfo.dfg.components.at(HO1facIdx+varSize));
  BOOST_CHECK( dfgInfo.dfg.components.at(HvarIdx) == dfgInfo.dfg.components.at(HO2facIdx+varSize));
  BOOST_CHECK( dfgInfo.dfg.components.at(O2varIdx) == dfgInfo.dfg.components.at(HO2facIdx+varSize));
  BOOST_CHECK( dfgInfo.dfg.components.at(O1varIdx) == dfgInfo.dfg.components.at(HO1facIdx+varSize));
  BOOST_CHECK( dfgInfo.dfg.components.at(O3varIdx) == dfgInfo.dfg.components.at(O3facIdx+varSize));
  BOOST_CHECK( dfgInfo.dfg.components.at(O4varIdx) == dfgInfo.dfg.components.at(O4facIdx+varSize));

  //Check different components
  BOOST_CHECK( dfgInfo.dfg.components.at(O3varIdx) != dfgInfo.dfg.components.at(HvarIdx) );
  BOOST_CHECK( dfgInfo.dfg.components.at(O3varIdx) != dfgInfo.dfg.components.at(O4varIdx) );
  BOOST_CHECK( dfgInfo.dfg.components.at(HvarIdx) != dfgInfo.dfg.components.at(O4varIdx) );

  //
}

BOOST_AUTO_TEST_CASE(DFG_simulation)
{
  string const inPrefix                = "./data/dfgSpecUnconnected/";
  string const stateMapsFile           = "stateMaps.txt";
  string const factorPotentialsFile    = "factorPotentials.txt";
  string const variablesFile           = "variables.txt";
  string const factorGraphFile         = "factorGraph.txt";

  //Read dfgInfo
  DfgInfo dfgInfo = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile, inPrefix + variablesFile, inPrefix + factorGraphFile);

  
  stateMaskVec_t stateMasks( dfgInfo.dfg.variables.size() );
  dfgInfo.dfg.runSumProduct(stateMasks);
  boost::mt19937 gen(std::time(0));

  vector<unsigned> acc(dfgInfo.dfg.variables.size());
  for(int i = 0; i < 100000; ++i){
    vector<unsigned> sample = dfgInfo.dfg.sample(gen);
    for(int j = 0; j < acc.size(); ++j)
      acc.at(j) += sample.at(j);
  }


  std::copy(acc.begin(), acc.end(), ostream_iterator<unsigned>(std::cout, "\t") );

  

}
