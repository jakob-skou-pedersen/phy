/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#define BOOST_TEST_DYN_LINK

// bost test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// tested code
// Note the hack to access protected data
#define private public
#define protected public

#include "phy/test/PhyTestDef.h"
#include "phy/TRCTMCTree.h"
#include "phy/PhyIO.h"
#include "phy/ama.h"
#include "phy/DiscreteFactorGraph.h"
#include "phy/DfgDataWrap.h"
#include "phy/PhyloTreeIO.h"

#undef protected
#undef private 

// input/output not included by above
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

// other stuff
#include <boost/foreach.hpp>
#include <utility>

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;


BOOST_AUTO_TEST_CASE(TECTMTTree_StateMap_1) 
{
  StateMap nucMap( mkNucStateMap() );
  unsigned n = nucMap.metaStateCount();
  for (unsigned i = 0; i < n; i++)
    BOOST_CHECK(nucMap.symbol2State( nucMap.state2Symbol(i) ) == i);
  
  BOOST_CHECK( nucMap.metaStateCount() == 4 );
  BOOST_CHECK( nucMap.stateCount() == 4 );
  BOOST_CHECK( nucMap.symbolSize() == 1 );
}


BOOST_AUTO_TEST_CASE(TECTMTTree_StateMap_2) 
{
  StateMap metaNucMap = mkMetaNucStateMap();
  unsigned n = metaNucMap.metaStateCount();
  for (unsigned i = 0; i < n; i++)
    BOOST_CHECK(metaNucMap.symbol2State( metaNucMap.state2Symbol(i) ) == i);

  BOOST_CHECK( metaNucMap.metaStateCount() == 4 + 14 );
  BOOST_CHECK( metaNucMap.stateCount() == 4 );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("-") ).size() == 4 );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("-") )[0] == metaNucMap.state2Symbol(0) );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("-") )[1] == metaNucMap.state2Symbol(1) );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("-") )[2] == metaNucMap.state2Symbol(2) );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("-") )[3] == metaNucMap.state2Symbol(3) );
}


BOOST_AUTO_TEST_CASE(TECTMTTree_StateMap_3) 
{
  StateMap metaNucMap = mkMultiMetaNucStateMap(2);
  unsigned n = metaNucMap.metaStateCount();
  for (unsigned i = 0; i < n; i++)
    BOOST_CHECK(metaNucMap.symbol2State( metaNucMap.state2Symbol(i) ) == i);

  BOOST_CHECK( metaNucMap.metaStateCount() == power(unsigned(4 + 14), unsigned(2) ) );
  BOOST_CHECK( metaNucMap.stateCount() == 16 );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("--") ).size() == 16 );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("--") )[0] == metaNucMap.state2Symbol(0) );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("--") )[1] == metaNucMap.state2Symbol(1) );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("--") )[2] == metaNucMap.state2Symbol(2) );
  BOOST_CHECK( metaNucMap.degeneracyVector( string("--") )[3] == metaNucMap.state2Symbol(3) );
}


BOOST_AUTO_TEST_CASE(TECTMTTree_StateMap_4) 
{
  StateMap nucMap = mkMultiNucStateMap(4);
  unsigned n = nucMap.metaStateCount();
  for (unsigned i = 0; i < n; i++) 
    BOOST_CHECK(nucMap.symbol2State( nucMap.state2Symbol(i) ) == i);

  BOOST_CHECK( nucMap.stateCount() == 256 );
}


BOOST_AUTO_TEST_CASE(TECTMTTree_StateMask_1) 
{
  StateMap metaNucMap = mkMetaNucStateMap();
  StateMaskMap stateMask(metaNucMap);
 
  // check that the bit vector for 'A' is correct
  ublas::vector<bool> const & v = stateMask.metaState2StateMask( metaNucMap.symbol2State("A") );
  BOOST_CHECK( v[ metaNucMap.symbol2State("A") ] == true  );
  BOOST_CHECK( v[ metaNucMap.symbol2State("C") ] == false );
  BOOST_CHECK( v[ metaNucMap.symbol2State("G") ] == false );
  BOOST_CHECK( v[ metaNucMap.symbol2State("T") ] == false );

  // check that the bit vector for 'N' is correct
  ublas::vector<bool> const & u = stateMask.metaState2StateMask( metaNucMap.symbol2State("N") );
  BOOST_CHECK( u[ metaNucMap.symbol2State("A") ] == true );
  BOOST_CHECK( u[ metaNucMap.symbol2State("C") ] == true );
  BOOST_CHECK( u[ metaNucMap.symbol2State("G") ] == true );
  BOOST_CHECK( u[ metaNucMap.symbol2State("T") ] == true );
}


BOOST_AUTO_TEST_CASE(jukesCantorRM_1) 
{
  ublas::matrix<number_t> Q = jukesCantorRM();
  BOOST_CHECK_CLOSE( Q(0,0), -1.0, 0.0000000001 );
  BOOST_CHECK_CLOSE( Q(0,1), 1 / 3.0, 0.0000000001 );
  // BOOST_TEST_MESSAGE( Q(0,0) );
}


BOOST_AUTO_TEST_CASE(diagonalizeReversibleRM_1) 
{
  ublas::matrix<number_t> Q = jukesCantorRM();
  ublas::matrix<number_t> V(4,4);
  ublas::matrix<number_t> VInv(4,4);
  ublas::vector<number_t> lambda(4);

  diagonalizeReversibleRM(Q, V, VInv, lambda);
  ublas::banded_matrix<number_t> Pi( vectorToDiagonalMatrix(lambda) );

 // BOOST_TEST_MESSAGE( lambda );
 // BOOST_TEST_MESSAGE( Q );
 // BOOST_TEST_MESSAGE( V );
 // BOOST_TEST_MESSAGE( VInv );

  ublas::matrix<number_t> tmp1 = prod(V, Pi);
  ublas::matrix<number_t> prd = prod(tmp1, VInv);

  // test that diagonalization worked
  // Test: V * Pi * VInv = Q
  unsigned dim = prd.size2();
  for (unsigned i = 0; i < dim; i++)
    for (unsigned j = 0; j < dim; j++)
      BOOST_CHECK_CLOSE( prd(i, j) , Q(i, j), EPS );

  // BOOST_TEST_MESSAGE( lambda );
}


BOOST_AUTO_TEST_CASE(diagonalizeReversibleRM_2) 
{
  // rate matrix
  StateMap metaNuc = mkMetaNucStateMap(); // mkMultiMetaNucStateMap(2);
  GTRRateMatrix rateMatrix(metaNuc);
  
  vector_t pi(4);
  pi[0] = 0.1;
  pi[1] = 0.1;
  pi[2] = 0.1;
  pi[3] = 0.7;
  rateMatrix.resetEquiFreqs(pi);

  vector_t rateParams(6);
  rateParams[0] = 1;
  rateParams[1] = 2;
  rateParams[2] = 3;
  rateParams[3] = 4;
  rateParams[4] = 5;
  rateParams[5] = 6;
  rateMatrix.resetRateParameters(rateParams);

  ublas::matrix<number_t> Q = rateMatrix.mkRateMatrix();
  ublas::matrix<number_t> V(4,4);
  ublas::matrix<number_t> VInv(4,4);
  ublas::vector<number_t> lambda(4);

  diagonalizeReversibleRM(Q, V, VInv, lambda);
  ublas::banded_matrix<number_t> Pi( vectorToDiagonalMatrix(lambda) );

 // BOOST_TEST_MESSAGE( lambda );
 // BOOST_TEST_MESSAGE( Q );
 // BOOST_TEST_MESSAGE( V );
 // BOOST_TEST_MESSAGE( VInv );

  ublas::matrix<number_t> tmp1 = prod(V, Pi);
  ublas::matrix<number_t> prd = prod(tmp1, VInv);

  // test that diagonalization worked
  // Test: V * Pi * VInv = Q
  unsigned dim = prd.size2();
  for (unsigned i = 0; i < dim; i++)
    for (unsigned j = 0; j < dim; j++)
      BOOST_CHECK_SMALL( prd(i, j) - Q(i, j), 5 * EPS ); // one gave a diff of 1.45e-10

  // BOOST_TEST_MESSAGE( lambda );
}


BOOST_AUTO_TEST_CASE(initTransitionMatrix_1) 
{
  unsigned dim = 4;
  ublas::matrix<number_t> Q = jukesCantorRM();
  ublas::matrix<number_t> V(dim, dim);
  ublas::matrix<number_t> VInv(dim,dim);
  ublas::vector<number_t> lambda(dim);

  diagonalizeReversibleRM(Q, V, VInv, lambda);

  // test init conditions
  ublas::matrix<number_t> T(dim, dim);
  number_t t = 0.0;
  initTransitionMatrix (V, VInv, lambda, T, t);
  ublas::identity_matrix<number_t> I(dim, dim);
  // BOOST_TEST_MESSAGE( T );
  for (unsigned i = 0; i < dim; i++)
    for (unsigned j = 0; j < dim; j++)
      BOOST_CHECK_SMALL( T(i, j) - I(i, j), EPS );

  // test transition matrix at time 0.75
  t = 0.75;
  initTransitionMatrix (V, VInv, lambda, T, t);
  // BOOST_TEST_MESSAGE( T );
  for (unsigned i = 0; i < dim; i++)
    for (unsigned j = 0; j < dim; j++) {
      if (i == j)
	BOOST_CHECK_CLOSE( T(i, j) , 1/4. + 3/4. * exp( 4/3. * t * Q(i,j) ), EPS );
      if (i != j)
	BOOST_CHECK_CLOSE( T(i, j) , 1/4. - 1/4. * exp( -4 * t * Q(i,j) ), EPS ); // explicit jukes-cantor formulas from Ziheng Yang's "Computational Molecular Evolution"
    }
 
  // test equilibrium conditions
  t = 1000.0;
  initTransitionMatrix (V, VInv, lambda, T, t);
  ublas::vector<number_t> pi = deriveEquiFreqForReversibleRM(Q);
  // BOOST_TEST_MESSAGE( T );
  for (unsigned i = 0; i < dim; i++)
    for (unsigned j = 0; j < dim; j++)
      BOOST_CHECK_CLOSE( T(i, j) , pi(i), EPS );
}


BOOST_AUTO_TEST_CASE(dev_1) 
{
  // phylo tree string
  string treeStr = readNewickStr("./data/vert28way.nwk");
  TreeInfo treeInfo = newickToTreeInfo(treeStr);

  // rate matrix
  StateMap metaNuc = mkMetaNucStateMap(); // mkMultiMetaNucStateMap(2);
  boost::shared_ptr<BaseTRRateMatrix> rateMatrixPtr( new GTRRateMatrix(metaNuc) );

  // dfg and factorSet
  TrctmcFactorSet factorSet = mkTrctmcFactorSet(treeInfo, rateMatrixPtr);
  DFG dfg = mkPhyloDfg(treeInfo, factorSet);

  // alg data
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );

  // stateMaskMapSet
  StateMaskMapSet stateMaskMapSet(metaNuc, dfg.variables.size());

//  for (unsigned i = 0; i < seqDataVec.size(); i++) {
//    stateMask2DVec_t stateMask2DVec = mkStateMask2DVec(seqDataVec[i], treeInfo.nodeMap, stateMaskMapSet, 1);
//    vector<xnumber_t> likVec = calcNormConsMultObs(stateMask2DVec, dfg);
//    vector<vector_t> varAccMar = calcVarAccMarMultObs(stateMask2DVec, dfg);
//
//    // max sum
//    vector<xnumber_t> maxProb;
//    vector<vector<unsigned> > maxVar;
//    calcMaxProbStatesMultObs(maxProb, maxVar, stateMask2DVec, dfg);
//
//    // debug output
//    //     // output likelihoods
//    //     cout << "ENTRY: " << seqDataVec[i].entryId << endl;
//    //     cout << "column likelihoods:" << endl;
//    //     for (unsigned i = 0; i < likVec.size(); i++)
//    //       cout << i << ": " << likVec[i] << endl;
//    //     cout << endl;
//    // 
//    //     // output var marginals
//    //     cout << "variable marginals:" << endl;
//    //     for (unsigned i = 0; i < varAccMar.size(); i++)
//    //       cout << i << ": " << varAccMar[i] << " sum: " << sum(varAccMar[i]) << endl;
//    //     cout << endl;;
//    // 
//    //     // output max probs and max states
//    //     cout << "max probs and states:" << endl;
//    //     for (unsigned i = 0; i < maxProb.size(); i++) {
//    //       cout << i << ": " << maxProb[i] << " Max states: ";
//    //       for (unsigned j = 0; j < maxVar[i].size(); j++)
//    // 	cout << metaNuc.state2Symbol(maxVar[i][j]);
//    //       cout << " root: " << metaNuc.state2Symbol( maxVar[i][treeInfo.root] );
//    //       cout << endl;
//    //     }
//    //     cout << endl;
//  }

  // test rateMatrixPtr->mkRateMatrix
  matrix_t Q = rateMatrixPtr->mkRateMatrix();
  for (unsigned i = 0; i < Q.size1(); i++) 
    BOOST_CHECK_SMALL(sum(row(Q, i)), EPS); // all rate matrix rows must sum to zero
    
  // test rateMatrixPtr->equiFreqs
  vector_t pi = rateMatrixPtr->equiFreqs();
  BOOST_CHECK_CLOSE(sum(pi), 1.0, EPS); // equiFreqs must sum to one

  // test TrctmcFactorSet
  BOOST_CHECK( treeInfo.branchMap.size() + 1 == factorSet.facCount() ); // a factor for each branch and one for the root prior
  vector<matrix_t> v = factorSet.mkFactorVec();
  for (unsigned i = 0; i < v.size(); i++) 
    for (unsigned j = 0; j < v[i].size1(); j++) 
      BOOST_CHECK_CLOSE(sum(row(v[i], j)), 1.0, EPS); // all factor rows must sum to one

//   // testing the EM algorithm
//   vector<stateMask2DVec_t> dataVec = mkStateMask2DVec(seqDataVec, treeInfo.nodeMap, stateMaskMapSet, 1);
//   vector<stateMask2DVec_t> dataVec1;
//   dataVec1.push_back(dataVec[4]);
//   dfgEm(factorSet, dataVec1, dfg, 1e-4, 1000, "output/emTestLog.txt");
//   
//   // check that optimization worked. Input had 50% A's:
//   number_t AFreq = rateMatrixPtr->equiFreqs()[ metaNuc.symbol2State("A") ];
//   BOOST_CHECK_CLOSE(AFreq, 0.50, 0.1); 
// 
//   // regression test: check that rate scale converges to previously observed value
//   BOOST_CHECK_CLOSE(factorSet.rateScale_,  0.0583, 0.1); 
// 

  // debug
  //  for (unsigned i = 0; i < stateMask2DVec.size(); i++) {
  //    for (unsigned j = 0; j < stateMask2DVec[i].size(); j++) {
  //      cout << j << " " << treeInfo.nodeMap[j] << ": " ;
  //      if (stateMask2DVec[i][j])
  //	cout << * stateMask2DVec[i][j] << endl;
  //      else
  //	cout << "NULL" << endl;
  //    }
  //    cout << endl;
  //  }
  //  

  // calc norm const
  //
  //  
  //  cout << endl << endl;
  //  vector<xvector_t> varAccMar = calcVarAccMarMultObs(stateMask2DVec, dfg);
  //  for (unsigned i = 0; i < varAccMar.size(); i++)
  //    cout << i << ": " << varAccMar[i] << " sum: " << sum(varAccMar[i]) << endl;
}


BOOST_AUTO_TEST_CASE(TrctmcFactorSetAndsumProductWithMissingData_2) 
{
  // test of specific set of parameters that behaved pathologically in optimization
  // note: bug fixed, example no longer pathological.

  // phylo tree string
  string treeStr = readNewickStr("./data/vert28way.nwk");
  TreeInfo treeInfo = newickToTreeInfo(treeStr);

  // rate matrix
  StateMap metaNuc = mkMetaNucStateMap(); // mkMultiMetaNucStateMap(2);

  vector_t pi(4);
  pi[0] = 0.8552592;
  pi[1] = 0.0482468;
  pi[2] = 0.0482471;
  pi[3] = 0.0482469;

  vector_t rateParams(6);
  rateParams[0] = 0.989851;
  rateParams[1] = 0.989877;
  rateParams[2] = 0.989859;
  rateParams[3] = 1.01024;
  rateParams[4] = 1.01026;
  rateParams[5] = 1.01025;

  boost::shared_ptr<BaseTRRateMatrix> rateMatrixPtr( new GTRRateMatrix(metaNuc, rateParams, pi) );

  // dfg and factorSet
  TrctmcFactorSet factorSet = mkTrctmcFactorSet(treeInfo, rateMatrixPtr, 0.517981);
  DFG dfg = mkPhyloDfg(treeInfo, factorSet);

  // alg data
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );

  // stateMaskMapSet
  StateMaskMapSet stateMaskMapSet(metaNuc, dfg.variables.size());

  stateMask2DVec_t stateMask2DVec = mkStateMask2DVec(seqDataVec[4], treeInfo.nodeMap, stateMaskMapSet, 1);
  stateMaskVec_t stateMaskVec = stateMask2DVec[6]; // all missing data

  BOOST_CHECK_CLOSE(toNumber( dfg.calcNormConst(stateMaskVec) ), 1.0, EPS); // since all missing data
}



//string writeDot(FG_t const &tree)
//{
//  stringstream ss;
//  write_graphviz(ss, tree);
//  return ss.str();
//}
//
//void writeDot(FG_t const &tree, string const &fileName)
//{
//  ofstream f(fileName.c_str(), ios::out);
//  if (!f)
//    errorAbort("Cannot open file: " + fileName + "\n");
//  
//  f << writeDot(tree);
//  f.close();
//}
//
//
//BOOST_AUTO_TEST_CASE(mkFactorGraph_1) 
//{
//  // setup
//  string const IN_TREE_STRING = "('Bovine':0.693950,('Hylobates':0.360790,('Pongo':0.336360,('G._Gorilla':0.171470,('P._paniscus':0.192680,'H._sapiens':0.119270):0.083860):0.061240):0.150570):0.549390,'Rodent':1.214600);\n";
//  phy::tree_t phyTree = phy::readNewick(IN_TREE_STRING);
//  FG_t factorGraph; 
//  std::vector<unsigned> nodeMap; 
//  std::vector<unsigned> branchMap; 
//  std::vector<FGVertex_t> factorFirstVarMap;
//  unsigned initDistrNode;
//
//
//  mkFactorGraph(phyTree, factorGraph, nodeMap, branchMap, factorFirstVarMap, initDistrNode);
//
//  // output
//  BOOST_TEST_MESSAGE( num_vertices(factorGraph) );
//  writeDot(factorGraph, "fg.dot");
//
////   cout << "nodeMap: ";
////   BOOST_FOREACH(unsigned i, nodeMap)
////     cout << i << ", ";
////   cout << endl;
//// 
////   cout << "branchMap: ";
////   BOOST_FOREACH(unsigned i, branchMap)
////     cout << i << ", ";
////   cout << endl;
//}
//
//// predeclaration
//namespace phy {
//  std::vector<FGVertex_t> getNeighbors(FG_t const & factorGraph, FGVertex_t v);
//  std::vector<FGVertex_t> getNeighbors(FG_t const & factorGraph, FGVertex_t v, FGVertex_t firstVar);
//}
//
//BOOST_AUTO_TEST_CASE(getNeighbors_1) 
//{
//  // setup
//  string const IN_TREE_STRING = "('Bovine':0.693950,('Hylobates':0.360790,('Pongo':0.336360,('G._Gorilla':0.171470,('P._paniscus':0.192680,'H._sapiens':0.119270):0.083860):0.061240):0.150570):0.549390,'Rodent':1.214600);\n";
//  phy::tree_t phyTree = phy::readNewick(IN_TREE_STRING);
//  FG_t factorGraph; 
//  std::vector<unsigned> nodeMap; 
//  std::vector<unsigned> branchMap; 
//  std::vector<FGVertex_t> factorFirstVarMap;
//  unsigned initDistrNode;
//
//  mkFactorGraph(phyTree, factorGraph, nodeMap, branchMap, factorFirstVarMap, initDistrNode);
//
//  // test neighbors of the root node
//  std::vector<FGVertex_t> v0 = phy::getNeighbors(factorGraph, 0);
//  unsigned a0[] = {2, 4, 22, 23};
//  for (unsigned i = 0; i < v0.size(); i++)
//    BOOST_CHECK( v0[i] == a0[i] );
//
//  // test neighbors node 2 (true order taken from dot file)
//  std::vector<FGVertex_t> v2 = phy::getNeighbors(factorGraph, 2);
//  unsigned a2[] = {0, 1};
//  for (unsigned i = 0; i < v2.size(); i++)
//    BOOST_CHECK( v2[i] == a2[i] );
//}
//
//
//BOOST_AUTO_TEST_CASE(getNeighbors_2) 
//{
//  // setup
//  string const IN_TREE_STRING = "('Bovine':0.693950,('Hylobates':0.360790,('Pongo':0.336360,('G._Gorilla':0.171470,('P._paniscus':0.192680,'H._sapiens':0.119270):0.083860):0.061240):0.150570):0.549390,'Rodent':1.214600);\n";
//  phy::tree_t phyTree = phy::readNewick(IN_TREE_STRING);
//  FG_t factorGraph; 
//  std::vector<unsigned> nodeMap; 
//  std::vector<unsigned> branchMap; 
//  std::vector<FGVertex_t> factorFirstVarMap;
//  unsigned initDistrNode;
//
//  mkFactorGraph(phyTree, factorGraph, nodeMap, branchMap, factorFirstVarMap, initDistrNode);
//
//  // check first var version of getNeighbors
//  std::vector<FGVertex_t> v4 = phy::getNeighbors(factorGraph, 4, 0);
//  unsigned a4[] = {0, 3};
//  for (unsigned i = 0; i < v4.size(); i++)
//    BOOST_CHECK( v4[i] == a4[i] );
//}
//
//
//BOOST_AUTO_TEST_CASE(mkFactorGraphNodes_1) 
//{
//  // setup
//  string const IN_TREE_STRING = "('Bovine':0.693950,('Hylobates':0.360790,('Pongo':0.336360,('G._Gorilla':0.171470,('P._paniscus':0.192680,'H._sapiens':0.119270):0.083860):0.061240):0.150570):0.549390,'Rodent':1.214600);\n";
//  phy::tree_t phyTree = phy::readNewick(IN_TREE_STRING);
//  FG_t factorGraph; 
//  std::vector<unsigned> nodeMap; 
//  std::vector<unsigned> branchMap; 
//  std::vector<FGVertex_t> factorFirstVarMap;
//  unsigned initDistrNode;
//
//  mkFactorGraph(phyTree, factorGraph, nodeMap, branchMap, factorFirstVarMap, initDistrNode);
//
//  unsigned dim = 4;
//  ublas::matrix<number_t> Q = jukesCantorRM();
//  ublas::vector<number_t> branchLengths = extractBranchLengths(phyTree);
//  std::vector< ublas::matrix<number_t> > TM( branchLengths.size(), ublas::matrix<number_t>(dim, dim) );
//  updateTransitionMatrices(Q, branchLengths, TM);
//  ublas::vector<number_t> equiFreq = deriveEquiFreqForReversibleRM(Q);
//
//  vector<FGBaseNode *> nodes = mkFactorGraphNodes(nodeMap, factorGraph, branchMap, factorFirstVarMap, initDistrNode, TM, equiFreq);
//
//  StateMap metaNucMap = mkMetaNucStateMap();
//  StateMaskMap stateMask( metaNucMap );
//
//  // init observables
//  dynamic_cast<DiscreteVariableNode *>(nodes[1])->stateMask = & stateMask.metaState2StateMask( metaNucMap.symbol2State("-") );
//  dynamic_cast<DiscreteVariableNode *>(nodes[1])->observed = true;
//  dynamic_cast<DiscreteVariableNode *>(nodes[5])->stateMask = & stateMask.metaState2StateMask( metaNucMap.symbol2State("-") );
//  dynamic_cast<DiscreteVariableNode *>(nodes[5])->observed = true;
//  dynamic_cast<DiscreteVariableNode *>(nodes[9])->stateMask = & stateMask.metaState2StateMask( metaNucMap.symbol2State("-") );
//  dynamic_cast<DiscreteVariableNode *>(nodes[9])->observed = true;
//  dynamic_cast<DiscreteVariableNode *>(nodes[13])->stateMask = & stateMask.metaState2StateMask( metaNucMap.symbol2State("-") );
//  dynamic_cast<DiscreteVariableNode *>(nodes[13])->observed = true;
//  dynamic_cast<DiscreteVariableNode *>(nodes[17])->stateMask = & stateMask.metaState2StateMask( metaNucMap.symbol2State("-") );
//  dynamic_cast<DiscreteVariableNode *>(nodes[17])->observed = true;
//  dynamic_cast<DiscreteVariableNode *>(nodes[19])->stateMask = & stateMask.metaState2StateMask( metaNucMap.symbol2State("-") );
//  dynamic_cast<DiscreteVariableNode *>(nodes[19])->observed = true;
//  dynamic_cast<DiscreteVariableNode *>(nodes[21])->stateMask = & stateMask.metaState2StateMask( metaNucMap.symbol2State("-") );
//  dynamic_cast<DiscreteVariableNode *>(nodes[21])->observed = true;
//
//  unsigned root = 0;
//  phy::calcAllSumMessages(nodes, root);
//
//  std::vector<vector_t const *> & inMessages = nodes[root]->inMessages;
//
//  // node 2 -- should be all 1.0
//  unsigned nbIdx = nodes[root]->getNeighborIndex(2);
//  for (unsigned i = 0; i < dim; i++)
//    BOOST_CHECK_CLOSE( (*inMessages[nbIdx])[i], 1.0, 1e-10);
//
//  // node 4 -- should be all 1.0
//  nbIdx = nodes[root]->getNeighborIndex(4);
//  for (unsigned i = 0; i < dim; i++)
//    BOOST_CHECK_CLOSE( (*inMessages[nbIdx])[i], 1.0, 1e-10);
//
//  // node 22 -- should be all 1.0
//  nbIdx = nodes[root]->getNeighborIndex(22);
//  for (unsigned i = 0; i < dim; i++)
//    BOOST_CHECK_CLOSE( (*inMessages[nbIdx])[i], 1.0, 1e-10);
//
//  // node 23 -- should be all 0.25 (== equiFreq[i])
//  nbIdx = nodes[root]->getNeighborIndex(23);
//  for (unsigned i = 0; i < dim; i++)
//    BOOST_CHECK_CLOSE( (*inMessages[nbIdx])[i], equiFreq[i], 1e-10);
//}
//


