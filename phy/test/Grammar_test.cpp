/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#define BOOST_TEST_DYN_LINK

// bost test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// phy test constants, etc.
#include "phy/test/PhyTestDef.h"

// tested code
// Note the hack to access protected data
#define private public
#define protected public

#include "phy/Grammar.h"
#include "phy/GrammarIO.h"
#include "phy/EmitWrapIO.h"

#define private public
#define protected public

// input/output not included by above
#include <boost/numeric/ublas/io.hpp>

// other stuff
#include <boost/foreach.hpp>

using namespace phy;
using namespace phy::grammar;
using boost::unit_test::test_suite;
using namespace std;
 
 
BOOST_AUTO_TEST_CASE(TransitionParser_1) 
{
  ifstream f;
  f.open("data/grammarOneState.txt");
 
  // skipping header for now
  string s;
  getline(f, s);
  getline(f, s);
 
  Transition trans;
  f >> trans;
 
  BOOST_CHECK(trans.name == "0.0");
  BOOST_CHECK(toNumber(trans.p) == 0.9);
}
 
 
BOOST_AUTO_TEST_CASE(GrammarParser_1) 
{
  // one state grammar
  Grammar g = readGrammar("data/grammarOneState.txt");
 
  BOOST_CHECK(g.stateCount == 1);
  BOOST_CHECK(g.outTransitions.size() == 1);
  BOOST_CHECK(g.outTransitions[0].size() == 2);
  BOOST_CHECK(g.inTransitions.size() == 1);
  BOOST_CHECK(g.inTransitions[0].size() == 1);
   
  // knudsen grammar
  g = readGrammar("data/grammarKnudsen99.txt");
 
  BOOST_CHECK_EQUAL( g.outTransitions.size(), unsigned(3) );
  BOOST_CHECK_EQUAL( g.inTransitions.size(), unsigned(3) );
  BOOST_CHECK_EQUAL( g.outTransitions[2].size(), unsigned(2) );  // Transitions of F state
  BOOST_CHECK_EQUAL( g.inTransitions[2].size(), unsigned(3) );
  BOOST_CHECK_EQUAL( g.stateCount, unsigned(3) );
 
  //  cout << g << endl;
}
 

xmatrix_t mkUniformPairEmitMatrix(unsigned L, xnumber_t const & p)
{
  xmatrix_t m(L, L);
  reset(m);
  for (unsigned i = 0; i < L; i++) 
    for (unsigned j = i + 1; j < L; j++) // min dist is 1 (i denotes beginning of left part of symbol and j denotes beginning of right part)
      m(i, j) = p;
  return m;
}

 
BOOST_AUTO_TEST_CASE(resetEmissions_2) 
{
  Grammar g = readGrammar("data/grammarKnudsen99.txt");
  
  unsigned L = 10;
  
  // emit models
  xvector_t probVec(L, 0.5);
  vector<xvector_t> sngEmit(1, probVec);
  vector<xmatrix_t> dblEmit(1, mkUniformPairEmitMatrix(L, 0.25) );
   
  g.resetEmissions(sngEmit, dblEmit);
 
  BOOST_CHECK_EQUAL(g.e(0, 0), 0.5);     // singleEmissions
  BOOST_CHECK_EQUAL(g.e(0, 0, 2), 0.25); // pairEmissions
  
  BOOST_CHECK_EQUAL(g.inLength_, L); // pairEmissions
}
 

BOOST_AUTO_TEST_CASE(inside_grammarOneState_1) 
{
  Grammar g = readGrammar("data/grammarOneState.txt");
 
  // emit models
  xvector_t missingData(4, 1);
  vector<xvector_t> sngEmit(1, missingData);
 
  BOOST_CHECK_EQUAL( g.inTransitions.size(), unsigned(1) );
  BOOST_CHECK_EQUAL(g.outTransitions.size(), unsigned(1) );
  BOOST_CHECK_EQUAL(g.stateCount, unsigned(1) );
 
  g.resetEmissions(sngEmit);
  BOOST_CHECK_EQUAL(g.inLength_, unsigned(4) );
 
  // regression tests
  BOOST_CHECK_SMALL( 0.06561 - toNumber( g.inside() ), EPS ); // inside
  BOOST_CHECK_SMALL( 0.06561 - toNumber( g.cyk() ), EPS );    // cyk
}
 

BOOST_AUTO_TEST_CASE(inside_grammarKnudsen99_1) 
{
  Grammar g = readGrammar("data/grammarKnudsen99.txt");
 
  unsigned L = 5;
  // emit models
  xvector_t probVec(L, 1.0);
  vector<xvector_t> sngEmit(1, probVec);
  vector<xmatrix_t> dblEmit(1, mkUniformPairEmitMatrix(L, 1.0) );
 
  g.resetEmissions(sngEmit, dblEmit);
 
  // regression tests
  BOOST_CHECK_SMALL( 0.0126953125 - toNumber( g.inside() ), EPS ); // inside
  BOOST_CHECK_SMALL( 0.00390625 - toNumber( g.cyk() ), EPS );      // cyk
}


BOOST_AUTO_TEST_CASE(phyloEmission_grammarKnudsen99_1) 
// test of phylo emission probs, not really of grammar
{
  EmitWrappers ew;
  readEmitWrappers("data/grammarKnudsen99EmitModels.txt", ew);
  Grammar g = readGrammar( "data/grammarEvoFoldStr.v1.txt", ew.vecWrapNames(), ew.matWrapNames() );
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  SeqData const & seqData = seqDataVec[5];
  vector<xmatrix_t> v = ew.matrixEmissions(seqData);

  // checking that the phylo models are working as they should (as much a check of the model as of the code)
  BOOST_CHECK_SMALL( toNumber( v[0](0,4) ) - 1.0, EPS ); // inside
}


 

BOOST_AUTO_TEST_CASE(inside_grammarOnlyPairs_1) 
{
  Grammar g = readGrammar("data/grammarOnlyPairs.txt");
 
  unsigned L = 4;
  // emit model
  vector<xmatrix_t> dblEmit(1, mkUniformPairEmitMatrix(L, 0.25) );
  g.resetEmissions(dblEmit);
 
  // regression tests
  BOOST_CHECK_SMALL( 0.0050625 - toNumber( g.inside() ), EPS ); // inside
  BOOST_CHECK_SMALL( 0.0050625 - toNumber( g.cyk() ), EPS );    // cyk
}
 


BOOST_AUTO_TEST_CASE(outside_grammarKnudsen99_1) 
{
  Grammar g = readGrammar("data/grammarKnudsen99.txt");

  unsigned L = 5;
  // emit models
  xvector_t probVec(L, 1.0);
  vector<xvector_t> sngEmit(1, probVec);
  vector<xmatrix_t> dblEmit(1, mkUniformPairEmitMatrix(L, 1.0) );

  g.resetEmissions(sngEmit, dblEmit);

  g.inside();
  g.outside();

  Grammar::vector2D_t tv;
  g.calcAccumulatedMarginals(tv);

  // add marginals for emitting states (*2 for pairs)
  ynumber_t emitMar = 2 * (tv[2][1] + tv[1][0]) + tv[2][0];
  BOOST_CHECK_CLOSE( toNumber(emitMar), toNumber(L), EPS ); // count accumulated marginal prob of all emitting transitions (weighing pairs by 2) should equal length

  //  for (unsigned v = 0; v < tv.size(); v++) {
  //    cout << endl;
  //    cout << "state:\t" << v << endl;
  //    for (unsigned u = 0; u < tv[v].size(); u++) {
  //      Transition const & t = g.outTransitions[v][u];
  //      cout << t.name << "\t" << tv[v][u] << endl;
  //    }
  //  }
}


BOOST_AUTO_TEST_CASE(maxParse_grammarKnudsen99_1) 
{
  Grammar g = readGrammar("data/grammarKnudsen99.txt");

  unsigned L = 5;
  // emit models
  xvector_t probVec(L, 0.1);
  vector<xvector_t> sngEmit(1, probVec);
  vector<xmatrix_t> dblEmit(1, mkUniformPairEmitMatrix(L, 1.0) );

  g.resetEmissions(sngEmit, dblEmit);
  g.cyk();
  vector<EmitInfo> emitInfoVec = g.maxParse();

  ynumber_t insideProb = g.insideOutside();
  // post probs of emissions in max parse
  vector_t ppv = g.postProbParse(emitInfoVec);

  // regression test
  BOOST_CHECK_SMALL( 1.0 - toNumber( ppv[2] ), EPS );    // 
  BOOST_CHECK_SMALL( 0.3330557868 - toNumber( ppv[4] ), EPS );    // 

  //  for (unsigned i = 0; i < ppv.size(); i++) 
  //    cout << i << "\t" << emitInfoVec[i].name << "\t" << ( (emitInfoVec[i].left) ? "left" : "right") << "\tpp\t" << ppv[i] << endl; 

  vector<Transition> tv = g.emitTransitions();
  for (unsigned i = 0; i < g.inLength_; i++) {
    ynumber_t sum = 0;
    for (unsigned j = 0; j < tv.size(); j++) {
      Transition & t = tv[j];
      if (t.type == PAIR) {
	sum += g.emitTransitionAccUseProb(t, i, true);
	sum += g.emitTransitionAccUseProb(t, i, false);
      }
      else
	sum += g.emitTransitionAccUseProb(t, i);
    }
    BOOST_CHECK_CLOSE(1.0, toNumber(sum / insideProb), EPS ); 
  }
 
  //  // output instead of test
  //  for (unsigned i = 0; i < g.inLength_; i++) {
  //    for (unsigned j = 0; j < tv.size(); j++) {
  //      Transition & t = tv[j];
  //      if (t.type == PAIR) {
  //	cout << i << "\t" << t.name + ".l"  << "\t" << g.emitTransitionAccUseProb(t, i, true) / insideProb << endl;
  //	cout << i << "\t" << t.name + ".r"  << "\t" << g.emitTransitionAccUseProb(t, i, false) / insideProb << endl;
  //      }
  //      else
  //	cout << i << "\t" << t.name << "\t" << g.emitTransitionAccUseProb(t, i) / insideProb << endl;
  //    }
  //    cout << endl;
  //  }
}


BOOST_AUTO_TEST_CASE(EmitAnno_1) 
{
  EmitAnnoMap eam     = readEmitAnnoMap("data/grammarKnudsen99AnnoMap.txt");
  // check single anno
  EmitAnno const & eaSingle = eam["single"];
  BOOST_CHECK( eaSingle.singleCompatible( "." ) );
  BOOST_CHECK( eaSingle.singleCompatible( "*" ) );
  BOOST_CHECK( not eaSingle.singleCompatible( ")" ) );

 // check pair anno
  EmitAnno const & eaPair = eam["pair"];
  BOOST_CHECK( eaPair.pairCompatible( "(", ")" ) );
  BOOST_CHECK( eaPair.pairCompatible( "(", "*" ) );
  BOOST_CHECK( eaPair.pairCompatible( "*", "*" ) );
  BOOST_CHECK( not eaPair.pairCompatible( ")", "*" ) );
  BOOST_CHECK( not eaPair.pairCompatible( "-", "j" ) );
}


BOOST_AUTO_TEST_CASE(AnnoMap_1) 
{
  // read grammar, load emmisions, run cyk
  Grammar g = readGrammar("data/grammarKnudsen99.txt");
  unsigned L = 5;
  // emit models
  xvector_t probVec(L, 0.1);
  vector<xvector_t> sngEmit(1, probVec);
  vector<xmatrix_t> dblEmit(1, mkUniformPairEmitMatrix(L, 1.0) );

  g.resetEmissions(sngEmit, dblEmit);
  g.cyk();

  // AnnoMap
  TransAnnoMap am = mkTransAnnoMap("data/grammarKnudsen99AnnoMap.txt", g);
  string anno = am.convertToString( g.maxParse() );

  // regression test
  BOOST_CHECK_EQUAL(anno, string(".(..)") );  
}


//
// 
//BOOST_AUTO_TEST_CASE(EmitWrap_1) 
//{
//  // read grammar, load emmisions, run cyk
//  EmitWrappers ew;
//  readEmitWrappers("data/grammarKnudsen99EmitModels.txt", ew);
//  Grammar g = readGrammar( "data/grammarKnudsen99.txt", ew.vecWrapNames(), ew.matWrapNames() );
//  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
//  SeqData const & seqData = seqDataVec[0];
//   
//  g.resetEmissions( ew.vectorEmissions(seqData), ew.matrixEmissions(seqData) );
//  g.cyk();
// 
//  // AnnoMap
//  AnnoMap am = readAnnoMap("data/grammarKnudsen99AnnoMap.txt");
//  string anno = am.convertToString( g.maxParse() );
// 
//  cout << anno << endl;
//}
// 

BOOST_AUTO_TEST_CASE(EvoFoldGrammarStr_v1_1) 
{
  // read grammar, load emmisions, run cyk
  EmitWrappers ew;
  readEmitWrappers("data/grammarKnudsen99EmitModels.txt", ew);
  Grammar g = readGrammar( "data/grammarEvoFoldStr.v1.txt", ew.vecWrapNames(), ew.matWrapNames() );
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  SeqData const & seqData = seqDataVec[0];
  
  g.resetEmissions( ew.vectorEmissions(seqData), ew.matrixEmissions(seqData) );
  g.cyk();

  // AnnoMap
  TransAnnoMap am = mkTransAnnoMap("data/grammarEvoFoldStrAnnoMap.v1.txt", g);
  string anno = am.convertToString( g.maxParse() );

  //  cout << anno << endl;
}


BOOST_AUTO_TEST_CASE(equivalentTransitions_v1_1) 
{
  EmitWrappers ew;
  readEmitWrappers("data/grammarKnudsen99EmitModels.txt", ew);
  Grammar g = readGrammar( "data/grammarEvoFoldStr.v1.txt", ew.vecWrapNames(), ew.matWrapNames() );

  TransAnnoMap am = mkTransAnnoMap("data/grammarEvoFoldStrAnnoMap.v1.txt", g);
  TransAnno ta;
  ta.anno = ".";

  // single emitting
  vector<string> v = am.equivalentTransitions(ta);
  BOOST_CHECK( hasElement(v, string("E.1") ) );  
  BOOST_CHECK( hasElement(v, string("E_NS.0") )  );  
  BOOST_CHECK( hasElement(v, string("E_NS.1") ) );  
  
  // pair emitting
  ta.anno = "";
  ta.annoLeft = "(";
  ta.annoRight = ")";
  v = am.equivalentTransitions(ta);
  BOOST_CHECK( hasElement(v, string("E.0") )  );  
  BOOST_CHECK( hasElement(v, string("P.0") ) );  
}


BOOST_AUTO_TEST_CASE(postProbParse_1) 
{
  EmitWrappers ew;
  readEmitWrappers("data/grammarKnudsen99EmitModels.txt", ew);
  Grammar g = readGrammar( "data/grammarEvoFoldStr.v1.txt", ew.vecWrapNames(), ew.matWrapNames() );
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  SeqData const & seqData = seqDataVec[4];
  
  g.resetEmissions( ew.vectorEmissions(seqData), ew.matrixEmissions(seqData) );
  g.cyk();
  g.insideOutside();
  vector<EmitInfo> vei = g.maxParse();
  
  // per position post probs
  xvector_t v1 = g.postProbParse(vei);

  // per position post probs summing over all transitions with equivalent annotations
  TransAnnoMap am = mkTransAnnoMap("data/grammarEvoFoldStrAnnoMap.v1.txt", g);
  xvector_t v2 = g.postProbParse( vei, mkEquivalenceMap(g, am) );

  // anno
  string anno = am.convertToString( g.maxParse() );

  BOOST_CHECK( v2[0] >= v1[0] );  
  BOOST_CHECK( v2[1] >= v1[1] );  
  BOOST_CHECK( v2[2] >= v1[2] );  

  // cout << ew.matrixEmissions(seqData)[0] << endl << endl;
  // cout << ew.vectorEmissions(seqData)[0] << endl << endl;
  // cout << anno << endl;
  // cout << toString(v1) << endl;
  // cout << toString(v2) << endl;
}


BOOST_AUTO_TEST_CASE(postProb_EvoFoldGrammarStr_1) 
{
  Grammar g = readGrammar("data/grammarEvoFoldStr.v1.txt");

  unsigned L = 5;
  // emit models
  xvector_t probVec(L, 0.1);
  vector<xvector_t> sngEmit(1, probVec);
  vector<xmatrix_t> dblEmit(1, mkUniformPairEmitMatrix(L, 1.0) );

  g.resetEmissions(sngEmit, dblEmit);

  ynumber_t insideProb = g.insideOutside();
  // post probs of emissions in max parse
  //  g.cyk();
  //  vector<EmitInfo> emitInfoVec = g.maxParse();
  //  vector_t ppv = g.postProbParse(emitInfoVec);
  //  
  //  cout << "in prob: " << insideProb << endl;
  //  cout << ppv << endl;
  //
  //  // print state order for out-transitions as well as in-transitions 
  //  BOOST_CHECK( g.stateCount == g.outTransitions.size() );
  //  BOOST_CHECK( g.outTransitions.size() == g.inTransitions.size() );
  //  cout << "out-transition state order" << endl;
  //  for (unsigned v = 0; v < g.stateCount; v++) 
  //    cout << "state index " << v << ": " << g.outTransitions[v][0].from_id << endl;
  //
  //  cout << endl << "in-transition state order" << endl;
  //  for (unsigned v = 0; v < g.stateCount; v++) 
  //    if (g.inTransitions[v].size() > 0)
  //      cout << "state index " << v << ": " << g.inTransitions[v][0].to_id << endl;

  vector<Transition> tv = g.emitTransitions();
  for (unsigned i = 0; i < g.inLength_; i++) {
    ynumber_t sum = 0;
    for (unsigned j = 0; j < tv.size(); j++) {
      Transition & t = tv[j];
      if (t.type == PAIR) {
	sum += g.emitTransitionAccUseProb(t, i, true);
	sum += g.emitTransitionAccUseProb(t, i, false);
      }
      else
	sum += g.emitTransitionAccUseProb(t, i);
    }
    //    cout << "pos " << i << ": " << sum / insideProb << endl;

    BOOST_CHECK_CLOSE(1.0, toNumber(sum / insideProb), EPS ); 
  }
}
