/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#define BOOST_TEST_DYN_LINK

// bost test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// convenient definitions 
#include "phy/test/PhyTestDef.h"

// tested code
// Note the hack to access protected data
#define private public
#define protected public

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


BOOST_AUTO_TEST_CASE(PhyloTreeIO_readPhyloTree_sngNuc_1) 
{
  PhyloTree phyloTree = readPhyloTree("./data/sngNucPhyloModel.txt");

  // read in some test data
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  vector<stateMask2DVec_t> stateMask3DVec = mkStateMask3DVec(seqDataVec, phyloTree);

  // do a basic calculation with the dfg (tree)
  stateMaskVec_t stateMaskVec = stateMask3DVec[4][6]; // all missing data
  BOOST_CHECK_CLOSE(toNumber( phyloTree.dfg.calcNormConst(stateMaskVec) ), 1.0, EPS); // since all missing data

  // ouput rate matrix and equiFreqs for comparison with input
  // cout << phyloTree.rmPtr->mkRateMatrix() << endl << endl;
  // cout << phyloTree.rmPtr->equiFreqs() << endl << endl;
}


BOOST_AUTO_TEST_CASE(PhyloTreeIO_readPhyloTree_sngNuc_2) 
{
  PhyloTree phyloTree = readPhyloTree("./data/sngNucPhyloModel.v2.txt");

  // read in some test data
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  vector<stateMask2DVec_t> stateMask3DVec = mkStateMask3DVec(seqDataVec, phyloTree);

  // do a basic calculation with the dfg (tree)
  stateMaskVec_t stateMaskVec = stateMask3DVec[4][6]; // all missing data
  BOOST_CHECK_CLOSE(toNumber( phyloTree.dfg.calcNormConst(stateMaskVec) ), 1.0, EPS); // since all missing data
}


BOOST_AUTO_TEST_CASE(PhyloTreeIO_readPhyloTree_sngNuc_3) 
{
  PhyloTree phyloTree = readPhyloTree("./data/sngNucPhyloModel.v3.txt");

  // read in some test data
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  vector<stateMask2DVec_t> stateMask3DVec = mkStateMask3DVec(seqDataVec, phyloTree);

  // do a basic calculation with the dfg (tree)
  stateMaskVec_t stateMaskVec = stateMask3DVec[4][6]; // all missing data
  BOOST_CHECK_CLOSE(toNumber( phyloTree.dfg.calcNormConst(stateMaskVec) ), 1.0, EPS); // since all missing data

  //  cout << phyloTree.rmPtr->mkRateMatrix() << endl << endl;
  //  cout << phyloTree.rmPtr->equiFreqs() << endl << endl;
}


BOOST_AUTO_TEST_CASE(PhyloTreeIO_readPhyloTree_diNuc_2) 
{
  PhyloTree phyloTree = readPhyloTree("./data/dinucPhyloModel.v2.txt");

  // read in some test data
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  vector<stateMask2DVec_t> stateMask3DVec = mkStateMask3DVec(seqDataVec, phyloTree);

  // do a basic calculation with the dfg (tree)
  stateMaskVec_t stateMaskVec = stateMask3DVec[4][10]; // all missing data
  BOOST_CHECK_CLOSE(toNumber( phyloTree.dfg.calcNormConst(stateMaskVec) ), 1.0, EPS); // since all missing data

  // ouput rate matrix and equiFreqs for comparison with input
  //  cout << phyloTree.rmPtr->mkRateMatrix() << endl << endl;
  //  cout << phyloTree.rmPtr->equiFreqs() << endl << endl;
}


BOOST_AUTO_TEST_CASE(PhyloTreeIO_readPhyloTree_diNuc_3) 
{
  PhyloTree phyloTree = readPhyloTree("./data/dinucPhyloModel.v3.txt");

  // read in some test data
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  vector<stateMask2DVec_t> stateMask3DVec = mkStateMask3DVec(seqDataVec, phyloTree);

  // do a basic calculation with the dfg (tree)
  stateMaskVec_t stateMaskVec = stateMask3DVec[4][10]; // all missing data
  BOOST_CHECK_CLOSE(toNumber( phyloTree.dfg.calcNormConst(stateMaskVec) ), 1.0, EPS); // since all missing data

  // ouput rate matrix and equiFreqs for comparison with input
  // cout << phyloTree.rmPtr->mkRateMatrix() << endl << endl;
  // cout << phyloTree.rmPtr->equiFreqs() << endl << endl;
}


xvector_t toXVector(vector<xnumber_t> const & v)
{
  xvector_t u( v.size() );
  for (unsigned i = 0; i < u.size(); i++)
    u[i] = v[i];
  return u;
}


BOOST_AUTO_TEST_CASE(PhyloTreeIO_calcLikelihoodVector_sngNuc_1) 
{
  PhyloTree phyloTree = readPhyloTree("./data/sngNucPhyloModel.v3.txt");

  // read in some test data
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  vector<stateMask2DVec_t> stateMask3DVec = mkStateMask3DVec(seqDataVec, phyloTree);
  xvector_t result1( symbolCount(seqDataVec[0].sequences, 1) );
  vector<xnumber_t> tmpResult2( symbolCount(seqDataVec[0].sequences, 1) );
  
  // check that calcLikelihoodVector gives the same results as calcNormConsMultObs
  calcLikelihoodVector(result1, seqDataVec[0], phyloTree);
  calcNormConsMultObs(tmpResult2, stateMask3DVec[0], phyloTree.dfg);

  xvector_t result2 = toXVector(tmpResult2);
  assert( result1.size() == result2.size() );
  for (unsigned i = 0; i < result1.size(); i++)
    BOOST_CHECK_CLOSE(toNumber(result1[i]), toNumber(result2[i]), EPS); // since all missing data

  // check that calcLikelihoodVector gives the same results when using hashing (currently a map)
  probMap_t probMap;
  calcLikelihoodVector(result2, seqDataVec[0], phyloTree, probMap, true);
  for (unsigned i = 0; i < result1.size(); i++)
    BOOST_CHECK_CLOSE(toNumber(result1[i]), toNumber(result2[i]), EPS); // since all missing data

  BOOST_CHECK(probMap.size() > 0); // since all missing data

  //  map<stateMaskVec_t, xnumber_t>::iterator iter;
  //  cout << probMap.size() << endl;
  //  for ( iter = probMap.begin(); iter != probMap.end(); iter++)
  //    cout << iter->second << ", ";
  //  cout << endl;
  //
  //  cout << result1 << endl;
  //  cout << result2 << endl;
}


BOOST_AUTO_TEST_CASE(PhyloTreeIO_calcLikelihoodMatrix_diNuc_1) 
{
  PhyloTree phyloTree = readPhyloTree("./data/dinucPhyloModel.v3.txt");

  // read in some test data
  vector<SeqData> seqDataVec = amaToSeqData( readAmaFile("data/test1.ama") );
  vector<stateMask2DVec_t> stateMask3DVec = mkStateMask3DVec(seqDataVec, phyloTree);

  xmatrix_t result = calcLikelihoodMatrix(seqDataVec[3], phyloTree);

  // cout << result << endl;
}

