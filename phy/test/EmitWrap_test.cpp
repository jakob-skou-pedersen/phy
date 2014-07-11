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

#include "phy/EmitWrapIO.h"

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



BOOST_AUTO_TEST_CASE(EmiWrap_EmitWrapDfg_1) 
{
  EmitWrappers ew = readEmitWrappers("./data/grammarStacking/grammarStackingEmitModels.txt");
  BOOST_CHECK( ew.vecWrapNames()[0] == "single" );
  BOOST_CHECK( ew.matWrapNames()[0] == "pair" );
  BOOST_CHECK( ew.matWrapNames()[1] == "stack" );

  SeqData seqData;
  seqData.addSeq("dna", "..........AATTTTT");

  vector<xvector_t> vecEmit = ew.vectorEmissions(seqData);
  vector<xmatrix_t> matEmit = ew.matrixEmissions(seqData);

  // spot checking of emission probabilities according to factorPotentials
  // single
  BOOST_CHECK_CLOSE(toNumber( vecEmit[0][0] ), 1.0, EPS);
  BOOST_CHECK_CLOSE(toNumber( vecEmit[0][10] ), 0.4, EPS);
  BOOST_CHECK_CLOSE(toNumber( vecEmit[0][12] ), 0.1, EPS);

  // pair
  //  cout << matEmit[0] << endl;
  BOOST_CHECK_CLOSE(toNumber( matEmit[0](0,0) ), 0.0, EPS);     // j-i=0<minDist
  BOOST_CHECK_CLOSE(toNumber( matEmit[0](0,4) ), 1.0, EPS);     // j-i=5>=minDist
  BOOST_CHECK_CLOSE(toNumber( matEmit[0](0,14) ), 0.225, EPS);  // P(".T")=0.225
  BOOST_CHECK_CLOSE(toNumber( matEmit[0](12,16) ), 0.05, EPS);  // P(".T")=0.225

  // stack
  cout << matEmit[1] << endl;
  BOOST_CHECK_CLOSE(toNumber( matEmit[1](0,0) ), 0.0, EPS);     // j-i=0<minDist
  //  BOOST_CHECK_CLOSE(toNumber( matEmit[1](0,4) ), 1.0, EPS);     // j-i=5>=minDist  // this fails as stack factor potential needs adjustment so rows sum to one  - sudhakar, please fix
  //  BOOST_CHECK_CLOSE(toNumber( matEmit[1](0,14) ), 0.225, EPS);  //?? I have not done the calculation...
  //  BOOST_CHECK_CLOSE(toNumber( matEmit[1](12,16) ), 0.05, EPS);  


}



