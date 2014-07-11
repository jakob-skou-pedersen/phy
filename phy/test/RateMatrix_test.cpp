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
#include "phy/RateMatrix.h"

// input/output not included by above
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

// other stuff
#include <boost/foreach.hpp>

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;

BOOST_AUTO_TEST_CASE(GTRRateMatrix_mkRateMatrix_1) 
{
  StateMap metaNucMap = mkMetaNucStateMap();

  vector_t equiFreq(4);
  equiFreq[0] = 0.1;
  equiFreq[1] = 0.5;
  equiFreq[2] = 0.1;
  equiFreq[3] = 0.3;

  vector_t params(6);
  params[0] = 1;
  params[1] = 1;
  params[2] = 1;
  params[3] = 1;
  params[4] = 1;
  params[5] = 1;

  GTRRateMatrix rmFac(metaNucMap, params, equiFreq);
  matrix_t Q =  rmFac.mkRateMatrix();
  BOOST_CHECK(Q(0,0) < 0);
  BOOST_CHECK_CLOSE(Q(0,2), Q(2,0), EPS); // since pi_0 == pi_2
  BOOST_CHECK_CLOSE(Q(0,1), Q(0,2) * 5, EPS);
  BOOST_CHECK_CLOSE(Q(0,2) * 3, Q(0,3), EPS);
}

BOOST_AUTO_TEST_CASE(GTRRateMatrix_mkRateMatrix_2) 
{
  // rate matrix
  StateMap metaNuc = mkMetaNucStateMap();
  
  vector_t pi(4);
  pi[0] = 0.1;
  pi[1] = 0.1;
  pi[2] = 0.1;
  pi[3] = 0.7;

  vector_t rateParams(6);
  rateParams[0] = 1;
  rateParams[1] = 2;
  rateParams[2] = 3;
  rateParams[3] = 4;
  rateParams[4] = 5;
  rateParams[5] = 6;

  StateMap metaNucMap = mkMetaNucStateMap();
  GTRRateMatrix rm(metaNucMap, rateParams, pi);

  matrix_t Q =  rm.mkRateMatrix();
  BOOST_CHECK(Q(0,0) < 0);
  //  BOOST_CHECK_CLOSE(Q(0,1), Q(2,0), EPS); // since pi_0 == pi_2
  //  BOOST_CHECK_CLOSE(Q(0,1), Q(0,2) * 5, EPS);
  //  BOOST_CHECK_CLOSE(Q(0,2) * 3, Q(0,3), EPS);
  //
  //  cout << rm.mkRateMatrix() << endl;
}


BOOST_AUTO_TEST_CASE(deriveEquiFreqAndGtrParamsForReversibleRM_1)
{
  StateMap metaNuc = mkMetaNucStateMap();
  
  vector_t equiFreqs(4);
  equiFreqs[0] = 0.1;
  equiFreqs[1] = 0.1;
  equiFreqs[2] = 0.1;
  equiFreqs[3] = 0.7;

  vector_t rateParams(6);
  rateParams[0] = 1;
  rateParams[1] = 2;
  rateParams[2] = 3;
  rateParams[3] = 4;
  rateParams[4] = 5;
  rateParams[5] = 6;

  StateMap metaNucMap = mkMetaNucStateMap();
  GTRRateMatrix rm(metaNucMap, rateParams, equiFreqs);

  vector_t pv = rm.rateParameters();
  vector_t pi = rm.equiFreqs();
  matrix_t Q =  rm.mkRateMatrix();
  vector_t pvDerived, piDerived; 

  deriveEquiFreqAndGtrParamsForReversibleRM(Q, piDerived, pvDerived);

  for (unsigned i = 0; i < pi.size(); i++)
    BOOST_CHECK_CLOSE(pi[i], piDerived[i], EPS);

  for (unsigned i = 1; i < pv.size(); i++)
    BOOST_CHECK_CLOSE(pv[0] / pvDerived[0], pv[i] / pvDerived[i], EPS);  // should be proportional, but not necessarily equal due to rate normalization

  //  cout << "pi         " << pi << endl;
  //  cout << "piDerived: " << piDerived << endl;
  //
  //  cout << "pv         " << pv << endl;
  //  cout << "pvDerived: " << pvDerived << endl;
}
