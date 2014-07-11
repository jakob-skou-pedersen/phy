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

#include "phy/StatDist.h"

#define private public
#define protected public

// input/output not included by above
#include <boost/numeric/ublas/io.hpp>

// other stuff
#include  "phy/PhyDef.h"
#include "phy/utils.h"
#include <boost/foreach.hpp>
#include <boost/math/distributions/binomial.hpp>

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;
 
BOOST_AUTO_TEST_CASE(multinomialDist_1) 
{
  vector<number_t> pv(2);
  vector<unsigned> cv(2);
  // define a simple case (binomial)
  pv[0] = 0.7;
  pv[1] = 0.3;
  cv[0] = 2;
  cv[1] = 2;
  
  //  cout << "multinomialCoefficient<number_t>(cv):\t" << multinomialCoefficient<number_t>(cv) << endl;
  //  cout << "multinomialDistribution<number_t>(pv, cv):\t" << multinomialDistribution<number_t>(pv, cv) << endl;

  boost::math::binomial_distribution<number_t> bnDist(2 + 2, 0.7);
  number_t result = boost::math::pdf(bnDist, 2);
  BOOST_CHECK_CLOSE(multinomialDistribution(pv, cv), result, EPS);
}


BOOST_AUTO_TEST_CASE(multinomialDist_2) 
{
  typedef NTL::xdoubleMod ntlNumber_t;

  vector<ntlNumber_t> pv(2);
  vector<unsigned> cv(2);
  // define a simple case (binomial)
  pv[0] = 0.7;
  pv[1] = 0.3;
  cv[0] = 2;
  cv[1] = 2;
  
  boost::math::binomial_distribution<number_t> bnDist(2 + 2, 0.7);
  number_t result = boost::math::pdf(bnDist, 2);
  BOOST_CHECK_CLOSE(NTL::to_double( multinomialDistribution(pv, cv) ), result, EPS);
}

//
//BOOST_AUTO_TEST_CASE(multinomialDist_3) 
//{
//  typedef NTL::xdoubleMod ntlNumber_t;
//
//  vector<ntlNumber_t> pv(2);
//  vector<unsigned> cv(2);
//  // define a simple case (binomial)
//  pv[0] = 0.5;
//  pv[1] = 0.5;
//  cv[0] = 530000;
//  cv[1] = 470000;
//  
//  cout << "multinomialCoefficient<ntlNumber_t>(cv):\t" << multinomialCoefficient<ntlNumber_t>(cv) << endl;
//  cout << "multinomialDistribution(pv, cv):\t" << multinomialDistribution(pv, cv) << endl;
//}
