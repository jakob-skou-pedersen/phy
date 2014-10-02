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

// tested code (the defines allow access to private data)
// Beware though to avoid testing implementation details but only interface
#define private public
#define protected public
#include "phy/Mixtures.h"
#include "phy/Mixtures.cpp"
#undef protected
#undef private 

// input/output not included by above
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <iostream>
#include <sstream>

// other stuff
#include <boost/foreach.hpp>

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;

// helper functions

BOOST_AUTO_TEST_CASE(Mixture_Normal_1) 
{
  vector_t means(3);
  vector_t vars(3);
  means <<= 1,2,3;
  vars <<= 1,1,4;

  NormalMixture nm( means, vars, 0, 5, 10);
  matrix_t m(3,10,0);
  nm.mkFactor(m);
  //R: pnorm(1,1,1)-pnorm(0,1,1)
  BOOST_CHECK_CLOSE( m(0,0), 0.3413447461, 0.0001);
  //R: pnorm(1,2,1)-pnorm(0,2,1)
  BOOST_CHECK_CLOSE( m(1,0), 0.1359051220, 0.0001);
  //R: pnorm(2,3,2)-pnorm(1,3,2)
  BOOST_CHECK_CLOSE( m(2,1), 0.149882284795, 0.0001);
}


BOOST_AUTO_TEST_CASE(Mixture_Normal_IO_1)
{
  std::stringstream ss;
  //When using IO functions the internal scaling of parameters should be handled automatically
  ss << "MEANS:	[3](1, 3, 2)" << std::endl << "VARS:	[3](1,2,1)" << std::endl;
  
  NormalMixture nm( ss, 1, 6, 10);
  matrix_t m(3,10,0);
  nm.mkFactor(m);
  //R: pnorm(1.5,1,1)-pnorm(1,1,1)
  BOOST_CHECK_CLOSE( m(0,0), 0.1914624613, 0.0001);
  //R: pnorm(2,2,1)-pnorm(1.5,2,1)
  BOOST_CHECK_CLOSE( m(2,1), 0.1914624613, 0.0001);
  //R: pnorm(2.5,3,sqrt(2))-pnorm(2,3,sqrt(2))
  BOOST_CHECK_CLOSE( m(1,2), 0.1220867438, 0.0001);
}


BOOST_AUTO_TEST_CASE(Mixture_Gamma_1) 
{
  vector_t alphas(3);
  vector_t betas(3);
  alphas <<= 1,2,3;
  betas <<= 1,1,2;

  GammaMixture gm( alphas, betas, 0, 5, 10);
  //R: > pgamma(1,1,1) - pgamma(0,1,1)
  //BOOST_CHECK_CLOSE( gm.binProb(0,0), 0.6321205588, 0.0001);
  //R: > pgamma(1,2,1) - pgamma(0,2,1)
  //BOOST_CHECK_CLOSE( gm.binProb(0,1), 0.2642411177, 0.0001);
  //R: > pgamma(2,3,2) - pgamma(1,3,2)
  //BOOST_CHECK_CLOSE( gm.binProb(1,2), 0.4385731106, 0.0001);
}


BOOST_AUTO_TEST_CASE(Mixture_Gamma_IO_1)
{
  std::stringstream ss;
  ss << "ALPHAS:	[3](1, 3, 2)" << std::endl << "BETAS:	[3](1,2,1)" << std::endl;
  
  GammaMixture gm( ss, 1, 6, 10);
  matrix_t m(3,10,0);
  //R: > pgamma(1.5,1,1)-pgamma(1,1,1)
  BOOST_CHECK_CLOSE( m(0,0), 0.144749281, 0.0001);
  //R: > pgamma(2,2,1)-pgamma(1.5,2,1)
  BOOST_CHECK_CLOSE( m(2,1), 0.1518195507, 0.0001);
  //R: > pgamma(2.5,3,2)-pgamma(2,3,2)
  BOOST_CHECK_CLOSE( m(1,2), 0.1134512861, 0.0001);
}

