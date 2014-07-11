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
#include "phy/PhyDef.h"

// input/output not included by above
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

// other stuff
#include <boost/foreach.hpp>

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;

BOOST_AUTO_TEST_CASE(toNumber_1) 
{
  xnumber_t x = 0.22e-22;
  number_t y = toNumber(x);
  
  BOOST_CHECK_CLOSE(y, 0.22e-22, EPS);
}

BOOST_AUTO_TEST_CASE(xnumber_t_1) 
{
//   // testing xnumber_t
//    xnumber_t a(2.1);
//    for (unsigned i = 0; i < 100; i++) {
//      cout << a << endl;
//      a *= 2;
//    }
//    cout << endl;
//  
//    xnumber_t b(2.1);
//    for (unsigned i = 0; i < 100; i++) {
//      cout << b << endl;
//      b *= b;
//    }
//    cout << endl;
//  
// 
//   // test boost matrix with xdouble -- works when double constructor is added
//   ublas::matrix<xnumber_t> m (3, 3);
//   for (unsigned i = 0; i < m.size1 (); ++ i)
//     for (unsigned j = 0; j < m.size2 (); ++ j) {
//       m (i, j) = (3.0 * i + j)/100;
//     }
//   std::cout << m << std::endl;
// 
//   m = m / 2.0;
//   std::cout << trans(m) << std::endl;
// 
//   std::cout << std::endl;
//   for (unsigned i = 0; i < 10; ++ i)
//     m = prod(m,m);
//   std::cout << m << std::endl;
// 
// 
//   matrix_t o(3, 3);
//   for (unsigned i = 0; i < o.size1 (); ++ i)
//     for (unsigned j = 0; j < o.size2 (); ++ j) {
//       o (i, j) = (3.0 * i + j);
//     }
//   cout << "o: " << o << endl;
// 
//   xvector_t v(3);
//   for (unsigned i = 0; i < v.size (); ++ i)
//     v(i) = i;
//   cout << "v: " << v << endl;
//   
//   xvector_t w = prod(o, v);
//   cout << "prod(o, v) " << w << endl;
//   cout << "power(w[0] / 0.9, 1000) " << power(w[0] / 0.9, 1000) << endl;
// 
//   // test xnumber_t
//   xnumber_t x = 0.5e200;
//   xnumber_t y = 0.25e200;
//   cout << "x + y = " << x + y << endl;
//   cout << "x * 10**1000 + y * 10**1000 = " << x * NTL::power(xnumber_t(10), 1000) + y * NTL::power(xnumber_t(10), 1000) << endl;
// 
//   xnumber_t z = NTL::power2_xdouble(100);
//   cout << "z=2**100 = " << z << endl;
// 
//   cout << "log(z) = " << NTL::log(z) << endl;
//   cout << "log(10**50) = " << NTL::log( NTL::power(xnumber_t(10), 50) ) << endl;
//   cout << "log_10(10**50) = " <<  NTL::log( NTL::power(xnumber_t(10), 50) ) / log(10) << endl;
}
