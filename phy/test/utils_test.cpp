/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#define BOOST_TEST_DYN_LINK

// bost test
#include <boost/test/unit_test.hpp>

// tested code
#include "phy/utils.h"

// input/output not included by above
#include <iostream>

// other stuff
#include <boost/foreach.hpp>

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;

// testing how simple tag-value pairs are being parsed
BOOST_AUTO_TEST_CASE(parserUtils_1) 
{
  ifstream f;
  f.open("data/parserTest.txt");

  int intVal = 0;
  string s;
  getline(f, s);
  getFeature(s, "PARSER_TEST", intVal);
  BOOST_CHECK(intVal == 1);
  
  skipWhiteSpaceAndComments(f);
  
  string tag;
  string strVal;
  f >> tag;
  f >> strVal;
  skipLine(f);
  BOOST_CHECK(tag == "NAME1:");
  BOOST_CHECK(strVal == "first");
  BOOST_CHECK( moreTags(f) );

  f >> tag;
  f >> strVal;
  skipLine(f);
  BOOST_CHECK(tag == "NAME2:");
  BOOST_CHECK(strVal == "second");
  BOOST_CHECK( not moreTags(f) );

  // next section
  skipWhiteSpaceAndComments(f);

  // test of comment at end of line
  f >> tag;
  f >> strVal;
  skipLine(f);
  BOOST_CHECK(tag == "NAME1:");
  BOOST_CHECK(strVal == "first.2");  
  BOOST_CHECK( moreTags(f) );

  f >> tag;
  f >> strVal;
  skipLine(f);
  BOOST_CHECK(tag == "NAME2:");
  BOOST_CHECK(strVal == "second.2");  
  BOOST_CHECK( not moreTags(f) );
}
