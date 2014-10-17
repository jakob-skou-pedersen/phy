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
#include <string>
#include <vector>

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

BOOST_AUTO_TEST_CASE(conversionUtils_1) {
  ConvertIndexNumber cin1(0,10,10);
  BOOST_CHECK_CLOSE( cin1.getToNumberBeta(), 0.5, 0.001);
  BOOST_CHECK_CLOSE( cin1.getToNumberAlpha(), 1, 0.001);
  BOOST_CHECK_CLOSE( cin1.indexToNumber( 1 ), 1.5, 0.001);
  BOOST_CHECK_EQUAL( cin1.numberToIndex( 1.5), 1);

  ConvertIndexNumber cin2(-0.10,0.10,20);
  BOOST_CHECK_CLOSE( cin2.getToNumberBeta(), -.095, 0.001);
  BOOST_CHECK_CLOSE( cin2.getToNumberAlpha(), 0.01, 0.001);
  BOOST_CHECK_CLOSE( cin2.indexToNumber( 10 ), 0.005, 0.001);
  BOOST_CHECK_EQUAL( cin2.numberToIndex( 0.005), 10);

  ConvertIndexNumber cin3(-1, 2, 201);
}

BOOST_AUTO_TEST_CASE(maps_1)
{
  int tmpFrom[] = {5,1,2,3};
  int tmpTo[]   = {4,3,2,1};
  vector<int> from(tmpFrom, tmpFrom+4);
  vector<int> to(tmpTo, tmpTo+4);

  vector<unsigned> map = mkMap( to, from);

  BOOST_CHECK_EQUAL( map.at(0), -1);
  BOOST_CHECK_EQUAL( map.at(1), 3);
  BOOST_CHECK_EQUAL( map.at(2), 2);
  BOOST_CHECK_EQUAL( map.at(3), 1);
}
