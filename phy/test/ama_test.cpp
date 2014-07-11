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
#include "phy/ama.h"

// input/output not included by above
#include <iostream>
#include <fstream>

// other stuff
#include <boost/foreach.hpp>

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;

BOOST_AUTO_TEST_CASE(ama_io_1) 
{
  vector<ama> amaVec;

  string inFileName("data/test1.ama");
  ifstream inFile(inFileName.c_str(), ios::in);
  if (!inFile)
    errorAbort("Cannot open file: " + inFileName + "\n");

  inFile >> amaVec;

  string outFileName("output/tmpTest1.ama");
  ofstream outFile(outFileName.c_str(), ios::out);
  
  outFile << amaVec;

  // to do: write test comparison of ama files
}

