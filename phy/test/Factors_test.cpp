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
#define private public
#define protected public
#include "phy/Factors.h"
#include "phy/FactorsIO.h"
#undef protected
#undef private 

// input/output not included by above
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

// other stuff
#include <boost/foreach.hpp>

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;

// helper functions
matrix_t mkAndFillSquareMatrixIncr(unsigned n)
{
  matrix_t m(n, n);
  for (unsigned i = 0; i < n; i++) 
    for (unsigned j = 0; j < n; j++) 
      m(i, j) = i * n + j;
  return m;
}


matrix_t mkAndFillSquareMatrixUnif(unsigned n)
{
  matrix_t m(n, n);
  for (unsigned i = 0; i < n; i++) 
    for (unsigned j = 0; j < n; j++) 
      m(i, j) = 1.0 / (n * n);
  return m;
}


BOOST_AUTO_TEST_CASE(GlobalNormFactor_mkFactor_1) 
{
  GlobalNormFactor fac("noName", mkAndFillSquareMatrixIncr(4));
  matrix_t n = fac.mkFactor();

  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      BOOST_CHECK(n(i, j) == i * 4 + j);
}

BOOST_AUTO_TEST_CASE(GlobalNormFactor_optimizeParameters_1) 
{
  matrix_t m = mkAndFillSquareMatrixIncr(4);
  matrix_t pseudo = mkAndFillSquareMatrixIncr(4);
  GlobalNormFactor fac("noName", m , pseudo); // m will be overridden by optimization below
  fac.submitCounts(m);
  fac.optimizeParameters();
  matrix_t n = fac.mkFactor();

  BOOST_CHECK_CLOSE(sumMatrix(n), 1.0, EPS);
  BOOST_CHECK_CLOSE(n(0, 0), 0.0, EPS);
  BOOST_CHECK_CLOSE(n(3, 3), 0.125, EPS);
}


BOOST_AUTO_TEST_CASE(GlobalNormFactor_constructor_1) 
{
  matrix_t m = mkAndFillSquareMatrixIncr(4);
  matrix_t pseudo = mkAndFillSquareMatrixIncr(4);
  GlobalNormFactor fac("noName", m, pseudo);

  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      BOOST_CHECK(fac.m_(i, j) == i * 4 + j);
}


BOOST_AUTO_TEST_CASE(GlobalNormFactor_submitCounts_1) 
{
  matrix_t m = mkAndFillSquareMatrixIncr(4);
  GlobalNormFactor fac("noName", m);

  // uniform count matrix
  matrix_t c(4, 4);
  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      c(i, j) = 1;

  fac.submitCounts(c);
  fac.optimizeParameters();
  matrix_t par = fac.m_;

  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      BOOST_CHECK_CLOSE(par(i, j), 1.0/16, EPS);
}  


BOOST_AUTO_TEST_CASE(GlobalNormFactor_submitCounts_2) 
{
  matrix_t m = mkAndFillSquareMatrixIncr(4);
  matrix_t pseudo = mkAndFillSquareMatrixIncr(4);
  GlobalNormFactor fac("noName", m, pseudo);

  // uniform count matrix
  matrix_t c(4, 4);
  reset(c);

  fac.submitCounts(c);
  fac.optimizeParameters();

  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      BOOST_CHECK_CLOSE(fac.m_(i, j), pseudo(i, j) / sumMatrix(pseudo), EPS);
}  


BOOST_AUTO_TEST_CASE(RowNormFactor_general_1) 
{
  matrix_t m(4, 4);
  reset(m);
  RowNormFactor fac("noName", m);

  matrix_t c = mkAndFillSquareMatrixIncr(4);
  fac.submitCounts(c);
  fac.optimizeParameters();

  matrix_t n = fac.mkFactor();
    BOOST_CHECK(n(0, 0) == 0.0);
  for (unsigned i = 0; i < 4; i++)
    BOOST_CHECK_CLOSE(sum( row(n,i) ), 1.0, EPS);
}


BOOST_AUTO_TEST_CASE(ColumnNormFactor_general_1) 
{
  matrix_t m(4, 4);
  reset(m);
  ColumnNormFactor fac("noName", m);

  matrix_t c = mkAndFillSquareMatrixIncr(4);
  fac.submitCounts(c);
  fac.optimizeParameters();

  matrix_t n = fac.mkFactor();
    BOOST_CHECK(n(0, 0) == 0.0);
    BOOST_CHECK_CLOSE(n(3, 0), 0.5, EPS);

  for (unsigned i = 0; i < 4; i++)
    BOOST_CHECK_CLOSE(sum( column(n,i) ), 1.0, EPS);

  //  cout << n << endl;
}

//BOOST_AUTO_TEST_CASE(DiscContFactor_general_1)
//{
//  vector_t means(1, 100);
//  vector_t vars(1, 2500);
//  DiscContFactor dcf("noName", means, vars, 0, 2, 1, 200);
//  
//  //TODO Add some data and run init
//  //  dcf.counts_(0,i) = -(i-3)*(i-3)+64;
// 
//  //  nf.init(); 
//  //  nf.optimizeParameters();
//
//  /* R-code
//     > bp <- c(0.6,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.4)
//     > counts <- c((55),(60),(63),(64),(63),(60),(55),(48),(39),(28),(15))
//     > mean(sweep(as.matrix(bp),1,counts,FUN=rep))
//     [1] 1.347273
//     > var(sweep(as.matrix(bp),1,counts,FUN=rep))
//     [1] 0.290621
//  */
//
//  //Test equality with some margin
//  //  BOOST_CHECK_CLOSE( nf.mean_ , 1.347273, 0.001);
//  //  BOOST_CHECK_CLOSE( nf.var_  , 0.290621, 0.001);
//}
//
BOOST_AUTO_TEST_CASE(DiscContFactor_general_1){
  unsigned bins = 10;
  double minv = -2;
  double maxv = 3;
  int states = 2;

  //Initialize means and variances
  vector_t means(states);
  vector_t vars(states, bins*bins/1.96/1.96);
  for(unsigned i = 0; i < states; ++i)
    means(i) = (double)bins/states*(0.5+i);
 
  DiscContFactor dcf("dcf", means, vars, minv, maxv, states, bins);

  matrix_t c(2, 10, 1);
  for( int i = 0; i < 2; ++i)
    for( int j = 0; j < 10; ++j)
      c(i,j) = i+j+3;

  dcf.submitCounts(c);
  dcf.optimizeParameters();


  dcf.optimizeParametersImpl();
}


BOOST_AUTO_TEST_CASE(CompositeFactorSet_mkFactor_1) 
{
  // define factor set
  matrix_t m0 = mkAndFillSquareMatrixIncr(4);
  matrix_t m1 = mkAndFillSquareMatrixUnif(4);

  AbsBasFacPtr_t fac1Ptr(new GlobalNormFactor("noName", m0) );
  AbsBasFacPtr_t fac2Ptr(new GlobalNormFactor("noName", m1) );

  // note that factor 0 and 2 are the same.
  std::vector<AbsBasFacPtr_t> factorPtrs;
  factorPtrs.push_back(fac1Ptr);
  factorPtrs.push_back(fac2Ptr);
  factorPtrs.push_back(fac1Ptr);

  CompositeFactorSet facSet(factorPtrs);

  // test mkFactor
  BOOST_CHECK(matrixEqual(facSet.mkFactor(0), m0) );
  BOOST_CHECK(matrixEqual(facSet.mkFactor(1), m1) );
  BOOST_CHECK(matrixEqual(facSet.mkFactor(2), m0) );

  BOOST_CHECK(matrixEqual(facSet.mkFactorVec()[0], m0) );
  BOOST_CHECK(matrixEqual(facSet.mkFactorVec()[1], m1) );
  BOOST_CHECK(matrixEqual(facSet.mkFactorVec()[2], m0) );

  // check submitCounts and optimizeParameters
  facSet.submitCounts(m1, 0);
  facSet.submitCounts(m1, 1);
  facSet.submitCounts(m1, 2);

  facSet.optimizeParameters();

  BOOST_CHECK(matrixEqual(facSet.mkFactor(0), m1) );
  BOOST_CHECK(matrixEqual(facSet.mkFactor(1), m1) );
  BOOST_CHECK(matrixEqual(facSet.mkFactor(2), m1) );

  // diff counts
  facSet.submitCounts(m1, 0);
  facSet.submitCounts(m1, 1);
  facSet.submitCounts(m0, 2);

  facSet.optimizeParameters();

  BOOST_CHECK(not matrixEqual(facSet.mkFactor(0), m1) );
  BOOST_CHECK(matrixEqual(facSet.mkFactor(1), m1) );
  BOOST_CHECK(not matrixEqual(facSet.mkFactor(2), m1) );

  BOOST_CHECK(matrixEqual( facSet.mkFactor(2), facSet.mkFactor(0) ) );

  // output factor matrices
  // BOOST_FOREACH(matrix_t const & m, facSet.mkFactor() )
  //   cout << m << endl << sumMatrix(m) << endl << endl;

}


BOOST_AUTO_TEST_CASE(readFactorFile_1) 
{
  map<string, AbsBasFacPtr_t> factorMap = readFactorFile("./data/dfgSpec/test1Potentials.txt");

  BOOST_CHECK(factorMap["prior"] != NULL);
  BOOST_CHECK_EQUAL(factorMap["prior"]->mkFactor().size1(), (unsigned) 1);
  BOOST_CHECK_EQUAL(factorMap["prior"]->mkFactor().size2(), (unsigned) 2);
  BOOST_CHECK_EQUAL(factorMap["prior"]->mkFactor()(0,1), 0.2);
  
  BOOST_CHECK_EQUAL(factorMap["inner"]->mkFactor().size1(), (unsigned) 2);
  BOOST_CHECK_EQUAL(factorMap["inner"]->mkFactor().size2(), (unsigned) 2);
  BOOST_CHECK_EQUAL(factorMap["inner"]->mkFactor()(1,0), 0.1);
}

BOOST_AUTO_TEST_CASE(readFactorFile_2)
{
  //Check DiscContFactor
}

/*
BOOST_AUTO_TEST_CASE(readFactorFile_2) 
{
  map<string, AbsBasFacPtr_t> factorMap = readFactorFile("./data/dfgSpec/test2Potentials.txt");

  BOOST_CHECK(factorMap["prior"] != NULL);
  BOOST_CHECK_EQUAL(factorMap["prior"]->mkFactor().size1(), (unsigned) 1);
  BOOST_CHECK_EQUAL(factorMap["prior"]->mkFactor().size2(), (unsigned) 12);
  // BOOST_CHECK_EQUAL(factorMap["prior"]->mkFactor()(0,1), 0.2);
}
*/


BOOST_AUTO_TEST_CASE(writeFactorMap_1) 
{
  string const fileName1 = "./data/dfgSpec/test1Potentials.txt";
  string const fileName2 = "./output/test1Potentials.out.txt";

  map<string, AbsBasFacPtr_t> factorMap1 = readFactorFile(fileName1);
  writeFactorMap(fileName2, factorMap1);
  map<string, AbsBasFacPtr_t> factorMap2 = readFactorFile(fileName2);

  BOOST_CHECK( matrixEqual( factorMap1["prior"]->mkFactor(), factorMap2["prior"]->mkFactor() ) );
  BOOST_CHECK( matrixEqual( factorMap1["inner"]->mkFactor(), factorMap2["inner"]->mkFactor() ) );
}
