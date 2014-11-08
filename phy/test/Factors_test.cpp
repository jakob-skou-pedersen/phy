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
#include "phy/Mixtures.h"
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
}

BOOST_AUTO_TEST_CASE(DiscContFactor_general_2){
  unsigned bins = 10;
  double minv = 0;
  double maxv = 10;
  int states = 2;

  //Initialize means and variances
  vector_t means(states,5);
  vector_t vars(states, 4);

  MixPtr_t mixDist( new NormalMixture(means, vars, minv, maxv, bins));
  DiscContFactor dcf( "dcf", minv, maxv, states, bins, mixDist );

  //Test mk_factor
  matrix_t m(2,10);
  dcf.mkFactor(m);
  BOOST_CHECK_CLOSE( m(0,4) , 0.1914625, 0.0001);
}

BOOST_AUTO_TEST_CASE(DiscContFactor_general_3)
{
  vector_t means(1, 7);
  vector_t vars(1, 7);
  DiscContFactor dcf("noName", means, vars, 0, 2, 1, 10);
  
  //TODO Add some data and run init
  matrix_t counts(1,10,0);
  for(int i = 0; i < 10; ++i)
    counts(0,i) = -(i-3)*(i-3)+37;

  //  dcf.hasChanged_ = true;
  dcf.submitCounts( counts );
  dcf.optimizeParameters();

  /* R-code
     > bp <- c(0.6,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.4)
     > counts <- c((55),(60),(63),(64),(63),(60),(55),(48),(39),(28),(15))
     > mean(sweep(as.matrix(bp),1,counts,FUN=rep))
     [1] 1.347273
     > var(sweep(as.matrix(bp),1,counts,FUN=rep))
     [1] 0.290621
  */

  //Get our hands on the underlying NormalMixture
  boost::shared_ptr<NormalMixture> nm = boost::dynamic_pointer_cast<NormalMixture>( dcf.mixDist_ );

  //Test equality with some margin
  BOOST_CHECK_CLOSE( nm->means_(0) , 3.56603773585, 0.001);
  BOOST_CHECK_CLOSE( nm->vars_(0)  , 5.40566037736, 0.001);
}

BOOST_AUTO_TEST_CASE(DiscContFactor_update_1)
{
  vector_t alphas(1,1);
  vector_t betas(1,1);

  MixPtr_t mixDist( new BetaMixture(alphas, betas, 0, 1, 10));
  DiscContFactor dcf( "dcf", 0, 1, 2, 10, mixDist );

  //Test update
  vector<string> obs;
  obs.push_back( "10" );
  obs.push_back( "4" );
  dcf.setUpdateType(1);
  dcf.update( obs);
  matrix_t m(1,10,0);
  dcf.mkFactor(m);

  //R: > pbeta(0.3, 5,7) - pbeta(0.2,5,7)
  BOOST_CHECK_CLOSE( m(0,2), 0.1598950429, 0.0001);
}

BOOST_AUTO_TEST_CASE(ContContFactor_general_1)
{
  number_t alpha = 1;
  number_t beta = 0;
  number_t var = 1;
  unsigned bins = 10;
  number_t minv = 0;
  number_t maxv = 10;
  ContContFactor ccf("ccf", alpha, beta, var, bins, bins, minv, maxv, minv, maxv);

  matrix_t m = ccf.AbstractBaseFactor::mkFactor();
  BOOST_CHECK_CLOSE( m(4,4), 0.3413447, 0.0001);
  BOOST_CHECK_CLOSE( m(4,5), 0.3413447, 0.0001);
  
}

BOOST_AUTO_TEST_CASE(ContContFactor_reflection_1)
{
  map<string, AbsBasFacPtr_t> factorMap1 = readFactorFile("./data/dfgSpecFiles/factorPotentials.txt");
  writeFactorMap("./output/out.factorPotentials.txt", factorMap1);
  map<string, AbsBasFacPtr_t> factorMap2 = readFactorFile("./output/out.factorPotentials.txt");
  
  BOOST_CHECK( factorMap1["lin_reg"] != NULL);
  BOOST_CHECK( factorMap2["lin_reg"] != NULL);

  boost::shared_ptr<ContContFactor> ccf1 = boost::dynamic_pointer_cast<ContContFactor>( factorMap1["lin_reg"] );
  boost::shared_ptr<ContContFactor> ccf2 = boost::dynamic_pointer_cast<ContContFactor>( factorMap2["lin_reg"] );
  
  
  BOOST_CHECK_CLOSE( ccf1->alpha_, ccf2->alpha_, 0.0001);
  BOOST_CHECK_CLOSE( ccf1->beta_, ccf2->beta_, 0.0001);
  BOOST_CHECK_CLOSE( ccf1->var_, ccf2->var_, 0.0001);

  matrix_t pot1 = factorMap1["lin_reg"]->mkFactor();
  matrix_t pot2 = factorMap2["lin_reg"]->mkFactor();
  BOOST_CHECK_EQUAL( pot1.size1(), pot2.size1());
  BOOST_CHECK_EQUAL( pot1.size2(), pot2.size2());

  //Commented out for flooding output
  //for(int i = 0; i < pot1.size1(); ++i)
  //    for(int j = 0; j < pot1.size2(); ++j)
  //      BOOST_CHECK_CLOSE(pot1(i,j), pot2(i,j), 0.0001);
}

BOOST_AUTO_TEST_CASE(BinomialFactor_1)
{
  BinomialFactor bf("noname", 0.4, 5, 0, 5);
  matrix_t m(1,6);
  bf.mkFactor(m);

  //Test equalities
  BOOST_CHECK_CLOSE( m(0,0), 0.07776, 0.001);
  BOOST_CHECK_CLOSE( m(0,2), 0.34560, 0.001);

  //Test optimization
  matrix_t counts(1,6,0);
  counts(0,1) = 2;
  counts(0,2) = 2;
  bf.submitCounts( counts);
  bf.optimizeParameters();

  BOOST_CHECK_CLOSE( bf.prob_ , 6./20 , 0.001);

  //Test update
  vector<string> param;
  param.push_back( "10" );
  param.push_back( "0.2" );
  bf.setUpdateType(12);
  bf.update( param);
  bf.mkFactor(m);

  BOOST_CHECK_CLOSE( m(0,1), 0.268435456, 0.001);
  BOOST_CHECK_CLOSE( m(0,2), 0.301989888, 0.001);
}

BOOST_AUTO_TEST_CASE(NormalMeanPostFactor_1)
{
  //Construct
  number_t var = 2;
  number_t minv = -5;
  number_t maxv = 5;
  unsigned bins = 100;
  vector<string> subs;
  subs.push_back("X1");
  subs.push_back("X3");
  NormalMeanPostFactor nmpf("", var, minv, maxv, bins, subs);
  
  //Check subscriptions
  BOOST_CHECK_EQUAL( nmpf.getSubscriptions().size(), 2);
  BOOST_CHECK_EQUAL( nmpf.getSubscriptions().at(0), "X1");
  BOOST_CHECK_EQUAL( nmpf.getSubscriptions().at(1), "X3");

  //Submit variables
  vector<string> vars;
  vars.push_back("2.4");
  vars.push_back("1.8");
  nmpf.update( vars);
  BOOST_CHECK_CLOSE( nmpf.postmean_, 2.1, 0.0001);
  BOOST_CHECK_CLOSE( nmpf.postvar_, 1.0, 0.0001);

  //Check mkFactor
  matrix_t m(1,100);
  nmpf.mkFactor(m);
  //R> pnorm(0.1,2.1)-pnorm(0,2.1)
  BOOST_CHECK_CLOSE( m(0,50), 0.004885711, 0.0001);
  //R> pnorm(0.2,2.1)-pnorm(0.1,2.1)
  BOOST_CHECK_CLOSE( m(0,51), 0.005966428, 0.0001);

  //Resubmit variables
  vars.at(0) = "1.3";
  vars.at(1) = "0.9";
  nmpf.update(vars);

  //Check mkFactor again
  nmpf.mkFactor(m);
  //R> pnorm(0.1,1.1)-pnorm(0,1.1)
  BOOST_CHECK_CLOSE( m(0,50), 0.02298919, 0.0001);
  //R> pnorm(0.2,1.1)-pnorm(0.1,1.1)
  BOOST_CHECK_CLOSE( m(0,51), 0.02540487, 0.0001);
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
  map<string, AbsBasFacPtr_t> factorMap = readFactorFile("./data/dfgSpecBeta/factorPotentials.txt");
  
  BOOST_CHECK(factorMap["beta"] != NULL);

  //Check update
  BOOST_CHECK_EQUAL( factorMap["beta"]->getSubscriptions().at(0), "N1");
  BOOST_CHECK_EQUAL( factorMap["beta"]->getSubscriptions().at(1), "X1");
}

BOOST_AUTO_TEST_CASE(readFactorFile_3)
{
  //Check BinomialFactor
  map<string, AbsBasFacPtr_t> factorMap = readFactorFile("./data/dfgSpecBinom/factorPotentials.txt");

  //Check Existence
  BOOST_CHECK(factorMap["binomial"] != NULL);

  //Check dimensions
  BOOST_CHECK_EQUAL(factorMap["binomial"]->mkFactor().size1(), (unsigned) 1);
  BOOST_CHECK_EQUAL(factorMap["binomial"]->mkFactor().size2(), (unsigned) 6);

  //Check entries
  BOOST_CHECK_CLOSE(factorMap["binomial"]->mkFactor()(0,0), 0.07776, 0.000001);
  BOOST_CHECK_CLOSE(factorMap["binomial"]->mkFactor()(0,1), 0.2592, 0.000001);
  BOOST_CHECK_CLOSE(factorMap["binomial"]->mkFactor()(0,2), 0.3456, 0.000001);
}

BOOST_AUTO_TEST_CASE(readFactorFile_4)
{
  //Check BinomialFactor
  map<string, AbsBasFacPtr_t> factorMap = readFactorFile("./data/dfgSpecBinom2/factorPotentials.txt");

  //Check Existence
  BOOST_CHECK(factorMap["binomial"] != NULL);

  //Check dimensions
  BOOST_CHECK_EQUAL(factorMap["binomial"]->mkFactor().size1(), (unsigned) 1);
  BOOST_CHECK_EQUAL(factorMap["binomial"]->mkFactor().size2(), (unsigned) 11);

  //Check subscriptions
  BOOST_CHECK_EQUAL(factorMap["binomial"]->getSubscriptions().at(0), "Ns");
}

BOOST_AUTO_TEST_CASE(readFactorFile_5)
{
  //Check BinomialFactor
  map<string, AbsBasFacPtr_t> factorMap = readFactorFile("./data/dfgSpecSubscribe/factorPotentials.txt");

  //Check Existence
  BOOST_CHECK(factorMap["binomial1"] != NULL);
  BOOST_CHECK(factorMap["binomial2"] != NULL);


  //Check dimensions
  BOOST_CHECK_EQUAL(factorMap["binomial1"]->mkFactor().size1(), (unsigned) 1);
  BOOST_CHECK_EQUAL(factorMap["binomial1"]->mkFactor().size2(), (unsigned) 6);

  //Check subscriptions
  BOOST_CHECK_EQUAL(factorMap["binomial1"]->getSubscriptions().at(0), "N");
  BOOST_CHECK_EQUAL(factorMap["binomial1"]->getSubscriptions().at(1), "P1");
  BOOST_CHECK_EQUAL(factorMap["binomial2"]->getSubscriptions().at(0), "N");
  BOOST_CHECK_EQUAL(factorMap["binomial2"]->getSubscriptions().at(1), "P2");
}

BOOST_AUTO_TEST_CASE(readFactorFile_6){
  //Check NormalPostMeanFactor
  map<string, AbsBasFacPtr_t> factorMap = readFactorFile("./data/dfgSpecNormObs/factorPotentials.txt");

  //Check Existence
  BOOST_CHECK(factorMap["normObs"] != NULL);

  //Check dimensions
  BOOST_CHECK_EQUAL(factorMap["normObs"]->mkFactor().size1(), (unsigned) 1);
  BOOST_CHECK_EQUAL(factorMap["normObs"]->mkFactor().size2(), (unsigned) 100);

  //Check subscriptions
  BOOST_CHECK_EQUAL(factorMap["normObs"]->getSubscriptions().at(0), "X1");
  BOOST_CHECK_EQUAL(factorMap["normObs"]->getSubscriptions().at(1), "X3");
}

BOOST_AUTO_TEST_CASE(readFactorFile_7)
{
  //Check ContContFactor

}

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

BOOST_AUTO_TEST_CASE(writeFactorMap_2)
{
  //Check reflection of new continuous factors
  //If default parameters are used the output should be with parameters
}

