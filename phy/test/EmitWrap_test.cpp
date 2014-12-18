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



BOOST_AUTO_TEST_CASE(EmitWrap_EmitWrapDfg_1) 
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
  // cout << matEmit[1] << endl;
  BOOST_CHECK_CLOSE(toNumber( matEmit[1](0,0) ), 0.0, EPS);     // j-i=0<minDist
  //  BOOST_CHECK_CLOSE(toNumber( matEmit[1](0,4) ), 1.0, EPS);     // j-i=5>=minDist  // this fails as stack factor potential needs adjustment so rows sum to one  - sudhakar, please fix
  //  BOOST_CHECK_CLOSE(toNumber( matEmit[1](0,14) ), 0.225, EPS);  //?? I have not done the calculation...
  //  BOOST_CHECK_CLOSE(toNumber( matEmit[1](12,16) ), 0.05, EPS);  
}

BOOST_AUTO_TEST_CASE(EmitWrap_EmitWrapDfg_2)
{
  //Read discrete factor graph
  string const inPrefix = "./data/probfold/dfgSpec/singleEmit/";
  
  string const stateMapsFile           = "stateMaps.txt";
  string const factorPotentialsFile    = "factorPotentials.txt";
  string const variablesFile           = "variables.txt";
  string const factorGraphFile         = "factorGraph.txt";
    
  DfgInfo dfgInfo = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile, inPrefix + variablesFile, inPrefix + factorGraphFile);

  //Make some seq data
  SeqData seqData;
  seqData.addSeq("X_1", "00020,00015,00001");
  seqData.addSeq("N_1", "00050,00025,00003");
  seqData.addSeq("P_1", "0.305,0.508,0.305");

  //Read stv
  string const stvFile = "seqToVarSymbol.txt";
  vector<SeqToVarSymbol> stvVec = readSeqToVarSymbols(inPrefix + stvFile);

  xvector_t res = calcLikelihoodVector(seqData, dfgInfo, stvVec, '.');

  //R:> pbeta(0.310,21,31)-pbeta(0.305,21,31)
  BOOST_CHECK_CLOSE( 0.01120509, toNumber(res(0)) , 0.0001);
  //R:> pbeta(0.510,16,11)-pbeta(0.505,16,11)
  BOOST_CHECK_CLOSE( 0.01361204, toNumber(res(1)), 0.0001);
  //R:> pbeta(0.310,2,3)-pbeta(0.305,2,3)
  BOOST_CHECK_CLOSE( 0.008847678, toNumber(res(2)), 0.0001);

}

BOOST_AUTO_TEST_CASE(EmitWrap_EmitWrapDfg_3){
  //Read two discrete factor graphs
  string const inPrefix = "./data/dfgSpecSahoo/";

  string const stateMapsFile           = "stateMaps.txt";
  string const factorPotentialsFile1   = "factorPotentials1.txt";
  string const factorPotentialsFile2   = "factorPotentials2.txt";
  string const variablesFile           = "variables.txt";
  string const factorGraphFile         = "factorGraph.txt";

  DfgInfo dfgInfo1 = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile1, inPrefix + variablesFile, inPrefix + factorGraphFile);
  DfgInfo dfgInfo2 = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile2, inPrefix + variablesFile, inPrefix + factorGraphFile);

  //Make some seq data
  SeqData seqData;
  seqData.addSeq("X_1", "00020,00015,00001");
  seqData.addSeq("N_1", "00050,00025,00003");

  //Read stv
  string const stvFile = "seqToVarSymbol.txt";
  vector<SeqToVarSymbol> stvVec = readSeqToVarSymbols(inPrefix + stvFile);
  
  //Calculate likelihood in both models
  xvector_t res1 = calcLikelihoodVector(seqData, dfgInfo1, stvVec, '.');
  xvector_t res2 = calcLikelihoodVector(seqData, dfgInfo2, stvVec, '.');

  //R:> library(VGAM)
  //R:> dbetabinom.ab(20,50,2,3)/dbetabinom.ab(20,50,8,12)
  BOOST_CHECK_CLOSE(0.5429574, toNumber(res1(0)/res2(0)) , 0.01);
  //R:> dbetabinom.ab(15,25,2,3)/dbetabinom.ab(15,25,8,12)
  BOOST_CHECK_CLOSE(1.041309, toNumber(res1(1)/res2(1)), 0.01);
  //R:> dbetabinom.ab(1,3,2,3)/dbetabinom.ab(1,3,8,12)
  BOOST_CHECK_CLOSE(0.8461538, toNumber(res1(2)/res2(2)), 0.01);
}

BOOST_AUTO_TEST_CASE(EmitWrap_EmitWrapDfg_4){
  //Testing pair emit and the likelihood matrix
  string const inPrefix = "./data/dfgSpecSahooPair/";

  string const stateMapsFile           = "stateMaps.txt";
  string const factorPotentialsFile    = "factorPotentials.txt";
  string const variablesFile           = "variables.txt";
  string const factorGraphFile         = "factorGraph.txt";
  
  DfgInfo dfgInfo = readDfgInfo(inPrefix + stateMapsFile, inPrefix + factorPotentialsFile, inPrefix + variablesFile, inPrefix + factorGraphFile);

  //Make some seq data
  SeqData seqData;
  seqData.addSeq("X_1", "00020,00015,00001");
  seqData.addSeq("N_1", "00050,00025,00003");

  //Read stv
  string const stvFile = "seqToVarSymbol.txt";
  vector<SeqToVarSymbol> stvVec = readSeqToVarSymbols(inPrefix + stvFile);

  //Calculate likelihood matrix
  xmatrix_t res = calcLikelihoodMatrix(seqData, dfgInfo, stvVec);
  std::cout << res << std::endl;
}


