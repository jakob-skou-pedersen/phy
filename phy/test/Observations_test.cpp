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
// Note the hack to access protected data
#define private public
#define protected public

#include "phy/test/PhyTestDef.h"
#include "phy/Observations.h"
#include "phy/ObservationsIO.h"

#undef protected
#undef private 

// input/output not included by above
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

// other stuff
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

using namespace phy;
using boost::unit_test::test_suite;
using namespace std;

BOOST_AUTO_TEST_CASE(Observations_mkSymbol_1) 
{
  //  inline symbol_t mkSymbol(unsigned symSize, string const & obsString, long obsStringSize, long symStartPos, char missingDataChar = '.')
  string const obsString = "ACGT";
  BOOST_CHECK( mkSymbol(2, obsString, obsString.size(), -1, '.') == string(".A") );
  BOOST_CHECK( mkSymbol(2, obsString, obsString.size(), 0, '.') == string("AC") );
  BOOST_CHECK( mkSymbol(2, obsString, obsString.size(), 2, '.') == string("GT") );
  BOOST_CHECK( mkSymbol(2, obsString, obsString.size(), 4, '.') == string("..") );
  BOOST_CHECK( mkSymbol(1, obsString, obsString.size(), 1, '.') == string("C") );
}



BOOST_AUTO_TEST_CASE(ObservationsIO_overloadedInputOperator_1) 
{
  ifstream fIn("./data/metaNuc.txt");
  ofstream fOut("./output/metaNuc.txt");

  StateMap sm("");

  // nuc
  fIn >> sm;
  fOut << sm << endl;
  BOOST_CHECK(sm.name() == "nuc");
  BOOST_CHECK(sm.stateCount() == 4);

  // metaNuc
  fIn >> sm;
  fOut << sm << endl;
  BOOST_CHECK(sm.name() == "metaNuc");
  BOOST_CHECK(sm.stateCount() == 4);

  // 2-metaNuc
  fIn >> sm;
  fOut << sm << endl;
  BOOST_CHECK(sm.name() == "2-metaNuc");
  BOOST_CHECK(sm.stateCount() == 16);
}


BOOST_AUTO_TEST_CASE(ObservationsIO_stateMap_1) 
{
  //Discrete type statemap

  ifstream fIn("./data/dfgSpec/test1StateMaps.txt");
  //  ofstream fOut("./output/metaNuc.txt");

  StateMap sm("");
  vector<state_t> v(2);
  v.at(0) = 1; v.at(1) = 0;

  // twoState
  fIn >> sm;
  BOOST_CHECK(sm.name() == "twoState");
  BOOST_CHECK(sm.stateCount() == 2);
  BOOST_CHECK(sm.symbol2State("A") == 0);
  BOOST_CHECK(sm.state2Symbol(0) == "A");
  BOOST_CHECK(sm.symbol2State("B") == 1);
  BOOST_CHECK(sm.state2Symbol(1) == "B");

  BOOST_CHECK_EQUAL(sm.state2Symbol(v).size(), 2);
  BOOST_CHECK_EQUAL(sm.state2Symbol(v)[0], "B");
  BOOST_CHECK_EQUAL(sm.state2Symbol(v)[1], "A");

  BOOST_CHECK(sm.metaStateCount() == 2 + 2);
  BOOST_CHECK(sm.degeneracyVector("*")[0] == "A");
  BOOST_CHECK(sm.degeneracyVector("*")[1] == "B");



  // 2-twoState  (meaning four two-long states)
  fIn >> sm;
  BOOST_CHECK(sm.name() == "2-twoState");
  BOOST_CHECK(sm.stateCount() == 4);
  BOOST_CHECK(sm.symbol2State("AA") == 0);
  BOOST_CHECK(sm.symbol2State("AB") == 1);
  BOOST_CHECK(sm.symbol2State("BA") == 2);
  BOOST_CHECK(sm.symbol2State("BB") == 3);
  BOOST_CHECK(sm.metaStateCount() == 4 + 12);
  BOOST_CHECK(sm.degeneracyVector("**")[0] == "AA");
  BOOST_CHECK(sm.degeneracyVector("..")[1] == "AB");
}

BOOST_AUTO_TEST_CASE(ObservationsIO_stateMap_2) 
{
  //Continuous type statemap
  ifstream fIn("./data/dfgSpecNorm2/stateMaps.txt");
  
  StateMap sm("");
  fIn >> sm;
  
  //Testing Parsing
  BOOST_CHECK_EQUAL( sm.symbol2State("0.9999"), 199);
  BOOST_CHECK_EQUAL( sm.symbol2State("-0.9999"), 0);
  BOOST_CHECK_EQUAL( sm.symbol2State("0.004"),100);
  BOOST_CHECK_EQUAL( sm.symbol2State("0.014"),101);
}


BOOST_AUTO_TEST_CASE(ObservationsIO_stateMap_3)
{
  //Count type statemap
  ifstream fIn("./data/dfgSpecBinom2/stateMaps.txt");

  StateMap sm("");
  fIn >> sm;

  //Test parsing
  BOOST_CHECK_EQUAL( sm.symbol2State("5"), 0);
  BOOST_CHECK_EQUAL( sm.symbol2State("7"), 2);
  BOOST_CHECK_EQUAL( sm.symbol2State("15"), 10);
}

BOOST_AUTO_TEST_CASE(Observations_mkDiSymbol_1) 
{
  //  inline symbol_t mkDiSymbol(unsigned const sizeLeft, unsigned const sizeRight, string const & seq, long const seqSize, long const leftStartPos, long const rightStartPos, char const missingDataChar = '.')
  string const seq = "ACGT";
  BOOST_CHECK( mkDiSymbol(1, 1, seq, seq.size(), -1, 0, '.') == string(".A") );
  BOOST_CHECK( mkDiSymbol(1, 1, seq, seq.size(), 0, 3, '.') == string("AT") );
  BOOST_CHECK( mkDiSymbol(2, 1, seq, seq.size(), 0, 3, '.') == string("ACT") );
}


BOOST_AUTO_TEST_CASE(ObservationsIO_readAndWriteStateMap_1) 
{
  string const fileName1 = "./data/dfgSpec/test1StateMaps.txt";
  string const fileName2 = "./output/test1StateMaps.out.txt";

  map<string, StateMapPtr_t> smMapMap = readStateMapFile(fileName1);
  writeStateMapFile2(fileName2, smMapMap);
  map<string, StateMapPtr_t> smMapMap2 = readStateMapFile(fileName2);
  BOOST_CHECK( smMapMap["twoState"]->stateCount() == smMapMap2["twoState"]->stateCount() );
  BOOST_CHECK( smMapMap["twoState"]->metaStateCount() == smMapMap2["twoState"]->metaStateCount() );
  BOOST_CHECK( smMapMap["twoState"]->degeneracyVector("*") == smMapMap2["twoState"]->degeneracyVector("*") );
}


BOOST_AUTO_TEST_CASE(Observations_readSeqToVarSymbols_1) 
{
  string const fileName1 = "./data/seqToVarSymbolTest.txt";
  string const fileName2 = "./output/seqToVarSymbolTest.out.txt";

  vector<SeqToVarSymbol> v = readSeqToVarSymbols(fileName1);
  writeSeqToVarSymbols(fileName2, v);
  vector<SeqToVarSymbol> u = readSeqToVarSymbols(fileName2);
  BOOST_CHECK( v.size() == u.size() );

  SeqToVarSymbol const & stv1 = v[0];
  SeqToVarSymbol const & stv2 = u[0];
  SeqToVarSymbol const & stv3 = v[1];
  SeqToVarSymbol const & stv4 = v[4];

  BOOST_CHECK( stv1.varName == "H");
  BOOST_CHECK( stv1.seqName == "dna");
  BOOST_CHECK( stv1.symbolOffset == 2);
  BOOST_CHECK( stv1.indexOffset == 3);
  BOOST_CHECK( stv1.optional == true);
  BOOST_CHECK( stv1.varName == stv2.varName);
  BOOST_CHECK( stv1.varName != stv3.varName);
  BOOST_CHECK( stv1.size == stv2.size);
  BOOST_CHECK( stv1.rSize == stv2.rSize);

  BOOST_CHECK( stv3.varName == "O");
  BOOST_CHECK( stv3.seqName == "dna");
  BOOST_CHECK( stv3.rIndexOffset == -2);

  BOOST_CHECK( stv4.varName == "O4");
  BOOST_CHECK( stv4.seqName == "numbers");
  BOOST_CHECK( stv4.indexOffset == -1);
  BOOST_CHECK( stv4.rIndexOffset == 1);
}


BOOST_AUTO_TEST_CASE(Observations_stvToSeqMap_1) 
{
  // setup basic data structures
  string const fileName1 = "./data/seqToVarSymbolTest.txt";
  vector<SeqToVarSymbol> stvVec = readSeqToVarSymbols(fileName1);
  SeqData seqData;
  seqData.addSeq("numbers", "0123456789");
  seqData.addSeq("dna", "ACTGGACTGGACTGGACTGG");

  // map from stv to seq indexes. If sequence is missing for stvVec[i], then stvToSeqMap[i] = -1.
  vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames);
  for (unsigned i = 0; i < stvToSeqMap.size(); i++)
    BOOST_CHECK( stvToSeqMap[i] >= 0 );

  // test map
  vector<int> v( stvToSeqMap.size() );
  fromString(split("1,1,0,0,0", ","), v);
  BOOST_CHECK( stvToSeqMap == v );
}


BOOST_AUTO_TEST_CASE(Observations_mkStvToVarMap_1) 
{
  // setup basic data structures
  string const fileName1 = "./data/seqToVarSymbolTest.txt";
  vector<SeqToVarSymbol> stvVec = readSeqToVarSymbols(fileName1);

  vector<string> varNames = split("O,O3,O2,H,O4", ",");
  vector<unsigned> stvToVarMap = mkStvToVarMap(stvVec, varNames);

  vector<unsigned> v( stvToVarMap.size() );
  fromString(split("3,0,2,1,4", ","), v);

  BOOST_CHECK( stvToVarMap == v );
}


BOOST_AUTO_TEST_CASE(Observations_mkSymbol_2) 
{
  // setup basic data structures
  string const fileName1 = "./data/seqToVarSymbolTest.txt";
  vector<SeqToVarSymbol> stvVec = readSeqToVarSymbols(fileName1);
  SeqData seqData;
  seqData.addSeq("dna", "ACTGGACTGGACTGGACTGG");
  seqData.addSeq("numbers", "0123456789");

  // map from stv to seq indexes. If sequence is missing for stvVec[i], then stvToSeqMap[i] = -1.
  vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames);

  // get sequence sizes and check all have same symbol-lengths / symbol-offset sizes
  vector<long> seqSizes = getSeqSizes(seqData);
  sameSymbolCountOrAbort(stvVec, stvToSeqMap, seqSizes, seqData.entryId);

  // init symbol vector
  vector<string> symVec = initSymVec(stvVec, false);

  // make some "single" symbols
  mkSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 0, '.');
  vector<string> symVec_sng_idx0 = symVec;
  mkSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 1, '.');
  vector<string> symVec_sng_idx1 = symVec;
  mkSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 11, 'N');
  vector<string> symVec_sng_idx11 = symVec;
  // check expected values
  BOOST_CHECK( symVec_sng_idx0 == split("GG,CT,0,2,.", ",") );
  BOOST_CHECK( symVec_sng_idx1 == split("AC,GG,1,3,0", ",") );
  BOOST_CHECK( symVec_sng_idx11 == split("NN,NN,N,N,N", ",") );


  // make some "pair" symbols
  symVec = initSymVec(stvVec, true); // need to adjust sizes of symVec
  mkDiSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 0, 5,'.');
  vector<string> symVec_pair_idx0 = symVec;
  mkDiSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 1, 8,'.');
  vector<string> symVec_pair_idx1 = symVec;
  mkDiSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 11, 20,'.');
  vector<string> symVec_pair_idx11 = symVec;
  // check expected values

  //  cout << toString(symVec_pair_idx0) << endl;
  BOOST_CHECK( symVec_pair_idx0 == split("GGTG,CTG,05,234,.6", ",") );
  BOOST_CHECK( symVec_pair_idx1 == split("ACGG,GGG,18,367,09", ",") );
  BOOST_CHECK( symVec_pair_idx11 == split("....,...,..,...,..", ",") );
}

BOOST_AUTO_TEST_CASE(Observations_mkSymbol_3)
{
  // setup basic data structures
  string const fileName1 = "./data/seqToVarSymbolTestSub.txt";
  vector<SeqToVarSymbol> stvVec = readSeqToVarSymbols(fileName1);
  SeqData seqData;

  seqData.addSeq("X_2","00042,00065,00015");
  seqData.addSeq("X_1","00031,00050,00020");
  seqData.addSeq("N_1","00150,00160,00080");

  vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames);
  vector<long> seqSizes = getSeqSizes(seqData);
  
  // Check mapping to subvardata
  vector<string> subNames;
  subNames.push_back("N_1"); subNames.push_back("X_2");  subNames.push_back("X_1");
  vector<unsigned> stvToSubVarInvMap = mkStvToVarInvMap(stvVec, subNames);
  BOOST_CHECK( stvToSubVarInvMap.at(0) == 1);
  BOOST_CHECK( stvToSubVarInvMap.at(1) == 2);
  BOOST_CHECK( stvToSubVarInvMap.at(2) == 0);

  // Check symbol generation
  vector<string> symVec = initSymVec(stvVec, false);

  mkSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 0, '.');
  vector<string> symVec1 = symVec;
  mkSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 1, '.');
  vector<string> symVec2 = symVec;
  mkSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 2, '.');
  vector<string> symVec3 = symVec;

  BOOST_CHECK( symVec1.at(0) == "00031");
  BOOST_CHECK( symVec1.at(1) == "00150");
  BOOST_CHECK( symVec2.at(0) == "00050");
  BOOST_CHECK( symVec3.at(1) == "00080");

  BOOST_CHECK( boost::lexical_cast<int>(symVec1.at(0)) == 31);
  BOOST_CHECK( boost::lexical_cast<int>(symVec1.at(1)) == 150);
  BOOST_CHECK( boost::lexical_cast<int>(symVec2.at(0)) == 50);
  BOOST_CHECK( boost::lexical_cast<int>(symVec3.at(1)) == 80);

  // Check diSymbols
  symVec = initSymVec(stvVec, true); // adjust sizes of symVec
  mkDiSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 0, 1, '.');
  vector<string> symVec4 = symVec;
  mkDiSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, 0, 2, '.');
  vector<string> symVec5 = symVec;

  BOOST_CHECK( boost::lexical_cast<int>(symVec4.at(3)) == 50);
  BOOST_CHECK( boost::lexical_cast<int>(symVec5.at(3)) == 20);
}
