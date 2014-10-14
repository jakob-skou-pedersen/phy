/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/FactorsIO.h"
#include "phy/Mixtures.h"
#include <stdio.h>
#include <cassert>

namespace phy {

  // helper functions
  AbsBasFacPtr_t readAbstractFullyParameterizedFactor(istream & str, string const & type, string const & name)
  {
    matrix_t potMat;
    getFeatureAndSkipLine(str, "POT_MAT:", potMat);

    matrix_t pcMat;
    if ( moreTags(str) ) {
      getFeatureAndSkipLine(str, "PC_MAT:", pcMat);
      if ( pcMat.size1() == 0 )
	errorAbort("From readAbstractFullyParameterizedFactor: PC_MAT must define a matrix of non-zero size or not be defined at all. Problem in specification of factor with name '" + name + "'.");
    }

    if (type == "rowNorm")
      return  AbsBasFacPtr_t(new RowNormFactor(name, potMat, pcMat) );
    if (type == "colNorm") 
      return  AbsBasFacPtr_t(new ColumnNormFactor(name, potMat, pcMat) );
    if (type == "globNorm")
      return  AbsBasFacPtr_t(new GlobalNormFactor(name, potMat, pcMat) );

    errorAbort("From readAbstractFullyParameterizedFactor: Unknown type ('" + type + "') in specification of factor with name '" + name + "'.");
    exit(1); // to satisfy compiler
  }

  AbsBasFacPtr_t readContContFactor(istream & str, string const & name){
    unsigned bins1, bins2;
    number_t min1v, max1v, min2v, max2v;
    getFeatureAndSkipLine(str, "BINS1:", bins1);
    getFeatureAndSkipLine(str, "BINS2:", bins2);
    getFeatureAndSkipLine(str, "MIN1:", min1v);
    getFeatureAndSkipLine(str, "MAX1:", max1v);
    getFeatureAndSkipLine(str, "MIN2:", min2v);
    getFeatureAndSkipLine(str, "MAX2:", max2v);
    
    return AbsBasFacPtr_t(new ContContFactor(name, bins1, bins2, min1v, max1v, min2v, max2v));
  }

  AbsBasFacPtr_t readDiscContFactor(istream & str, string const & name){
    int bins = 0;
    double minv = 0;
    double maxv = 0;
    int states = 0;
    string dist;

    getFeatureAndSkipLine(str, "BINS:", bins);
    getFeatureAndSkipLine(str, "MIN:", minv);
    getFeatureAndSkipLine(str, "MAX:", maxv);
    if(maxv < minv )
      errorAbort("readDiscContFactor: Normal factor max < min");

    getFeatureAndSkipLine(str, "STATES:", states);
    if( moreTags(str) ){
      getFeatureAndSkipLine(str, "DIST:", dist);
      MixPtr_t mixDist = readMixture(str, dist, minv, maxv, bins, states);
      return AbsBasFacPtr_t(new DiscContFactor(name, minv, maxv, states, bins, mixDist) );
    }
    else{
      //Use default(normal distribution)
      return AbsBasFacPtr_t(new DiscContFactor(name, minv, maxv, states, bins) );
    }
  }

  AbsBasFacPtr_t readBinomFactor(istream & str, string const & name){
    int minv, maxv;
    double prob = 0;
    int N = 0;
    getFeatureAndSkipLine(str, "MIN:", minv);
    getFeatureAndSkipLine(str, "MAX:", maxv);
    getFeatureAndSkipLine(str, "PROB:", prob);
    getFeatureAndSkipLine(str, "N:", N);

    return AbsBasFacPtr_t( new BinomialFactor(name, prob, N, minv, maxv) );
  }

  void writeAbstractFullyParameterizedFactor(ostream & str, AbsBasFacPtr_t const & factorPtr)
  {
    factorPtr->serialize(str);
  }

  // definition of header declared functions

  // read factor specifications. 
  // Reading of parameters and pseudo-counts done by helper functions
  pair<string, AbsBasFacPtr_t> readFactor(istream & str)
  {
    skipWhiteSpaceAndComments(str);

    string name;
    getFeatureAndSkipLine(str, "NAME:", name);

    string type;
    getFeatureAndSkipLine(str, "TYPE:", type);

    // dispatch the correct parser depending on type
    if (type == "rowNorm" or type == "colNorm" or type == "globNorm")
      return pair<string, AbsBasFacPtr_t>(name, readAbstractFullyParameterizedFactor(str, type, name) );
    if (type == "discCont")
      return pair<string, AbsBasFacPtr_t>(name, readDiscContFactor(str, name) );
    if (type == "contCont")
      return pair<string, AbsBasFacPtr_t>(name, readContContFactor(str, name) );
    if (type == "binom" or type == "binomial")
      return pair<string, AbsBasFacPtr_t>(name, readBinomFactor(str, name) );
    else {
      errorAbort("From readFactor: Unknown type ('" + type + "') in specification of factor with name '" + name + "'.");
      exit(1); // to satisfy compiler
    }
  }


  map<string, AbsBasFacPtr_t> readFactorFile(string const & file)
  {
    ifstream f;
    openInFile(f, file);
    return readFactorFile(f);
  }


  map<string, AbsBasFacPtr_t> readFactorFile(istream & str)
  {
    map<string, AbsBasFacPtr_t> factorMap;

    skipWhiteSpaceAndComments(str);
    while ( str.good() ) {
      while ( moreTags(str) ) {
	factorMap.insert( readFactor(str) );
	skipWhiteSpaceAndComments(str);
      }
    }

    return factorMap;
  }


  void writeFactor(ostream & str, AbsBasFacPtr_t const & factorPtr, string const & name)
  {
    str << "NAME:\t" << name << endl;
    string const type = factorPtr->type();
    str << "TYPE:\t" << type << endl;
    factorPtr->serialize(str);

    str << endl;
  }


  void writeFactorMap(string const & file, map<string, AbsBasFacPtr_t> const & factorMap)
  {
    ofstream f;
    f.open( file.c_str() );
    writeFactorMap(f, factorMap);
    f.close();
  }


  void writeFactorMap(ostream & str, map<string, AbsBasFacPtr_t> const & factorMap)
  {
    for (map<string, AbsBasFacPtr_t>::const_iterator iter = factorMap.begin(); 
	 iter != factorMap.end(); 
	 iter++)
      writeFactor(str, iter->second, iter->first);
  }


  void writeFactorVec(string const & file, vector<AbsBasFacPtr_t> const & facVec)
  {
    ofstream f;
    f.open( file.c_str() );
    writeFactorVec(f, facVec);
  }


  // helper function for writeFactorVec
  map<string, AbsBasFacPtr_t> facVecToFacMap(vector<AbsBasFacPtr_t> const & facVec)
  {
    map<string, AbsBasFacPtr_t> facMap;
    for (unsigned i = 0; i < facVec.size(); i++)
      facMap[ facVec[i]->name() ]  = facVec[i];
    return facMap;
  }


  void writeFactorVec(ostream & str, vector<AbsBasFacPtr_t> const & facVec)
  {
    writeFactorMap( str, facVecToFacMap(facVec) );
  }


  vector<AbsBasFacPtr_t> mkFactorVector(vector<string> const & facNames, map<string const, AbsBasFacPtr_t> const & factorMap)
  {
    vector<AbsBasFacPtr_t> facVec;
    map<string const, AbsBasFacPtr_t>::const_iterator iter;
    BOOST_FOREACH(string const & name, facNames) {
      iter = factorMap.find(name);
      if ( iter == factorMap.end() )
	errorAbort("From mkFactorVector: FactorMap holds no factor with name '" + name + "'.");
      facVec.push_back(iter->second);
    }
    return facVec;
  }


} // end namespace phy
