/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/FactorsIO.h"
#include "phy/Mixtures.h"
#include "boost/lexical_cast.hpp"
#include <stdio.h>
#include <cassert>
#include <sstream>

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
    number_t minv1, maxv1, minv2, maxv2;
    number_t alpha, beta, var;
    getFeatureAndSkipLine(str, "BINS1:", bins1);
    getFeatureAndSkipLine(str, "BINS2:", bins2);
    getFeatureAndSkipLine(str, "MIN1:", minv1);
    getFeatureAndSkipLine(str, "MAX1:", maxv1);
    getFeatureAndSkipLine(str, "MIN2:", minv2);
    getFeatureAndSkipLine(str, "MAX2:", maxv2);
    
    if( moreTags(str)){
      getFeatureAndSkipLine(str, "ALPHA:", alpha);
      getFeatureAndSkipLine(str, "BETA:", beta);
      getFeatureAndSkipLine(str, "VAR:", var);
      
      //Convert alphas and beta to "index-space"
      //Reflection: See serialize method in ContContFactor
      ConvertIndexNumber is1(minv1, maxv1, bins1);
      ConvertIndexNumber is2(minv2, maxv2, bins2);
      alpha = alpha*is1.getToNumberAlpha()/is2.getToNumberAlpha();
      beta =  (beta - is2.getToNumberBeta() + alpha*is2.getToNumberAlpha()*is1.getToNumberBeta()/is1.getToNumberAlpha())/is2.getToNumberAlpha();
      var = var/is2.getToNumberAlpha()/is2.getToNumberAlpha();
      return AbsBasFacPtr_t(new ContContFactor(name, alpha, beta, var, bins1, bins2, minv1, maxv1, minv2, maxv2));
    }
    return AbsBasFacPtr_t(new ContContFactor(name, bins1, bins2, minv1, maxv1, minv2, maxv2));
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
    string sProb, sN;
    bool probSubscribe, NSubscribe; 
    double prob = 0;
    int N = 0;
    getFeatureAndSkipLine(str, "MIN:", minv);
    getFeatureAndSkipLine(str, "MAX:", maxv);
    getFeatureAndSkipLine(str, "PROB:", sProb);
    getFeatureAndSkipLine(str, "N:", sN);

    //Check if prob and N are references to variables and set factor subscriptions
    try{
      prob = boost::lexical_cast<double>( sProb);
      probSubscribe = false;
    }
    catch(boost::bad_lexical_cast &){
      prob = 0.5;
      probSubscribe = true;
    }

    try{
      N = boost::lexical_cast<int>( sN);
      NSubscribe = false;
    }
    catch(boost::bad_lexical_cast &){
      N = 1;
      NSubscribe = true;
    }
    
    AbsBasFacPtr_t bf(new BinomialFactor(name, prob, N, minv, maxv));
    if( probSubscribe && NSubscribe){
      bf->setUpdateType(12);
      bf->addSubscription( sN );
      bf->addSubscription( sProb );      
    }
    else if( NSubscribe){
      bf->setUpdateType(1);
      bf->addSubscription( sN );
    }
    else if( probSubscribe){
      bf->setUpdateType(2);
      bf->addSubscription( sProb );
    }

    return bf;
  }

  AbsBasFacPtr_t readNormalMeanPostFactor(istream & str, string const & name){
    number_t minv, maxv;
    number_t var;
    unsigned bins;
    vector<string> subs;

    getFeatureAndSkipLine(str, "BINS:", bins);
    getFeatureAndSkipLine(str, "MIN:", minv);
    getFeatureAndSkipLine(str, "MAX:", maxv);
    getFeatureAndSkipLine(str, "VAR:", var);
    
    string tag;
    str >> tag;
    if(tag != "X:")
      errorAbort("FactorsIO.cpp::readNormalMeanPostFactor: Next tag should be X:");
    
    string sub_string;
    getline(str, sub_string);
    std::stringstream ss(sub_string);
    istream_iterator<string> begin(ss);
    istream_iterator<string> end;
    subs = vector<string>(begin,end);

    return AbsBasFacPtr_t(new NormalMeanPostFactor(name, var, minv, maxv, bins, subs));
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
    if (type == "normalMeanPost")
      return pair<string, AbsBasFacPtr_t>(name, readNormalMeanPostFactor(str, name) );
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
