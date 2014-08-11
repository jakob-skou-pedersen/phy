/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/FactorsIO.h"
#include <stdio.h>
#include <cassert>

namespace phy {

  // helper functions
  AbsBasFacPtr_t readAbstractFullyParameterizedFactor(istream & str, string const & type, string const & name)
  {
    if(type == "rowNorm" or type == "colNorm" or type == "globNorm" ){
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

    }
    if (type == "normal"){
      //TODO Read all sorts of tags
      int no_bp = 0;
      getFeatureAndSkipLine(str, "BREAKPOINTS:", no_bp); std::cout << "Breakpoints read from file: " << no_bp << std::endl;

      double min = 0;
      double max = 0;
      getFeatureAndSkipLine(str, "MIN:", min);
      getFeatureAndSkipLine(str, "MAX:", max);
      printf("Range [%f;%f]\n", min, max);

      //TODO Assert min < max
      assert(min < max );

      vector_t breakpoints(no_bp);
      for(int i = 0; i < no_bp; ++i)
	breakpoints(i) = min + i* (max-min)/(no_bp-1);

      return AbsBasFacPtr_t(new NormalFactor(name, 1, 1, breakpoints));
    }

    errorAbort("From readAbstractFullyParameterizedFactor: Unknown type ('" + type + "') in specification of factor with name '" + name + "'.");
    exit(1); // to satisfy compiler
  }


  void writeAbstractFullyParameterizedFactor(ostream & str, AbsBasFacPtr_t const & factorPtr)
  {
    AbstractFullyParameterizedFactor * derivedFctorPtr = dynamic_cast<AbstractFullyParameterizedFactor *>( factorPtr.get() );

    str << "POT_MAT:\t" << derivedFctorPtr->m_ << endl;
    if ( derivedFctorPtr->pseudoCounts_.size1() != 0 )
      str << "PC_MAT:\t"  << derivedFctorPtr->pseudoCounts_ << endl;
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
    if (type == "rowNorm" or type == "colNorm" or type == "globNorm" or type == "normal")
      return pair<string, AbsBasFacPtr_t>(name, readAbstractFullyParameterizedFactor(str, type, name) );
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

    if (type == "rowNorm" or type == "colNorm" or type == "globNorm")
      writeAbstractFullyParameterizedFactor(str, factorPtr);
    else 
      errorAbort("From writeFactor: Unknown type ('" + type + "') in write request for factor with name '" + name + "'.");
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
