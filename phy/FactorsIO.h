/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __FactorsIO_h
#define __FactorsIO_h

#include "phy/Factors.h"
#include <boost/numeric/ublas/io.hpp>
#include <map>

/** IO functions for factors (in this use perhaps more aptly called
    potentials). */


namespace phy {

  using namespace std;
  
  // read functions
  pair<string, AbsBasFacPtr_t> readFactor(istream & str);
  map<string, AbsBasFacPtr_t> readFactorFile(string const & file);
  map<string, AbsBasFacPtr_t> readFactorFile(istream & str);

  // write functions
  void writeFactor(ostream & str, AbsBasFacPtr_t const & factorPtr, string const & name);
  void writeFactorMap(string const & file, map<string, AbsBasFacPtr_t> const & factorMap);
  void writeFactorMap(ostream & str, map<string, AbsBasFacPtr_t> const & factorMap);

  // note that the facNames should hold the potential names (may be different from the factor names in a given implementation). 
  // The facVec entries should correspond to these.
  void writeFactorVec(string const & file, vector<AbsBasFacPtr_t> const & facVec);
  void writeFactorVec(ostream & str, vector<AbsBasFacPtr_t> const & facVec);

  // conversion
  vector<AbsBasFacPtr_t> mkFactorVector(vector<string> const & facNames, map<string const, AbsBasFacPtr_t> const & factorMap);



} // end namespace phy

#endif  // __FactorsIO_h
