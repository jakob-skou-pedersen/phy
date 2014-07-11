/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __StatDist_h
#define __StatDist_h

#include "phy/utils.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/unordered_map.hpp> 

namespace phy {

  using namespace std;

  /** return probability of the counts cv under the multinomial distribution defined by the parameters pv (should sum to one). */
  // note that the implementation uses the non-standard name power for the pow function
  template<typename T>
  T multinomialDistribution(vector<T> const & pv, vector<unsigned> const & cv);
  
  template<typename T>
  T multinomialCoefficient(vector<unsigned> const & cv);


  // helper functions
  double hashedLGamma(unsigned x);
  

  // template and inline function definitions
  template<typename T>
  T multinomialDistribution(vector<T> const & pv, vector<unsigned> const & cv)
  {
    assert(pv.size() == cv.size() );
    T p = multinomialCoefficient<T>(cv);
    for (unsigned i = 0; i < cv.size(); i++)
      p *= power(pv[i], cv[i]); 
    return p;
  }


  double hashedLGamma(unsigned x) 
  {
    static boost::unordered_map<unsigned, double> hashMap;
    boost::unordered_map<unsigned, double>::const_iterator iter = hashMap.find(x);
    if ( iter != hashMap.end() )
      return iter->second;
    else {
      double result = boost::math::lgamma(x);
      hashMap[x] = result;
      return result;
    }
  }
  

  template<typename T>
  T multinomialCoefficient(vector<unsigned> const & cv)
  {
    unsigned n = 0;
    for (unsigned i = 0; i < cv.size(); i++)
      n += cv[i];
    T d = 0;
    for (unsigned i = 0; i < cv.size(); i++)
      d += hashedLGamma(cv[i] + 1);
    return exp(hashedLGamma(n + 1) - d);
  }



} // end namespace phy

#endif  // __StatDist_h
