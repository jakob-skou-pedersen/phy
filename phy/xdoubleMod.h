/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __xdoubleMod_h
#define __xdoubleMod_h

#include <NTL/xdouble.h>
#include <algorithm>

// xdouble with double constructor, which is needed by ublas. 
// need to add some other constructors to allow all needed conversions
// to do: investigate performance of xdoubleMod compared to xdouble. Alternative is to add double constructor directly to xdouble header, however, this also have NTL version compatiblility issues.

namespace NTL {

  class xdoubleMod : public xdouble {
  public:
    
    inline xdoubleMod() : xdouble() {}
    inline xdoubleMod(xdouble const & a) : xdouble(a) {}
    inline xdoubleMod(double a) : xdouble( to_xdouble(a) ) {}
    inline xdoubleMod& operator=(double a) {*this = to_xdouble(a); return *this;}
  };

  /** fabs renamed to abs to conform with conventions (and ublas expectations). */
  inline xdouble abs(const xdouble& a) {return fabs(a);}

  /** max and min definitions for mixed xdouble and xdoubleMod arguments, since automatic conversion does not take place */
  template<class T>
  inline xdoubleMod max(T const & x, xdoubleMod const & y) {return std::max(static_cast<xdoubleMod>(x), y);}

  template<class T>
  inline xdoubleMod max(xdoubleMod const & x, T const & y) {return std::max(x, static_cast<xdoubleMod>(y) );}

  template<class T>
  inline xdoubleMod min(T const & x, xdoubleMod const & y) {return std::min(static_cast<xdoubleMod>(x), y);}

  template<class T>
  inline xdoubleMod min(xdoubleMod const & x, T const & y) {return std::min(x, static_cast<xdoubleMod>(y) );}

  inline xdoubleMod exp(xdoubleMod const & x) {return xexp( to_double(x) );}  // will fail if xdouble can be converted to double...


  ////////////////////////////////////////////////////////////////
  // inline function implementations

}


#endif  //__xdoubleMod_h
