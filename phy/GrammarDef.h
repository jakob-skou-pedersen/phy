/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __GrammarDef_h
#define __GrammarDef_h

#include "phy/PhyDef.h"
#include "phy/xdoubleMod.h"

////////////////////////////////////////////////////////////////
// Compilation constants that define what float type to use
#define YNUMBER_IS_XNUMBER  // alternatively it is of type xdouble, although xnumber may not be

namespace phy {
  namespace grammar {

    // Definition of ynumber_t. Note that it must be possible to
    // convert from a xnumber_t to a ynumber_t without loss of
    // precision
    /** Float type used in Grammar calculations */
#ifdef YNUMBER_IS_XNUMBER
    typedef xnumber_t ynumber_t ;
#else
    typedef NTL::xdoubleMod ynumber_t ;
#endif

    // conversion from ynumber_t to number_t (only define if type of ynumber_t != xnumber_t)
#ifndef YNUMBER_IS_XNUMBER
    inline void toNumber(number_t & r, ynumber_t const & a) {r = NTL::to_double(a);}
    inline number_t toNumber(ynumber_t const & a) {return NTL::to_double(a);}

    inline void toNumber(number_t & r, ynumber_t const & a) {r = static_cast<double>(a);}
    inline number_t toNumber(ynumber_t const & a) {return static_cast<double>(a);}
#endif

  /** Matrix types holding number_t and xnumber_t types */
    typedef ublas::vector<ynumber_t> yvector_t;
    typedef ublas::matrix<ynumber_t> ymatrix_t;

    // conversion from number_t to ynumber_t (only define if type of ynumber_t != xnumber_t)
    template<class T>
    void toYNumber(yvector_t & r, ublas::vector<T> const & v) {for (unsigned i = 0; i < v.size(); i++) r[i] = v[i];}

    template<class T>
    yvector_t toYNumber(ublas::vector<T> const & v) {yvector_t r(v.size()); toYNumber(r, v); return r;}
    
    template<class T>
    void toYNumber(ymatrix_t & r, ublas::matrix<T> const & m) {for (unsigned i = 0; i < m.size1(); i++) for (unsigned j = 0; j < m.size2(); j++) r(i, j) = m(i, j);}

    template<class T>
    ymatrix_t toYNumber(ublas::matrix<T> const & m) {ymatrix_t r(m.size1(), m.size2()); toYNumber(r, m); return r;}

  } // end namespace grammar
} // end namespace phy


#endif  //__GrammarDef_h
