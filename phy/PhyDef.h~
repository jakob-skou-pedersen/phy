#ifndef __phyDef_h
#define __phyDef_h

#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "phy/xdoubleMod.h"

////////////////////////////////////////////////////////////////
//                  DEFINITION OF FLOAT TYPES
//
// The phy library defines two potentially different floating point
// types: number_t and xnumber_t. number_t must be a standard float
// type such as float or double. xnumber_t must have at least the same
// exponent range as number_t. If the calculations won't lead to
// underflow (or overflow) with long double, which often has exponent
// ranges of 1.0e-2466 to 1.0e2466, then that is the preferred
// xnumber_t type, since it is often implemented in hardware and thus
// fast. If larger exponent ranges are needed, then the the NTL
// xdouble type can be used (setting XNUMBER_IS_XDOUBLE), which stores
// the exponent in a long type (+ the original exponent range of
// double).
//
// Compilation constants that define what float types to use.  At most
// one should be uncommented (xnumber_t is a long double by default)

//#define XNUMBER_IS_NUMBER   // xnumber_t is the same type as number_t
#define XNUMBER_IS_XDOUBLE  // xnumber_t is a (modified) NTL::xdouble

namespace phy {

  namespace ublas = boost::numeric::ublas;
  using namespace std;

  /** Standard float type */
  typedef double number_t ;

  /** Float type with (potentially) extended exponent range */
#if   defined(XNUMBER_IS_NUMBER)
  typedef number_t xnumber_t ;
#elif defined(XNUMBER_IS_XDOUBLE)
  typedef NTL::xdoubleMod xnumber_t ;
#else
  typedef long double xnumber_t ;
#endif 

  /** Vector types holding number_t and xnumber_t types */
  typedef ublas::vector<number_t> vector_t;
  typedef ublas::vector<xnumber_t> xvector_t;
  
  /** Matrix types holding number_t and xnumber_t types */
  typedef ublas::matrix<number_t> matrix_t;
  typedef ublas::matrix<xnumber_t> xmatrix_t;
  
  /** Type used for observed data */
  typedef std::string symbol_t;

  /** Type used to represent states of random variables */
  typedef unsigned state_t;

  /** Type used for observed random variables */
  typedef ublas::vector<bool> stateMask_t;

  /** Type used for an enumerated set of stateMasks (e.g., input to factor graph) */
  typedef std::vector<stateMask_t const *> stateMaskVec_t;

  /** Type used for sequential data sets translated to stateMasks */
  typedef std::vector<stateMaskVec_t> stateMask2DVec_t;

  // conversion functions involving above defined types
  /** Convert extended-range-type (xnumber_t) to normal float type (number_t) */
#ifdef XNUMBER_IS_XDOUBLE
  inline void toNumber(number_t & r, xnumber_t const & a) {r = NTL::to_double(a);}
  inline number_t toNumber(xnumber_t const & a) {return NTL::to_double(a);}
#else
  inline void toNumber(number_t & r, xnumber_t const & a) {r = a;}
  inline number_t toNumber(xnumber_t const & a) {return a;}
#endif 

  inline void toNumber(vector_t & r, xvector_t const & v) {for (unsigned i = 0; i < v.size(); i++) r[i] = toNumber(v[i]);}
  inline vector_t toNumber(xvector_t const & v) {vector_t r(v.size()); toNumber(r, v); return r;}

  inline void toNumber(matrix_t & r, xmatrix_t const & m) {for (unsigned i = 0; i < m.size1(); i++) for (unsigned j = 0; j < m.size2(); j++) r(i, j) = toNumber( m(i, j) );}
  inline matrix_t toNumber(xmatrix_t const & m) {matrix_t r(m.size1(), m.size2()); toNumber(r, m); return r;}

  inline void toXNumber(xvector_t & r, vector_t const & v) {for (unsigned i = 0; i < v.size(); i++) r[i] = v[i];}
  inline xvector_t toXNumber(vector_t const & v) {xvector_t r(v.size()); toXNumber(r, v); return r;}

  inline void toXNumber(xmatrix_t & r, matrix_t const & m) {for (unsigned i = 0; i < m.size1(); i++) for (unsigned j = 0; j < m.size2(); j++) r(i, j) = m(i, j);}
  inline xmatrix_t toXNumber(matrix_t const & m) {xmatrix_t r(m.size1(), m.size2()); toXNumber(r, m); return r;}

  /** Convert between boost::numeric::ublas::vector and std::vector */
  template<class T>
  std::vector<T> toStdVector(ublas::vector<T> v) {std::vector<T> u(v.size()); for (unsigned i = 0; i < v.size(); i++) u[i] = v[i]; return u;}

  template<class T>
  ublas::vector<T> toNumVector(std::vector<T> v) {ublas::vector<T> u(v.size()); for (unsigned i = 0; i < v.size(); i++) u[i] = v[i]; return u;}




} // end namespace phy


#endif  //__phyDef_h
